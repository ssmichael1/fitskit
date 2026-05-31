//! Tiled image compression (the BINTABLE-based compressed-image convention).
//!
//! A tile-compressed image is stored physically as a `BINTABLE` extension but is
//! logically an image. The original image is divided into a rectangular grid of
//! *tiles*; each tile is compressed independently and stored as a byte stream in a
//! variable-length-array column (typically `COMPRESSED_DATA`). One table row holds
//! one tile, in row-major (axis-1-fastest) order.
//!
//! See the *Tiled Image Convention for Storing Compressed Images in FITS Binary
//! Tables* (v2.3) and the cfitsio implementation.
//!
//! # Status
//!
//! This module is being built up in phases (see `COMPRESSION_PLAN.md`):
//!
//! - Phase 0 (current): scaffolding + the [`rice_decompress_i32`] family of
//!   decoders with unit tests.
//! - Phase 1: RICE_1 integer read path ([`CompressedImage::decompress`]).
//! - Later: GZIP_1/GZIP_2, float quantization/dithering, PLIO_1, HCOMPRESS_1, and a
//!   compression (write) path.
//!
//! Most of the high-level pipeline is currently stubbed with `todo!()`.

use crate::bintable::BinTable;
use crate::error::{Error, Result};
use crate::header::Header;

/// Compression algorithm named by the `ZCMPTYPE` keyword.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionType {
    /// `RICE_1` (also `RICE_ONE`) — Rice coding of pixel differences. Most common.
    Rice1,
    /// `GZIP_1` — DEFLATE over the raw tile bytes.
    Gzip1,
    /// `GZIP_2` — DEFLATE over byte-shuffled tile bytes.
    Gzip2,
    /// `HCOMPRESS_1` — H-transform + quantization (lossy/lossless).
    Hcompress1,
    /// `PLIO_1` — IRAF pixel-list run-length encoding (integer masks).
    Plio1,
    /// `NOCOMPRESS` — tiles stored uncompressed.
    NoCompress,
}

impl CompressionType {
    /// Parse a `ZCMPTYPE` keyword value.
    pub fn parse(s: &str) -> Result<Self> {
        match s.trim().to_ascii_uppercase().as_str() {
            "RICE_1" | "RICE_ONE" => Ok(CompressionType::Rice1),
            "GZIP_1" => Ok(CompressionType::Gzip1),
            "GZIP_2" => Ok(CompressionType::Gzip2),
            "HCOMPRESS_1" => Ok(CompressionType::Hcompress1),
            "PLIO_1" => Ok(CompressionType::Plio1),
            "NOCOMPRESS" => Ok(CompressionType::NoCompress),
            other => Err(Error::UnsupportedCompression(other.to_string())),
        }
    }
}

/// Float quantization / dithering method named by `ZQUANTIZ`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantize {
    /// No dithering (or keyword absent).
    None,
    /// `SUBTRACTIVE_DITHER_1`.
    SubtractiveDither1,
    /// `SUBTRACTIVE_DITHER_2` (preserves exact zeros).
    SubtractiveDither2,
}

impl Quantize {
    /// Parse a `ZQUANTIZ` keyword value (absent ⇒ [`Quantize::None`]).
    pub fn parse(s: Option<&str>) -> Result<Self> {
        match s.map(|v| v.trim().to_ascii_uppercase()) {
            None => Ok(Quantize::None),
            Some(v) => match v.as_str() {
                "NO_DITHER" => Ok(Quantize::None),
                "SUBTRACTIVE_DITHER_1" => Ok(Quantize::SubtractiveDither1),
                "SUBTRACTIVE_DITHER_2" => Ok(Quantize::SubtractiveDither2),
                other => Err(Error::UnsupportedCompression(format!(
                    "ZQUANTIZ={other}"
                ))),
            },
        }
    }
}

/// Tile-grid geometry derived from `ZNAXIS`/`ZNAXISn`/`ZTILEn`.
#[derive(Debug, Clone)]
pub struct TileGeometry {
    /// Original image dimensions (`ZNAXISn`), axis-1 first.
    pub znaxis: Vec<usize>,
    /// Tile dimensions (`ZTILEn`), axis-1 first.
    pub ztile: Vec<usize>,
}

impl TileGeometry {
    /// Build geometry from the compressed-image header.
    ///
    /// Defaults to row-by-row tiling (`ZTILE1 = ZNAXIS1`, others = 1) when the
    /// `ZTILEn` keywords are absent.
    pub fn from_header(header: &Header) -> Result<Self> {
        let znaxis = header.require_int("ZNAXIS")? as usize;
        let mut dims = Vec::with_capacity(znaxis);
        for i in 1..=znaxis {
            dims.push(header.require_int(&format!("ZNAXIS{i}"))? as usize);
        }

        let mut tile = Vec::with_capacity(znaxis);
        for i in 1..=znaxis {
            let default = if i == 1 { dims[0] } else { 1 };
            let t = header
                .get_int(&format!("ZTILE{i}"))
                .map(|v| v as usize)
                .unwrap_or(default);
            tile.push(t);
        }

        Ok(TileGeometry {
            znaxis: dims,
            ztile: tile,
        })
    }

    /// Number of tiles along each axis.
    pub fn tiles_per_axis(&self) -> Vec<usize> {
        self.znaxis
            .iter()
            .zip(&self.ztile)
            .map(|(&n, &t)| n.div_ceil(t.max(1)))
            .collect()
    }

    /// Total number of tiles (== expected number of table rows).
    pub fn num_tiles(&self) -> usize {
        self.tiles_per_axis().iter().product()
    }
}

/// A compressed image read from a BINTABLE, ready to be decompressed.
///
/// Borrows the already-parsed [`BinTable`]; the heavy lifting happens in
/// [`CompressedImage::decompress`].
#[derive(Debug)]
pub struct CompressedImage<'a> {
    /// Compression algorithm.
    pub ctype: CompressionType,
    /// BITPIX of the original (uncompressed) image.
    pub zbitpix: i64,
    /// Tile-grid geometry.
    pub geometry: TileGeometry,
    /// Float quantization/dither method.
    pub quantize: Quantize,
    /// Dither seed (`ZDITHER0`), if present.
    pub zdither0: Option<i64>,
    /// Integer blank sentinel (`ZBLANK` keyword), if present.
    pub blank: Option<i64>,
    /// RICE `BLOCKSIZE` parameter (default 32).
    pub blocksize: usize,
    /// RICE `BYTEPIX` parameter (default 4).
    pub bytepix: usize,
    /// The borrowed binary table holding the compressed tiles + heap.
    /// Read once [`CompressedImage::decompress`] is implemented (phase 1).
    #[allow(dead_code)]
    pub(crate) table: &'a BinTable,
}

impl<'a> CompressedImage<'a> {
    /// True if `header` describes a tile-compressed image (`ZIMAGE = T`).
    pub fn detect(header: &Header) -> bool {
        header.get_bool("ZIMAGE").unwrap_or(false)
    }

    /// Build a [`CompressedImage`] from a compressed-image BINTABLE header and the
    /// parsed table. Parses the `Z*` driver keywords; does not yet decompress.
    pub fn from_bintable(header: &Header, table: &'a BinTable) -> Result<Self> {
        if !Self::detect(header) {
            return Err(Error::CompressionError(
                "header is not a tile-compressed image (ZIMAGE != T)".into(),
            ));
        }

        let ctype = CompressionType::parse(
            header
                .get_string("ZCMPTYPE")
                .ok_or_else(|| Error::MissingKeyword("ZCMPTYPE".into()))?,
        )?;
        let zbitpix = header.require_int("ZBITPIX")?;
        let geometry = TileGeometry::from_header(header)?;
        let quantize = Quantize::parse(header.get_string("ZQUANTIZ"))?;
        let zdither0 = header.get_int("ZDITHER0");
        let blank = header.get_int("ZBLANK");

        // Default RICE parameters; overridden by ZNAMEi='BLOCKSIZE'/'BYTEPIX'.
        let mut blocksize = 32usize;
        let mut bytepix = 4usize;
        let nzparams = count_zname_params(header);
        for i in 1..=nzparams {
            if let Some(name) = header.get_string(&format!("ZNAME{i}")) {
                let val = header.get_int(&format!("ZVAL{i}")).unwrap_or(0) as usize;
                match name.trim().to_ascii_uppercase().as_str() {
                    "BLOCKSIZE" => blocksize = val,
                    "BYTEPIX" => bytepix = val,
                    _ => {}
                }
            }
        }

        Ok(CompressedImage {
            ctype,
            zbitpix,
            geometry,
            quantize,
            zdither0,
            blank,
            blocksize,
            bytepix,
            table,
        })
    }

    /// The compression algorithm (`ZCMPTYPE`). Cheap inspector.
    pub fn compression(&self) -> CompressionType {
        self.ctype
    }

    /// The tile-grid geometry (`ZNAXISn`/`ZTILEn`). Cheap inspector.
    pub fn geometry(&self) -> &TileGeometry {
        &self.geometry
    }

    /// Decompress all tiles and reassemble the full image into an [`ImageData`].
    ///
    /// Implemented for integer `ZBITPIX` (8/16/32/64) with `RICE_1` (and, behind the
    /// `gzip` feature, `GZIP_1`/`GZIP_2`). Each tile is read from the
    /// `COMPRESSED_DATA` variable-length-array column, decoded to a flat array of
    /// integers, and scattered into the correct sub-rectangle of the output buffer
    /// per the tile grid (axis-1 fastest). `ZBLANK` integer sentinels are carried
    /// through unchanged into the output array.
    ///
    /// # Not yet supported (see `COMPRESSION_PLAN.md`)
    ///
    /// - Float `ZBITPIX` (-32/-64): requires quantization inverse (`ZSCALE`/`ZZERO`)
    ///   and subtractive dithering — phase 3.
    /// - `PLIO_1` — phase 4.
    /// - `HCOMPRESS_1` — phase 5.
    /// - `GZIP_1`/`GZIP_2` without the `gzip` feature return
    ///   [`Error::UnsupportedCompression`].
    pub fn decompress(&self) -> Result<crate::image_data::ImageData> {
        use crate::image_data::{ImageData, PixelData};
        use crate::types::Bitpix;

        // Float originals require quantization inverse + dithering (phase 3).
        if self.zbitpix < 0 {
            return Err(Error::UnsupportedCompression(format!(
                "float ZBITPIX={} (quantization/dithering decode is phase 3, see COMPRESSION_PLAN.md)",
                self.zbitpix
            )));
        }

        let bitpix = Bitpix::from_i64(self.zbitpix)?;
        let npix: usize = self.geometry.znaxis.iter().product();

        // Decode every tile into i64 values, in tile order (axis-1 fastest), then
        // scatter into the full-image buffer.
        let mut full = vec![0i64; npix];
        let tiles_per_axis = self.geometry.tiles_per_axis();
        let num_tiles = self.geometry.num_tiles();

        if self.table.nrows < num_tiles {
            return Err(Error::CompressionError(format!(
                "compressed image expects {} tiles but BINTABLE has {} rows",
                num_tiles, self.table.nrows
            )));
        }

        let cdata_col = self.compressed_data_column()?;

        for tile_index in 0..num_tiles {
            let raw = self.tile_bytes(tile_index, cdata_col)?;
            // Number of pixels in this (possibly edge-truncated) tile.
            let coords = unravel(tile_index, &tiles_per_axis);
            let tile_dims = self.tile_dims_at(&coords);
            let tile_npix: usize = tile_dims.iter().product();

            let values = self.decode_tile(&raw, tile_npix)?;
            if values.len() != tile_npix {
                return Err(Error::CompressionError(format!(
                    "tile {tile_index} decoded {} values, expected {tile_npix}",
                    values.len()
                )));
            }
            scatter_tile(
                &mut full,
                &self.geometry.znaxis,
                &self.geometry.ztile,
                &tile_dims,
                &coords,
                &values,
            );
        }

        // Narrow i64 -> the target integer storage.
        let pixels = match bitpix {
            Bitpix::U8 => PixelData::U8(full.iter().map(|&v| v as u8).collect()),
            Bitpix::I16 => PixelData::I16(full.iter().map(|&v| v as i16).collect()),
            Bitpix::I32 => PixelData::I32(full.iter().map(|&v| v as i32).collect()),
            Bitpix::I64 => PixelData::I64(full),
            Bitpix::F32 | Bitpix::F64 => unreachable!("float handled above"),
        };

        Ok(ImageData::new(self.geometry.znaxis.clone(), pixels))
    }

    /// Resolve the index of the `COMPRESSED_DATA` column.
    fn compressed_data_column(&self) -> Result<usize> {
        self.table
            .columns
            .iter()
            .position(|c| c.name.trim().eq_ignore_ascii_case("COMPRESSED_DATA"))
            .ok_or_else(|| {
                Error::CompressionError("compressed image has no COMPRESSED_DATA column".into())
            })
    }

    /// Read the raw compressed bytes for one tile (one table row) from the VLA column.
    fn tile_bytes(&self, row: usize, col: usize) -> Result<Vec<u8>> {
        match self.table.get_cell(row, col)? {
            crate::bintable::BinCellValue::Bytes(b) => Ok(b),
            other => Err(Error::CompressionError(format!(
                "COMPRESSED_DATA cell is not a byte VLA: {other:?}"
            ))),
        }
    }

    /// Tile dimensions at tile-grid coordinates `coords`, accounting for edge tiles
    /// that are truncated when `ZNAXISn` is not a multiple of `ZTILEn`.
    fn tile_dims_at(&self, coords: &[usize]) -> Vec<usize> {
        let mut dims = Vec::with_capacity(coords.len());
        for (axis, &c) in coords.iter().enumerate() {
            let n = self.geometry.znaxis[axis];
            let t = self.geometry.ztile[axis].max(1);
            let start = c * t;
            dims.push((n - start).min(t));
        }
        dims
    }

    /// Decode one tile's bytes into a flat array of `tile_npix` integer values
    /// (axis-1 fastest) according to the compression type.
    fn decode_tile(&self, raw: &[u8], tile_npix: usize) -> Result<Vec<i64>> {
        match self.ctype {
            CompressionType::Rice1 => match self.bytepix {
                1 => Ok(rice_decompress_i8(raw, tile_npix, self.blocksize)?
                    .into_iter()
                    .map(|v| v as i64)
                    .collect()),
                2 => Ok(rice_decompress_i16(raw, tile_npix, self.blocksize)?
                    .into_iter()
                    .map(|v| v as i64)
                    .collect()),
                4 => Ok(rice_decompress_i32(raw, tile_npix, self.blocksize)?
                    .into_iter()
                    .map(|v| v as i64)
                    .collect()),
                other => Err(Error::CompressionError(format!(
                    "RICE_1 BYTEPIX={other} not supported (expected 1, 2, or 4)"
                ))),
            },
            CompressionType::NoCompress => self.bytes_to_ints(raw, tile_npix),
            CompressionType::Gzip1 => {
                let inflated = gzip_inflate(raw)?;
                self.bytes_to_ints(&inflated, tile_npix)
            }
            CompressionType::Gzip2 => {
                let inflated = gzip_inflate(raw)?;
                let unshuffled = gzip2_unshuffle(&inflated, self.zbitpix_bytes());
                self.bytes_to_ints(&unshuffled, tile_npix)
            }
            CompressionType::Plio1 => Err(Error::UnsupportedCompression(
                "PLIO_1 decode is phase 4 (see COMPRESSION_PLAN.md)".into(),
            )),
            CompressionType::Hcompress1 => Err(Error::UnsupportedCompression(
                "HCOMPRESS_1 decode is phase 5 (see COMPRESSION_PLAN.md)".into(),
            )),
        }
    }

    /// Bytes-per-pixel of the original (uncompressed) integer image.
    fn zbitpix_bytes(&self) -> usize {
        (self.zbitpix.unsigned_abs() as usize) / 8
    }

    /// Interpret a big-endian byte buffer as `tile_npix` integers of width
    /// `ZBITPIX` (used by NOCOMPRESS and the GZIP codecs, which store the raw image
    /// integers, not RICE-coded differences).
    fn bytes_to_ints(&self, bytes: &[u8], tile_npix: usize) -> Result<Vec<i64>> {
        let width = self.zbitpix_bytes();
        if width == 0 {
            return Err(Error::CompressionError("invalid ZBITPIX width".into()));
        }
        let avail = bytes.len() / width;
        if avail < tile_npix {
            return Err(Error::CompressionError(format!(
                "tile has {avail} integers, expected {tile_npix}"
            )));
        }
        let mut out = Vec::with_capacity(tile_npix);
        for c in bytes.chunks_exact(width).take(tile_npix) {
            let v = match width {
                1 => c[0] as i64, // ZBITPIX=8 is unsigned bytes
                2 => i16::from_be_bytes([c[0], c[1]]) as i64,
                4 => i32::from_be_bytes([c[0], c[1], c[2], c[3]]) as i64,
                8 => i64::from_be_bytes(c.try_into().unwrap()),
                _ => return Err(Error::CompressionError("invalid ZBITPIX width".into())),
            };
            out.push(v);
        }
        Ok(out)
    }
}

/// Convert a linear tile index into per-axis tile coordinates (axis-1 fastest).
fn unravel(mut index: usize, tiles_per_axis: &[usize]) -> Vec<usize> {
    let mut coords = vec![0usize; tiles_per_axis.len()];
    for axis in 0..tiles_per_axis.len() {
        let n = tiles_per_axis[axis].max(1);
        coords[axis] = index % n;
        index /= n;
    }
    coords
}

/// Scatter a decoded tile's flat `values` (axis-1 fastest) into the full-image
/// buffer `full` (also axis-1 fastest), placing it at tile-grid `coords`.
///
/// - `image_dims` are the full `ZNAXISn` (axis-1 first).
/// - `ztile` are the nominal tile sizes (`ZTILEn`); used to locate the tile origin.
/// - `tile_dims` are the (possibly edge-truncated) dimensions of *this* tile; used
///   to size the copy. They equal `ztile` except for tiles on a high edge where
///   `ZNAXISn` is not a multiple of `ZTILEn`.
fn scatter_tile(
    full: &mut [i64],
    image_dims: &[usize],
    ztile: &[usize],
    tile_dims: &[usize],
    coords: &[usize],
    values: &[i64],
) {
    let ndim = image_dims.len();
    // Per-axis stride into the full image buffer (axis-1 fastest => stride 1).
    let mut img_stride = vec![1usize; ndim];
    for axis in 1..ndim {
        img_stride[axis] = img_stride[axis - 1] * image_dims[axis - 1];
    }
    // Origin pixel of this tile in the full image: coords stepped by the *nominal*
    // tile size along each axis.
    let mut origin = 0usize;
    for axis in 0..ndim {
        origin += coords[axis] * ztile[axis].max(1) * img_stride[axis];
    }

    // Iterate over every pixel of the tile via its multi-index (axis-1 fastest).
    let tile_npix: usize = tile_dims.iter().product();
    let mut tcoord = vec![0usize; ndim];
    for &val in values.iter().take(tile_npix) {
        let mut dst = origin;
        for axis in 0..ndim {
            dst += tcoord[axis] * img_stride[axis];
        }
        full[dst] = val;
        // Increment the tile multi-index (axis-1 fastest).
        for axis in 0..ndim {
            tcoord[axis] += 1;
            if tcoord[axis] < tile_dims[axis] {
                break;
            }
            tcoord[axis] = 0;
        }
    }
}

/// Count the contiguous `ZNAMEi`/`ZVALi` parameter pairs present in the header.
fn count_zname_params(header: &Header) -> usize {
    let mut n = 0;
    while header.find(&format!("ZNAME{}", n + 1)).is_some() {
        n += 1;
    }
    n
}

// ---------------------------------------------------------------------------
// RICE_1 decompression
// ---------------------------------------------------------------------------
//
// Port of cfitsio `fits_rdecomp` / `_short` / `_byte` (R. White, STScI). The first
// pixel is stored verbatim; subsequent pixels are stored as Rice-coded, zigzag-mapped
// differences from the previous pixel, in blocks of `blocksize` pixels. Each block is
// prefixed by `fsbits` bits giving `fs+1`:
//   * fs < 0      -> all differences in the block are zero
//   * fs == fsmax -> each difference is stored verbatim in `bbits` bits
//   * else        -> Rice code: unary leading zeros give the high bits, then `fs`
//                    low bits; combine as `(nzero << fs) | low`.

/// A simple big-endian-ish MSB-first bit reader over a byte slice.
struct BitReader<'a> {
    data: &'a [u8],
    byte_pos: usize,
    /// Number of valid bits remaining in the current `buffer` (0..=8).
    bits_in_buf: u32,
    buffer: u32,
}

impl<'a> BitReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        BitReader {
            data,
            byte_pos: 0,
            bits_in_buf: 0,
            buffer: 0,
        }
    }

    /// Read `n` bits (0..=32) MSB-first as an unsigned value.
    fn read_bits(&mut self, n: u32) -> Result<u32> {
        let mut result: u32 = 0;
        let mut need = n;
        while need > 0 {
            if self.bits_in_buf == 0 {
                let b = *self
                    .data
                    .get(self.byte_pos)
                    .ok_or_else(|| Error::CompressionError("RICE: unexpected end of stream".into()))?;
                self.byte_pos += 1;
                self.buffer = b as u32;
                self.bits_in_buf = 8;
            }
            let take = need.min(self.bits_in_buf);
            let shift = self.bits_in_buf - take;
            let mask = if take == 32 { u32::MAX } else { (1u32 << take) - 1 };
            let bits = (self.buffer >> shift) & mask;
            result = (result << take) | bits;
            self.bits_in_buf -= take;
            need -= take;
        }
        Ok(result)
    }

    /// Count and consume leading zero bits, then consume the terminating one-bit.
    /// Returns the number of zero bits seen.
    fn count_leading_zeros(&mut self) -> Result<u32> {
        let mut count = 0u32;
        loop {
            if self.bits_in_buf == 0 {
                let b = *self
                    .data
                    .get(self.byte_pos)
                    .ok_or_else(|| Error::CompressionError("RICE: unexpected end of stream".into()))?;
                self.byte_pos += 1;
                self.buffer = b as u32;
                self.bits_in_buf = 8;
            }
            // Inspect the top valid bit.
            let top = (self.buffer >> (self.bits_in_buf - 1)) & 1;
            self.bits_in_buf -= 1;
            if top == 1 {
                return Ok(count);
            }
            count += 1;
        }
    }
}

/// Map a zigzag-encoded unsigned difference back to a signed difference.
#[inline]
fn unzigzag(u: u64) -> i64 {
    if u & 1 == 0 {
        (u >> 1) as i64
    } else {
        !((u >> 1) as i64)
    }
}

/// Core RICE_1 decode shared by the 8/16/32-bit variants.
///
/// `fsbits`/`fsmax`/`bbits` are the per-variant constants. Returns reconstructed
/// values as `i64` (caller narrows to the target type).
fn rice_decompress_core(
    src: &[u8],
    nvals: usize,
    blocksize: usize,
    fsbits: u32,
    fsmax: u32,
    bbits: u32,
) -> Result<Vec<i64>> {
    if nvals == 0 {
        return Ok(Vec::new());
    }
    let blocksize = blocksize.max(1);
    let mut out = Vec::with_capacity(nvals);
    let mut reader = BitReader::new(src);

    // First value: bbits bits, verbatim (unsigned), interpreted as the raw integer.
    //
    // Note: the verbatim first value is the *seed* (`lastpix`), not itself a stored
    // pixel. cfitsio computes `array[0] = lastpix + diff[0]` and decodes `nvals`
    // differences in total (see `fits_rdecomp` in cfitsio `ricecomp.c`). Storing the
    // seed as the first output pixel — and then decoding only `nvals - 1` diffs —
    // shifts every value one position later, which is exactly wrong.
    let first = reader.read_bits(bbits)? as u64;
    let mut lastpix: i64 = sign_extend(first, bbits);

    while out.len() < nvals {
        let fs_plus_1 = reader.read_bits(fsbits)?;
        let fs = fs_plus_1 as i64 - 1;
        let remaining = nvals - out.len();
        let block_n = remaining.min(blocksize);

        if fs < 0 {
            // All differences zero.
            for _ in 0..block_n {
                out.push(lastpix);
            }
        } else if fs as u32 == fsmax {
            // Verbatim differences in bbits bits each.
            for _ in 0..block_n {
                let raw = reader.read_bits(bbits)? as u64;
                let diff = unzigzag(raw);
                lastpix = lastpix.wrapping_add(diff);
                out.push(lastpix);
            }
        } else {
            let fs = fs as u32;
            for _ in 0..block_n {
                let high = reader.count_leading_zeros()? as u64;
                let low = if fs > 0 { reader.read_bits(fs)? as u64 } else { 0 };
                let mapped = (high << fs) | low;
                let diff = unzigzag(mapped);
                lastpix = lastpix.wrapping_add(diff);
                out.push(lastpix);
            }
        }
    }

    Ok(out)
}

/// Sign-extend the low `bits` bits of `v` into an `i64`.
#[inline]
fn sign_extend(v: u64, bits: u32) -> i64 {
    if bits == 0 || bits >= 64 {
        return v as i64;
    }
    let shift = 64 - bits;
    ((v << shift) as i64) >> shift
}

/// Decompress a RICE_1 stream of 32-bit integers (`BYTEPIX = 4`).
pub fn rice_decompress_i32(src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<i32>> {
    let v = rice_decompress_core(src, nvals, blocksize, 5, 25, 32)?;
    Ok(v.into_iter().map(|x| x as i32).collect())
}

/// Decompress a RICE_1 stream of 16-bit integers (`BYTEPIX = 2`).
pub fn rice_decompress_i16(src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<i16>> {
    let v = rice_decompress_core(src, nvals, blocksize, 4, 14, 16)?;
    Ok(v.into_iter().map(|x| x as i16).collect())
}

/// Decompress a RICE_1 stream of 8-bit integers (`BYTEPIX = 1`).
pub fn rice_decompress_i8(src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<u8>> {
    let v = rice_decompress_core(src, nvals, blocksize, 3, 6, 8)?;
    Ok(v.into_iter().map(|x| x as u8).collect())
}

// ---------------------------------------------------------------------------
// Other algorithms (stubs)
// ---------------------------------------------------------------------------

/// Inflate a `GZIP_1` tile.
///
/// In the FITS Tiled Image Compression convention, `GZIP_1` tiles are full gzip
/// members (RFC 1952: gzip header + DEFLATE body + CRC32/ISIZE trailer), exactly as
/// produced by zlib's `gzip` routines and cfitsio. This routine strips the gzip
/// wrapper and runs a raw DEFLATE inflate over the body.
///
/// Available only with the `gzip` feature enabled; otherwise this returns
/// [`Error::UnsupportedCompression`].
#[cfg(feature = "gzip")]
pub fn gzip_inflate(src: &[u8]) -> Result<Vec<u8>> {
    let body = strip_gzip_wrapper(src)?;
    miniz_oxide::inflate::decompress_to_vec(body)
        .map_err(|e| Error::CompressionError(format!("DEFLATE inflate failed: {e:?}")))
}

/// Stub when the `gzip` feature is disabled: GZIP_1/GZIP_2 tiles cannot be decoded.
#[cfg(not(feature = "gzip"))]
pub fn gzip_inflate(_src: &[u8]) -> Result<Vec<u8>> {
    Err(Error::UnsupportedCompression(
        "GZIP_1/GZIP_2 require the `gzip` feature (miniz_oxide); rebuild with --features gzip".into(),
    ))
}

/// Parse and skip an RFC 1952 gzip member header, returning the DEFLATE body
/// (excluding the 8-byte CRC32+ISIZE trailer). Used by [`gzip_inflate`].
#[cfg(feature = "gzip")]
fn strip_gzip_wrapper(src: &[u8]) -> Result<&[u8]> {
    const FHCRC: u8 = 1 << 1;
    const FEXTRA: u8 = 1 << 2;
    const FNAME: u8 = 1 << 3;
    const FCOMMENT: u8 = 1 << 4;

    if src.len() < 18 || src[0] != 0x1f || src[1] != 0x8b {
        return Err(Error::CompressionError(
            "GZIP_1 tile is not a valid gzip member (bad magic)".into(),
        ));
    }
    if src[2] != 8 {
        return Err(Error::CompressionError(format!(
            "GZIP_1 unsupported compression method {}",
            src[2]
        )));
    }
    let flags = src[3];
    let mut pos = 10usize; // fixed header: magic(2)+CM(1)+FLG(1)+MTIME(4)+XFL(1)+OS(1)

    let end = || Error::CompressionError("GZIP_1 truncated gzip header".to_string());

    if flags & FEXTRA != 0 {
        if pos + 2 > src.len() {
            return Err(end());
        }
        let xlen = u16::from_le_bytes([src[pos], src[pos + 1]]) as usize;
        pos += 2 + xlen;
    }
    if flags & FNAME != 0 {
        while pos < src.len() && src[pos] != 0 {
            pos += 1;
        }
        pos += 1; // skip NUL
    }
    if flags & FCOMMENT != 0 {
        while pos < src.len() && src[pos] != 0 {
            pos += 1;
        }
        pos += 1;
    }
    if flags & FHCRC != 0 {
        pos += 2;
    }
    if pos + 8 > src.len() {
        return Err(end());
    }
    // Body is everything between the header and the 8-byte trailer.
    Ok(&src[pos..src.len() - 8])
}

/// `GZIP_2` stores the tile's integers with bytes shuffled into planes (all the
/// most-significant bytes first, then the next, etc.) before gzip. After inflating,
/// undo the shuffle to recover the original big-endian integer byte stream.
///
/// `bytepix` is the bytes-per-pixel of the original integers (1/2/4/8). For
/// `bytepix <= 1` the shuffle is a no-op.
#[cfg(feature = "gzip")]
fn gzip2_unshuffle(shuffled: &[u8], bytepix: usize) -> Vec<u8> {
    if bytepix <= 1 {
        return shuffled.to_vec();
    }
    let n = shuffled.len() / bytepix;
    let mut out = vec![0u8; n * bytepix];
    for (i, chunk) in out.chunks_exact_mut(bytepix).enumerate().take(n) {
        for (b, slot) in chunk.iter_mut().enumerate() {
            *slot = shuffled[b * n + i];
        }
    }
    out
}

/// `gzip2_unshuffle` is only reachable when the `gzip` feature is on (the GZIP_2
/// branch errors out in `gzip_inflate` otherwise), so it is feature-gated to avoid a
/// dead-code warning in the default build.
#[cfg(not(feature = "gzip"))]
#[allow(dead_code)]
fn gzip2_unshuffle(shuffled: &[u8], _bytepix: usize) -> Vec<u8> {
    shuffled.to_vec()
}

/// Decompress a PLIO_1 (IRAF pixel-list RLE) tile into i32 mask values. TODO(phase 4).
pub fn plio_decompress(_src: &[u8], _nvals: usize) -> Result<Vec<i32>> {
    todo!("phase 4: PLIO_1 pixel-list decode")
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal MSB-first bit writer mirroring [`BitReader`], for building test inputs.
    struct BitWriter {
        bytes: Vec<u8>,
        cur: u8,
        nbits: u32,
    }
    impl BitWriter {
        fn new() -> Self {
            BitWriter {
                bytes: Vec::new(),
                cur: 0,
                nbits: 0,
            }
        }
        fn put_bits(&mut self, val: u32, n: u32) {
            for i in (0..n).rev() {
                let bit = ((val >> i) & 1) as u8;
                self.cur = (self.cur << 1) | bit;
                self.nbits += 1;
                if self.nbits == 8 {
                    self.bytes.push(self.cur);
                    self.cur = 0;
                    self.nbits = 0;
                }
            }
        }
        fn finish(mut self) -> Vec<u8> {
            if self.nbits > 0 {
                self.cur <<= 8 - self.nbits;
                self.bytes.push(self.cur);
            }
            self.bytes
        }
    }

    fn zigzag(v: i64) -> u64 {
        ((v << 1) ^ (v >> 63)) as u64
    }

    #[test]
    fn unzigzag_round_trip() {
        for v in [-5i64, -1, 0, 1, 2, 100, -100, i32::MIN as i64, i32::MAX as i64] {
            assert_eq!(unzigzag(zigzag(v)), v);
        }
    }

    #[test]
    fn bit_reader_basic() {
        // 0b1011_0010, 0b1100_0000
        let data = [0b1011_0010u8, 0b1100_0000u8];
        let mut r = BitReader::new(&data);
        assert_eq!(r.read_bits(4).unwrap(), 0b1011);
        assert_eq!(r.read_bits(4).unwrap(), 0b0010);
        assert_eq!(r.read_bits(2).unwrap(), 0b11);
    }

    /// Hand-encode a small i32 RICE_1 stream (one block, fs known) and decode it.
    ///
    /// Mirrors cfitsio `fits_rcomp`: the first value is written verbatim as the *seed*
    /// (`lastpix = vals[0]`), and `nvals` differences are written — the first of which
    /// is `vals[0] - lastpix == 0`. The decoder then reconstructs `array[0] = seed + 0`.
    #[test]
    fn rice_i32_single_block_roundtrip() {
        // Original values; encode diffs with a chosen fs so we control the layout.
        let vals: Vec<i32> = vec![1000, 1003, 1001, 1005, 1002];
        let fs: u32 = 2; // low bits per code
        let encoded = rice_encode_i32(&vals, fs);

        let decoded = rice_decompress_i32(&encoded, vals.len(), 32).unwrap();
        assert_eq!(decoded, vals);
    }

    /// Zero-difference (constant) block: fs encoded as 0 (fs_plus_1 = 0 -> fs = -1).
    #[test]
    fn rice_i32_zero_block() {
        let vals: Vec<i32> = vec![42, 42, 42, 42];
        let mut w = BitWriter::new();
        w.put_bits(vals[0] as u32, 32);
        w.put_bits(0, 5); // fs_plus_1 = 0 => fs = -1 => all-zero block
        let encoded = w.finish();
        let decoded = rice_decompress_i32(&encoded, vals.len(), 32).unwrap();
        assert_eq!(decoded, vals);
    }

    #[test]
    fn tile_geometry_default_row_by_row() {
        let mut h = Header::new();
        h.set("ZNAXIS", HeaderValue::Integer(2), None);
        h.set("ZNAXIS1", HeaderValue::Integer(100), None);
        h.set("ZNAXIS2", HeaderValue::Integer(50), None);
        let g = TileGeometry::from_header(&h).unwrap();
        assert_eq!(g.znaxis, vec![100, 50]);
        assert_eq!(g.ztile, vec![100, 1]); // row-by-row default
        assert_eq!(g.num_tiles(), 50);
    }

    use crate::keyword::HeaderValue;

    // --- tile reassembly geometry --------------------------------------------

    #[test]
    fn unravel_axis1_fastest() {
        // 3x2 grid of tiles: index increments axis-1 (x) fastest.
        let tpa = [3usize, 2];
        assert_eq!(unravel(0, &tpa), vec![0, 0]);
        assert_eq!(unravel(1, &tpa), vec![1, 0]);
        assert_eq!(unravel(2, &tpa), vec![2, 0]);
        assert_eq!(unravel(3, &tpa), vec![0, 1]);
        assert_eq!(unravel(5, &tpa), vec![2, 1]);
    }

    #[test]
    fn scatter_full_tiles_2d() {
        // 4x4 image, 2x2 tiles => 2x2 grid. Each tile filled with its index*10+local.
        let image_dims = [4usize, 4];
        let ztile = [2usize, 2];
        let tiles_per_axis = [2usize, 2];
        let mut full = vec![-1i64; 16];
        for tile in 0..4 {
            let coords = unravel(tile, &tiles_per_axis);
            // tile-local values 0..4 (axis-1 fastest within the tile)
            let vals: Vec<i64> = (0..4).map(|l| (tile as i64) * 100 + l as i64).collect();
            scatter_tile(&mut full, &image_dims, &ztile, &[2, 2], &coords, &vals);
        }
        // Expected layout (row-major, axis-1 fastest):
        // row0: tile0[0] tile0[1] tile1[0] tile1[1]
        // row1: tile0[2] tile0[3] tile1[2] tile1[3]
        // row2: tile2[0] tile2[1] tile3[0] tile3[1]
        // row3: tile2[2] tile2[3] tile3[2] tile3[3]
        let expected = vec![
            0, 1, 100, 101, // row0
            2, 3, 102, 103, // row1
            200, 201, 300, 301, // row2
            202, 203, 302, 303, // row3
        ];
        assert_eq!(full, expected);
    }

    #[test]
    fn scatter_edge_truncated_tiles() {
        // 3x3 image with 2x2 tiles => 2x2 grid, but edge tiles are 1 wide/tall.
        let image_dims = [3usize, 3];
        let ztile = [2usize, 2];
        let tiles_per_axis = [2usize, 2];
        // tile dims: (0,0)=2x2, (1,0)=1x2, (0,1)=2x1, (1,1)=1x1
        let tile_dims = [vec![2, 2], vec![1, 2], vec![2, 1], vec![1, 1]];
        let mut full = vec![-1i64; 9];
        for tile in 0..4 {
            let coords = unravel(tile, &tiles_per_axis);
            let td = &tile_dims[tile];
            let n: usize = td.iter().product();
            let vals: Vec<i64> = (0..n).map(|l| (tile as i64) * 100 + l as i64).collect();
            scatter_tile(&mut full, &image_dims, &ztile, td, &coords, &vals);
        }
        // No -1 should remain; every pixel covered exactly once.
        assert!(!full.contains(&-1));
        // Spot check origins: tile0 at (0,0)->index0; tile1 (x-tile 1) origin x=2 row0 => index2.
        assert_eq!(full[0], 0); // tile0 local0
        assert_eq!(full[2], 100); // tile1 local0
        assert_eq!(full[6], 200); // tile2 local0 (row2, x0 => index 6)
        assert_eq!(full[8], 300); // tile3 local0 (row2, x2 => index 8)
    }

    // --- end-to-end RICE_1 decompress ----------------------------------------

    /// Build a one-block (<=32 values) i32 RICE_1 stream for `vals` with low-bits `fs`.
    ///
    /// Mirrors cfitsio `fits_rcomp`: writes `vals[0]` verbatim as the seed, sets
    /// `lastpix = vals[0]`, then emits exactly `vals.len()` differences — the first
    /// being `vals[0] - lastpix == 0`.
    fn rice_encode_i32(vals: &[i32], fs: u32) -> Vec<u8> {
        let mut w = BitWriter::new();
        w.put_bits(vals[0] as u32, 32);
        w.put_bits(fs + 1, 5);
        let mut last = vals[0];
        for &v in vals {
            let diff = (v - last) as i64;
            last = v;
            let mapped = zigzag(diff);
            let high = (mapped >> fs) as u32;
            let low = (mapped & ((1 << fs) - 1)) as u32;
            for _ in 0..high {
                w.put_bits(0, 1);
            }
            w.put_bits(1, 1);
            w.put_bits(low, fs);
        }
        w.finish()
    }

    #[test]
    fn end_to_end_rice1_i32_two_tiles() {
        use crate::bintable::{BinColumnType, BinTableBuilder};

        // 4x1 image, tile = 2x1 => two row tiles of 2 pixels each.
        let tile0 = [1000i32, 1003];
        let tile1 = [50i32, 47];
        let enc0 = rice_encode_i32(&tile0, 2);
        let enc1 = rice_encode_i32(&tile1, 2);

        let table = BinTableBuilder::new()
            .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
            .push_row(|r| r.write_var_p(enc0.len() as i32, |heap| heap.extend_from_slice(&enc0)))
            .push_row(|r| r.write_var_p(enc1.len() as i32, |heap| heap.extend_from_slice(&enc1)))
            .build();

        let mut h = Header::new();
        h.set("ZIMAGE", HeaderValue::Logical(true), None);
        h.set("ZCMPTYPE", HeaderValue::String("RICE_1".into()), None);
        h.set("ZBITPIX", HeaderValue::Integer(32), None);
        h.set("ZNAXIS", HeaderValue::Integer(1), None);
        h.set("ZNAXIS1", HeaderValue::Integer(4), None);
        h.set("ZTILE1", HeaderValue::Integer(2), None);
        // BYTEPIX = 4 for i32 RICE.
        h.set("ZNAME1", HeaderValue::String("BYTEPIX".into()), None);
        h.set("ZVAL1", HeaderValue::Integer(4), None);

        let cimg = CompressedImage::from_bintable(&h, &table).unwrap();
        assert_eq!(cimg.compression(), CompressionType::Rice1);
        assert_eq!(cimg.geometry().num_tiles(), 2);

        let img = cimg.decompress().unwrap();
        assert_eq!(img.axes, vec![4]);
        match img.pixels {
            crate::image_data::PixelData::I32(v) => {
                assert_eq!(v, vec![1000, 1003, 50, 47]);
            }
            other => panic!("expected I32, got {other:?}"),
        }
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn end_to_end_gzip1_i16_one_tile() {
        use crate::bintable::{BinColumnType, BinTableBuilder};

        // 4x1 i16 image, one row tile. GZIP_1 stores raw big-endian image integers.
        let vals: [i16; 4] = [100, -200, 300, -400];
        let mut raw = Vec::new();
        for v in vals {
            raw.extend_from_slice(&v.to_be_bytes());
        }
        // Wrap raw bytes as a gzip member (header + raw DEFLATE + trailer).
        let body = miniz_oxide::deflate::compress_to_vec(&raw, 6);
        let mut gz = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
        gz.extend_from_slice(&body);
        let crc = crc32(&raw);
        gz.extend_from_slice(&crc.to_le_bytes());
        gz.extend_from_slice(&(raw.len() as u32).to_le_bytes());

        let table = BinTableBuilder::new()
            .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
            .push_row(|r| r.write_var_p(gz.len() as i32, |heap| heap.extend_from_slice(&gz)))
            .build();

        let mut h = Header::new();
        h.set("ZIMAGE", HeaderValue::Logical(true), None);
        h.set("ZCMPTYPE", HeaderValue::String("GZIP_1".into()), None);
        h.set("ZBITPIX", HeaderValue::Integer(16), None);
        h.set("ZNAXIS", HeaderValue::Integer(1), None);
        h.set("ZNAXIS1", HeaderValue::Integer(4), None);

        let cimg = CompressedImage::from_bintable(&h, &table).unwrap();
        let img = cimg.decompress().unwrap();
        match img.pixels {
            crate::image_data::PixelData::I16(v) => assert_eq!(v, vec![100, -200, 300, -400]),
            other => panic!("expected I16, got {other:?}"),
        }
    }

    /// Minimal CRC32 (gzip/PNG polynomial) for building gzip test fixtures.
    #[cfg(feature = "gzip")]
    fn crc32(data: &[u8]) -> u32 {
        let mut crc: u32 = 0xffff_ffff;
        for &byte in data {
            crc ^= byte as u32;
            for _ in 0..8 {
                let mask = (crc & 1).wrapping_neg();
                crc = (crc >> 1) ^ (0xedb8_8320 & mask);
            }
        }
        !crc
    }

    #[cfg(not(feature = "gzip"))]
    #[test]
    fn gzip_disabled_errors() {
        assert!(matches!(
            gzip_inflate(&[0x1f, 0x8b]),
            Err(Error::UnsupportedCompression(_))
        ));
    }

    #[test]
    fn float_zbitpix_is_unsupported() {
        use crate::bintable::{BinColumnType, BinTableBuilder};
        let table = BinTableBuilder::new()
            .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
            .push_row(|r| r.write_var_p(0, |_| {}))
            .build();
        let mut h = Header::new();
        h.set("ZIMAGE", HeaderValue::Logical(true), None);
        h.set("ZCMPTYPE", HeaderValue::String("RICE_1".into()), None);
        h.set("ZBITPIX", HeaderValue::Integer(-32), None);
        h.set("ZNAXIS", HeaderValue::Integer(1), None);
        h.set("ZNAXIS1", HeaderValue::Integer(2), None);
        let cimg = CompressedImage::from_bintable(&h, &table).unwrap();
        assert!(matches!(
            cimg.decompress(),
            Err(Error::UnsupportedCompression(_))
        ));
    }
}
