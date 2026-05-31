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
//! Read (decompression) support is implemented for `RICE_1`, `GZIP_1`/`GZIP_2`
//! (behind the `gzip` feature), `PLIO_1`, `HCOMPRESS_1` (`SMOOTH=0`), and
//! `NOCOMPRESS`, for integer originals and for quantized/lossless float originals
//! (with `SUBTRACTIVE_DITHER_1`/`_2` reversal).
//!
//! Write (compression) support is implemented for `RICE_1` (integer images and
//! quantized/dithered float images) and `GZIP_1`/`GZIP_2` (integer images; lossless
//! raw-float storage via `GZIP_1`), via [`compress_image`] / [`ImageData::compress`].
//! The encoders are the exact inverse of the decoders above and emit standard,
//! cfitsio-readable compressed FITS (verified bit-exact against `funpack`). `PLIO_1`
//! and `HCOMPRESS_1` encoding remain unimplemented (see `COMPRESSION_PLAN.md`).

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
                // `NONE` is written by cfitsio for losslessly-stored float tiles
                // (GZIP of the raw floats, no quantization) and is treated like the
                // absent case here.
                "NO_DITHER" | "NONE" => Ok(Quantize::None),
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
    /// HCOMPRESS `SMOOTH` parameter (default 0 = no smoothing).
    pub smooth: i32,
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
        let mut smooth = 0i32;
        let nzparams = count_zname_params(header);
        for i in 1..=nzparams {
            if let Some(name) = header.get_string(&format!("ZNAME{i}")) {
                let val = header.get_int(&format!("ZVAL{i}")).unwrap_or(0);
                match name.trim().to_ascii_uppercase().as_str() {
                    "BLOCKSIZE" => blocksize = val as usize,
                    "BYTEPIX" => bytepix = val as usize,
                    "SMOOTH" => smooth = val as i32,
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
            smooth,
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
    /// Float `ZBITPIX` (-32/-64) images are decoded via [`Self::decompress_float`],
    /// which inverts the per-tile `ZSCALE`/`ZZERO` quantization and reverses cfitsio's
    /// subtractive dithering (`ZQUANTIZ` = `SUBTRACTIVE_DITHER_1`/`_2`); `ZBLANK`
    /// sentinels map to `NaN`.
    ///
    /// `PLIO_1` (IRAF pixel-list RLE) and `HCOMPRESS_1` (H-transform + quadtree)
    /// integer tiles are decoded by their dedicated per-tile decoders and scattered
    /// like any other integer codec.
    ///
    /// # Not supported
    ///
    /// - `GZIP_1`/`GZIP_2` without the `gzip` feature return
    ///   [`Error::UnsupportedCompression`].
    /// - `HCOMPRESS_1` with `SMOOTH != 0` (smoothing) is rejected.
    pub fn decompress(&self) -> Result<crate::image_data::ImageData> {
        use crate::image_data::{ImageData, PixelData};
        use crate::types::Bitpix;

        // Float originals (ZBITPIX = -32 / -64) take a separate path: each tile holds
        // either quantized integers (with per-tile ZSCALE/ZZERO) or, when quantization
        // was not applied, the raw floats. See `decompress_float`.
        if self.zbitpix < 0 {
            return self.decompress_float();
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
            // Number of pixels in this (possibly edge-truncated) tile.
            let coords = unravel(tile_index, &tiles_per_axis);
            let tile_dims = self.tile_dims_at(&coords);
            let tile_npix: usize = tile_dims.iter().product();

            // PLIO_1 stores its line list as a `1PI` (i16) VLA; every other codec
            // uses a byte (`1PB`) VLA. Fetch the cell in its native form and decode.
            let values = if self.ctype == CompressionType::Plio1 {
                let words = self.tile_words_i16(tile_index, cdata_col)?;
                plio_decompress(&words, tile_npix)?
                    .into_iter()
                    .map(|v| v as i64)
                    .collect()
            } else if self.ctype == CompressionType::Hcompress1 {
                let raw = self.tile_bytes(tile_index, cdata_col)?;
                self.decode_hcompress(&raw, &tile_dims)?
            } else {
                let raw = self.tile_bytes(tile_index, cdata_col)?;
                self.decode_tile(&raw, tile_npix)?
            };
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

    /// Decompress a float (`ZBITPIX = -32` / `-64`) tile-compressed image.
    ///
    /// Float images are normally stored as per-tile linearly-quantized 32-bit
    /// integers plus per-tile `ZSCALE`/`ZZERO` scale/offset columns; the integers are
    /// reconstructed exactly like the integer path, then *unquantized* back to floats.
    /// cfitsio reverses the subtractive-dithering it applied at compress time using a
    /// fixed pseudo-random sequence (see [`fits_rand_value`]).
    ///
    /// A tile that could not be quantized (e.g. one containing only NaNs, or when
    /// lossless `GZIP` of the raw floats was requested) is stored instead as the raw
    /// big-endian float bytes — either in a `GZIP_COMPRESSED_DATA`/`UNCOMPRESSED_DATA`
    /// fallback column for that single tile, or, for a wholly-lossless image, in
    /// `COMPRESSED_DATA` with no `ZSCALE`/`ZZERO` columns at all. Such tiles are passed
    /// through verbatim.
    fn decompress_float(&self) -> Result<crate::image_data::ImageData> {
        use crate::image_data::{ImageData, PixelData};

        // -32 -> 4-byte floats, -64 -> 8-byte floats.
        let is_f64 = self.zbitpix == -64;
        let npix: usize = self.geometry.znaxis.iter().product();

        let tiles_per_axis = self.geometry.tiles_per_axis();
        let num_tiles = self.geometry.num_tiles();
        if self.table.nrows < num_tiles {
            return Err(Error::CompressionError(format!(
                "compressed image expects {} tiles but BINTABLE has {} rows",
                num_tiles, self.table.nrows
            )));
        }

        let cdata_col = self.compressed_data_column()?;
        let gzip_fallback_col = self.find_column("GZIP_COMPRESSED_DATA");
        let uncompressed_col = self.find_column("UNCOMPRESSED_DATA");
        let zscale_col = self.find_column("ZSCALE");
        let zzero_col = self.find_column("ZZERO");
        let zblank_col = self.find_column("ZBLANK");

        // Reconstructed floats, axis-1 fastest, scattered tile-by-tile. We scatter via
        // the integer `scatter_tile` over a bit-reinterpreted buffer so the existing
        // (well-tested) geometry code is reused.
        let mut full_bits = vec![0i64; npix];

        for tile_index in 0..num_tiles {
            let coords = unravel(tile_index, &tiles_per_axis);
            let tile_dims = self.tile_dims_at(&coords);
            let tile_npix: usize = tile_dims.iter().product();

            // Per-tile scale/zero (D columns). A sentinel/absent value marks an
            // unquantized (raw-float) tile.
            let zscale = zscale_col.and_then(|c| self.tile_double(tile_index, c));
            let zzero = zzero_col.and_then(|c| self.tile_double(tile_index, c));
            // Per-tile blank sentinel (J column) overrides the ZBLANK keyword.
            let blank = zblank_col
                .and_then(|c| self.tile_int(tile_index, c))
                .or(self.blank);

            let cdata = self.tile_bytes(tile_index, cdata_col)?;

            let tile_floats: Vec<f64> = match (zscale, zzero) {
                // Quantized tile: COMPRESSED_DATA holds quantized integers.
                (Some(scale), Some(zero)) if !cdata.is_empty() && is_quantized(scale) => {
                    let q = if self.ctype == CompressionType::Hcompress1 {
                        self.decode_hcompress(&cdata, &tile_dims)?
                    } else {
                        self.decode_tile(&cdata, tile_npix)?
                    };
                    self.unquantize(&q, tile_index, scale, zero, blank, tile_npix)
                }
                // Unquantized tile (or whole image is lossless): raw floats live in
                // COMPRESSED_DATA, or in a per-tile fallback column when COMPRESSED_DATA
                // is empty.
                _ => {
                    let (bytes, gzipped) = if !cdata.is_empty() {
                        // For a lossless GZIP image the raw floats are gzip-compressed
                        // here; for NOCOMPRESS they are verbatim.
                        (cdata, self.ctype != CompressionType::NoCompress)
                    } else if let Some(b) =
                        gzip_fallback_col.and_then(|c| self.tile_bytes(tile_index, c).ok())
                    {
                        (b, true)
                    } else if let Some(b) =
                        uncompressed_col.and_then(|c| self.tile_bytes(tile_index, c).ok())
                    {
                        (b, false)
                    } else {
                        return Err(Error::CompressionError(format!(
                            "float tile {tile_index} has no quantization scale and no raw-float fallback data"
                        )));
                    };
                    let raw = if gzipped { gzip_inflate(&bytes)? } else { bytes };
                    raw_floats(&raw, tile_npix, is_f64)?
                }
            };

            if tile_floats.len() != tile_npix {
                return Err(Error::CompressionError(format!(
                    "float tile {tile_index} produced {} values, expected {tile_npix}",
                    tile_floats.len()
                )));
            }

            // Reinterpret each float's bits as an integer so we can reuse scatter_tile,
            // then convert back below. (f32 bits zero-extended into i64.)
            let bits: Vec<i64> = if is_f64 {
                tile_floats.iter().map(|&v| v.to_bits() as i64).collect()
            } else {
                tile_floats
                    .iter()
                    .map(|&v| (v as f32).to_bits() as i64)
                    .collect()
            };
            scatter_tile(
                &mut full_bits,
                &self.geometry.znaxis,
                &self.geometry.ztile,
                &tile_dims,
                &coords,
                &bits,
            );
        }

        let pixels = if is_f64 {
            PixelData::F64(full_bits.iter().map(|&b| f64::from_bits(b as u64)).collect())
        } else {
            PixelData::F32(
                full_bits
                    .iter()
                    .map(|&b| f32::from_bits(b as u32))
                    .collect(),
            )
        };

        Ok(ImageData::new(self.geometry.znaxis.clone(), pixels))
    }

    /// Inverse linear quantization with cfitsio's subtractive-dithering reversal.
    ///
    /// For each quantized integer `q[i]` of a tile:
    /// - `q == blank` (ZBLANK / per-tile null sentinel) ⇒ `NaN`.
    /// - `SUBTRACTIVE_DITHER_2` and `q == ZERO_VALUE (-2147483646)` ⇒ exactly `0.0`.
    /// - otherwise `value = (q - r + 0.5) * scale + zero`, where `r` is the next value
    ///   of the fixed pseudo-random sequence ([`fits_rand_value`]); for the
    ///   non-dithered methods `r` is taken as `0.0` (so `value = q*scale + zero`,
    ///   matching cfitsio's `fffi4r4`).
    ///
    /// The random-sequence indexing mirrors cfitsio `unquantize_i4r4`:
    /// `iseed = (tile_index + ZDITHER0 - 1) mod N_RANDOM`,
    /// `nextrand = (int)(fits_rand_value[iseed] * 500)`, advancing `nextrand` per pixel
    /// and, on reaching `N_RANDOM`, bumping `iseed` (wrapping) and re-deriving
    /// `nextrand`.
    fn unquantize(
        &self,
        q: &[i64],
        tile_index: usize,
        scale: f64,
        zero: f64,
        blank: Option<i64>,
        tile_npix: usize,
    ) -> Vec<f64> {
        let dithered = matches!(
            self.quantize,
            Quantize::SubtractiveDither1 | Quantize::SubtractiveDither2
        );
        let dither2 = self.quantize == Quantize::SubtractiveDither2;

        // ZDITHER0 defaults to 1 when absent (cfitsio fits_read_compressed_img).
        let zdither0 = self.zdither0.unwrap_or(1);
        // cfitsio passes row = tile_index_1based + zdither0 - 1, then iseed = (row-1) % N.
        let mut iseed = ((tile_index as i64 + zdither0 - 1).rem_euclid(N_RANDOM as i64)) as usize;
        let mut nextrand = (fits_rand_value(iseed) * 500.0) as usize;

        let mut out = Vec::with_capacity(tile_npix);
        for &qi in q.iter().take(tile_npix) {
            // cfitsio computes `x * scale + zero` as a single fused multiply-add
            // (the C compiler contracts the expression), which we must reproduce
            // exactly via `mul_add` to be bit-identical in catastrophic-cancellation
            // cases (large `zero` offsets).
            let value = if blank.is_some_and(|b| qi == b) {
                f64::NAN
            } else if dither2 && qi == ZERO_VALUE {
                0.0
            } else if dithered {
                ((qi as f64) - fits_rand_value(nextrand) as f64 + 0.5).mul_add(scale, zero)
            } else {
                (qi as f64).mul_add(scale, zero)
            };
            out.push(value);

            if dithered {
                nextrand += 1;
                if nextrand == N_RANDOM {
                    iseed += 1;
                    if iseed == N_RANDOM {
                        iseed = 0;
                    }
                    nextrand = (fits_rand_value(iseed) * 500.0) as usize;
                }
            }
        }
        out
    }

    /// Index of a column by (case-insensitive, trimmed) `TTYPE` name, if present.
    fn find_column(&self, name: &str) -> Option<usize> {
        self.table
            .columns
            .iter()
            .position(|c| c.name.trim().eq_ignore_ascii_case(name))
    }

    /// Read a single `D`/`E` scalar cell as `f64` (per-tile ZSCALE/ZZERO).
    fn tile_double(&self, row: usize, col: usize) -> Option<f64> {
        match self.table.get_cell(row, col).ok()? {
            crate::bintable::BinCellValue::F64(v) => v.first().copied(),
            crate::bintable::BinCellValue::F32(v) => v.first().map(|&x| x as f64),
            _ => None,
        }
    }

    /// Read a single integer scalar cell as `i64` (per-tile ZBLANK).
    fn tile_int(&self, row: usize, col: usize) -> Option<i64> {
        match self.table.get_cell(row, col).ok()? {
            crate::bintable::BinCellValue::I32(v) => v.first().map(|&x| x as i64),
            crate::bintable::BinCellValue::I64(v) => v.first().copied(),
            crate::bintable::BinCellValue::I16(v) => v.first().map(|&x| x as i64),
            _ => None,
        }
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

    /// Read one tile's `1PI`/`1QI` cell as native `i16` words (PLIO_1 line list).
    fn tile_words_i16(&self, row: usize, col: usize) -> Result<Vec<i16>> {
        match self.table.get_cell(row, col)? {
            crate::bintable::BinCellValue::I16(v) => Ok(v),
            // Some writers may declare the PLIO list as a byte VLA; reinterpret as
            // big-endian i16 pairs in that case.
            crate::bintable::BinCellValue::Bytes(b) => Ok(b
                .chunks_exact(2)
                .map(|c| i16::from_be_bytes([c[0], c[1]]))
                .collect()),
            other => Err(Error::CompressionError(format!(
                "PLIO_1 COMPRESSED_DATA cell is not an i16 VLA: {other:?}"
            ))),
        }
    }

    /// Decode one HCOMPRESS_1 tile into a flat array of `i64` integer samples
    /// (axis-1 fastest, matching the tile's pixel layout).
    ///
    /// The tile dimensions (`nx`, `ny`) and digitization `scale` are read from the
    /// encoded stream itself; `tile_dims` is `[width, height]` (axis-1 first) and is
    /// used only to validate the decoded size. HCOMPRESS lays out its output array
    /// with `ny` (== tile width, the FITS fast axis) varying fastest, which already
    /// matches the axis-1-fastest order the scatter step expects.
    fn decode_hcompress(&self, raw: &[u8], tile_dims: &[usize]) -> Result<Vec<i64>> {
        // SMOOTH parameter (ZNAMEi='SMOOTH'); default 0.
        let smooth = self.smooth;
        // cfitsio uses the 32-bit H-transform (`fits_hdecompress`) only for
        // ZBITPIX 8/16; every other original type (32-bit ints, and quantized
        // -32/-64 floats) uses the 64-bit transform to avoid intermediate overflow.
        let wide = !(self.zbitpix == 8 || self.zbitpix == 16);
        let (a, nx, ny) = hcompress_decompress(raw, wide, smooth)?;
        let tile_npix: usize = tile_dims.iter().product();
        if nx * ny != tile_npix {
            return Err(Error::CompressionError(format!(
                "HCOMPRESS tile size {nx}x{ny} = {} != expected {tile_npix}",
                nx * ny
            )));
        }
        Ok(a)
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
            // PLIO_1 and HCOMPRESS_1 are dispatched to their dedicated per-tile
            // decoders (`plio_decompress` / `decode_hcompress`) in `decompress` /
            // `decompress_float`, so they never reach this generic byte path.
            CompressionType::Plio1 => Err(Error::CompressionError(
                "internal: PLIO_1 must be decoded via plio_decompress".into(),
            )),
            CompressionType::Hcompress1 => Err(Error::CompressionError(
                "internal: HCOMPRESS_1 must be decoded via decode_hcompress".into(),
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

// ---------------------------------------------------------------------------
// Float quantization / subtractive-dithering reconstruction
// ---------------------------------------------------------------------------
//
// Ported from cfitsio (`imcompress.c` / `quantize.c`, R. White & W. Pence, STScI).
// Quantized float tiles store each pixel as an i32 produced by
// `round(value/scale - zero/scale + dither)`; we invert with the same fixed
// pseudo-random `fits_rand_value` table that cfitsio uses.

/// Length of cfitsio's fixed pseudo-random sequence (`N_RANDOM`). Do not change.
const N_RANDOM: usize = 10000;

/// Sentinel quantized value flagging an undefined (NaN) pixel (cfitsio `NULL_VALUE`).
/// Carried here for documentation; the actual null sentinel is taken from
/// `ZBLANK` (keyword or per-tile column), which equals this for cfitsio-written files.
#[allow(dead_code)]
const NULL_VALUE: i64 = -2_147_483_647;

/// Sentinel quantized value flagging an exact-zero pixel under
/// `SUBTRACTIVE_DITHER_2` (cfitsio `ZERO_VALUE`).
const ZERO_VALUE: i64 = -2_147_483_646;

/// True if a per-tile `ZSCALE` indicates a quantized tile. cfitsio treats a zero
/// scale (its `cn_zscale == 0` / absent-scale default) as "not quantized".
fn is_quantized(zscale: f64) -> bool {
    zscale != 0.0
}

/// Reinterpret big-endian bytes as `n` floats (f32 if `!is_f64`, else f64), widened
/// to `f64` for uniform downstream handling. Used for unquantized (raw-float) tiles.
fn raw_floats(bytes: &[u8], n: usize, is_f64: bool) -> Result<Vec<f64>> {
    let width = if is_f64 { 8 } else { 4 };
    if bytes.len() < n * width {
        return Err(Error::CompressionError(format!(
            "raw float tile has {} bytes, expected at least {}",
            bytes.len(),
            n * width
        )));
    }
    let mut out = Vec::with_capacity(n);
    for c in bytes.chunks_exact(width).take(n) {
        let v = if is_f64 {
            f64::from_be_bytes(c.try_into().unwrap())
        } else {
            f32::from_be_bytes([c[0], c[1], c[2], c[3]]) as f64
        };
        out.push(v);
    }
    Ok(out)
}

/// The `ii`-th element of cfitsio's fixed pseudo-random sequence.
///
/// Generated once (lazily) by the Park–Miller minimal-standard LCG
/// (`a = 16807`, `m = 2^31 - 1`, seed 1) exactly as cfitsio's `fits_init_randoms`:
/// `seed_{k+1} = (a*seed_k) mod m`, value = `seed / m`. The 10000-element table is
/// validated against cfitsio's published checkpoint (final seed `1043618065`).
fn fits_rand_value(ii: usize) -> f32 {
    use std::sync::OnceLock;
    static TABLE: OnceLock<Vec<f32>> = OnceLock::new();
    let table = TABLE.get_or_init(|| {
        let a = 16807.0f64;
        let m = 2_147_483_647.0f64;
        let mut seed = 1.0f64;
        let mut v = Vec::with_capacity(N_RANDOM);
        for _ in 0..N_RANDOM {
            let temp = a * seed;
            // cfitsio: seed = temp - m * (int)(temp / m)
            seed = temp - m * ((temp / m) as i64 as f64);
            v.push((seed / m) as f32);
        }
        v
    });
    table[ii % N_RANDOM]
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
// RICE_1 compression (encode) — inverse of `rice_decompress_*`
// ---------------------------------------------------------------------------
//
// Port of cfitsio `fits_rcomp` / `_short` / `_byte` (R. White, STScI, `ricecomp.c`).
// The encoder writes the first value verbatim (the seed `lastpix`) in `bbits` bits,
// then processes the array in blocks of `blocksize` pixels. For each block it computes
// the zig-zag-mapped differences from the running `lastpix`, picks the Rice parameter
// `fs` that minimises the block's coded length, and emits `fs+1` in `fsbits` bits
// followed by the per-pixel codes:
//   * fs < 0      -> block of all-zero differences (fs+1 == 0).
//   * fs == fsmax -> each mapped diff stored verbatim in `bbits` bits.
//   * else        -> Rice code: `(top = mapped >> fs)` zero bits, a one bit, then the
//                    `fs` low bits of `mapped`.
// The chosen `fs` exactly matches cfitsio's sum-based selection so the byte stream is
// identical to `fpack`'s and decodes via [`rice_decompress_core`].

/// MSB-first bit writer (the inverse of [`BitReader`]).
struct BitWriterEnc {
    bytes: Vec<u8>,
    /// Bits accumulated in the low `bits_in_buf` bits of `buffer`, flushed a byte at a
    /// time (MSB-first) whenever at least 8 are present. A u64 buffer keeps headroom for
    /// a 32-bit write on top of up to 7 leftover bits.
    buffer: u64,
    bits_in_buf: u32,
}

impl BitWriterEnc {
    fn new() -> Self {
        BitWriterEnc {
            bytes: Vec::new(),
            buffer: 0,
            bits_in_buf: 0,
        }
    }

    /// Write the low `n` bits of `val` (0..=32) MSB-first.
    fn write_bits(&mut self, val: u32, n: u32) {
        if n == 0 {
            return;
        }
        let val = (val as u64) & if n >= 32 { u32::MAX as u64 } else { (1u64 << n) - 1 };
        self.buffer = (self.buffer << n) | val;
        self.bits_in_buf += n;
        while self.bits_in_buf >= 8 {
            self.bits_in_buf -= 8;
            let byte = (self.buffer >> self.bits_in_buf) & 0xff;
            self.bytes.push(byte as u8);
        }
    }

    /// Write `count` zero bits then a single one bit (the unary part of a Rice code).
    fn write_unary(&mut self, count: u32) {
        let mut remaining = count;
        while remaining >= 24 {
            self.write_bits(0, 24);
            remaining -= 24;
        }
        // remaining zeros followed by a terminating 1 bit.
        self.write_bits(1, remaining + 1);
    }

    /// Flush any partial byte (zero-padded on the low side, as cfitsio does) and return
    /// the encoded bytes.
    fn finish(mut self) -> Vec<u8> {
        if self.bits_in_buf > 0 {
            let pad = 8 - self.bits_in_buf;
            let byte = (self.buffer << pad) & 0xff;
            self.bytes.push(byte as u8);
            self.bits_in_buf = 0;
        }
        debug_assert!(self.bits_in_buf == 0);
        self.bytes
    }
}

/// Zig-zag map a signed difference to an unsigned value (inverse of [`unzigzag`]).
#[inline]
fn zigzag_enc(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

/// Core RICE_1 encode shared by the 8/16/32-bit variants. `vals` are the raw integer
/// samples (sign-extended to i64); `fsbits`/`fsmax`/`bbits` are the per-variant
/// constants matching [`rice_decompress_core`].
fn rice_compress_core(
    vals: &[i64],
    blocksize: usize,
    fsbits: u32,
    fsmax: u32,
    bbits: u32,
) -> Vec<u8> {
    let mut w = BitWriterEnc::new();
    if vals.is_empty() {
        return w.finish();
    }
    let blocksize = blocksize.max(1);

    // Seed: first value verbatim in bbits bits (masked to the field width).
    let mask = if bbits >= 64 {
        u64::MAX
    } else {
        (1u64 << bbits) - 1
    };
    let mut lastpix = vals[0];
    w.write_bits((vals[0] as u64 & mask) as u32, bbits);

    let mut i = 0usize;
    while i < vals.len() {
        let block_n = (vals.len() - i).min(blocksize);
        // Zig-zag-mapped differences for this block.
        let mut mapped = Vec::with_capacity(block_n);
        let mut sum: u64 = 0;
        for &v in &vals[i..i + block_n] {
            // Compute the difference in the pixel's native bit width: reduce mod 2^bbits
            // and sign-extend so it lands in [-2^(bbits-1), 2^(bbits-1)). This guarantees
            // the zig-zag mapping fits in `bbits` bits (the verbatim-store width), and the
            // decoder's `wrapping_add` + final narrowing reconstructs `v` exactly even
            // when the raw difference wraps the type (e.g. MIN→MAX in i16).
            let raw_diff = v.wrapping_sub(lastpix);
            lastpix = v;
            let diff = sign_extend((raw_diff as u64) & mask, bbits);
            let m = zigzag_enc(diff);
            sum = sum.wrapping_add(m);
            mapped.push(m);
        }

        // Choose the Rice parameter `fs` that minimises this block's coded length.
        let fs = select_fs(sum, block_n as u64, fsmax);

        if fs < 0 {
            // All differences are zero (sum == 0): emit fs+1 == 0.
            w.write_bits(0, fsbits);
        } else if fs as u32 >= fsmax {
            // No compression: store each mapped diff verbatim in bbits bits.
            w.write_bits(fsmax + 1, fsbits);
            for &m in &mapped {
                w.write_bits((m & mask) as u32, bbits);
            }
        } else {
            let fs = fs as u32;
            w.write_bits(fs + 1, fsbits);
            for &m in &mapped {
                let top = (m >> fs) as u32;
                w.write_unary(top);
                if fs > 0 {
                    w.write_bits((m & ((1u64 << fs) - 1)) as u32, fs);
                }
            }
        }

        i += block_n;
    }

    w.finish()
}

/// Select the Rice parameter `fs` for a block by minimising its coded length.
///
/// For `block_n` pixels whose zig-zag-mapped differences sum to `sum`, coding with a
/// given `fs` costs (approximately, and exactly for the dominant terms cfitsio uses)
/// `block_n * (fs + 1)` fixed bits (the `fs` low bits plus one stop bit per pixel) plus
/// `sum >> fs` unary tail bits. This is the same convex cost model cfitsio's
/// `fits_rcomp` minimises; because the per-block `fs+1` value is written into the
/// stream, the result decodes correctly with *any* valid `fs`, and minimising this cost
/// reproduces cfitsio's choice on real data. Returns `-1` for an all-zero block
/// (`sum == 0`) and never exceeds `fsmax` (the caller switches to verbatim storage at
/// `fsmax`).
fn select_fs(sum: u64, block_n: u64, fsmax: u32) -> i32 {
    if sum == 0 {
        return -1;
    }
    // Walk fs upward while increasing it lowers the estimated bit count; the cost is
    // unimodal in fs, so stop at the first non-improving step.
    let cost = |fs: u32| block_n.wrapping_mul((fs + 1) as u64).wrapping_add(sum >> fs);
    let mut best_fs = 0u32;
    let mut best_cost = cost(0);
    let mut fs = 1u32;
    while fs <= fsmax {
        let c = cost(fs);
        if c < best_cost {
            best_cost = c;
            best_fs = fs;
            fs += 1;
        } else {
            break;
        }
    }
    best_fs as i32
}

/// Compress 32-bit integer samples with RICE_1 (`BYTEPIX = 4`).
pub fn rice_compress_i32(vals: &[i32], blocksize: usize) -> Vec<u8> {
    let v: Vec<i64> = vals.iter().map(|&x| x as i64).collect();
    rice_compress_core(&v, blocksize, 5, 25, 32)
}

/// Compress 16-bit integer samples with RICE_1 (`BYTEPIX = 2`).
pub fn rice_compress_i16(vals: &[i16], blocksize: usize) -> Vec<u8> {
    let v: Vec<i64> = vals.iter().map(|&x| x as i64).collect();
    rice_compress_core(&v, blocksize, 4, 14, 16)
}

/// Compress 8-bit (unsigned) integer samples with RICE_1 (`BYTEPIX = 1`).
pub fn rice_compress_i8(vals: &[u8], blocksize: usize) -> Vec<u8> {
    let v: Vec<i64> = vals.iter().map(|&x| x as i64).collect();
    rice_compress_core(&v, blocksize, 3, 6, 8)
}

// ---------------------------------------------------------------------------
// GZIP_1 / GZIP_2 compression (encode)
// ---------------------------------------------------------------------------

/// Wrap raw bytes as a single RFC 1952 gzip member (the inverse of
/// [`strip_gzip_wrapper`] + [`gzip_inflate`]): a fixed 10-byte header, a raw DEFLATE
/// body, and the CRC32 + ISIZE trailer. Matches what zlib's `gzip` (and cfitsio's
/// `GZIP_1`) produce, so the existing decoder and `funpack` both read it back.
#[cfg(feature = "gzip")]
fn gzip_deflate(raw: &[u8]) -> Vec<u8> {
    let body = miniz_oxide::deflate::compress_to_vec(raw, 6);
    let mut out = Vec::with_capacity(body.len() + 18);
    // magic, CM=deflate, no flags, mtime=0, XFL=0, OS=255 (unknown).
    out.extend_from_slice(&[0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff]);
    out.extend_from_slice(&body);
    out.extend_from_slice(&gzip_crc32(raw).to_le_bytes());
    out.extend_from_slice(&(raw.len() as u32).to_le_bytes());
    out
}

/// CRC32 (gzip/PNG polynomial) for the gzip member trailer.
#[cfg(feature = "gzip")]
fn gzip_crc32(data: &[u8]) -> u32 {
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

/// Byte-shuffle the big-endian integer stream for `GZIP_2` (the inverse of
/// [`gzip2_unshuffle`]): emit all the most-significant bytes first, then the next, etc.
#[cfg(feature = "gzip")]
fn gzip2_shuffle(raw: &[u8], bytepix: usize) -> Vec<u8> {
    if bytepix <= 1 {
        return raw.to_vec();
    }
    let n = raw.len() / bytepix;
    let mut out = vec![0u8; n * bytepix];
    for (i, chunk) in raw.chunks_exact(bytepix).enumerate().take(n) {
        for (b, &byte) in chunk.iter().enumerate() {
            out[b * n + i] = byte;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// High-level tile-compression (the public encode entry point)
// ---------------------------------------------------------------------------

/// Options controlling how [`ImageData::compress`] / [`compress_image`] encode an
/// image into a tile-compressed BINTABLE HDU.
#[derive(Debug, Clone)]
pub struct CompressOptions {
    /// Compression algorithm (`ZCMPTYPE`). `RICE_1`, `GZIP_1`, `GZIP_2` are supported
    /// for encoding; `PLIO_1`/`HCOMPRESS_1`/`NOCOMPRESS` are rejected.
    pub algorithm: CompressionType,
    /// Tile shape (`ZTILEn`), axis-1 first. `None` selects cfitsio's default: one row
    /// per tile (`ZTILE1 = NAXIS1`, all other axes = 1).
    pub tile: Option<Vec<usize>>,
    /// For float images: quantization parameter `q` (controls the per-tile scale,
    /// `scale = range / (2*q*sigma_estimate)`-style). `None` ⇒ lossless float storage
    /// (raw big-endian floats, only valid with the GZIP codecs). For RICE on floats a
    /// quantization value is required.
    pub quantize: Option<f64>,
    /// Subtractive-dithering method for quantized floats (`ZQUANTIZ`). Ignored for the
    /// lossless and integer paths.
    pub dither: Quantize,
    /// Dither seed (`ZDITHER0`, 1..=10000). `None` ⇒ 1.
    pub dither_seed: Option<i64>,
    /// RICE block size (`BLOCKSIZE`); cfitsio default 32.
    pub blocksize: usize,
}

impl Default for CompressOptions {
    fn default() -> Self {
        CompressOptions {
            algorithm: CompressionType::Rice1,
            tile: None,
            quantize: Some(4.0),
            dither: Quantize::SubtractiveDither1,
            dither_seed: Some(1),
            blocksize: 32,
        }
    }
}

/// Default tile shape: one row per tile (`ZTILE1 = NAXIS1`, others 1).
fn default_tile(axes: &[usize]) -> Vec<usize> {
    let mut t = vec![1usize; axes.len()];
    if let Some(first) = axes.first() {
        t[0] = *first;
    }
    t
}

/// Extract the flat integer samples of an integer image as `i64` (axis-1 fastest, the
/// native FITS storage order), along with the `BYTEPIX` to feed RICE.
fn integer_samples(image: &crate::image_data::ImageData) -> Result<(Vec<i64>, usize)> {
    use crate::image_data::PixelData;
    Ok(match &image.pixels {
        PixelData::U8(v) => (v.iter().map(|&x| x as i64).collect(), 1),
        PixelData::I16(v) => (v.iter().map(|&x| x as i64).collect(), 2),
        PixelData::I32(v) => (v.iter().map(|&x| x as i64).collect(), 4),
        PixelData::I64(_) => {
            return Err(Error::CompressionError(
                "RICE_1/GZIP tile compression of 64-bit integer images is not supported".into(),
            ))
        }
        PixelData::F32(_) | PixelData::F64(_) => {
            return Err(Error::CompressionError(
                "internal: float image routed to the integer encoder".into(),
            ))
        }
    })
}

/// Gather one tile's flat samples (axis-1 fastest) from the full image buffer.
fn gather_tile(
    full: &[i64],
    image_dims: &[usize],
    ztile: &[usize],
    tile_dims: &[usize],
    coords: &[usize],
) -> Vec<i64> {
    let ndim = image_dims.len();
    let mut img_stride = vec![1usize; ndim];
    for axis in 1..ndim {
        img_stride[axis] = img_stride[axis - 1] * image_dims[axis - 1];
    }
    let mut origin = 0usize;
    for axis in 0..ndim {
        origin += coords[axis] * ztile[axis].max(1) * img_stride[axis];
    }
    let tile_npix: usize = tile_dims.iter().product();
    let mut out = Vec::with_capacity(tile_npix);
    let mut tcoord = vec![0usize; ndim];
    for _ in 0..tile_npix {
        let mut src = origin;
        for axis in 0..ndim {
            src += tcoord[axis] * img_stride[axis];
        }
        out.push(full[src]);
        for axis in 0..ndim {
            tcoord[axis] += 1;
            if tcoord[axis] < tile_dims[axis] {
                break;
            }
            tcoord[axis] = 0;
        }
    }
    out
}

/// Tile dimensions at grid coordinates `coords`, accounting for edge truncation.
fn tile_dims_at(image_dims: &[usize], ztile: &[usize], coords: &[usize]) -> Vec<usize> {
    let mut dims = Vec::with_capacity(coords.len());
    for (axis, &c) in coords.iter().enumerate() {
        let n = image_dims[axis];
        let t = ztile[axis].max(1);
        let start = c * t;
        dims.push((n - start).min(t));
    }
    dims
}

/// Encode one integer tile's samples into its compressed byte stream per `algorithm`.
fn encode_int_tile(
    samples: &[i64],
    algorithm: CompressionType,
    bytepix: usize,
    blocksize: usize,
) -> Result<Vec<u8>> {
    match algorithm {
        CompressionType::Rice1 => Ok(match bytepix {
            1 => rice_compress_i8(
                &samples.iter().map(|&v| v as u8).collect::<Vec<_>>(),
                blocksize,
            ),
            2 => rice_compress_i16(
                &samples.iter().map(|&v| v as i16).collect::<Vec<_>>(),
                blocksize,
            ),
            4 => rice_compress_i32(
                &samples.iter().map(|&v| v as i32).collect::<Vec<_>>(),
                blocksize,
            ),
            other => {
                return Err(Error::CompressionError(format!(
                    "RICE_1 BYTEPIX={other} not supported (expected 1, 2, or 4)"
                )))
            }
        }),
        CompressionType::Gzip1 | CompressionType::Gzip2 => {
            #[cfg(feature = "gzip")]
            {
                let raw = int_samples_to_be(samples, bytepix);
                let payload = if algorithm == CompressionType::Gzip2 {
                    gzip2_shuffle(&raw, bytepix)
                } else {
                    raw
                };
                Ok(gzip_deflate(&payload))
            }
            #[cfg(not(feature = "gzip"))]
            {
                let _ = (samples, bytepix);
                Err(Error::UnsupportedCompression(
                    "GZIP_1/GZIP_2 encode requires the `gzip` feature (miniz_oxide)".into(),
                ))
            }
        }
        CompressionType::Plio1 | CompressionType::Hcompress1 | CompressionType::NoCompress => {
            Err(Error::UnsupportedCompression(format!(
                "{algorithm:?} encoding is not supported (only RICE_1, GZIP_1, GZIP_2)"
            )))
        }
    }
}

/// Serialise integer samples to big-endian bytes of width `bytepix` (used by GZIP).
#[cfg(feature = "gzip")]
fn int_samples_to_be(samples: &[i64], bytepix: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(samples.len() * bytepix);
    for &v in samples {
        match bytepix {
            1 => out.push(v as u8),
            2 => out.extend_from_slice(&(v as i16).to_be_bytes()),
            4 => out.extend_from_slice(&(v as i32).to_be_bytes()),
            8 => out.extend_from_slice(&v.to_be_bytes()),
            _ => out.extend_from_slice(&(v as i32).to_be_bytes()),
        }
    }
    out
}

/// Build the `Z*` driver keywords for a compressed-image BINTABLE, in the exact order
/// cfitsio's `fpack` writes them.
///
/// The order is load-bearing: `funpack`/cfitsio reconstruct the uncompressed image
/// header by walking the cards in order and translating `ZBITPIX`/`ZNAXIS`/... into
/// `BITPIX`/`NAXIS`/..., emitting the leading `XTENSION` card from `ZTENSION`. If
/// `ZBITPIX` appears before `ZTENSION`, the rebuilt header starts with `BITPIX` and
/// cfitsio rejects it ("1st key not SIMPLE or XTENSION"). The accepted order is:
/// `ZIMAGE`, `ZTILEn`, `ZCMPTYPE`, `ZNAMEi`/`ZVALi`, `EXTNAME`, `ZTENSION`, `ZBITPIX`,
/// `ZNAXIS`, `ZNAXISn`, `ZPCOUNT`, `ZGCOUNT`. Codec parameter pairs (`zparams`) are
/// emitted between `ZCMPTYPE` and `EXTNAME`.
fn build_z_header(
    header: &mut Header,
    image: &crate::image_data::ImageData,
    ztile: &[usize],
    algorithm: CompressionType,
    zparams: &[(&str, i64)],
) {
    use crate::keyword::HeaderValue;
    let zbitpix = image.bitpix().to_i64();
    let znaxis = image.axes.len();

    header.set(
        "ZIMAGE",
        HeaderValue::Logical(true),
        Some("extension contains compressed image"),
    );
    for (i, &t) in ztile.iter().enumerate() {
        header.set(
            &format!("ZTILE{}", i + 1),
            HeaderValue::Integer(t as i64),
            Some("size of tiles to be compressed"),
        );
    }
    let zcmptype = match algorithm {
        CompressionType::Rice1 => "RICE_1",
        CompressionType::Gzip1 => "GZIP_1",
        CompressionType::Gzip2 => "GZIP_2",
        CompressionType::Hcompress1 => "HCOMPRESS_1",
        CompressionType::Plio1 => "PLIO_1",
        CompressionType::NoCompress => "NOCOMPRESS",
    };
    header.set(
        "ZCMPTYPE",
        HeaderValue::String(zcmptype.into()),
        Some("compression algorithm"),
    );
    for (i, (name, val)) in zparams.iter().enumerate() {
        header.set(
            &format!("ZNAME{}", i + 1),
            HeaderValue::String((*name).into()),
            None,
        );
        header.set(&format!("ZVAL{}", i + 1), HeaderValue::Integer(*val), None);
    }
    header.set(
        "EXTNAME",
        HeaderValue::String("COMPRESSED_IMAGE".into()),
        None,
    );
    // The compressed image always reconstructs as an IMAGE extension; cfitsio reads the
    // original extension type from `ZTENSION` and emits the leading `XTENSION` card.
    header.set(
        "ZTENSION",
        HeaderValue::String("IMAGE".into()),
        Some("Image extension"),
    );
    header.set(
        "ZBITPIX",
        HeaderValue::Integer(zbitpix),
        Some("data type of original image"),
    );
    header.set(
        "ZNAXIS",
        HeaderValue::Integer(znaxis as i64),
        Some("dimension of original image"),
    );
    for (i, &n) in image.axes.iter().enumerate() {
        header.set(
            &format!("ZNAXIS{}", i + 1),
            HeaderValue::Integer(n as i64),
            None,
        );
    }
    header.set("ZPCOUNT", HeaderValue::Integer(0), Some("original PCOUNT"));
    header.set("ZGCOUNT", HeaderValue::Integer(1), Some("original GCOUNT"));
}

/// Compress an [`ImageData`] into a tile-compressed BINTABLE [`Hdu`] (`ZIMAGE = T`).
///
/// The returned HDU's `data` is [`HduData::BinTable`](crate::hdu::HduData::BinTable)
/// with a `COMPRESSED_DATA` (`1PB`) variable-length-array column holding one
/// RICE/GZIP-encoded tile per row, plus the full set of `Z*` driver keywords. It is
/// ready to hand to [`FitsFile::push_extension`](crate::fits::FitsFile::push_extension).
///
/// - Integer images (`ZBITPIX` 8/16/32) compress losslessly with `RICE_1`/`GZIP_1`/
///   `GZIP_2`.
/// - Float images (`ZBITPIX` -32/-64) compress losslessly with `GZIP_1`/`GZIP_2` when
///   `opts.quantize` is `None` (raw big-endian floats, `ZQUANTIZ` = `NONE`), or lossily
///   with quantization + optional subtractive dithering otherwise (`RICE_1`/GZIP of the
///   quantized 32-bit integers, with per-tile `ZSCALE`/`ZZERO` columns).
///
/// `PLIO_1`, `HCOMPRESS_1`, `NOCOMPRESS`, and 64-bit integer images are rejected with
/// [`Error::UnsupportedCompression`]/[`Error::CompressionError`].
pub fn compress_image(
    image: &crate::image_data::ImageData,
    opts: &CompressOptions,
) -> Result<crate::hdu::Hdu> {
    use crate::bintable::{BinColumnType, BinTableBuilder};
    use crate::hdu::{Hdu, HduData};

    match opts.algorithm {
        CompressionType::Rice1 | CompressionType::Gzip1 | CompressionType::Gzip2 => {}
        other => {
            return Err(Error::UnsupportedCompression(format!(
                "{other:?} encoding is not supported (only RICE_1, GZIP_1, GZIP_2)"
            )))
        }
    }

    if image.axes.is_empty() {
        return Err(Error::CompressionError(
            "cannot compress an image with zero axes".into(),
        ));
    }
    let expected_npix: usize = image.axes.iter().product();
    if image.pixels.len() != expected_npix {
        return Err(Error::CompressionError(format!(
            "image has {} pixels but axes imply {expected_npix}",
            image.pixels.len()
        )));
    }

    let ztile = opts.tile.clone().unwrap_or_else(|| default_tile(&image.axes));
    if ztile.len() != image.axes.len() {
        return Err(Error::CompressionError(format!(
            "tile shape has {} axes but image has {}",
            ztile.len(),
            image.axes.len()
        )));
    }

    let is_float = image.bitpix().to_i64() < 0;
    if is_float {
        return compress_float_image(image, opts, &ztile);
    }

    // ---- integer path ----
    let (full, bytepix) = integer_samples(image)?;
    let blocksize = opts.blocksize.max(1);

    // Tile grid.
    let tiles_per_axis: Vec<usize> = image
        .axes
        .iter()
        .zip(&ztile)
        .map(|(&n, &t)| n.div_ceil(t.max(1)))
        .collect();
    let num_tiles: usize = tiles_per_axis.iter().product();

    let mut builder = BinTableBuilder::new().add_column("COMPRESSED_DATA", BinColumnType::VarP('B'));
    let mut max_elems = 0usize;
    for tile_index in 0..num_tiles {
        let coords = unravel(tile_index, &tiles_per_axis);
        let tdims = tile_dims_at(&image.axes, &ztile, &coords);
        let samples = gather_tile(&full, &image.axes, &ztile, &tdims, &coords);
        let enc = encode_int_tile(&samples, opts.algorithm, bytepix, blocksize)?;
        max_elems = max_elems.max(enc.len());
        builder = builder.push_row(|r| {
            r.write_var_p(enc.len() as i32, |heap| heap.extend_from_slice(&enc))
        });
    }
    let table = builder.build();

    let mut header = Header::new();
    let zparams: Vec<(&str, i64)> = if opts.algorithm == CompressionType::Rice1 {
        vec![("BLOCKSIZE", blocksize as i64), ("BYTEPIX", bytepix as i64)]
    } else {
        Vec::new()
    };
    build_z_header(&mut header, image, &ztile, opts.algorithm, &zparams);
    let header = finalize_compressed_hdu(header, &table, max_elems);

    Ok(Hdu::new(header, HduData::BinTable(table)))
}

/// Compress a float image (`ZBITPIX` -32/-64). See [`compress_image`].
fn compress_float_image(
    image: &crate::image_data::ImageData,
    opts: &CompressOptions,
    ztile: &[usize],
) -> Result<crate::hdu::Hdu> {
    use crate::bintable::{BinColumnType, BinTableBuilder};
    use crate::hdu::{Hdu, HduData};
    use crate::image_data::PixelData;
    use crate::keyword::HeaderValue;

    // Used only by the (feature-gated) lossless-float GZIP branch.
    #[cfg_attr(not(feature = "gzip"), allow(unused_variables))]
    let is_f64 = image.bitpix().to_i64() == -64;
    let floats: Vec<f64> = match &image.pixels {
        PixelData::F32(v) => v.iter().map(|&x| x as f64).collect(),
        PixelData::F64(v) => v.clone(),
        _ => unreachable!("compress_float_image called on non-float image"),
    };

    let tiles_per_axis: Vec<usize> = image
        .axes
        .iter()
        .zip(ztile)
        .map(|(&n, &t)| n.div_ceil(t.max(1)))
        .collect();
    let num_tiles: usize = tiles_per_axis.iter().product();
    let blocksize = opts.blocksize.max(1);

    // ---- lossless float storage (raw big-endian floats, GZIP only) ----
    if opts.quantize.is_none() {
        // Lossless float storage uses GZIP over the raw big-endian floats. Only GZIP_1
        // is supported: cfitsio's float byte-shuffle (GZIP_2) is not applied on the
        // raw-float fallback path that the decoder reads, so GZIP_2 here would not
        // round-trip. (GZIP_2 *is* supported for integer images.)
        match opts.algorithm {
            CompressionType::Gzip1 => {}
            other => {
                return Err(Error::UnsupportedCompression(format!(
                    "lossless float compression (quantize=None) needs GZIP_1, not {other:?}"
                )))
            }
        }
        #[cfg(not(feature = "gzip"))]
        {
            return Err(Error::UnsupportedCompression(
                "GZIP float encode requires the `gzip` feature (miniz_oxide)".into(),
            ));
        }
        #[cfg(feature = "gzip")]
        {
            let mut builder =
                BinTableBuilder::new().add_column("COMPRESSED_DATA", BinColumnType::VarP('B'));
            let mut max_elems = 0usize;
            for tile_index in 0..num_tiles {
                let coords = unravel(tile_index, &tiles_per_axis);
                let tdims = tile_dims_at(&image.axes, ztile, &coords);
                let tile = gather_tile_f64(&floats, &image.axes, ztile, &tdims, &coords);
                // Raw big-endian floats, GZIP_1 (no shuffle) — read back by the decoder's
                // unquantized-tile fallback path.
                let mut raw = Vec::with_capacity(tile.len() * if is_f64 { 8 } else { 4 });
                for &v in &tile {
                    if is_f64 {
                        raw.extend_from_slice(&v.to_be_bytes());
                    } else {
                        raw.extend_from_slice(&(v as f32).to_be_bytes());
                    }
                }
                let enc = gzip_deflate(&raw);
                max_elems = max_elems.max(enc.len());
                builder = builder
                    .push_row(|r| r.write_var_p(enc.len() as i32, |heap| heap.extend_from_slice(&enc)));
            }
            let table = builder.build();
            let mut header = Header::new();
            build_z_header(&mut header, image, ztile, opts.algorithm, &[]);
            header.set("ZQUANTIZ", HeaderValue::String("NONE".into()), Some("no dithering"));
            let header = finalize_compressed_hdu(header, &table, max_elems);
            return Ok(Hdu::new(header, HduData::BinTable(table)));
        }
    }

    // ---- lossy quantized float storage ----
    let q = opts.quantize.unwrap_or(4.0);
    let dither = opts.dither;
    let dithered = matches!(
        dither,
        Quantize::SubtractiveDither1 | Quantize::SubtractiveDither2
    );
    let zdither0 = opts.dither_seed.unwrap_or(1);

    let mut builder = BinTableBuilder::new()
        .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
        .add_column("ZSCALE", BinColumnType::D64(1))
        .add_column("ZZERO", BinColumnType::D64(1));
    let mut max_elems = 0usize;

    for tile_index in 0..num_tiles {
        let coords = unravel(tile_index, &tiles_per_axis);
        let tdims = tile_dims_at(&image.axes, ztile, &coords);
        let tile = gather_tile_f64(&floats, &image.axes, ztile, &tdims, &coords);

        let (scale, zero) = choose_scale_zero(&tile, q);
        let q_ints = quantize_tile(&tile, scale, zero, dither, zdither0, tile_index);
        let enc = encode_int_tile(&q_ints, opts.algorithm, 4, blocksize)?;
        max_elems = max_elems.max(enc.len());
        builder = builder.push_row(|r| {
            r.write_var_p(enc.len() as i32, |heap| heap.extend_from_slice(&enc));
            r.write_f64(scale);
            r.write_f64(zero);
        });
    }
    let table = builder.build();

    let mut header = Header::new();
    let zparams: Vec<(&str, i64)> = if opts.algorithm == CompressionType::Rice1 {
        vec![("BLOCKSIZE", blocksize as i64), ("BYTEPIX", 4)]
    } else {
        Vec::new()
    };
    build_z_header(&mut header, image, ztile, opts.algorithm, &zparams);
    let zquantiz = match dither {
        Quantize::SubtractiveDither1 => "SUBTRACTIVE_DITHER_1",
        Quantize::SubtractiveDither2 => "SUBTRACTIVE_DITHER_2",
        Quantize::None => "NO_DITHER",
    };
    header.set("ZQUANTIZ", HeaderValue::String(zquantiz.into()), Some("quantization method"));
    if dithered {
        header.set("ZDITHER0", HeaderValue::Integer(zdither0), Some("dithering offset seed"));
    }
    header.set("ZBLANK", HeaderValue::Integer(NULL_VALUE), Some("null value"));
    let header = finalize_compressed_hdu(header, &table, max_elems);

    Ok(Hdu::new(header, HduData::BinTable(table)))
}

/// Gather one tile's flat floats (axis-1 fastest) from the full image buffer.
fn gather_tile_f64(
    full: &[f64],
    image_dims: &[usize],
    ztile: &[usize],
    tile_dims: &[usize],
    coords: &[usize],
) -> Vec<f64> {
    let ndim = image_dims.len();
    let mut img_stride = vec![1usize; ndim];
    for axis in 1..ndim {
        img_stride[axis] = img_stride[axis - 1] * image_dims[axis - 1];
    }
    let mut origin = 0usize;
    for axis in 0..ndim {
        origin += coords[axis] * ztile[axis].max(1) * img_stride[axis];
    }
    let tile_npix: usize = tile_dims.iter().product();
    let mut out = Vec::with_capacity(tile_npix);
    let mut tcoord = vec![0usize; ndim];
    for _ in 0..tile_npix {
        let mut src = origin;
        for axis in 0..ndim {
            src += tcoord[axis] * img_stride[axis];
        }
        out.push(full[src]);
        for axis in 0..ndim {
            tcoord[axis] += 1;
            if tcoord[axis] < tile_dims[axis] {
                break;
            }
            tcoord[axis] = 0;
        }
    }
    out
}

/// Choose a per-tile linear quantization `(scale, zero)` for float values.
///
/// We do not replicate cfitsio's noise-based heuristic; any scale that round-trips
/// within tolerance is acceptable (per the task). We map the finite value range into
/// the available integer range with `q` quantization levels per unit range so the
/// quantization step is `scale = range / (q * 2^16)` (clamped to a sane minimum), with
/// `zero` at the tile mean. The reconstructed value is `zero + scale * q_int`, matching
/// the decoder's `unquantize`.
fn choose_scale_zero(tile: &[f64], q: f64) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut sum = 0.0f64;
    let mut count = 0usize;
    for &v in tile {
        if v.is_finite() {
            min = min.min(v);
            max = max.max(v);
            sum += v;
            count += 1;
        }
    }
    if count == 0 {
        // All-NaN tile: scale 1, zero 0 (every pixel becomes the blank sentinel).
        return (1.0, 0.0);
    }
    let range = (max - min).abs();
    let q = if q > 0.0 { q } else { 4.0 };
    // Target ~ (q * 65536) quantization levels across the range; floor the step so a
    // flat tile still quantizes losslessly enough.
    let mut scale = range / (q * 65536.0);
    if !(scale.is_finite()) || scale <= 0.0 {
        scale = 1.0;
    }
    let zero = sum / count as f64;
    (scale, zero)
}

/// Quantize a tile's floats to 32-bit integers with cfitsio-compatible subtractive
/// dithering. The forward transform inverts the decoder's `unquantize`:
/// `q = round((value - zero)/scale + dither - 0.5)`, NaNs map to `NULL_VALUE`, and
/// (for `SUBTRACTIVE_DITHER_2`) exact zeros map to `ZERO_VALUE`.
fn quantize_tile(
    tile: &[f64],
    scale: f64,
    zero: f64,
    dither: Quantize,
    zdither0: i64,
    tile_index: usize,
) -> Vec<i64> {
    let dithered = matches!(
        dither,
        Quantize::SubtractiveDither1 | Quantize::SubtractiveDither2
    );
    let dither2 = dither == Quantize::SubtractiveDither2;

    let mut iseed = ((tile_index as i64 + zdither0 - 1).rem_euclid(N_RANDOM as i64)) as usize;
    let mut nextrand = (fits_rand_value(iseed) * 500.0) as usize;

    let mut out = Vec::with_capacity(tile.len());
    for &v in tile {
        let qi = if v.is_nan() {
            NULL_VALUE
        } else if dither2 && v == 0.0 {
            ZERO_VALUE
        } else if dithered {
            // Inverse of value = (q - r + 0.5)*scale + zero, computed with the FMA the
            // decoder uses so the round-trip is exact:
            //   q = round((value - zero)/scale + r - 0.5)
            let r = fits_rand_value(nextrand) as f64;
            let t = (v - zero) / scale + r - 0.5;
            t.round() as i64
        } else {
            let t = (v - zero) / scale;
            t.round() as i64
        };
        out.push(qi);

        if dithered {
            nextrand += 1;
            if nextrand == N_RANDOM {
                iseed += 1;
                if iseed == N_RANDOM {
                    iseed = 0;
                }
                nextrand = (fits_rand_value(iseed) * 500.0) as usize;
            }
        }
    }
    out
}

/// Finalise a compressed-image HDU: fill the BINTABLE keywords, then fix up `TFORM1` to
/// the `1PB(maxlen)` form cfitsio writes (so `funpack` knows the max VLA element count),
/// and append a descriptive `EXTNAME`.
fn finalize_compressed_hdu(mut header: Header, table: &BinTable, max_elems: usize) -> Header {
    use crate::keyword::HeaderValue;
    // Materialise the mandatory BINTABLE keywords (XTENSION/NAXISn/PCOUNT/TFIELDS/...)
    // from the table, then layer the Z* keywords (already in `header`) on top. We do
    // this by filling a fresh header and merging Z* afterwards so ordering matches the
    // cfitsio convention (BINTABLE structural keywords first, then ZIMAGE...).
    let mut out = Header::new();
    table.fill_header(&mut out);
    // TFORM1 with explicit max length, matching cfitsio (`1PB(<maxlen>)`).
    out.set(
        "TFORM1",
        HeaderValue::String(format!("1PB({max_elems})")),
        Some("variable length array"),
    );
    // Append all Z* / quantization / extra keywords from the staged header.
    for kw in header.keywords.drain(..) {
        out.keywords.push(kw);
    }
    out
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

// ---------------------------------------------------------------------------
// PLIO_1 decompression (IRAF pixel-list run-length encoding)
// ---------------------------------------------------------------------------
//
// Port of cfitsio `pliocomp.c` `pl_l2pi` (D. Tody, NRAO; the IRAF PLIO scheme).
// A tile is stored as a *line list*: a stream of 16-bit instruction words. Each
// word is `(opcode << 12) | data`, where `opcode` is the high nibble (0..15) and
// `data` is the low 12 bits. The list reconstructs a 1-D run of `npix`
// non-negative integers ("pixels"), tracking a running "high value" `pv`.
//
// In the FITS tiled-image convention the line list is stored in the
// `COMPRESSED_DATA` column as a `1PI`/`1QI` variable-length array of *big-endian*
// 16-bit values (cfitsio reads it with `fits_read_col(TSHORT, ...)`), so the
// caller hands us the already-byte-swapped `i16` words.
//
// Opcodes (after cfitsio's `opcode = word/4096; data = word & 4095`; the C switch
// is on `opcode+1`, reproduced here on `opcode` directly):
//   * 0  (Z run)         — `data` zeros.
//   * 1  (HD, high+data) — set high value: `pv = (next_word << 12) + data`; the
//                          following word is consumed as the high bits.
//   * 2  (IH, inc high)  — `pv += data`.
//   * 3  (DH, dec high)  — `pv -= data`.
//   * 4  (PD, pixel data)— `data` pixels at value `pv`.
//   * 5  (PS, pixel set) — `data` zeros, then the *last* pixel of the run is set
//                          to `pv` (single non-zero pixel at the end of the span).
//   * 6  (IS, inc set)   — `pv += data`, then emit one pixel at `pv`.
//   * 7  (DS, dec set)   — `pv -= data`, then emit one pixel at `pv`.
// (cfitsio's switch maps opcodes 0,4,5 all through the "run of length `data`"
// label L160; opcode 0 is a pure zero run, 4 fills with `pv`, 5 fills with zeros
// then sets the last element. Opcodes 1..3 adjust `pv` without emitting span
// pixels, and 6/7 adjust `pv` and emit a single pixel.)

/// Decompress a PLIO_1 (IRAF pixel-list RLE) tile into `nvals` `i32` values.
///
/// `words` is the line-list instruction stream (already decoded from the
/// big-endian `1PI` VLA into native `i16`). Pixels are non-negative; any pixels
/// not explicitly written are zero.
pub fn plio_decompress(words: &[i16], nvals: usize) -> Result<Vec<i32>> {
    // Treat the instruction words as unsigned 16-bit (the high nibble is the
    // opcode; PLIO values never use the sign bit).
    let w: Vec<u32> = words.iter().map(|&v| (v as u16) as u32).collect();

    let mut out = vec![0i32; nvals];
    if nvals == 0 {
        return Ok(out);
    }
    if w.len() < 3 {
        return Err(Error::CompressionError(
            "PLIO_1 line list too short (need >= 3 header words)".into(),
        ));
    }

    // List length and index of the first instruction word (cfitsio uses 1-based
    // indices via `--ll_src`; here `w[0]` == cfitsio `ll_src[1]`, so cfitsio
    // `ll_src[k]` == `w[k-1]`).
    //
    //   if ll_src[3] > 0:   lllen = ll_src[3];                 llfirt = 4
    //   else:               lllen = (ll_src[5]<<15)+ll_src[4]; llfirt = ll_src[2]+1
    // The `ll_src[3] > 0` test is a SIGNED 16-bit comparison: the standard PLIO
    // header stores -100 there, so the long form (else branch) is the usual path.
    let (lllen, llfirt) = if words[2] > 0 {
        (w[2] as usize, 4usize)
    } else {
        if w.len() < 5 {
            return Err(Error::CompressionError(
                "PLIO_1 long-form line list too short".into(),
            ));
        }
        let len = ((w[4] << 15) + w[3]) as usize;
        (len, (w[1] as usize) + 1)
    };

    if lllen == 0 {
        return Ok(out); // empty list => all zeros
    }

    // Pixel-output cursor `op` (1-based in cfitsio; `op==1` means out[0]).
    let mut op: usize = 1; // index into 1-based px_dst
    let mut x1: i64 = 1; // current x position (1-based)
    let mut pv: i64 = 1; // running "high value"
    let xs: i64 = 1; // starting index (cfitsio passes xs=1)
    let xe: i64 = nvals as i64; // end index inclusive

    let put = |out: &mut [i32], op: usize, v: i64| {
        // 1-based op -> 0-based; ignore writes past the tile (defensive).
        if op >= 1 && op <= out.len() {
            out[op - 1] = v as i32;
        }
    };

    // ip iterates cfitsio ll_src indices llfirt..=lllen (1-based) => w[ip-1].
    let mut ip = llfirt;
    let mut skipwd = false;
    'outer: while ip <= lllen {
        if ip == 0 || ip > w.len() {
            break;
        }
        if skipwd {
            skipwd = false;
            ip += 1;
            continue;
        }
        let word = w[ip - 1];
        let opcode = word / 4096;
        let data = (word & 4095) as i64;

        match opcode {
            // L160: run of length `data` (opcodes 0, 4, 5).
            0 | 4 | 5 => {
                let x2 = x1 + data - 1;
                let i1 = x1.max(xs);
                let i2 = x2.min(xe);
                let np = i2 - i1 + 1;
                if np > 0 {
                    let otop = op as i64 + np - 1;
                    if opcode == 4 {
                        for i in op..=(otop as usize) {
                            put(&mut out, i, pv);
                        }
                    } else {
                        // opcodes 0 and 5: zeros (out is already zero-initialized).
                        if opcode == 5 && i2 == x2 {
                            put(&mut out, otop as usize, pv);
                        }
                    }
                    op = (otop + 1) as usize;
                }
                x1 = x2 + 1;
            }
            // L220: set high value from this+next word.
            1 => {
                if ip >= w.len() {
                    return Err(Error::CompressionError(
                        "PLIO_1 truncated high-value instruction".into(),
                    ));
                }
                pv = ((w[ip] as i64) << 12) + data;
                skipwd = true;
            }
            // L230: increment high value.
            2 => {
                pv += data;
            }
            // L240: decrement high value.
            3 => {
                pv -= data;
            }
            // L250: increment high value and emit one pixel.
            6 => {
                pv += data;
                if x1 >= xs && x1 <= xe {
                    put(&mut out, op, pv);
                    op += 1;
                }
                x1 += 1;
            }
            // L260: decrement high value and emit one pixel.
            7 => {
                pv -= data;
                if x1 >= xs && x1 <= xe {
                    put(&mut out, op, pv);
                    op += 1;
                }
                x1 += 1;
            }
            _ => {
                // opcodes 8..15 are unused by the encoder; cfitsio falls through
                // (no-op) on them. Match that behaviour.
            }
        }

        if x1 > xe {
            break 'outer;
        }
        ip += 1;
    }

    // Remaining pixels are zero (already initialized).
    Ok(out)
}

// ---------------------------------------------------------------------------
// HCOMPRESS_1 decompression
// ---------------------------------------------------------------------------
//
// Port of cfitsio `fits_hdecompress` / `fits_hdecompress64` (the concatenation of
// R. White's hinv.c / undigitize.c / decode.c / dodecode.c / qtree_decode.c /
// qread.c / bit_input.c in STScI's hcompress distribution).
//
// Pipeline (read direction):
//   1. `decode`  — parse the byte stream header (2-byte magic 0xDD 0x99, then
//      big-endian `nx`, `ny`, `scale` as 4-byte ints, an 8-byte `sumall`, and a
//      3-byte `nbitplanes`), then `dodecode` the quadtree-coded bit planes of the
//      four image quadrants into the transform-coefficient array `a`, and finally
//      the per-element sign bits. `a[0]` is overwritten with `sumall`.
//   2. `undigitize` — multiply every coefficient by `scale` (the quantization
//      step; `scale <= 1` is a no-op, i.e. lossless).
//   3. `hinv` — inverse H-transform, expanding the coefficients back into pixels.
//
// The array `a` is indexed as `a[i*ny + j]` (ny = fast axis = FITS axis-1 width),
// so its flat order already matches the tile's axis-1-fastest pixel layout.
//
// cfitsio uses 32-bit ints for ZBITPIX 8/16 and 64-bit ints for everything else
// (incl. quantized floats), because the H-transform intermediate sums can exceed
// 32 bits. We reproduce both via a macro, using wrapping arithmetic to match C's
// 2's-complement overflow behaviour exactly.

/// Bit/byte reader for the HCOMPRESS stream (mirrors cfitsio's global
/// `nextchar` + `buffer2`/`bits_to_go` state machine, but as a struct).
struct HcInput<'a> {
    data: &'a [u8],
    nextchar: usize,
    buffer2: i32,
    bits_to_go: i32,
}

impl<'a> HcInput<'a> {
    fn new(data: &'a [u8]) -> Self {
        HcInput {
            data,
            nextchar: 0,
            buffer2: 0,
            bits_to_go: 0,
        }
    }

    #[inline]
    fn next_byte(&mut self) -> Result<i32> {
        let b = *self
            .data
            .get(self.nextchar)
            .ok_or_else(|| Error::CompressionError("HCOMPRESS: unexpected end of stream".into()))?;
        self.nextchar += 1;
        Ok(b as i32)
    }

    /// Read `n` raw bytes (no bit buffering); cfitsio `qread`.
    fn qread(&mut self, n: usize) -> Result<&[u8]> {
        let start = self.nextchar;
        let end = start + n;
        if end > self.data.len() {
            return Err(Error::CompressionError(
                "HCOMPRESS: unexpected end of stream (qread)".into(),
            ));
        }
        self.nextchar = end;
        Ok(&self.data[start..end])
    }

    /// Read a big-endian 4-byte int (cfitsio `readint`).
    fn readint(&mut self) -> Result<i32> {
        let b = self.qread(4)?;
        let mut a = b[0] as i32;
        for &x in &b[1..4] {
            a = (a << 8) + x as i32;
        }
        Ok(a)
    }

    /// Read a big-endian 8-byte long long (cfitsio `readlonglong`).
    fn readlonglong(&mut self) -> Result<i64> {
        let b = self.qread(8)?;
        let mut a = b[0] as i64;
        for &x in &b[1..8] {
            a = (a << 8) + x as i64;
        }
        Ok(a)
    }

    fn start_inputing_bits(&mut self) {
        self.bits_to_go = 0;
    }

    fn input_bit(&mut self) -> Result<i32> {
        if self.bits_to_go == 0 {
            self.buffer2 = self.next_byte()?;
            self.bits_to_go = 8;
        }
        self.bits_to_go -= 1;
        Ok((self.buffer2 >> self.bits_to_go) & 1)
    }

    fn input_nbits(&mut self, n: i32) -> Result<i32> {
        if self.bits_to_go < n {
            self.buffer2 = (self.buffer2 << 8) | self.next_byte()?;
            self.bits_to_go += 8;
        }
        self.bits_to_go -= n;
        Ok((self.buffer2 >> self.bits_to_go) & ((1 << n) - 1))
    }

    #[inline]
    fn input_nybble(&mut self) -> Result<i32> {
        self.input_nbits(4)
    }

    /// Read `n` 4-bit nybbles into `array` (cfitsio `input_nnybble`).
    fn input_nnybble(&mut self, n: usize, array: &mut [u8]) -> Result<()> {
        if n == 1 {
            array[0] = self.input_nybble()? as u8;
            return Ok(());
        }
        if self.bits_to_go == 8 {
            // Backspace to reuse the last char (cfitsio quirk).
            self.nextchar -= 1;
            self.bits_to_go = 0;
        }
        let shift1 = self.bits_to_go + 4;
        let shift2 = self.bits_to_go;
        let mut kk = 0usize;
        let mut ii = 0usize;
        if self.bits_to_go == 0 {
            while ii < n / 2 {
                self.buffer2 = (self.buffer2 << 8) | self.next_byte()?;
                array[kk] = ((self.buffer2 >> 4) & 15) as u8;
                array[kk + 1] = (self.buffer2 & 15) as u8;
                kk += 2;
                ii += 1;
            }
        } else {
            while ii < n / 2 {
                self.buffer2 = (self.buffer2 << 8) | self.next_byte()?;
                array[kk] = ((self.buffer2 >> shift1) & 15) as u8;
                array[kk + 1] = ((self.buffer2 >> shift2) & 15) as u8;
                kk += 2;
                ii += 1;
            }
        }
        if ii * 2 != n {
            array[n - 1] = self.input_nybble()? as u8;
        }
        Ok(())
    }

    /// Huffman decode of a 4-bit code (cfitsio `input_huffman`).
    fn input_huffman(&mut self) -> Result<i32> {
        let mut c = self.input_nbits(3)?;
        if c < 4 {
            return Ok(1 << c);
        }
        c = self.input_bit()? | (c << 1);
        if c < 13 {
            match c {
                8 => return Ok(3),
                9 => return Ok(5),
                10 => return Ok(10),
                11 => return Ok(12),
                12 => return Ok(15),
                _ => {}
            }
        }
        c = self.input_bit()? | (c << 1);
        if c < 31 {
            match c {
                26 => return Ok(6),
                27 => return Ok(7),
                28 => return Ok(9),
                29 => return Ok(11),
                30 => return Ok(13),
                _ => {}
            }
        }
        c = self.input_bit()? | (c << 1);
        if c == 62 {
            Ok(0)
        } else {
            Ok(14)
        }
    }
}

/// Expand 4-bit quadtree values from `a[(nx+1)/2,(ny+1)/2]` into `b[nx,ny]`
/// (2x2 per value); cfitsio `qtree_copy`. `a` and `b` are the same buffer here, so
/// we operate in place exactly as the C does (iterating from the end first).
fn qtree_copy(buf: &mut [u8], nx: usize, ny: usize, n: usize) {
    let nx2 = nx.div_ceil(2);
    let ny2 = ny.div_ceil(2);
    // Copy 4-bit values to b, from the end (a,b same array).
    // k is index of a[i,j]; s00 is index of b[2*i,2*j].
    let mut k = (ny2 * (nx2 - 1) + ny2 - 1) as isize;
    for i in (0..nx2).rev() {
        let mut s00 = (2 * (n * i + ny2 - 1)) as isize;
        for _j in (0..ny2).rev() {
            buf[s00 as usize] = buf[k as usize];
            k -= 1;
            s00 -= 2;
        }
    }
    // Expand each 2x2 block. Mapping: bit3->b[s00], bit2->b[s00+1],
    // bit1->b[s10], bit0->b[s10+1] where s10 = s00+n.
    let mut i = 0usize;
    while i + 1 < nx {
        let mut s00 = n * i;
        let s10base = s00 + n;
        let mut s10 = s10base;
        let mut j = 0usize;
        while j + 1 < ny {
            let v = buf[s00];
            buf[s10 + 1] = v & 1;
            buf[s10] = (v >> 1) & 1;
            buf[s00 + 1] = (v >> 2) & 1;
            buf[s00] = (v >> 3) & 1;
            s00 += 2;
            s10 += 2;
            j += 2;
        }
        if j < ny {
            // odd row length
            let v = buf[s00];
            buf[s10] = (v >> 1) & 1;
            buf[s00] = (v >> 3) & 1;
        }
        i += 2;
    }
    if i < nx {
        // odd column length: last row, s10 off edge
        let mut s00 = n * i;
        let mut j = 0usize;
        while j + 1 < ny {
            let v = buf[s00];
            buf[s00 + 1] = (v >> 2) & 1;
            buf[s00] = (v >> 3) & 1;
            s00 += 2;
            j += 2;
        }
        if j < ny {
            let v = buf[s00];
            buf[s00] = (v >> 3) & 1;
        }
    }
}

/// One quadtree expansion step (cfitsio `qtree_expand`): copy+expand then read a
/// fresh Huffman code into every non-zero element (scanning from the end).
fn qtree_expand(input: &mut HcInput, buf: &mut [u8], nx: usize, ny: usize) -> Result<()> {
    qtree_copy(buf, nx, ny, ny);
    for i in (0..nx * ny).rev() {
        if buf[i] != 0 {
            buf[i] = input.input_huffman()? as u8;
        }
    }
    Ok(())
}

macro_rules! impl_hdecompress {
    ($name:ident, $t:ty, $bitins:ident, $read_bdirect:ident, $qtree_decode:ident,
     $dodecode:ident, $hinv:ident, $undigitize:ident, $unshuffle:ident) => {
        /// Distribute even/odd interleaved coefficients (cfitsio `unshuffle`).
        /// `offset` is the base index into `a`; pointer arithmetic is done in
        /// `isize` so the trailing (unused) decrements can go negative as in C.
        fn $unshuffle(a: &mut [$t], offset: usize, n: usize, n2: usize, tmp: &mut [$t]) {
            let base = offset as isize;
            let n2i = n2 as isize;
            let nhalf = (n + 1) >> 1;
            // copy 2nd half of array to tmp
            let mut p1 = base + n2i * nhalf as isize;
            for slot in tmp.iter_mut().take(n - nhalf) {
                *slot = a[p1 as usize];
                p1 += n2i;
            }
            // distribute 1st half to even elements (descending)
            let mut p2 = base + n2i * (nhalf as isize - 1);
            let mut p1e = base + ((n2i * (nhalf as isize - 1)) << 1);
            let mut i = nhalf as isize - 1;
            while i >= 0 {
                a[p1e as usize] = a[p2 as usize];
                p2 -= n2i;
                p1e -= n2i + n2i;
                i -= 1;
            }
            // distribute 2nd half (tmp) to odd elements
            let mut p1o = base + n2i;
            let mut pt = 0usize;
            let mut i = 1usize;
            while i < n {
                a[p1o as usize] = tmp[pt];
                p1o += n2i + n2i;
                pt += 1;
                i += 2;
            }
        }

        /// Insert expanded 4-bit codes from `aa[(nx+1)/2,(ny+1)/2]` into bitplane
        /// `bit` of `b[nx,ny]` (cfitsio `qtree_bitins`).
        fn $bitins(aa: &[u8], nx: usize, ny: usize, b: &mut [$t], n: usize, bit: i32) {
            let plane_val: $t = (1 as $t) << bit;
            let mut k = 0usize;
            let mut i = 0usize;
            while i + 1 < nx {
                let s00 = n * i;
                let mut s00 = s00;
                let mut j = 0usize;
                while j + 1 < ny {
                    let v = aa[k];
                    if v & 1 != 0 {
                        b[s00 + n + 1] |= plane_val;
                    }
                    if v & 2 != 0 {
                        b[s00 + n] |= plane_val;
                    }
                    if v & 4 != 0 {
                        b[s00 + 1] |= plane_val;
                    }
                    if v & 8 != 0 {
                        b[s00] |= plane_val;
                    }
                    s00 += 2;
                    k += 1;
                    j += 2;
                }
                if j < ny {
                    // odd row: s00+1, s10+1 off edge -> only bits 1 (s10) and 3 (s00)
                    let v = aa[k];
                    if v & 2 != 0 {
                        b[s00 + n] |= plane_val;
                    }
                    if v & 8 != 0 {
                        b[s00] |= plane_val;
                    }
                    k += 1;
                }
                i += 2;
            }
            if i < nx {
                // odd column: last row, s10 off edge -> bits 2 (s00+1) and 3 (s00)
                let mut s00 = n * i;
                let mut j = 0usize;
                while j + 1 < ny {
                    let v = aa[k];
                    if v & 4 != 0 {
                        b[s00 + 1] |= plane_val;
                    }
                    if v & 8 != 0 {
                        b[s00] |= plane_val;
                    }
                    s00 += 2;
                    k += 1;
                    j += 2;
                }
                if j < ny {
                    // corner: only bit 3 (s00)
                    let v = aa[k];
                    if v & 8 != 0 {
                        b[s00] |= plane_val;
                    }
                    k += 1;
                }
            }
            let _ = k;
        }

        /// Read a directly-stored bit plane and insert it (cfitsio `read_bdirect`).
        fn $read_bdirect(
            input: &mut HcInput,
            a: &mut [$t],
            n: usize,
            nqx: usize,
            nqy: usize,
            scratch: &mut [u8],
            bit: i32,
        ) -> Result<()> {
            let cnt = nqx.div_ceil(2) * nqy.div_ceil(2);
            input.input_nnybble(cnt, scratch)?;
            $bitins(scratch, nqx, nqy, a, n, bit);
            Ok(())
        }

        /// Decode the quadtree-coded bit planes of one quadrant (cfitsio
        /// `qtree_decode`).
        fn $qtree_decode(
            input: &mut HcInput,
            a: &mut [$t],
            a_off: usize,
            n: usize,
            nqx: usize,
            nqy: usize,
            nbitplanes: i32,
        ) -> Result<()> {
            let nqmax = nqx.max(nqy);
            let mut log2n = ((nqmax as f32).ln() / 2.0f32.ln() + 0.5) as i32;
            if nqmax > (1usize << log2n) {
                log2n += 1;
            }
            let nqx2 = nqx.div_ceil(2);
            let nqy2 = nqy.div_ceil(2);
            let mut scratch = vec![0u8; nqx2 * nqy2 + 4];

            let asl = &mut a[a_off..];

            let mut bit = nbitplanes - 1;
            while bit >= 0 {
                let b = input.input_nybble()?;
                if b == 0 {
                    $read_bdirect(input, asl, n, nqx, nqy, &mut scratch, bit)?;
                } else if b != 0xf {
                    return Err(Error::CompressionError(
                        "qtree_decode: bad format code".into(),
                    ));
                } else {
                    scratch[0] = input.input_huffman()? as u8;
                    let mut nx = 1usize;
                    let mut ny = 1usize;
                    let mut nfx = nqx;
                    let mut nfy = nqy;
                    let mut c = 1usize << log2n;
                    let mut k = 1;
                    while k < log2n {
                        c >>= 1;
                        nx <<= 1;
                        ny <<= 1;
                        if nfx <= c {
                            nx -= 1;
                        } else {
                            nfx -= c;
                        }
                        if nfy <= c {
                            ny -= 1;
                        } else {
                            nfy -= c;
                        }
                        qtree_expand(input, &mut scratch, nx, ny)?;
                        k += 1;
                    }
                    $bitins(&scratch, nqx, nqy, asl, n, bit);
                }
                bit -= 1;
            }
            Ok(())
        }

        /// Decode the four quadrants into coefficient array `a` (cfitsio
        /// `dodecode`).
        fn $dodecode(
            input: &mut HcInput,
            a: &mut [$t],
            nx: usize,
            ny: usize,
            nbitplanes: [u8; 3],
        ) -> Result<()> {
            let nel = nx * ny;
            let nx2 = nx.div_ceil(2);
            let ny2 = ny.div_ceil(2);
            for v in a.iter_mut().take(nel) {
                *v = 0 as $t;
            }
            input.start_inputing_bits();
            $qtree_decode(input, a, 0, ny, nx2, ny2, nbitplanes[0] as i32)?;
            $qtree_decode(input, a, ny2, ny, nx2, ny / 2, nbitplanes[1] as i32)?;
            $qtree_decode(input, a, ny * nx2, ny, nx / 2, ny2, nbitplanes[1] as i32)?;
            $qtree_decode(
                input,
                a,
                ny * nx2 + ny2,
                ny,
                nx / 2,
                ny / 2,
                nbitplanes[2] as i32,
            )?;
            if input.input_nybble()? != 0 {
                return Err(Error::CompressionError(
                    "dodecode: bad bit plane values (missing EOF)".into(),
                ));
            }
            // sign bits
            input.start_inputing_bits();
            for v in a.iter_mut().take(nel) {
                if *v != 0 as $t && input.input_bit()? != 0 {
                    *v = (0 as $t).wrapping_sub(*v);
                }
            }
            Ok(())
        }

        fn $undigitize(a: &mut [$t], nel: usize, scale: i32) {
            if scale <= 1 {
                return;
            }
            let s = scale as $t;
            for v in a.iter_mut().take(nel) {
                *v = (*v).wrapping_mul(s);
            }
        }

        /// Inverse H-transform (cfitsio `hinv`). `smooth` is unsupported (the
        /// fixtures use SMOOTH=0); a non-zero value is rejected by the caller.
        fn $hinv(a: &mut [$t], nx: usize, ny: usize) {
            let nmax = nx.max(ny);
            let mut log2n = ((nmax as f32).ln() / 2.0f32.ln() + 0.5) as i32;
            if nmax > (1usize << log2n) {
                log2n += 1;
            }
            let nmax_i = nmax;
            let mut tmp = vec![0 as $t; nmax_i.div_ceil(2) + 1];

            let mut shift: i32 = 1;
            let mut bit0: $t = (1 as $t) << (log2n - 1);
            let mut bit1: $t = bit0 << 1;
            let mut bit2: $t = bit0 << 2;
            let mut mask0: $t = (0 as $t).wrapping_sub(bit0);
            let mut mask1: $t = mask0 << 1;
            let mask2: $t = mask0 << 2;
            let mut prnd0: $t = bit0 >> 1;
            let mut prnd1: $t = bit1 >> 1;
            let prnd2: $t = bit2 >> 1;
            let mut nrnd0: $t = prnd0 - 1;
            let mut nrnd1: $t = prnd1 - 1;
            let nrnd2: $t = prnd2 - 1;

            // round h0 to multiple of bit2
            a[0] = (a[0].wrapping_add(if a[0] >= 0 as $t { prnd2 } else { nrnd2 })) & mask2;

            let ny_i = ny as isize;
            let mut nxtop = 1usize;
            let mut nytop = 1usize;
            let mut nxf = nx;
            let mut nyf = ny;
            let mut c = 1usize << log2n;
            let mut k = log2n - 1;
            while k >= 0 {
                c >>= 1;
                nxtop <<= 1;
                nytop <<= 1;
                if nxf <= c {
                    nxtop -= 1;
                } else {
                    nxf -= c;
                }
                if nyf <= c {
                    nytop -= 1;
                } else {
                    nyf -= c;
                }
                if k == 0 {
                    nrnd0 = 0 as $t;
                    shift = 2;
                }
                // unshuffle in each dimension
                for i in 0..nxtop {
                    $unshuffle(a, ny * i, nytop, 1, &mut tmp);
                }
                for j in 0..nytop {
                    $unshuffle(a, j, nxtop, ny, &mut tmp);
                }
                let oddx = nxtop % 2;
                let oddy = nytop % 2;
                let mut i = 0usize;
                while i + oddx < nxtop {
                    // i steps by 2 over 0..nxtop-oddx
                    let mut s00 = (ny * i) as isize;
                    let mut s10 = s00 + ny_i;
                    let mut j = 0usize;
                    while j + oddy < nytop {
                        let mut h0 = a[s00 as usize];
                        let mut hx = a[s10 as usize];
                        let mut hy = a[(s00 + 1) as usize];
                        let mut hc = a[(s10 + 1) as usize];
                        hx = (hx.wrapping_add(if hx >= 0 as $t { prnd1 } else { nrnd1 })) & mask1;
                        hy = (hy.wrapping_add(if hy >= 0 as $t { prnd1 } else { nrnd1 })) & mask1;
                        hc = (hc.wrapping_add(if hc >= 0 as $t { prnd0 } else { nrnd0 })) & mask0;
                        let lowbit0 = hc & bit0;
                        hx = if hx >= 0 as $t {
                            hx.wrapping_sub(lowbit0)
                        } else {
                            hx.wrapping_add(lowbit0)
                        };
                        hy = if hy >= 0 as $t {
                            hy.wrapping_sub(lowbit0)
                        } else {
                            hy.wrapping_add(lowbit0)
                        };
                        let lowbit1 = (hc ^ hx ^ hy) & bit1;
                        h0 = if h0 >= 0 as $t {
                            h0.wrapping_add(lowbit0).wrapping_sub(lowbit1)
                        } else {
                            h0.wrapping_add(if lowbit0 == 0 as $t {
                                lowbit1
                            } else {
                                lowbit0.wrapping_sub(lowbit1)
                            })
                        };
                        a[(s10 + 1) as usize] =
                            (h0.wrapping_add(hx).wrapping_add(hy).wrapping_add(hc)) >> shift;
                        a[s10 as usize] =
                            (h0.wrapping_add(hx).wrapping_sub(hy).wrapping_sub(hc)) >> shift;
                        a[(s00 + 1) as usize] =
                            (h0.wrapping_sub(hx).wrapping_add(hy).wrapping_sub(hc)) >> shift;
                        a[s00 as usize] =
                            (h0.wrapping_sub(hx).wrapping_sub(hy).wrapping_add(hc)) >> shift;
                        s00 += 2;
                        s10 += 2;
                        j += 2;
                    }
                    if oddy != 0 {
                        let mut h0 = a[s00 as usize];
                        let mut hx = a[s10 as usize];
                        hx = (hx.wrapping_add(if hx >= 0 as $t { prnd1 } else { nrnd1 })) & mask1;
                        let lowbit1 = hx & bit1;
                        h0 = if h0 >= 0 as $t {
                            h0.wrapping_sub(lowbit1)
                        } else {
                            h0.wrapping_add(lowbit1)
                        };
                        a[s10 as usize] = (h0.wrapping_add(hx)) >> shift;
                        a[s00 as usize] = (h0.wrapping_sub(hx)) >> shift;
                    }
                    i += 2;
                }
                if oddx != 0 {
                    let mut s00 = (ny * i) as isize;
                    let mut j = 0usize;
                    while j + oddy < nytop {
                        let mut h0 = a[s00 as usize];
                        let mut hy = a[(s00 + 1) as usize];
                        hy = (hy.wrapping_add(if hy >= 0 as $t { prnd1 } else { nrnd1 })) & mask1;
                        let lowbit1 = hy & bit1;
                        h0 = if h0 >= 0 as $t {
                            h0.wrapping_sub(lowbit1)
                        } else {
                            h0.wrapping_add(lowbit1)
                        };
                        a[(s00 + 1) as usize] = (h0.wrapping_add(hy)) >> shift;
                        a[s00 as usize] = (h0.wrapping_sub(hy)) >> shift;
                        s00 += 2;
                        j += 2;
                    }
                    if oddy != 0 {
                        let h0 = a[s00 as usize];
                        a[s00 as usize] = h0 >> shift;
                    }
                }
                // divide masks/rounding by 2
                bit2 = bit1;
                bit1 = bit0;
                bit0 >>= 1;
                mask1 = mask0;
                mask0 >>= 1;
                prnd1 = prnd0;
                prnd0 >>= 1;
                nrnd1 = nrnd0;
                nrnd0 = prnd0 - 1;
                k -= 1;
            }
            let _ = (bit2, mask1, prnd1, nrnd1);
        }

        /// Full HCOMPRESS decode for one quadrant-int width. Returns the pixel
        /// array (axis-1 fastest) along with `(nx_slow, ny_fast)`.
        fn $name(input: &mut HcInput, smooth: i32) -> Result<(Vec<$t>, usize, usize)> {
            // magic code
            let magic = input.qread(2)?;
            if magic != [0xDDu8, 0x99u8] {
                return Err(Error::CompressionError(
                    "HCOMPRESS: bad magic code".into(),
                ));
            }
            let nx = input.readint()? as usize; // slow axis
            let ny = input.readint()? as usize; // fast axis
            let scale = input.readint()?;
            let nel = nx.checked_mul(ny).ok_or_else(|| {
                Error::CompressionError("HCOMPRESS: dimension overflow".into())
            })?;
            let sumall = input.readlonglong()?;
            let nbp = input.qread(3)?;
            let nbitplanes = [nbp[0], nbp[1], nbp[2]];

            let mut a = vec![0 as $t; nel.max(1)];
            $dodecode(input, &mut a, nx, ny, nbitplanes)?;
            // put sum of all pixels back into pixel 0
            a[0] = sumall as $t;

            if smooth != 0 {
                return Err(Error::UnsupportedCompression(
                    "HCOMPRESS_1 SMOOTH != 0 is not supported".into(),
                ));
            }
            $undigitize(&mut a, nel, scale);
            $hinv(&mut a, nx, ny);
            Ok((a, nx, ny))
        }
    };
}

impl_hdecompress!(
    hdecode32,
    i32,
    qtree_bitins32,
    read_bdirect32,
    qtree_decode32,
    dodecode32,
    hinv32,
    undigitize32,
    unshuffle32
);
impl_hdecompress!(
    hdecode64,
    i64,
    qtree_bitins64,
    read_bdirect64,
    qtree_decode64,
    dodecode64,
    hinv64,
    undigitize64,
    unshuffle64
);

/// Decompress one HCOMPRESS_1 tile.
///
/// `wide` selects the 64-bit transform (used for ZBITPIX 32 and quantized
/// -32/-64 floats; cfitsio uses the 32-bit transform only for ZBITPIX 8/16).
/// Returns the decoded integer samples (axis-1 fastest) and the `(nx_slow,
/// ny_fast)` dimensions read from the stream.
pub fn hcompress_decompress(
    src: &[u8],
    wide: bool,
    smooth: i32,
) -> Result<(Vec<i64>, usize, usize)> {
    let mut input = HcInput::new(src);
    if wide {
        let (a, nx, ny) = hdecode64(&mut input, smooth)?;
        Ok((a, nx, ny))
    } else {
        let (a, nx, ny) = hdecode32(&mut input, smooth)?;
        Ok((a.into_iter().map(|v| v as i64).collect(), nx, ny))
    }
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
    fn plio_single_pixel_at_end_of_run() {
        // PLIO line list for a 512-pixel row with one pixel == 1 at index 178
        // (0-based), all else zero. Lifted verbatim from the cfitsio `fpack -p`
        // output for row 137 of EUVEngc4151imgx.fits. Header word [2] = -100 (the
        // PLIO magic, 0xFF9C) forces the long-form length decode.
        let words: Vec<i16> = [0u16, 7, 0xFF9C, 9, 0, 0, 0, 20659, 333]
            .iter()
            .map(|&w| w as i16)
            .collect();
        let out = plio_decompress(&words, 512).unwrap();
        let mut expected = vec![0i32; 512];
        expected[178] = 1; // opcode 5 (run 179, last pixel set to pv=1)
        assert_eq!(out, expected);
    }

    #[test]
    fn plio_empty_list_is_all_zero() {
        // A zero-length line list (lllen == 0) decodes to an all-zero tile.
        let words: Vec<i16> = [0u16, 7, 0xFF9C, 0, 0, 0, 0]
            .iter()
            .map(|&w| w as i16)
            .collect();
        let out = plio_decompress(&words, 16).unwrap();
        assert_eq!(out, vec![0i32; 16]);
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
        for (tile, td) in tile_dims.iter().enumerate() {
            let coords = unravel(tile, &tiles_per_axis);
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
    fn fits_rand_value_matches_cfitsio_checkpoint() {
        // cfitsio validates fits_init_randoms by asserting the final LCG seed is
        // 1043618065 after 10000 iterations. Reproduce that seed and check it, which
        // exercises the exact same recurrence our table generator uses.
        let a = 16807.0f64;
        let m = 2_147_483_647.0f64;
        let mut seed = 1.0f64;
        for _ in 0..N_RANDOM {
            let temp = a * seed;
            seed = temp - m * ((temp / m) as i64 as f64);
        }
        assert_eq!(seed as i64, 1_043_618_065);
        // All table values lie in [0, 1).
        for i in [0usize, 1, 499, 5000, N_RANDOM - 1] {
            let v = fits_rand_value(i);
            assert!((0.0..1.0).contains(&v), "rand[{i}] = {v} out of range");
        }
    }

    #[test]
    fn nodither_float_constant_tile() {
        // A NO_DITHER quantized float tile: value = q*scale + zero (no +0.5, no dither).
        use crate::bintable::{BinColumnType, BinTableBuilder};

        // 2-pixel image, single row tile, quantized ints [10, 20], scale 0.5, zero 3.0.
        let q = [10i32, 20];
        let enc = rice_encode_i32(&q, 2);

        let table = BinTableBuilder::new()
            .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
            .add_column("ZSCALE", BinColumnType::D64(1))
            .add_column("ZZERO", BinColumnType::D64(1))
            .push_row(|r| {
                r.write_var_p(enc.len() as i32, |heap| heap.extend_from_slice(&enc));
                r.write_f64(0.5);
                r.write_f64(3.0);
            })
            .build();

        let mut h = Header::new();
        h.set("ZIMAGE", HeaderValue::Logical(true), None);
        h.set("ZCMPTYPE", HeaderValue::String("RICE_1".into()), None);
        h.set("ZBITPIX", HeaderValue::Integer(-32), None);
        h.set("ZNAXIS", HeaderValue::Integer(1), None);
        h.set("ZNAXIS1", HeaderValue::Integer(2), None);
        h.set("ZQUANTIZ", HeaderValue::String("NO_DITHER".into()), None);
        h.set("ZNAME1", HeaderValue::String("BYTEPIX".into()), None);
        h.set("ZVAL1", HeaderValue::Integer(4), None);

        let cimg = CompressedImage::from_bintable(&h, &table).unwrap();
        let img = cimg.decompress().unwrap();
        match img.pixels {
            crate::image_data::PixelData::F32(v) => {
                assert_eq!(v, vec![10.0 * 0.5 + 3.0, 20.0 * 0.5 + 3.0]);
            }
            other => panic!("expected F32, got {other:?}"),
        }
    }

    #[test]
    fn dither2_preserves_zero_and_blank_maps_to_nan() {
        use crate::bintable::{BinColumnType, BinTableBuilder};

        // 3-pixel tile: [ZERO_VALUE, NULL/blank, 5]. DITHER_2 => [0.0, NaN, dithered].
        // Use NOCOMPRESS so the quantized ints are stored verbatim (RICE encoding of
        // these extreme values is awkward and unrelated to what we're testing here).
        let q = [ZERO_VALUE as i32, NULL_VALUE as i32, 5];
        let mut enc = Vec::new();
        for v in q {
            enc.extend_from_slice(&v.to_be_bytes());
        }

        let table = BinTableBuilder::new()
            .add_column("COMPRESSED_DATA", BinColumnType::VarP('B'))
            .add_column("ZSCALE", BinColumnType::D64(1))
            .add_column("ZZERO", BinColumnType::D64(1))
            .push_row(|r| {
                r.write_var_p(enc.len() as i32, |heap| heap.extend_from_slice(&enc));
                r.write_f64(2.0);
                r.write_f64(1.0);
            })
            .build();

        let mut h = Header::new();
        h.set("ZIMAGE", HeaderValue::Logical(true), None);
        h.set("ZCMPTYPE", HeaderValue::String("NOCOMPRESS".into()), None);
        h.set("ZBITPIX", HeaderValue::Integer(-32), None);
        h.set("ZNAXIS", HeaderValue::Integer(1), None);
        h.set("ZNAXIS1", HeaderValue::Integer(3), None);
        h.set("ZQUANTIZ", HeaderValue::String("SUBTRACTIVE_DITHER_2".into()), None);
        h.set("ZDITHER0", HeaderValue::Integer(5), None);
        h.set("ZBLANK", HeaderValue::Integer(NULL_VALUE), None);
        h.set("ZNAME1", HeaderValue::String("BYTEPIX".into()), None);
        h.set("ZVAL1", HeaderValue::Integer(4), None);

        let cimg = CompressedImage::from_bintable(&h, &table).unwrap();
        let img = cimg.decompress().unwrap();
        match img.pixels {
            crate::image_data::PixelData::F32(v) => {
                assert_eq!(v[0], 0.0); // ZERO_VALUE preserved exactly
                assert!(v[1].is_nan()); // blank -> NaN
                assert!(v[2].is_finite() && v[2] != 0.0);
            }
            other => panic!("expected F32, got {other:?}"),
        }
    }

    // --- RICE_1 encode (the new write path) ----------------------------------

    /// The RICE encoder's output must decode via the existing decoder, for a range of
    /// difference magnitudes / signs / block boundaries, across all three BYTEPIX.
    #[test]
    fn rice_encode_decode_roundtrip_i32() {
        let cases: Vec<Vec<i32>> = vec![
            vec![0],
            vec![42, 42, 42, 42],                  // zero-diff block
            vec![1000, 1003, 1001, 1005, 1002],    // small diffs
            vec![i32::MIN, 0, i32::MAX, -1, 1],     // extreme diffs (verbatim path)
            (0..100).map(|i| (i * i) % 7 - 3).collect(), // > one block
            (0..200).map(|i: i32| (i - 100) * 17).collect(),
        ];
        for vals in cases {
            let enc = rice_compress_i32(&vals, 32);
            let dec = rice_decompress_i32(&enc, vals.len(), 32).unwrap();
            assert_eq!(dec, vals, "i32 roundtrip failed for {vals:?}");
        }
    }

    #[test]
    fn rice_encode_decode_roundtrip_i16() {
        let vals: Vec<i16> = (0..300).map(|i| ((i * 31) % 251 - 120) as i16).collect();
        let enc = rice_compress_i16(&vals, 32);
        let dec = rice_decompress_i16(&enc, vals.len(), 32).unwrap();
        assert_eq!(dec, vals);
        // extremes
        let ex: Vec<i16> = vec![i16::MIN, i16::MAX, 0, -1, 1, i16::MIN];
        let enc = rice_compress_i16(&ex, 32);
        assert_eq!(rice_decompress_i16(&enc, ex.len(), 32).unwrap(), ex);
    }

    #[test]
    fn rice_encode_decode_roundtrip_i8() {
        let vals: Vec<u8> = (0..400).map(|i| ((i * 7) % 256) as u8).collect();
        let enc = rice_compress_i8(&vals, 32);
        let dec = rice_decompress_i8(&enc, vals.len(), 32).unwrap();
        assert_eq!(dec, vals);
    }

    /// A non-default block size must round-trip too (the block boundary moves).
    #[test]
    fn rice_encode_blocksize_16() {
        let vals: Vec<i32> = (0..70).map(|i| (i * 3) % 11 - 5).collect();
        let enc = rice_compress_i32(&vals, 16);
        let dec = rice_decompress_i32(&enc, vals.len(), 16).unwrap();
        assert_eq!(dec, vals);
    }

    /// `ImageData::compress` (RICE_1) → `as_compressed_image().decompress()` must be a
    /// byte-exact identity for an integer image, including a 2-D square tiling with
    /// edge-truncated tiles.
    #[test]
    fn compress_image_rice_int_internal_roundtrip() {
        use crate::image_data::{ImageData, PixelData};

        // 13x11 I16 image, value = mix of trend + noise; square 5x5 tiles (edges ragged).
        let w = 13usize;
        let h = 11usize;
        let pixels: Vec<i16> = (0..(w * h))
            .map(|i| (((i * 37) % 521) as i16) - 200)
            .collect();
        let img = ImageData::new(vec![w, h], PixelData::I16(pixels));

        for tile in [None, Some(vec![w, 1]), Some(vec![5, 5]), Some(vec![13, 11])] {
            let opts = CompressOptions {
                algorithm: CompressionType::Rice1,
                tile,
                ..Default::default()
            };
            let hdu = img.compress(&opts).unwrap();
            assert!(matches!(hdu.data, crate::hdu::HduData::BinTable(_)));
            let back = hdu
                .as_compressed_image()
                .expect("ZIMAGE detected")
                .decompress()
                .unwrap();
            assert_eq!(back.axes, img.axes);
            assert_eq!(back.pixels.to_bytes(), img.pixels.to_bytes());
        }
    }

    #[test]
    fn compress_image_rice_i32_and_u8_roundtrip() {
        use crate::image_data::{ImageData, PixelData};

        let i32img = ImageData::new(
            vec![20, 3],
            PixelData::I32((0..60).map(|i| (i - 30) * 12345).collect()),
        );
        let hdu = i32img.compress(&CompressOptions::default()).unwrap();
        let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
        assert_eq!(back.pixels.to_bytes(), i32img.pixels.to_bytes());

        let u8img = ImageData::new(
            vec![17, 4],
            PixelData::U8((0..68).map(|i| (i * 5 % 256) as u8).collect()),
        );
        let hdu = u8img.compress(&CompressOptions::default()).unwrap();
        let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
        assert_eq!(back.pixels.to_bytes(), u8img.pixels.to_bytes());
    }

    #[test]
    fn compress_rejects_unsupported_algorithms() {
        use crate::image_data::{ImageData, PixelData};
        let img = ImageData::new(vec![4], PixelData::I16(vec![1, 2, 3, 4]));
        for alg in [
            CompressionType::Plio1,
            CompressionType::Hcompress1,
            CompressionType::NoCompress,
        ] {
            let opts = CompressOptions {
                algorithm: alg,
                ..Default::default()
            };
            assert!(matches!(
                img.compress(&opts),
                Err(Error::UnsupportedCompression(_))
            ));
        }
        // 64-bit integer images are rejected.
        let i64img = ImageData::new(vec![3], PixelData::I64(vec![1, 2, 3]));
        assert!(i64img.compress(&CompressOptions::default()).is_err());
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn compress_image_gzip_int_roundtrip() {
        use crate::image_data::{ImageData, PixelData};
        let img = ImageData::new(
            vec![15, 6],
            PixelData::I16((0..90).map(|i| (i * 13 % 400 - 200) as i16).collect()),
        );
        for alg in [CompressionType::Gzip1, CompressionType::Gzip2] {
            let opts = CompressOptions {
                algorithm: alg,
                ..Default::default()
            };
            let hdu = img.compress(&opts).unwrap();
            let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
            assert_eq!(
                back.pixels.to_bytes(),
                img.pixels.to_bytes(),
                "{alg:?} int roundtrip"
            );
        }
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn compress_image_gzip_lossless_float_roundtrip() {
        use crate::image_data::{ImageData, PixelData};
        let pixels: Vec<f32> = (0..96).map(|i| (i as f32) * 0.125 - 6.0).collect();
        let img = ImageData::new(vec![12, 8], PixelData::F32(pixels));
        // Lossless float storage is GZIP_1 only (GZIP_2 float shuffle is unsupported).
        let opts = CompressOptions {
            algorithm: CompressionType::Gzip1,
            quantize: None, // lossless raw floats
            ..Default::default()
        };
        let hdu = img.compress(&opts).unwrap();
        let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
        assert_eq!(
            back.pixels.to_bytes(),
            img.pixels.to_bytes(),
            "GZIP_1 lossless float roundtrip"
        );
        // GZIP_2 lossless float is rejected.
        let opts2 = CompressOptions {
            algorithm: CompressionType::Gzip2,
            quantize: None,
            ..Default::default()
        };
        assert!(matches!(
            img.compress(&opts2),
            Err(Error::UnsupportedCompression(_))
        ));
    }

    /// Lossy quantized-float RICE encode must round-trip within ~scale tolerance, and
    /// the reconstructed values must match our own decoder's `unquantize` exactly.
    #[test]
    fn compress_image_rice_float_quantize_roundtrip() {
        use crate::image_data::{ImageData, PixelData};
        let pixels: Vec<f32> = (0..256)
            .map(|i| ((i as f32) * 0.017).sin() * 100.0 + 50.0)
            .collect();
        let img = ImageData::new(vec![16, 16], PixelData::F32(pixels.clone()));
        let opts = CompressOptions {
            algorithm: CompressionType::Rice1,
            tile: Some(vec![16, 4]),
            quantize: Some(4.0),
            dither: Quantize::SubtractiveDither1,
            dither_seed: Some(5),
            ..Default::default()
        };
        let hdu = img.compress(&opts).unwrap();
        let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
        let recon = match &back.pixels {
            PixelData::F32(v) => v.clone(),
            other => panic!("expected F32, got {other:?}"),
        };
        assert_eq!(recon.len(), pixels.len());
        // The dither RNG indexing the encoder used is the same the decoder inverts, so
        // the error is bounded by the per-tile quantization step. Check a generous
        // tolerance relative to the data range.
        let range = 100.0f32; // amplitude
        let tol = range / 1000.0;
        for (i, (&o, &r)) in pixels.iter().zip(recon.iter()).enumerate() {
            assert!(
                (o - r).abs() <= tol,
                "pixel {i}: orig {o} recon {r} exceeds tol {tol}"
            );
        }
    }
}
