# FITS Tiled Image Compression — Implementation Plan

Status: planning + initial scaffold. This document describes how to add support for
the FITS *Tiled Image Compression* convention (the BINTABLE-based compressed-image
format) to `fits4`, currently the main feature gap versus `fitsrs`.

References:
- Tiled Image Convention for Storing Compressed Images in FITS Binary Tables,
  v2.3 — https://fits.gsfc.nasa.gov/registry/tilecompression/tilecompression2.3.pdf
  (also arXiv:1201.1336)
- HEASARC FITSIO compression docs — https://heasarc.gsfc.nasa.gov/docs/software/fitsio/compression.html
- cfitsio `ricecomp.c` — https://github.com/HEASARC/cfitsio/blob/develop/ricecomp.c
- astropy compressed codecs / `_tiled_compression.py`

## 1. What the convention is

A tile-compressed image is physically a **BINTABLE extension** (`XTENSION = 'BINTABLE'`)
but is logically an image. The original image is divided into a rectangular grid of
*tiles*; each tile is compressed independently and stored as a byte stream in a
variable-length-array (VLA) cell of a binary-table column (typically `COMPRESSED_DATA`,
a `1PB` or `1QB` column). One table row == one tile, in row-major (FORTRAN/NAXIS1-fastest)
tile order.

### Driving header keywords

Mandatory:
- `ZIMAGE = T` — marks this BINTABLE as a compressed image.
- `ZCMPTYPE` — algorithm name: `RICE_1`, `GZIP_1`, `GZIP_2`, `HCOMPRESS_1`, `PLIO_1`
  (also `RICE_ONE` alias, `NOCOMPRESS`).
- `ZBITPIX` — BITPIX of the *original* (uncompressed) image (8/16/32/64/-32/-64).
- `ZNAXIS` — NAXIS of the original image.
- `ZNAXISn` — NAXISn of the original image (n = 1..ZNAXIS).

Optional / conditional:
- `ZTILEn` — tile size along axis n. Default = row-by-row: `ZTILE1 = ZNAXIS1`, all
  others = 1.
- `ZNAMEi` / `ZVALi` — algorithm parameter name/value pairs:
  - RICE_1: `BLOCKSIZE` (default 32), `BYTEPIX` (default 4 — bytes per pixel of the
    *integers* fed to Rice; note this is the post-quantization integer size).
  - HCOMPRESS_1: `SCALE`, `SMOOTH`.
- `ZQUANTIZ` — float quantization method: `NO_DITHER`, `SUBTRACTIVE_DITHER_1`,
  `SUBTRACTIVE_DITHER_2`. Absent ⇒ no dithering (older files quantized without dither).
- `ZDITHER0` — integer dither seed (1..10000) for subtractive dithering.
- `ZBLANK` — integer value used in the integer stream to flag undefined/NaN pixels
  (when given as a keyword it is constant for all tiles; can also be a per-tile column).
- `ZSIMPLE`, `ZEXTEND`, `ZBLOCKED`, `ZPCOUNT`, `ZGCOUNT`, `ZHECKSUM`, `ZDATASUM` —
  preserve the original primary/extension header context (mostly carried through;
  `ZSIMPLE`/`ZTENSION` distinguish original primary vs extension).
- `ZMASKCMP` — algorithm used to compress the null-pixel mask (rarely needed for read).
- The original image's other header keywords are copied verbatim into the BINTABLE
  header (so WCS etc. survive); on decompression we reconstruct an image header from
  the `Z*` keywords plus the carried-through cards.

### Binary-table columns

- `COMPRESSED_DATA` (required) — VLA (`1PB`/`1QB`) holding the compressed tile bytes.
  When a tile failed to compress/quantize, its descriptor count is 0 and the data
  lives in the fallback column instead.
- `GZIP_COMPRESSED_DATA` (optional) — VLA fallback: a gzip-compressed copy of the
  *un-quantized* tile, used when lossy quantization was rejected for that tile.
- `UNCOMPRESSED_DATA` (optional, older variant) — raw uncompressed tile fallback.
- `ZSCALE`, `ZZERO` (optional, columns of D) — per-tile linear scale/zero for float
  quantization: `float_value = ZZERO + ZSCALE * integer_value`. If absent, use the
  `ZSCALE`/`ZZERO` *keywords*, or 1.0/0.0.
- `ZBLANK` (optional column of J) — per-tile blank sentinel.

### Reassembly

For each row r (tile index): decode the tile's compressed bytes with the ZCMPTYPE
algorithm into a flat array of `ZBITPIX`-typed integers (or, for the gzip fallback,
the original dtype). For float originals, apply inverse quantization
(`ZZERO + ZSCALE * q`, plus dither reversal). Then scatter the tile's pixels into the
correct sub-rectangle of the full image buffer using the tile grid geometry derived
from `ZNAXISn` / `ZTILEn`. Tiles iterate with axis 1 fastest.

## 2. Float quantization & dithering

Float (`ZBITPIX < 0`) images are usually stored as quantized 32-bit integers
(`BYTEPIX` typically 4) plus per-tile `ZSCALE`/`ZZERO`. Decode:

1. Decompress the integer stream `q[i]`.
2. If `q[i] == ZBLANK` ⇒ output `NaN`.
3. Otherwise `value = ZZERO + ZSCALE * q[i]` — with dither correction:
   - `NO_DITHER` / absent: `value = ZZERO + ZSCALE * q[i]`.
   - `SUBTRACTIVE_DITHER_1`: a per-pixel uniform dither `r ∈ [0,1)` was *added* to
     `q` before rounding at compress time, so decode is
     `value = ZZERO + ZSCALE * (q[i] - r + 0.5)` (cfitsio uses `(q - dither + 0.5)`).
   - `SUBTRACTIVE_DITHER_2`: same as DITHER_1 except exact-zero pixels are preserved
     (encoded with the special integer `-2147483647` ⇒ decode to exactly 0.0).
   The dither value `r` comes from cfitsio's fixed pseudo-random sequence
   (`fits_rand_value`, a 10000-element table generated by a specific LCG), indexed by
   `(tile_row + ZDITHER0 - 1)` to pick a starting offset, advancing per pixel and
   wrapping at 10000. This table/sequence must be reproduced exactly to round-trip
   dithered files — see cfitsio `quantize.c` `fits_init_randoms`.

This is the trickiest correctness item; for the first read-only pass we can support
`NO_DITHER` and the absent case, and treat dithered files as a follow-up (decoding
without dither reversal yields values off by up to ~0.5 ZSCALE — acceptable only as a
clearly-flagged approximation, so better to error until the RNG table is implemented).

## 3. Module design

New module `src/tile_compress.rs`, gated behind nothing (pure Rust; gzip needs an
inflate implementation — see open questions). Public surface:

```rust
pub enum CompressionType { Rice1, Gzip1, Gzip2, Hcompress1, Plio1, NoCompress }
pub enum Quantize { None, SubtractiveDither1, SubtractiveDither2 }

pub struct TileGeometry { znaxis: Vec<usize>, ztile: Vec<usize>, /* grids */ }

pub struct CompressedImage<'a> {
    pub ctype: CompressionType,
    pub zbitpix: i64,
    pub geometry: TileGeometry,
    pub quantize: Quantize,
    pub zdither0: i64,
    pub blank: Option<i64>,
    table: &'a BinTable,           // borrows the parsed BINTABLE
    // resolved column indices: compressed_data, gzip_fallback, zscale, zzero, zblank
}

impl CompressedImage<'_> {
    pub fn detect(header: &Header) -> bool;          // ZIMAGE == T
    pub fn from_bintable(header: &Header, table: &BinTable) -> Result<Self>;
    pub fn decompress(&self) -> Result<ImageData>;   // full image
    fn decompress_tile(&self, row: usize) -> Result<Vec<i64-ish>>;
}

// algorithm primitives (pure functions, unit-testable in isolation):
pub fn rice_decompress_i32(src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<i32>>;
pub fn rice_decompress_i16(src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<i16>>;
pub fn rice_decompress_i8 (src: &[u8], nvals: usize, blocksize: usize) -> Result<Vec<u8>>;
pub fn gzip_inflate(src: &[u8]) -> Result<Vec<u8>>;          // GZIP_1
// GZIP_2 = byte-shuffled gzip (de-shuffle after inflate)
pub fn plio_decompress(src: &[u8], nvals: usize) -> Result<Vec<i32>>;
pub fn hcompress_decompress(...) -> Result<Vec<i32>>;        // last
```

### Read-pipeline hookup

In `hdu.rs::read_from`, after the BINTABLE branch parses a `BinTable`, check
`CompressedImage::detect(&header)`. Options for surfacing it (see open questions):

- **Option A (recommended):** add a new `HduData::Image(ImageData)` produced eagerly —
  i.e. when `ZIMAGE=T`, decompress and store the result as `HduData::Image`, so callers
  treat it like any other image. Keep the original BINTABLE accessible only if needed.
  Downside: loses the raw table; re-writing won't round-trip the compressed bytes.
- **Option B:** add a new variant `HduData::CompressedImage(BinTable)` (or keep the
  `BinTable` and add a method `Hdu::as_image()` that lazily decompresses). Preserves
  round-trip and is lazy, at the cost of a new public enum variant.
- **Option C:** keep `HduData::BinTable`, and add `BinTable::as_compressed_image()`
  / `Hdu::decompressed_image()` returning `Option<Result<ImageData>>`. Smallest API
  change, fully backward compatible, lazy. Likely the best first step.

The plan adopts **Option C** for the initial implementation (no enum change, no
behavior change for existing users), with a possible convenience accessor on `FitsFile`.

### Error handling

Add `Error` variants: `UnsupportedCompression(String)`, `CompressionError(String)`.

## 4. Staged implementation

**Phase 0 (this pass — scaffold):** module skeleton, enums, `detect`, `from_bintable`
header parsing, function stubs with `todo!()`, RICE_1 32/16/8-bit decoder + unit test.
Wire `mod tile_compress;` into `lib.rs`. Build must stay green.

**Phase 1 — RICE_1 read path (integers):** finish `decompress()` for integer
`ZBITPIX` with `NO_DITHER`; tile scatter/reassembly; `ZBLANK` handling. Validate
against an `fpack`-generated test file.

**Phase 2 — GZIP_1 / GZIP_2:** integer + float read. Needs a DEFLATE decoder
(see open questions re: dependency). GZIP_2 adds byte de-shuffle.

**Phase 3 — float quantization + dithering:** ZSCALE/ZZERO inverse, NaN/ZBLANK,
and the cfitsio random-number table for SUBTRACTIVE_DITHER_1/2.

**Phase 4 — PLIO_1:** IRAF pixel-list run-length decoder (fully specified in the
convention doc, simplest of the lossless coders; mainly used for integer masks).

**Phase 5 — HCOMPRESS_1:** wavelet-ish H-transform decode + quantization. Most
complex; lowest priority.

**Phase 6 — writing/compression:** mirror each decoder with an encoder, build the
BINTABLE + `Z*` keywords, VLA heap packing. RICE_1 + GZIP first; HCOMPRESS optional.
Hook into `FitsFile` builder (e.g. `add_compressed_image(img, CompressionType::Rice1)`).

## 5. Test files

None of the existing `samp/*.fits` files are tile-compressed (verified: they are
plain IMAGE / TABLE / BINTABLE HDUs). We need test fixtures:

- Generate with cfitsio `fpack`:
  - `fpack -r image.fits`  → RICE_1
  - `fpack -g image.fits`  → GZIP_1
  - `fpack -g2 image.fits` → GZIP_2
  - `fpack -h image.fits`  → HCOMPRESS_1
  - `fpack -p image.fits`  → PLIO_1 (integer images / masks)
  Cover at least: I16 image, I32 image, F32 image (quantized, both dithered and
  `-q 0`/no-dither), and a small 3-D cube to exercise multi-axis tiling.
- Cross-check decompressed output against the original uncompressed FITS (byte-exact
  for lossless integer codecs; within quantization tolerance for float).
- Store fixtures alongside the existing samples (GCS bucket `fits4_samples`, per
  MEMORY.md) so CI can fetch them; add a `tests/tile_compress.rs` integration test.

## 6. Open questions for the human

1. **API shape (most important):** Option C (lazy `Hdu::decompressed_image()` accessor,
   no enum change) vs Option B (new `HduData::CompressedImage` variant) vs Option A
   (eager — present compressed images transparently as `HduData::Image`). Transparent
   presentation is friendliest for consumers like `fitsview`/`startracker` but breaks
   round-trip writing. Recommendation: start with C, consider a transparent
   `FitsFile::images()` helper later.
2. **GZIP dependency:** GZIP_1/GZIP_2 require DEFLATE. The crate currently has *zero*
   dependencies for core. Choices: (a) put gzip behind an optional `flate2`/`miniz_oxide`
   feature; (b) write a small pure-Rust inflate (more work, keeps zero-dep promise);
   (c) ship RICE_1/PLIO only in core and gate gzip/hcompress behind a feature.
   Recommendation: optional `flate` feature using `miniz_oxide` (pure Rust, no C).
3. **Dithering scope:** is byte-exact round-trip of `SUBTRACTIVE_DITHER_1/2` float
   files required, or is read-only (and erroring on dithered files until Phase 3) okay
   for the first release? This decides how soon we must port the cfitsio RNG table.
4. **Surfacing scaled vs raw:** for float quantized images, do we return the
   reconstructed `PixelData::F32`/`F64` directly (applying ZSCALE/ZZERO + dither), or
   also expose the raw quantized integers + scale like the existing BSCALE/BZERO split
   on `ImageData`? Recommend returning reconstructed floats (matches user expectation);
   raw access can be a later addition.
5. **Writing priority:** is the immediate goal read parity with `fitsrs`, or also
   write/compress? Read-only RICE+GZIP closes most of the practical gap.
