# fitskit

**Pure-Rust, zero-dependency reader and writer for [FITS](https://fits.gsfc.nasa.gov/fits_standard.html) v4.0** — the Flexible Image Transport System format used throughout astronomy.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

fitskit reads *and* writes the full FITS v4.0 standard — primary and extension HDUs, images, ASCII tables, and binary tables (including variable-length arrays) — with **no external dependencies** in the default build and **no C toolchain** required. It decodes every tile-compressed image type in the standard (RICE, GZIP, PLIO, HCOMPRESS) and — uniquely among pure-Rust FITS crates — **writes** RICE- and GZIP-compressed images that cfitsio's `funpack` reads back byte-for-byte.

## Features

- **All standard HDU types** — Primary, IMAGE, ASCII `TABLE`, `BINTABLE`
- **Full read *and* write** — round-trip files, bytes, or any `Read`/`Write`
- **All BITPIX types** — 8, 16, 32, 64, -32, -64
- **BSCALE/BZERO** — raw and scaled access, with the unsigned-integer convention
- **Variable-length arrays** — `P` and `Q` heap descriptors in binary tables
- **CHECKSUM/DATASUM** — ones-complement integrity computation and verification
- **Tile-compressed images (read)** — RICE_1, GZIP_1, GZIP_2, PLIO_1, and HCOMPRESS_1, including float quantization with subtractive dithering; decoded bit-exactly against cfitsio's `funpack`
- **Tile-compressed images (write)** — encode RICE_1 and GZIP_1/GZIP_2 (integer and float); output verified byte-exact through `funpack`
- **Zero dependencies by default** — optional `image` and `gzip` features pull in pure-Rust crates only

## Installation

```toml
[dependencies]
fitskit = "0.2"
```

## Usage

### Reading a FITS file

```rust
use fitskit::{FitsFile, HduData, PixelData};

let fits = FitsFile::from_file("image.fits")?;

let primary = fits.primary();
println!("BITPIX = {}", primary.header.get_int("BITPIX").unwrap());

if let HduData::Image(img) = &primary.data {
    println!("{}x{}", img.width().unwrap(), img.height().unwrap());
    if let PixelData::F32(pixels) = &img.pixels {
        println!("first pixel = {}", pixels[0]);
    }
}

for hdu in fits.extensions() {
    match &hdu.data {
        HduData::Image(img) => println!("IMAGE: {:?}", img.axes),
        HduData::BinTable(t) => println!("BINTABLE: {} rows, {} cols", t.nrows, t.columns.len()),
        HduData::AsciiTable(t) => println!("TABLE: {} rows", t.nrows),
        HduData::Empty => {}
    }
}
# Ok::<(), fitskit::Error>(())
```

### Writing a FITS file

```rust
use fitskit::{FitsFile, ImageData, PixelData, HeaderValue};

let pixels: Vec<i16> = (0..10000).map(|i| (i % 1000) as i16).collect();
let img = ImageData::new(vec![100, 100], PixelData::I16(pixels));

let mut fits = FitsFile::with_primary_image(img);
fits.primary_mut().header.set("OBJECT", HeaderValue::String("M31".into()), None);

fits.to_file("output.fits")?;
# Ok::<(), fitskit::Error>(())
```

### Building a binary table

```rust
use fitskit::{FitsFile, Hdu, BinTableBuilder, BinColumnType};

let table = BinTableBuilder::new()
    .add_column("RA", BinColumnType::D64(1))
    .add_column("DEC", BinColumnType::D64(1))
    .add_column("MAG", BinColumnType::E32(1))
    .push_row(|row| {
        row.write_f64(180.0);
        row.write_f64(45.0);
        row.write_f32(12.5);
    })
    .push_row(|row| {
        row.write_f64(90.0);
        row.write_f64(-30.0);
        row.write_f32(8.2);
    })
    .build();

let mut fits = FitsFile::with_empty_primary();
fits.push_extension(Hdu::bintable_extension(table));
# let _ = fits.to_bytes();
```

### Reading a tile-compressed image

Tile-compressed images are stored on disk as a `BINTABLE` extension (the FITS
Tiled Image Compression convention). `Hdu::as_compressed_image` cheaply detects
them; `decompress` reassembles the full image lazily.

```rust
use fitskit::FitsFile;

let fits = FitsFile::from_file("compressed.fits")?;

for hdu in fits.extensions() {
    if let Some(cimg) = hdu.as_compressed_image() {
        println!("compression = {:?}", cimg.compression());
        let image = cimg.decompress()?;
        println!("decompressed to {:?}", image.axes);
    }
}
# Ok::<(), fitskit::Error>(())
```

### Writing a tile-compressed image

`ImageData::compress` produces a compressed-image `BINTABLE` HDU ready to push
into a file. The output is read back byte-for-byte by cfitsio's `funpack`.

```rust
use fitskit::{FitsFile, ImageData, PixelData, CompressOptions};

let pixels: Vec<i16> = (0..10000).map(|i| (i % 1000) as i16).collect();
let img = ImageData::new(vec![100, 100], PixelData::I16(pixels));

// Default options: RICE_1, one tile per row, lossless for integers.
let hdu = img.compress(&CompressOptions::default())?;

let mut fits = FitsFile::with_empty_primary();
fits.push_extension(hdu);
fits.to_file("compressed.fits")?;
# Ok::<(), fitskit::Error>(())
```

> **`gzip` feature:** the `GZIP_1`/`GZIP_2` algorithms (e.g. `CompressOptions { algorithm: CompressionType::Gzip1, .. }` when writing, or decoding a GZIP-compressed tile on read) require the `gzip` feature; `RICE_1` works without it. Build with `--features gzip` or `fitskit = { version = "0.2", features = ["gzip"] }`.

### Converting to/from the `image` crate (feature `image`)

With the `image` feature, `ImageData` converts to and from the [`image`](https://crates.io/crates/image) crate's `DynamicImage` — e.g. to save a FITS image as a PNG, or to ingest a raster as FITS. Build with `fitskit = { version = "0.2", features = ["image"] }`.

```rust
use fitskit::{FitsFile, HduData, ImageData};

let fits = FitsFile::from_file("image.fits")?;
if let HduData::Image(img) = &fits.primary().data {
    // FITS → image crate (BSCALE/BZERO applied; 1.0/0.0 = identity)
    let dynamic = img.to_dynamic_image(1.0, 0.0)?;
    dynamic.save("image.png").unwrap();

    // image crate → FITS, returning (ImageData, bscale, bzero)
    let (restored, _bscale, _bzero) = ImageData::from_dynamic_image(&dynamic)?;
    assert_eq!(restored.axes, img.axes);
}
# Ok::<(), fitskit::Error>(())
```

### Checksums

```rust
use fitskit::{FitsFile, ImageData, PixelData};

let img = ImageData::new(vec![4], PixelData::U8(vec![1, 2, 3, 4]));
let fits = FitsFile::with_primary_image(img);

// Write with CHECKSUM/DATASUM keywords
let bytes = fits.to_bytes_with_checksum()?;

// Verify on read
let fits2 = FitsFile::from_bytes(&bytes)?;
fits2.primary().verify_datasum()?;
# Ok::<(), fitskit::Error>(())
```

### World Coordinate System (WCS)

With the `wcs` feature, parse a two-axis celestial WCS straight from a header and
map between **1-based FITS pixel** coordinates and **world coordinates in
degrees**. The spherical-projection math is backed by the zero-dependency
[`mapproj`](https://crates.io/crates/mapproj) crate.

```rust,no_run
use fitskit::FitsFile;

let fits = FitsFile::from_file("image.fits")?;
let wcs = fits.primary().header.wcs()?;

// Pixel (1-based) -> world (lon, lat) in degrees
let (ra, dec) = wcs.pixel_to_world(256.5, 256.5).unwrap();

// ...and back to pixel
let (x, y) = wcs.world_to_pixel(ra, dec).unwrap();
# Ok::<(), fitskit::Error>(())
```

Supported: the common 2-axis celestial case — `CTYPEn = xxx--CCC` for any
projection `mapproj` implements (`TAN`, `SIN`, `ARC`, `ZEA`, `STG`, `AIT`, `MER`,
`CAR`, `CEA`, `SFL`, `MOL`, conics, …), with the linear transform from `CDi_j` or
`PCi_j` + `CDELTi` (CD takes precedence). SIP distortions and 3+-axis / spectral
WCS are out of scope.

## Core types

A FITS file is an ordered sequence of **Header-Data Units (HDUs)**. fitskit mirrors
that structure directly:

| Type | Role |
|------|------|
| `FitsFile` | Top-level container — an ordered `Vec<Hdu>`. Reads/writes files, byte slices, or any `Read`/`Write`; builder methods (`with_primary_image`, `with_empty_primary`, `push_extension`, `to_file`, `to_bytes`, `to_bytes_with_checksum`). |
| `Hdu` | One Header-Data Unit: a `Header` plus an `HduData` payload. `as_compressed_image()` exposes a tile-compressed image stored in a `BINTABLE`. |
| `HduData` | The payload enum: `Empty`, `Image(ImageData)`, `AsciiTable(AsciiTable)`, `BinTable(BinTable)`. |
| `Header` | Ordered list of `Keyword`s with typed accessors (`get_int`, `get_float`, `get_string`, `get_bool`) and `set`. |
| `Keyword` / `HeaderValue` | An 80-byte header card and its typed value (with `CONTINUE` long-string handling). |
| `ImageData` / `PixelData` | Image `axes` plus a typed pixel buffer (`U8`/`I16`/`I32`/`I64`/`F32`/`F64`); raw access and BSCALE/BZERO-scaled access (`scaled_values`). |
| `BinTable` / `BinColumn` / `BinColumnType` / `BinCellValue` | Binary-table model — columns, rows, and heap (variable-length arrays); built with `BinTableBuilder`. |
| `AsciiTable` | ASCII `TABLE` model with `TFORMn` column parsing and `TSCALn`/`TZEROn` scaling. |
| `CompressedImage` / `CompressionType` / `TileGeometry` / `Quantize` | Read-side view over a tile-compressed image: detect the algorithm and geometry, then `decompress()` to an `ImageData`. |
| `CompressOptions` | Write-side encode options (algorithm, tiling, quantization, dithering) for `ImageData::compress` / `compress_image`. |
| `Bitpix` | The `BITPIX` data type (`8`/`16`/`32`/`64`/`-32`/`-64`). |
| `Error` / `Result` | Crate error type and result alias. |

Decoding is done eagerly at read time (big-endian → native), except tile-compressed
images, which are decoded lazily on `decompress()` so the original compressed tiles
are preserved for a lossless round-trip write.

## Feature flags

| Flag | Default | Description |
|------|---------|-------------|
| *(none)* | ✓ | Core read/write, all HDU types, RICE_1 / PLIO_1 / HCOMPRESS_1 decompression, and RICE_1 compression — zero dependencies |
| `image` | | Conversion between `ImageData` and the [`image`](https://crates.io/crates/image) crate's `DynamicImage` |
| `gzip` | | `GZIP_1`/`GZIP_2` tile compression and decompression via the pure-Rust [`miniz_oxide`](https://crates.io/crates/miniz_oxide) crate |
| `wcs` | | Two-axis celestial World Coordinate System pixel ⇄ world transforms (`Wcs`, `Header::wcs`) via the zero-dependency [`mapproj`](https://crates.io/crates/mapproj) crate |

The default build stays dependency-free; `RICE_1`, `PLIO_1`, and `HCOMPRESS_1` decompression and `RICE_1` compression all work without any feature (only `GZIP_1`/`GZIP_2` need the `gzip` feature).

## Why fitskit?

Among Rust FITS crates, fitskit fills a specific niche — **pure Rust, zero default dependencies, full read + write including tables, complete compressed-image reading, and compressed-image writing**, with no C toolchain to install. No other pure-Rust crate writes compressed FITS.

| Crate | Pure Rust | Write | Tables | Compressed read | Compressed write | Notes |
|-------|-----------|-------|--------|-----------------|------------------|-------|
| [`fitsio`](https://crates.io/crates/fitsio) | ✗ | ✓ | ✓ | ✓ | ✓ | Wraps the cfitsio C library; needs a C toolchain |
| [`fitsrs`](https://crates.io/crates/fitsrs) | ✓ | ✗ | partial | — | ✗ | Read-only |
| [`fitrs`](https://crates.io/crates/fitrs) | ✓ | ✓ | ✗ | ✗ | ✗ | Dormant; no table support |
| **fitskit** | ✓ | ✓ | ✓ | ✓ (all types) | ✓ (RICE/GZIP) | Zero default deps; no C dependency |

## Supported / not supported

**Supported**

- Primary + extension HDUs; images, ASCII tables, binary tables (read + write)
- All BITPIX types; BSCALE/BZERO scaling and the unsigned-integer convention
- Variable-length arrays (`P`/`Q` descriptors) in binary tables
- CHECKSUM/DATASUM computation and verification
- Tile-compressed image **reading**: RICE_1, GZIP_1, GZIP_2, PLIO_1, and
  HCOMPRESS_1 — for integer and floating-point images, including quantization
  with subtractive dithering (`NO_DITHER` / `SUBTRACTIVE_DITHER_1` /
  `SUBTRACTIVE_DITHER_2`)
- Tile-compressed image **writing**: RICE_1 and GZIP_1/GZIP_2 — integer
  (lossless) and float (lossless via GZIP, or lossy quantized + dithered);
  output verified byte-exact through cfitsio's `funpack`
- Two-axis celestial WCS pixel ⇄ world transforms (`wcs` feature) for the common
  `CTYPEn = xxx--CCC` case, linear transform from `CDi_j` or `PCi_j` + `CDELTi`

**Not supported**

- Encoding `PLIO_1` or `HCOMPRESS_1` (these decode only)
- HCOMPRESS image smoothing (`SMOOTH` ≠ 0) on decode
- Random groups (deprecated in FITS v4.0)
- WCS: SIP distortions, 3+-axis / spectral WCS, `PVi_m` projection parameters,
  and non-degree `CUNIT` (the `wcs` feature handles the 2-axis celestial linear
  + projection case only)

## License

[MIT](LICENSE)
</content>
</invoke>
