# fits4

**Pure-Rust, zero-dependency reader and writer for [FITS](https://fits.gsfc.nasa.gov/fits_standard.html) v4.0** — the Flexible Image Transport System format used throughout astronomy.

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

fits4 reads *and* writes the full FITS v4.0 standard — primary and extension HDUs, images, ASCII tables, and binary tables (including variable-length arrays) — with **no external dependencies** in the default build and **no C toolchain** required. It also decodes every tile-compressed image type in the standard (RICE, GZIP, PLIO, HCOMPRESS) on the read path.

## Features

- **All standard HDU types** — Primary, IMAGE, ASCII `TABLE`, `BINTABLE`
- **Full read *and* write** — round-trip files, bytes, or any `Read`/`Write`
- **All BITPIX types** — 8, 16, 32, 64, -32, -64
- **BSCALE/BZERO** — raw and scaled access, with the unsigned-integer convention
- **Variable-length arrays** — `P` and `Q` heap descriptors in binary tables
- **CHECKSUM/DATASUM** — ones-complement integrity computation and verification
- **Tile-compressed images (read)** — RICE_1, GZIP_1, GZIP_2, PLIO_1, and HCOMPRESS_1, including float quantization with subtractive dithering; decoded bit-exactly against cfitsio's `funpack`
- **Zero dependencies by default** — optional `image` and `gzip` features pull in pure-Rust crates only

## Installation

```toml
[dependencies]
fits4 = "0.1"
```

## Usage

### Reading a FITS file

```rust
use fits4::{FitsFile, HduData, PixelData};

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
# Ok::<(), fits4::Error>(())
```

### Writing a FITS file

```rust
use fits4::{FitsFile, ImageData, PixelData, HeaderValue};

let pixels: Vec<i16> = (0..10000).map(|i| (i % 1000) as i16).collect();
let img = ImageData::new(vec![100, 100], PixelData::I16(pixels));

let mut fits = FitsFile::with_primary_image(img);
fits.primary_mut().header.set("OBJECT", HeaderValue::String("M31".into()), None);

fits.to_file("output.fits")?;
# Ok::<(), fits4::Error>(())
```

### Building a binary table

```rust
use fits4::{FitsFile, Hdu, BinTableBuilder, BinColumnType};

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
use fits4::FitsFile;

let fits = FitsFile::from_file("compressed.fits")?;

for hdu in fits.extensions() {
    if let Some(cimg) = hdu.as_compressed_image() {
        println!("compression = {:?}", cimg.compression());
        let image = cimg.decompress()?;
        println!("decompressed to {:?}", image.axes);
    }
}
# Ok::<(), fits4::Error>(())
```

### Checksums

```rust
use fits4::{FitsFile, ImageData, PixelData};

let img = ImageData::new(vec![4], PixelData::U8(vec![1, 2, 3, 4]));
let fits = FitsFile::with_primary_image(img);

// Write with CHECKSUM/DATASUM keywords
let bytes = fits.to_bytes_with_checksum()?;

// Verify on read
let fits2 = FitsFile::from_bytes(&bytes)?;
fits2.primary().verify_datasum()?;
# Ok::<(), fits4::Error>(())
```

## Feature flags

| Flag | Default | Description |
|------|---------|-------------|
| *(none)* | ✓ | Core read/write, all HDU types, RICE_1 / PLIO_1 / HCOMPRESS_1 tile decompression — zero dependencies |
| `image` | | Conversion between `ImageData` and the [`image`](https://crates.io/crates/image) crate's `DynamicImage` |
| `gzip` | | Decoding of `GZIP_1`/`GZIP_2` tile-compressed images via the pure-Rust [`miniz_oxide`](https://crates.io/crates/miniz_oxide) crate |

The default build stays dependency-free; `RICE_1`, `PLIO_1`, and `HCOMPRESS_1` decompression all work without any feature (only `GZIP_1`/`GZIP_2` need the `gzip` feature).

## Why fits4?

Among Rust FITS crates, fits4 fills a specific niche — **pure Rust, zero default dependencies, full read + write including tables, and complete compressed-image read support**, with no C toolchain to install.

| Crate | Pure Rust | Write | Tables | Compressed read | Notes |
|-------|-----------|-------|--------|-----------------|-------|
| [`fitsio`](https://crates.io/crates/fitsio) | ✗ | ✓ | ✓ | ✓ | Wraps the cfitsio C library; needs a C toolchain |
| [`fitsrs`](https://crates.io/crates/fitsrs) | ✓ | ✗ | partial | — | Read-only |
| [`fitrs`](https://crates.io/crates/fitrs) | ✓ | ✓ | ✗ | ✗ | Dormant; no table support |
| **fits4** | ✓ | ✓ | ✓ | ✓ (all types) | Zero default deps; no C dependency |

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

**Not supported**

- Writing/encoding tile-compressed images (read/decode only)
- HCOMPRESS image smoothing (`SMOOTH` ≠ 0) on decode
- Random groups (deprecated in FITS v4.0)

## License

[MIT](LICENSE)
</content>
</invoke>
