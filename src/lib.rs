//! # fits4 — Pure Rust FITS v4.0 reader/writer
//!
//! A zero-dependency implementation of the FITS (Flexible Image Transport System)
//! standard v4.0 for reading and writing astronomical data files.
//!
//! ## Supported HDU types
//!
//! - **Primary HDU** — image data or header-only
//! - **IMAGE extension** — all BITPIX types (8, 16, 32, 64, -32, -64)
//! - **ASCII TABLE extension** — Aw, Iw, Fw.d, Ew.d, Dw.d column formats
//! - **BINTABLE extension** — all type codes including variable-length arrays (P/Q descriptors)
//!
//! Tile-compressed images (`ZIMAGE`) are supported for **reading** integer images
//! (`RICE_1`; `GZIP_1`/`GZIP_2` behind the `gzip` feature) via
//! [`Hdu::as_compressed_image`](hdu::Hdu::as_compressed_image) +
//! [`CompressedImage::decompress`](tile_compress::CompressedImage::decompress).
//! Float quantization/dithering, `PLIO_1`, `HCOMPRESS_1`, and the write/compress
//! path are not yet implemented. Random groups are not supported.
//!
//! ## Quick start — reading
//!
//! ```no_run
//! use fits4::{FitsFile, HduData, PixelData};
//!
//! let fits = FitsFile::from_file("image.fits").unwrap();
//!
//! // Access the primary HDU
//! let primary = fits.primary();
//! println!("BITPIX = {}", primary.header.get_int("BITPIX").unwrap());
//!
//! if let HduData::Image(img) = &primary.data {
//!     println!("{}x{}", img.width().unwrap(), img.height().unwrap());
//!     if let PixelData::F32(pixels) = &img.pixels {
//!         println!("first pixel = {}", pixels[0]);
//!     }
//! }
//!
//! // Iterate over extensions
//! for hdu in fits.extensions() {
//!     match &hdu.data {
//!         HduData::Image(_) => println!("IMAGE extension"),
//!         HduData::BinTable(t) => println!("BINTABLE: {} rows", t.nrows),
//!         HduData::AsciiTable(t) => println!("TABLE: {} rows", t.nrows),
//!         HduData::Empty => println!("empty"),
//!     }
//! }
//! ```
//!
//! ## Quick start — writing
//!
//! ```
//! use fits4::{FitsFile, Hdu, ImageData, PixelData, HeaderValue};
//!
//! // Create a 100x100 16-bit image
//! let pixels: Vec<i16> = (0..10000).map(|i| (i % 1000) as i16).collect();
//! let img = ImageData::new(vec![100, 100], PixelData::I16(pixels));
//!
//! let mut fits = FitsFile::with_primary_image(img);
//! fits.primary_mut().header.set("OBJECT", HeaderValue::String("M31".into()), None);
//!
//! let bytes = fits.to_bytes().unwrap();
//! assert_eq!(bytes.len() % 2880, 0); // block-aligned
//! ```
//!
//! ## BSCALE/BZERO
//!
//! Physical values are computed as `BZERO + BSCALE * array_value`. The unsigned
//! integer convention stores unsigned values in signed storage:
//!
//! ```
//! use fits4::{ImageData, PixelData};
//!
//! // Unsigned u16 via BZERO=32768
//! let img = ImageData::new(vec![3], PixelData::I16(vec![-32768, 0, 32767]));
//! let physical = img.scaled_values(1.0, 32768.0);
//! assert_eq!(physical, vec![0.0, 32768.0, 65535.0]);
//! ```
//!
//! ## Checksums
//!
//! Write with CHECKSUM/DATASUM integrity keywords:
//!
//! ```
//! use fits4::{FitsFile, ImageData, PixelData};
//!
//! let img = ImageData::new(vec![4], PixelData::U8(vec![1, 2, 3, 4]));
//! let fits = FitsFile::with_primary_image(img);
//! let bytes = fits.to_bytes_with_checksum().unwrap();
//!
//! // Read back and verify
//! let fits2 = FitsFile::from_bytes(&bytes).unwrap();
//! fits2.primary().verify_datasum().unwrap();
//! ```
//!
//! ## Feature flags
//!
//! - **`image`** — enables conversion between [`ImageData`] and the
//!   [`image`](https://crates.io/crates/image) crate's `DynamicImage`.
//! - **`gzip`** — enables decoding of `GZIP_1`/`GZIP_2` tile-compressed images via
//!   the pure-Rust [`miniz_oxide`](https://crates.io/crates/miniz_oxide) crate. The
//!   default build stays dependency-free; `RICE_1` works without this feature.

pub mod ascii_table;
pub mod bintable;
pub mod checksum;
pub mod error;
pub mod fits;
pub mod hdu;
pub mod header;
pub mod image_data;
pub mod io_utils;
pub mod keyword;
pub mod tile_compress;
pub mod types;

#[cfg(feature = "image")]
pub mod image_conv;

pub use ascii_table::AsciiTable;
pub use bintable::{BinCellValue, BinColumn, BinColumnType, BinTable, BinTableBuilder};
pub use error::{Error, Result};
pub use fits::FitsFile;
pub use hdu::{Hdu, HduData};
pub use header::Header;
pub use image_data::{ImageData, PixelData};
pub use keyword::{HeaderValue, Keyword};
pub use tile_compress::{CompressedImage, CompressionType, Quantize, TileGeometry};
pub use types::Bitpix;
