//! # fitskit — Pure Rust FITS v4.0 reader/writer
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
//! Tile-compressed images (`ZIMAGE`) are supported for **reading** integer and float
//! images (`RICE_1`, `PLIO_1`, `HCOMPRESS_1`; `GZIP_1`/`GZIP_2` behind the `gzip`
//! feature), including float quantization with `SUBTRACTIVE_DITHER_1`/`_2`, via
//! [`Hdu::as_compressed_image`](hdu::Hdu::as_compressed_image) +
//! [`CompressedImage::decompress`](tile_compress::CompressedImage::decompress).
//!
//! **Writing** tile-compressed images is supported for `RICE_1` (integer + quantized
//! float) and `GZIP_1`/`GZIP_2` (integer; lossless float via `GZIP_1`) via
//! [`ImageData::compress`](image_data::ImageData::compress) /
//! [`compress_image`](tile_compress::compress_image), producing standard, cfitsio-
//! readable (`funpack`-decodable) compressed FITS. `PLIO_1`/`HCOMPRESS_1` encoding and
//! random groups are not supported.
//!
//! ## Quick start — reading
//!
//! ```no_run
//! use fitskit::{FitsFile, HduData, PixelData};
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
//! use fitskit::{FitsFile, Hdu, ImageData, PixelData, HeaderValue};
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
//! ## Binary tables
//!
//! Build a `BINTABLE` extension column-by-column, then push it onto a file:
//!
//! ```
//! use fitskit::{FitsFile, Hdu, BinTableBuilder, BinColumnType};
//!
//! let table = BinTableBuilder::new()
//!     .add_column("RA", BinColumnType::D64(1))
//!     .add_column("DEC", BinColumnType::D64(1))
//!     .add_column("MAG", BinColumnType::E32(1))
//!     .push_row(|row| {
//!         row.write_f64(180.0);
//!         row.write_f64(45.0);
//!         row.write_f32(12.5);
//!     })
//!     .build();
//!
//! let mut fits = FitsFile::with_empty_primary();
//! fits.push_extension(Hdu::bintable_extension(table));
//! let _ = fits.to_bytes().unwrap();
//! ```
//!
//! ## Tile-compressed images
//!
//! [`ImageData::compress`](image_data::ImageData::compress) produces a compressed-image
//! `BINTABLE` HDU; [`Hdu::as_compressed_image`](hdu::Hdu::as_compressed_image) +
//! [`CompressedImage::decompress`](tile_compress::CompressedImage::decompress) read it
//! back. RICE_1 is lossless for integers, so this round-trips exactly:
//!
//! ```
//! use fitskit::{FitsFile, ImageData, PixelData, CompressOptions, HduData};
//!
//! let pixels: Vec<i16> = (0..10000).map(|i| (i % 1000) as i16).collect();
//! let img = ImageData::new(vec![100, 100], PixelData::I16(pixels.clone()));
//!
//! // Compress (default: RICE_1, one tile per row) into a BINTABLE HDU and write it.
//! let mut fits = FitsFile::with_empty_primary();
//! fits.push_extension(img.compress(&CompressOptions::default()).unwrap());
//! let bytes = fits.to_bytes().unwrap();
//!
//! // Read back and decompress.
//! let fits2 = FitsFile::from_bytes(&bytes).unwrap();
//! let cimg = fits2.extensions()[0].as_compressed_image().unwrap();
//! let restored = cimg.decompress().unwrap();
//! assert!(matches!(restored.pixels, PixelData::I16(ref v) if *v == pixels));
//! ```
//!
//! `RICE_1`, `PLIO_1`, and `HCOMPRESS_1` work in the default build. The `GZIP_1`/
//! `GZIP_2` algorithms (selected via `CompressOptions { algorithm:
//! CompressionType::Gzip1, .. }`, or encountered when decoding a GZIP-compressed
//! tile) require the `gzip` feature.
//!
//! ## BSCALE/BZERO
//!
//! Physical values are computed as `BZERO + BSCALE * array_value`. The unsigned
//! integer convention stores unsigned values in signed storage:
//!
//! ```
//! use fitskit::{ImageData, PixelData};
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
//! use fitskit::{FitsFile, ImageData, PixelData};
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
//! ## Core types
//!
//! A FITS file is an ordered sequence of Header-Data Units (HDUs):
//!
//! - [`FitsFile`] — top-level container ([`Vec<Hdu>`](Hdu)); reads/writes files,
//!   bytes, or any `Read`/`Write`, plus builder methods.
//! - [`Hdu`] — one HDU: a [`Header`] plus an [`HduData`] payload.
//! - [`HduData`] — payload enum: `Empty`, `Image`, `AsciiTable`, `BinTable`.
//! - [`Header`] / [`Keyword`] / [`HeaderValue`] — ordered cards with typed accessors.
//! - [`ImageData`] / [`PixelData`] — axes plus a typed pixel buffer; raw and
//!   BSCALE/BZERO-scaled access.
//! - [`BinTable`] / [`BinColumn`] / [`BinColumnType`] / [`BinCellValue`] — binary-table
//!   model (with VLA heap), built via [`BinTableBuilder`].
//! - [`AsciiTable`] — ASCII `TABLE` model.
//! - [`CompressedImage`] / [`CompressionType`] / [`CompressOptions`] — tile-compression
//!   read view and write options.
//! - [`Bitpix`] — the `BITPIX` data type; [`Error`] / [`Result`] — error handling.
//!
//! ## Image-crate interop (feature `image`)
//!
//! With the `image` feature, [`ImageData`] converts to and from the
//! [`image`](https://crates.io/crates/image) crate's `DynamicImage` — e.g. to
//! save a FITS image as a PNG, or to ingest a raster as FITS.
//!
//! ```no_run
//! # #[cfg(feature = "image")] {
//! use fitskit::{FitsFile, HduData, ImageData};
//!
//! let fits = FitsFile::from_file("image.fits").unwrap();
//! if let HduData::Image(img) = &fits.primary().data {
//!     // FITS -> image crate (BSCALE/BZERO applied; 1.0/0.0 = identity)
//!     let dynamic = img.to_dynamic_image(1.0, 0.0).unwrap();
//!     // image crate -> FITS, returning (ImageData, bscale, bzero)
//!     let (restored, _bscale, _bzero) = ImageData::from_dynamic_image(&dynamic).unwrap();
//!     assert_eq!(restored.axes, img.axes);
//! }
//! # }
//! ```
//!
//! ## World Coordinate System (WCS)
//!
//! With the `wcs` feature, parse the celestial WCS from a header and map between
//! 1-based FITS pixel coordinates and world coordinates (degrees). The spherical
//! projection math is backed by the zero-dependency
//! [`mapproj`](https://crates.io/crates/mapproj) crate.
//!
//! ```no_run
//! # #[cfg(feature = "wcs")] {
//! use fitskit::FitsFile;
//!
//! let fits = FitsFile::from_file("image.fits").unwrap();
//! let wcs = fits.primary().header.wcs().unwrap();
//!
//! // Pixel (1-based) -> world (lon, lat) in degrees
//! let (ra, dec) = wcs.pixel_to_world(256.5, 256.5).unwrap();
//!
//! // ...and back
//! let (x, y) = wcs.world_to_pixel(ra, dec).unwrap();
//! # }
//! ```
//!
//! Only the common two-axis celestial case (e.g. `RA---TAN`/`DEC--TAN`) is
//! supported; SIP distortions and 3+-axis/spectral WCS are out of scope. See the
//! [`wcs`] module docs for details.
//!
//! ## Feature flags
//!
//! - **`image`** — enables conversion between [`ImageData`] and the
//!   [`image`](https://crates.io/crates/image) crate's `DynamicImage`.
//! - **`gzip`** — enables encoding *and* decoding of `GZIP_1`/`GZIP_2` tile-compressed
//!   images via the pure-Rust [`miniz_oxide`](https://crates.io/crates/miniz_oxide)
//!   crate. The default build stays dependency-free; `RICE_1` (read and write),
//!   `PLIO_1`, and `HCOMPRESS_1` decoding all work without this feature.
//! - **`wcs`** — enables two-axis celestial World Coordinate System pixel <-> world
//!   transforms ([`Wcs`], [`Header::wcs`](header::Header::wcs)) backed by the
//!   zero-dependency [`mapproj`](https://crates.io/crates/mapproj) crate. The default
//!   build stays dependency-free.

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

#[cfg(feature = "wcs")]
pub mod wcs;

pub use ascii_table::AsciiTable;
pub use bintable::{BinCellValue, BinColumn, BinColumnType, BinTable, BinTableBuilder};
pub use error::{Error, Result};
pub use fits::FitsFile;
pub use hdu::{Hdu, HduData};
pub use header::Header;
pub use image_data::{ImageData, PixelData};
pub use keyword::{HeaderValue, Keyword};
pub use tile_compress::{
    compress_image, CompressOptions, CompressedImage, CompressionType, Quantize, TileGeometry,
};
pub use types::Bitpix;

#[cfg(feature = "wcs")]
pub use wcs::Wcs;
