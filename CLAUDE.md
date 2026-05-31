# fits4 — Pure Rust FITS v4.0 Library

## Build Commands

```bash
cargo build
cargo test
cargo test --features image   # with image crate integration
cargo clippy
```

## Architecture

Zero external dependencies for core functionality. Optional `image` crate behind feature flag.

### Module Structure

| Module | Purpose |
|---|---|
| `error.rs` | `Error` enum, `Result` type alias |
| `types.rs` | `Bitpix` enum, constants (`BLOCK_SIZE=2880`, `RECORD_SIZE=80`) |
| `keyword.rs` | `Keyword` struct, `HeaderValue` enum, 80-byte card parse/serialize, CONTINUE handling |
| `header.rs` | `Header` (ordered keyword vec), typed accessors, block I/O |
| `io_utils.rs` | Block-aligned reading/writing, padding helpers |
| `image_data.rs` | `ImageData`, `PixelData` enum, BSCALE/BZERO scaling |
| `ascii_table.rs` | `AsciiTable`: TFORMn parsing (Aw/Iw/Fw.d/Ew.d/Dw.d), column access, TSCALn/TZEROn |
| `bintable.rs` | `BinTable`: all type codes (L,X,B,I,J,K,A,E,D,C,M,P,Q), heap/VLA |
| `checksum.rs` | CHECKSUM/DATASUM ones-complement computation |
| `tile_compress.rs` | Tiled-image compression. Decode: RICE_1, GZIP_1/2 (feature `gzip`), PLIO_1, HCOMPRESS_1; quantization + subtractive dithering for floats. Encode (`compress_image`/`ImageData::compress`): RICE_1 (int + quantized float) and GZIP_1/2 (int; lossless float via GZIP_1) |
| `hdu.rs` | `Hdu` struct, `HduData` enum (Empty/Image/AsciiTable/BinTable); `as_compressed_image()` accessor |
| `fits.rs` | `FitsFile`: top-level read/write, HDU iteration, builder API |
| `image_conv.rs` | (feature="image") `DynamicImage` <-> `ImageData` conversion |

### Key Types

- **`FitsFile`** — `Vec<Hdu>`, read/write from files, bytes, or readers
- **`Hdu`** — header + data payload (`HduData` enum)
- **`Header`** — ordered `Vec<Keyword>` with typed accessors (`get_int`, `get_float`, `get_string`, `get_bool`)
- **`ImageData`** — axes + `PixelData` enum (U8/I16/I32/I64/F32/F64)
- **`BinTable`** — columns + main data + heap (for VLAs)
- **`AsciiTable`** — columns + raw ASCII data

### Reading Pipeline

1. Parse 2880-byte header blocks → extract 36×80-byte keywords → stop at END
2. Compute data size from header keywords
3. Read padded data block
4. Dispatch by HDU type (SIMPLE → primary, XTENSION → extension)
5. Decode big-endian bytes to native at read time

### Writing Pipeline

1. Build header with mandatory keywords → serialize to 80-byte cards → pad to 2880
2. Convert native values to big-endian bytes
3. Pad data to 2880-byte boundary
4. Write header blocks then data blocks per HDU

## Test Files

NASA sample FITS files in `samp/`:
- `EUVEngc4151imgx.fits` — NAXIS=0 primary + IMAGE(I16) + BINTABLE extensions
- `FGSf64y0106m_a1f.fits` — I32 89688×7 image + ASCII TABLE (6 cols, 7 rows)
- `FOCx38i0101t_c0f.fits` — F32 1024×1024 image + ASCII TABLE (18 cols, 1 row)
- `IUElwp25637mxlo.fits` — Header-only (NAXIS=0, no extensions)
- `WFPC2u5780205r_c0fx.fits` — F32 200×200×4 cube + ASCII TABLE (49 cols, 4 rows)

## Design Decisions

- **No external dependencies** for core — no `byteorder`, no `thiserror`. Only optional, feature-gated deps: `image` and `gzip` (`miniz_oxide`)
- **Random groups** skipped (deprecated per standard)
- BSCALE/BZERO: raw vs scaled access modes on `ImageData`
- Unsigned integer convention: BZERO offset (32768 for u16, etc.)
- **Tile compression**: decode (read) for RICE_1/PLIO_1/HCOMPRESS_1 in the zero-dep core, GZIP_1/2 behind the `gzip` feature. Lazy `hdu.as_compressed_image()?.decompress()?`; `HduData` stays `BinTable` so the compressed tiles survive for lossless round-trip. Float decode reproduces cfitsio's fused multiply-add to stay bit-exact vs `funpack` (see memory). Encode (write) via `image.compress(&CompressOptions{..})? -> Hdu` (then `fits.push_extension(..)`): RICE_1 (int lossless + quantized/dithered float lossy) and GZIP_1/2 (int lossless; lossless raw-float storage via GZIP_1). Encoders are byte-exact inverses of the decoders and emit `funpack`-readable FITS. **Z\* keyword order is load-bearing**: `funpack` rebuilds the image header by walking cards, so `ZTENSION` must precede `ZBITPIX`/`ZNAXIS` (else "1st key not SIMPLE or XTENSION"); `build_z_header` emits fpack's order. PLIO_1/HCOMPRESS_1 encode and HCOMPRESS `SMOOTH≠0` are not implemented
- **Compressed fixtures**: `scripts/gen_compressed_fixtures.sh` builds fpack `.fz` test files into `samp/` (gitignored, served from GCS bucket `fits4_samples`); compression tests skip when fixtures/`funpack` are absent
