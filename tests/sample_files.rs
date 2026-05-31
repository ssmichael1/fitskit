//! Tests reading NASA sample FITS files from the samp/ directory.
//! These tests are skipped if the samp/ directory is not present.

use fits4::*;
use std::path::Path;

const SAMP_DIR: &str = "samp";

fn samp(name: &str) -> std::path::PathBuf {
    Path::new(SAMP_DIR).join(name)
}

macro_rules! require_samples {
    () => {
        if !Path::new(SAMP_DIR).is_dir() {
            eprintln!("skipping: samp/ directory not present");
            return;
        }
    };
}

#[test]
fn read_euv_image() {
    require_samples!();
    let fits = FitsFile::from_file(samp("EUVEngc4151imgx.fits")).unwrap();

    let primary = fits.primary();
    assert_eq!(primary.header.get_bool("SIMPLE"), Some(true));
    assert_eq!(primary.header.get_int("BITPIX"), Some(8));
    assert_eq!(primary.header.get_int("NAXIS"), Some(0));
    assert!(matches!(primary.data, HduData::Empty));

    assert!(
        fits.len() > 1,
        "expected extensions, got {} HDUs",
        fits.len()
    );

    let mut image_count = 0;
    let mut bintable_count = 0;
    for hdu in fits.extensions() {
        match &hdu.data {
            HduData::Image(_) => image_count += 1,
            HduData::BinTable(_) => bintable_count += 1,
            _ => {}
        }
    }
    assert!(image_count > 0, "expected IMAGE extensions");
    assert!(bintable_count > 0, "expected BINTABLE extensions");

    for hdu in fits.extensions() {
        if let HduData::Image(img) = &hdu.data {
            assert_eq!(img.bitpix(), Bitpix::I16);
            assert_eq!(img.axes.len(), 2);
            assert!(img.width().unwrap() > 0);
            assert!(img.height().unwrap() > 0);
            break;
        }
    }
}

#[test]
fn read_fgs_with_ascii_table() {
    require_samples!();
    let fits = FitsFile::from_file(samp("FGSf64y0106m_a1f.fits")).unwrap();

    let primary = fits.primary();
    assert_eq!(primary.header.get_bool("SIMPLE"), Some(true));
    assert_eq!(primary.header.get_int("BITPIX"), Some(32));
    assert_eq!(primary.header.get_int("NAXIS"), Some(2));
    assert_eq!(primary.header.get_int("NAXIS1"), Some(89688));
    assert_eq!(primary.header.get_int("NAXIS2"), Some(7));

    if let HduData::Image(img) = &primary.data {
        assert_eq!(img.bitpix(), Bitpix::I32);
        assert_eq!(img.axes, vec![89688, 7]);
        assert_eq!(img.num_pixels(), 89688 * 7);
    } else {
        panic!("expected image data in primary HDU");
    }

    assert_eq!(fits.len(), 2, "expected primary + 1 extension");
    if let HduData::AsciiTable(table) = &fits.extensions()[0].data {
        assert_eq!(table.nrows, 7);
        assert_eq!(table.columns.len(), 6);
    } else {
        panic!("expected ASCII TABLE extension");
    }
}

#[test]
fn read_foc_float_image() {
    require_samples!();
    let fits = FitsFile::from_file(samp("FOCx38i0101t_c0f.fits")).unwrap();

    let primary = fits.primary();
    assert_eq!(primary.header.get_int("BITPIX"), Some(-32));
    assert_eq!(primary.header.get_int("NAXIS"), Some(2));
    assert_eq!(primary.header.get_int("NAXIS1"), Some(1024));
    assert_eq!(primary.header.get_int("NAXIS2"), Some(1024));

    if let HduData::Image(img) = &primary.data {
        assert_eq!(img.bitpix(), Bitpix::F32);
        assert_eq!(img.num_pixels(), 1024 * 1024);

        if let PixelData::F32(data) = &img.pixels {
            assert!(!data.is_empty());
            let non_zero = data.iter().filter(|&&v| v != 0.0).count();
            assert!(non_zero > 0, "expected some non-zero pixels");
        }
    } else {
        panic!("expected float image data");
    }

    assert_eq!(fits.len(), 2);
    if let HduData::AsciiTable(table) = &fits.extensions()[0].data {
        assert_eq!(table.nrows, 1);
        assert_eq!(table.columns.len(), 18);
    } else {
        panic!("expected ASCII TABLE extension");
    }
}

#[test]
fn read_iue_header_only() {
    require_samples!();
    let fits = FitsFile::from_file(samp("IUElwp25637mxlo.fits")).unwrap();

    let primary = fits.primary();
    assert_eq!(primary.header.get_bool("SIMPLE"), Some(true));
    assert_eq!(primary.header.get_int("BITPIX"), Some(8));
    assert_eq!(primary.header.get_int("NAXIS"), Some(0));
    assert!(matches!(primary.data, HduData::Empty));

    assert!(
        primary.header.find("TELESCOP").is_some() || primary.header.find("INSTRUME").is_some(),
        "expected instrument metadata"
    );
}

#[test]
fn read_wfpc2_3d_cube() {
    require_samples!();
    let fits = FitsFile::from_file(samp("WFPC2u5780205r_c0fx.fits")).unwrap();

    let primary = fits.primary();
    assert_eq!(primary.header.get_int("BITPIX"), Some(-32));
    assert_eq!(primary.header.get_int("NAXIS"), Some(3));
    assert_eq!(primary.header.get_int("NAXIS1"), Some(200));
    assert_eq!(primary.header.get_int("NAXIS2"), Some(200));
    assert_eq!(primary.header.get_int("NAXIS3"), Some(4));

    if let HduData::Image(img) = &primary.data {
        assert_eq!(img.bitpix(), Bitpix::F32);
        assert_eq!(img.axes, vec![200, 200, 4]);
        assert_eq!(img.num_pixels(), 200 * 200 * 4);
    } else {
        panic!("expected 3D image data");
    }

    assert_eq!(fits.len(), 2);
    if let HduData::AsciiTable(table) = &fits.extensions()[0].data {
        assert_eq!(table.nrows, 4);
        assert_eq!(table.columns.len(), 49);
    } else {
        panic!("expected ASCII TABLE extension");
    }
}

#[test]
fn all_sample_files_readable() {
    require_samples!();
    let fits_files = [
        "EUVEngc4151imgx.fits",
        "FGSf64y0106m_a1f.fits",
        "FOCx38i0101t_c0f.fits",
        "IUElwp25637mxlo.fits",
        "WFPC2u5780205r_c0fx.fits",
    ];

    for name in &fits_files {
        let path = samp(name);
        let result = FitsFile::from_file(&path);
        assert!(result.is_ok(), "failed to read {name}: {:?}", result.err());

        let fits = result.unwrap();
        assert!(!fits.is_empty(), "{name}: expected at least one HDU");
        assert_eq!(
            fits.primary().header.get_bool("SIMPLE"),
            Some(true),
            "{name}: missing SIMPLE=T"
        );
    }
}

// --- Tiled-image compression: RICE_1 / GZIP_1 round-trip vs originals --------
//
// Each compressed `.fz` fixture has the SAME HDU ordering as its uncompressed
// source: every `HduData::Image` HDU in the source appears, in order, as a
// compressed-image `HduData::BinTable` (ZIMAGE=T) in the `.fz`. We decompress
// each and assert byte-exact equality of the pixel buffers.
//
// These are guarded per-fixture with `Path::exists()` so they skip cleanly on a
// machine without the (gitignored, GCS-hosted) fixtures.

/// Assert that decompressing every compressed-image HDU in `fz_name` reproduces,
/// byte-for-byte, the corresponding `HduData::Image` HDUs of `src_name`.
fn assert_rice_roundtrip(src_name: &str, fz_name: &str) {
    let fz_path = samp(fz_name);
    if !fz_path.exists() {
        eprintln!("skipping: fixture {fz_name} not present");
        return;
    }
    let src = FitsFile::from_file(samp(src_name)).expect("read source");
    let fz = FitsFile::from_file(&fz_path).expect("read .fz");

    // Source image HDUs, in order.
    let src_images: Vec<&ImageData> = src
        .hdus
        .iter()
        .filter_map(|h| match &h.data {
            HduData::Image(im) if !im.pixels.to_bytes().is_empty() => Some(im),
            _ => None,
        })
        .collect();

    // Compressed-image HDUs in the .fz, in order.
    let mut matched = 0usize;
    let mut src_iter = src_images.iter();
    for hdu in &fz.hdus {
        if let Some(cimg) = hdu.as_compressed_image() {
            let orig = src_iter
                .next()
                .expect("more compressed images than source images");
            let dec = cimg.decompress().expect("decompress");
            assert_eq!(
                dec.axes, orig.axes,
                "{fz_name}: axes mismatch on compressed HDU #{matched}"
            );
            assert_eq!(
                dec.pixels.to_bytes(),
                orig.pixels.to_bytes(),
                "{fz_name}: pixel bytes differ on compressed HDU #{matched}"
            );
            matched += 1;
        }
    }
    assert!(matched > 0, "{fz_name}: found no compressed-image HDUs");
}

#[test]
fn rice_roundtrip_euv_row_tiled() {
    // 512x512 I16, RICE row tiling (512x1 = one tile per row).
    assert_rice_roundtrip("EUVEngc4151imgx.fits", "EUVEngc4151imgx.rice.fits.fz");
}

#[test]
fn rice_roundtrip_euv_square_tiled() {
    // Same source, forced 100x100 tiles: exercises 2-D edge-truncated tiling.
    assert_rice_roundtrip("EUVEngc4151imgx.fits", "EUVEngc4151imgx.rice_t100.fits.fz");
}

#[test]
fn rice_roundtrip_fgs_i32() {
    // 89688x7 I32 image.
    assert_rice_roundtrip("FGSf64y0106m_a1f.fits", "FGSf64y0106m_a1f.rice.fits.fz");
}

#[test]
fn plio_roundtrip_euv() {
    // PLIO_1 (IRAF pixel-list RLE) over the I16 EUV image extensions. Lossless,
    // so the decode must reproduce the source pixels byte-for-byte.
    assert_rice_roundtrip("EUVEngc4151imgx.fits", "EUVEngc4151imgx.plio.fits.fz");
}

#[test]
fn hcompress_int_roundtrip_euv() {
    // HCOMPRESS_1 with SCALE=0 is LOSSLESS for integer images, so the decode of
    // the I16 EUV extensions must be byte-exact vs the source. Exercises both the
    // 512x16 and 2048x16 tile geometries.
    assert_rice_roundtrip(
        "EUVEngc4151imgx.fits",
        "EUVEngc4151imgx.hcomp_int.fits.fz",
    );
}

// --- float tile-compression: bit-exact vs funpack -------------------------
//
// Quantized-float decompression is *lossy*, so the reconstructed floats will not
// equal the original uncompressed sample. The authoritative oracle is funpack's own
// reconstruction: fits4 and funpack must apply the identical dither table + per-tile
// scaling, so fits4's output must be BIT-IDENTICAL to funpack's. We invoke the
// `funpack` CLI at test time and compare; the test skips cleanly when either the
// fixture or the `funpack` binary is absent (so CI without cfitsio stays green).

/// True if a `funpack` binary is on PATH.
fn funpack_available() -> bool {
    std::process::Command::new("funpack")
        .arg("-V")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Decompress `fz_name` with the `funpack` CLI into a temp `.fits` and read it back.
fn funpack_reference(fz_name: &str) -> Option<FitsFile> {
    let fz_path = samp(fz_name);
    let mut out = std::env::temp_dir();
    out.push(format!(
        "fits4_funpack_ref_{}_{}.fits",
        std::process::id(),
        fz_name.replace(['/', '.'], "_")
    ));
    let _ = std::fs::remove_file(&out);
    let status = std::process::Command::new("funpack")
        .arg("-O")
        .arg(&out)
        .arg("-C") // don't update checksums
        .arg(&fz_path)
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let f = FitsFile::from_file(&out).ok();
    let _ = std::fs::remove_file(&out);
    f
}

/// Assert that fits4's float decompression of every compressed-image HDU in `fz_name`
/// is byte-exact against funpack's reconstruction of the same file.
fn assert_float_matches_funpack(fz_name: &str) {
    let fz_path = samp(fz_name);
    if !fz_path.exists() {
        eprintln!("skipping: fixture {fz_name} not present");
        return;
    }
    if !funpack_available() {
        eprintln!("skipping: funpack binary not available");
        return;
    }
    let reference = match funpack_reference(fz_name) {
        Some(r) => r,
        None => {
            eprintln!("skipping: funpack failed on {fz_name}");
            return;
        }
    };

    // funpack reference image HDUs (the reconstructed floats), in order.
    let ref_images: Vec<&ImageData> = reference
        .hdus
        .iter()
        .filter_map(|h| match &h.data {
            HduData::Image(im) if !im.pixels.to_bytes().is_empty() => Some(im),
            _ => None,
        })
        .collect();

    let fz = FitsFile::from_file(&fz_path).expect("read .fz");
    let mut ref_iter = ref_images.iter();
    let mut matched = 0usize;
    for hdu in &fz.hdus {
        if let Some(cimg) = hdu.as_compressed_image() {
            let reference_img = ref_iter
                .next()
                .expect("more compressed images than funpack reference images");
            let dec = cimg.decompress().expect("fits4 float decompress");
            assert_eq!(
                dec.axes, reference_img.axes,
                "{fz_name}: axes mismatch on compressed HDU #{matched}"
            );
            compare_floats_vs_funpack(fz_name, matched, &dec.pixels, &reference_img.pixels);
            matched += 1;
        }
    }
    assert!(matched > 0, "{fz_name}: found no compressed-image HDUs");
}

/// Compare fits4's reconstructed floats against funpack's, element-by-element.
///
/// fits4 reproduces cfitsio's unquantization *exactly*, including the fused
/// multiply-add (`x*scale + zero` contracted to a single FMA) that the cfitsio C code
/// relies on, so every reconstructed pixel must be bit-identical to funpack's. NaNs
/// (from `ZBLANK`) only have to match as NaN — IEEE leaves the payload unspecified.
fn compare_floats_vs_funpack(fz_name: &str, hdu: usize, got: &PixelData, want: &PixelData) {
    fn cmp(fz_name: &str, hdu: usize, a: f64, b: f64) {
        if a.is_nan() && b.is_nan() {
            return;
        }
        assert!(
            a.to_bits() == b.to_bits(),
            "{fz_name} HDU#{hdu}: fits4={a} ({:#x}) != funpack={b} ({:#x})",
            a.to_bits(),
            b.to_bits()
        );
    }
    match (got, want) {
        (PixelData::F32(g), PixelData::F32(w)) => {
            assert_eq!(g.len(), w.len(), "{fz_name} HDU#{hdu}: length mismatch");
            for (x, y) in g.iter().zip(w.iter()) {
                cmp(fz_name, hdu, *x as f64, *y as f64);
            }
        }
        (PixelData::F64(g), PixelData::F64(w)) => {
            assert_eq!(g.len(), w.len(), "{fz_name} HDU#{hdu}: length mismatch");
            for (x, y) in g.iter().zip(w.iter()) {
                cmp(fz_name, hdu, *x, *y);
            }
        }
        _ => panic!("{fz_name} HDU#{hdu}: pixel type mismatch vs funpack reference"),
    }
}

#[test]
fn float_rice_nodither_matches_funpack() {
    // 1024x1024 F32, RICE_1, ZQUANTIZ=NO_DITHER.
    assert_float_matches_funpack("FOCx38i0101t_c0f.rice_nodith.fits.fz");
}

#[cfg(feature = "gzip")]
#[test]
fn float_gzip_lossless_matches_funpack() {
    // 1024x1024 F32, GZIP_1, ZQUANTIZ=NONE (raw floats, no quantization). Needs gzip.
    assert_float_matches_funpack("FOCx38i0101t_c0f.gzip_lossless.fits.fz");
}

#[test]
fn float_rice_dither1_matches_funpack() {
    // 1024x1024 F32, RICE_1, ZQUANTIZ=SUBTRACTIVE_DITHER_1.
    assert_float_matches_funpack("FOCx38i0101t_c0f.rice_dith.fits.fz");
}

#[test]
fn float_rice_dither2_matches_funpack() {
    // 1024x1024 F32, RICE_1, ZQUANTIZ=SUBTRACTIVE_DITHER_2 (preserves exact zeros),
    // generated by scripts/gen_compressed_fixtures.sh with `fpack -r -qz5 16`.
    assert_float_matches_funpack("FOCx38i0101t_c0f.rice_dith2.fits.fz");
}

#[test]
fn float_cube_rice_dither1_matches_funpack() {
    // 200x200x4 F32 cube, RICE_1, SUBTRACTIVE_DITHER_1.
    assert_float_matches_funpack("WFPC2u5780205r_c0fx.rice_dith.fits.fz");
}

#[test]
fn float_cube_rice_nodither_matches_funpack() {
    // 200x200x4 F32 cube, RICE_1, NO_DITHER.
    assert_float_matches_funpack("WFPC2u5780205r_c0fx.rice_nodith.fits.fz");
}

#[test]
fn float_hcompress_dither1_matches_funpack() {
    // 1024x1024 F32, HCOMPRESS_1 (SCALE=0 on the quantized ints),
    // ZQUANTIZ=SUBTRACTIVE_DITHER_1. Lossy (quantized), so validate bit-exact
    // against funpack's own reconstruction. Tiles are 1024x16.
    assert_float_matches_funpack("FOCx38i0101t_c0f.hcomp.fits.fz");
}

#[cfg(feature = "gzip")]
#[test]
fn gzip1_roundtrip_euv() {
    // GZIP_1 stores raw big-endian image integers per tile; needs the `gzip` feature.
    assert_rice_roundtrip("EUVEngc4151imgx.fits", "EUVEngc4151imgx.gzip1.fits.fz");
}
