//! World Coordinate System (WCS) tests (feature `wcs`).
//!
//! Two layers:
//!  * Always-runnable: a hand-built TAN WCS with analytic checks, plus
//!    pixel->world->pixel round-trips. These need no external data.
//!  * Oracle: compare `Wcs` against `astropy.wcs.WCS` on the NASA sample files
//!    that carry a full 2-axis celestial TAN WCS. Skipped cleanly when either
//!    the `samp/` files or a Python `astropy` are absent (mirroring the
//!    `funpack`-oracle skip pattern in `tests/sample_files.rs`).
#![cfg(feature = "wcs")]

use fitskit::{FitsFile, Header, HeaderValue, Wcs};
use std::path::Path;

const SAMP_DIR: &str = "samp";

fn samp(name: &str) -> std::path::PathBuf {
    Path::new(SAMP_DIR).join(name)
}

// ---------------------------------------------------------------------------
// Always-runnable: hand-built TAN WCS.
// ---------------------------------------------------------------------------

fn tan_header() -> Header {
    let mut h = Header::new();
    h.set("CTYPE1", HeaderValue::String("RA---TAN".into()), None);
    h.set("CTYPE2", HeaderValue::String("DEC--TAN".into()), None);
    h.set("CRVAL1", HeaderValue::Float(150.0), None);
    h.set("CRVAL2", HeaderValue::Float(30.0), None);
    h.set("CRPIX1", HeaderValue::Float(50.5), None);
    h.set("CRPIX2", HeaderValue::Float(50.5), None);
    h.set("CD1_1", HeaderValue::Float(-1.0 / 3600.0), None);
    h.set("CD1_2", HeaderValue::Float(0.0), None);
    h.set("CD2_1", HeaderValue::Float(0.0), None);
    h.set("CD2_2", HeaderValue::Float(1.0 / 3600.0), None);
    h
}

#[test]
fn tan_reference_pixel_is_crval() {
    let wcs = tan_header().wcs().unwrap();
    let (lon, lat) = wcs.pixel_to_world(50.5, 50.5).unwrap();
    assert!((lon - 150.0).abs() < 1e-9);
    assert!((lat - 30.0).abs() < 1e-9);
}

#[test]
fn tan_round_trip_identity() {
    let wcs = tan_header().wcs().unwrap();
    for &(x, y) in &[(1.0, 1.0), (100.0, 100.0), (1.0, 100.0), (33.0, 78.0)] {
        let (lon, lat) = wcs.pixel_to_world(x, y).unwrap();
        let (x2, y2) = wcs.world_to_pixel(lon, lat).unwrap();
        assert!((x - x2).abs() < 1e-7, "x {x} != {x2}");
        assert!((y - y2).abs() < 1e-7, "y {y} != {y2}");
    }
}

// ---------------------------------------------------------------------------
// Oracle: astropy comparison on the celestial-TAN sample files.
// ---------------------------------------------------------------------------

/// True if `python3 -c "import astropy"` succeeds.
fn astropy_available() -> bool {
    std::process::Command::new("python3")
        .args(["-c", "import astropy.wcs"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Ask astropy for (ra, dec) in degrees at 1-based pixel (px, py) for HDU `hdu`
/// of `file`. Returns None if python/astropy can't be run.
fn astropy_pix2world(file: &Path, hdu: usize, pts: &[(f64, f64)]) -> Option<Vec<(f64, f64)>> {
    let coords: String = pts
        .iter()
        .map(|(x, y)| format!("({x},{y})"))
        .collect::<Vec<_>>()
        .join(",");
    let script = format!(
        r#"
import sys, warnings
warnings.simplefilter('ignore')
from astropy.io import fits
from astropy.wcs import WCS
w = WCS(fits.open(r"{}")[{}].header, naxis=2)
pts = [{}]
for px, py in pts:
    ra, dec = w.all_pix2world(px, py, 1)
    print(float(ra) % 360.0, float(dec))
"#,
        file.display(),
        hdu,
        coords
    );
    let out = std::process::Command::new("python3")
        .args(["-c", &script])
        .output()
        .ok()?;
    if !out.status.success() {
        eprintln!(
            "astropy script failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        return None;
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let mut v = Vec::new();
    for line in text.lines() {
        let mut it = line.split_whitespace();
        let ra: f64 = it.next()?.parse().ok()?;
        let dec: f64 = it.next()?.parse().ok()?;
        v.push((ra, dec));
    }
    if v.len() == pts.len() {
        Some(v)
    } else {
        None
    }
}

/// Load the `Wcs` for HDU `hdu` of a sample file.
fn sample_wcs(file: &str, hdu: usize) -> Wcs {
    let fits = FitsFile::from_file(samp(file)).expect("read sample");
    fits.hdus[hdu].header.wcs().expect("parse WCS")
}

/// Compare `Wcs::pixel_to_world` against astropy at the given 1-based pixels,
/// then confirm `world_to_pixel` inverts `pixel_to_world` to tight tolerance.
fn assert_matches_astropy(file: &str, hdu: usize, pts: &[(f64, f64)]) {
    let path = samp(file);
    if !path.exists() {
        eprintln!("skipping: sample {file} not present");
        return;
    }
    if !astropy_available() {
        eprintln!("skipping: astropy not available");
        return;
    }
    let oracle = match astropy_pix2world(&path, hdu, pts) {
        Some(o) => o,
        None => {
            eprintln!("skipping: astropy oracle failed for {file}");
            return;
        }
    };

    let wcs = sample_wcs(file, hdu);
    for (&(px, py), &(ra_a, dec_a)) in pts.iter().zip(oracle.iter()) {
        let (ra, dec) = wcs
            .pixel_to_world(px, py)
            .unwrap_or_else(|| panic!("{file}: pixel_to_world None at ({px},{py})"));
        assert!(
            (ra - ra_a).abs() < 1e-6 && (dec - dec_a).abs() < 1e-6,
            "{file} HDU#{hdu} pix ({px},{py}): fitskit ({ra},{dec}) vs astropy ({ra_a},{dec_a})"
        );
        // Inverse must return to the input pixel.
        let (x2, y2) = wcs
            .world_to_pixel(ra, dec)
            .unwrap_or_else(|| panic!("{file}: world_to_pixel None"));
        assert!(
            (x2 - px).abs() < 1e-6 && (y2 - py).abs() < 1e-6,
            "{file}: round-trip pixel ({px},{py}) -> ({x2},{y2})"
        );
    }
}

#[test]
fn astropy_euv_tan_cdelt() {
    // EUVEngc4151imgx HDU 1: 512x512 RA---TAN/DEC--TAN with CDELT (no CD matrix).
    assert_matches_astropy(
        "EUVEngc4151imgx.fits",
        1,
        &[(1.0, 1.0), (256.5, 256.5), (512.0, 512.0), (1.0, 512.0)],
    );
}

#[test]
fn astropy_foc_tan_cd() {
    // FOCx38i0101t_c0f HDU 0: 1024x1024 RA---TAN/DEC--TAN with a full CD matrix.
    assert_matches_astropy(
        "FOCx38i0101t_c0f.fits",
        0,
        &[(1.0, 1.0), (512.0, 512.0), (1024.0, 1024.0), (1.0, 1024.0)],
    );
}

#[test]
fn astropy_wfpc2_tan_cd_first_two_axes() {
    // WFPC2u5780205r_c0fx HDU 0: NAXIS=3 cube, but the first two axes carry a
    // RA---TAN/DEC--TAN CD-matrix WCS. We use only those two axes.
    assert_matches_astropy(
        "WFPC2u5780205r_c0fx.fits",
        0,
        &[(1.0, 1.0), (100.0, 100.0), (200.0, 200.0)],
    );
}

#[test]
fn euv_spectral_hdu_is_unsupported() {
    // EUV HDU 2 is WAVELENGTH/ANGLE (not celestial) -> unsupported projection.
    let path = samp("EUVEngc4151imgx.fits");
    if !path.exists() {
        eprintln!("skipping: sample not present");
        return;
    }
    let fits = FitsFile::from_file(&path).unwrap();
    assert!(fits.hdus[2].header.wcs().is_err());
}
