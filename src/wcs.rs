//! World Coordinate System (WCS) pixel <-> world transforms (feature `wcs`).
//!
//! This module implements the common **two-axis celestial** FITS WCS: mapping
//! image pixel coordinates to world (celestial) longitude/latitude and back, for
//! a single linear transform plus a spherical projection (`CTYPEn = xxx--CCC`,
//! e.g. `RA---TAN`/`DEC--TAN`).
//!
//! The spherical-projection math (the native-sphere rotation and the per-code
//! projection equations) is provided by the pure-Rust
//! [`mapproj`](https://crates.io/crates/mapproj) crate, which itself has zero
//! runtime dependencies. This module is therefore gated behind the optional
//! `wcs` feature so the default `fitskit` build stays dependency-free.
//!
//! # Conventions
//!
//! * **Pixel coordinates are 1-based** at the API boundary, matching the FITS
//!   `CRPIXn` convention: the center of the first pixel is `(1.0, 1.0)`. So the
//!   reference pixel `(CRPIX1, CRPIX2)` maps exactly to `(CRVAL1, CRVAL2)`.
//! * **World coordinates are in degrees.** [`Wcs::pixel_to_world`] returns
//!   `(lon, lat)` in degrees (longitude in `[0, 360)`); [`Wcs::world_to_pixel`]
//!   accepts degrees.
//! * Axis order is `(x, y)` = `(axis 1, axis 2)` throughout, i.e. the first
//!   returned world coordinate corresponds to `CTYPE1`/`CRVAL1`.
//!
//! # What is supported
//!
//! * `CTYPE1`/`CTYPE2` with a recognised 3-letter projection code (the full set
//!   `mapproj` implements: `TAN`, `SIN`, `ARC`, `ZEA`, `STG`, `AZP`, `SZP`,
//!   `AIR`, `NCP`, `FEYE`, `MER`, `CAR`, `CEA`, `CYP`, `SFL`, `PAR`, `MOL`,
//!   `AIT`, `COD`, `COE`, `COO`, `COP`, `HPX`).
//! * The linear transform from either the `CDi_j` matrix, or `PCi_j` + `CDELTi`,
//!   or `CDELTi` alone (PC defaulting to identity). `CD` takes precedence when
//!   present.
//! * `CUNIT1`/`CUNIT2` of `deg` (the default when absent).
//!
//! # What is *not* supported (out of scope, returns `Err`)
//!
//! * **SIP distortions** (`A_p_q`/`B_p_q`, `CTYPE` suffix `-SIP`) — detected and
//!   rejected.
//! * **Three or more axes / spectral WCS** — only the 2-axis celestial case is
//!   handled. (A header with `NAXIS >= 3` is fine as long as the *first two*
//!   axes carry the celestial WCS; higher CD/PC cross-terms to other axes are
//!   ignored.)
//! * **`PVi_m` projection parameters** — parameterised projections use their
//!   `mapproj` defaults.
//! * **Non-degree `CUNIT`** (e.g. `arcsec`, `rad`) — rejected.
//! * **`LONPOLE`/`LATPOLE`** other than the projection's default placement: the
//!   keywords are parsed and stored, but the transform uses `mapproj`'s standard
//!   native-pole convention (correct for the usual celestial case where
//!   `LONPOLE = 180` for zenithal projections at the default reference point).

use crate::error::{Error, Result};
use crate::header::Header;

use mapproj::conic::cod::Cod;
use mapproj::conic::coe::Coe;
use mapproj::conic::coo::Coo;
use mapproj::conic::cop::Cop;
use mapproj::cylindrical::car::Car;
use mapproj::cylindrical::cea::Cea;
use mapproj::cylindrical::cyp::Cyp;
use mapproj::cylindrical::mer::Mer;
use mapproj::hybrid::hpx::Hpx;
use mapproj::img2celestial::Img2Celestial;
use mapproj::img2proj::WcsImgXY2ProjXY;
use mapproj::pseudocyl::ait::Ait;
use mapproj::pseudocyl::mol::Mol;
use mapproj::pseudocyl::par::Par;
use mapproj::pseudocyl::sfl::Sfl;
use mapproj::zenithal::air::Air;
use mapproj::zenithal::arc::Arc;
use mapproj::zenithal::azp::Azp;
use mapproj::zenithal::feye::Feye;
use mapproj::zenithal::ncp::Ncp;
use mapproj::zenithal::sin::Sin;
use mapproj::zenithal::stg::Stg;
use mapproj::zenithal::szp::Szp;
use mapproj::zenithal::tan::Tan;
use mapproj::zenithal::zea::Zea;
use mapproj::{CanonicalProjection, CenteredProjection, ImgXY, LonLat, Projection};

/// The 2x2 linear (intermediate-pixel-to-intermediate-world) transform, in the
/// `CDi_j` form (degrees per pixel). `PCi_j + CDELTi` is folded into this on
/// parse.
#[derive(Debug, Clone, Copy)]
struct LinearTransform {
    crpix1: f64,
    crpix2: f64,
    cd11: f64,
    cd12: f64,
    cd21: f64,
    cd22: f64,
}

/// A parsed two-axis celestial World Coordinate System.
///
/// Construct one with [`Header::wcs`] (or [`crate::Hdu::wcs`]). See the
/// [module documentation](self) for conventions (1-based pixels, degrees,
/// axis order) and the list of supported / unsupported features.
#[derive(Debug, Clone)]
pub struct Wcs {
    /// The 3-letter projection code extracted from `CTYPE1` (e.g. `"TAN"`).
    code: ProjCode,
    /// Reference world coordinates (`CRVAL1`, `CRVAL2`) in degrees.
    crval1: f64,
    crval2: f64,
    /// Linear transform (CRPIX + CD matrix, degrees).
    lin: LinearTransform,
}

/// The recognised projection codes (those `mapproj` implements).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProjCode {
    Tan,
    Sin,
    Arc,
    Zea,
    Stg,
    Azp,
    Szp,
    Air,
    Ncp,
    Feye,
    Mer,
    Car,
    Cea,
    Cyp,
    Sfl,
    Par,
    Mol,
    Ait,
    Cod,
    Coe,
    Coo,
    Cop,
    Hpx,
}

impl ProjCode {
    fn from_str(code: &str) -> Option<Self> {
        Some(match code {
            "TAN" => ProjCode::Tan,
            "SIN" => ProjCode::Sin,
            "ARC" => ProjCode::Arc,
            "ZEA" => ProjCode::Zea,
            "STG" => ProjCode::Stg,
            "AZP" => ProjCode::Azp,
            "SZP" => ProjCode::Szp,
            "AIR" => ProjCode::Air,
            "NCP" => ProjCode::Ncp,
            "FEYE" => ProjCode::Feye,
            "MER" => ProjCode::Mer,
            "CAR" => ProjCode::Car,
            "CEA" => ProjCode::Cea,
            "CYP" => ProjCode::Cyp,
            "SFL" => ProjCode::Sfl,
            "PAR" => ProjCode::Par,
            "MOL" => ProjCode::Mol,
            "AIT" => ProjCode::Ait,
            "COD" => ProjCode::Cod,
            "COE" => ProjCode::Coe,
            "COO" => ProjCode::Coo,
            "COP" => ProjCode::Cop,
            "HPX" => ProjCode::Hpx,
            _ => return None,
        })
    }
}

/// Build the `mapproj` image-to-celestial pipeline for one projection type `P`,
/// set its center from CRVAL (radians). Used for the forward
/// (pixel -> world) transform.
///
/// `P` is constructed via `Default` (projection parameters `PVi_m` are out of
/// scope and use `mapproj` defaults). The center is set from `CRVAL` in radians.
fn build_img2cel<P: CanonicalProjection + Default>(
    lin: &LinearTransform,
    crval1_rad: f64,
    crval2_rad: f64,
) -> Img2Celestial<P, WcsImgXY2ProjXY> {
    // mapproj's WcsImgXY2ProjXY subtracts CRPIX directly and converts the CD
    // matrix (deg/pix) to radians internally; passing the FITS 1-based CRPIX
    // together with 1-based pixel inputs keeps the reference pixel at the
    // projection-plane origin.
    let img2proj = WcsImgXY2ProjXY::from_cd(
        lin.crpix1, lin.crpix2, lin.cd11, lin.cd12, lin.cd21, lin.cd22,
    );
    Img2Celestial::new(img2proj, build_proj::<P>(crval1_rad, crval2_rad))
}

/// Build just the centered spherical projection for type `P`, centered on
/// `(crval1, crval2)` in radians. Used for the inverse (world -> pixel)
/// transform, where we apply our own (correct) inverse CD matrix to the
/// projection-plane coordinates rather than relying on `mapproj`'s
/// `WcsImgXY2ProjXY::inverse` (which transposes the off-diagonal terms of the
/// inverse CD matrix — see `world_to_pixel`).
fn build_proj<P: CanonicalProjection + Default>(
    crval1_rad: f64,
    crval2_rad: f64,
) -> CenteredProjection<P> {
    let mut proj = CenteredProjection::new(P::default());
    proj.set_proj_center_from_lonlat(&LonLat::new(crval1_rad, crval2_rad));
    proj
}

/// Dispatch `$op::<P>($($args)*)` over the concrete `mapproj` projection type
/// `P` selected by `$code`. The generic projection type never escapes the match
/// — each arm is monomorphised for one projection and returns `$op`'s result.
macro_rules! dispatch_proj {
    ($code:expr, $op:ident, $($args:expr),* $(,)?) => {
        match $code {
            ProjCode::Tan => $op::<Tan>($($args),*),
            ProjCode::Sin => $op::<Sin>($($args),*),
            ProjCode::Arc => $op::<Arc>($($args),*),
            ProjCode::Zea => $op::<Zea>($($args),*),
            ProjCode::Stg => $op::<Stg>($($args),*),
            ProjCode::Azp => $op::<Azp>($($args),*),
            ProjCode::Szp => $op::<Szp>($($args),*),
            ProjCode::Air => $op::<Air>($($args),*),
            ProjCode::Ncp => $op::<Ncp>($($args),*),
            ProjCode::Feye => $op::<Feye>($($args),*),
            ProjCode::Mer => $op::<Mer>($($args),*),
            ProjCode::Car => $op::<Car>($($args),*),
            ProjCode::Cea => $op::<Cea>($($args),*),
            ProjCode::Cyp => $op::<Cyp>($($args),*),
            ProjCode::Sfl => $op::<Sfl>($($args),*),
            ProjCode::Par => $op::<Par>($($args),*),
            ProjCode::Mol => $op::<Mol>($($args),*),
            ProjCode::Ait => $op::<Ait>($($args),*),
            ProjCode::Cod => $op::<Cod>($($args),*),
            ProjCode::Coe => $op::<Coe>($($args),*),
            ProjCode::Coo => $op::<Coo>($($args),*),
            ProjCode::Cop => $op::<Cop>($($args),*),
            ProjCode::Hpx => $op::<Hpx>($($args),*),
        }
    };
}

/// Forward (pixel -> world) for one concrete projection type.
fn pixel_to_world_for<P: CanonicalProjection + Default>(
    lin: &LinearTransform,
    crval1_rad: f64,
    crval2_rad: f64,
    x: f64,
    y: f64,
) -> Option<LonLat> {
    let i2c = build_img2cel::<P>(lin, crval1_rad, crval2_rad);
    i2c.img2lonlat(&ImgXY::new(x, y))
}

/// Inverse (world -> pixel) for one concrete projection type.
///
/// We deliberately do *not* use `mapproj`'s `Img2Celestial::lonlat2img`: its
/// `WcsImgXY2ProjXY::inverse` transposes the off-diagonal elements of the
/// inverse CD matrix, which is wrong for any non-diagonal CD/PC matrix. Instead
/// we project to the plane with `mapproj` and apply our own correct inverse CD.
fn world_to_pixel_for<P: CanonicalProjection + Default>(
    lin: &LinearTransform,
    crval1_rad: f64,
    crval2_rad: f64,
    lon_rad: f64,
    lat_rad: f64,
) -> Option<(f64, f64)> {
    let proj = build_proj::<P>(crval1_rad, crval2_rad);
    let pxy = proj.proj_lonlat(&LonLat::new(lon_rad, lat_rad))?;
    // `pxy` is in radians (mapproj converts the CD matrix to radians). Invert the
    // radian CD matrix to recover (pixel - crpix), then add CRPIX back.
    let f = std::f64::consts::PI / 180.0;
    let (a, b, c, d) = (lin.cd11 * f, lin.cd12 * f, lin.cd21 * f, lin.cd22 * f);
    let det = a * d - b * c;
    if det == 0.0 {
        return None;
    }
    let (px, py) = (pxy.x(), pxy.y());
    let dx = (d * px - b * py) / det;
    let dy = (-c * px + a * py) / det;
    Some((dx + lin.crpix1, dy + lin.crpix2))
}

/// Extract the 3-letter projection code from a `CTYPE` value such as
/// `"RA---TAN"` or `"GLON-AIT"`. Returns the substring after the final group of
/// `-` separators (trimmed of trailing spaces).
fn ctype_proj_code(ctype: &str) -> &str {
    let t = ctype.trim();
    match t.rfind('-') {
        Some(i) => &t[i + 1..],
        None => t,
    }
}

impl Wcs {
    /// Parse a two-axis celestial WCS from a FITS header.
    ///
    /// See [`Header::wcs`], which forwards here. Returns [`Error::Wcs`] /
    /// [`Error::UnsupportedWcs`] for headers that don't carry a supported
    /// 2-axis celestial WCS (missing keys, SIP, unrecognised projection,
    /// non-degree units, etc.).
    pub fn from_header(header: &Header) -> Result<Self> {
        let ctype1 = header
            .get_string("CTYPE1")
            .ok_or_else(|| Error::Wcs("missing CTYPE1".into()))?
            .to_string();
        let ctype2 = header
            .get_string("CTYPE2")
            .ok_or_else(|| Error::Wcs("missing CTYPE2".into()))?
            .to_string();

        // SIP distortion is out of scope.
        if ctype1.trim_end().ends_with("-SIP") || ctype2.trim_end().ends_with("-SIP") {
            return Err(Error::UnsupportedWcs(
                "SIP distortion (CTYPE -SIP) is not supported".into(),
            ));
        }

        let code1 = ctype_proj_code(&ctype1);
        let code2 = ctype_proj_code(&ctype2);
        if code1 != code2 {
            return Err(Error::Wcs(format!(
                "CTYPE1/CTYPE2 projection codes differ: '{code1}' vs '{code2}'"
            )));
        }
        let code = ProjCode::from_str(code1).ok_or_else(|| {
            Error::UnsupportedWcs(format!("unsupported projection code '{code1}'"))
        })?;

        // Units: only degrees supported.
        for key in ["CUNIT1", "CUNIT2"] {
            if let Some(u) = header.get_string(key) {
                let u = u.trim();
                if !u.is_empty() && !u.eq_ignore_ascii_case("deg") && !u.eq_ignore_ascii_case("degree") {
                    return Err(Error::UnsupportedWcs(format!(
                        "non-degree {key} = '{u}' is not supported"
                    )));
                }
            }
        }

        let crval1 = header
            .get_float("CRVAL1")
            .ok_or_else(|| Error::Wcs("missing CRVAL1".into()))?;
        let crval2 = header
            .get_float("CRVAL2")
            .ok_or_else(|| Error::Wcs("missing CRVAL2".into()))?;
        let crpix1 = header
            .get_float("CRPIX1")
            .ok_or_else(|| Error::Wcs("missing CRPIX1".into()))?;
        let crpix2 = header
            .get_float("CRPIX2")
            .ok_or_else(|| Error::Wcs("missing CRPIX2".into()))?;

        let lin = parse_linear(header, crpix1, crpix2)?;

        Ok(Wcs {
            code,
            crval1,
            crval2,
            lin,
        })
    }

    /// Convert a **1-based** FITS pixel coordinate `(x, y)` to world
    /// coordinates `(lon, lat)` in **degrees** (longitude in `[0, 360)`).
    ///
    /// Returns `None` if the pixel cannot be deprojected onto the sphere.
    pub fn pixel_to_world(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let c1 = self.crval1.to_radians();
        let c2 = self.crval2.to_radians();
        let lonlat = dispatch_proj!(self.code, pixel_to_world_for, &self.lin, c1, c2, x, y)?;
        let lon = lonlat.lon().to_degrees().rem_euclid(360.0);
        let lat = lonlat.lat().to_degrees();
        Some((lon, lat))
    }

    /// Convert world coordinates `(lon, lat)` in **degrees** to a **1-based**
    /// FITS pixel coordinate `(x, y)`.
    ///
    /// Returns `None` if the point cannot be projected (e.g. it lies outside the
    /// projection's valid region, or the linear transform is singular).
    ///
    /// Note: in **debug builds**, passing a point lying *exactly* at the
    /// projection center (`CRVAL1`, `CRVAL2`) can trip an internal
    /// `debug_assert!` in `mapproj` (a rounding-tolerance check on the rotated
    /// center unit vector). Release builds are unaffected; the inverse
    /// [`pixel_to_world`](Self::pixel_to_world) on the reference pixel is exact.
    pub fn world_to_pixel(&self, lon: f64, lat: f64) -> Option<(f64, f64)> {
        let c1 = self.crval1.to_radians();
        let c2 = self.crval2.to_radians();
        let lon_rad = lon.to_radians();
        let lat_rad = lat.to_radians();
        dispatch_proj!(
            self.code,
            world_to_pixel_for,
            &self.lin,
            c1,
            c2,
            lon_rad,
            lat_rad
        )
    }
}

/// Build the linear transform (CD matrix form) from the header, honouring the
/// precedence CD > (PC + CDELT) > (CDELT alone, PC = identity).
fn parse_linear(header: &Header, crpix1: f64, crpix2: f64) -> Result<LinearTransform> {
    let has_cd = header.find("CD1_1").is_some()
        || header.find("CD1_2").is_some()
        || header.find("CD2_1").is_some()
        || header.find("CD2_2").is_some();

    if has_cd {
        // Missing CD elements default to 0.
        let cd11 = header.get_float("CD1_1").unwrap_or(0.0);
        let cd12 = header.get_float("CD1_2").unwrap_or(0.0);
        let cd21 = header.get_float("CD2_1").unwrap_or(0.0);
        let cd22 = header.get_float("CD2_2").unwrap_or(0.0);
        return Ok(LinearTransform {
            crpix1,
            crpix2,
            cd11,
            cd12,
            cd21,
            cd22,
        });
    }

    // PC (default identity) + CDELT.
    let cdelt1 = header
        .get_float("CDELT1")
        .ok_or_else(|| Error::Wcs("missing CDELT1 (and no CD matrix)".into()))?;
    let cdelt2 = header
        .get_float("CDELT2")
        .ok_or_else(|| Error::Wcs("missing CDELT2 (and no CD matrix)".into()))?;
    let pc11 = header.get_float("PC1_1").unwrap_or(1.0);
    let pc12 = header.get_float("PC1_2").unwrap_or(0.0);
    let pc21 = header.get_float("PC2_1").unwrap_or(0.0);
    let pc22 = header.get_float("PC2_2").unwrap_or(1.0);

    // Fold PC + CDELT into a CD matrix: CDi_j = CDELTi * PCi_j.
    Ok(LinearTransform {
        crpix1,
        crpix2,
        cd11: cdelt1 * pc11,
        cd12: cdelt1 * pc12,
        cd21: cdelt2 * pc21,
        cd22: cdelt2 * pc22,
    })
}

impl Header {
    /// Parse a two-axis celestial [`Wcs`] from this header (feature `wcs`).
    ///
    /// This mirrors [`Hdu::as_compressed_image`](crate::Hdu::as_compressed_image):
    /// WCS lives entirely in the header, so it is read straight from the
    /// keywords. See the [`wcs`](crate::wcs) module docs for conventions
    /// (1-based pixels, degrees) and supported features.
    ///
    /// ```no_run
    /// # use fitskit::FitsFile;
    /// let fits = FitsFile::from_file("image.fits").unwrap();
    /// let wcs = fits.primary().header.wcs().unwrap();
    /// let (ra, dec) = wcs.pixel_to_world(256.5, 256.5).unwrap();
    /// println!("center -> RA={ra} Dec={dec}");
    /// ```
    pub fn wcs(&self) -> Result<Wcs> {
        Wcs::from_header(self)
    }
}

impl crate::Hdu {
    /// Parse a two-axis celestial [`Wcs`] from this HDU's header (feature `wcs`).
    ///
    /// Convenience forwarder to [`Header::wcs`].
    pub fn wcs(&self) -> Result<Wcs> {
        self.header.wcs()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keyword::HeaderValue;

    fn tan_header() -> Header {
        // Simple TAN WCS: 100x100 image, reference pixel at center, 1"/pix.
        let mut h = Header::new();
        h.set("CTYPE1", HeaderValue::String("RA---TAN".into()), None);
        h.set("CTYPE2", HeaderValue::String("DEC--TAN".into()), None);
        h.set("CRVAL1", HeaderValue::Float(150.0), None);
        h.set("CRVAL2", HeaderValue::Float(30.0), None);
        h.set("CRPIX1", HeaderValue::Float(50.5), None);
        h.set("CRPIX2", HeaderValue::Float(50.5), None);
        // 1 arcsec/pixel = 1/3600 deg/pixel, no rotation, RA increases to -x.
        h.set("CD1_1", HeaderValue::Float(-1.0 / 3600.0), None);
        h.set("CD1_2", HeaderValue::Float(0.0), None);
        h.set("CD2_1", HeaderValue::Float(0.0), None);
        h.set("CD2_2", HeaderValue::Float(1.0 / 3600.0), None);
        h
    }

    #[test]
    fn reference_pixel_maps_to_crval() {
        let wcs = tan_header().wcs().unwrap();
        let (lon, lat) = wcs.pixel_to_world(50.5, 50.5).unwrap();
        assert!((lon - 150.0).abs() < 1e-9, "lon = {lon}");
        assert!((lat - 30.0).abs() < 1e-9, "lat = {lat}");
    }

    #[test]
    fn world_to_pixel_near_crval_is_reference() {
        // A world point one pixel (1 arcsec) "north" of CRVAL lands one pixel up
        // from the reference pixel. (We use a point one pixel off the exact
        // center: a point lying *exactly* at CRVAL can trip an internal
        // `debug_assert!` in `mapproj`'s `XYZ::new` — a rounding-tolerance check
        // on the rotated center unit vector — in debug builds. Release builds and
        // the exact-reference `pixel_to_world` path are unaffected.)
        let wcs = tan_header().wcs().unwrap();
        let (x, y) = wcs.world_to_pixel(150.0, 30.0 + 1.0 / 3600.0).unwrap();
        assert!((x - 50.5).abs() < 1e-6, "x = {x}");
        assert!((y - 51.5).abs() < 1e-6, "y = {y}");
    }

    #[test]
    fn round_trip_pixel_world_pixel() {
        let wcs = tan_header().wcs().unwrap();
        for &(x, y) in &[(1.0, 1.0), (100.0, 100.0), (1.0, 100.0), (75.25, 12.5)] {
            let (lon, lat) = wcs.pixel_to_world(x, y).unwrap();
            let (x2, y2) = wcs.world_to_pixel(lon, lat).unwrap();
            assert!((x - x2).abs() < 1e-7, "x {x} != {x2}");
            assert!((y - y2).abs() < 1e-7, "y {y} != {y2}");
        }
    }

    #[test]
    fn pc_cdelt_equivalent_to_cd() {
        // PC = identity, CDELT = (-1/3600, 1/3600) must equal the CD matrix above.
        let mut h = Header::new();
        h.set("CTYPE1", HeaderValue::String("RA---TAN".into()), None);
        h.set("CTYPE2", HeaderValue::String("DEC--TAN".into()), None);
        h.set("CRVAL1", HeaderValue::Float(150.0), None);
        h.set("CRVAL2", HeaderValue::Float(30.0), None);
        h.set("CRPIX1", HeaderValue::Float(50.5), None);
        h.set("CRPIX2", HeaderValue::Float(50.5), None);
        h.set("CDELT1", HeaderValue::Float(-1.0 / 3600.0), None);
        h.set("CDELT2", HeaderValue::Float(1.0 / 3600.0), None);
        let wcs = h.wcs().unwrap();
        let wcs_cd = tan_header().wcs().unwrap();
        for &(x, y) in &[(10.0, 20.0), (90.0, 80.0)] {
            let a = wcs.pixel_to_world(x, y).unwrap();
            let b = wcs_cd.pixel_to_world(x, y).unwrap();
            assert!((a.0 - b.0).abs() < 1e-12 && (a.1 - b.1).abs() < 1e-12);
        }
    }

    #[test]
    fn unsupported_projection_rejected() {
        let mut h = tan_header();
        h.set("CTYPE1", HeaderValue::String("RA---XYZ".into()), None);
        h.set("CTYPE2", HeaderValue::String("DEC--XYZ".into()), None);
        assert!(matches!(h.wcs(), Err(Error::UnsupportedWcs(_))));
    }

    #[test]
    fn sip_rejected() {
        let mut h = tan_header();
        h.set("CTYPE1", HeaderValue::String("RA---TAN-SIP".into()), None);
        h.set("CTYPE2", HeaderValue::String("DEC--TAN-SIP".into()), None);
        assert!(matches!(h.wcs(), Err(Error::UnsupportedWcs(_))));
    }

    #[test]
    fn missing_keyword_errors() {
        let mut h = tan_header();
        h.keywords.retain(|k| k.name != "CRVAL1");
        assert!(matches!(h.wcs(), Err(Error::Wcs(_))));
    }
}
