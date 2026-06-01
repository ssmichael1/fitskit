#[cfg(feature = "image")]
use crate::error::{Error, Result};
#[cfg(feature = "image")]
use crate::image_data::{ImageData, PixelData};
#[cfg(feature = "image")]
use image::{DynamicImage, GrayImage, ImageBuffer, Luma};
type Gray16Image = ImageBuffer<Luma<u16>, Vec<u16>>;

#[cfg(feature = "image")]
impl ImageData {
    /// Convert to a `DynamicImage` from the `image` crate.
    ///
    /// - BITPIX=8 → `GrayImage` (Luma8)
    /// - BITPIX=16, BZERO=32768 → `ImageLuma16` (unsigned u16)
    /// - BITPIX=16, signed → `ImageLuma16` with BSCALE/BZERO applied, clamped to 0..65535
    /// - All others → normalize min..max to `ImageLuma16`
    pub fn to_dynamic_image(&self, bscale: f64, bzero: f64) -> Result<DynamicImage> {
        let width = self
            .width()
            .ok_or_else(|| Error::InvalidFormat("need at least 1 axis".into()))?
            as u32;
        let height = self.height().unwrap_or(1) as u32;

        match &self.pixels {
            PixelData::U8(data) => {
                let img = GrayImage::from_raw(width, height, data.clone())
                    .ok_or_else(|| Error::InvalidFormat("image dimensions mismatch".into()))?;
                Ok(DynamicImage::ImageLuma8(img))
            }
            PixelData::I16(data) => {
                let pixels: Vec<u16> =
                    if (bzero - 32768.0).abs() < 0.5 && (bscale - 1.0).abs() < 1e-10 {
                        // Unsigned u16 convention
                        data.iter().map(|&v| (v as i32 + 32768) as u16).collect()
                    } else {
                        data.iter()
                            .map(|&v| (bzero + bscale * v as f64).clamp(0.0, 65535.0) as u16)
                            .collect()
                    };
                let img = Gray16Image::from_raw(width, height, pixels)
                    .ok_or_else(|| Error::InvalidFormat("image dimensions mismatch".into()))?;
                Ok(DynamicImage::ImageLuma16(img))
            }
            _ => {
                // Normalize to u16 range
                let scaled = self.scaled_values(bscale, bzero);
                let (min, max) = scaled
                    .iter()
                    .fold((f64::MAX, f64::MIN), |(mn, mx), &v| (mn.min(v), mx.max(v)));
                let range = if (max - min).abs() < 1e-30 {
                    1.0
                } else {
                    max - min
                };

                let pixels: Vec<u16> = scaled
                    .iter()
                    .map(|&v| ((v - min) / range * 65535.0) as u16)
                    .collect();
                let img = Gray16Image::from_raw(width, height, pixels)
                    .ok_or_else(|| Error::InvalidFormat("image dimensions mismatch".into()))?;
                Ok(DynamicImage::ImageLuma16(img))
            }
        }
    }

    /// Create `ImageData` from a `DynamicImage`.
    ///
    /// Returns `(image_data, bscale, bzero)`:
    /// - Luma8 → BITPIX=8, BSCALE=1, BZERO=0
    /// - Luma16 → BITPIX=16, BSCALE=1, BZERO=32768 (unsigned convention)
    /// - Others → converted to Luma16 first
    pub fn from_dynamic_image(img: &DynamicImage) -> Result<(Self, f64, f64)> {
        match img {
            DynamicImage::ImageLuma8(gray) => {
                let (w, h) = gray.dimensions();
                let data = gray.as_raw().clone();
                Ok((
                    ImageData::new(vec![w as usize, h as usize], PixelData::U8(data)),
                    1.0,
                    0.0,
                ))
            }
            DynamicImage::ImageLuma16(gray16) => Ok(luma16_to_image(gray16)),
            other => Ok(luma16_to_image(&other.to_luma16())),
        }
    }
}

/// Convert a 16-bit grayscale image to `ImageData` using the unsigned-16 FITS
/// convention (`BZERO=32768`). Returns `(image_data, bscale, bzero)`.
fn luma16_to_image(gray16: &Gray16Image) -> (ImageData, f64, f64) {
    let (w, h) = gray16.dimensions();
    let pixels: Vec<i16> = gray16
        .as_raw()
        .iter()
        .map(|&v| (v as i32 - 32768) as i16)
        .collect();
    (
        ImageData::new(vec![w as usize, h as usize], PixelData::I16(pixels)),
        1.0,
        32768.0,
    )
}
