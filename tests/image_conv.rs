#![cfg(feature = "image")]

use fits4::*;
use image::{DynamicImage, GrayImage, ImageBuffer, Luma};
type Gray16Image = ImageBuffer<Luma<u16>, Vec<u16>>;

#[test]
fn u8_image_round_trip() {
    let gray = GrayImage::from_fn(4, 3, |x, y| Luma([(x * 30 + y * 60) as u8]));
    let dyn_img = DynamicImage::ImageLuma8(gray.clone());

    let (img_data, bscale, bzero) = ImageData::from_dynamic_image(&dyn_img).unwrap();
    assert_eq!(img_data.bitpix(), Bitpix::U8);
    assert_eq!(img_data.axes, vec![4, 3]);
    assert!((bscale - 1.0).abs() < 1e-10);
    assert!((bzero - 0.0).abs() < 1e-10);

    let back = img_data.to_dynamic_image(bscale, bzero).unwrap();
    match back {
        DynamicImage::ImageLuma8(g) => {
            assert_eq!(g.dimensions(), (4, 3));
            assert_eq!(g.as_raw(), gray.as_raw());
        }
        _ => panic!("expected Luma8"),
    }
}

#[test]
fn u16_image_round_trip() {
    let gray16 = Gray16Image::from_fn(5, 4, |x, y| Luma([(x * 1000 + y * 2000) as u16]));
    let dyn_img = DynamicImage::ImageLuma16(gray16.clone());

    let (img_data, bscale, bzero) = ImageData::from_dynamic_image(&dyn_img).unwrap();
    assert_eq!(img_data.bitpix(), Bitpix::I16);
    assert!((bzero - 32768.0).abs() < 1e-10);

    let back = img_data.to_dynamic_image(bscale, bzero).unwrap();
    match back {
        DynamicImage::ImageLuma16(g) => {
            assert_eq!(g.dimensions(), (5, 4));
            assert_eq!(g.as_raw(), gray16.as_raw());
        }
        _ => panic!("expected Luma16"),
    }
}

#[test]
fn u16_full_range() {
    // Test u16 min/max values survive the round trip
    let pixels = vec![0u16, 1, 32767, 32768, 65534, 65535];
    let gray16 = Gray16Image::from_raw(6, 1, pixels.clone()).unwrap();
    let dyn_img = DynamicImage::ImageLuma16(gray16);

    let (img_data, bscale, bzero) = ImageData::from_dynamic_image(&dyn_img).unwrap();
    let back = img_data.to_dynamic_image(bscale, bzero).unwrap();

    match back {
        DynamicImage::ImageLuma16(g) => {
            assert_eq!(g.as_raw().as_slice(), &pixels);
        }
        _ => panic!("expected Luma16"),
    }
}

#[test]
fn f32_normalizes_to_u16() {
    let img_data = ImageData::new(vec![4], PixelData::F32(vec![0.0, 1.0, 0.5, 0.25]));

    let back = img_data.to_dynamic_image(1.0, 0.0).unwrap();
    match back {
        DynamicImage::ImageLuma16(g) => {
            let raw = g.as_raw();
            // 0.0 → 0, 1.0 → 65535, 0.5 → ~32767, 0.25 → ~16383
            assert_eq!(raw[0], 0);
            assert_eq!(raw[1], 65535);
            assert!((raw[2] as i32 - 32767).abs() < 2);
            assert!((raw[3] as i32 - 16383).abs() < 2);
        }
        _ => panic!("expected Luma16"),
    }
}

#[test]
fn fits_file_with_image_crate() {
    // Full pipeline: image crate → FITS → bytes → FITS → image crate
    let gray = GrayImage::from_fn(10, 8, |x, y| Luma([((x + y) % 256) as u8]));
    let dyn_img = DynamicImage::ImageLuma8(gray.clone());

    let (img_data, bscale, bzero) = ImageData::from_dynamic_image(&dyn_img).unwrap();
    let mut fits = FitsFile::with_primary_image(img_data);
    fits.primary_mut()
        .header
        .set("BSCALE", HeaderValue::Float(bscale), None);
    fits.primary_mut()
        .header
        .set("BZERO", HeaderValue::Float(bzero), None);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    let bs = fits2.primary().header.get_float("BSCALE").unwrap_or(1.0);
    let bz = fits2.primary().header.get_float("BZERO").unwrap_or(0.0);

    if let HduData::Image(img) = &fits2.primary().data {
        let back = img.to_dynamic_image(bs, bz).unwrap();
        match back {
            DynamicImage::ImageLuma8(g) => {
                assert_eq!(g.as_raw(), gray.as_raw());
            }
            _ => panic!("expected Luma8"),
        }
    } else {
        panic!("expected image");
    }
}
