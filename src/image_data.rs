use crate::error::{Error, Result};
use crate::header::Header;
use crate::keyword::HeaderValue;
use crate::types::Bitpix;

/// Raw pixel data storage.
#[derive(Debug, Clone)]
pub enum PixelData {
    U8(Vec<u8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
}

impl PixelData {
    pub fn len(&self) -> usize {
        match self {
            PixelData::U8(v) => v.len(),
            PixelData::I16(v) => v.len(),
            PixelData::I32(v) => v.len(),
            PixelData::I64(v) => v.len(),
            PixelData::F32(v) => v.len(),
            PixelData::F64(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn bitpix(&self) -> Bitpix {
        match self {
            PixelData::U8(_) => Bitpix::U8,
            PixelData::I16(_) => Bitpix::I16,
            PixelData::I32(_) => Bitpix::I32,
            PixelData::I64(_) => Bitpix::I64,
            PixelData::F32(_) => Bitpix::F32,
            PixelData::F64(_) => Bitpix::F64,
        }
    }

    /// Convert pixel data to big-endian bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        match self {
            PixelData::U8(v) => v.clone(),
            PixelData::I16(v) => v.iter().flat_map(|x| x.to_be_bytes()).collect(),
            PixelData::I32(v) => v.iter().flat_map(|x| x.to_be_bytes()).collect(),
            PixelData::I64(v) => v.iter().flat_map(|x| x.to_be_bytes()).collect(),
            PixelData::F32(v) => v.iter().flat_map(|x| x.to_be_bytes()).collect(),
            PixelData::F64(v) => v.iter().flat_map(|x| x.to_be_bytes()).collect(),
        }
    }

    /// Decode big-endian bytes into pixel data.
    pub fn from_bytes(bitpix: Bitpix, data: &[u8]) -> Result<Self> {
        let bpv = bitpix.bytes_per_value();
        if !data.len().is_multiple_of(bpv) {
            return Err(Error::DataSizeMismatch {
                expected: (data.len() / bpv + 1) * bpv,
                actual: data.len(),
            });
        }

        Ok(match bitpix {
            Bitpix::U8 => PixelData::U8(data.to_vec()),
            Bitpix::I16 => PixelData::I16(
                data.chunks_exact(2)
                    .map(|c| i16::from_be_bytes([c[0], c[1]]))
                    .collect(),
            ),
            Bitpix::I32 => PixelData::I32(
                data.chunks_exact(4)
                    .map(|c| i32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect(),
            ),
            Bitpix::I64 => PixelData::I64(
                data.chunks_exact(8)
                    .map(|c| i64::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            ),
            Bitpix::F32 => PixelData::F32(
                data.chunks_exact(4)
                    .map(|c| f32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect(),
            ),
            Bitpix::F64 => PixelData::F64(
                data.chunks_exact(8)
                    .map(|c| f64::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            ),
        })
    }
}

/// Image data with axes and scaling info.
#[derive(Debug, Clone)]
pub struct ImageData {
    pub axes: Vec<usize>,
    pub pixels: PixelData,
}

impl ImageData {
    pub fn new(axes: Vec<usize>, pixels: PixelData) -> Self {
        ImageData { axes, pixels }
    }

    pub fn bitpix(&self) -> Bitpix {
        self.pixels.bitpix()
    }

    pub fn num_pixels(&self) -> usize {
        self.axes.iter().product()
    }

    /// Width (NAXIS1) if 2D.
    pub fn width(&self) -> Option<usize> {
        self.axes.first().copied()
    }

    /// Height (NAXIS2) if 2D.
    pub fn height(&self) -> Option<usize> {
        self.axes.get(1).copied()
    }

    /// Get scaled pixel values as f64 using BSCALE and BZERO.
    pub fn scaled_values(&self, bscale: f64, bzero: f64) -> Vec<f64> {
        match &self.pixels {
            PixelData::U8(v) => v.iter().map(|&x| bzero + bscale * x as f64).collect(),
            PixelData::I16(v) => v.iter().map(|&x| bzero + bscale * x as f64).collect(),
            PixelData::I32(v) => v.iter().map(|&x| bzero + bscale * x as f64).collect(),
            PixelData::I64(v) => v.iter().map(|&x| bzero + bscale * x as f64).collect(),
            PixelData::F32(v) => v.iter().map(|&x| bzero + bscale * x as f64).collect(),
            PixelData::F64(v) => v.iter().map(|&x| bzero + bscale * x).collect(),
        }
    }

    /// Read image data from a header and raw bytes.
    pub fn from_header_and_data(header: &Header, data: &[u8]) -> Result<Self> {
        let bitpix_val = header.require_int("BITPIX")?;
        let bitpix = Bitpix::from_i64(bitpix_val)?;
        let naxis = header.require_int("NAXIS")? as usize;

        let mut axes = Vec::with_capacity(naxis);
        for i in 1..=naxis {
            let key = format!("NAXIS{i}");
            axes.push(header.require_int(&key)? as usize);
        }

        let pixels = PixelData::from_bytes(bitpix, data)?;

        Ok(ImageData { axes, pixels })
    }

    /// Populate header keywords for this image data.
    pub fn fill_header(&self, header: &mut Header) {
        header.set("BITPIX", HeaderValue::Integer(self.bitpix().to_i64()), None);
        header.set("NAXIS", HeaderValue::Integer(self.axes.len() as i64), None);
        for (i, &ax) in self.axes.iter().enumerate() {
            header.set(
                &format!("NAXIS{}", i + 1),
                HeaderValue::Integer(ax as i64),
                None,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pixel_data_round_trip_u8() {
        let data = vec![1u8, 2, 3, 4, 5];
        let pixels = PixelData::U8(data.clone());
        let bytes = pixels.to_bytes();
        let back = PixelData::from_bytes(Bitpix::U8, &bytes).unwrap();
        if let PixelData::U8(v) = back {
            assert_eq!(v, data);
        } else {
            panic!("wrong type");
        }
    }

    #[test]
    fn pixel_data_round_trip_i16() {
        let data = vec![-1000i16, 0, 1000, i16::MIN, i16::MAX];
        let pixels = PixelData::I16(data.clone());
        let bytes = pixels.to_bytes();
        let back = PixelData::from_bytes(Bitpix::I16, &bytes).unwrap();
        if let PixelData::I16(v) = back {
            assert_eq!(v, data);
        } else {
            panic!("wrong type");
        }
    }

    #[test]
    fn pixel_data_round_trip_f32() {
        let data = vec![1.5f32, -3.125, 0.0, f32::MAX];
        let pixels = PixelData::F32(data.clone());
        let bytes = pixels.to_bytes();
        let back = PixelData::from_bytes(Bitpix::F32, &bytes).unwrap();
        if let PixelData::F32(v) = back {
            assert_eq!(v, data);
        } else {
            panic!("wrong type");
        }
    }

    #[test]
    fn scaled_values() {
        let img = ImageData::new(vec![3], PixelData::I16(vec![0, 1, 2]));
        let scaled = img.scaled_values(2.0, 100.0);
        assert_eq!(scaled, vec![100.0, 102.0, 104.0]);
    }

    #[test]
    fn unsigned_u16_via_bzero() {
        // BITPIX=16, BZERO=32768 → unsigned u16
        let img = ImageData::new(vec![2], PixelData::I16(vec![-32768, 32767]));
        let scaled = img.scaled_values(1.0, 32768.0);
        assert_eq!(scaled, vec![0.0, 65535.0]);
    }
}
