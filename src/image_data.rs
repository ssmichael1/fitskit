use crate::error::{Error, Result};
use crate::header::Header;
use crate::keyword::HeaderValue;
use crate::types::Bitpix;
use std::io::{Read, Write};

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

/// Chunk size (bytes) for streaming conversion between native and big-endian
/// pixel storage. A multiple of 8 so no value straddles a chunk boundary, and
/// of 4 so checksums can be accumulated chunk by chunk. Large enough that
/// file I/O happens in few syscalls.
const IO_CHUNK: usize = 1 << 20;

/// Fixed-width big-endian (de)serialization for pixel element types.
trait BigEndian: Copy {
    const SIZE: usize;
    /// Decode `bytes` (length a multiple of `SIZE`) and append to `out`.
    fn decode_be(bytes: &[u8], out: &mut Vec<Self>);
    /// Encode `vals` into `out` (length exactly `vals.len() * SIZE`).
    fn encode_be(vals: &[Self], out: &mut [u8]);
}

macro_rules! impl_big_endian {
    ($($t:ty => $n:literal),* $(,)?) => {$(
        impl BigEndian for $t {
            const SIZE: usize = $n;

            #[inline]
            fn decode_be(bytes: &[u8], out: &mut Vec<Self>) {
                let (words, rem) = bytes.as_chunks::<$n>();
                debug_assert!(rem.is_empty());
                out.extend(words.iter().map(|&w| <$t>::from_be_bytes(w)));
            }

            #[inline]
            fn encode_be(vals: &[Self], out: &mut [u8]) {
                let (words, rem) = out.as_chunks_mut::<$n>();
                debug_assert!(rem.is_empty() && words.len() == vals.len());
                for (dst, &v) in words.iter_mut().zip(vals) {
                    *dst = v.to_be_bytes();
                }
            }
        }
    )*};
}

impl_big_endian!(i16 => 2, i32 => 4, i64 => 8, f32 => 4, f64 => 8);

/// Decode a whole big-endian byte slice into a new vector.
fn decode_vec<T: BigEndian>(bytes: &[u8]) -> Vec<T> {
    let mut out = Vec::with_capacity(bytes.len() / T::SIZE);
    T::decode_be(bytes, &mut out);
    out
}

/// Read `n` big-endian values from `reader`, decoding through a bounded
/// scratch buffer so peak memory is the output vector plus one chunk.
fn read_vec<T: BigEndian, R: Read>(reader: &mut R, n: usize) -> Result<Vec<T>> {
    let total = n * T::SIZE;
    let mut out = Vec::with_capacity(n);
    let mut scratch = vec![0u8; IO_CHUNK.min(total)];
    let mut remaining = total;
    while remaining > 0 {
        let take = remaining.min(scratch.len());
        reader.read_exact(&mut scratch[..take])?;
        T::decode_be(&scratch[..take], &mut out);
        remaining -= take;
    }
    Ok(out)
}

/// Encode `vals` as big-endian bytes, handing each chunk to `f`.
fn for_each_be_chunk<T: BigEndian, E>(
    vals: &[T],
    mut f: impl FnMut(&[u8]) -> std::result::Result<(), E>,
) -> std::result::Result<(), E> {
    let mut scratch = vec![0u8; IO_CHUNK.min(vals.len() * T::SIZE)];
    for chunk in vals.chunks(IO_CHUNK / T::SIZE) {
        let buf = &mut scratch[..chunk.len() * T::SIZE];
        T::encode_be(chunk, buf);
        f(buf)?;
    }
    Ok(())
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

    /// Size of the big-endian on-disk encoding in bytes (unpadded).
    pub fn byte_len(&self) -> usize {
        self.len() * self.bitpix().bytes_per_value()
    }

    /// Convert pixel data to big-endian bytes.
    ///
    /// (The `flat_map` form compiles to a tight vectorized loop and measures
    /// slightly faster than encoding into a pre-zeroed buffer.)
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

    /// Visit the big-endian on-disk encoding in bounded chunks without
    /// materializing the whole byte buffer. Every chunk but the last is a
    /// multiple of 8 bytes.
    pub(crate) fn for_each_be_chunk<E>(
        &self,
        mut f: impl FnMut(&[u8]) -> std::result::Result<(), E>,
    ) -> std::result::Result<(), E> {
        match self {
            PixelData::U8(v) => {
                if v.is_empty() {
                    Ok(())
                } else {
                    f(v)
                }
            }
            PixelData::I16(v) => for_each_be_chunk(v, f),
            PixelData::I32(v) => for_each_be_chunk(v, f),
            PixelData::I64(v) => for_each_be_chunk(v, f),
            PixelData::F32(v) => for_each_be_chunk(v, f),
            PixelData::F64(v) => for_each_be_chunk(v, f),
        }
    }

    /// Write the big-endian encoding to `writer` (unpadded).
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.for_each_be_chunk(|chunk| writer.write_all(chunk))
            .map_err(Error::Io)
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
            Bitpix::I16 => PixelData::I16(decode_vec(data)),
            Bitpix::I32 => PixelData::I32(decode_vec(data)),
            Bitpix::I64 => PixelData::I64(decode_vec(data)),
            Bitpix::F32 => PixelData::F32(decode_vec(data)),
            Bitpix::F64 => PixelData::F64(decode_vec(data)),
        })
    }

    /// Read `n` big-endian values of type `bitpix` from `reader`.
    ///
    /// Streams through a bounded scratch buffer, so peak memory is the decoded
    /// pixel vector plus about 1 MiB rather than raw bytes plus pixels.
    pub fn read_from<R: Read>(reader: &mut R, bitpix: Bitpix, n: usize) -> Result<Self> {
        Ok(match bitpix {
            Bitpix::U8 => {
                let mut v = vec![0u8; n];
                reader.read_exact(&mut v)?;
                PixelData::U8(v)
            }
            Bitpix::I16 => PixelData::I16(read_vec(reader, n)?),
            Bitpix::I32 => PixelData::I32(read_vec(reader, n)?),
            Bitpix::I64 => PixelData::I64(read_vec(reader, n)?),
            Bitpix::F32 => PixelData::F32(read_vec(reader, n)?),
            Bitpix::F64 => PixelData::F64(read_vec(reader, n)?),
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

    /// Parse BITPIX and NAXISn from a header.
    fn shape_from_header(header: &Header) -> Result<(Bitpix, Vec<usize>)> {
        let bitpix = Bitpix::from_i64(header.require_int("BITPIX")?)?;
        let naxis = header.require_int("NAXIS")? as usize;

        let mut axes = Vec::with_capacity(naxis);
        for i in 1..=naxis {
            let key = format!("NAXIS{i}");
            axes.push(header.require_int(&key)? as usize);
        }
        Ok((bitpix, axes))
    }

    /// Read image data from a header and raw bytes.
    pub fn from_header_and_data(header: &Header, data: &[u8]) -> Result<Self> {
        let (bitpix, axes) = Self::shape_from_header(header)?;
        let pixels = PixelData::from_bytes(bitpix, data)?;
        Ok(ImageData { axes, pixels })
    }

    /// Read image data described by `header` directly from `reader`
    /// (unpadded; the caller skips the block padding). Decodes in a
    /// streaming fashion, avoiding an intermediate copy of the raw bytes.
    pub fn read_from<R: Read>(header: &Header, reader: &mut R) -> Result<Self> {
        let (bitpix, axes) = Self::shape_from_header(header)?;
        let n = axes.iter().product();
        let pixels = PixelData::read_from(reader, bitpix, n)?;
        Ok(ImageData { axes, pixels })
    }

    /// Compress this image into a tile-compressed BINTABLE [`Hdu`](crate::hdu::Hdu)
    /// (`ZIMAGE = T`), ready for [`FitsFile::push_extension`](crate::fits::FitsFile::push_extension).
    ///
    /// Thin wrapper over [`compress_image`](crate::tile_compress::compress_image). See
    /// [`CompressOptions`](crate::tile_compress::CompressOptions) for the algorithm,
    /// tiling, and float-quantization knobs. The round-trip inverse is
    /// [`Hdu::as_compressed_image`](crate::hdu::Hdu::as_compressed_image) +
    /// [`CompressedImage::decompress`](crate::tile_compress::CompressedImage::decompress).
    ///
    /// ```
    /// use fitskit::{ImageData, PixelData};
    /// use fitskit::tile_compress::CompressOptions;
    ///
    /// let img = ImageData::new(vec![8, 4], PixelData::I16((0..32).collect()));
    /// let hdu = img.compress(&CompressOptions::default()).unwrap();
    /// let back = hdu.as_compressed_image().unwrap().decompress().unwrap();
    /// assert_eq!(back.pixels.to_bytes(), img.pixels.to_bytes());
    /// ```
    pub fn compress(
        &self,
        opts: &crate::tile_compress::CompressOptions,
    ) -> Result<crate::hdu::Hdu> {
        crate::tile_compress::compress_image(self, opts)
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
    fn stream_round_trip_all_types_across_chunks() {
        // More than one IO_CHUNK of bytes for every element size, so the
        // chunk loops take several iterations.
        let n = IO_CHUNK / 2 + 123;
        let cases = vec![
            PixelData::U8((0..n).map(|i| i as u8).collect()),
            PixelData::I16((0..n).map(|i| i as i16).collect()),
            PixelData::I32((0..n).map(|i| (i as i32).wrapping_mul(7919)).collect()),
            PixelData::I64(
                (0..n)
                    .map(|i| (i as i64).wrapping_mul(-1_000_003))
                    .collect(),
            ),
            PixelData::F32((0..n).map(|i| i as f32 * 0.5).collect()),
            PixelData::F64((0..n).map(|i| i as f64 * -0.25).collect()),
        ];
        for px in cases {
            let direct = px.to_bytes();
            let mut streamed = Vec::new();
            px.write_to(&mut streamed).unwrap();
            assert_eq!(streamed, direct);
            assert_eq!(streamed.len(), px.byte_len());

            let mut cursor = std::io::Cursor::new(&streamed);
            let back = PixelData::read_from(&mut cursor, px.bitpix(), n).unwrap();
            assert_eq!(back.to_bytes(), direct);
            let back2 = PixelData::from_bytes(px.bitpix(), &streamed).unwrap();
            assert_eq!(back2.to_bytes(), direct);
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
