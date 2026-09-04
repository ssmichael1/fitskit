use crate::ascii_table::AsciiTable;
use crate::bintable::BinTable;
use crate::checksum;
use crate::error::{Error, Result};
use crate::header::Header;
use crate::image_data::ImageData;
use crate::io_utils;
use crate::keyword::HeaderValue;
use std::io::{Read, Seek, Write};

/// The data payload of an HDU.
#[derive(Debug, Clone)]
pub enum HduData {
    Empty,
    Image(ImageData),
    AsciiTable(AsciiTable),
    BinTable(BinTable),
}

/// A single Header Data Unit.
#[derive(Debug, Clone)]
pub struct Hdu {
    pub header: Header,
    pub data: HduData,
}

impl Hdu {
    pub fn new(header: Header, data: HduData) -> Self {
        Hdu { header, data }
    }

    /// Create a primary HDU with image data.
    pub fn primary_image(image: ImageData) -> Self {
        let mut header = Header::new();
        header.set(
            "SIMPLE",
            HeaderValue::Logical(true),
            Some("conforms to FITS standard"),
        );
        image.fill_header(&mut header);
        Hdu {
            header,
            data: HduData::Image(image),
        }
    }

    /// Create a primary HDU with no data.
    pub fn primary_empty() -> Self {
        let mut header = Header::new();
        header.set(
            "SIMPLE",
            HeaderValue::Logical(true),
            Some("conforms to FITS standard"),
        );
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(0), None);
        Hdu {
            header,
            data: HduData::Empty,
        }
    }

    /// Create an IMAGE extension HDU.
    pub fn image_extension(image: ImageData) -> Self {
        let mut header = Header::new();
        header.set(
            "XTENSION",
            HeaderValue::String("IMAGE".into()),
            Some("image extension"),
        );
        image.fill_header(&mut header);
        header.set("PCOUNT", HeaderValue::Integer(0), None);
        header.set("GCOUNT", HeaderValue::Integer(1), None);
        Hdu {
            header,
            data: HduData::Image(image),
        }
    }

    /// Create an ASCII TABLE extension HDU.
    pub fn ascii_table_extension(table: AsciiTable) -> Self {
        let mut header = Header::new();
        table.fill_header(&mut header);
        Hdu {
            header,
            data: HduData::AsciiTable(table),
        }
    }

    /// Create a BINTABLE extension HDU.
    pub fn bintable_extension(table: BinTable) -> Self {
        let mut header = Header::new();
        table.fill_header(&mut header);
        Hdu {
            header,
            data: HduData::BinTable(table),
        }
    }

    /// View this HDU as a tile-compressed image, if it is one.
    ///
    /// A tile-compressed image is stored on disk as a `BINTABLE` extension with
    /// `ZIMAGE = T` (the FITS Tiled Image Compression convention). This is a *cheap*
    /// detection step: it inspects the header and, when it matches, returns a
    /// [`CompressedImage`] view that borrows the underlying [`BinTable`] — preserving
    /// the original compressed tiles for lossless round-trip writing.
    ///
    /// Returns `None` when the HDU is not a compressed-image BINTABLE. The actual
    /// decoding happens in [`CompressedImage::decompress`]:
    ///
    /// ```no_run
    /// # use fitskit::FitsFile;
    /// # let fits = FitsFile::from_file("compressed.fits").unwrap();
    /// for hdu in fits.extensions() {
    ///     if let Some(cimg) = hdu.as_compressed_image() {
    ///         let image = cimg.decompress().unwrap();
    ///         println!("{:?}", image.axes);
    ///     }
    /// }
    /// ```
    pub fn as_compressed_image(&self) -> Option<crate::tile_compress::CompressedImage<'_>> {
        if !crate::tile_compress::CompressedImage::detect(&self.header) {
            return None;
        }
        match &self.data {
            HduData::BinTable(table) => {
                crate::tile_compress::CompressedImage::from_bintable(&self.header, table).ok()
            }
            _ => None,
        }
    }

    /// Read an HDU from a reader.
    pub fn read_from<R: Read + Seek>(reader: &mut R) -> Result<Self> {
        let header = Header::read_from(reader)?;
        let data_bytes = header.data_byte_count()?;

        // Determine HDU type
        let is_primary = header.find("SIMPLE").is_some();
        let xtension = header.get_string("XTENSION").map(|s| s.to_string());
        let is_image = is_primary || xtension.as_deref() == Some("IMAGE");

        if data_bytes == 0 {
            return Ok(Hdu {
                header,
                data: HduData::Empty,
            });
        }

        // For extensions the data size is |BITPIX|/8 * GCOUNT * (PCOUNT + product(NAXISn)).
        // For BINTABLE (BITPIX=8, GCOUNT=1) that is main table + heap; for a
        // plain image PCOUNT is 0.
        let pcount = header.get_int("PCOUNT").unwrap_or(0) as usize;

        if is_image && pcount == 0 {
            // Fast path: decode pixels straight from the reader (no raw copy),
            // then step over the block padding.
            let img = ImageData::read_from(&header, reader)?;
            io_utils::skip_padding(reader, data_bytes)?;
            return Ok(Hdu {
                header,
                data: HduData::Image(img),
            });
        }

        let raw = io_utils::read_data_block(reader, data_bytes)?;

        let data = if is_image {
            // Image data — exclude pcount bytes
            let img = ImageData::from_header_and_data(&header, &raw[..raw.len() - pcount])?;
            HduData::Image(img)
        } else if xtension.as_deref() == Some("TABLE") {
            let table = AsciiTable::from_header_and_vec(&header, raw)?;
            HduData::AsciiTable(table)
        } else if xtension.as_deref() == Some("BINTABLE") {
            let table = BinTable::from_header_and_vec(&header, raw)?;
            HduData::BinTable(table)
        } else if let Some(ext) = &xtension {
            return Err(Error::UnsupportedExtension(ext.clone()));
        } else {
            HduData::Empty
        };

        Ok(Hdu { header, data })
    }

    /// Write this HDU to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.write_impl(writer, false)
    }

    /// Write this HDU with CHECKSUM and DATASUM keywords computed and inserted.
    pub fn write_with_checksum<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.write_impl(writer, true)
    }

    /// Size of this HDU's data payload on disk, before block padding.
    pub fn data_byte_len(&self) -> usize {
        match &self.data {
            HduData::Empty => 0,
            HduData::Image(img) => img.pixels.byte_len(),
            HduData::AsciiTable(table) => table.raw_data.len(),
            HduData::BinTable(table) => table.main_data.len() + table.heap.len(),
        }
    }

    /// Byte value used to pad this HDU's data unit to a block boundary:
    /// ASCII blanks for ASCII tables, zeros otherwise (the FITS standard requires blank fill for ASCII tables).
    fn fill_byte(&self) -> u8 {
        match &self.data {
            HduData::AsciiTable(_) => b' ',
            _ => 0,
        }
    }

    /// Visit this HDU's on-disk (unpadded) data bytes in bounded chunks,
    /// without materializing a copy of the whole payload.
    fn for_each_data_chunk<E>(
        &self,
        mut f: impl FnMut(&[u8]) -> std::result::Result<(), E>,
    ) -> std::result::Result<(), E> {
        match &self.data {
            HduData::Empty => Ok(()),
            HduData::Image(img) => img.pixels.for_each_be_chunk(f),
            HduData::AsciiTable(table) => f(&table.raw_data),
            HduData::BinTable(table) => {
                f(&table.main_data)?;
                f(&table.heap)
            }
        }
    }

    /// Compute the DATASUM of this HDU's block-padded data unit.
    ///
    /// Zero padding contributes nothing to the sum, so only the payload is
    /// visited; the blank fill of ASCII tables is included explicitly.
    pub fn datasum(&self) -> u32 {
        let mut acc = checksum::Checksum::new();
        let _ = self.for_each_data_chunk(|chunk| {
            acc.update(chunk);
            Ok::<(), std::convert::Infallible>(())
        });
        let fill = self.fill_byte();
        if fill != 0 {
            let len = self.data_byte_len();
            let padding = io_utils::padded_size(len) - len;
            acc.update(&[fill; crate::types::BLOCK_SIZE][..padding]);
        }
        acc.finish()
    }

    fn write_impl<W: Write>(&self, writer: &mut W, with_checksum: bool) -> Result<()> {
        let mut header = self.header.clone();

        match &self.data {
            HduData::Empty => {}
            HduData::Image(img) => img.fill_header(&mut header),
            HduData::AsciiTable(table) => table.fill_header(&mut header),
            HduData::BinTable(table) => table.fill_header(&mut header),
        }

        if with_checksum {
            let header_bytes = checksum::stamp_hdu_with_datasum(&mut header, self.datasum())?;
            writer.write_all(&header_bytes)?;
        } else {
            header.write_to(writer)?;
        }

        self.for_each_data_chunk(|chunk| writer.write_all(chunk))?;
        io_utils::write_padding(writer, self.data_byte_len(), self.fill_byte())?;

        Ok(())
    }

    /// Verify the DATASUM of this HDU (if the keyword is present).
    pub fn verify_datasum(&self) -> Result<()> {
        if self.header.find("DATASUM").is_none() {
            return Ok(());
        }
        checksum::verify_datasum_value(&self.header, self.datasum())
    }
}
