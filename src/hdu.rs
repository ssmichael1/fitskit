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
    /// # use fits4::FitsFile;
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

        if data_bytes == 0 {
            // Skip padding if any
            io_utils::skip_data_block(reader, 0)?;
            return Ok(Hdu {
                header,
                data: HduData::Empty,
            });
        }

        // For bintable, we need to read main data + heap (PCOUNT bytes)
        let pcount = header.get_int("PCOUNT").unwrap_or(0) as usize;
        // data_byte_count already includes PCOUNT for extensions with the formula:
        // |BITPIX|/8 * GCOUNT * (PCOUNT + product(NAXISn))
        // But for BINTABLE, PCOUNT is heap size, separate from main table.
        // The standard formula: data_bytes = NAXIS1 * NAXIS2 (main) + PCOUNT (heap) for BINTABLE
        // Actually the standard formula is: |BITPIX|/8 * GCOUNT * (PCOUNT + NAXIS1*NAXIS2)
        // For BINTABLE: BITPIX=8, GCOUNT=1, so total = PCOUNT + NAXIS1*NAXIS2
        // data_byte_count() already computes this correctly.

        let raw = io_utils::read_data_block(reader, data_bytes)?;

        let data = if is_primary || xtension.as_deref() == Some("IMAGE") {
            // Image data — exclude pcount bytes (should be 0 for images)
            let naxis = header.get_int("NAXIS").unwrap_or(0) as usize;
            if naxis == 0 {
                HduData::Empty
            } else {
                let image_bytes = if pcount > 0 {
                    &raw[..raw.len() - pcount]
                } else {
                    &raw
                };
                let img = ImageData::from_header_and_data(&header, image_bytes)?;
                HduData::Image(img)
            }
        } else if xtension.as_deref() == Some("TABLE") {
            let table = AsciiTable::from_header_and_data(&header, &raw)?;
            HduData::AsciiTable(table)
        } else if xtension.as_deref() == Some("BINTABLE") {
            let table = BinTable::from_header_and_data(&header, &raw)?;
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

    fn write_impl<W: Write>(&self, writer: &mut W, with_checksum: bool) -> Result<()> {
        let mut header = self.header.clone();

        match &self.data {
            HduData::Empty => {}
            HduData::Image(img) => img.fill_header(&mut header),
            HduData::AsciiTable(table) => table.fill_header(&mut header),
            HduData::BinTable(table) => table.fill_header(&mut header),
        }

        let data_bytes = match &self.data {
            HduData::Empty => Vec::new(),
            HduData::Image(img) => img.pixels.to_bytes(),
            HduData::AsciiTable(table) => table.raw_data.clone(),
            HduData::BinTable(table) => table.to_bytes(),
        };

        let padded_data = io_utils::pad_to_block(&data_bytes);

        if with_checksum {
            let header_bytes = checksum::stamp_hdu(&mut header, &padded_data)?;
            writer.write_all(&header_bytes)?;
        } else {
            header.write_to(writer)?;
        }

        io_utils::write_data_block(writer, &data_bytes)?;

        Ok(())
    }

    /// Verify the DATASUM of this HDU (if the keyword is present).
    pub fn verify_datasum(&self) -> Result<()> {
        let data_bytes = match &self.data {
            HduData::Empty => Vec::new(),
            HduData::Image(img) => img.pixels.to_bytes(),
            HduData::AsciiTable(table) => table.raw_data.clone(),
            HduData::BinTable(table) => table.to_bytes(),
        };
        let padded = io_utils::pad_to_block(&data_bytes);
        checksum::verify_from_header(&self.header, &padded)
    }
}
