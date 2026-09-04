use crate::error::Result;
use crate::hdu::Hdu;
use crate::image_data::ImageData;
use std::fs::File;
use std::io::{BufReader, BufWriter, Cursor, Read, Seek, Write};
use std::path::Path;

/// Top-level FITS file containing one or more HDUs.
#[derive(Debug, Clone)]
pub struct FitsFile {
    pub hdus: Vec<Hdu>,
}

impl FitsFile {
    pub fn new() -> Self {
        FitsFile { hdus: Vec::new() }
    }

    /// Read a FITS file from a seekable reader.
    pub fn from_reader<R: Read + Seek>(reader: &mut R) -> Result<Self> {
        let mut hdus = Vec::new();

        // Read primary HDU
        hdus.push(Hdu::read_from(reader)?);

        // Read extension HDUs until EOF
        loop {
            match Hdu::read_from(reader) {
                Ok(hdu) => hdus.push(hdu),
                Err(crate::error::Error::Io(ref e))
                    if e.kind() == std::io::ErrorKind::UnexpectedEof =>
                {
                    break;
                }
                Err(e) => return Err(e),
            }
        }

        Ok(FitsFile { hdus })
    }

    /// Read a FITS file from a file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::from_reader(&mut reader)
    }

    /// Read a FITS file from a byte slice.
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut cursor = Cursor::new(data);
        Self::from_reader(&mut cursor)
    }

    /// Write to a writer.
    pub fn to_writer<W: Write>(&self, writer: &mut W) -> Result<()> {
        for hdu in &self.hdus {
            hdu.write_to(writer)?;
        }
        Ok(())
    }

    /// Write to a writer with CHECKSUM/DATASUM keywords in every HDU.
    pub fn to_writer_with_checksum<W: Write>(&self, writer: &mut W) -> Result<()> {
        for hdu in &self.hdus {
            hdu.write_with_checksum(writer)?;
        }
        Ok(())
    }

    /// Write to a byte vector with CHECKSUM/DATASUM keywords.
    pub fn to_bytes_with_checksum(&self) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(self.estimated_byte_len());
        self.to_writer_with_checksum(&mut buf)?;
        Ok(buf)
    }

    /// Write to a file path.
    pub fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        self.to_writer(&mut writer)
    }

    /// Write to a byte vector.
    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(self.estimated_byte_len());
        self.to_writer(&mut buf)?;
        Ok(buf)
    }

    /// Estimate of the serialized size for pre-allocation: padded data units
    /// plus header blocks sized from the keyword count (with slack for the
    /// mandatory keywords `write` may add). Long-string CONTINUE cards can push
    /// a header past this, in which case the buffer simply grows once.
    fn estimated_byte_len(&self) -> usize {
        use crate::types::{BLOCK_SIZE, RECORDS_PER_BLOCK};
        self.hdus
            .iter()
            .map(|h| {
                let cards = h.header.keywords.len() + 8;
                crate::io_utils::padded_size(h.data_byte_len())
                    + cards.div_ceil(RECORDS_PER_BLOCK) * BLOCK_SIZE
            })
            .sum()
    }

    /// Get the primary HDU.
    pub fn primary(&self) -> &Hdu {
        &self.hdus[0]
    }

    /// Get the primary HDU mutably.
    pub fn primary_mut(&mut self) -> &mut Hdu {
        &mut self.hdus[0]
    }

    /// Get extension HDUs (everything after the primary).
    pub fn extensions(&self) -> &[Hdu] {
        if self.hdus.len() > 1 {
            &self.hdus[1..]
        } else {
            &[]
        }
    }

    /// Number of HDUs.
    pub fn len(&self) -> usize {
        self.hdus.len()
    }

    /// Whether the file has no HDUs.
    pub fn is_empty(&self) -> bool {
        self.hdus.is_empty()
    }

    /// Create a FitsFile with a single primary image HDU.
    pub fn with_primary_image(image: ImageData) -> Self {
        FitsFile {
            hdus: vec![Hdu::primary_image(image)],
        }
    }

    /// Create a FitsFile with an empty primary HDU.
    pub fn with_empty_primary() -> Self {
        FitsFile {
            hdus: vec![Hdu::primary_empty()],
        }
    }

    /// Add an extension HDU.
    pub fn push_extension(&mut self, hdu: Hdu) {
        self.hdus.push(hdu);
    }

    /// Find the first extension HDU with a given EXTNAME.
    pub fn find_extension(&self, extname: &str) -> Option<&Hdu> {
        self.extensions()
            .iter()
            .find(|hdu| hdu.header.get_string("EXTNAME") == Some(extname))
    }

    /// Iterate over all HDUs.
    pub fn iter(&self) -> std::slice::Iter<'_, Hdu> {
        self.hdus.iter()
    }
}

impl Default for FitsFile {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> IntoIterator for &'a FitsFile {
    type Item = &'a Hdu;
    type IntoIter = std::slice::Iter<'a, Hdu>;

    fn into_iter(self) -> Self::IntoIter {
        self.hdus.iter()
    }
}
