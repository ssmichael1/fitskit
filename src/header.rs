use crate::error::{Error, Result};
use crate::keyword::{HeaderValue, Keyword};
use crate::types::{BLOCK_SIZE, RECORDS_PER_BLOCK, RECORD_SIZE};
use std::io::{Read, Write};

/// An ordered collection of FITS header keywords.
#[derive(Debug, Clone)]
pub struct Header {
    pub keywords: Vec<Keyword>,
}

impl Header {
    pub fn new() -> Self {
        Header {
            keywords: Vec::new(),
        }
    }

    /// Find the first keyword with the given name.
    pub fn find(&self, name: &str) -> Option<&Keyword> {
        self.keywords
            .iter()
            .find(|k| k.name.eq_ignore_ascii_case(name))
    }

    /// Get an integer value by keyword name.
    pub fn get_int(&self, name: &str) -> Option<i64> {
        self.find(name)?.value.as_ref()?.as_int()
    }

    /// Get a float value by keyword name.
    pub fn get_float(&self, name: &str) -> Option<f64> {
        self.find(name)?.value.as_ref()?.as_float()
    }

    /// Get a string value by keyword name.
    pub fn get_string(&self, name: &str) -> Option<&str> {
        self.find(name)?.value.as_ref()?.as_str()
    }

    /// Get a boolean value by keyword name.
    pub fn get_bool(&self, name: &str) -> Option<bool> {
        self.find(name)?.value.as_ref()?.as_bool()
    }

    /// Get a required integer, returning error if missing.
    pub fn require_int(&self, name: &str) -> Result<i64> {
        self.get_int(name)
            .ok_or_else(|| Error::MissingKeyword(name.into()))
    }

    /// Set or update a keyword. If name exists, update in place; otherwise append.
    pub fn set(&mut self, name: &str, value: HeaderValue, comment: Option<&str>) {
        if let Some(kw) = self
            .keywords
            .iter_mut()
            .find(|k| k.name.eq_ignore_ascii_case(name))
        {
            kw.value = Some(value);
            if let Some(c) = comment {
                kw.comment = Some(c.to_string());
            }
        } else {
            self.keywords
                .push(Keyword::with_value(name, value, comment));
        }
    }

    /// Add a keyword without checking for duplicates.
    pub fn push(&mut self, keyword: Keyword) {
        self.keywords.push(keyword);
    }

    /// Iterate over keywords.
    pub fn iter(&self) -> std::slice::Iter<'_, Keyword> {
        self.keywords.iter()
    }

    /// Read a header from a block-aligned reader.
    /// Handles CONTINUE long-string assembly.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let mut keywords = Vec::new();
        let mut buf = [0u8; BLOCK_SIZE];

        'outer: loop {
            reader.read_exact(&mut buf).map_err(Error::Io)?;

            for i in 0..RECORDS_PER_BLOCK {
                let start = i * RECORD_SIZE;
                let record: &[u8; RECORD_SIZE] =
                    buf[start..start + RECORD_SIZE].try_into().unwrap();

                // Check for END keyword
                if record[..3] == *b"END" && record[3..8].iter().all(|&b| b == b' ') {
                    break 'outer;
                }

                let kw = Keyword::parse(record)?;

                // Handle CONTINUE - append to previous string value
                if kw.name == "CONTINUE" && !keywords.is_empty() {
                    let prev: &mut Keyword = keywords.last_mut().unwrap();
                    if let Some(HeaderValue::String(ref mut s)) = prev.value {
                        if let Some(HeaderValue::String(ref cont)) = kw.value {
                            s.push_str(cont);
                        }
                        // Update comment if the CONTINUE card has one
                        if kw.comment.is_some() {
                            prev.comment = kw.comment;
                        }
                    }
                } else {
                    keywords.push(kw);
                }
            }
        }

        Ok(Header { keywords })
    }

    /// Serialize the header to block-aligned bytes (a multiple of 2880).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out =
            Vec::with_capacity((self.keywords.len() + 1).div_ceil(RECORDS_PER_BLOCK) * BLOCK_SIZE);

        for kw in &self.keywords {
            for card in kw.to_cards() {
                out.extend_from_slice(&card);
            }
        }

        // Append END card
        let mut end_card = [b' '; RECORD_SIZE];
        end_card[..3].copy_from_slice(b"END");
        out.extend_from_slice(&end_card);

        // Pad with blank cards to fill the last block
        let padded = out.len().div_ceil(BLOCK_SIZE) * BLOCK_SIZE;
        out.resize(padded, b' ');
        out
    }

    /// Write the header as block-aligned 2880-byte blocks.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(&self.to_bytes())?;
        Ok(())
    }

    /// Compute the number of data bytes described by this header.
    /// For primary/image: |BITPIX|/8 * product(NAXISn)
    /// For extensions: |BITPIX|/8 * GCOUNT * (PCOUNT + product(NAXISn))
    pub fn data_byte_count(&self) -> Result<usize> {
        let naxis = self.require_int("NAXIS")? as usize;
        if naxis == 0 {
            return Ok(0);
        }

        let bitpix = self.require_int("BITPIX")?;
        let bytes_per_val = (bitpix.unsigned_abs() / 8) as usize;

        let mut product: usize = 1;
        for i in 1..=naxis {
            let key = format!("NAXIS{i}");
            let n = self.require_int(&key)? as usize;
            product = product
                .checked_mul(n)
                .ok_or_else(|| Error::InvalidFormat("axis dimensions overflow".into()))?;
        }

        let pcount = self.get_int("PCOUNT").unwrap_or(0) as usize;
        let gcount = self.get_int("GCOUNT").unwrap_or(1) as usize;

        let total = bytes_per_val
            .checked_mul(gcount)
            .and_then(|v| v.checked_mul(pcount + product))
            .ok_or_else(|| Error::InvalidFormat("data size overflow".into()))?;

        Ok(total)
    }
}

impl Default for Header {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_accessors() {
        let mut h = Header::new();
        h.set("BITPIX", HeaderValue::Integer(16), None);
        h.set("BSCALE", HeaderValue::Float(1.0), None);
        h.set("OBJECT", HeaderValue::String("M31".into()), None);
        h.set("SIMPLE", HeaderValue::Logical(true), None);

        assert_eq!(h.get_int("BITPIX"), Some(16));
        assert_eq!(h.get_float("BSCALE"), Some(1.0));
        assert_eq!(h.get_string("OBJECT"), Some("M31"));
        assert_eq!(h.get_bool("SIMPLE"), Some(true));
    }

    #[test]
    fn set_updates_existing() {
        let mut h = Header::new();
        h.set("NAXIS", HeaderValue::Integer(0), None);
        h.set("NAXIS", HeaderValue::Integer(2), None);
        assert_eq!(h.get_int("NAXIS"), Some(2));
        assert_eq!(h.keywords.len(), 1);
    }

    #[test]
    fn data_byte_count_image() {
        let mut h = Header::new();
        h.set("BITPIX", HeaderValue::Integer(16), None);
        h.set("NAXIS", HeaderValue::Integer(2), None);
        h.set("NAXIS1", HeaderValue::Integer(100), None);
        h.set("NAXIS2", HeaderValue::Integer(200), None);
        assert_eq!(h.data_byte_count().unwrap(), 2 * 100 * 200);
    }

    #[test]
    fn data_byte_count_zero_naxis() {
        let mut h = Header::new();
        h.set("BITPIX", HeaderValue::Integer(8), None);
        h.set("NAXIS", HeaderValue::Integer(0), None);
        assert_eq!(h.data_byte_count().unwrap(), 0);
    }

    #[test]
    fn header_round_trip() {
        let mut h = Header::new();
        h.set(
            "SIMPLE",
            HeaderValue::Logical(true),
            Some("conforms to standard"),
        );
        h.set("BITPIX", HeaderValue::Integer(16), None);
        h.set("NAXIS", HeaderValue::Integer(0), None);

        let mut buf = Vec::new();
        h.write_to(&mut buf).unwrap();

        assert_eq!(buf.len() % BLOCK_SIZE, 0);

        let mut cursor = std::io::Cursor::new(&buf);
        let h2 = Header::read_from(&mut cursor).unwrap();

        assert_eq!(h2.get_bool("SIMPLE"), Some(true));
        assert_eq!(h2.get_int("BITPIX"), Some(16));
        assert_eq!(h2.get_int("NAXIS"), Some(0));
    }
}
