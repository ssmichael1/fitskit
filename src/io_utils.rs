use crate::error::Result;
use crate::types::BLOCK_SIZE;
use std::io::{Read, Seek, SeekFrom, Write};

/// Read exactly `n` bytes of data, then skip padding to next block boundary.
pub fn read_data_block<R: Read>(reader: &mut R, n: usize) -> Result<Vec<u8>> {
    if n == 0 {
        return Ok(Vec::new());
    }

    let padded = padded_size(n);
    let mut buf = vec![0u8; padded];
    reader.read_exact(&mut buf)?;
    buf.truncate(n);
    Ok(buf)
}

/// Write data bytes, then pad with zeros to next block boundary.
pub fn write_data_block<W: Write>(writer: &mut W, data: &[u8]) -> Result<()> {
    if data.is_empty() {
        return Ok(());
    }

    writer.write_all(data)?;
    write_padding(writer, data.len(), 0)
}

/// Write the padding that follows `data_len` bytes of data so the total
/// reaches the next block boundary. `fill` is 0 for images and binary
/// tables, and ASCII blank (0x20) for ASCII tables, per the standard.
pub fn write_padding<W: Write>(writer: &mut W, data_len: usize, fill: u8) -> Result<()> {
    static ZEROS: [u8; BLOCK_SIZE] = [0u8; BLOCK_SIZE];
    static BLANKS: [u8; BLOCK_SIZE] = [b' '; BLOCK_SIZE];
    let padding = padded_size(data_len) - data_len;
    if padding != 0 {
        let block = if fill == b' ' { &BLANKS } else { &ZEROS };
        debug_assert!(fill == 0 || fill == b' ');
        writer.write_all(&block[..padding])?;
    }
    Ok(())
}

/// Skip the block padding that follows `data_len` bytes of data.
pub fn skip_padding<R: Read + Seek>(reader: &mut R, data_len: usize) -> Result<()> {
    let padding = padded_size(data_len) - data_len;
    if padding != 0 {
        reader.seek(SeekFrom::Current(padding as i64))?;
    }
    Ok(())
}

/// Skip over `n` data bytes plus padding to next block boundary.
pub fn skip_data_block<R: Read + Seek>(reader: &mut R, n: usize) -> Result<()> {
    if n == 0 {
        return Ok(());
    }
    let padded = padded_size(n);
    reader.seek(SeekFrom::Current(padded as i64))?;
    Ok(())
}

/// Pad data with zeros to the next BLOCK_SIZE boundary. Returns the data as-is
/// if already aligned, or an empty vec if input is empty.
pub fn pad_to_block(data: &[u8]) -> Vec<u8> {
    if data.is_empty() {
        return Vec::new();
    }
    let ps = padded_size(data.len());
    let mut out = Vec::with_capacity(ps);
    out.extend_from_slice(data);
    out.resize(ps, 0);
    out
}

/// Round up to the next multiple of BLOCK_SIZE.
pub fn padded_size(n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    let remainder = n % BLOCK_SIZE;
    if remainder == 0 {
        n
    } else {
        n + BLOCK_SIZE - remainder
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn padded_size_exact() {
        assert_eq!(padded_size(2880), 2880);
        assert_eq!(padded_size(5760), 5760);
    }

    #[test]
    fn padded_size_rounded() {
        assert_eq!(padded_size(1), 2880);
        assert_eq!(padded_size(2881), 5760);
    }

    #[test]
    fn padded_size_zero() {
        assert_eq!(padded_size(0), 0);
    }

    #[test]
    fn write_read_round_trip() {
        let data = vec![42u8; 1000];
        let mut buf = Vec::new();
        write_data_block(&mut buf, &data).unwrap();
        assert_eq!(buf.len(), 2880);

        let mut cursor = std::io::Cursor::new(&buf);
        let read_back = read_data_block(&mut cursor, 1000).unwrap();
        assert_eq!(read_back, data);
    }
}
