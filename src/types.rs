use crate::error::{Error, Result};

/// Size of one FITS block in bytes.
pub const BLOCK_SIZE: usize = 2880;

/// Size of one header record (keyword card) in bytes.
pub const RECORD_SIZE: usize = 80;

/// Number of records per header block.
pub const RECORDS_PER_BLOCK: usize = BLOCK_SIZE / RECORD_SIZE;

/// BITPIX values per FITS standard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bitpix {
    U8,
    I16,
    I32,
    I64,
    F32,
    F64,
}

impl Bitpix {
    pub fn to_i64(self) -> i64 {
        match self {
            Bitpix::U8 => 8,
            Bitpix::I16 => 16,
            Bitpix::I32 => 32,
            Bitpix::I64 => 64,
            Bitpix::F32 => -32,
            Bitpix::F64 => -64,
        }
    }

    pub fn from_i64(val: i64) -> Result<Self> {
        match val {
            8 => Ok(Bitpix::U8),
            16 => Ok(Bitpix::I16),
            32 => Ok(Bitpix::I32),
            64 => Ok(Bitpix::I64),
            -32 => Ok(Bitpix::F32),
            -64 => Ok(Bitpix::F64),
            _ => Err(Error::UnsupportedBitpix(val)),
        }
    }

    pub fn bytes_per_value(self) -> usize {
        match self {
            Bitpix::U8 => 1,
            Bitpix::I16 => 2,
            Bitpix::I32 => 4,
            Bitpix::I64 => 8,
            Bitpix::F32 => 4,
            Bitpix::F64 => 8,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bitpix_round_trip() {
        for bp in [
            Bitpix::U8,
            Bitpix::I16,
            Bitpix::I32,
            Bitpix::I64,
            Bitpix::F32,
            Bitpix::F64,
        ] {
            assert_eq!(Bitpix::from_i64(bp.to_i64()).unwrap(), bp);
        }
    }

    #[test]
    fn bitpix_invalid() {
        assert!(Bitpix::from_i64(0).is_err());
        assert!(Bitpix::from_i64(128).is_err());
    }

    #[test]
    fn bytes_per_value() {
        assert_eq!(Bitpix::U8.bytes_per_value(), 1);
        assert_eq!(Bitpix::I16.bytes_per_value(), 2);
        assert_eq!(Bitpix::F64.bytes_per_value(), 8);
    }
}
