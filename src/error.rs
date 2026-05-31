use std::fmt;
use std::io;

/// All errors produced by this crate.
#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    InvalidFormat(String),
    InvalidKeyword(String),
    UnsupportedBitpix(i64),
    MissingKeyword(String),
    KeywordTypeMismatch(String),
    DataSizeMismatch { expected: usize, actual: usize },
    ChecksumMismatch { expected: u32, actual: u32 },
    UnsupportedExtension(String),
    InvalidTableFormat(String),
    UnsupportedCompression(String),
    CompressionError(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::InvalidFormat(s) => write!(f, "invalid FITS format: {s}"),
            Error::InvalidKeyword(s) => write!(f, "invalid keyword: {s}"),
            Error::UnsupportedBitpix(b) => write!(f, "unsupported BITPIX value: {b}"),
            Error::MissingKeyword(s) => write!(f, "missing required keyword: {s}"),
            Error::KeywordTypeMismatch(s) => write!(f, "keyword type mismatch: {s}"),
            Error::DataSizeMismatch { expected, actual } => {
                write!(f, "data size mismatch: expected {expected}, got {actual}")
            }
            Error::ChecksumMismatch { expected, actual } => {
                write!(
                    f,
                    "checksum mismatch: expected {expected:#010x}, got {actual:#010x}"
                )
            }
            Error::UnsupportedExtension(s) => write!(f, "unsupported extension: {s}"),
            Error::InvalidTableFormat(s) => write!(f, "invalid table format: {s}"),
            Error::UnsupportedCompression(s) => write!(f, "unsupported compression: {s}"),
            Error::CompressionError(s) => write!(f, "compression error: {s}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Error::Io(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;
