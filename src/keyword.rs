use crate::error::{Error, Result};
use crate::types::RECORD_SIZE;
use std::fmt;

/// A value stored in a FITS header keyword.
#[derive(Debug, Clone, PartialEq)]
pub enum HeaderValue {
    Logical(bool),
    Integer(i64),
    Float(f64),
    String(String),
    ComplexInteger(i64, i64),
    ComplexFloat(f64, f64),
    Undefined,
}

impl HeaderValue {
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            HeaderValue::Logical(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_int(&self) -> Option<i64> {
        match self {
            HeaderValue::Integer(i) => Some(*i),
            _ => None,
        }
    }

    pub fn as_float(&self) -> Option<f64> {
        match self {
            HeaderValue::Float(f) => Some(*f),
            HeaderValue::Integer(i) => Some(*i as f64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            HeaderValue::String(s) => Some(s),
            _ => None,
        }
    }
}

impl fmt::Display for HeaderValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HeaderValue::Logical(b) => write!(f, "{}", if *b { "T" } else { "F" }),
            HeaderValue::Integer(i) => write!(f, "{i}"),
            HeaderValue::Float(v) => write!(f, "{v}"),
            HeaderValue::String(s) => write!(f, "'{s}'"),
            HeaderValue::ComplexInteger(r, i) => write!(f, "({r}, {i})"),
            HeaderValue::ComplexFloat(r, i) => write!(f, "({r}, {i})"),
            HeaderValue::Undefined => write!(f, ""),
        }
    }
}

/// A single FITS header keyword (card image).
#[derive(Debug, Clone)]
pub struct Keyword {
    pub name: String,
    pub value: Option<HeaderValue>,
    pub comment: Option<String>,
}

impl Keyword {
    pub fn new(name: &str, value: Option<HeaderValue>, comment: Option<&str>) -> Self {
        Keyword {
            name: name.to_uppercase(),
            value,
            comment: comment.map(|s| s.to_string()),
        }
    }

    /// Create a valued keyword (has `= ` indicator).
    pub fn with_value(name: &str, value: HeaderValue, comment: Option<&str>) -> Self {
        Keyword::new(name, Some(value), comment)
    }

    /// Create a commentary keyword (COMMENT, HISTORY, or blank).
    pub fn commentary(name: &str, text: &str) -> Self {
        Keyword {
            name: name.to_uppercase(),
            value: None,
            comment: Some(text.to_string()),
        }
    }

    /// Parse a single 80-byte card image into a Keyword.
    pub fn parse(record: &[u8; RECORD_SIZE]) -> Result<Self> {
        let card = std::str::from_utf8(record)
            .map_err(|_| Error::InvalidKeyword("non-ASCII card image".into()))?;

        let name = card[..8].trim_end().to_string();

        // Check for value indicator "= "
        if card.len() >= 10 && &card[8..10] == "= " {
            let value_comment = &card[10..];
            let (value, comment) = parse_value_comment(value_comment)?;
            Ok(Keyword {
                name,
                value: Some(value),
                comment,
            })
        } else {
            // Commentary or blank keyword
            let text = card[8..].trim_end();
            let comment = if text.is_empty() {
                None
            } else {
                Some(text.to_string())
            };
            Ok(Keyword {
                name,
                value: None,
                comment,
            })
        }
    }

    /// Serialize this keyword into one or more 80-byte card images.
    /// Returns multiple cards if CONTINUE is needed for long strings.
    pub fn to_cards(&self) -> Vec<[u8; RECORD_SIZE]> {
        let mut cards = Vec::new();

        if self.name == "END" {
            let mut card = [b' '; RECORD_SIZE];
            card[..3].copy_from_slice(b"END");
            cards.push(card);
            return cards;
        }

        if self.value.is_none() {
            // Commentary keyword
            let mut card = [b' '; RECORD_SIZE];
            let name_bytes = self.name.as_bytes();
            let len = name_bytes.len().min(8);
            card[..len].copy_from_slice(&name_bytes[..len]);
            if let Some(ref comment) = self.comment {
                let cmt = comment.as_bytes();
                let clen = cmt.len().min(RECORD_SIZE - 8);
                card[8..8 + clen].copy_from_slice(&cmt[..clen]);
            }
            cards.push(card);
            return cards;
        }

        let value = self.value.as_ref().unwrap();

        // For string values, handle potential CONTINUE
        if let HeaderValue::String(ref s) = value {
            let cards_out = serialize_string_keyword(&self.name, s, self.comment.as_deref());
            return cards_out;
        }

        // Non-string valued keyword
        let mut card = [b' '; RECORD_SIZE];
        let name_bytes = self.name.as_bytes();
        let len = name_bytes.len().min(8);
        card[..len].copy_from_slice(&name_bytes[..len]);
        card[8] = b'=';
        card[9] = b' ';

        let val_str = format_value(value);
        let val_bytes = val_str.as_bytes();

        if let Some(ref comment) = self.comment {
            // Value right-justified in columns 11-30, then " / comment"
            let val_start = if val_bytes.len() < 20 {
                30 - val_bytes.len()
            } else {
                10
            };
            card[val_start..val_start + val_bytes.len()].copy_from_slice(val_bytes);

            let cmt_start = val_start + val_bytes.len() + 1;
            if cmt_start + 2 < RECORD_SIZE {
                card[cmt_start] = b'/';
                card[cmt_start + 1] = b' ';
                let cmt = comment.as_bytes();
                let avail = RECORD_SIZE - cmt_start - 2;
                let clen = cmt.len().min(avail);
                card[cmt_start + 2..cmt_start + 2 + clen].copy_from_slice(&cmt[..clen]);
            }
        } else {
            // Right-justify in columns 11-30 for numeric/logical
            if val_bytes.len() <= 20 {
                let start = 30 - val_bytes.len();
                card[start..30].copy_from_slice(val_bytes);
            } else {
                card[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
            }
        }

        cards.push(card);
        cards
    }
}

fn format_value(value: &HeaderValue) -> String {
    match value {
        HeaderValue::Logical(b) => {
            if *b {
                "T".to_string()
            } else {
                "F".to_string()
            }
        }
        HeaderValue::Integer(i) => format!("{i}"),
        HeaderValue::Float(f) => format_float(*f),
        HeaderValue::String(s) => format!("'{}'", pad_string_value(s)),
        HeaderValue::ComplexInteger(r, i) => format!("({r}, {i})"),
        HeaderValue::ComplexFloat(r, i) => format!("({}, {})", format_float(*r), format_float(*i)),
        HeaderValue::Undefined => String::new(),
    }
}

fn format_float(f: f64) -> String {
    // Use scientific notation with enough precision
    let s = format!("{:.15E}", f);
    s
}

fn pad_string_value(s: &str) -> String {
    // String values must be at least 8 characters inside the quotes
    if s.len() < 8 {
        format!("{:<8}", s)
    } else {
        s.to_string()
    }
}

fn serialize_string_keyword(
    name: &str,
    value: &str,
    comment: Option<&str>,
) -> Vec<[u8; RECORD_SIZE]> {
    let mut cards = Vec::new();

    // Escape single quotes by doubling them
    let escaped = value.replace('\'', "''");

    // Max string content in first card: columns 11-80 = 70 chars
    // Format: '<string padded to >=8>' or with & for continuation
    // Need room for opening quote, content, closing quote
    // First card: 70 chars available after "= "
    // Need: ' + content + ' + optional_comment = 70 chars max
    // Available content space in first card (cols 11-80 = 70 bytes)
    // '...' takes at minimum '        ' (8 chars) + 2 quotes = 10
    // With comment: '<content>' / comment
    let first_avail = if comment.is_some() { 55 } else { 67 }; // conservative

    if escaped.len() <= first_avail {
        // Fits in one card
        let mut card = [b' '; RECORD_SIZE];
        let name_bytes = name.as_bytes();
        let len = name_bytes.len().min(8);
        card[..len].copy_from_slice(&name_bytes[..len]);
        card[8] = b'=';
        card[9] = b' ';

        let padded = pad_string_value(&escaped);
        let val_str = format!("'{padded}'");
        let val_bytes = val_str.as_bytes();
        card[10..10 + val_bytes.len()].copy_from_slice(val_bytes);

        if let Some(cmt) = comment {
            let cmt_start = 10 + val_bytes.len() + 1;
            if cmt_start + 2 < RECORD_SIZE {
                card[cmt_start] = b'/';
                card[cmt_start + 1] = b' ';
                let cmt_bytes = cmt.as_bytes();
                let avail = RECORD_SIZE - cmt_start - 2;
                let clen = cmt_bytes.len().min(avail);
                card[cmt_start + 2..cmt_start + 2 + clen].copy_from_slice(&cmt_bytes[..clen]);
            }
        }

        cards.push(card);
    } else {
        // Need CONTINUE cards
        let mut remaining = escaped.as_str();

        // First card: can fit 67 chars of string content (70 - quote - quote - &)
        let first_chunk_len = 67.min(remaining.len());
        let chunk = &remaining[..first_chunk_len];
        remaining = &remaining[first_chunk_len..];

        let mut card = [b' '; RECORD_SIZE];
        let name_bytes = name.as_bytes();
        let nlen = name_bytes.len().min(8);
        card[..nlen].copy_from_slice(&name_bytes[..nlen]);
        card[8] = b'=';
        card[9] = b' ';
        let val_str = format!("'{chunk}&'");
        let val_bytes = val_str.as_bytes();
        card[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
        cards.push(card);

        // CONTINUE cards
        while !remaining.is_empty() {
            let is_last = remaining.len() <= 67;
            let chunk_len = if is_last { remaining.len() } else { 67 };
            let chunk = &remaining[..chunk_len];
            remaining = &remaining[chunk_len..];

            let mut card = [b' '; RECORD_SIZE];
            card[..8].copy_from_slice(b"CONTINUE");
            card[8] = b' ';
            card[9] = b' ';

            if is_last {
                let padded = pad_string_value(chunk);
                let val_str = format!("'{padded}'");
                let val_bytes = val_str.as_bytes();
                card[10..10 + val_bytes.len()].copy_from_slice(val_bytes);

                if let Some(cmt) = comment {
                    let cmt_start = 10 + val_bytes.len() + 1;
                    if cmt_start + 2 < RECORD_SIZE {
                        card[cmt_start] = b'/';
                        card[cmt_start + 1] = b' ';
                        let cmt_bytes = cmt.as_bytes();
                        let avail = RECORD_SIZE - cmt_start - 2;
                        let clen = cmt_bytes.len().min(avail);
                        card[cmt_start + 2..cmt_start + 2 + clen]
                            .copy_from_slice(&cmt_bytes[..clen]);
                    }
                }
            } else {
                let val_str = format!("'{chunk}&'");
                let val_bytes = val_str.as_bytes();
                card[10..10 + val_bytes.len()].copy_from_slice(val_bytes);
            }

            cards.push(card);
        }
    }

    cards
}

/// Parse the value+comment portion of a valued keyword (columns 11-80).
fn parse_value_comment(s: &str) -> Result<(HeaderValue, Option<String>)> {
    let trimmed = s.trim_start();

    if trimmed.is_empty() {
        return Ok((HeaderValue::Undefined, None));
    }

    // String value: starts with single quote
    if trimmed.starts_with('\'') {
        return parse_string_value(trimmed);
    }

    // Complex value: starts with (
    if trimmed.starts_with('(') {
        return parse_complex_value(trimmed);
    }

    // Logical: T or F at column 30 (index 20 from value start), but be lenient
    // Split on '/' to separate value from comment
    let (val_part, comment) = split_comment(trimmed);
    let val_trimmed = val_part.trim();

    if val_trimmed == "T" {
        return Ok((HeaderValue::Logical(true), comment));
    }
    if val_trimmed == "F" {
        return Ok((HeaderValue::Logical(false), comment));
    }

    // Try integer
    if let Ok(i) = val_trimmed.parse::<i64>() {
        return Ok((HeaderValue::Integer(i), comment));
    }

    // Try float (handle D exponent notation)
    let float_str = val_trimmed.replace('D', "E").replace('d', "e");
    if let Ok(f) = float_str.parse::<f64>() {
        return Ok((HeaderValue::Float(f), comment));
    }

    // Undefined
    Ok((HeaderValue::Undefined, comment))
}

fn parse_string_value(s: &str) -> Result<(HeaderValue, Option<String>)> {
    // Find the string content between quotes, handling doubled quotes
    let bytes = s.as_bytes();
    if bytes[0] != b'\'' {
        return Err(Error::InvalidKeyword("expected opening quote".into()));
    }

    let mut i = 1;
    let mut value = String::new();
    let mut closed = false;

    while i < bytes.len() {
        if bytes[i] == b'\'' {
            if i + 1 < bytes.len() && bytes[i + 1] == b'\'' {
                // Doubled quote -> literal quote
                value.push('\'');
                i += 2;
            } else {
                // End of string
                closed = true;
                i += 1;
                break;
            }
        } else {
            value.push(bytes[i] as char);
            i += 1;
        }
    }

    if !closed {
        return Err(Error::InvalidKeyword("unterminated string value".into()));
    }

    // Trim trailing spaces from string value
    let value = value.trim_end().to_string();

    // Check for continuation marker (trailing &)
    let continues = value.ends_with('&');
    let value = if continues {
        value[..value.len() - 1].to_string()
    } else {
        value
    };

    // Remaining after closing quote is potential comment
    let rest = &s[i..];
    let comment = parse_trailing_comment(rest);

    Ok((HeaderValue::String(value), comment))
}

fn parse_complex_value(s: &str) -> Result<(HeaderValue, Option<String>)> {
    let close = s
        .find(')')
        .ok_or_else(|| Error::InvalidKeyword("unterminated complex value".into()))?;
    let inner = &s[1..close];
    let parts: Vec<&str> = inner.split(',').collect();
    if parts.len() != 2 {
        return Err(Error::InvalidKeyword(
            "complex value must have two components".into(),
        ));
    }

    let rest = &s[close + 1..];
    let comment = parse_trailing_comment(rest);

    let a_str = parts[0].trim().replace('D', "E").replace('d', "e");
    let b_str = parts[1].trim().replace('D', "E").replace('d', "e");

    // Try integer complex first
    if let (Ok(a), Ok(b)) = (a_str.parse::<i64>(), b_str.parse::<i64>()) {
        return Ok((HeaderValue::ComplexInteger(a, b), comment));
    }

    let a: f64 = a_str
        .parse()
        .map_err(|_| Error::InvalidKeyword("invalid complex component".into()))?;
    let b: f64 = b_str
        .parse()
        .map_err(|_| Error::InvalidKeyword("invalid complex component".into()))?;

    Ok((HeaderValue::ComplexFloat(a, b), comment))
}

fn split_comment(s: &str) -> (&str, Option<String>) {
    if let Some(pos) = s.find('/') {
        let val = &s[..pos];
        let cmt = s[pos + 1..].trim();
        let comment = if cmt.is_empty() {
            None
        } else {
            Some(cmt.to_string())
        };
        (val, comment)
    } else {
        (s, None)
    }
}

fn parse_trailing_comment(s: &str) -> Option<String> {
    let trimmed = s.trim_start();
    if let Some(rest) = trimmed.strip_prefix('/') {
        let cmt = rest.trim();
        if cmt.is_empty() {
            None
        } else {
            Some(cmt.to_string())
        }
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_card(s: &str) -> [u8; RECORD_SIZE] {
        let mut card = [b' '; RECORD_SIZE];
        let bytes = s.as_bytes();
        let len = bytes.len().min(RECORD_SIZE);
        card[..len].copy_from_slice(&bytes[..len]);
        card
    }

    #[test]
    fn parse_integer() {
        let card = make_card("BITPIX  =                   16 / bits per data value");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.name, "BITPIX");
        assert_eq!(kw.value.unwrap().as_int(), Some(16));
        assert_eq!(kw.comment.unwrap(), "bits per data value");
    }

    #[test]
    fn parse_negative_integer() {
        let card = make_card("BITPIX  =                  -32 / IEEE float");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.value.unwrap().as_int(), Some(-32));
    }

    #[test]
    fn parse_float() {
        let card = make_card("BSCALE  =   1.000000000000E+00 / scale");
        let kw = Keyword::parse(&card).unwrap();
        let v = kw.value.unwrap().as_float().unwrap();
        assert!((v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn parse_logical() {
        let card = make_card("SIMPLE  =                    T / Standard FITS");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.value.unwrap().as_bool(), Some(true));
    }

    #[test]
    fn parse_string() {
        let card = make_card("TELESCOP= 'Hubble  '           / telescope name");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.value.unwrap().as_str(), Some("Hubble"));
    }

    #[test]
    fn parse_string_with_embedded_quote() {
        let card = make_card("COMMENT = 'It''s OK'           / note");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.value.unwrap().as_str(), Some("It's OK"));
    }

    #[test]
    fn parse_commentary() {
        let card = make_card("COMMENT This is a comment");
        let kw = Keyword::parse(&card).unwrap();
        assert_eq!(kw.name, "COMMENT");
        assert!(kw.value.is_none());
        assert!(kw.comment.unwrap().contains("This is a comment"));
    }

    #[test]
    fn round_trip_integer() {
        let kw = Keyword::with_value("NAXIS", HeaderValue::Integer(2), Some("number of axes"));
        let cards = kw.to_cards();
        assert_eq!(cards.len(), 1);
        let parsed = Keyword::parse(&cards[0]).unwrap();
        assert_eq!(parsed.name, "NAXIS");
        assert_eq!(parsed.value.unwrap().as_int(), Some(2));
    }

    #[test]
    fn round_trip_logical() {
        let kw = Keyword::with_value("SIMPLE", HeaderValue::Logical(true), None);
        let cards = kw.to_cards();
        let parsed = Keyword::parse(&cards[0]).unwrap();
        assert_eq!(parsed.value.unwrap().as_bool(), Some(true));
    }

    #[test]
    fn round_trip_string() {
        let kw = Keyword::with_value("OBJECT", HeaderValue::String("NGC 1234".into()), None);
        let cards = kw.to_cards();
        let parsed = Keyword::parse(&cards[0]).unwrap();
        assert_eq!(parsed.value.unwrap().as_str(), Some("NGC 1234"));
    }

    #[test]
    fn round_trip_float() {
        let kw = Keyword::with_value("CRVAL1", HeaderValue::Float(123.456), None);
        let cards = kw.to_cards();
        let parsed = Keyword::parse(&cards[0]).unwrap();
        let v = parsed.value.unwrap().as_float().unwrap();
        assert!((v - 123.456).abs() < 1e-10);
    }

    #[test]
    fn long_string_continue() {
        let long = "A".repeat(100);
        let kw = Keyword::with_value("LONGSTR", HeaderValue::String(long.clone()), None);
        let cards = kw.to_cards();
        assert!(cards.len() > 1);
        // Verify first card has CONTINUE indicator
        let first = std::str::from_utf8(&cards[0]).unwrap();
        assert!(first.contains('&'));
        // Verify second card starts with CONTINUE
        let second = std::str::from_utf8(&cards[1]).unwrap();
        assert!(second.starts_with("CONTINUE"));
    }

    #[test]
    fn end_card() {
        let kw = Keyword::new("END", None, None);
        let cards = kw.to_cards();
        assert_eq!(cards.len(), 1);
        let s = std::str::from_utf8(&cards[0]).unwrap();
        assert!(s.starts_with("END"));
        assert!(s[3..].chars().all(|c| c == ' '));
    }
}
