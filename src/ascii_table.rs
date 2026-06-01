use crate::error::{Error, Result};
use crate::header::Header;
use crate::keyword::HeaderValue;

/// ASCII table column format types.
#[derive(Debug, Clone)]
pub enum AsciiFormat {
    /// Character string of width w
    Aw(usize),
    /// Integer of width w
    Iw(usize),
    /// Fixed-point float of width w with d decimal places
    Fwd(usize, usize),
    /// Exponential float of width w with d decimal places
    Ewd(usize, usize),
    /// Double exponential of width w with d decimal places
    Dwd(usize, usize),
}

impl AsciiFormat {
    pub fn parse(s: &str) -> Result<Self> {
        let s = s.trim();
        if s.is_empty() {
            return Err(Error::InvalidTableFormat("empty TFORM".into()));
        }

        let code = s.as_bytes()[0];
        let rest = &s[1..];

        match code {
            b'A' => {
                let w: usize = rest
                    .parse()
                    .map_err(|_| Error::InvalidTableFormat(format!("invalid Aw format: {s}")))?;
                Ok(AsciiFormat::Aw(w))
            }
            b'I' => {
                let w: usize = rest
                    .parse()
                    .map_err(|_| Error::InvalidTableFormat(format!("invalid Iw format: {s}")))?;
                Ok(AsciiFormat::Iw(w))
            }
            b'F' | b'E' | b'D' => {
                let (w, d) = parse_wd(rest, s)?;
                match code {
                    b'F' => Ok(AsciiFormat::Fwd(w, d)),
                    b'E' => Ok(AsciiFormat::Ewd(w, d)),
                    b'D' => Ok(AsciiFormat::Dwd(w, d)),
                    _ => unreachable!(),
                }
            }
            _ => Err(Error::InvalidTableFormat(format!(
                "unknown ASCII TFORM code: {s}"
            ))),
        }
    }

    pub fn width(&self) -> usize {
        match self {
            AsciiFormat::Aw(w) | AsciiFormat::Iw(w) => *w,
            AsciiFormat::Fwd(w, _) | AsciiFormat::Ewd(w, _) | AsciiFormat::Dwd(w, _) => *w,
        }
    }
}

fn parse_wd(rest: &str, original: &str) -> Result<(usize, usize)> {
    let parts: Vec<&str> = rest.split('.').collect();
    if parts.len() != 2 {
        return Err(Error::InvalidTableFormat(format!(
            "expected w.d format: {original}"
        )));
    }
    let w: usize = parts[0]
        .parse()
        .map_err(|_| Error::InvalidTableFormat(format!("invalid width: {original}")))?;
    let d: usize = parts[1]
        .parse()
        .map_err(|_| Error::InvalidTableFormat(format!("invalid decimals: {original}")))?;
    Ok((w, d))
}

/// Column descriptor for an ASCII table.
#[derive(Debug, Clone)]
pub struct AsciiColumn {
    pub name: String,
    pub format: AsciiFormat,
    pub tbcol: usize, // 1-based start column
    pub tscal: f64,
    pub tzero: f64,
    pub tunit: Option<String>,
}

/// An ASCII table extension.
#[derive(Debug, Clone)]
pub struct AsciiTable {
    pub columns: Vec<AsciiColumn>,
    pub nrows: usize,
    pub row_len: usize, // NAXIS1
    pub raw_data: Vec<u8>,
}

impl AsciiTable {
    /// Read column definitions from header and raw data.
    pub fn from_header_and_data(header: &Header, data: &[u8]) -> Result<Self> {
        let nrows = header.require_int("NAXIS2")? as usize;
        let row_len = header.require_int("NAXIS1")? as usize;
        let tfields = header.require_int("TFIELDS")? as usize;

        let mut columns = Vec::with_capacity(tfields);
        for i in 1..=tfields {
            let name = header
                .get_string(&format!("TTYPE{i}"))
                .unwrap_or("")
                .to_string();

            let fmt_str = header
                .get_string(&format!("TFORM{i}"))
                .ok_or_else(|| Error::MissingKeyword(format!("TFORM{i}")))?;
            let format = AsciiFormat::parse(fmt_str)?;

            let tbcol = header.require_int(&format!("TBCOL{i}"))? as usize;

            let tscal = header.get_float(&format!("TSCAL{i}")).unwrap_or(1.0);
            let tzero = header.get_float(&format!("TZERO{i}")).unwrap_or(0.0);
            let tunit = header
                .get_string(&format!("TUNIT{i}"))
                .map(|s| s.to_string());

            columns.push(AsciiColumn {
                name,
                format,
                tbcol,
                tscal,
                tzero,
                tunit,
            });
        }

        Ok(AsciiTable {
            columns,
            nrows,
            row_len,
            raw_data: data.to_vec(),
        })
    }

    /// Get a string cell value (raw, before scaling).
    pub fn get_cell_raw(&self, row: usize, col: usize) -> Result<&str> {
        if row >= self.nrows || col >= self.columns.len() {
            return Err(Error::InvalidFormat("cell index out of bounds".into()));
        }
        let column = &self.columns[col];
        let start = row * self.row_len + (column.tbcol - 1);
        let end = start + column.format.width();
        if end > self.raw_data.len() {
            return Err(Error::DataSizeMismatch {
                expected: end,
                actual: self.raw_data.len(),
            });
        }
        std::str::from_utf8(&self.raw_data[start..end])
            .map_err(|_| Error::InvalidFormat("non-UTF8 table data".into()))
    }

    /// Get a float value with TSCAL/TZERO applied.
    pub fn get_float(&self, row: usize, col: usize) -> Result<f64> {
        let raw = self.get_cell_raw(row, col)?.trim();
        let column = &self.columns[col];

        // Handle Fortran D exponent
        let s = crate::keyword::fortran_exp(raw);
        let val: f64 = s
            .parse()
            .map_err(|_| Error::InvalidTableFormat(format!("cannot parse float: '{raw}'")))?;

        Ok(column.tzero + column.tscal * val)
    }

    /// Get an integer value with TSCAL/TZERO applied.
    pub fn get_int(&self, row: usize, col: usize) -> Result<i64> {
        let raw = self.get_cell_raw(row, col)?.trim();
        let column = &self.columns[col];

        let val: i64 = raw
            .parse()
            .map_err(|_| Error::InvalidTableFormat(format!("cannot parse integer: '{raw}'")))?;

        Ok((column.tzero + column.tscal * val as f64) as i64)
    }

    /// Get a string value (no scaling).
    pub fn get_string(&self, row: usize, col: usize) -> Result<String> {
        Ok(self.get_cell_raw(row, col)?.trim_end().to_string())
    }

    /// Populate header keywords for this table.
    pub fn fill_header(&self, header: &mut Header) {
        header.set(
            "XTENSION",
            HeaderValue::String("TABLE".into()),
            Some("ASCII table extension"),
        );
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(2), None);
        header.set("NAXIS1", HeaderValue::Integer(self.row_len as i64), None);
        header.set("NAXIS2", HeaderValue::Integer(self.nrows as i64), None);
        header.set("PCOUNT", HeaderValue::Integer(0), None);
        header.set("GCOUNT", HeaderValue::Integer(1), None);
        header.set(
            "TFIELDS",
            HeaderValue::Integer(self.columns.len() as i64),
            None,
        );

        for (i, col) in self.columns.iter().enumerate() {
            let idx = i + 1;
            header.set(
                &format!("TTYPE{idx}"),
                HeaderValue::String(col.name.clone()),
                None,
            );
            header.set(
                &format!("TBCOL{idx}"),
                HeaderValue::Integer(col.tbcol as i64),
                None,
            );

            let fmt_str = match &col.format {
                AsciiFormat::Aw(w) => format!("A{w}"),
                AsciiFormat::Iw(w) => format!("I{w}"),
                AsciiFormat::Fwd(w, d) => format!("F{w}.{d}"),
                AsciiFormat::Ewd(w, d) => format!("E{w}.{d}"),
                AsciiFormat::Dwd(w, d) => format!("D{w}.{d}"),
            };
            header.set(&format!("TFORM{idx}"), HeaderValue::String(fmt_str), None);

            if col.tscal != 1.0 {
                header.set(&format!("TSCAL{idx}"), HeaderValue::Float(col.tscal), None);
            }
            if col.tzero != 0.0 {
                header.set(&format!("TZERO{idx}"), HeaderValue::Float(col.tzero), None);
            }
            if let Some(ref unit) = col.tunit {
                header.set(
                    &format!("TUNIT{idx}"),
                    HeaderValue::String(unit.clone()),
                    None,
                );
            }
        }
    }

    /// Build an AsciiTable from column definitions and row data.
    pub fn build(columns: Vec<AsciiColumn>, nrows: usize, row_data: Vec<u8>) -> Self {
        let row_len = if let Some(last) = columns.last() {
            last.tbcol - 1 + last.format.width()
        } else {
            0
        };
        AsciiTable {
            columns,
            nrows,
            row_len,
            raw_data: row_data,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_ascii_formats() {
        let f = AsciiFormat::parse("A10").unwrap();
        assert!(matches!(f, AsciiFormat::Aw(10)));

        let f = AsciiFormat::parse("I5").unwrap();
        assert!(matches!(f, AsciiFormat::Iw(5)));

        let f = AsciiFormat::parse("F10.3").unwrap();
        assert!(matches!(f, AsciiFormat::Fwd(10, 3)));

        let f = AsciiFormat::parse("E15.7").unwrap();
        assert!(matches!(f, AsciiFormat::Ewd(15, 7)));

        let f = AsciiFormat::parse("D25.17").unwrap();
        assert!(matches!(f, AsciiFormat::Dwd(25, 17)));
    }

    #[test]
    fn invalid_format() {
        assert!(AsciiFormat::parse("Z10").is_err());
        assert!(AsciiFormat::parse("").is_err());
    }
}
