use crate::error::{Error, Result};
use crate::header::Header;
use crate::keyword::HeaderValue;

/// Binary table column type codes per FITS standard.
#[derive(Debug, Clone)]
pub enum BinColumnType {
    /// L: Logical (1 byte)
    Logical(usize),
    /// X: Bit array (ceil(n/8) bytes)
    Bit(usize),
    /// B: Unsigned byte
    Byte(usize),
    /// I: 16-bit integer
    I16(usize),
    /// J: 32-bit integer
    J32(usize),
    /// K: 64-bit integer
    K64(usize),
    /// A: Character
    Char(usize),
    /// E: 32-bit float
    E32(usize),
    /// D: 64-bit float
    D64(usize),
    /// C: 32-bit complex float (pair of E)
    C64(usize),
    /// M: 64-bit complex float (pair of D)
    M128(usize),
    /// P: 32-bit variable-length array descriptor (2×i32)
    VarP(char),
    /// Q: 64-bit variable-length array descriptor (2×i64)
    VarQ(char),
}

impl BinColumnType {
    /// Parse a TFORMn value like "1J", "10E", "8A", "1PJ(1000)", "1QD(2000)" etc.
    pub fn parse(s: &str) -> Result<Self> {
        let s = s.trim();
        if s.is_empty() {
            return Err(Error::InvalidTableFormat("empty TFORM".into()));
        }

        // Check for variable-length array: rPt or rQt (with optional max in parens)
        // Pattern: optional_repeat + P/Q + type_code + optional(maxlen)
        let bytes = s.as_bytes();

        // Find the type code letter position
        let mut code_pos = 0;
        while code_pos < bytes.len() && bytes[code_pos].is_ascii_digit() {
            code_pos += 1;
        }

        if code_pos >= bytes.len() {
            return Err(Error::InvalidTableFormat(format!(
                "no type code in TFORM: {s}"
            )));
        }

        let repeat: usize = if code_pos == 0 {
            1
        } else {
            s[..code_pos]
                .parse()
                .map_err(|_| Error::InvalidTableFormat(format!("invalid repeat count: {s}")))?
        };

        let code = bytes[code_pos];

        // Variable-length array
        if code == b'P' || code == b'Q' {
            if code_pos + 1 >= bytes.len() {
                return Err(Error::InvalidTableFormat(format!(
                    "P/Q format needs element type: {s}"
                )));
            }
            let elem_type = bytes[code_pos + 1] as char;
            if code == b'P' {
                return Ok(BinColumnType::VarP(elem_type));
            } else {
                return Ok(BinColumnType::VarQ(elem_type));
            }
        }

        match code {
            b'L' => Ok(BinColumnType::Logical(repeat)),
            b'X' => Ok(BinColumnType::Bit(repeat)),
            b'B' => Ok(BinColumnType::Byte(repeat)),
            b'I' => Ok(BinColumnType::I16(repeat)),
            b'J' => Ok(BinColumnType::J32(repeat)),
            b'K' => Ok(BinColumnType::K64(repeat)),
            b'A' => Ok(BinColumnType::Char(repeat)),
            b'E' => Ok(BinColumnType::E32(repeat)),
            b'D' => Ok(BinColumnType::D64(repeat)),
            b'C' => Ok(BinColumnType::C64(repeat)),
            b'M' => Ok(BinColumnType::M128(repeat)),
            _ => Err(Error::InvalidTableFormat(format!(
                "unknown BINTABLE type code: {}",
                code as char
            ))),
        }
    }

    /// Number of bytes this column occupies per row in the main table.
    pub fn byte_size(&self) -> usize {
        match self {
            BinColumnType::Logical(n) => *n,
            BinColumnType::Bit(n) => (*n).div_ceil(8),
            BinColumnType::Byte(n) => *n,
            BinColumnType::I16(n) => n * 2,
            BinColumnType::J32(n) => n * 4,
            BinColumnType::K64(n) => n * 8,
            BinColumnType::Char(n) => *n,
            BinColumnType::E32(n) => n * 4,
            BinColumnType::D64(n) => n * 8,
            BinColumnType::C64(n) => n * 8,
            BinColumnType::M128(n) => n * 16,
            BinColumnType::VarP(_) => 8,  // 2 × i32
            BinColumnType::VarQ(_) => 16, // 2 × i64
        }
    }

    /// Format string for writing TFORM.
    pub fn to_tform_string(&self) -> String {
        match self {
            BinColumnType::Logical(n) => format!("{n}L"),
            BinColumnType::Bit(n) => format!("{n}X"),
            BinColumnType::Byte(n) => format!("{n}B"),
            BinColumnType::I16(n) => format!("{n}I"),
            BinColumnType::J32(n) => format!("{n}J"),
            BinColumnType::K64(n) => format!("{n}K"),
            BinColumnType::Char(n) => format!("{n}A"),
            BinColumnType::E32(n) => format!("{n}E"),
            BinColumnType::D64(n) => format!("{n}D"),
            BinColumnType::C64(n) => format!("{n}C"),
            BinColumnType::M128(n) => format!("{n}M"),
            BinColumnType::VarP(t) => format!("1P{t}"),
            BinColumnType::VarQ(t) => format!("1Q{t}"),
        }
    }
}

/// A column in a binary table.
#[derive(Debug, Clone)]
pub struct BinColumn {
    pub name: String,
    pub format: BinColumnType,
    pub tscal: f64,
    pub tzero: f64,
    pub tunit: Option<String>,
}

/// Cell value from a binary table.
#[derive(Debug, Clone)]
pub enum BinCellValue {
    Logical(Vec<bool>),
    Bytes(Vec<u8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
    F32(Vec<f32>),
    F64(Vec<f64>),
    String(String),
    ComplexF32(Vec<(f32, f32)>),
    ComplexF64(Vec<(f64, f64)>),
    Bits(Vec<u8>, usize), // bytes + bit count
}

/// A binary table extension.
#[derive(Debug, Clone)]
pub struct BinTable {
    pub columns: Vec<BinColumn>,
    pub nrows: usize,
    pub row_len: usize, // NAXIS1
    /// Main data area (nrows × row_len bytes).
    pub main_data: Vec<u8>,
    /// Heap area (variable-length array storage).
    pub heap: Vec<u8>,
}

impl BinTable {
    /// Read from header and raw data (including heap).
    pub fn from_header_and_data(header: &Header, data: &[u8]) -> Result<Self> {
        let nrows = header.require_int("NAXIS2")? as usize;
        let row_len = header.require_int("NAXIS1")? as usize;
        let tfields = header.require_int("TFIELDS")? as usize;
        let pcount = header.get_int("PCOUNT").unwrap_or(0) as usize;

        let main_size = nrows * row_len;
        if data.len() < main_size {
            return Err(Error::DataSizeMismatch {
                expected: main_size,
                actual: data.len(),
            });
        }

        let main_data = data[..main_size].to_vec();
        let heap = if pcount > 0 && data.len() >= main_size + pcount {
            data[main_size..main_size + pcount].to_vec()
        } else {
            Vec::new()
        };

        let mut columns = Vec::with_capacity(tfields);
        for i in 1..=tfields {
            let name = header
                .get_string(&format!("TTYPE{i}"))
                .unwrap_or("")
                .to_string();

            let fmt_str = header
                .get_string(&format!("TFORM{i}"))
                .ok_or_else(|| Error::MissingKeyword(format!("TFORM{i}")))?;
            let format = BinColumnType::parse(fmt_str)?;

            let tscal = header.get_float(&format!("TSCAL{i}")).unwrap_or(1.0);
            let tzero = header.get_float(&format!("TZERO{i}")).unwrap_or(0.0);
            let tunit = header
                .get_string(&format!("TUNIT{i}"))
                .map(|s| s.to_string());

            columns.push(BinColumn {
                name,
                format,
                tscal,
                tzero,
                tunit,
            });
        }

        Ok(BinTable {
            columns,
            nrows,
            row_len,
            main_data,
            heap,
        })
    }

    /// Get the byte offset of a column within a row.
    fn column_offset(&self, col: usize) -> usize {
        self.columns[..col]
            .iter()
            .map(|c| c.format.byte_size())
            .sum()
    }

    /// Get the raw bytes for a cell.
    pub fn cell_bytes(&self, row: usize, col: usize) -> Result<&[u8]> {
        if row >= self.nrows || col >= self.columns.len() {
            return Err(Error::InvalidFormat("cell index out of bounds".into()));
        }
        let col_offset = self.column_offset(col);
        let col_size = self.columns[col].format.byte_size();
        let start = row * self.row_len + col_offset;
        let end = start + col_size;
        if end > self.main_data.len() {
            return Err(Error::DataSizeMismatch {
                expected: end,
                actual: self.main_data.len(),
            });
        }
        Ok(&self.main_data[start..end])
    }

    /// Read a cell and decode it into a typed value.
    pub fn get_cell(&self, row: usize, col: usize) -> Result<BinCellValue> {
        let bytes = self.cell_bytes(row, col)?;
        let fmt = &self.columns[col].format;

        match fmt {
            BinColumnType::Logical(n) => {
                let vals: Vec<bool> = bytes[..*n].iter().map(|&b| b == b'T').collect();
                Ok(BinCellValue::Logical(vals))
            }
            BinColumnType::Bit(n) => Ok(BinCellValue::Bits(bytes.to_vec(), *n)),
            BinColumnType::Byte(n) => Ok(BinCellValue::Bytes(bytes[..*n].to_vec())),
            BinColumnType::I16(n) => {
                let vals: Vec<i16> = bytes
                    .chunks_exact(2)
                    .take(*n)
                    .map(|c| i16::from_be_bytes([c[0], c[1]]))
                    .collect();
                Ok(BinCellValue::I16(vals))
            }
            BinColumnType::J32(n) => {
                let vals: Vec<i32> = bytes
                    .chunks_exact(4)
                    .take(*n)
                    .map(|c| i32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect();
                Ok(BinCellValue::I32(vals))
            }
            BinColumnType::K64(n) => {
                let vals: Vec<i64> = bytes
                    .chunks_exact(8)
                    .take(*n)
                    .map(|c| i64::from_be_bytes(c.try_into().unwrap()))
                    .collect();
                Ok(BinCellValue::I64(vals))
            }
            BinColumnType::Char(n) => {
                let s = std::str::from_utf8(&bytes[..*n])
                    .map_err(|_| Error::InvalidFormat("non-UTF8 in char column".into()))?
                    .trim_end()
                    .to_string();
                Ok(BinCellValue::String(s))
            }
            BinColumnType::E32(n) => {
                let vals: Vec<f32> = bytes
                    .chunks_exact(4)
                    .take(*n)
                    .map(|c| f32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect();
                Ok(BinCellValue::F32(vals))
            }
            BinColumnType::D64(n) => {
                let vals: Vec<f64> = bytes
                    .chunks_exact(8)
                    .take(*n)
                    .map(|c| f64::from_be_bytes(c.try_into().unwrap()))
                    .collect();
                Ok(BinCellValue::F64(vals))
            }
            BinColumnType::C64(n) => {
                let vals: Vec<(f32, f32)> = bytes
                    .chunks_exact(8)
                    .take(*n)
                    .map(|c| {
                        let r = f32::from_be_bytes([c[0], c[1], c[2], c[3]]);
                        let i = f32::from_be_bytes([c[4], c[5], c[6], c[7]]);
                        (r, i)
                    })
                    .collect();
                Ok(BinCellValue::ComplexF32(vals))
            }
            BinColumnType::M128(n) => {
                let vals: Vec<(f64, f64)> = bytes
                    .chunks_exact(16)
                    .take(*n)
                    .map(|c| {
                        let r = f64::from_be_bytes(c[..8].try_into().unwrap());
                        let i = f64::from_be_bytes(c[8..16].try_into().unwrap());
                        (r, i)
                    })
                    .collect();
                Ok(BinCellValue::ComplexF64(vals))
            }
            BinColumnType::VarP(elem_type) => {
                // P descriptor: 2 × i32 (count, offset)
                let count = i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
                let offset = i32::from_be_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize;
                self.read_heap_array(*elem_type, count, offset)
            }
            BinColumnType::VarQ(elem_type) => {
                // Q descriptor: 2 × i64 (count, offset)
                let count = i64::from_be_bytes(bytes[..8].try_into().unwrap()) as usize;
                let offset = i64::from_be_bytes(bytes[8..16].try_into().unwrap()) as usize;
                self.read_heap_array(*elem_type, count, offset)
            }
        }
    }

    /// Read variable-length array data from the heap.
    fn read_heap_array(
        &self,
        elem_type: char,
        count: usize,
        offset: usize,
    ) -> Result<BinCellValue> {
        if count == 0 {
            return match elem_type {
                'B' => Ok(BinCellValue::Bytes(Vec::new())),
                'I' => Ok(BinCellValue::I16(Vec::new())),
                'J' => Ok(BinCellValue::I32(Vec::new())),
                'K' => Ok(BinCellValue::I64(Vec::new())),
                'E' => Ok(BinCellValue::F32(Vec::new())),
                'D' => Ok(BinCellValue::F64(Vec::new())),
                'A' => Ok(BinCellValue::String(String::new())),
                'L' => Ok(BinCellValue::Logical(Vec::new())),
                _ => Err(Error::InvalidTableFormat(format!(
                    "unsupported VLA element type: {elem_type}"
                ))),
            };
        }

        let elem_size = match elem_type {
            'L' | 'B' | 'A' => 1,
            'I' => 2,
            'J' | 'E' => 4,
            'K' | 'D' => 8,
            'C' => 8,
            'M' => 16,
            _ => {
                return Err(Error::InvalidTableFormat(format!(
                    "unsupported VLA element type: {elem_type}"
                )))
            }
        };

        let end = offset + count * elem_size;
        if end > self.heap.len() {
            return Err(Error::DataSizeMismatch {
                expected: end,
                actual: self.heap.len(),
            });
        }

        let data = &self.heap[offset..end];

        match elem_type {
            'L' => Ok(BinCellValue::Logical(
                data.iter().map(|&b| b == b'T').collect(),
            )),
            'B' => Ok(BinCellValue::Bytes(data.to_vec())),
            'A' => {
                let s = std::str::from_utf8(data)
                    .map_err(|_| Error::InvalidFormat("non-UTF8 in VLA char".into()))?
                    .trim_end()
                    .to_string();
                Ok(BinCellValue::String(s))
            }
            'I' => Ok(BinCellValue::I16(
                data.chunks_exact(2)
                    .map(|c| i16::from_be_bytes([c[0], c[1]]))
                    .collect(),
            )),
            'J' => Ok(BinCellValue::I32(
                data.chunks_exact(4)
                    .map(|c| i32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect(),
            )),
            'K' => Ok(BinCellValue::I64(
                data.chunks_exact(8)
                    .map(|c| i64::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            )),
            'E' => Ok(BinCellValue::F32(
                data.chunks_exact(4)
                    .map(|c| f32::from_be_bytes([c[0], c[1], c[2], c[3]]))
                    .collect(),
            )),
            'D' => Ok(BinCellValue::F64(
                data.chunks_exact(8)
                    .map(|c| f64::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            )),
            _ => Err(Error::InvalidTableFormat(format!(
                "unsupported VLA type: {elem_type}"
            ))),
        }
    }

    /// Populate header keywords for this binary table.
    pub fn fill_header(&self, header: &mut Header) {
        header.set(
            "XTENSION",
            HeaderValue::String("BINTABLE".into()),
            Some("binary table extension"),
        );
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(2), None);
        header.set("NAXIS1", HeaderValue::Integer(self.row_len as i64), None);
        header.set("NAXIS2", HeaderValue::Integer(self.nrows as i64), None);
        header.set("PCOUNT", HeaderValue::Integer(self.heap.len() as i64), None);
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
                &format!("TFORM{idx}"),
                HeaderValue::String(col.format.to_tform_string()),
                None,
            );

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

    /// Get the serialized data bytes (main table + heap).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = self.main_data.clone();
        out.extend_from_slice(&self.heap);
        out
    }
}

/// Builder for constructing a `BinTable` row by row.
///
/// ```
/// use fitskit::{BinTableBuilder, BinColumnType};
///
/// let table = BinTableBuilder::new()
///     .add_column("RA", BinColumnType::D64(1))
///     .add_column("DEC", BinColumnType::D64(1))
///     .add_column("MAG", BinColumnType::E32(1))
///     .push_row(|row| {
///         row.write_f64(180.0);
///         row.write_f64(45.0);
///         row.write_f32(12.5);
///     })
///     .push_row(|row| {
///         row.write_f64(90.0);
///         row.write_f64(-30.0);
///         row.write_f32(8.2);
///     })
///     .build();
///
/// assert_eq!(table.nrows, 2);
/// assert_eq!(table.columns.len(), 3);
/// ```
pub struct BinTableBuilder {
    columns: Vec<BinColumn>,
    row_len: usize,
    main_data: Vec<u8>,
    heap: Vec<u8>,
    nrows: usize,
}

impl BinTableBuilder {
    pub fn new() -> Self {
        BinTableBuilder {
            columns: Vec::new(),
            row_len: 0,
            main_data: Vec::new(),
            heap: Vec::new(),
            nrows: 0,
        }
    }

    /// Add a column with default TSCAL=1, TZERO=0, no unit.
    pub fn add_column(mut self, name: &str, format: BinColumnType) -> Self {
        self.row_len += format.byte_size();
        self.columns.push(BinColumn {
            name: name.to_string(),
            format,
            tscal: 1.0,
            tzero: 0.0,
            tunit: None,
        });
        self
    }

    /// Add a column with full metadata.
    pub fn add_column_full(
        mut self,
        name: &str,
        format: BinColumnType,
        tscal: f64,
        tzero: f64,
        tunit: Option<&str>,
    ) -> Self {
        self.row_len += format.byte_size();
        self.columns.push(BinColumn {
            name: name.to_string(),
            format,
            tscal,
            tzero,
            tunit: tunit.map(|s| s.to_string()),
        });
        self
    }

    /// Append a row using a closure that writes cell values via a [`RowWriter`].
    pub fn push_row<F: FnOnce(&mut RowWriter)>(mut self, f: F) -> Self {
        let mut writer = RowWriter {
            data: Vec::with_capacity(self.row_len),
            heap: &mut self.heap,
        };
        f(&mut writer);
        // Pad row to row_len if the closure wrote fewer bytes
        writer.data.resize(self.row_len, 0);
        self.main_data.extend_from_slice(&writer.data);
        self.nrows += 1;
        self
    }

    /// Consume the builder and produce a `BinTable`.
    pub fn build(self) -> BinTable {
        BinTable {
            columns: self.columns,
            nrows: self.nrows,
            row_len: self.row_len,
            main_data: self.main_data,
            heap: self.heap,
        }
    }
}

impl Default for BinTableBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Writes cell values into a single binary table row.
pub struct RowWriter<'a> {
    data: Vec<u8>,
    heap: &'a mut Vec<u8>,
}

impl<'a> RowWriter<'a> {
    pub fn write_bool(&mut self, val: bool) {
        self.data.push(if val { b'T' } else { b'F' });
    }

    pub fn write_u8(&mut self, val: u8) {
        self.data.push(val);
    }

    pub fn write_i16(&mut self, val: i16) {
        self.data.extend_from_slice(&val.to_be_bytes());
    }

    pub fn write_i32(&mut self, val: i32) {
        self.data.extend_from_slice(&val.to_be_bytes());
    }

    pub fn write_i64(&mut self, val: i64) {
        self.data.extend_from_slice(&val.to_be_bytes());
    }

    pub fn write_f32(&mut self, val: f32) {
        self.data.extend_from_slice(&val.to_be_bytes());
    }

    pub fn write_f64(&mut self, val: f64) {
        self.data.extend_from_slice(&val.to_be_bytes());
    }

    /// Write a fixed-length ASCII string, padded with spaces to `len` bytes.
    pub fn write_string(&mut self, val: &str, len: usize) {
        let bytes = val.as_bytes();
        let n = bytes.len().min(len);
        self.data.extend_from_slice(&bytes[..n]);
        for _ in n..len {
            self.data.push(b' ');
        }
    }

    /// Write a variable-length array using a P descriptor (2×i32: count, heap offset).
    /// The data is appended to the heap. `data_fn` writes the element bytes.
    pub fn write_var_p<F: FnOnce(&mut Vec<u8>)>(&mut self, count: i32, data_fn: F) {
        let offset = self.heap.len() as i32;
        self.data.extend_from_slice(&count.to_be_bytes());
        self.data.extend_from_slice(&offset.to_be_bytes());
        data_fn(self.heap);
    }

    /// Write raw bytes directly.
    pub fn write_bytes(&mut self, bytes: &[u8]) {
        self.data.extend_from_slice(bytes);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tform_codes() {
        assert!(matches!(
            BinColumnType::parse("1J").unwrap(),
            BinColumnType::J32(1)
        ));
        assert!(matches!(
            BinColumnType::parse("10E").unwrap(),
            BinColumnType::E32(10)
        ));
        assert!(matches!(
            BinColumnType::parse("8A").unwrap(),
            BinColumnType::Char(8)
        ));
        assert!(matches!(
            BinColumnType::parse("D").unwrap(),
            BinColumnType::D64(1)
        ));
        assert!(matches!(
            BinColumnType::parse("1PJ").unwrap(),
            BinColumnType::VarP('J')
        ));
        assert!(matches!(
            BinColumnType::parse("1QD").unwrap(),
            BinColumnType::VarQ('D')
        ));
    }

    #[test]
    fn byte_sizes() {
        assert_eq!(BinColumnType::Logical(1).byte_size(), 1);
        assert_eq!(BinColumnType::I16(3).byte_size(), 6);
        assert_eq!(BinColumnType::J32(1).byte_size(), 4);
        assert_eq!(BinColumnType::D64(2).byte_size(), 16);
        assert_eq!(BinColumnType::VarP('J').byte_size(), 8);
        assert_eq!(BinColumnType::VarQ('D').byte_size(), 16);
        assert_eq!(BinColumnType::Bit(10).byte_size(), 2);
    }

    #[test]
    fn tform_round_trip() {
        for s in &[
            "1J", "10E", "8A", "1D", "3I", "1L", "20X", "1C", "1M", "1PJ", "1QD",
        ] {
            let ct = BinColumnType::parse(s).unwrap();
            let out = ct.to_tform_string();
            let ct2 = BinColumnType::parse(&out).unwrap();
            assert_eq!(ct.byte_size(), ct2.byte_size());
        }
    }
}
