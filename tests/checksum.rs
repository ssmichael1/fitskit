use fitskit::ascii_table::{AsciiColumn, AsciiFormat};
use fitskit::checksum;
use fitskit::*;

#[test]
fn write_with_checksum_and_verify() {
    let pixels = PixelData::I16(vec![100, 200, 300, 400, 500, 600]);
    let img = ImageData::new(vec![3, 2], pixels);
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes_with_checksum().unwrap();

    // Read back and check DATASUM/CHECKSUM keywords are present
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();
    let hdr = &fits2.primary().header;

    let datasum_str = hdr.get_string("DATASUM").expect("DATASUM missing");
    assert!(!datasum_str.is_empty());
    let _datasum: u64 = datasum_str.parse().expect("DATASUM not a number");

    let checksum_str = hdr.get_string("CHECKSUM").expect("CHECKSUM missing");
    assert_eq!(checksum_str.len(), 16);
    assert!(checksum_str.bytes().all(|b| b.is_ascii_alphanumeric()));

    // Verify HDU checksum is valid (sums to all-ones)
    // We need the raw bytes for this
    assert!(checksum::verify_hdu(
        &bytes[..2880], // header is one block for this small HDU
        &bytes[2880..], // rest is data
    ));
}

#[test]
fn checksum_with_extensions() {
    use fitskit::bintable::{BinColumn, BinColumnType};

    let mut fits = FitsFile::with_empty_primary();

    let col = BinColumn {
        name: "X".into(),
        format: BinColumnType::D64(1),
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };
    let row_len = 8;
    let mut main_data = vec![0u8; 16];
    main_data[0..8].copy_from_slice(&3.125f64.to_be_bytes());
    main_data[8..16].copy_from_slice(&2.5f64.to_be_bytes());

    let table = fitskit::BinTable {
        columns: vec![col],
        nrows: 2,
        row_len,
        main_data,
        heap: Vec::new(),
    };
    fits.push_extension(Hdu::bintable_extension(table));

    let bytes = fits.to_bytes_with_checksum().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    // Both HDUs should have DATASUM
    for (i, hdu) in fits2.iter().enumerate() {
        assert!(
            hdu.header.get_string("DATASUM").is_some(),
            "HDU {i} missing DATASUM"
        );
        assert!(
            hdu.header.get_string("CHECKSUM").is_some(),
            "HDU {i} missing CHECKSUM"
        );
    }
}

#[test]
fn verify_datasum_on_read() {
    let pixels = PixelData::U8(vec![10, 20, 30, 40]);
    let img = ImageData::new(vec![4], pixels);
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes_with_checksum().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    // Verify DATASUM matches
    fits2.primary().verify_datasum().unwrap();
}

#[test]
fn write_without_checksum_has_no_keywords() {
    let fits = FitsFile::with_empty_primary();
    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert!(fits2.primary().header.get_string("DATASUM").is_none());
    assert!(fits2.primary().header.get_string("CHECKSUM").is_none());
}

#[test]
fn checksum_encode_decode_all_byte_values() {
    // Test that every possible byte value in a checksum can be encoded/decoded
    for byte_val in 0..=255u32 {
        let val = byte_val | (byte_val << 8) | (byte_val << 16) | (byte_val << 24);
        let encoded = checksum::encode_checksum(val, false);
        assert_eq!(encoded.len(), 16);
        assert!(
            encoded.bytes().all(|b| b.is_ascii_alphanumeric()),
            "non-alphanumeric for byte {byte_val:#x}: {encoded:?}"
        );
        let decoded = checksum::decode_checksum(&encoded, false);
        assert_eq!(decoded, val, "round-trip failed for {val:#010x}");
    }
}

/// Locate each HDU in `bytes` and check that the ones-complement sum of its
/// full header+data bytes is all-ones (the CHECKSUM invariant), and that the
/// stored DATASUM matches the padded data unit. This is the check `fitsverify`
/// performs, done on the raw bytes rather than through fitskit's own parser.
fn assert_raw_checksums_valid(bytes: &[u8]) -> usize {
    let mut cursor = std::io::Cursor::new(bytes);
    let mut n = 0;
    while (cursor.position() as usize) < bytes.len() {
        let start = cursor.position() as usize;
        let header = Header::read_from(&mut cursor).unwrap();
        let hdr_end = cursor.position() as usize;
        let data_len = header.data_byte_count().unwrap();
        let padded = data_len.div_ceil(2880) * 2880;
        let data = &bytes[hdr_end..hdr_end + padded];
        let stored: u32 = header.get_string("DATASUM").unwrap().parse().unwrap();
        assert_eq!(checksum::datasum(data), stored, "DATASUM of HDU {n}");
        assert!(
            checksum::verify_hdu(&bytes[start..hdr_end], data),
            "CHECKSUM of HDU {n} does not sum to all-ones"
        );
        cursor.set_position((hdr_end + padded) as u64);
        n += 1;
    }
    n
}

#[test]
fn raw_hdu_checksums_all_types() {
    // Image sizes chosen so data units need non-trivial padding and exceed the
    // ~256 KiB at which a 32-bit accumulator would overflow.
    let n = 700usize;
    let mut fits = FitsFile::with_primary_image(ImageData::new(
        vec![n, n],
        PixelData::F32((0..n * n).map(|i| -(i as f32) * 1e30).collect()),
    ));
    fits.push_extension(Hdu::image_extension(ImageData::new(
        vec![301, 5],
        PixelData::I16((0..1505).map(|i| (i as i16).wrapping_mul(-7)).collect()),
    )));
    let mut b = BinTableBuilder::new()
        .add_column("X", BinColumnType::D64(1))
        .add_column("S", BinColumnType::Char(3));
    for i in 0..1001 {
        b = b.push_row(|r| {
            r.write_f64(i as f64 * 1e300);
            r.write_string("abc", 3);
        });
    }
    fits.push_extension(Hdu::bintable_extension(b.build()));
    // ASCII table: 3 rows of 11 chars -> 33 bytes of data, 2847 bytes of blank fill
    let cols = vec![AsciiColumn {
        name: "V".into(),
        format: AsciiFormat::parse("F10.3").unwrap(),
        tbcol: 1,
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    }];
    let rows = b"   1.500   -2.250   3.125 ".to_vec();
    fits.push_extension(Hdu::ascii_table_extension(AsciiTable::build(cols, 3, rows)));

    let bytes = fits.to_bytes_with_checksum().unwrap();
    assert_eq!(assert_raw_checksums_valid(&bytes), 4);

    // ASCII table fill must be blanks, everything else zeros.
    let reread = FitsFile::from_bytes(&bytes).unwrap();
    let mut off = 0;
    for (i, hdu) in reread.iter().enumerate() {
        let mut hb = Vec::new();
        hdu.header.write_to(&mut hb).unwrap();
        off += hb.len();
        let len = hdu.data_byte_len();
        let padded = len.div_ceil(2880) * 2880;
        let fill = &bytes[off + len..off + padded];
        let expect = if i == 3 { b' ' } else { 0 };
        assert!(fill.iter().all(|&b| b == expect), "HDU {i} fill byte");
        off += padded;
        hdu.verify_datasum().unwrap();
    }
}
