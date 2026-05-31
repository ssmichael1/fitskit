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
