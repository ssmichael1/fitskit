use fitskit::*;

#[test]
fn round_trip_empty_primary() {
    let fits = FitsFile::with_empty_primary();
    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.len(), 1);
    assert!(matches!(fits2.primary().data, HduData::Empty));
    assert_eq!(fits2.primary().header.get_bool("SIMPLE"), Some(true));
}

#[test]
fn round_trip_u8_image() {
    let pixels = PixelData::U8(vec![0, 128, 255, 1, 2, 3]);
    let img = ImageData::new(vec![3, 2], pixels);
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    match &fits2.primary().data {
        HduData::Image(img) => {
            assert_eq!(img.axes, vec![3, 2]);
            assert_eq!(img.bitpix(), Bitpix::U8);
            if let PixelData::U8(data) = &img.pixels {
                assert_eq!(data, &[0, 128, 255, 1, 2, 3]);
            } else {
                panic!("wrong pixel type");
            }
        }
        _ => panic!("expected image data"),
    }
}

#[test]
fn round_trip_i16_image() {
    let pixels: Vec<i16> = vec![-32768, -1, 0, 1, 32767];
    let img = ImageData::new(vec![5], PixelData::I16(pixels.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        if let PixelData::I16(data) = &img.pixels {
            assert_eq!(data, &pixels);
        } else {
            panic!("wrong pixel type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_i32_image() {
    let pixels: Vec<i32> = vec![i32::MIN, -1, 0, 1, i32::MAX];
    let img = ImageData::new(vec![5], PixelData::I32(pixels.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        if let PixelData::I32(data) = &img.pixels {
            assert_eq!(data, &pixels);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_i64_image() {
    let pixels: Vec<i64> = vec![i64::MIN, 0, i64::MAX];
    let img = ImageData::new(vec![3], PixelData::I64(pixels.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        if let PixelData::I64(data) = &img.pixels {
            assert_eq!(data, &pixels);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_f32_image() {
    let pixels: Vec<f32> = vec![0.0, 1.0, -1.0, 3.125, f32::MAX];
    let img = ImageData::new(vec![5], PixelData::F32(pixels.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        if let PixelData::F32(data) = &img.pixels {
            assert_eq!(data, &pixels);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_f64_image() {
    let pixels: Vec<f64> = vec![0.0, std::f64::consts::PI, -1e300, 1e300];
    let img = ImageData::new(vec![4], PixelData::F64(pixels.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        if let PixelData::F64(data) = &img.pixels {
            assert_eq!(data, &pixels);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_2d_image() {
    // 4x3 image
    let mut px = Vec::new();
    for y in 0..3u16 {
        for x in 0..4u16 {
            px.push((y * 4 + x) as i16);
        }
    }
    let img = ImageData::new(vec![4, 3], PixelData::I16(px.clone()));
    let fits = FitsFile::with_primary_image(img);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::Image(img) = &fits2.primary().data {
        assert_eq!(img.axes, vec![4, 3]);
        assert_eq!(img.width(), Some(4));
        assert_eq!(img.height(), Some(3));
        if let PixelData::I16(data) = &img.pixels {
            assert_eq!(data, &px);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_image_extension() {
    let primary_img = ImageData::new(vec![2, 2], PixelData::U8(vec![1, 2, 3, 4]));
    let mut fits = FitsFile::with_primary_image(primary_img);

    let ext_img = ImageData::new(vec![3], PixelData::F32(vec![1.0, 2.0, 3.0]));
    fits.push_extension(Hdu::image_extension(ext_img));

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.len(), 2);

    // Check primary
    if let HduData::Image(img) = &fits2.primary().data {
        assert_eq!(img.bitpix(), Bitpix::U8);
    } else {
        panic!("expected primary image");
    }

    // Check extension
    if let HduData::Image(img) = &fits2.extensions()[0].data {
        assert_eq!(img.bitpix(), Bitpix::F32);
        if let PixelData::F32(data) = &img.pixels {
            assert_eq!(data, &[1.0, 2.0, 3.0]);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected image extension");
    }
}

#[test]
fn round_trip_bscale_bzero() {
    // Unsigned u16 via BZERO=32768
    let raw_pixels: Vec<i16> = vec![-32768, 0, 32767]; // physical: 0, 32768, 65535
    let img = ImageData::new(vec![3], PixelData::I16(raw_pixels));

    let mut fits = FitsFile::with_primary_image(img);
    fits.primary_mut()
        .header
        .set("BSCALE", HeaderValue::Float(1.0), None);
    fits.primary_mut()
        .header
        .set("BZERO", HeaderValue::Float(32768.0), None);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    let bscale = fits2.primary().header.get_float("BSCALE").unwrap();
    let bzero = fits2.primary().header.get_float("BZERO").unwrap();
    assert!((bscale - 1.0).abs() < 1e-10);
    assert!((bzero - 32768.0).abs() < 1e-10);

    if let HduData::Image(img) = &fits2.primary().data {
        let scaled = img.scaled_values(bscale, bzero);
        assert_eq!(scaled, vec![0.0, 32768.0, 65535.0]);
    } else {
        panic!("expected image");
    }
}

#[test]
fn round_trip_bintable() {
    use fitskit::bintable::{BinCellValue, BinColumn, BinColumnType};

    // Build a simple bintable with 2 rows, 3 columns: i32, f64, 8-char string
    let col_j = BinColumn {
        name: "COUNT".into(),
        format: BinColumnType::J32(1),
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };
    let col_d = BinColumn {
        name: "VALUE".into(),
        format: BinColumnType::D64(1),
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };
    let col_a = BinColumn {
        name: "LABEL".into(),
        format: BinColumnType::Char(8),
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };

    let columns = vec![col_j, col_d, col_a];
    let row_len: usize = columns.iter().map(|c| c.format.byte_size()).sum();
    assert_eq!(row_len, 4 + 8 + 8); // 20

    let nrows = 2;
    let mut main_data = vec![0u8; row_len * nrows];

    // Row 0: COUNT=42, VALUE=3.14, LABEL="hello   "
    main_data[0..4].copy_from_slice(&42i32.to_be_bytes());
    main_data[4..12].copy_from_slice(&3.125f64.to_be_bytes());
    main_data[12..20].copy_from_slice(b"hello   ");

    // Row 1: COUNT=-1, VALUE=2.718, LABEL="world   "
    main_data[20..24].copy_from_slice(&(-1i32).to_be_bytes());
    main_data[24..32].copy_from_slice(&2.5f64.to_be_bytes());
    main_data[32..40].copy_from_slice(b"world   ");

    let table = BinTable {
        columns,
        nrows,
        row_len,
        main_data,
        heap: Vec::new(),
    };

    let mut fits = FitsFile::with_empty_primary();
    fits.push_extension(Hdu::bintable_extension(table));

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.len(), 2);

    if let HduData::BinTable(t) = &fits2.extensions()[0].data {
        assert_eq!(t.nrows, 2);
        assert_eq!(t.columns.len(), 3);

        // Check row 0
        if let BinCellValue::I32(v) = t.get_cell(0, 0).unwrap() {
            assert_eq!(v, vec![42]);
        } else {
            panic!("wrong cell type");
        }
        if let BinCellValue::F64(v) = t.get_cell(0, 1).unwrap() {
            assert!((v[0] - 3.125).abs() < 1e-10);
        } else {
            panic!("wrong cell type");
        }
        if let BinCellValue::String(s) = t.get_cell(0, 2).unwrap() {
            assert_eq!(s, "hello");
        } else {
            panic!("wrong cell type");
        }

        // Check row 1
        if let BinCellValue::I32(v) = t.get_cell(1, 0).unwrap() {
            assert_eq!(v, vec![-1]);
        } else {
            panic!("wrong cell type");
        }
        if let BinCellValue::String(s) = t.get_cell(1, 2).unwrap() {
            assert_eq!(s, "world");
        } else {
            panic!("wrong cell type");
        }
    } else {
        panic!("expected bintable");
    }
}

#[test]
fn round_trip_bintable_vla() {
    use fitskit::bintable::{BinCellValue, BinColumn, BinColumnType};

    // Variable-length array column using P descriptor
    let col = BinColumn {
        name: "DATA".into(),
        format: BinColumnType::VarP('J'),
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };

    let row_len = col.format.byte_size(); // 8 bytes per P descriptor
    let nrows = 2;
    let mut main_data = vec![0u8; row_len * nrows];

    // Heap: row 0 has 3 i32s at offset 0, row 1 has 2 i32s at offset 12
    let mut heap = Vec::new();
    // Row 0 data: [10, 20, 30]
    heap.extend_from_slice(&10i32.to_be_bytes());
    heap.extend_from_slice(&20i32.to_be_bytes());
    heap.extend_from_slice(&30i32.to_be_bytes());
    // Row 1 data: [40, 50]
    heap.extend_from_slice(&40i32.to_be_bytes());
    heap.extend_from_slice(&50i32.to_be_bytes());

    // Row 0 descriptor: count=3, offset=0
    main_data[0..4].copy_from_slice(&3i32.to_be_bytes());
    main_data[4..8].copy_from_slice(&0i32.to_be_bytes());
    // Row 1 descriptor: count=2, offset=12
    main_data[8..12].copy_from_slice(&2i32.to_be_bytes());
    main_data[12..16].copy_from_slice(&12i32.to_be_bytes());

    let table = BinTable {
        columns: vec![col],
        nrows,
        row_len,
        main_data,
        heap,
    };

    let mut fits = FitsFile::with_empty_primary();
    fits.push_extension(Hdu::bintable_extension(table));

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    if let HduData::BinTable(t) = &fits2.extensions()[0].data {
        if let BinCellValue::I32(v) = t.get_cell(0, 0).unwrap() {
            assert_eq!(v, vec![10, 20, 30]);
        } else {
            panic!("wrong type");
        }
        if let BinCellValue::I32(v) = t.get_cell(1, 0).unwrap() {
            assert_eq!(v, vec![40, 50]);
        } else {
            panic!("wrong type");
        }
    } else {
        panic!("expected bintable");
    }
}

#[test]
fn round_trip_ascii_table() {
    use fitskit::ascii_table::{AsciiColumn, AsciiFormat, AsciiTable};

    // 2 rows, 2 columns: name (A10) and value (F10.3)
    let col_name = AsciiColumn {
        name: "NAME".into(),
        format: AsciiFormat::Aw(10),
        tbcol: 1,
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };
    let col_val = AsciiColumn {
        name: "VALUE".into(),
        format: AsciiFormat::Fwd(10, 3),
        tbcol: 11,
        tscal: 1.0,
        tzero: 0.0,
        tunit: None,
    };

    let row_len = 20;
    let nrows = 2;
    let mut raw = vec![b' '; row_len * nrows];

    // Row 0: "alpha     " + "    12.345"
    raw[0..10].copy_from_slice(b"alpha     ");
    raw[10..20].copy_from_slice(b"    12.345");

    // Row 1: "beta      " + "   -99.000"
    raw[20..30].copy_from_slice(b"beta      ");
    raw[30..40].copy_from_slice(b"   -99.000");

    let table = AsciiTable {
        columns: vec![col_name, col_val],
        nrows,
        row_len,
        raw_data: raw,
    };

    let mut fits = FitsFile::with_empty_primary();
    fits.push_extension(Hdu::ascii_table_extension(table));

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.len(), 2);

    if let HduData::AsciiTable(t) = &fits2.extensions()[0].data {
        assert_eq!(t.nrows, 2);
        assert_eq!(t.get_string(0, 0).unwrap(), "alpha");
        assert!((t.get_float(0, 1).unwrap() - 12.345).abs() < 1e-6);
        assert_eq!(t.get_string(1, 0).unwrap(), "beta");
        assert!((t.get_float(1, 1).unwrap() - (-99.0)).abs() < 1e-6);
    } else {
        panic!("expected ascii table");
    }
}

#[test]
fn multiple_extensions() {
    let mut fits = FitsFile::with_empty_primary();

    // Add image extension
    let img = ImageData::new(vec![2], PixelData::F32(vec![1.0, 2.0]));
    let mut img_hdu = Hdu::image_extension(img);
    img_hdu
        .header
        .set("EXTNAME", HeaderValue::String("MYIMAGE".into()), None);
    fits.push_extension(img_hdu);

    // Add another image extension
    let img2 = ImageData::new(vec![3], PixelData::U8(vec![10, 20, 30]));
    let mut img_hdu2 = Hdu::image_extension(img2);
    img_hdu2
        .header
        .set("EXTNAME", HeaderValue::String("OTHER".into()), None);
    fits.push_extension(img_hdu2);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.len(), 3);
    assert!(fits2.find_extension("MYIMAGE").is_some());
    assert!(fits2.find_extension("OTHER").is_some());
    assert!(fits2.find_extension("NOTHERE").is_none());
}

#[test]
fn block_alignment() {
    // Every FITS file should be a multiple of 2880 bytes
    let fits = FitsFile::with_empty_primary();
    let bytes = fits.to_bytes().unwrap();
    assert_eq!(bytes.len() % 2880, 0);

    let img = ImageData::new(vec![100, 100], PixelData::I16(vec![0; 10000]));
    let fits = FitsFile::with_primary_image(img);
    let bytes = fits.to_bytes().unwrap();
    assert_eq!(bytes.len() % 2880, 0);
}

#[test]
fn header_keyword_preservation() {
    let mut fits = FitsFile::with_empty_primary();
    fits.primary_mut()
        .header
        .push(Keyword::commentary("COMMENT", "Test comment line"));
    fits.primary_mut()
        .header
        .push(Keyword::commentary("HISTORY", "Created by fitskit tests"));
    fits.primary_mut()
        .header
        .set("AUTHOR", HeaderValue::String("test".into()), None);

    let bytes = fits.to_bytes().unwrap();
    let fits2 = FitsFile::from_bytes(&bytes).unwrap();

    assert_eq!(fits2.primary().header.get_string("AUTHOR"), Some("test"));

    let comments: Vec<_> = fits2
        .primary()
        .header
        .iter()
        .filter(|k| k.name == "COMMENT")
        .collect();
    assert_eq!(comments.len(), 1);

    let history: Vec<_> = fits2
        .primary()
        .header
        .iter()
        .filter(|k| k.name == "HISTORY")
        .collect();
    assert_eq!(history.len(), 1);
}
