//! FITS checksum computation per the FITS Checksum Proposal (Seaman, Pence, Rots 2012).
//!
//! Uses ones-complement arithmetic on 2880-byte blocks. The DATASUM keyword stores
//! the data checksum as a decimal string. The CHECKSUM keyword stores an ASCII-encoded
//! complement such that the checksum of the entire HDU (header + data) is zero.

use crate::error::{Error, Result};
use crate::header::Header;
use crate::keyword::HeaderValue;

/// Compute the ones-complement checksum of a byte buffer.
///
/// Treats the data as a sequence of 32-bit big-endian words summed with
/// end-around carry (equivalently, 16-bit hi/lo halves with cross carries as
/// in the checksum proposal). A trailing partial word is zero-padded on the
/// right, which is exactly what the FITS block padding would contribute.
///
/// The sum is accumulated in 64-bit and folded periodically, so it is exact
/// for arbitrarily large buffers.
pub fn checksum(data: &[u8]) -> u32 {
    let mut acc = Checksum::new();
    acc.update(data);
    acc.finish()
}

/// Incremental FITS checksum over a byte stream.
///
/// Feed bytes in any chunking with [`update`](Self::update); the result of
/// [`finish`](Self::finish) is identical to [`checksum`] over the concatenated
/// input. Partial 32-bit words at chunk boundaries are carried across calls.
///
/// ```
/// use fitskit::checksum::{checksum, Checksum};
///
/// let data: Vec<u8> = (0..10_000u32).map(|i| (i * 7) as u8).collect();
/// let mut acc = Checksum::new();
/// for chunk in data.chunks(333) {
///     acc.update(chunk);
/// }
/// assert_eq!(acc.finish(), checksum(&data));
/// ```
#[derive(Debug, Clone, Default)]
pub struct Checksum {
    /// Running sum, kept below 2^32 between chunks by [`fold32`].
    sum: u64,
    /// Bytes of an incomplete trailing word from the previous `update`.
    pending: [u8; 4],
    npending: usize,
}

impl Checksum {
    /// Fold at least every 2^32 words so the u64 accumulator cannot overflow.
    /// A 1 MiB chunk (a multiple of 4) keeps the folds negligible and the inner
    /// loop a plain vectorizable reduction.
    const CHUNK: usize = 1 << 20;

    pub fn new() -> Self {
        Self::default()
    }

    /// Add `data` to the running sum.
    pub fn update(&mut self, mut data: &[u8]) {
        if self.npending > 0 {
            let take = (4 - self.npending).min(data.len());
            self.pending[self.npending..self.npending + take].copy_from_slice(&data[..take]);
            self.npending += take;
            data = &data[take..];
            if self.npending < 4 {
                return;
            }
            self.sum = fold32(self.sum + u32::from_be_bytes(self.pending) as u64);
            self.npending = 0;
        }

        for chunk in data.chunks(Self::CHUNK) {
            let (words, rem) = chunk.as_chunks::<4>();
            let chunk_sum: u64 = words.iter().map(|&w| u32::from_be_bytes(w) as u64).sum();
            self.sum = fold32(self.sum + chunk_sum);
            // Only the final chunk can have a remainder (CHUNK is a multiple of 4).
            self.pending[..rem.len()].copy_from_slice(rem);
            self.npending = rem.len();
        }
    }

    /// Finish the sum, zero-padding any trailing partial word.
    pub fn finish(&self) -> u32 {
        let mut sum = self.sum;
        if self.npending > 0 {
            let mut word = [0u8; 4];
            word[..self.npending].copy_from_slice(&self.pending[..self.npending]);
            sum = fold32(sum + u32::from_be_bytes(word) as u64);
        }
        sum as u32
    }
}

/// Reduce a 64-bit sum to a 32-bit ones-complement value by end-around carry.
///
/// The result is `sum mod (2^32 - 1)` in the canonical ones-complement form
/// where a non-zero multiple of `2^32 - 1` maps to `0xFFFF_FFFF` rather than
/// `0` (matching the hi/lo fold in the FITS checksum proposal and cfitsio).
fn fold32(mut sum: u64) -> u64 {
    while sum >> 32 != 0 {
        sum = (sum & 0xFFFF_FFFF) + (sum >> 32);
    }
    sum
}

/// Fold end-around carries from the hi/lo accumulators until both are 16-bit.
fn fold_carries(mut hi: u32, mut lo: u32) -> (u32, u32) {
    loop {
        let hicarry = hi >> 16;
        let locarry = lo >> 16;
        if hicarry == 0 && locarry == 0 {
            break;
        }
        hi = (hi & 0xFFFF) + locarry;
        lo = (lo & 0xFFFF) + hicarry;
    }
    (hi, lo)
}

/// Accumulate a checksum: ones-complement add `new_data` checksum into `existing`.
pub fn checksum_accumulate(existing: u32, new_data: &[u8]) -> u32 {
    let new_sum = checksum(new_data);
    ones_complement_add(existing, new_sum)
}

/// Ones-complement addition of two 32-bit values.
fn ones_complement_add(a: u32, b: u32) -> u32 {
    let hi = (a >> 16) + (b >> 16);
    let lo = (a & 0xFFFF) + (b & 0xFFFF);

    let (hi, lo) = fold_carries(hi, lo);
    (hi << 16) | lo
}

/// Encode a 32-bit checksum as a 16-character ASCII string.
///
/// If `complement` is true, encodes `!sum` (used for CHECKSUM keyword).
/// If false, encodes `sum` directly.
///
/// Algorithm per cfitsio reference implementation (Seaman, Pence, Rots).
pub fn encode_checksum(sum: u32, complement: bool) -> String {
    let exclude: [u8; 13] = [
        0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40, // :;<=>?@
        0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60, // [\]^_`
    ];
    let masks: [u32; 4] = [0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff];
    let offset: u32 = 0x30;

    let value = if complement { !sum } else { sum };

    let mut asc = [0u8; 16];

    for i in 0..4 {
        let byte = (value & masks[i]) >> (24 - 8 * i);
        let quotient = byte / 4 + offset;
        let remainder = byte % 4;

        let mut ch = [quotient as u8; 4];
        ch[0] = (quotient + remainder) as u8;

        // Adjust excluded characters by incrementing even, decrementing odd
        let mut check = true;
        while check {
            check = false;
            for &ex in &exclude {
                for j in (0..4).step_by(2) {
                    if ch[j] == ex || ch[j + 1] == ex {
                        ch[j] = ch[j].wrapping_add(1);
                        ch[j + 1] = ch[j + 1].wrapping_sub(1);
                        check = true;
                    }
                }
            }
        }

        // Distribute into columns (transposed layout)
        for j in 0..4 {
            asc[4 * j + i] = ch[j];
        }
    }

    // Rotate by 1 position to the right (cyclic permutation for FITS alignment)
    let mut result = [0u8; 16];
    for i in 0..16 {
        result[i] = asc[(i + 15) % 16];
    }

    String::from_utf8(result.to_vec()).expect("checksum encoding produced non-UTF8")
}

/// Decode a 16-character ASCII checksum string back to a 32-bit value.
///
/// If `complement` is true, returns the complement of the decoded value.
pub fn decode_checksum(ascii: &str, complement: bool) -> u32 {
    let bytes = ascii.as_bytes();
    assert!(bytes.len() >= 16, "checksum string must be 16 characters");

    // Undo rotation
    let mut cbuf = [0u8; 16];
    for i in 0..16 {
        cbuf[i] = bytes[(i + 1) % 16].wrapping_sub(0x30);
    }

    let mut hi: u32 = 0;
    let mut lo: u32 = 0;

    for i in (0..16).step_by(4) {
        hi += ((cbuf[i] as u32) << 8) + cbuf[i + 1] as u32;
        lo += ((cbuf[i + 2] as u32) << 8) + cbuf[i + 3] as u32;
    }

    let (hi, lo) = fold_carries(hi, lo);

    let sum = (hi << 16) | lo;
    if complement {
        !sum
    } else {
        sum
    }
}

/// Compute DATASUM for a data byte buffer (should be block-padded).
pub fn datasum(data: &[u8]) -> u32 {
    checksum(data)
}

/// Verify the CHECKSUM of an HDU given its serialized header and data bytes.
///
/// Returns true if the ones-complement sum of all bytes is 0xFFFFFFFF (all ones).
pub fn verify_hdu(header_bytes: &[u8], data_bytes: &[u8]) -> bool {
    let sum = checksum_accumulate(checksum(header_bytes), data_bytes);
    sum == 0xFFFF_FFFF || sum == 0
}

/// Compute and insert DATASUM and CHECKSUM keywords into a header.
///
/// `data_bytes` should be the block-padded data for this HDU.
/// Returns the serialized header bytes (for use by the caller).
pub fn stamp_hdu(header: &mut Header, data_bytes: &[u8]) -> Result<Vec<u8>> {
    stamp_hdu_with_datasum(header, datasum(data_bytes))
}

/// Insert DATASUM and CHECKSUM keywords into a header given a precomputed
/// data-unit checksum (see [`datasum`] / [`Checksum`]).
///
/// Returns the serialized header bytes (for use by the caller).
pub fn stamp_hdu_with_datasum(header: &mut Header, dsum: u32) -> Result<Vec<u8>> {
    header.set(
        "DATASUM",
        HeaderValue::String(dsum.to_string()),
        Some("data unit checksum"),
    );

    // Set CHECKSUM to zeros initially
    header.set(
        "CHECKSUM",
        HeaderValue::String("0000000000000000".into()),
        Some("HDU checksum"),
    );

    // Serialize header to get its bytes
    let mut header_bytes = Vec::new();
    header.write_to(&mut header_bytes)?;

    // Compute header checksum
    let hsum = checksum(&header_bytes);

    // Total HDU checksum
    let total = ones_complement_add(hsum, dsum);

    // Encode complement so that re-checksum of HDU yields all-ones
    let encoded = encode_checksum(total, true);
    header.set(
        "CHECKSUM",
        HeaderValue::String(encoded),
        Some("HDU checksum"),
    );

    // Re-serialize with final CHECKSUM value
    header_bytes.clear();
    header.write_to(&mut header_bytes)?;

    Ok(header_bytes)
}

/// Verify CHECKSUM/DATASUM on an HDU that was read from a file.
///
/// Reconstructs the header + data bytes and checks the ones-complement sum.
/// Returns Ok(()) if valid or no checksum keywords present, Err on mismatch.
pub fn verify_from_header(header: &Header, data_bytes: &[u8]) -> Result<()> {
    verify_datasum_value(header, datasum(data_bytes))
}

/// Verify a precomputed data-unit checksum against the header's DATASUM
/// keyword. Returns Ok(()) if it matches, if DATASUM is absent/unparseable,
/// or if the stored value is 0.
pub fn verify_datasum_value(header: &Header, computed: u32) -> Result<()> {
    if let Some(stored_datasum_str) = header.get_string("DATASUM") {
        if let Ok(stored) = stored_datasum_str.parse::<u64>() {
            let stored = stored as u32;
            if stored != 0 && computed != stored {
                return Err(Error::ChecksumMismatch {
                    expected: stored,
                    actual: computed,
                });
            }
        }
    }

    // For full CHECKSUM verification, we'd need the original header bytes
    // (before parsing). Since we don't retain those, we skip CHECKSUM verification
    // on read. The DATASUM check above catches data corruption.
    // Full verification is possible via verify_hdu() if the caller has raw bytes.

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io_utils;

    #[test]
    fn checksum_empty() {
        assert_eq!(checksum(&[]), 0);
    }

    #[test]
    fn checksum_zeros() {
        let data = vec![0u8; 2880];
        assert_eq!(checksum(&data), 0);
    }

    #[test]
    fn checksum_chunked_matches_whole() {
        // Pseudo-random bytes, fed through the accumulator in odd-sized pieces
        // so partial words straddle every boundary.
        let mut x: u32 = 0x1234_5678;
        let data: Vec<u8> = (0..100_003)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 17;
                x ^= x << 5;
                x as u8
            })
            .collect();
        let whole = checksum(&data);
        for size in [1usize, 3, 4, 7, 4096, 65_537] {
            let mut acc = Checksum::new();
            for c in data.chunks(size) {
                acc.update(c);
            }
            assert_eq!(acc.finish(), whole, "chunk size {size}");
        }
        // ones_complement_add on 4-aligned pieces agrees too
        let (a, b) = data.split_at(40_000);
        assert_eq!(checksum_accumulate(checksum(a), b), whole);
    }

    #[test]
    fn checksum_large_no_overflow() {
        // 1 MiB of 0xFF: every 16-bit word is 0xFFFF, so the sum of either
        // half is 0xFFFF * 262144 = 0x3FFF_C000_0 which overflows a u32 accumulator.
        // Ones-complement: each word contributes 0xFFFF_FFFF == 0 (mod 2^32-1),
        // and a non-zero multiple maps to all-ones.
        let data = vec![0xFFu8; 1 << 20];
        assert_eq!(checksum(&data), 0xFFFF_FFFF);
        // 1 MiB of 0xAB: 262144 words of 0xABAB_ABAB, reduced mod 2^32-1.
        let data = vec![0xABu8; 1 << 20];
        let expected = ((0xABAB_ABABu64 * 262_144) % 0xFFFF_FFFF) as u32;
        assert_eq!(checksum(&data), expected);
    }

    #[test]
    fn checksum_deterministic() {
        let data: Vec<u8> = (0..2880).map(|i| (i % 256) as u8).collect();
        let c1 = checksum(&data);
        let c2 = checksum(&data);
        assert_eq!(c1, c2);
        assert_ne!(c1, 0);
    }

    #[test]
    fn encode_produces_valid_ascii() {
        for val in [0u32, 1, 0xDEADBEEF, 0xFFFFFFFF, 0x12345678] {
            let s = encode_checksum(val, false);
            assert_eq!(s.len(), 16);
            for b in s.bytes() {
                assert!(
                    b.is_ascii_alphanumeric(),
                    "non-alphanumeric char {b:#x} in encoded checksum"
                );
            }
        }
    }

    #[test]
    fn encode_complement_produces_valid_ascii() {
        for val in [0u32, 1, 0xDEADBEEF, 0xFFFFFFFF] {
            let s = encode_checksum(val, true);
            assert_eq!(s.len(), 16);
            for b in s.bytes() {
                assert!(b.is_ascii_alphanumeric());
            }
        }
    }

    #[test]
    fn encode_decode_round_trip() {
        for val in [0u32, 1, 42, 0xDEADBEEF, 0xFFFFFFFF, 0x12345678, 2503531142] {
            let encoded = encode_checksum(val, false);
            let decoded = decode_checksum(&encoded, false);
            assert_eq!(decoded, val, "round-trip failed for {val:#x}");
        }
    }

    #[test]
    fn encode_decode_complement_round_trip() {
        for val in [0u32, 1, 0xDEADBEEF, 0xFFFFFFFF] {
            let encoded = encode_checksum(val, true);
            let decoded = decode_checksum(&encoded, true);
            assert_eq!(decoded, val, "complement round-trip failed for {val:#x}");
        }
    }

    #[test]
    fn stamp_and_verify() {
        let mut header = Header::new();
        header.set("SIMPLE", HeaderValue::Logical(true), Some("standard"));
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(0), None);

        let data_bytes = vec![0u8; 0]; // no data
        let padded_data = io_utils::pad_to_block(&data_bytes);

        let header_bytes = stamp_hdu(&mut header, &padded_data).unwrap();

        // Verify the HDU
        assert!(verify_hdu(&header_bytes, &padded_data));
    }

    #[test]
    fn stamp_and_verify_with_data() {
        let mut header = Header::new();
        header.set("SIMPLE", HeaderValue::Logical(true), Some("standard"));
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(1), None);
        header.set("NAXIS1", HeaderValue::Integer(100), None);

        let data_bytes: Vec<u8> = (0..100).map(|i| (i * 3) as u8).collect();
        let padded_data = io_utils::pad_to_block(&data_bytes);

        let header_bytes = stamp_hdu(&mut header, &padded_data).unwrap();

        assert!(verify_hdu(&header_bytes, &padded_data));

        // Verify DATASUM is present and correct
        let dsum_str = header.get_string("DATASUM").unwrap();
        let dsum: u32 = dsum_str.parse().unwrap();
        assert_eq!(dsum, datasum(&padded_data));
    }

    #[test]
    fn ones_complement_add_identity() {
        assert_eq!(ones_complement_add(0, 0), 0);
        assert_eq!(ones_complement_add(42, 0), 42);
        assert_eq!(ones_complement_add(0, 42), 42);
    }

    #[test]
    fn ones_complement_add_complement() {
        // x + !x should give all ones
        let x = 0x12345678u32;
        let result = ones_complement_add(x, !x);
        assert_eq!(result, 0xFFFFFFFF);
    }

    #[test]
    fn verify_corruption_detected() {
        let mut header = Header::new();
        header.set("SIMPLE", HeaderValue::Logical(true), Some("standard"));
        header.set("BITPIX", HeaderValue::Integer(8), None);
        header.set("NAXIS", HeaderValue::Integer(1), None);
        header.set("NAXIS1", HeaderValue::Integer(100), None);

        let data_bytes: Vec<u8> = (0..100).collect();
        let padded_data = io_utils::pad_to_block(&data_bytes);

        let header_bytes = stamp_hdu(&mut header, &padded_data).unwrap();

        // Corrupt the data
        let mut corrupted = padded_data.clone();
        corrupted[0] ^= 0xFF;

        // Verification should fail
        assert!(!verify_hdu(&header_bytes, &corrupted));
    }
}
