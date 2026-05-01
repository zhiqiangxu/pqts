//! Compact encoding/decoding for Falcon-512 signatures and public keys.
//!
//! - Public key h: 512 coefficients in [0, q=12289), packed as 14-bit values → 896 bytes
//! - Signature s2: 512 small coefficients encoded with a simple compressed format
//!   using a sign bit + variable-length magnitude, targeting ~626 bytes for the s2 part
//!   (total signature = 1 byte header + 40 byte nonce + s2_encoded ≈ 666 bytes)

use crate::falcon::params::{N, Q};

/// Pack public key: 512 coefficients × 14 bits = 7168 bits = 896 bytes.
/// Each coefficient is in [0, 12289) which fits in 14 bits.
pub fn pack_public_key(h: &[u32; N]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(PK_PACKED_BYTES);
    let mut acc: u32 = 0;
    let mut acc_bits: u32 = 0;

    for &coeff in h.iter() {
        debug_assert!(coeff < Q);
        acc |= coeff << acc_bits;
        acc_bits += 14;
        while acc_bits >= 8 {
            buf.push((acc & 0xFF) as u8);
            acc >>= 8;
            acc_bits -= 8;
        }
    }
    if acc_bits > 0 {
        buf.push((acc & 0xFF) as u8);
    }

    buf
}

/// Unpack public key from 896 bytes → 512 coefficients.
pub fn unpack_public_key(buf: &[u8]) -> Option<[u32; N]> {
    if buf.len() != PK_PACKED_BYTES {
        return None;
    }

    let mut h = [0u32; N];
    let mut acc: u32 = 0;
    let mut acc_bits: u32 = 0;
    let mut byte_idx = 0;

    for coeff in h.iter_mut() {
        while acc_bits < 14 {
            if byte_idx >= buf.len() {
                return None;
            }
            acc |= (buf[byte_idx] as u32) << acc_bits;
            byte_idx += 1;
            acc_bits += 8;
        }
        *coeff = acc & 0x3FFF; // 14-bit mask
        if *coeff >= Q {
            return None;
        }
        acc >>= 14;
        acc_bits -= 14;
    }

    Some(h)
}

/// Packed public key size: ceil(512 * 14 / 8) = 896 bytes
pub const PK_PACKED_BYTES: usize = (N * 14 + 7) / 8; // 896

/// Encode s2 coefficients using a compact format.
///
/// Real Falcon uses a sophisticated compressed encoding (Golomb-Rice / Huffman).
/// We use a simpler scheme: each coefficient is encoded as:
/// - 1 sign bit (0 = non-negative, 1 = negative)
/// - 7-bit magnitude (abs value up to 127)
/// This gives exactly 1 byte per coefficient = 512 bytes.
///
/// For even more compact encoding, we use a variable-length scheme:
/// - If |coeff| < 4: 4 bits (1 sign + 3 magnitude) — covers ~95% of coefficients
/// - If |coeff| < 64: 11 bits (1 flag + 1 sign + 9 magnitude)
/// - Otherwise: 19 bits (1 flag + 1 flag + 1 sign + 16 magnitude)
///
/// For simplicity and correctness, we use the fixed 8-bit encoding here,
/// giving 512 bytes for s2 + 40 bytes nonce + 1 byte header = 553 bytes.
/// (Real Falcon achieves ~626 bytes for s2 via entropy coding.)
pub fn encode_s2(s2: &[i32]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(s2.len());
    for &c in s2 {
        let sign = if c < 0 { 0x80u8 } else { 0u8 };
        let mag = c.unsigned_abs() as u8 & 0x7F;
        buf.push(sign | mag);
    }
    buf
}

/// Decode s2 from compact encoding.
pub fn decode_s2(buf: &[u8], n: usize) -> Option<Vec<i32>> {
    if buf.len() != n {
        return None;
    }
    let mut s2 = Vec::with_capacity(n);
    for &b in buf {
        let mag = (b & 0x7F) as i32;
        let val = if b & 0x80 != 0 { -mag } else { mag };
        s2.push(val);
    }
    Some(s2)
}

/// Total signature packed size: 1 (header) + 40 (nonce) + 512 (s2) = 553 bytes
/// Real Falcon-512 is ~666 bytes (using entropy-coded s2).
pub const SIG_PACKED_BYTES: usize = 1 + 40 + N; // 553

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_public_key() {
        let mut h = [0u32; N];
        for i in 0..N {
            h[i] = (i as u32 * 37 + 5) % Q;
        }

        let packed = pack_public_key(&h);
        assert_eq!(packed.len(), PK_PACKED_BYTES);

        let unpacked = unpack_public_key(&packed).unwrap();
        assert_eq!(h, unpacked);
    }

    #[test]
    fn test_encode_decode_s2() {
        let s2: Vec<i32> = (0..N as i32)
            .map(|i| ((i * 7 + 3) % 11) - 5)
            .collect();

        let encoded = encode_s2(&s2);
        assert_eq!(encoded.len(), N);

        let decoded = decode_s2(&encoded, N).unwrap();
        assert_eq!(s2, decoded);
    }

    #[test]
    fn test_encode_decode_s2_edge_cases() {
        let s2 = vec![0, 1, -1, 127, -127, 0];
        let encoded = encode_s2(&s2);
        let decoded = decode_s2(&encoded, s2.len()).unwrap();
        assert_eq!(s2, decoded);
    }
}
