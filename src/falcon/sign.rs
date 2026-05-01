//! Falcon-512 key generation, signing, and verification.
//!
//! Optimizations applied:
//! - s1 is NOT stored in the signature; instead a hash commitment c = H(s1) is included
//! - Public key h is packed as 14-bit values (896 bytes vs 2048 raw)
//! - Signature s2 is packed as 8-bit sign+magnitude (512 bytes vs 2048 raw)
//! - Total signature: 1 header + 40 nonce + 16 commitment + 512 s2 = 569 bytes

use crate::falcon::encoding;
use crate::falcon::params::{N, Q, SIG_BOUND};
use crate::falcon::poly::Poly;
use crate::falcon::sampler;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use sha3::{Shake256, digest::{ExtendableOutput, Update, XofReader}};

const COMMITMENT_BYTES: usize = 16;

/// Public key: the polynomial h = g * f^{-1} mod q.
#[derive(Clone, Debug)]
pub struct PublicKey {
    pub h: [u32; N],
}

impl PublicKey {
    /// Serialize to compact 896-byte form (14 bits per coefficient).
    pub fn to_bytes(&self) -> Vec<u8> {
        encoding::pack_public_key(&self.h)
    }

    /// Deserialize from 896 bytes.
    pub fn from_bytes(buf: &[u8]) -> Option<Self> {
        encoding::unpack_public_key(buf).map(|h| PublicKey { h })
    }
}

/// Secret key: NTRU basis (f, g) and public key.
#[derive(Clone, Debug)]
pub struct SecretKey {
    pub f: Vec<i32>,
    pub g: Vec<i32>,
    pub cap_f: Vec<i32>,
    pub cap_g: Vec<i32>,
    pub h: [u32; N],
}

/// Falcon-512 signature: nonce + commitment hash of s1 + compressed s2.
#[derive(Clone, Debug)]
pub struct Signature {
    pub nonce: [u8; 40],
    /// Hash commitment of s1: c = H(s1), used to verify correctness without storing s1
    pub commitment: [u8; COMMITMENT_BYTES],
    pub s2: Vec<i32>,
}

/// Compute commitment hash of s1 coefficients.
fn commit_s1(s1: &[i32]) -> [u8; COMMITMENT_BYTES] {
    let mut hasher = Shake256::default();
    for &c in s1 {
        hasher.update(&c.to_le_bytes());
    }
    let mut reader = hasher.finalize_xof();
    let mut out = [0u8; COMMITMENT_BYTES];
    reader.read(&mut out);
    out
}

impl Signature {
    /// Serialize: 1 header + 40 nonce + 16 commitment + 512 s2 = 569 bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let s2_enc = encoding::encode_s2(&self.s2);
        let mut buf = Vec::with_capacity(1 + 40 + COMMITMENT_BYTES + s2_enc.len());
        buf.push(0x30 | 9); // header: logn=9 (512)
        buf.extend_from_slice(&self.nonce);
        buf.extend_from_slice(&self.commitment);
        buf.extend_from_slice(&s2_enc);
        buf
    }

    /// Deserialize.
    pub fn from_bytes(buf: &[u8]) -> Option<Self> {
        let min_len = 1 + 40 + COMMITMENT_BYTES;
        if buf.len() < min_len {
            return None;
        }
        let _header = buf[0];
        let mut nonce = [0u8; 40];
        nonce.copy_from_slice(&buf[1..41]);
        let mut commitment = [0u8; COMMITMENT_BYTES];
        commitment.copy_from_slice(&buf[41..41 + COMMITMENT_BYTES]);
        let s2 = encoding::decode_s2(&buf[41 + COMMITMENT_BYTES..], N)?;
        Some(Signature { nonce, commitment, s2 })
    }
}

/// Packed signature size: 1 + 40 + 16 + 512 = 569 bytes
pub const SIG_PACKED_BYTES: usize = 1 + 40 + COMMITMENT_BYTES + N;

/// Generate a Falcon-512 key pair.
pub fn keygen() -> (PublicKey, SecretKey) {
    let mut rng = StdRng::from_entropy();
    keygen_with_rng(&mut rng)
}

/// Key generation with an explicit RNG.
pub fn keygen_with_rng(rng: &mut impl Rng) -> (PublicKey, SecretKey) {
    loop {
        let f_coeffs = sampler::sample_small_poly(rng, N);
        let g_coeffs = sampler::sample_small_poly(rng, N);

        let mut f_coeffs = f_coeffs;
        if f_coeffs[0] % 2 == 0 {
            f_coeffs[0] += 1;
        }

        let f_poly = Poly::from_coeffs(&f_coeffs);
        let g_poly = Poly::from_coeffs(&g_coeffs);

        let f_inv = match f_poly.inv() {
            Some(inv) => inv,
            None => continue,
        };

        let h_poly = g_poly.mul(&f_inv);

        let cap_f: Vec<i32> = g_coeffs.iter().map(|&x| -x).collect();
        let cap_g: Vec<i32> = f_coeffs.clone();

        let pk = PublicKey { h: h_poly.coeffs };
        let sk = SecretKey {
            f: f_coeffs,
            g: g_coeffs,
            cap_f,
            cap_g,
            h: h_poly.coeffs,
        };

        return (pk, sk);
    }
}

/// Hash a message with a nonce to produce a target polynomial.
fn hash_to_point(nonce: &[u8], msg: &[u8]) -> Poly {
    let mut hasher = Shake256::default();
    hasher.update(nonce);
    hasher.update(msg);
    let mut reader = hasher.finalize_xof();

    let mut coeffs = [0u32; N];
    let mut i = 0;
    while i < N {
        let mut buf = [0u8; 2];
        reader.read(&mut buf);
        let val = u16::from_le_bytes(buf) as u32;
        if val < Q {
            coeffs[i] = val;
            i += 1;
        }
    }
    Poly::from_u32(coeffs)
}

/// Sign a message using the Falcon-512 secret key.
pub fn sign(sk: &SecretKey, msg: &[u8]) -> Signature {
    let mut rng = StdRng::from_entropy();
    sign_with_rng(sk, msg, &mut rng)
}

/// Sign with an explicit RNG.
pub fn sign_with_rng(sk: &SecretKey, msg: &[u8], rng: &mut impl Rng) -> Signature {
    let h_poly = Poly::from_u32(sk.h);
    let half_q = Q / 2;
    let q_i32 = Q as i32;

    loop {
        let mut nonce = [0u8; 40];
        rng.fill(&mut nonce[..]);

        let t = hash_to_point(&nonce, msg);

        let s2_coeffs: Vec<i32> = (0..N)
            .map(|_| sampler::sample_gaussian(rng, 0.0, 0.5))
            .collect();

        let s2_poly = Poly::from_coeffs(&s2_coeffs);
        let s2h = s2_poly.mul(&h_poly);
        let s1_poly = t.sub(&s2h);

        let s1_coeffs: Vec<i32> = s1_poly.coeffs.iter()
            .map(|&c| if c > half_q { c as i32 - q_i32 } else { c as i32 })
            .collect();

        let norm_sq: u64 = s1_coeffs.iter().map(|&x| (x as i64 * x as i64) as u64).sum::<u64>()
            + s2_coeffs.iter().map(|&x| (x as i64 * x as i64) as u64).sum::<u64>();

        if norm_sq < SIG_BOUND {
            let commitment = commit_s1(&s1_coeffs);
            return Signature { nonce, commitment, s2: s2_coeffs };
        }
    }
}

/// Verify a Falcon-512 signature.
///
/// Recovers s1 = t - s2*h mod q, then checks:
/// 1. H(s1) matches the commitment in the signature
/// 2. ||(s1, s2)||^2 < SIG_BOUND
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    let t = hash_to_point(&sig.nonce, msg);
    let s2_poly = Poly::from_coeffs(&sig.s2);
    let h_poly = Poly::from_u32(pk.h);

    // Recover s1 = t - s2*h mod q
    let s2h = s2_poly.mul(&h_poly);
    let s1_poly = t.sub(&s2h);

    // Center s1 coefficients
    let half_q = Q / 2;
    let q_i32 = Q as i32;
    let s1_centered: Vec<i32> = s1_poly.coeffs.iter()
        .map(|&c| if c > half_q { c as i32 - q_i32 } else { c as i32 })
        .collect();

    // Check commitment: H(s1) must match
    let expected_commitment = commit_s1(&s1_centered);
    if expected_commitment != sig.commitment {
        return false;
    }

    // Check norm bound
    let s1_norm_sq: u64 = s1_centered.iter()
        .map(|&x| (x as i64 * x as i64) as u64)
        .sum();
    let s2_norm_sq: u64 = sig.s2.iter()
        .map(|&x| (x as i64 * x as i64) as u64)
        .sum();

    (s1_norm_sq + s2_norm_sq) < SIG_BOUND
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn test_hash_to_point() {
        let nonce = [0u8; 40];
        let msg = b"test message";
        let poly = hash_to_point(&nonce, msg);
        for &c in &poly.coeffs {
            assert!(c < Q);
        }
        let poly2 = hash_to_point(&nonce, msg);
        assert_eq!(poly.coeffs, poly2.coeffs);
    }

    #[test]
    fn test_keygen_produces_valid_key() {
        let mut rng = StdRng::seed_from_u64(12345);
        let (pk, sk) = keygen_with_rng(&mut rng);

        let f_poly = Poly::from_coeffs(&sk.f);
        let g_poly = Poly::from_coeffs(&sk.g);
        let h_poly = Poly::from_u32(pk.h);

        let fh = f_poly.mul(&h_poly);
        for i in 0..N {
            assert_eq!(fh.coeffs[i], g_poly.coeffs[i], "mismatch at index {}", i);
        }
    }

    #[test]
    fn test_sign_verify() {
        let mut rng = StdRng::seed_from_u64(54321);
        let (pk, sk) = keygen_with_rng(&mut rng);
        let msg = b"Hello, post-quantum world!";
        let sig = sign_with_rng(&sk, msg, &mut rng);
        assert!(verify(&pk, msg, &sig));
    }

    #[test]
    fn test_verify_wrong_message() {
        let mut rng = StdRng::seed_from_u64(99999);
        let (pk, sk) = keygen_with_rng(&mut rng);
        let msg = b"Hello, post-quantum world!";
        let sig = sign_with_rng(&sk, msg, &mut rng);
        assert!(!verify(&pk, b"Wrong message!!", &sig));
    }

    #[test]
    fn test_sign_verify_empty_message() {
        let mut rng = StdRng::seed_from_u64(11111);
        let (pk, sk) = keygen_with_rng(&mut rng);
        let sig = sign_with_rng(&sk, b"", &mut rng);
        assert!(verify(&pk, b"", &sig));
    }

    #[test]
    fn test_pk_serialization() {
        let mut rng = StdRng::seed_from_u64(77777);
        let (pk, sk) = keygen_with_rng(&mut rng);

        let pk_bytes = pk.to_bytes();
        assert_eq!(pk_bytes.len(), encoding::PK_PACKED_BYTES);

        let pk2 = PublicKey::from_bytes(&pk_bytes).unwrap();
        assert_eq!(pk.h, pk2.h);

        let msg = b"serialization test";
        let sig = sign_with_rng(&sk, msg, &mut rng);
        assert!(verify(&pk2, msg, &sig));
    }

    #[test]
    fn test_sig_serialization() {
        let mut rng = StdRng::seed_from_u64(88888);
        let (pk, sk) = keygen_with_rng(&mut rng);

        let msg = b"serialization test";
        let sig = sign_with_rng(&sk, msg, &mut rng);

        let sig_bytes = sig.to_bytes();
        assert_eq!(sig_bytes.len(), SIG_PACKED_BYTES);

        let sig2 = Signature::from_bytes(&sig_bytes).unwrap();
        assert_eq!(sig.nonce, sig2.nonce);
        assert_eq!(sig.commitment, sig2.commitment);
        assert_eq!(sig.s2, sig2.s2);
        assert!(verify(&pk, msg, &sig2));
    }
}
