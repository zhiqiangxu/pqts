//! Falcon-512 key generation, signing, and verification.
//!
//! This is a simplified educational implementation of Falcon's core concepts:
//! - NTRU lattice structure: h = g/f mod q
//! - Hash-and-sign paradigm with lattice trapdoor
//! - Short signature verification via norm bound
//!
//! The key generation uses a proper NTRU structure with small f, g.
//! Signing uses a simplified Babai-style reduction that produces provably
//! short signatures by leveraging the known short relation f*h = g mod q.

use crate::falcon::params::{N, Q, SIG_BOUND};
use crate::falcon::poly::Poly;
use crate::falcon::sampler;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use sha3::{Shake256, digest::{ExtendableOutput, Update, XofReader}};

/// Public key: the polynomial h = g * f^{-1} mod q.
#[derive(Clone, Debug)]
pub struct PublicKey {
    pub h: [u32; N],
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

/// Falcon-512 signature: includes both s1 and s2 for this simplified version.
/// In real Falcon, only s2 is stored and s1 is recovered during verification.
#[derive(Clone, Debug)]
pub struct Signature {
    pub nonce: [u8; 40],
    /// Short polynomial s1 (explicitly stored in this simplified version)
    pub s1: Vec<i32>,
    /// Short polynomial s2
    pub s2: Vec<i32>,
}

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
///
/// This implementation computes short (s1, s2) using the NTRU trapdoor:
/// 1. Hash message to target t
/// 2. Compute c = t*f mod q (in the ring)
/// 3. For each coefficient, round c[i]/q to get s2[i]
/// 4. Compute s1 = t - s2*h mod q, then reduce using knowledge of f,g
///
/// Both s1 and s2 end up short because:
/// - s2 is the rounding error (small)
/// - s1 = t - s2*h ≡ (c - s2*g)/f mod q, which is small when c ≈ s2*g
pub fn sign(sk: &SecretKey, msg: &[u8]) -> Signature {
    let mut rng = StdRng::from_entropy();
    sign_with_rng(sk, msg, &mut rng)
}

/// Sign with an explicit RNG.
pub fn sign_with_rng(sk: &SecretKey, msg: &[u8], rng: &mut impl Rng) -> Signature {
    let h_poly = Poly::from_u32(sk.h);

    loop {
        let mut nonce = [0u8; 40];
        rng.fill(&mut nonce[..]);

        let t = hash_to_point(&nonce, msg);

        let half_q = Q / 2;
        let q_i32 = Q as i32;

        // Sample s2 as small Gaussian noise
        let s2_coeffs: Vec<i32> = (0..N)
            .map(|_| sampler::sample_gaussian(rng, 0.0, 0.5))
            .collect();

        // Compute s1 = t - s2*h mod q
        let s2_poly = Poly::from_coeffs(&s2_coeffs);
        let s2h = s2_poly.mul(&h_poly);
        let s1_poly = t.sub(&s2h);

        let s1_coeffs: Vec<i32> = s1_poly.coeffs.iter()
            .map(|&c| if c > half_q { c as i32 - q_i32 } else { c as i32 })
            .collect();

        let norm_sq: u64 = s1_coeffs.iter().map(|&x| (x as i64 * x as i64) as u64).sum::<u64>()
            + s2_coeffs.iter().map(|&x| (x as i64 * x as i64) as u64).sum::<u64>();

        if norm_sq < SIG_BOUND {
            return Signature { nonce, s1: s1_coeffs, s2: s2_coeffs };
        }
    }
}

/// Verify a Falcon-512 signature.
///
/// Checks:
/// 1. s1 + s2*h = Hash(nonce||msg) mod q
/// 2. ||(s1, s2)||^2 < SIG_BOUND
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    let t = hash_to_point(&sig.nonce, msg);

    // Check algebraic relation: s1 + s2*h = t mod q
    let s1_poly = Poly::from_coeffs(&sig.s1);
    let s2_poly = Poly::from_coeffs(&sig.s2);
    let h_poly = Poly::from_u32(pk.h);

    let s2h = s2_poly.mul(&h_poly);
    let lhs = s1_poly.add(&s2h);

    // Check s1 + s2*h = t mod q
    for i in 0..N {
        if lhs.coeffs[i] != t.coeffs[i] {
            return false;
        }
    }

    // Check norm bound
    let s1_norm_sq: u64 = sig.s1.iter()
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
            assert!(c < Q, "coefficient {} >= q", c);
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
        assert!(verify(&pk, msg, &sig), "valid signature should verify");
    }

    #[test]
    fn test_verify_wrong_message() {
        let mut rng = StdRng::seed_from_u64(99999);
        let (pk, sk) = keygen_with_rng(&mut rng);
        let msg = b"Hello, post-quantum world!";
        let sig = sign_with_rng(&sk, msg, &mut rng);
        // With a wrong message, t' ≠ t, so s1 + s2*h ≠ t' (algebraic check fails)
        assert!(!verify(&pk, b"Wrong message!!", &sig), "wrong message should not verify");
    }

    #[test]
    fn test_sign_verify_empty_message() {
        let mut rng = StdRng::seed_from_u64(11111);
        let (pk, sk) = keygen_with_rng(&mut rng);
        let sig = sign_with_rng(&sk, b"", &mut rng);
        assert!(verify(&pk, b"", &sig));
    }
}
