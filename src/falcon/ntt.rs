//! Number Theoretic Transform (NTT) over Z_q[X]/(X^512 + 1) with q = 12289.
//!
//! Uses the negacyclic NTT: multiply by psi^i, apply standard NTT, and reverse.
//! PSI is a primitive 1024th root of unity mod q (PSI^512 = -1).

use crate::falcon::params::{N, PSI, Q, N_INV};

/// Modular exponentiation: base^exp mod q
fn mod_pow(mut base: u32, mut exp: u32, modulus: u32) -> u32 {
    let mut result: u64 = 1;
    base %= modulus;
    let m = modulus as u64;
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base as u64) % m;
        }
        exp >>= 1;
        base = ((base as u64 * base as u64) % m) as u32;
    }
    result as u32
}

/// Modular inverse via Fermat's little theorem: a^(q-2) mod q
pub fn mod_inv(a: u32, q: u32) -> u32 {
    mod_pow(a, q - 2, q)
}

/// Bit-reverse an index of `log_n` bits.
fn bit_reverse(x: usize, log_n: u32) -> usize {
    let mut result = 0usize;
    let mut val = x;
    for _ in 0..log_n {
        result = (result << 1) | (val & 1);
        val >>= 1;
    }
    result
}

/// Bit-reverse permutation of an array.
fn bit_reverse_copy(a: &[u32; N]) -> [u32; N] {
    let mut out = [0u32; N];
    for i in 0..N {
        out[bit_reverse(i, 9)] = a[i];
    }
    out
}

/// Forward negacyclic NTT in-place.
/// Computes the NTT of a polynomial in Z_q[X]/(X^N+1).
pub fn ntt(a: &mut [u32; N]) {
    let q64 = Q as u64;
    let omega = ((PSI as u64 * PSI as u64) % q64) as u32; // primitive N-th root of unity

    // Step 1: Multiply by psi^i (negacyclic twist)
    let mut psi_power: u64 = 1;
    for i in 0..N {
        a[i] = ((a[i] as u64 * psi_power) % q64) as u32;
        psi_power = (psi_power * PSI as u64) % q64;
    }

    // Step 2: Bit-reverse permutation
    *a = bit_reverse_copy(a);

    // Step 3: Cooley-Tukey DIT butterfly
    let mut m = 2usize;
    while m <= N {
        let half_m = m >> 1;
        let w_m = mod_pow(omega, (N / m) as u32, Q);
        let mut k = 0;
        while k < N {
            let mut w: u64 = 1;
            for j in 0..half_m {
                let u = a[k + j] as u64;
                let v = (a[k + j + half_m] as u64 * w) % q64;
                a[k + j] = ((u + v) % q64) as u32;
                a[k + j + half_m] = ((u + q64 - v) % q64) as u32;
                w = (w * w_m as u64) % q64;
            }
            k += m;
        }
        m <<= 1;
    }
}

/// Inverse negacyclic NTT in-place.
pub fn inv_ntt(a: &mut [u32; N]) {
    let q64 = Q as u64;
    let omega = ((PSI as u64 * PSI as u64) % q64) as u32;
    let omega_inv = mod_inv(omega, Q);

    // Step 1: Inverse NTT = NTT with omega^{-1}, then scale by N^{-1}
    *a = bit_reverse_copy(a);

    let mut m = 2usize;
    while m <= N {
        let half_m = m >> 1;
        let w_m = mod_pow(omega_inv, (N / m) as u32, Q);
        let mut k = 0;
        while k < N {
            let mut w: u64 = 1;
            for j in 0..half_m {
                let u = a[k + j] as u64;
                let v = (a[k + j + half_m] as u64 * w) % q64;
                a[k + j] = ((u + v) % q64) as u32;
                a[k + j + half_m] = ((u + q64 - v) % q64) as u32;
                w = (w * w_m as u64) % q64;
            }
            k += m;
        }
        m <<= 1;
    }

    // Step 2: Scale by N^{-1}
    let n_inv = N_INV as u64;
    for i in 0..N {
        a[i] = ((a[i] as u64 * n_inv) % q64) as u32;
    }

    // Step 3: Multiply by psi^{-i} (undo negacyclic twist)
    let psi_inv = mod_inv(PSI, Q) as u64;
    let mut psi_inv_power: u64 = 1;
    for i in 0..N {
        a[i] = ((a[i] as u64 * psi_inv_power) % q64) as u32;
        psi_inv_power = (psi_inv_power * psi_inv) % q64;
    }
}

/// Pointwise multiplication in NTT domain: c[i] = a[i] * b[i] mod q.
pub fn ntt_mul(a: &[u32; N], b: &[u32; N]) -> [u32; N] {
    let mut c = [0u32; N];
    let q64 = Q as u64;
    for i in 0..N {
        c[i] = ((a[i] as u64 * b[i] as u64) % q64) as u32;
    }
    c
}

/// Pointwise inverse in NTT domain: out[i] = 1/a[i] mod q.
/// Returns None if any element is zero.
pub fn ntt_inv(a: &[u32; N]) -> Option<[u32; N]> {
    let mut out = [0u32; N];
    for i in 0..N {
        if a[i] == 0 {
            return None;
        }
        out[i] = mod_inv(a[i], Q);
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_roundtrip() {
        let mut a = [0u32; N];
        a[0] = 1;
        a[1] = 2;
        a[2] = 3;
        let original = a;

        ntt(&mut a);
        inv_ntt(&mut a);

        for i in 0..N {
            assert_eq!(a[i], original[i], "mismatch at index {}", i);
        }
    }

    #[test]
    fn test_ntt_multiplication() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2 in Z_q[X]/(X^N+1)
        let mut a = [0u32; N];
        let mut b = [0u32; N];
        a[0] = 1; a[1] = 1;
        b[0] = 1; b[1] = 1;

        ntt(&mut a);
        ntt(&mut b);
        let mut c = ntt_mul(&a, &b);
        inv_ntt(&mut c);

        assert_eq!(c[0], 1, "c[0]");
        assert_eq!(c[1], 2, "c[1]");
        assert_eq!(c[2], 1, "c[2]");
        for i in 3..N {
            assert_eq!(c[i], 0, "nonzero at index {}", i);
        }
    }

    #[test]
    fn test_mod_inv() {
        let a = 3u32;
        let inv = mod_inv(a, Q);
        assert_eq!((a as u64 * inv as u64) % Q as u64, 1);
    }
}
