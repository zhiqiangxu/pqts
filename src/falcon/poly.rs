//! Polynomial operations in Z_q[X]/(X^N + 1).

use crate::falcon::ntt::{inv_ntt, ntt, ntt_inv, ntt_mul};
use crate::falcon::params::{N, Q};

/// A polynomial in Z_q[X]/(X^N+1), stored as coefficients in [0, q).
#[derive(Clone, Debug)]
pub struct Poly {
    pub coeffs: [u32; N],
}

impl Poly {
    /// Zero polynomial.
    pub fn zero() -> Self {
        Poly { coeffs: [0u32; N] }
    }

    /// Create from a coefficient slice. Values are reduced mod q.
    pub fn from_coeffs(c: &[i32]) -> Self {
        let mut coeffs = [0u32; N];
        for (i, &v) in c.iter().enumerate().take(N) {
            coeffs[i] = ((v % Q as i32 + Q as i32) % Q as i32) as u32;
        }
        Poly { coeffs }
    }

    /// Create from a u32 coefficient array (already reduced mod q).
    pub fn from_u32(c: [u32; N]) -> Self {
        Poly { coeffs: c }
    }

    /// Add two polynomials mod q.
    pub fn add(&self, other: &Poly) -> Poly {
        let mut r = [0u32; N];
        for i in 0..N {
            r[i] = (self.coeffs[i] + other.coeffs[i]) % Q;
        }
        Poly { coeffs: r }
    }

    /// Subtract two polynomials mod q: self - other.
    pub fn sub(&self, other: &Poly) -> Poly {
        let mut r = [0u32; N];
        for i in 0..N {
            r[i] = (self.coeffs[i] + Q - other.coeffs[i]) % Q;
        }
        Poly { coeffs: r }
    }

    /// Negate polynomial mod q.
    pub fn neg(&self) -> Poly {
        let mut r = [0u32; N];
        for i in 0..N {
            r[i] = if self.coeffs[i] == 0 { 0 } else { Q - self.coeffs[i] };
        }
        Poly { coeffs: r }
    }

    /// Multiply two polynomials in Z_q[X]/(X^N+1) using NTT.
    pub fn mul(&self, other: &Poly) -> Poly {
        let mut a = self.coeffs;
        let mut b = other.coeffs;
        ntt(&mut a);
        ntt(&mut b);
        let mut c = ntt_mul(&a, &b);
        inv_ntt(&mut c);
        Poly { coeffs: c }
    }

    /// Compute the inverse of this polynomial in Z_q[X]/(X^N+1) using NTT.
    /// Returns None if the polynomial is not invertible.
    pub fn inv(&self) -> Option<Poly> {
        let mut a = self.coeffs;
        ntt(&mut a);
        let a_inv = ntt_inv(&a)?;
        let mut result = a_inv;
        inv_ntt(&mut result);
        Some(Poly { coeffs: result })
    }

    /// Compute the NTT representation and return it.
    pub fn to_ntt(&self) -> [u32; N] {
        let mut a = self.coeffs;
        ntt(&mut a);
        a
    }

    /// Create polynomial from NTT representation.
    pub fn from_ntt(mut ntt_repr: [u32; N]) -> Self {
        inv_ntt(&mut ntt_repr);
        Poly { coeffs: ntt_repr }
    }

    /// Compute the squared L2 norm, treating coefficients as centered in (-q/2, q/2].
    pub fn sq_norm(&self) -> u64 {
        let half_q = Q / 2;
        let mut norm: u64 = 0;
        for &c in &self.coeffs {
            let centered = if c > half_q {
                (Q - c) as i64
            } else {
                c as i64
            };
            norm += (centered * centered) as u64;
        }
        norm
    }

    /// Scalar multiply: multiply every coefficient by a scalar mod q.
    pub fn scalar_mul(&self, s: u32) -> Poly {
        let mut r = [0u32; N];
        let s64 = s as u64;
        for i in 0..N {
            r[i] = ((self.coeffs[i] as u64 * s64) % Q as u64) as u32;
        }
        Poly { coeffs: r }
    }

    /// Divide this polynomial by another in Z_q[X]/(X^N+1).
    /// Returns None if divisor is not invertible.
    pub fn div(&self, other: &Poly) -> Option<Poly> {
        let mut a = self.coeffs;
        let mut b = other.coeffs;
        ntt(&mut a);
        ntt(&mut b);
        let b_inv = ntt_inv(&b)?;
        let mut c = ntt_mul(&a, &b_inv);
        inv_ntt(&mut c);
        Some(Poly { coeffs: c })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_mul_identity() {
        let a = Poly::from_coeffs(&[1, 2, 3]);
        let one = Poly::from_coeffs(&[1]);
        let result = a.mul(&one);
        assert_eq!(result.coeffs[0], 1);
        assert_eq!(result.coeffs[1], 2);
        assert_eq!(result.coeffs[2], 3);
    }

    #[test]
    fn test_poly_add_sub() {
        let a = Poly::from_coeffs(&[100, 200]);
        let b = Poly::from_coeffs(&[50, 300]);
        let sum = a.add(&b);
        assert_eq!(sum.coeffs[0], 150);
        assert_eq!(sum.coeffs[1], 500);
        let diff = a.sub(&b);
        assert_eq!(diff.coeffs[0], 50);
        assert_eq!(diff.coeffs[1], (Q + 200 - 300) % Q);
    }

    #[test]
    fn test_poly_inv() {
        // Inverse of a constant polynomial
        let a = Poly::from_coeffs(&[3]);
        let a_inv = a.inv().unwrap();
        let product = a.mul(&a_inv);
        assert_eq!(product.coeffs[0], 1);
        for i in 1..N {
            assert_eq!(product.coeffs[i], 0, "nonzero at {}", i);
        }
    }

    #[test]
    fn test_poly_div() {
        let a = Poly::from_coeffs(&[6]);
        let b = Poly::from_coeffs(&[3]);
        let c = a.div(&b).unwrap();
        assert_eq!(c.coeffs[0], 2);
    }
}
