#![allow(dead_code)]

use crate::dilithium::ntt;
use crate::dilithium::params::*;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::{Shake128, Shake256};

/// A polynomial in Z_q[X]/(X^256 + 1)
#[derive(Clone, Debug)]
pub struct Poly {
    pub coeffs: [i32; N],
}

impl Default for Poly {
    fn default() -> Self {
        Self { coeffs: [0i32; N] }
    }
}

impl Poly {
    pub fn new() -> Self {
        Self::default()
    }

    /// Reduce all coefficients mod q
    pub fn reduce(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = ntt::reduce32(*c);
        }
    }

    /// Conditional add q to negative coefficients
    pub fn caddq(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = ntt::caddq(*c);
        }
    }

    /// Freeze: fully reduce to [0, q)
    pub fn freeze(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = ntt::freeze(*c);
        }
    }

    /// Forward NTT in place
    pub fn ntt(&mut self) {
        ntt::ntt(&mut self.coeffs);
    }

    /// Inverse NTT in place
    pub fn invntt(&mut self) {
        ntt::invntt(&mut self.coeffs);
    }

    /// Add two polynomials: self = self + other
    pub fn add_assign(&mut self, other: &Poly) {
        for i in 0..N {
            self.coeffs[i] += other.coeffs[i];
        }
    }

    /// Subtract: self = self - other
    pub fn sub_assign(&mut self, other: &Poly) {
        for i in 0..N {
            self.coeffs[i] -= other.coeffs[i];
        }
    }

    /// Shift left by d bits (multiply by 2^d)
    pub fn shiftl(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c <<= D;
        }
    }

    /// Check infinity norm <= bound
    pub fn chknorm(&self, bound: i32) -> bool {
        for &c in &self.coeffs {
            // Reduce to centered representation
            let mut t = ntt::reduce32(c);
            // Absolute value
            t = t - ((t >> 31) & (2 * t));
            if t >= bound {
                return false;
            }
        }
        true
    }

    /// Pointwise Montgomery multiplication: self = a * b (both in NTT domain)
    pub fn pointwise_montgomery(a: &Poly, b: &Poly) -> Poly {
        let mut c = Poly::new();
        ntt::pointwise_montgomery(&mut c.coeffs, &a.coeffs, &b.coeffs);
        c
    }

    // === Sampling ===

    /// Sample polynomial with uniformly random coefficients in [0, q) from SHAKE128 stream
    /// Used for expanding the matrix A from seed rho
    pub fn uniform(seed: &[u8; SEEDBYTES], nonce_i: u8, nonce_j: u8) -> Poly {
        let mut poly = Poly::new();
        let mut hasher = Shake128::default();
        hasher.update(seed);
        hasher.update(&[nonce_j, nonce_i]); // Note: column-major order per spec
        let mut reader = hasher.finalize_xof();

        let mut ctr = 0usize;
        let mut buf = [0u8; 3];
        while ctr < N {
            reader.read(&mut buf);
            // Rejection sampling from 3-byte (23-bit) chunks
            let val = ((buf[0] as u32)
                | ((buf[1] as u32) << 8)
                | ((buf[2] as u32) << 16))
                & 0x7FFFFF;
            if val < Q as u32 {
                poly.coeffs[ctr] = val as i32;
                ctr += 1;
            }
        }
        poly
    }

    /// Sample short polynomial with coefficients in [-eta, eta]
    /// using SHAKE256 with seed || nonce
    pub fn sample_eta(seed: &[u8], nonce: u16) -> Poly {
        let mut poly = Poly::new();
        let mut hasher = Shake256::default();
        hasher.update(seed);
        hasher.update(&nonce.to_le_bytes());
        let mut reader = hasher.finalize_xof();

        // For eta=2, sample from {0,1,2,3,4} and map to {-2,-1,0,1,2}
        // Use rejection from 4 bits: values 0..14 are valid (mod 5)
        // Actually the spec uses CBD (centered binomial) for eta=2:
        // sample 2*eta bits, count ones in first eta bits minus count of ones in second eta bits
        // But Dilithium uses uniform sampling with rejection for eta

        // Per FIPS 204 / Dilithium: for eta=2, we need uniform sampling
        // Each coefficient: sample a nibble (4 bits), if < 15, coeff = nibble mod 5 - 2...
        // Actually the reference uses: sample byte, split into two nibbles, each < 15 maps to eta range
        // Let me follow the reference implementation exactly:

        let buflen = N / 2; // 128 bytes for eta=2 (2 coefficients per byte from nibbles)
        let mut buf = vec![0u8; buflen];
        reader.read(&mut buf);

        let mut ctr = 0;
        for &b in &buf {
            if ctr >= N {
                break;
            }
            // Low nibble
            let t0 = (b & 0x0F) as i32;
            // High nibble
            let t1 = (b >> 4) as i32;

            // For eta=2: each nibble encodes coefficient via:
            // t mod 5 gives value in {0,1,2,3,4}, then subtract eta to get {-2,-1,0,1,2}
            // But we need rejection: only accept if nibble < 15
            // Actually for eta=2, per reference: coeff = (t0 & 0x07) - ((t0 >> 2) & 0x03) style...
            // Let me use the exact reference approach for eta=2:
            // Sample uniform byte, split nibbles, each nibble -> coefficient via:
            // bits b0,b1,b2,b3 -> (b0+b1) - (b2+b3) which gives range [-2, 2]
            if ctr < N {
                let b0 = t0 & 1;
                let b1 = (t0 >> 1) & 1;
                let b2 = (t0 >> 2) & 1;
                let b3 = (t0 >> 3) & 1;
                poly.coeffs[ctr] = (b0 + b1) - (b2 + b3);
                ctr += 1;
            }
            if ctr < N {
                let b0 = t1 & 1;
                let b1 = (t1 >> 1) & 1;
                let b2 = (t1 >> 2) & 1;
                let b3 = (t1 >> 3) & 1;
                poly.coeffs[ctr] = (b0 + b1) - (b2 + b3);
                ctr += 1;
            }
        }

        // If we still need more (shouldn't happen for eta=2 with 128 bytes -> 256 coeffs)
        while ctr < N {
            let mut extra = [0u8; 1];
            reader.read(&mut extra);
            let t0 = (extra[0] & 0x0F) as i32;
            let t1 = (extra[0] >> 4) as i32;
            if ctr < N {
                poly.coeffs[ctr] = (t0 & 1) + ((t0 >> 1) & 1) - ((t0 >> 2) & 1) - ((t0 >> 3) & 1);
                ctr += 1;
            }
            if ctr < N {
                poly.coeffs[ctr] = (t1 & 1) + ((t1 >> 1) & 1) - ((t1 >> 2) & 1) - ((t1 >> 3) & 1);
                ctr += 1;
            }
        }

        poly
    }

    /// Sample mask polynomial with coefficients in [-gamma1+1, gamma1]
    /// using SHAKE256 with seed || nonce
    pub fn sample_gamma1(seed: &[u8], nonce: u16) -> Poly {
        let mut poly = Poly::new();
        let mut hasher = Shake256::default();
        hasher.update(seed);
        hasher.update(&nonce.to_le_bytes());
        let mut reader = hasher.finalize_xof();

        // For gamma1 = 2^17, we need 18 bits per coefficient
        // Pack: 9 bytes -> 4 coefficients (9*8 = 72 = 4*18)
        let buflen = N * 18 / 8; // 576 bytes
        let mut buf = vec![0u8; buflen];
        reader.read(&mut buf);

        for i in 0..(N / 4) {
            // Read 9 bytes = 72 bits = 4 * 18 bits
            let base = i * 9;
            let b = &buf[base..base + 9];

            let mut vals = [0u32; 4];
            vals[0] = (b[0] as u32)
                | ((b[1] as u32) << 8)
                | (((b[2] as u32) & 0x03) << 16);
            vals[1] = ((b[2] as u32) >> 2)
                | ((b[3] as u32) << 6)
                | (((b[4] as u32) & 0x0F) << 14);
            vals[2] = ((b[4] as u32) >> 4)
                | ((b[5] as u32) << 4)
                | (((b[6] as u32) & 0x3F) << 12);
            vals[3] = ((b[6] as u32) >> 6)
                | ((b[7] as u32) << 2)
                | ((b[8] as u32) << 10);

            for j in 0..4 {
                let v = vals[j] & 0x3FFFF; // 18 bits
                // Map: gamma1 - v gives coefficient in centered range
                poly.coeffs[4 * i + j] = GAMMA1 - (v as i32);
            }
        }

        poly
    }

    /// Sample challenge polynomial c with exactly tau coefficients in {-1, +1}
    /// and the rest zero, from a seed
    pub fn sample_challenge(seed: &[u8; CTILDEBYTES]) -> Poly {
        let mut poly = Poly::new();
        let mut hasher = Shake256::default();
        hasher.update(seed);
        let mut reader = hasher.finalize_xof();

        // Read 8 bytes for sign bits
        let mut signs_bytes = [0u8; 8];
        reader.read(&mut signs_bytes);
        let signs = u64::from_le_bytes(signs_bytes);

        let mut idx_buf = [0u8; 1];
        for i in (N - TAU)..N {
            // Rejection sample an index j in [0, i]
            loop {
                reader.read(&mut idx_buf);
                let j = idx_buf[0] as usize;
                if j <= i {
                    poly.coeffs[i] = poly.coeffs[j];
                    let sign_bit = (signs >> (i - (N - TAU))) & 1;
                    poly.coeffs[j] = if sign_bit == 0 { 1 } else { -1 };
                    break;
                }
            }
        }

        poly
    }

    // === Decomposition ===

    /// Power2Round: decompose a into (a1, a0) such that a = a1*2^d + a0
    /// with -2^{d-1} < a0 <= 2^{d-1}
    pub fn power2round(a: i32) -> (i32, i32) {
        let a = ntt::freeze(a);
        let a0 = a - ((a + (1 << (D - 1)) - 1) >> D << D);
        let a1 = (a - a0) >> D;
        (a1, a0)
    }

    /// Decompose: decompose a into (a1, a0) such that a = a1*2*gamma2 + a0
    /// with -gamma2 < a0 <= gamma2
    pub fn decompose(a: i32) -> (i32, i32) {
        let a = ntt::freeze(a);
        let mut a0 = a % (2 * GAMMA2);
        if a0 > GAMMA2 {
            a0 -= 2 * GAMMA2;
        }
        let a1 = if a - a0 == Q - 1 {
            // Special case
            0
        } else {
            (a - a0) / (2 * GAMMA2)
        };
        let a0 = if a - a0 == Q - 1 { a0 - 1 } else { a0 };
        (a1, a0)
    }

    /// HighBits: extract high bits
    pub fn highbits(a: i32) -> i32 {
        Self::decompose(a).0
    }

    /// LowBits: extract low bits
    pub fn lowbits(a: i32) -> i32 {
        Self::decompose(a).1
    }

    /// MakeHint: returns 1 if highbits would change when adding a0 to a1*alpha
    pub fn make_hint(a0: i32, a1: i32) -> i32 {
        if a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0) {
            1
        } else {
            0
        }
    }

    /// UseHint: use hint to recover high bits
    pub fn use_hint(a: i32, hint: i32) -> i32 {
        let (a1, a0) = Self::decompose(a);
        if hint == 0 {
            return a1;
        }
        // m = (q-1) / (2*gamma2)
        let m = (Q - 1) / (2 * GAMMA2);
        if a0 > 0 {
            (a1 + 1) % m
        } else {
            (a1 - 1 + m) % m
        }
    }

    // === Packing ===

    /// Pack t1 polynomial (10 bits per coefficient)
    pub fn pack_t1(&self) -> Vec<u8> {
        let mut buf = vec![0u8; POLYT1_PACKEDBYTES];
        for i in 0..(N / 4) {
            buf[5 * i + 0] = (self.coeffs[4 * i + 0] & 0x3FF) as u8;
            buf[5 * i + 1] = ((self.coeffs[4 * i + 0] >> 8) | (self.coeffs[4 * i + 1] << 2)) as u8;
            buf[5 * i + 2] = ((self.coeffs[4 * i + 1] >> 6) | (self.coeffs[4 * i + 2] << 4)) as u8;
            buf[5 * i + 3] = ((self.coeffs[4 * i + 2] >> 4) | (self.coeffs[4 * i + 3] << 6)) as u8;
            buf[5 * i + 4] = (self.coeffs[4 * i + 3] >> 2) as u8;
        }
        buf
    }

    /// Unpack t1 polynomial
    pub fn unpack_t1(buf: &[u8]) -> Poly {
        let mut poly = Poly::new();
        for i in 0..(N / 4) {
            poly.coeffs[4 * i + 0] = ((buf[5 * i + 0] as i32) | ((buf[5 * i + 1] as i32) << 8)) & 0x3FF;
            poly.coeffs[4 * i + 1] = ((buf[5 * i + 1] as i32 >> 2) | ((buf[5 * i + 2] as i32) << 6)) & 0x3FF;
            poly.coeffs[4 * i + 2] = ((buf[5 * i + 2] as i32 >> 4) | ((buf[5 * i + 3] as i32) << 4)) & 0x3FF;
            poly.coeffs[4 * i + 3] = ((buf[5 * i + 3] as i32 >> 6) | ((buf[5 * i + 4] as i32) << 2)) & 0x3FF;
        }
        poly
    }

    /// Pack t0 polynomial (13 bits per coefficient, centered around 2^{d-1})
    pub fn pack_t0(&self) -> Vec<u8> {
        let mut buf = vec![0u8; POLYT0_PACKEDBYTES];
        for i in 0..(N / 8) {
            let mut t = [0i32; 8];
            for j in 0..8 {
                t[j] = (1 << (D - 1)) - self.coeffs[8 * i + j];
            }
            buf[13 * i + 0] = t[0] as u8;
            buf[13 * i + 1] = (t[0] >> 8) as u8;
            buf[13 * i + 1] |= (t[1] << 5) as u8;
            buf[13 * i + 2] = (t[1] >> 3) as u8;
            buf[13 * i + 3] = (t[1] >> 11) as u8;
            buf[13 * i + 3] |= (t[2] << 2) as u8;
            buf[13 * i + 4] = (t[2] >> 6) as u8;
            buf[13 * i + 4] |= (t[3] << 7) as u8;
            buf[13 * i + 5] = (t[3] >> 1) as u8;
            buf[13 * i + 6] = (t[3] >> 9) as u8;
            buf[13 * i + 6] |= (t[4] << 4) as u8;
            buf[13 * i + 7] = (t[4] >> 4) as u8;
            buf[13 * i + 8] = (t[4] >> 12) as u8;
            buf[13 * i + 8] |= (t[5] << 1) as u8;
            buf[13 * i + 9] = (t[5] >> 7) as u8;
            buf[13 * i + 9] |= (t[6] << 6) as u8;
            buf[13 * i + 10] = (t[6] >> 2) as u8;
            buf[13 * i + 11] = (t[6] >> 10) as u8;
            buf[13 * i + 11] |= (t[7] << 3) as u8;
            buf[13 * i + 12] = (t[7] >> 5) as u8;
        }
        buf
    }

    /// Unpack t0 polynomial
    pub fn unpack_t0(buf: &[u8]) -> Poly {
        let mut poly = Poly::new();
        for i in 0..(N / 8) {
            poly.coeffs[8 * i + 0] = (buf[13 * i + 0] as i32)
                | ((buf[13 * i + 1] as i32 & 0x1F) << 8);
            poly.coeffs[8 * i + 1] = ((buf[13 * i + 1] as i32) >> 5)
                | ((buf[13 * i + 2] as i32) << 3)
                | ((buf[13 * i + 3] as i32 & 0x03) << 11);
            poly.coeffs[8 * i + 2] = ((buf[13 * i + 3] as i32) >> 2)
                | ((buf[13 * i + 4] as i32 & 0x7F) << 6);
            poly.coeffs[8 * i + 3] = ((buf[13 * i + 4] as i32) >> 7)
                | ((buf[13 * i + 5] as i32) << 1)
                | ((buf[13 * i + 6] as i32 & 0x0F) << 9);
            poly.coeffs[8 * i + 4] = ((buf[13 * i + 6] as i32) >> 4)
                | ((buf[13 * i + 7] as i32) << 4)
                | ((buf[13 * i + 8] as i32 & 0x01) << 12);
            poly.coeffs[8 * i + 5] = ((buf[13 * i + 8] as i32) >> 1)
                | ((buf[13 * i + 9] as i32 & 0x3F) << 7);
            poly.coeffs[8 * i + 6] = ((buf[13 * i + 9] as i32) >> 6)
                | ((buf[13 * i + 10] as i32) << 2)
                | ((buf[13 * i + 11] as i32 & 0x07) << 10);
            poly.coeffs[8 * i + 7] = ((buf[13 * i + 11] as i32) >> 3)
                | ((buf[13 * i + 12] as i32) << 5);

            for j in 0..8 {
                poly.coeffs[8 * i + j] = (1 << (D - 1)) - poly.coeffs[8 * i + j];
            }
        }
        poly
    }

    /// Pack polynomial with eta-range coefficients (eta=2, 3 bits per coeff)
    pub fn pack_eta(&self) -> Vec<u8> {
        let mut buf = vec![0u8; POLYETA_PACKEDBYTES];
        // eta=2: 8 coefficients per 3 bytes
        for i in 0..(N / 8) {
            let mut t = [0u8; 8];
            for j in 0..8 {
                t[j] = (ETA as i32 - self.coeffs[8 * i + j]) as u8;
            }
            buf[3 * i + 0] = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
            buf[3 * i + 1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
            buf[3 * i + 2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
        }
        buf
    }

    /// Unpack polynomial with eta-range coefficients
    pub fn unpack_eta(buf: &[u8]) -> Poly {
        let mut poly = Poly::new();
        for i in 0..(N / 8) {
            poly.coeffs[8 * i + 0] = ((buf[3 * i + 0] >> 0) & 7) as i32;
            poly.coeffs[8 * i + 1] = ((buf[3 * i + 0] >> 3) & 7) as i32;
            poly.coeffs[8 * i + 2] = (((buf[3 * i + 0] >> 6) | (buf[3 * i + 1] << 2)) & 7) as i32;
            poly.coeffs[8 * i + 3] = ((buf[3 * i + 1] >> 1) & 7) as i32;
            poly.coeffs[8 * i + 4] = ((buf[3 * i + 1] >> 4) & 7) as i32;
            poly.coeffs[8 * i + 5] = (((buf[3 * i + 1] >> 7) | (buf[3 * i + 2] << 1)) & 7) as i32;
            poly.coeffs[8 * i + 6] = ((buf[3 * i + 2] >> 2) & 7) as i32;
            poly.coeffs[8 * i + 7] = ((buf[3 * i + 2] >> 5) & 7) as i32;

            for j in 0..8 {
                poly.coeffs[8 * i + j] = ETA as i32 - poly.coeffs[8 * i + j];
            }
        }
        poly
    }

    /// Pack z polynomial (gamma1 = 2^17, 18 bits per coeff)
    pub fn pack_z(&self) -> Vec<u8> {
        let mut buf = vec![0u8; POLYZ_PACKEDBYTES];
        for i in 0..(N / 4) {
            let mut t = [0u32; 4];
            for j in 0..4 {
                t[j] = (GAMMA1 - self.coeffs[4 * i + j]) as u32;
            }
            buf[9 * i + 0] = t[0] as u8;
            buf[9 * i + 1] = (t[0] >> 8) as u8;
            buf[9 * i + 2] = (t[0] >> 16) as u8 | (t[1] << 2) as u8;
            buf[9 * i + 3] = (t[1] >> 6) as u8;
            buf[9 * i + 4] = (t[1] >> 14) as u8 | (t[2] << 4) as u8;
            buf[9 * i + 5] = (t[2] >> 4) as u8;
            buf[9 * i + 6] = (t[2] >> 12) as u8 | (t[3] << 6) as u8;
            buf[9 * i + 7] = (t[3] >> 2) as u8;
            buf[9 * i + 8] = (t[3] >> 10) as u8;
        }
        buf
    }

    /// Unpack z polynomial
    pub fn unpack_z(buf: &[u8]) -> Poly {
        let mut poly = Poly::new();
        for i in 0..(N / 4) {
            poly.coeffs[4 * i + 0] = (buf[9 * i + 0] as i32)
                | ((buf[9 * i + 1] as i32) << 8)
                | (((buf[9 * i + 2] as i32) & 0x03) << 16);
            poly.coeffs[4 * i + 1] = ((buf[9 * i + 2] as i32) >> 2)
                | ((buf[9 * i + 3] as i32) << 6)
                | (((buf[9 * i + 4] as i32) & 0x0F) << 14);
            poly.coeffs[4 * i + 2] = ((buf[9 * i + 4] as i32) >> 4)
                | ((buf[9 * i + 5] as i32) << 4)
                | (((buf[9 * i + 6] as i32) & 0x3F) << 12);
            poly.coeffs[4 * i + 3] = ((buf[9 * i + 6] as i32) >> 6)
                | ((buf[9 * i + 7] as i32) << 2)
                | ((buf[9 * i + 8] as i32) << 10);

            for j in 0..4 {
                poly.coeffs[4 * i + j] = GAMMA1 - poly.coeffs[4 * i + j];
            }
        }
        poly
    }

    /// Pack w1 polynomial (6 bits per coefficient for gamma2 = (q-1)/88)
    pub fn pack_w1(&self) -> Vec<u8> {
        let mut buf = vec![0u8; POLYW1_PACKEDBYTES];
        for i in 0..(N / 4) {
            buf[3 * i + 0] = (self.coeffs[4 * i + 0]) as u8
                | ((self.coeffs[4 * i + 1]) << 6) as u8;
            buf[3 * i + 1] = (self.coeffs[4 * i + 1] >> 2) as u8
                | ((self.coeffs[4 * i + 2]) << 4) as u8;
            buf[3 * i + 2] = (self.coeffs[4 * i + 2] >> 4) as u8
                | ((self.coeffs[4 * i + 3]) << 2) as u8;
        }
        buf
    }
}

/// A vector of K polynomials
#[derive(Clone, Debug)]
pub struct PolyVecK {
    pub polys: [Poly; K],
}

impl Default for PolyVecK {
    fn default() -> Self {
        Self {
            polys: std::array::from_fn(|_| Poly::new()),
        }
    }
}

impl PolyVecK {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn ntt(&mut self) {
        for p in self.polys.iter_mut() {
            p.ntt();
        }
    }

    pub fn invntt(&mut self) {
        for p in self.polys.iter_mut() {
            p.invntt();
        }
    }

    pub fn reduce(&mut self) {
        for p in self.polys.iter_mut() {
            p.reduce();
        }
    }

    pub fn caddq(&mut self) {
        for p in self.polys.iter_mut() {
            p.caddq();
        }
    }

    pub fn freeze(&mut self) {
        for p in self.polys.iter_mut() {
            p.freeze();
        }
    }

    pub fn add_assign(&mut self, other: &PolyVecK) {
        for i in 0..K {
            self.polys[i].add_assign(&other.polys[i]);
        }
    }

    pub fn sub_assign(&mut self, other: &PolyVecK) {
        for i in 0..K {
            self.polys[i].sub_assign(&other.polys[i]);
        }
    }

    /// Check that all polynomials have infinity norm < bound
    pub fn chknorm(&self, bound: i32) -> bool {
        self.polys.iter().all(|p| p.chknorm(bound))
    }

    /// Power2Round on each coefficient: returns (v1, v0)
    pub fn power2round(&self) -> (PolyVecK, PolyVecK) {
        let mut v1 = PolyVecK::new();
        let mut v0 = PolyVecK::new();
        for i in 0..K {
            for j in 0..N {
                let (a1, a0) = Poly::power2round(self.polys[i].coeffs[j]);
                v1.polys[i].coeffs[j] = a1;
                v0.polys[i].coeffs[j] = a0;
            }
        }
        (v1, v0)
    }

    /// Decompose on each coefficient
    pub fn decompose(&self) -> (PolyVecK, PolyVecK) {
        let mut v1 = PolyVecK::new();
        let mut v0 = PolyVecK::new();
        for i in 0..K {
            for j in 0..N {
                let (a1, a0) = Poly::decompose(self.polys[i].coeffs[j]);
                v1.polys[i].coeffs[j] = a1;
                v0.polys[i].coeffs[j] = a0;
            }
        }
        (v1, v0)
    }

    /// Make hint vector
    pub fn make_hint(v0: &PolyVecK, v1: &PolyVecK) -> (PolyVecK, usize) {
        let mut h = PolyVecK::new();
        let mut count = 0;
        for i in 0..K {
            for j in 0..N {
                h.polys[i].coeffs[j] = Poly::make_hint(v0.polys[i].coeffs[j], v1.polys[i].coeffs[j]);
                count += h.polys[i].coeffs[j] as usize;
            }
        }
        (h, count)
    }

    /// Use hint vector
    pub fn use_hint(u: &PolyVecK, h: &PolyVecK) -> PolyVecK {
        let mut r = PolyVecK::new();
        for i in 0..K {
            for j in 0..N {
                r.polys[i].coeffs[j] = Poly::use_hint(u.polys[i].coeffs[j], h.polys[i].coeffs[j]);
            }
        }
        r
    }

    /// Pack w1 for hashing
    pub fn pack_w1(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(K * POLYW1_PACKEDBYTES);
        for i in 0..K {
            buf.extend_from_slice(&self.polys[i].pack_w1());
        }
        buf
    }
}

/// A vector of L polynomials
#[derive(Clone, Debug)]
pub struct PolyVecL {
    pub polys: [Poly; L],
}

impl Default for PolyVecL {
    fn default() -> Self {
        Self {
            polys: std::array::from_fn(|_| Poly::new()),
        }
    }
}

impl PolyVecL {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn ntt(&mut self) {
        for p in self.polys.iter_mut() {
            p.ntt();
        }
    }

    pub fn invntt(&mut self) {
        for p in self.polys.iter_mut() {
            p.invntt();
        }
    }

    pub fn reduce(&mut self) {
        for p in self.polys.iter_mut() {
            p.reduce();
        }
    }

    pub fn caddq(&mut self) {
        for p in self.polys.iter_mut() {
            p.caddq();
        }
    }

    pub fn freeze(&mut self) {
        for p in self.polys.iter_mut() {
            p.freeze();
        }
    }

    pub fn add_assign(&mut self, other: &PolyVecL) {
        for i in 0..L {
            self.polys[i].add_assign(&other.polys[i]);
        }
    }

    /// Check norm bound
    pub fn chknorm(&self, bound: i32) -> bool {
        self.polys.iter().all(|p| p.chknorm(bound))
    }

    /// Pointwise multiply accumulate: result[i] = sum_j(row[j] * v[j]) for each polynomial
    /// row is one row of the matrix (L polynomials), v is the vector
    /// All in NTT domain
    pub fn pointwise_acc(row: &[Poly; L], v: &PolyVecL) -> Poly {
        let mut acc = Poly::pointwise_montgomery(&row[0], &v.polys[0]);
        for j in 1..L {
            let t = Poly::pointwise_montgomery(&row[j], &v.polys[j]);
            acc.add_assign(&t);
        }
        acc.reduce();
        acc
    }
}

/// Matrix A represented as K rows of L polynomials (in NTT domain)
#[derive(Clone, Debug)]
pub struct MatrixA {
    pub rows: [[Poly; L]; K],
}

impl MatrixA {
    /// Expand matrix A from seed rho using SHAKE128
    pub fn expand(rho: &[u8; SEEDBYTES]) -> Self {
        let rows: [[Poly; L]; K] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let mut p = Poly::uniform(rho, i as u8, j as u8);
                p.ntt();
                p
            })
        });
        MatrixA { rows }
    }

    /// Multiply matrix A (in NTT domain) by vector v (in NTT domain): result = A * v
    pub fn mul_ntt(&self, v: &PolyVecL) -> PolyVecK {
        let mut result = PolyVecK::new();
        for i in 0..K {
            result.polys[i] = PolyVecL::pointwise_acc(&self.rows[i], v);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pack_unpack_t1() {
        let mut p = Poly::new();
        for i in 0..N {
            p.coeffs[i] = (i as i32 * 3 + 7) % 1024; // t1 is 10 bits
        }
        let packed = p.pack_t1();
        let unpacked = Poly::unpack_t1(&packed);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i], "t1 mismatch at {i}");
        }
    }

    #[test]
    fn test_pack_unpack_eta() {
        let p = Poly::sample_eta(&[0u8; 64], 0);
        for c in &p.coeffs {
            assert!(*c >= -(ETA as i32) && *c <= ETA as i32);
        }
        let packed = p.pack_eta();
        let unpacked = Poly::unpack_eta(&packed);
        for i in 0..N {
            assert_eq!(p.coeffs[i], unpacked.coeffs[i], "eta mismatch at {i}");
        }
    }

    #[test]
    fn test_sample_challenge() {
        let seed = [42u8; CTILDEBYTES];
        let c = Poly::sample_challenge(&seed);
        let nonzero: usize = c.coeffs.iter().filter(|&&x| x != 0).count();
        assert_eq!(nonzero, TAU);
        for &coeff in &c.coeffs {
            assert!(coeff == -1 || coeff == 0 || coeff == 1);
        }
    }

    #[test]
    fn test_decompose() {
        for a in [0, 1, Q / 2, Q - 1, 95232, 95233, 190463, 190464] {
            let (a1, a0) = Poly::decompose(a);
            let recomposed = ntt::freeze(a1 * 2 * GAMMA2 + a0);
            let orig = ntt::freeze(a);
            assert_eq!(recomposed, orig, "decompose failed for a={a}: a1={a1}, a0={a0}");
        }
    }
}
