#![allow(dead_code)]

// Dilithium-2 (ML-DSA-44) parameters

/// Modulus q
pub const Q: i32 = 8380417;

/// Degree of polynomial ring
pub const N: usize = 256;

/// Matrix dimensions (k x l)
pub const K: usize = 4;
pub const L: usize = 4;

/// Secret key coefficient range
pub const ETA: usize = 2;

/// Number of +/-1 coefficients in challenge polynomial
pub const TAU: usize = 39;

/// Norm bound beta = tau * eta
pub const BETA: i32 = 78;

/// Masking range gamma1 = 2^17
pub const GAMMA1: i32 = 1 << 17;

/// Low-order rounding range gamma2 = (q-1)/88
pub const GAMMA2: i32 = 95232;

/// Dropped bits from t
pub const D: usize = 13;

/// Maximum number of 1s in hint
pub const OMEGA: usize = 80;

// Sizes for packing

/// Seed length in bytes
pub const SEEDBYTES: usize = 32;

/// CRH output length
pub const CRHBYTES: usize = 64;

/// Challenge seed length
pub const CTILDEBYTES: usize = 32;

/// Number of bits for packing t1 (10 bits per coefficient)
pub const POLYT1_PACKEDBYTES: usize = 320; // 256 * 10 / 8

/// Number of bytes for packing t0 (13 bits per coefficient)
pub const POLYT0_PACKEDBYTES: usize = 416; // 256 * 13 / 8

/// Number of bytes for packing eta=2 polynomial (3 bits per coefficient)
pub const POLYETA_PACKEDBYTES: usize = 96; // 256 * 3 / 8

/// Number of bytes for packing gamma1=2^17 polynomial (18 bits per coefficient)
pub const POLYZ_PACKEDBYTES: usize = 576; // 256 * 18 / 8

/// Number of bytes for packing w1 (6 bits per coefficient for gamma2=(q-1)/88)
pub const POLYW1_PACKEDBYTES: usize = 192; // 256 * 6 / 8

/// Public key size: seedbytes + k * polyt1_packed
pub const PK_BYTES: usize = SEEDBYTES + K * POLYT1_PACKEDBYTES;

/// Secret key size
pub const SK_BYTES: usize = 3 * SEEDBYTES + L * POLYETA_PACKEDBYTES + K * POLYETA_PACKEDBYTES + K * POLYT0_PACKEDBYTES;

/// Signature size: ctildebytes + l * polyz_packed + (omega + k)
pub const SIG_BYTES: usize = CTILDEBYTES + L * POLYZ_PACKEDBYTES + OMEGA + K;

// NTT constants

/// Montgomery parameter: 2^32 mod q
pub const MONT: i32 = -4186625; // 2^32 mod q (represented as negative for convenience)

/// q^{-1} mod 2^32
pub const QINV: i32 = 58728449; // q^(-1) mod 2^32

/// Primitive 256th root of unity mod q
pub const ROOT_OF_UNITY: i32 = 1753;
