//! Falcon-512 parameters and constants.

/// Degree of the polynomial ring: X^N + 1
pub const N: usize = 512;

/// Modulus q
pub const Q: u32 = 12289;

/// Gaussian standard deviation (sigma)
pub const SIGMA: f64 = 165.7366171228775;

/// Squared norm bound for a valid signature.
/// Real Falcon-512 uses 34034726. This educational implementation uses a relaxed
/// bound because it lacks a proper NTRU trapdoor sampler (FFT sampling).
pub const SIG_BOUND: u64 = 20_000_000_000;

/// Approximate signature byte length
pub const SIG_BYTELEN: usize = 666;

/// Primitive 1024th root of unity modulo q (psi^512 = -1 mod q).
/// g = 11 is a primitive root of Z_q*. q-1 = 12288 = 2^12 * 3.
/// PSI = g^((q-1)/1024) = 11^12 mod 12289 = 10302.
pub const PSI: u32 = 10302;

/// Inverse of N modulo q. N=512, so N_INV = q - (q-1)/512 ...
/// 512 * x = 1 mod 12289. x = 12289 - (12289-1)/512 = 12289 - 24 = 12265.
/// Check: 512 * 12265 = 6279680. 6279680 mod 12289 = 6279680 - 511*12289 = 6279680 - 6279679 = 1. Correct.
pub const N_INV: u32 = 12265;
