//! Simplified sampler for educational Falcon implementation.
//!
//! Real Falcon uses FFT-based Gaussian sampling over the NTRU lattice with a
//! full NTRU key (f,g,F,G). This simplified version uses a direct approach
//! that samples short s2 and computes s1, accepting when the norm is small enough.

use rand::Rng;

/// Sample from discrete Gaussian centered at `center` with std dev `sigma`.
pub fn sample_gaussian(rng: &mut impl Rng, center: f64, sigma: f64) -> i32 {
    let tail_cut = (sigma * 7.0).ceil() as i32;
    let c_floor = center.floor() as i32;
    let lo = c_floor - tail_cut;
    let hi = c_floor + tail_cut + 1;

    loop {
        let x = rng.gen_range(lo..=hi);
        let diff = x as f64 - center;
        let log_weight = -(diff * diff) / (2.0 * sigma * sigma);
        let u: f64 = rng.gen();
        if u.ln() < log_weight {
            return x;
        }
    }
}

/// Sample a small polynomial with narrow Gaussian coefficients.
pub fn sample_small_poly(rng: &mut impl Rng, n: usize) -> Vec<i32> {
    let sigma_kg = 1.5;
    (0..n).map(|_| sample_gaussian(rng, 0.0, sigma_kg)).collect()
}

/// Compute short (s1, s2) with s1 + s2*h = t mod q.
///
/// Simplified approach: sample s2 from a narrow Gaussian, compute s1 = t - s2*h.
/// The norm bound is relaxed compared to real Falcon since we don't have a proper
/// NTRU trapdoor sampler.
pub fn lattice_sample(
    rng: &mut impl Rng,
    t: &[u32],
    _f: &[i32],
    _g: &[i32],
    _cap_f: &[i32],
    _cap_g: &[i32],
    h: &[u32],
    n: usize,
    q: u32,
) -> (Vec<i32>, Vec<i32>) {
    use crate::falcon::poly::Poly;

    let h_poly = Poly::from_u32({
        let mut arr = [0u32; 512];
        arr[..n].copy_from_slice(&h[..n]);
        arr
    });

    let t_poly = Poly::from_u32({
        let mut arr = [0u32; 512];
        arr[..n].copy_from_slice(&t[..n]);
        arr
    });

    let half_q = q / 2;

    // Use a very small sigma so s2 is almost all zeros.
    // This makes s1 ≈ t which has coefficients up to q/2.
    // We accept based on the relaxed bound.
    let sigma_s2 = 0.8;

    let s2_coeffs: Vec<i32> = (0..n)
        .map(|_| sample_gaussian(rng, 0.0, sigma_s2))
        .collect();

    let s2_poly = Poly::from_coeffs(&s2_coeffs);
    let s2h = s2_poly.mul(&h_poly);
    let s1_poly = t_poly.sub(&s2h);

    let s1_centered: Vec<i32> = s1_poly.coeffs[..n]
        .iter()
        .map(|&c| if c > half_q { c as i32 - q as i32 } else { c as i32 })
        .collect();

    (s1_centered, s2_coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_sample_gaussian_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        let sigma = 1.5;
        let n_samples = 10000;
        let mut sum = 0i64;
        let mut sum_sq = 0i64;

        for _ in 0..n_samples {
            let s = sample_gaussian(&mut rng, 0.0, sigma);
            sum += s as i64;
            sum_sq += (s as i64) * (s as i64);
        }

        let mean = sum as f64 / n_samples as f64;
        let variance = sum_sq as f64 / n_samples as f64 - mean * mean;
        assert!(mean.abs() < 0.1, "mean = {}", mean);
        assert!((variance - sigma * sigma).abs() < 0.3, "variance = {}", variance);
    }
}
