#![allow(dead_code)]

use crate::dilithium::params::{N, Q, QINV};

/// Precomputed zetas (powers of root of unity in Montgomery form)
/// zetas[i] = mont * root^{brv(i)} mod q where mont = 2^32 mod q, root = 1753
/// zetas[0] = 0 (unused, slot 0 is never accessed in the butterfly)
const ZETAS: [i32; N] = [
           0,    25847, -2608894,  -518909,   237124,  -777960,  -876248,   466468,
     1826347,  2353451,  -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
     2725464,  1024112, -1079900,  3585928,  -549488, -1119584,  2619752, -2108549,
    -2118186, -3859737, -1399561, -3277672,  1757237,   -19422,  4010497,   280005,
     2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
    -3861115, -3043716,  3574422, -2867647,  3539968,  -300467,  2348700,  -539299,
    -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
      811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
    -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
    -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
     3412210,  -983419,  2147896,  2715295, -2967645, -3693493,  -411027, -2477047,
     -671102, -1228525,   -22981, -1308169,  -381987,  1349076,  1852771, -1430430,
    -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
    -3724342,    -8578,  1653064, -3249728,  2389356,  -210977,   759969, -1316856,
      189548, -3553272,  3159746, -1851402, -2409325,  -177440,  1315589,  1341330,
     1285669, -1584928,  -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
     2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
      266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
      900702,  1859098,   909542,   819034,   495491, -1613174,   -43260,  -522500,
     -655327, -3122442,  2031748,  3207046, -3556995,  -525098,  -768622, -3595838,
      342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
     2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
    -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
    -1333058,  1237275, -3318210, -1430225,  -451100,  1312455,  3306115, -1962642,
    -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
     -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
    -2013608,  2432395,  2454455,  -164721,  1957272,  3369112,   185531, -1207385,
    -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
    -3038916,  3523897,  3866901,   269760,  2213111,  -975884,  1717735,   472078,
     -426683,  1723600, -1803090,  1910376, -1667432, -1104333,  -260646, -3833893,
    -2939036, -2235985,  -420899, -2286327,   183443,  -976891,  1612842, -3545687,
     -554416,  3919660,   -48306, -1362209,  3937738,  1400424,  -846154,  1976782,
];

/// Montgomery reduction: given a (i64), compute a * 2^{-32} mod q
/// Returns value in (-q, q)
#[inline(always)]
pub fn montgomery_reduce(a: i64) -> i32 {
    let t = (a as i32).wrapping_mul(QINV);
    ((a - (t as i64) * (Q as i64)) >> 32) as i32
}

/// Reduce a mod q into range roughly (-q/2, q/2)
#[inline(always)]
pub fn reduce32(a: i32) -> i32 {
    let t = (a + (1 << 22)) >> 23;
    a - t * Q
}

/// Conditional addition of q: if a is negative, add q
#[inline(always)]
pub fn caddq(a: i32) -> i32 {
    a + ((a >> 31) & Q)
}

/// Freeze: fully reduce to [0, q)
#[inline(always)]
pub fn freeze(a: i32) -> i32 {
    caddq(reduce32(a))
}

/// Forward NTT: in-place NTT on a polynomial of degree < 256
pub fn ntt(a: &mut [i32; N]) {
    let mut k = 0usize;
    let mut len = 128;
    while len >= 1 {
        let mut start = 0;
        while start < N {
            k += 1;
            let zeta = ZETAS[k];
            for j in start..(start + len) {
                let t = montgomery_reduce(zeta as i64 * a[j + len] as i64);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
            start += 2 * len;
        }
        len >>= 1;
    }
}

/// Inverse NTT.
///
/// The scaling factor f = mont^2 / 256 mod q = 41978 is chosen so that
/// invntt(pointwise_montgomery(ntt(a), ntt(b))) = a * b mod (X^256+1).
/// Note: invntt(ntt(a)) != a; it returns a * mont instead. This is by design —
/// the NTT is only used together with pointwise_montgomery multiplication.
pub fn invntt(a: &mut [i32; N]) {
    let mut k = N;
    let mut len = 1;
    while len < N {
        let mut start = 0;
        while start < N {
            k -= 1;
            let zeta = -ZETAS[k];
            for j in start..(start + len) {
                let t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = montgomery_reduce(zeta as i64 * a[j + len] as i64);
            }
            start += 2 * len;
        }
        len <<= 1;
    }

    // f = mont^2 / 256 mod q = 41978
    let f: i64 = 41978;
    for coeff in a.iter_mut() {
        *coeff = montgomery_reduce(f * (*coeff as i64));
    }
}

/// Pointwise multiplication of two polynomials in NTT domain (Montgomery)
pub fn pointwise_montgomery(c: &mut [i32; N], a: &[i32; N], b: &[i32; N]) {
    for i in 0..N {
        c[i] = montgomery_reduce(a[i] as i64 * b[i] as i64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_invntt_roundtrip() {
        // invntt(ntt(a)) = a * mont (by design, since f = mont^2/n).
        // This is correct: the NTT is designed for use with pointwise_montgomery,
        // where the extra mont factor gets cancelled.
        let a_copy: [i32; N] = std::array::from_fn(|i| ((i as i32) * 17 + 3) % Q);
        let mut a = a_copy;

        ntt(&mut a);
        invntt(&mut a);

        // Result should be a_copy * mont mod q, where mont = 2^32 mod q = 4193792
        const MONT: i64 = 4193792;
        for i in 0..N {
            let val = freeze(a[i]);
            let expected = freeze(((a_copy[i] as i64 * MONT) % Q as i64) as i32);
            assert_eq!(val, expected, "mismatch at index {i}");
        }
    }

    #[test]
    fn test_ntt_multiplication() {
        // Test that ntt-based multiplication matches schoolbook multiplication
        let a: [i32; N] = std::array::from_fn(|i| ((i * 3 + 1) % 100) as i32);
        let b: [i32; N] = std::array::from_fn(|i| ((i * 7 + 2) % 100) as i32);

        // Schoolbook multiplication in Z_q[X]/(X^256+1)
        let mut c_ref = [0i64; N];
        for i in 0..N {
            for j in 0..N {
                if i + j < N {
                    c_ref[i + j] += a[i] as i64 * b[j] as i64;
                } else {
                    c_ref[i + j - N] -= a[i] as i64 * b[j] as i64;
                }
            }
        }

        // NTT multiplication
        let mut a_ntt = a;
        let mut b_ntt = b;
        ntt(&mut a_ntt);
        ntt(&mut b_ntt);
        let mut c_ntt = [0i32; N];
        pointwise_montgomery(&mut c_ntt, &a_ntt, &b_ntt);
        invntt(&mut c_ntt);

        for i in 0..N {
            let val = freeze(c_ntt[i]);
            let expected = freeze(((c_ref[i] % Q as i64) + Q as i64) as i32);
            assert_eq!(val, expected, "ntt mul mismatch at index {i}");
        }
    }
}
