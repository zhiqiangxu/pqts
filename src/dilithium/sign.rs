#![allow(dead_code)]

use crate::dilithium::params::*;
use crate::dilithium::poly::*;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

/// Dilithium-2 public key
#[derive(Clone, Debug)]
pub struct PublicKey {
    /// Seed for regenerating matrix A
    pub rho: [u8; SEEDBYTES],
    /// Packed t1 vector
    pub t1: PolyVecK,
}

/// Dilithium-2 secret key
#[derive(Clone, Debug)]
pub struct SecretKey {
    /// Seed for matrix A
    pub rho: [u8; SEEDBYTES],
    /// Key for signing
    pub key: [u8; SEEDBYTES],
    /// Hash of public key: H(rho || pack(t1))
    pub tr: [u8; SEEDBYTES],
    /// Secret vectors
    pub s1: PolyVecL,
    pub s2: PolyVecK,
    /// Low bits of t
    pub t0: PolyVecK,
}

/// Dilithium-2 signature
#[derive(Clone, Debug)]
pub struct Signature {
    /// Challenge hash seed
    pub c_tilde: [u8; CTILDEBYTES],
    /// Response vector z
    pub z: PolyVecL,
    /// Hint vector h
    pub h: PolyVecK,
}

impl PublicKey {
    /// Serialize public key
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(PK_BYTES);
        buf.extend_from_slice(&self.rho);
        for i in 0..K {
            buf.extend_from_slice(&self.t1.polys[i].pack_t1());
        }
        buf
    }

    /// Deserialize public key
    pub fn from_bytes(buf: &[u8]) -> Option<Self> {
        if buf.len() != PK_BYTES {
            return None;
        }
        let mut rho = [0u8; SEEDBYTES];
        rho.copy_from_slice(&buf[..SEEDBYTES]);
        let mut t1 = PolyVecK::new();
        for i in 0..K {
            let start = SEEDBYTES + i * POLYT1_PACKEDBYTES;
            t1.polys[i] = Poly::unpack_t1(&buf[start..start + POLYT1_PACKEDBYTES]);
        }
        Some(PublicKey { rho, t1 })
    }
}

impl Signature {
    /// Serialize signature
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(SIG_BYTES);
        buf.extend_from_slice(&self.c_tilde);
        for i in 0..L {
            buf.extend_from_slice(&self.z.polys[i].pack_z());
        }
        // Encode hints
        self.encode_hints(&mut buf);
        buf
    }

    /// Encode hint vector h into the signature
    fn encode_hints(&self, buf: &mut Vec<u8>) {
        // Hints are encoded as: for each polynomial, list the positions of 1s,
        // then store the end index for each polynomial.
        let mut hints_buf = vec![0u8; OMEGA + K];
        let mut k_pos = 0usize;
        for i in 0..K {
            for j in 0..N {
                if self.h.polys[i].coeffs[j] != 0 {
                    if k_pos < OMEGA {
                        hints_buf[k_pos] = j as u8;
                        k_pos += 1;
                    }
                }
            }
            hints_buf[OMEGA + i] = k_pos as u8;
        }
        buf.extend_from_slice(&hints_buf);
    }

    /// Deserialize signature
    pub fn from_bytes(buf: &[u8]) -> Option<Self> {
        if buf.len() != SIG_BYTES {
            return None;
        }

        let mut c_tilde = [0u8; CTILDEBYTES];
        c_tilde.copy_from_slice(&buf[..CTILDEBYTES]);

        let mut z = PolyVecL::new();
        for i in 0..L {
            let start = CTILDEBYTES + i * POLYZ_PACKEDBYTES;
            z.polys[i] = Poly::unpack_z(&buf[start..start + POLYZ_PACKEDBYTES]);
        }

        // Decode hints
        let hints_start = CTILDEBYTES + L * POLYZ_PACKEDBYTES;
        let hints_buf = &buf[hints_start..];
        let mut h = PolyVecK::new();
        let mut k_pos = 0usize;
        for i in 0..K {
            let end = hints_buf[OMEGA + i] as usize;
            if end < k_pos || end > OMEGA {
                return None;
            }
            for idx in k_pos..end {
                let j = hints_buf[idx] as usize;
                if j >= N {
                    return None;
                }
                // Check that indices are in ascending order
                if idx > k_pos && hints_buf[idx] <= hints_buf[idx - 1] {
                    return None;
                }
                h.polys[i].coeffs[j] = 1;
            }
            k_pos = end;
        }
        // Remaining positions should be zero
        for idx in k_pos..OMEGA {
            if hints_buf[idx] != 0 {
                return None;
            }
        }

        Some(Signature { c_tilde, z, h })
    }
}

/// Compute CRH (collision-resistant hash) using SHAKE256
fn crh(input: &[u8], output_len: usize) -> Vec<u8> {
    let mut hasher = Shake256::default();
    hasher.update(input);
    let mut reader = hasher.finalize_xof();
    let mut out = vec![0u8; output_len];
    reader.read(&mut out);
    out
}

/// Generate Dilithium-2 keypair
pub fn keygen() -> (PublicKey, SecretKey) {
    use rand::RngCore;
    let mut rng = rand::thread_rng();

    // 1. Generate random seed
    let mut zeta = [0u8; SEEDBYTES];
    rng.fill_bytes(&mut zeta);

    // 2. Expand seed: (rho, rho', K) = H(zeta)
    let mut hasher = Shake256::default();
    hasher.update(&zeta);
    let mut reader = hasher.finalize_xof();
    let mut rho = [0u8; SEEDBYTES];
    let mut rho_prime = [0u8; CRHBYTES];
    let mut key = [0u8; SEEDBYTES];
    reader.read(&mut rho);
    reader.read(&mut rho_prime);
    reader.read(&mut key);

    // 3. Expand matrix A from rho
    let a_hat = MatrixA::expand(&rho);

    // 4. Sample secret vectors s1, s2
    let mut s1 = PolyVecL::new();
    let mut s2 = PolyVecK::new();
    for i in 0..L {
        s1.polys[i] = Poly::sample_eta(&rho_prime, i as u16);
    }
    for i in 0..K {
        s2.polys[i] = Poly::sample_eta(&rho_prime, (L + i) as u16);
    }

    // 5. Compute t = A * NTT(s1) + s2
    let mut s1_hat = s1.clone();
    s1_hat.ntt();
    let mut t = a_hat.mul_ntt(&s1_hat);
    t.invntt();
    t.add_assign(&s2);
    t.caddq();

    // 6. Power2Round: t = t1 * 2^d + t0
    let (t1, t0) = t.power2round();

    // 7. Compute tr = H(pk)
    let pk = PublicKey {
        rho,
        t1: t1.clone(),
    };
    let pk_bytes = pk.to_bytes();
    let tr_vec = crh(&pk_bytes, SEEDBYTES);
    let mut tr = [0u8; SEEDBYTES];
    tr.copy_from_slice(&tr_vec);

    let sk = SecretKey {
        rho,
        key,
        tr,
        s1,
        s2,
        t0,
    };

    (pk, sk)
}

/// Sign a message with the secret key
pub fn sign(sk: &SecretKey, msg: &[u8]) -> Signature {
    // 1. Expand matrix A
    let a_hat = MatrixA::expand(&sk.rho);

    // 2. Compute mu = CRH(tr || msg)
    let mut mu_input = Vec::with_capacity(SEEDBYTES + msg.len());
    mu_input.extend_from_slice(&sk.tr);
    mu_input.extend_from_slice(msg);
    let mu = crh(&mu_input, CRHBYTES);

    // 3. Compute rho' = CRH(key || mu) (deterministic signing)
    let mut rho_prime_input = Vec::with_capacity(SEEDBYTES + CRHBYTES);
    rho_prime_input.extend_from_slice(&sk.key);
    rho_prime_input.extend_from_slice(&mu);
    let rho_prime = crh(&rho_prime_input, CRHBYTES);

    // 4. Pre-compute NTT of s1, s2, t0
    let mut s1_hat = sk.s1.clone();
    s1_hat.ntt();
    let mut s2_hat = sk.s2.clone();
    s2_hat.ntt();
    let mut t0_hat = sk.t0.clone();
    t0_hat.ntt();

    // 5. Rejection sampling loop
    let mut kappa: u16 = 0;
    loop {
        // Sample masking vector y
        let mut y = PolyVecL::new();
        for i in 0..L {
            y.polys[i] = Poly::sample_gamma1(&rho_prime, kappa + i as u16);
        }
        kappa += L as u16;

        // w = A * NTT(y)
        let mut y_hat = y.clone();
        y_hat.ntt();
        let mut w = a_hat.mul_ntt(&y_hat);
        w.invntt();
        w.caddq();

        // Decompose w into w1, w0
        let (w1, w0) = w.decompose();

        // Hash w1: c_tilde = H(mu || pack(w1))
        let w1_packed = w1.pack_w1();
        let mut c_tilde_input = Vec::with_capacity(CRHBYTES + w1_packed.len());
        c_tilde_input.extend_from_slice(&mu);
        c_tilde_input.extend_from_slice(&w1_packed);
        let c_tilde_vec = crh(&c_tilde_input, CTILDEBYTES);
        let mut c_tilde = [0u8; CTILDEBYTES];
        c_tilde.copy_from_slice(&c_tilde_vec);

        // Sample challenge c from c_tilde
        let cp = Poly::sample_challenge(&c_tilde);
        let mut cp_hat = cp.clone();
        cp_hat.ntt();

        // Compute z = y + c*s1
        let mut z = PolyVecL::new();
        for i in 0..L {
            z.polys[i] = Poly::pointwise_montgomery(&cp_hat, &s1_hat.polys[i]);
            z.polys[i].invntt();
            z.polys[i].add_assign(&y.polys[i]);
            z.polys[i].reduce();
        }

        // Check ||z||_inf < gamma1 - beta
        if !z.chknorm(GAMMA1 - BETA) {
            continue;
        }

        // Compute w0 - c*s2
        let mut cs2 = PolyVecK::new();
        for i in 0..K {
            cs2.polys[i] = Poly::pointwise_montgomery(&cp_hat, &s2_hat.polys[i]);
            cs2.polys[i].invntt();
        }
        let mut r0 = w0.clone();
        r0.sub_assign(&cs2);
        r0.reduce();

        // Check ||w0 - c*s2||_inf < gamma2 - beta
        if !r0.chknorm(GAMMA2 - BETA) {
            continue;
        }

        // Compute c*t0
        let mut ct0 = PolyVecK::new();
        for i in 0..K {
            ct0.polys[i] = Poly::pointwise_montgomery(&cp_hat, &t0_hat.polys[i]);
            ct0.polys[i].invntt();
        }
        ct0.caddq();

        // Check ||c*t0||_inf < gamma2
        if !ct0.chknorm(GAMMA2) {
            continue;
        }

        // Compute hints: MakeHint checks if adding ct0 changes HighBits of (w - cs2)
        // w_minus_cs2 = w - cs2, w_minus_cs2_plus_ct0 = w - cs2 + ct0
        let mut w_minus_cs2 = w.clone();
        w_minus_cs2.sub_assign(&cs2);
        w_minus_cs2.caddq();

        let mut w_minus_cs2_plus_ct0 = w_minus_cs2.clone();
        w_minus_cs2_plus_ct0.add_assign(&ct0);
        w_minus_cs2_plus_ct0.caddq();

        let mut h = PolyVecK::new();
        let mut hint_count = 0usize;
        for i in 0..K {
            for j in 0..N {
                let hi_with = Poly::highbits(w_minus_cs2_plus_ct0.polys[i].coeffs[j]);
                let hi_without = Poly::highbits(w_minus_cs2.polys[i].coeffs[j]);
                if hi_with != hi_without {
                    h.polys[i].coeffs[j] = 1;
                    hint_count += 1;
                }
            }
        }

        if hint_count > OMEGA {
            continue;
        }

        return Signature { c_tilde, z, h };
    }
}

/// Verify a signature against a public key and message
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    // 1. Check z norm
    if !sig.z.chknorm(GAMMA1 - BETA) {
        return false;
    }

    // 2. Expand matrix A
    let a_hat = MatrixA::expand(&pk.rho);

    // 3. Compute mu = CRH(tr || msg) where tr = CRH(pk)
    let pk_bytes = pk.to_bytes();
    let tr = crh(&pk_bytes, SEEDBYTES);
    let mut mu_input = Vec::with_capacity(SEEDBYTES + msg.len());
    mu_input.extend_from_slice(&tr);
    mu_input.extend_from_slice(msg);
    let mu = crh(&mu_input, CRHBYTES);

    // 4. Recover challenge c
    let cp = Poly::sample_challenge(&sig.c_tilde);
    let mut cp_hat = cp.clone();
    cp_hat.ntt();

    // 5. Compute A*z - c*t1*2^d in NTT domain
    let mut z_hat = sig.z.clone();
    z_hat.ntt();
    let mut az = a_hat.mul_ntt(&z_hat);

    // Compute c * t1 * 2^d
    let mut t1_shifted = pk.t1.clone();
    for i in 0..K {
        t1_shifted.polys[i].shiftl(); // multiply by 2^d
    }
    t1_shifted.ntt();
    let mut ct1 = PolyVecK::new();
    for i in 0..K {
        ct1.polys[i] = Poly::pointwise_montgomery(&cp_hat, &t1_shifted.polys[i]);
    }

    // w' = A*z - c*t1*2^d
    az.sub_assign(&ct1);
    az.invntt();
    az.caddq();

    // 6. UseHint to recover w1'
    let w1_prime = PolyVecK::use_hint(&az, &sig.h);

    // 7. Recompute c_tilde' = H(mu || pack(w1'))
    let w1_packed = w1_prime.pack_w1();
    let mut c_tilde_input = Vec::with_capacity(CRHBYTES + w1_packed.len());
    c_tilde_input.extend_from_slice(&mu);
    c_tilde_input.extend_from_slice(&w1_packed);
    let c_tilde_prime = crh(&c_tilde_input, CTILDEBYTES);

    // 8. Check c_tilde == c_tilde'
    use subtle::ConstantTimeEq;
    let eq = sig.c_tilde[..].ct_eq(&c_tilde_prime[..CTILDEBYTES]);
    eq.into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_keygen_sign_verify() {
        let (pk, sk) = keygen();
        let msg = b"Hello, post-quantum world!";
        let sig = sign(&sk, msg);
        assert!(verify(&pk, msg, &sig), "valid signature should verify");
    }

    #[test]
    fn test_wrong_message_fails() {
        let (pk, sk) = keygen();
        let msg = b"Hello, post-quantum world!";
        let sig = sign(&sk, msg);
        let wrong_msg = b"Wrong message";
        assert!(!verify(&pk, wrong_msg, &sig), "wrong message should not verify");
    }

    #[test]
    fn test_wrong_key_fails() {
        let (_pk1, sk1) = keygen();
        let (pk2, _sk2) = keygen();
        let msg = b"Test message";
        let sig = sign(&sk1, msg);
        assert!(!verify(&pk2, msg, &sig), "wrong public key should not verify");
    }

    #[test]
    fn test_signature_serialization() {
        let (pk, sk) = keygen();
        let msg = b"Serialization test";
        let sig = sign(&sk, msg);
        let sig_bytes = sig.to_bytes();
        assert_eq!(sig_bytes.len(), SIG_BYTES);
        let sig2 = Signature::from_bytes(&sig_bytes).expect("deserialization should succeed");
        assert!(verify(&pk, msg, &sig2), "deserialized signature should verify");
    }

    #[test]
    fn test_pk_serialization() {
        let (pk, sk) = keygen();
        let pk_bytes = pk.to_bytes();
        assert_eq!(pk_bytes.len(), PK_BYTES);
        let pk2 = PublicKey::from_bytes(&pk_bytes).expect("pk deserialization should succeed");
        let msg = b"PK serialization test";
        let sig = sign(&sk, msg);
        assert!(verify(&pk2, msg, &sig), "deserialized pk should verify");
    }
}
