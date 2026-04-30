/// SHA-256-based hash functions for SPHINCS+ simple variant.
///
/// In the "simple" instantiation, tweakable hash functions are implemented as:
///   F(PK.seed, ADRS, M) = SHA-256(PK.seed || ADRS || M), truncated to n bytes
///   H(PK.seed, ADRS, M) = SHA-256(PK.seed || ADRS || M), truncated to n bytes
///   T(PK.seed, ADRS, M) = SHA-256(PK.seed || ADRS || M), truncated to n bytes
///   PRF(PK.seed, SK.seed, ADRS) = SHA-256(PK.seed || ADRS || SK.seed), truncated to n
///   PRFmsg(SK.prf, opt_rand, msg) = HMAC-SHA-256(SK.prf, opt_rand || msg), truncated to n
///   Hmsg(R, PK.seed, PK.root, M) = MGF1-SHA-256 based

use sha2::{Sha256, Digest};
use crate::sphincs::address::Address;
use crate::sphincs::params::*;

/// Truncate a SHA-256 digest to n bytes
fn truncate(hash: &[u8]) -> [u8; N] {
    let mut out = [0u8; N];
    out.copy_from_slice(&hash[..N]);
    out
}

/// F / H / T_l: tweakable hash (simple variant)
/// SHA-256(PK.seed || ADRS || M), truncated to n bytes
pub fn hash_f(pk_seed: &[u8; N], adrs: &Address, msg: &[u8]) -> [u8; N] {
    let mut hasher = Sha256::new();
    // Pad pk_seed to 32 bytes (SHA-256 block alignment for simple variant)
    let mut padded_seed = [0u8; 32];
    padded_seed[32 - N..].copy_from_slice(pk_seed);
    hasher.update(&padded_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(msg);
    let result = hasher.finalize();
    truncate(&result)
}

/// PRF: pseudo-random function
/// SHA-256(PK.seed || ADRS || SK.seed), truncated to n bytes
pub fn prf(pk_seed: &[u8; N], sk_seed: &[u8; N], adrs: &Address) -> [u8; N] {
    let mut hasher = Sha256::new();
    let mut padded_seed = [0u8; 32];
    padded_seed[32 - N..].copy_from_slice(pk_seed);
    hasher.update(&padded_seed);
    hasher.update(adrs.as_bytes());
    hasher.update(sk_seed);
    let result = hasher.finalize();
    truncate(&result)
}

/// PRF_msg: message PRF for randomized hashing
/// HMAC-SHA-256(SK.prf, opt_rand || msg), truncated
pub fn prf_msg(sk_prf: &[u8; N], opt_rand: &[u8; N], msg: &[u8]) -> [u8; N] {
    // HMAC-SHA-256 with sk_prf as key
    let block_size = 64usize;
    let mut key_pad = [0u8; 64];
    // sk_prf is N=16 bytes, which is less than block size, so just copy
    key_pad[..N].copy_from_slice(sk_prf);

    // ipad
    let mut ipad = [0x36u8; 64];
    for i in 0..block_size {
        ipad[i] ^= key_pad[i];
    }

    // opad
    let mut opad = [0x5cu8; 64];
    for i in 0..block_size {
        opad[i] ^= key_pad[i];
    }

    // inner hash: SHA-256(ipad || opt_rand || msg)
    let mut inner = Sha256::new();
    inner.update(&ipad);
    inner.update(opt_rand);
    inner.update(msg);
    let inner_hash = inner.finalize();

    // outer hash: SHA-256(opad || inner_hash)
    let mut outer = Sha256::new();
    outer.update(&opad);
    outer.update(&inner_hash);
    let result = outer.finalize();

    truncate(&result)
}

/// H_msg: hash message to produce FORS message digest + tree/leaf indices.
/// Uses MGF1-SHA-256 to produce enough bytes.
/// H_msg(R, PK.seed, PK.root, M) -> DIGEST_BYTES
pub fn hash_msg(
    r: &[u8; N],
    pk_seed: &[u8; N],
    pk_root: &[u8; N],
    msg: &[u8],
) -> Vec<u8> {
    // We need DIGEST_BYTES of output
    // Use SHA-256(R || PK.seed || PK.root || M) as the seed for MGF1
    let mut hasher = Sha256::new();
    hasher.update(r);
    hasher.update(pk_seed);
    hasher.update(pk_root);
    hasher.update(msg);
    let seed = hasher.finalize();

    mgf1_sha256(&seed, DIGEST_BYTES)
}

/// MGF1 with SHA-256
fn mgf1_sha256(seed: &[u8], length: usize) -> Vec<u8> {
    let mut output = Vec::new();
    let mut counter: u32 = 0;

    while output.len() < length {
        let mut hasher = Sha256::new();
        hasher.update(seed);
        hasher.update(&counter.to_be_bytes());
        let block = hasher.finalize();
        output.extend_from_slice(&block);
        counter += 1;
    }

    output.truncate(length);
    output
}
