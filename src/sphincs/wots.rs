/// WOTS+ one-time signature scheme.
///
/// WOTS+ with Winternitz parameter w=16 signs messages by
/// applying chains of hash function calls.

use crate::sphincs::address::{Address, WOTS_PK};
use crate::sphincs::hash::{hash_f, prf};
use crate::sphincs::params::*;

/// Compute one step of the WOTS+ chain.
/// chain(pk_seed, adrs, x, start, steps)
/// Applies hash_f `steps` times starting from index `start`.
fn chain(
    pk_seed: &[u8; N],
    adrs: &mut Address,
    input: &[u8; N],
    start: u32,
    steps: u32,
) -> [u8; N] {
    let mut tmp = *input;
    for i in start..(start + steps) {
        adrs.set_hash_address(i);
        tmp = hash_f(pk_seed, adrs, &tmp);
    }
    tmp
}

/// Convert message bytes to base-w representation.
/// Each base-w digit is log2(w)=4 bits.
fn base_w(msg: &[u8], out_len: usize) -> Vec<u32> {
    let mut result = Vec::with_capacity(out_len);
    let mut bits = 0u32;
    let mut total = 0u32;
    let mut idx = 0;

    for _ in 0..out_len {
        if bits == 0 {
            total = msg[idx] as u32;
            idx += 1;
            bits = 8;
        }
        bits -= WOTS_LOG_W as u32;
        result.push((total >> bits) & ((WOTS_W as u32) - 1));
    }
    result
}

/// Compute WOTS+ checksum and append to base-w message.
fn wots_checksum(msg_base_w: &[u32]) -> Vec<u32> {
    let mut csum: u32 = 0;
    for &val in msg_base_w {
        csum += (WOTS_W as u32 - 1) - val;
    }
    // csum <<= (8 - ((WOTS_LEN2 * WOTS_LOG_W) % 8)) % 8
    let shift = (8 - ((WOTS_LEN2 * WOTS_LOG_W) % 8)) % 8;
    csum <<= shift;

    // Convert checksum to base-w (need WOTS_LEN2 digits)
    let csum_bytes = csum.to_be_bytes();
    // We need ceil(WOTS_LEN2 * WOTS_LOG_W / 8) bytes
    let needed_bytes = (WOTS_LEN2 * WOTS_LOG_W + 7) / 8;
    let start = csum_bytes.len() - needed_bytes;
    let csum_base_w = base_w(&csum_bytes[start..], WOTS_LEN2);

    let mut full = msg_base_w.to_vec();
    full.extend_from_slice(&csum_base_w);
    full
}

/// Generate a WOTS+ public key.
pub fn wots_pk_gen(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    adrs: &mut Address,
) -> [u8; N] {
    let mut pk_adrs = adrs.clone();
    pk_adrs.set_type(WOTS_PK);
    pk_adrs.set_keypair_address(adrs.get_keypair_address());

    let mut tmp = Vec::with_capacity(WOTS_LEN * N);

    for i in 0..WOTS_LEN {
        adrs.set_chain_address(i as u32);
        adrs.set_hash_address(0);
        // sk = PRF(pk_seed, sk_seed, adrs)
        let sk = prf(pk_seed, sk_seed, adrs);
        // chain to the top: w-1 steps from 0
        let pk_i = chain(pk_seed, adrs, &sk, 0, (WOTS_W - 1) as u32);
        tmp.extend_from_slice(&pk_i);
    }

    // Compress: T_l(pk_seed, pk_adrs, tmp)
    hash_f(pk_seed, &pk_adrs, &tmp)
}

/// Sign a message digest (n bytes) using WOTS+.
pub fn wots_sign(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    msg: &[u8; N],
    adrs: &mut Address,
) -> Vec<u8> {
    let msg_base_w = base_w(msg, WOTS_LEN1);
    let lengths = wots_checksum(&msg_base_w);

    let mut sig = Vec::with_capacity(WOTS_SIG_BYTES);

    for i in 0..WOTS_LEN {
        adrs.set_chain_address(i as u32);
        adrs.set_hash_address(0);
        let sk = prf(pk_seed, sk_seed, adrs);
        let sig_i = chain(pk_seed, adrs, &sk, 0, lengths[i]);
        sig.extend_from_slice(&sig_i);
    }

    sig
}

/// Compute WOTS+ public key from a signature.
pub fn wots_pk_from_sig(
    pk_seed: &[u8; N],
    sig: &[u8],
    msg: &[u8; N],
    adrs: &mut Address,
) -> [u8; N] {
    let msg_base_w = base_w(msg, WOTS_LEN1);
    let lengths = wots_checksum(&msg_base_w);

    let mut pk_adrs = adrs.clone();
    pk_adrs.set_type(WOTS_PK);
    pk_adrs.set_keypair_address(adrs.get_keypair_address());

    let mut tmp = Vec::with_capacity(WOTS_LEN * N);

    for i in 0..WOTS_LEN {
        adrs.set_chain_address(i as u32);
        adrs.set_hash_address(0);
        let mut sig_i = [0u8; N];
        sig_i.copy_from_slice(&sig[i * N..(i + 1) * N]);
        let pk_i = chain(
            pk_seed,
            adrs,
            &sig_i,
            lengths[i],
            (WOTS_W as u32 - 1) - lengths[i],
        );
        tmp.extend_from_slice(&pk_i);
    }

    hash_f(pk_seed, &pk_adrs, &tmp)
}
