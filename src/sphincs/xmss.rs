/// XMSS (eXtended Merkle Signature Scheme) tree.
///
/// Each XMSS tree has 2^TREE_HEIGHT leaves. Each leaf is a WOTS+ public key.

use crate::sphincs::address::{Address, TREE, WOTS_HASH};
use crate::sphincs::hash::hash_f;
use crate::sphincs::params::*;
use crate::sphincs::wots;

/// Compute the root of a subtree of WOTS+ public keys.
fn treehash(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    leaf_idx: u32,
    height: u32,
    adrs: &Address,
) -> [u8; N] {
    if height == 0 {
        let mut wots_adrs = adrs.clone();
        wots_adrs.set_type(WOTS_HASH);
        wots_adrs.set_keypair_address(leaf_idx);
        return wots::wots_pk_gen(pk_seed, sk_seed, &mut wots_adrs);
    }

    let left = treehash(pk_seed, sk_seed, leaf_idx, height - 1, adrs);
    let right = treehash(
        pk_seed,
        sk_seed,
        leaf_idx + (1 << (height - 1)),
        height - 1,
        adrs,
    );

    let mut tree_adrs = adrs.clone();
    tree_adrs.set_type(TREE);
    tree_adrs.set_tree_height(height);
    tree_adrs.set_tree_index(leaf_idx >> height);

    let mut combined = [0u8; 2 * N];
    combined[..N].copy_from_slice(&left);
    combined[N..].copy_from_slice(&right);

    hash_f(pk_seed, &tree_adrs, &combined)
}

/// Compute the XMSS root node.
pub fn xmss_root(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    adrs: &mut Address,
) -> [u8; N] {
    treehash(pk_seed, sk_seed, 0, TREE_HEIGHT as u32, adrs)
}

/// Generate an XMSS signature for a given message at a specific leaf index.
pub fn xmss_sign(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    msg: &[u8; N],
    idx: u32,
    adrs: &mut Address,
) -> Vec<u8> {
    // 1. Generate WOTS+ signature
    let mut wots_adrs = adrs.clone();
    wots_adrs.set_type(WOTS_HASH);
    wots_adrs.set_keypair_address(idx);
    let wots_sig = wots::wots_sign(pk_seed, sk_seed, msg, &mut wots_adrs);

    // 2. Compute authentication path
    let auth = compute_auth_path(pk_seed, sk_seed, idx, adrs);

    let mut sig = wots_sig;
    sig.extend_from_slice(&auth);
    sig
}

/// Compute the authentication path for leaf `idx`.
fn compute_auth_path(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    idx: u32,
    adrs: &Address,
) -> Vec<u8> {
    let mut auth = Vec::with_capacity(TREE_HEIGHT * N);

    for j in 0..TREE_HEIGHT as u32 {
        let sibling = (idx >> j) ^ 1;
        let node = treehash(pk_seed, sk_seed, sibling << j, j, adrs);
        auth.extend_from_slice(&node);
    }

    auth
}

/// Compute XMSS root from a signature (for verification).
pub fn xmss_root_from_sig(
    pk_seed: &[u8; N],
    sig: &[u8],
    msg: &[u8; N],
    idx: u32,
    adrs: &mut Address,
) -> [u8; N] {
    let wots_sig = &sig[..WOTS_SIG_BYTES];
    let auth = &sig[WOTS_SIG_BYTES..];

    // Recover WOTS+ public key from signature
    let mut wots_adrs = adrs.clone();
    wots_adrs.set_type(WOTS_HASH);
    wots_adrs.set_keypair_address(idx);
    let mut node = wots::wots_pk_from_sig(pk_seed, wots_sig, msg, &mut wots_adrs);

    // Walk up the tree using the authentication path
    for j in 0..TREE_HEIGHT as u32 {
        let mut tree_adrs = adrs.clone();
        tree_adrs.set_type(TREE);
        tree_adrs.set_tree_height(j + 1);
        tree_adrs.set_tree_index(idx >> (j + 1));

        let mut combined = [0u8; 2 * N];
        let auth_node = &auth[(j as usize) * N..(j as usize + 1) * N];

        if (idx >> j) & 1 == 0 {
            combined[..N].copy_from_slice(&node);
            combined[N..].copy_from_slice(auth_node);
        } else {
            combined[..N].copy_from_slice(auth_node);
            combined[N..].copy_from_slice(&node);
        }

        node = hash_f(pk_seed, &tree_adrs, &combined);
    }

    node
}
