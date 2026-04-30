/// FORS (Forest of Random Subsets) few-time signature scheme.
///
/// FORS signs a message by revealing leaves from k binary trees,
/// each of height a. The message selects one leaf per tree.

use crate::sphincs::address::{Address, FORS_ROOTS};
use crate::sphincs::hash::{hash_f, prf};
use crate::sphincs::params::*;

/// Compute the value of a FORS leaf (secret value).
fn fors_sk_gen(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    adrs: &Address,
    idx: u32,
) -> [u8; N] {
    let mut sk_adrs = adrs.clone();
    sk_adrs.set_tree_height(0);
    sk_adrs.set_tree_index(idx);
    prf(pk_seed, sk_seed, &sk_adrs)
}

/// Compute a FORS tree node at given height and leaf offset.
/// leaf_idx is the index of the leftmost leaf of this subtree.
fn fors_treehash(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    leaf_idx: u32,
    height: u32,
    adrs: &Address,
) -> [u8; N] {
    if height == 0 {
        let sk = fors_sk_gen(pk_seed, sk_seed, adrs, leaf_idx);
        let mut leaf_adrs = adrs.clone();
        leaf_adrs.set_tree_height(0);
        leaf_adrs.set_tree_index(leaf_idx);
        return hash_f(pk_seed, &leaf_adrs, &sk);
    }

    let left = fors_treehash(pk_seed, sk_seed, leaf_idx, height - 1, adrs);
    let right = fors_treehash(
        pk_seed,
        sk_seed,
        leaf_idx + (1 << (height - 1)),
        height - 1,
        adrs,
    );

    let mut node_adrs = adrs.clone();
    node_adrs.set_tree_height(height);
    // At height h, the node index = leaf_idx / 2^h = leaf_idx >> h
    node_adrs.set_tree_index(leaf_idx >> height);

    let mut combined = [0u8; 2 * N];
    combined[..N].copy_from_slice(&left);
    combined[N..].copy_from_slice(&right);

    hash_f(pk_seed, &node_adrs, &combined)
}

/// Extract FORS message indices from the message digest.
/// Returns k indices, each in [0, 2^a).
fn message_to_indices(msg: &[u8]) -> Vec<u32> {
    let mut indices = Vec::with_capacity(FORS_TREES);
    let mut bits_consumed = 0usize;

    for _ in 0..FORS_TREES {
        let mut idx: u32 = 0;
        for _ in 0..FORS_HEIGHT {
            let byte_idx = bits_consumed / 8;
            let bit_idx = bits_consumed % 8;
            let bit = ((msg[byte_idx] >> (7 - bit_idx)) & 1) as u32;
            idx = (idx << 1) | bit;
            bits_consumed += 1;
        }
        indices.push(idx);
    }

    indices
}

/// Sign a message digest using FORS.
pub fn fors_sign(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    msg: &[u8],
    adrs: &mut Address,
) -> Vec<u8> {
    let indices = message_to_indices(msg);
    let mut sig = Vec::with_capacity(FORS_SIG_BYTES);

    for (i, &idx) in indices.iter().enumerate() {
        let tree_offset = (i as u32) * (1u32 << FORS_HEIGHT);

        // Include the secret leaf value
        let sk = fors_sk_gen(pk_seed, sk_seed, adrs, tree_offset + idx);
        sig.extend_from_slice(&sk);

        // Include authentication path
        for j in 0..FORS_HEIGHT as u32 {
            let sibling = (idx >> j) ^ 1;
            let node = fors_treehash(
                pk_seed,
                sk_seed,
                (tree_offset + (sibling << j)) as u32,
                j,
                adrs,
            );
            sig.extend_from_slice(&node);
        }
    }

    sig
}

/// Compute FORS public key from a signature.
pub fn fors_pk_from_sig(
    pk_seed: &[u8; N],
    sig: &[u8],
    msg: &[u8],
    adrs: &mut Address,
) -> [u8; N] {
    let indices = message_to_indices(msg);
    let mut roots = Vec::with_capacity(FORS_TREES * N);
    let mut sig_offset = 0;
    let keypair_addr = adrs.get_keypair_address();

    for (i, &idx) in indices.iter().enumerate() {
        let tree_offset = (i as u32) * (1u32 << FORS_HEIGHT);

        // Get secret leaf value from signature
        let mut sk = [0u8; N];
        sk.copy_from_slice(&sig[sig_offset..sig_offset + N]);
        sig_offset += N;

        // Hash the leaf value
        let mut leaf_adrs = adrs.clone();
        leaf_adrs.set_tree_height(0);
        leaf_adrs.set_tree_index(tree_offset + idx);
        let mut node = hash_f(pk_seed, &leaf_adrs, &sk);

        // Walk up using the authentication path
        for j in 0..FORS_HEIGHT as u32 {
            let mut auth_node = [0u8; N];
            auth_node.copy_from_slice(&sig[sig_offset..sig_offset + N]);
            sig_offset += N;

            let mut node_adrs = adrs.clone();
            node_adrs.set_tree_height(j + 1);
            // Tree index at height j+1: (tree_offset + idx) >> (j+1)
            node_adrs.set_tree_index((tree_offset + idx) >> (j + 1));

            let mut combined = [0u8; 2 * N];
            if (idx >> j) & 1 == 0 {
                combined[..N].copy_from_slice(&node);
                combined[N..].copy_from_slice(&auth_node);
            } else {
                combined[..N].copy_from_slice(&auth_node);
                combined[N..].copy_from_slice(&node);
            }

            node = hash_f(pk_seed, &node_adrs, &combined);
        }

        roots.extend_from_slice(&node);
    }

    // Compress all roots into the FORS public key
    let mut pk_adrs = adrs.clone();
    pk_adrs.set_type(FORS_ROOTS);
    pk_adrs.set_keypair_address(keypair_addr);

    hash_f(pk_seed, &pk_adrs, &roots)
}
