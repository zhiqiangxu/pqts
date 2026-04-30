/// SPHINCS+ key generation, signing, and verification.

use crate::sphincs::address::{Address, FORS_TREE};
use crate::sphincs::hash::{hash_msg, prf_msg};
use crate::sphincs::params::*;
use crate::sphincs::{fors, xmss};
use rand::RngCore;

/// SPHINCS+ public key: PK.seed || PK.root
#[derive(Clone, Debug)]
pub struct PublicKey {
    pub seed: [u8; N],
    pub root: [u8; N],
}

impl PublicKey {
    pub fn to_bytes(&self) -> [u8; PK_BYTES] {
        let mut out = [0u8; PK_BYTES];
        out[..N].copy_from_slice(&self.seed);
        out[N..].copy_from_slice(&self.root);
        out
    }
}

/// SPHINCS+ secret key
#[derive(Clone, Debug)]
pub struct SecretKey {
    pub sk_seed: [u8; N],
    pub sk_prf: [u8; N],
    pub pk: PublicKey,
}

impl SecretKey {
    pub fn to_bytes(&self) -> [u8; SK_BYTES] {
        let mut out = [0u8; SK_BYTES];
        out[..N].copy_from_slice(&self.sk_seed);
        out[N..2 * N].copy_from_slice(&self.sk_prf);
        out[2 * N..3 * N].copy_from_slice(&self.pk.seed);
        out[3 * N..].copy_from_slice(&self.pk.root);
        out
    }
}

/// SPHINCS+ signature
#[derive(Clone, Debug)]
pub struct Signature {
    pub randomness: [u8; N],
    pub fors_sig: Vec<u8>,
    pub ht_sig: Vec<u8>,
}

impl Signature {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(SIG_BYTES);
        out.extend_from_slice(&self.randomness);
        out.extend_from_slice(&self.fors_sig);
        out.extend_from_slice(&self.ht_sig);
        out
    }

    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != SIG_BYTES {
            return None;
        }
        let mut randomness = [0u8; N];
        randomness.copy_from_slice(&bytes[..N]);
        let fors_sig = bytes[N..N + FORS_SIG_BYTES].to_vec();
        let ht_sig = bytes[N + FORS_SIG_BYTES..].to_vec();
        Some(Signature { randomness, fors_sig, ht_sig })
    }
}

fn extract_tree_and_leaf(digest: &[u8]) -> (u64, u32) {
    let fors_end = FORS_MSG_BYTES;

    let mut tree: u64 = 0;
    for i in 0..TREE_BYTES {
        tree = (tree << 8) | (digest[fors_end + i] as u64);
    }
    tree &= (1u64 << TREE_BITS) - 1;

    let leaf_start = fors_end + TREE_BYTES;
    let mut leaf: u32 = 0;
    for i in 0..LEAF_BYTES {
        leaf = (leaf << 8) | (digest[leaf_start + i] as u32);
    }
    leaf &= (1u32 << LEAF_BITS) - 1;

    (tree, leaf)
}

/// Create a fresh Address for a given layer and tree.
fn make_adrs(layer: u32, tree: u64) -> Address {
    let mut adrs = Address::new();
    adrs.set_layer_address(layer);
    adrs.set_tree_address(tree);
    adrs
}

/// Generate a SPHINCS+ keypair.
pub fn keygen() -> (PublicKey, SecretKey) {
    let mut rng = rand::thread_rng();

    let mut sk_seed = [0u8; N];
    let mut sk_prf = [0u8; N];
    let mut pk_seed = [0u8; N];

    rng.fill_bytes(&mut sk_seed);
    rng.fill_bytes(&mut sk_prf);
    rng.fill_bytes(&mut pk_seed);

    let root = ht_root(&pk_seed, &sk_seed);

    let pk = PublicKey { seed: pk_seed, root };
    let sk = SecretKey { sk_seed, sk_prf, pk: pk.clone() };

    (pk, sk)
}

/// Compute the hypertree root.
fn ht_root(pk_seed: &[u8; N], sk_seed: &[u8; N]) -> [u8; N] {
    let mut adrs = make_adrs((D - 1) as u32, 0);
    xmss::xmss_root(pk_seed, sk_seed, &mut adrs)
}

/// Sign a message using SPHINCS+.
pub fn sign(sk: &SecretKey, msg: &[u8]) -> Signature {
    let pk_seed = &sk.pk.seed;
    let sk_seed = &sk.sk_seed;

    let opt_rand = sk.pk.seed;
    let r = prf_msg(&sk.sk_prf, &opt_rand, msg);
    let digest = hash_msg(&r, pk_seed, &sk.pk.root, msg);

    let fors_msg = &digest[..FORS_MSG_BYTES];
    let (tree_idx, leaf_idx) = extract_tree_and_leaf(&digest);

    // FORS signature
    let mut fors_adrs = Address::new();
    fors_adrs.set_tree_address(tree_idx);
    fors_adrs.set_type(FORS_TREE);
    fors_adrs.set_keypair_address(leaf_idx);
    let fors_sig = fors::fors_sign(pk_seed, sk_seed, fors_msg, &mut fors_adrs);

    // FORS public key
    let mut fors_adrs2 = Address::new();
    fors_adrs2.set_tree_address(tree_idx);
    fors_adrs2.set_type(FORS_TREE);
    fors_adrs2.set_keypair_address(leaf_idx);
    let fors_pk = fors::fors_pk_from_sig(pk_seed, &fors_sig, fors_msg, &mut fors_adrs2);

    // Hypertree signature
    let ht_sig = ht_sign(pk_seed, sk_seed, &fors_pk, tree_idx, leaf_idx);

    Signature { randomness: r, fors_sig, ht_sig }
}

/// Verify a SPHINCS+ signature.
pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    let pk_seed = &pk.seed;

    let digest = hash_msg(&sig.randomness, pk_seed, &pk.root, msg);
    let fors_msg = &digest[..FORS_MSG_BYTES];
    let (tree_idx, leaf_idx) = extract_tree_and_leaf(&digest);

    let mut fors_adrs = Address::new();
    fors_adrs.set_tree_address(tree_idx);
    fors_adrs.set_type(FORS_TREE);
    fors_adrs.set_keypair_address(leaf_idx);
    let fors_pk = fors::fors_pk_from_sig(pk_seed, &sig.fors_sig, fors_msg, &mut fors_adrs);

    ht_verify(pk_seed, &pk.root, &fors_pk, tree_idx, leaf_idx, &sig.ht_sig)
}

/// Sign with the hypertree (d layers of XMSS).
fn ht_sign(
    pk_seed: &[u8; N],
    sk_seed: &[u8; N],
    msg: &[u8; N],
    tree: u64,
    leaf: u32,
) -> Vec<u8> {
    let mut sig = Vec::with_capacity(HT_SIG_BYTES);

    // Sign at layer 0
    let mut adrs = make_adrs(0, tree);
    let xmss_sig = xmss::xmss_sign(pk_seed, sk_seed, msg, leaf, &mut adrs);
    sig.extend_from_slice(&xmss_sig);

    let mut adrs = make_adrs(0, tree);
    let mut root = xmss::xmss_root_from_sig(pk_seed, &xmss_sig, msg, leaf, &mut adrs);

    // Sign at layers 1..D-1
    let mut current_tree = tree;
    for layer in 1..D {
        let current_leaf = (current_tree & ((1u64 << TREE_HEIGHT) - 1)) as u32;
        current_tree >>= TREE_HEIGHT;

        let mut adrs = make_adrs(layer as u32, current_tree);
        let xmss_sig = xmss::xmss_sign(pk_seed, sk_seed, &root, current_leaf, &mut adrs);
        sig.extend_from_slice(&xmss_sig);

        if layer < D - 1 {
            let mut adrs = make_adrs(layer as u32, current_tree);
            root = xmss::xmss_root_from_sig(
                pk_seed, &xmss_sig, &root, current_leaf, &mut adrs,
            );
        }
    }

    sig
}

/// Verify a hypertree signature.
fn ht_verify(
    pk_seed: &[u8; N],
    pk_root: &[u8; N],
    msg: &[u8; N],
    tree: u64,
    leaf: u32,
    sig: &[u8],
) -> bool {
    let mut adrs = make_adrs(0, tree);
    let xmss_sig = &sig[..XMSS_SIG_BYTES];
    let mut root = xmss::xmss_root_from_sig(pk_seed, xmss_sig, msg, leaf, &mut adrs);

    let mut current_tree = tree;
    for layer in 1..D {
        let current_leaf = (current_tree & ((1u64 << TREE_HEIGHT) - 1)) as u32;
        current_tree >>= TREE_HEIGHT;

        let mut adrs = make_adrs(layer as u32, current_tree);
        let xmss_sig = &sig[layer * XMSS_SIG_BYTES..(layer + 1) * XMSS_SIG_BYTES];
        root = xmss::xmss_root_from_sig(pk_seed, xmss_sig, &root, current_leaf, &mut adrs);
    }

    root == *pk_root
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xmss_sign_verify_single() {
        use crate::sphincs::xmss;

        let pk_seed = [1u8; N];
        let sk_seed = [2u8; N];
        let msg = [3u8; N];

        let mut adrs = make_adrs(0, 0);
        let root = xmss::xmss_root(&pk_seed, &sk_seed, &mut adrs);

        let mut adrs = make_adrs(0, 0);
        let sig = xmss::xmss_sign(&pk_seed, &sk_seed, &msg, 0, &mut adrs);

        let mut adrs = make_adrs(0, 0);
        let recovered = xmss::xmss_root_from_sig(&pk_seed, &sig, &msg, 0, &mut adrs);
        assert_eq!(root, recovered, "XMSS root mismatch");
    }

    #[test]
    fn test_fors_sign_verify() {
        use crate::sphincs::fors;

        let pk_seed = [1u8; N];
        let sk_seed = [2u8; N];
        let msg = [0x42u8; FORS_MSG_BYTES];

        let mut adrs = Address::new();
        adrs.set_type(FORS_TREE);
        adrs.set_keypair_address(0);

        let sig = fors::fors_sign(&pk_seed, &sk_seed, &msg, &mut adrs);

        let mut adrs2 = Address::new();
        adrs2.set_type(FORS_TREE);
        adrs2.set_keypair_address(0);
        let pk1 = fors::fors_pk_from_sig(&pk_seed, &sig, &msg, &mut adrs2);

        let mut adrs3 = Address::new();
        adrs3.set_type(FORS_TREE);
        adrs3.set_keypair_address(0);
        let sig2 = fors::fors_sign(&pk_seed, &sk_seed, &msg, &mut adrs3);

        let mut adrs4 = Address::new();
        adrs4.set_type(FORS_TREE);
        adrs4.set_keypair_address(0);
        let pk2 = fors::fors_pk_from_sig(&pk_seed, &sig2, &msg, &mut adrs4);

        assert_eq!(pk1, pk2, "FORS pk mismatch");
    }

    #[test]
    fn test_ht_sign_verify() {
        let pk_seed = [1u8; N];
        let sk_seed = [2u8; N];
        let msg = [3u8; N];

        let root = ht_root(&pk_seed, &sk_seed);

        for &(tree_idx, leaf_idx) in &[(0u64, 0u32), (0, 5), (1, 3)] {
            let ht_sig = ht_sign(&pk_seed, &sk_seed, &msg, tree_idx, leaf_idx);
            let ok = ht_verify(&pk_seed, &root, &msg, tree_idx, leaf_idx, &ht_sig);
            assert!(ok, "hypertree verify failed for tree={}, leaf={}", tree_idx, leaf_idx);
        }
    }

    #[test]
    fn test_xmss_sign_verify_nonzero_leaf() {
        use crate::sphincs::xmss;

        let pk_seed = [1u8; N];
        let sk_seed = [2u8; N];
        let msg = [3u8; N];

        // Test with leaf=370, layer=6, tree=0 (same as the failing case)
        let mut adrs = make_adrs(6, 0);
        let root = xmss::xmss_root(&pk_seed, &sk_seed, &mut adrs);

        let mut adrs = make_adrs(6, 0);
        let sig = xmss::xmss_sign(&pk_seed, &sk_seed, &msg, 370, &mut adrs);

        let mut adrs = make_adrs(6, 0);
        let recovered = xmss::xmss_root_from_sig(&pk_seed, &sig, &msg, 370, &mut adrs);
        assert_eq!(root, recovered, "XMSS root mismatch for leaf=370");
    }

    #[test]
    fn test_param_sizes() {
        assert_eq!(TREE_HEIGHT, 9);
        assert_eq!(WOTS_LEN, 35);
        assert_eq!(WOTS_SIG_BYTES, 35 * 16);
        assert_eq!(XMSS_SIG_BYTES, 560 + 9 * 16);
        assert_eq!(HT_SIG_BYTES, 7 * 704);
        assert_eq!(FORS_SIG_BYTES, 14 * 13 * 16);
        assert_eq!(SIG_BYTES, 16 + 2912 + 4928);
        assert_eq!(PK_BYTES, 32);
    }
}
