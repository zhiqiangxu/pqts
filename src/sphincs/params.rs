/// SPHINCS+-128s-simple (SHA-256) parameters

/// Security parameter / hash output length in bytes
pub const N: usize = 16;

/// Total tree height
pub const FULL_HEIGHT: usize = 63;

/// Number of hypertree layers
pub const D: usize = 7;

/// Height of each XMSS tree (FULL_HEIGHT / D)
pub const TREE_HEIGHT: usize = FULL_HEIGHT / D; // 9

/// FORS: number of trees
pub const FORS_TREES: usize = 14; // k

/// FORS: height of each tree
pub const FORS_HEIGHT: usize = 12; // a

/// Winternitz parameter
pub const WOTS_W: usize = 16;

/// WOTS+ log2(w)
pub const WOTS_LOG_W: usize = 4;

/// WOTS+ chain length (len1)
/// len1 = ceil(8n / log2(w)) = ceil(128 / 4) = 32
pub const WOTS_LEN1: usize = 2 * N; // 32

/// WOTS+ checksum length (len2)
/// len2 = floor(log2(len1 * (w-1)) / log2(w)) + 1
/// = floor(log2(32 * 15) / 4) + 1 = floor(log2(480)/4)+1 = floor(8.9/4)+1 = 3
pub const WOTS_LEN2: usize = 3;

/// WOTS+ total length
pub const WOTS_LEN: usize = WOTS_LEN1 + WOTS_LEN2; // 35

/// WOTS+ signature size in bytes
pub const WOTS_SIG_BYTES: usize = WOTS_LEN * N;

/// FORS signature size: k * (a * n + n) = k * (a+1) * n
pub const FORS_SIG_BYTES: usize = FORS_TREES * (FORS_HEIGHT + 1) * N;

/// Number of FORS message bits: k * a
pub const FORS_MSG_BYTES: usize = (FORS_TREES * FORS_HEIGHT + 7) / 8;

/// XMSS signature size: WOTS sig + h' * n auth path
pub const XMSS_SIG_BYTES: usize = WOTS_SIG_BYTES + TREE_HEIGHT * N;

/// Hypertree signature size: d * XMSS sig
pub const HT_SIG_BYTES: usize = D * XMSS_SIG_BYTES;

/// Total signature size: randomness(n) + FORS sig + HT sig
pub const SIG_BYTES: usize = N + FORS_SIG_BYTES + HT_SIG_BYTES;

/// Public key size: 2 * n (seed + root)
pub const PK_BYTES: usize = 2 * N;

/// Secret key size: 2 * n (seed + prf key) + PK
pub const SK_BYTES: usize = 2 * N + PK_BYTES;

/// Number of bytes needed for tree index in message digest
pub const TREE_BITS: usize = FULL_HEIGHT - TREE_HEIGHT; // 54
pub const TREE_BYTES: usize = (TREE_BITS + 7) / 8; // 7

/// Number of bytes for leaf index
pub const LEAF_BITS: usize = TREE_HEIGHT; // 9
pub const LEAF_BYTES: usize = (LEAF_BITS + 7) / 8; // 2

/// Digest bytes needed: FORS_MSG_BYTES + TREE_BYTES + LEAF_BYTES
pub const DIGEST_BYTES: usize = FORS_MSG_BYTES + TREE_BYTES + LEAF_BYTES;
