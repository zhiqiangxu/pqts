/// SPHINCS+-128s-simple (SHA-256) post-quantum hash-based signature scheme.

pub mod params;
pub mod address;
pub mod hash;
pub mod wots;
pub mod xmss;
pub mod fors;
pub mod sign;

pub use sign::{PublicKey, SecretKey, Signature};

pub fn keygen() -> (PublicKey, SecretKey) {
    sign::keygen()
}

pub fn sign(sk: &SecretKey, msg: &[u8]) -> Signature {
    sign::sign(sk, msg)
}

pub fn verify(pk: &PublicKey, msg: &[u8], sig: &Signature) -> bool {
    sign::verify(pk, msg, sig)
}
