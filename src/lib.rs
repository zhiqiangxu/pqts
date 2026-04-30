#![allow(clippy::needless_range_loop)]

pub mod dilithium;
pub mod falcon;
pub mod sphincs;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dilithium_sign_verify() {
        let (pk, sk) = dilithium::keygen();
        let msg = b"Hello, post-quantum world!";
        let sig = dilithium::sign(&sk, msg);
        assert!(dilithium::verify(&pk, msg, &sig));
        // Verify with wrong message fails
        assert!(!dilithium::verify(&pk, b"wrong message", &sig));
    }

    #[test]
    fn test_falcon_sign_verify() {
        let (pk, sk) = falcon::keygen();
        let msg = b"Hello, post-quantum world!";
        let sig = falcon::sign(&sk, msg);
        assert!(falcon::verify(&pk, msg, &sig));
        assert!(!falcon::verify(&pk, b"wrong message", &sig));
    }

    #[test]
    fn test_sphincs_sign_verify() {
        let (pk, sk) = sphincs::keygen();
        let msg = b"Hello, post-quantum world!";
        let sig = sphincs::sign(&sk, msg);
        assert!(sphincs::verify(&pk, msg, &sig));
        assert!(!sphincs::verify(&pk, b"wrong message", &sig));
    }
}
