//! Falcon-512 post-quantum digital signature scheme.
//!
//! Falcon is a lattice-based signature scheme built on NTRU lattices,
//! using the "hash-and-sign" paradigm with a trapdoor sampler.
//!
//! This is a simplified educational implementation. It demonstrates the
//! core concepts of Falcon but should NOT be used in production:
//! - Key generation uses a simplified NTRU structure
//! - Signing uses rejection sampling instead of fast Fourier sampling
//! - No constant-time guarantees
//!
//! # Example
//!
//! ```no_run
//! use pqts::falcon::{keygen, sign, verify};
//!
//! let (pk, sk) = keygen();
//! let msg = b"Hello, post-quantum world!";
//! let sig = sign(&sk, msg);
//! assert!(verify(&pk, msg, &sig));
//! ```

pub mod ntt;
pub mod params;
pub mod poly;
pub mod sampler;

// The sign module contains keygen, sign, and verify implementations.
// We re-export the public API items for convenience.
// Note: `mod sign` (type namespace) and `use sign::sign` (value namespace) coexist.
pub mod sign;

pub use sign::{PublicKey, SecretKey, Signature};
pub use sign::{keygen, sign, verify};
