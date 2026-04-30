#![allow(dead_code)]

pub mod ntt;
pub mod params;
pub mod poly;
pub mod sign;

pub use sign::{keygen, sign, verify, PublicKey, SecretKey, Signature};
