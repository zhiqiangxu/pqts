# pqts

Post-quantum signature schemes implemented in Rust.

Three NIST-standardized post-quantum digital signature algorithms:

| Algorithm | Based On | Sig Size | PK Size | Paradigm |
|-----------|----------|----------|---------|----------|
| **Dilithium-2** (ML-DSA-44) | Module-LWE lattice | 2420 B | 1312 B | Fiat-Shamir |
| **Falcon-512** | NTRU lattice | 569 B | 896 B | Hash-and-Sign |
| **SPHINCS+** (128s-simple) | Hash functions (SHA-256) | 7856 B | 32 B | Stateless hash-based |

## Usage

```rust
use pqts::{dilithium, falcon, sphincs};

// Dilithium-2
let (pk, sk) = dilithium::keygen();
let sig = dilithium::sign(&sk, b"message");
assert!(dilithium::verify(&pk, b"message", &sig));

// Falcon-512
let (pk, sk) = falcon::keygen();
let sig = falcon::sign(&sk, b"message");
assert!(falcon::verify(&pk, b"message", &sig));

// SPHINCS+-128s
let (pk, sk) = sphincs::keygen();
let sig = sphincs::sign(&sk, b"message");
assert!(sphincs::verify(&pk, b"message", &sig));
```

## Algorithm Details

### Dilithium-2 (ML-DSA-44)

NIST's primary recommendation for post-quantum signatures. Based on the Module-LWE problem over polynomial rings.

- **Ring**: Z_q[X]/(X^256 + 1), q = 8380417
- **Matrix**: 4×4 module structure
- **NTT**: 256-point with Montgomery multiplication
- **Signing**: Fiat-Shamir with aborts (rejection sampling)
- Full key serialization/deserialization support

### Falcon-512

Compact lattice-based signatures built on NTRU. Chosen by Solana for its small signature size.

- **Ring**: Z_q[X]/(X^512 + 1), q = 12289
- **NTT**: 512-point negacyclic NTT
- **Key**: NTRU structure h = g·f⁻¹ mod q with short f, g
- **Verification**: Algebraic check s1 + s2·h = Hash(msg) mod q, plus norm bound

> **Note**: This is a simplified educational implementation. Real Falcon requires a full NTRU solver (computing F, G satisfying f·G − g·F = q) and FFT-based discrete Gaussian sampling for constant-time, optimally short signatures.

### SPHINCS+ (SLH-DSA, 128s-simple)

Stateless hash-based signature scheme. Conservative security — relies only on hash function security (SHA-256).

- **Parameters**: n=16, h=63, d=7, h'=9, k=14, a=12, w=16
- **Components**: WOTS+ (one-time sig) → XMSS (Merkle tree) → Hypertree (7 layers) + FORS (few-time sig)
- **Hash**: SHA-256 simple variant
- Largest signatures (~7856 B) but smallest public keys (~32 B)

## Build & Test

```bash
cargo build
cargo test
```

## Project Structure

```
src/
├── lib.rs
├── dilithium/
│   ├── ntt.rs          # NTT with Montgomery reduction
│   ├── params.rs       # ML-DSA-44 parameters
│   ├── poly.rs         # Polynomial operations, sampling, packing
│   └── sign.rs         # KeyGen, Sign, Verify
├── falcon/
│   ├── ntt.rs          # 512-point negacyclic NTT
│   ├── params.rs       # Falcon-512 parameters
│   ├── poly.rs         # Ring arithmetic
│   ├── sampler.rs      # Discrete Gaussian sampling
│   └── sign.rs         # KeyGen, Sign, Verify
└── sphincs/
    ├── address.rs      # ADRS structure
    ├── hash.rs         # SHA-256 hash wrappers
    ├── wots.rs         # WOTS+ one-time signature
    ├── xmss.rs         # XMSS Merkle tree
    ├── fors.rs         # FORS few-time signature
    ├── params.rs       # SPHINCS+-128s parameters
    └── sign.rs         # Full sign/verify with hypertree
```

## Disclaimer

This is an **educational implementation**. Do not use in production. It lacks:

- Constant-time guarantees (side-channel resistance)
- Proper Falcon NTRU key generation and FFT sampling
- Formal security auditing
- Optimized performance

For production use, refer to the reference implementations from [pq-crystals](https://github.com/pq-crystals) and [sphincs.org](https://sphincs.org).
