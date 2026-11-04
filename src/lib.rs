//! Rust port of `c-fn-dsa-alt` (FN-DSA, i.e. the Falcon-based FIPS 206 draft
//! signature scheme; RAM-efficient alternate implementation).
//!
//! Current scope: **signature verification and signing** (including the
//! "weak" degrees 4..=256), plus the building blocks they need (SHAKE256,
//! codecs, arithmetic modulo q = 12289 with a runtime-selected AVX2 path,
//! floating-point FFT and Gaussian sampling), and **key generation**
//! (fixed-point FFT, big-integer NTRU solver).
//!
//! The public key / signature formats are identical to the C implementation.

mod codec;
mod fpoly;
mod fxp;
mod fxp_tables;
mod keccak_gen;
mod kgen;
mod kpoly;
mod mp31;
mod mq;
mod mqsign;
mod ntru;
mod primes;
mod sampler;
mod sha3;
mod sign;
mod sign_core;
mod sign_tables;
mod tables;
mod util;
mod vrfy;
mod zint31;

pub use util::{compute_mu, hashed_vrfykey_from_vrfykey};
pub use kgen::{keygen, keygen_seeded};
pub use sign::{sign, sign_seeded, sign_weak, sign_weak_seeded};
pub use vrfy::{verify, verify_weak};

/// Hash identifier for a raw (not pre-hashed) message.
pub const HASH_ID_RAW: &[u8] = b"";
/// DER-encoded OID for SHA-256.
pub const HASH_ID_SHA256: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x01";
/// DER-encoded OID for SHA-384.
pub const HASH_ID_SHA384: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x02";
/// DER-encoded OID for SHA-512.
pub const HASH_ID_SHA512: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x03";
/// DER-encoded OID for SHA-512/256.
pub const HASH_ID_SHA512_256: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x06";
/// DER-encoded OID for SHA3-256.
pub const HASH_ID_SHA3_256: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x08";
/// DER-encoded OID for SHA3-384.
pub const HASH_ID_SHA3_384: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x09";
/// DER-encoded OID for SHA3-512.
pub const HASH_ID_SHA3_512: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0A";
/// DER-encoded OID for SHAKE128.
pub const HASH_ID_SHAKE128: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0B";
/// DER-encoded OID for SHAKE256.
pub const HASH_ID_SHAKE256: &[u8] = b"\x06\x09\x60\x86\x48\x01\x65\x03\x04\x02\x0C";
/// "External mu": the `hv` argument is the 64-byte message representative.
pub const HASH_ID_EXTMU: &[u8] = b"\xFE";

/// Signing key size in bytes (2 <= logn <= 10).
pub const fn sign_key_size(logn: u32) -> usize {
    65 + ((12 - (logn >= 6) as usize - (logn >= 8) as usize - (logn >= 10) as usize) << (logn - 2))
}

/// Number of bits per coefficient of f and g in an encoded signing key.
pub(crate) const fn fg_nbits(logn: u32) -> u32 {
    match logn {
        2..=5 => 8,
        6 | 7 => 7,
        8 | 9 => 6,
        _ => 5,
    }
}

/// Verifying key size in bytes (2 <= logn <= 10).
pub const fn vrfy_key_size(logn: u32) -> usize {
    1 + (7usize << (logn - 2))
}

/// Signature size in bytes (2 <= logn <= 10).
pub const fn signature_size(logn: u32) -> usize {
    let s = 10 - logn;
    44 + 3 * (256usize >> s) + 2 * (128usize >> s) + 3 * (64usize >> s) + 2 * (16usize >> s)
        - 2 * (2usize >> s)
        - 8 * (1usize >> s)
}

/// Whether the AVX2 code paths are in use on this machine.
pub fn using_avx2() -> bool {
    mq::has_avx2()
}

#[doc(hidden)]
pub mod internals {
    //! Exposed only for benchmarks and tests; not a stable API.
    pub use crate::codec::{comp_decode, mqpoly_decode};
    pub use crate::mq::{avx2 as avx2_kernels, best as best_kernels, Kernels, PORTABLE};
    pub use crate::vrfy::verify_with;
    pub use crate::util::hash_to_point;
}
