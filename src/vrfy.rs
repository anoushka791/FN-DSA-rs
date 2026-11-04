//! Signature verification.

use crate::codec::{comp_decode, mqpoly_decode};
use crate::mq::{self, Kernels};
use crate::util::{compute_mu, hash_to_point, hashed_vrfykey_from_vrfykey};
use crate::{signature_size, vrfy_key_size, HASH_ID_EXTMU};

/// Largest buffer needed: two polynomials of 1024 16-bit values.
const MAX_N: usize = 1024;


fn inner_verify(
    k: &Kernels,
    logn_min: u32,
    logn_max: u32,
    sig: &[u8],
    vrfy_key: &[u8],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
) -> bool {
    if sig.is_empty() || vrfy_key.is_empty() {
        return false;
    }
    let logn = vrfy_key[0] as u32;
    if logn < logn_min || logn > logn_max || sig[0] as u32 != 0x30 + logn {
        return false;
    }
    if sig.len() != signature_size(logn) || vrfy_key.len() != vrfy_key_size(logn) {
        return false;
    }
    let n = 1usize << logn;

    // Message representative mu.
    let mu: [u8; 64] = if id.first() == Some(&HASH_ID_EXTMU[0]) {
        match hv.try_into() {
            Ok(m) => m,
            Err(_) => return false,
        }
    } else {
        let hvk = hashed_vrfykey_from_vrfykey(vrfy_key);
        match compute_mu(&hvk, ctx, id, hv) {
            Some(m) => m,
            None => return false,
        }
    };

    let mut buf = [0u16; 2 * MAX_N];
    let (t1, t2) = buf.split_at_mut(n);
    let t2 = &mut t2[..n];

    // t1 <- h (verifying key, decoded); h is in ntt representation.
    if mqpoly_decode(logn, &vrfy_key[1..], t1) != vrfy_key.len() - 1 {
        return false;
    }
    (k.ext_to_int)(logn, t1);

    // t2 <- s2 (signature, decoded, converted to ntt); also get its norm.
    // comp_decode() checks the L-infinity norm of s2.
    if !comp_decode(logn, &sig[41..], t2) {
        return false;
    }
    let norm2 = (k.sqnorm_signed)(logn, t2);
    (k.signed_to_int)(logn, t2);
    (k.int_to_ntt)(logn, t2);

    // t2 <- s2*h (converted to int)
    (k.mul_ntt)(logn, t2, t1);
    (k.ntt_to_int)(logn, t2);

    // Hash message into polynomial c (into t1, converted to int).
    let nonce: &[u8; 40] = sig[1..41].try_into().unwrap();
    hash_to_point(logn, nonce, &mu, t1);
    (k.ext_to_int)(logn, t1);

    // t1 <- s1 = c - s2*h (converted to ext), and compute its norm; this
    // also verifies the L-infinity norm of s1.
    (k.sub)(logn, t1, t2);
    (k.int_to_ext)(logn, t1);
    let norm1 = (k.sqnorm_binf_ext)(logn, t1);

    // Valid if the total squared norm of (s1, s2) is small enough. Beware
    // of overflows.
    if norm2 != 0 && norm1 >= norm2.wrapping_neg() {
        return false;
    }
    mq::sqnorm_is_acceptable(logn, norm1.wrapping_add(norm2))
}

/// Verify a signature (standard degrees 512 and 1024).
///
/// `id` is the hash identifier (empty for a raw message, a DER-encoded OID
/// for a pre-hashed message, or [`HASH_ID_EXTMU`] with `hv` being the 64-byte
/// message representative; `ctx` is then ignored).
pub fn verify(sig: &[u8], vrfy_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8]) -> bool {
    inner_verify(mq::best(), 9, 10, sig, vrfy_key, ctx, id, hv)
}

/// Verify a signature for the weak degrees (4 to 256).
pub fn verify_weak(sig: &[u8], vrfy_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8]) -> bool {
    inner_verify(mq::best(), 2, 8, sig, vrfy_key, ctx, id, hv)
}

#[doc(hidden)]
pub fn verify_with(
    k: &Kernels,
    weak: bool,
    sig: &[u8],
    vrfy_key: &[u8],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
) -> bool {
    if weak {
        inner_verify(k, 2, 8, sig, vrfy_key, ctx, id, hv)
    } else {
        inner_verify(k, 9, 10, sig, vrfy_key, ctx, id, hv)
    }
}
