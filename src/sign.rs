//! Top-level signature generation.

use crate::sha3::Shake;
use crate::sign_core::sign_core;
use crate::util::{compute_mu, sysrng};
use crate::{fg_nbits, sign_key_size, signature_size, HASH_ID_EXTMU};

const SEEDBUF_LEN: usize = 40;

/// Make a 40-byte seed from the provided seed (or from the system RNG if
/// `seed` is `None`). The private key (f||g) and the message representative
/// are mixed in ("hedging"), so that a poor system RNG is tolerated.
#[inline(never)]
fn make_seedbuf(sign_key: &[u8], mu: &[u8; 64], seed: Option<&[u8]>) -> Option<[u8; SEEDBUF_LEN]> {
    let logn = (sign_key[0] & 0x0F) as u32;
    let nbits = fg_nbits(logn) as usize;

    // Hash of the private key (f||g).
    let mut sc = Shake::new(256);
    sc.inject(&sign_key[1..1 + (nbits << (logn - 2))]);
    sc.flip();
    let mut seedbuf = [0u8; SEEDBUF_LEN];
    sc.extract(&mut seedbuf);

    // Derived seed:
    //   SHAKE256(impl_id || SHAKE256(f||g)[40] || mu || seed)[40]
    //
    // `impl_id` characterizes the implementation variant: two implementations
    // that do not perform exactly the same floating-point operations in the
    // same order may diverge on the same key and seed, and that divergence
    // could leak information about the private key. This port performs
    // exactly the same operations as the C code (c-fn-dsa-alt), hence the
    // same identifier. If this code is ever changed so that the sequence of
    // floating-point operations changes, this identifier MUST be changed.
    let mut sc = Shake::new(256);
    sc.inject(b"20260830");
    sc.inject(&seedbuf);
    sc.inject(mu);
    match seed {
        Some(s) => sc.inject(s),
        None => {
            let mut rnd = [0u8; SEEDBUF_LEN];
            if !sysrng(&mut rnd) {
                return None;
            }
            sc.inject(&rnd);
        }
    }
    sc.flip();
    sc.extract(&mut seedbuf);
    Some(seedbuf)
}

/// Temporary area, aligned on a 32-byte boundary.
#[repr(C, align(32))]
struct Tmp<const N: usize>([u8; N]);

fn sign_step1(
    logn: u32,
    sign_key: &[u8],
    mu: &[u8; 64],
    seed: Option<&[u8]>,
    sig: &mut [u8],
    tmp: *mut u8,
) -> usize {
    let seedbuf = match make_seedbuf(sign_key, mu, seed) {
        Some(s) => s,
        None => return 0,
    };
    // SAFETY: `tmp` points to a 32-byte aligned area of 20*n bytes (see the
    // wrappers below); the caller checked the sizes of the key and of `sig`.
    unsafe { sign_core(logn, &sign_key[1..], mu, &seedbuf, sig, tmp) }
}

// Wrappers allocating the temporary area on the stack; several sizes are
// defined so that the stack use is not always the worst case.
macro_rules! sign_wrap {
    ($name:ident, $sz:expr) => {
        #[inline(never)]
        fn $name(
            logn: u32,
            sign_key: &[u8],
            mu: &[u8; 64],
            seed: Option<&[u8]>,
            sig: &mut [u8],
        ) -> usize {
            let mut tmp = Tmp::<{ $sz * 20 }>([0u8; $sz * 20]);
            sign_step1(logn, sign_key, mu, seed, sig, tmp.0.as_mut_ptr())
        }
    };
}

sign_wrap!(sign_32, 32);
sign_wrap!(sign_64, 64);
sign_wrap!(sign_128, 128);
sign_wrap!(sign_256, 256);
sign_wrap!(sign_512, 512);
sign_wrap!(sign_1024, 1024);

fn sign_wrapper(
    weak: bool,
    sign_key: &[u8],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
    seed: Option<&[u8]>,
    sig: &mut [u8],
) -> Option<usize> {
    if sign_key.is_empty() {
        return None;
    }
    let head = sign_key[0] as u32;
    if head & 0xF0 != 0x50 {
        return None;
    }
    let logn = head & 0x0F;
    let range_ok = if weak { (2..=8).contains(&logn) } else { (9..=10).contains(&logn) };
    if !range_ok || sign_key.len() != sign_key_size(logn) {
        return None;
    }
    let sig_len = signature_size(logn);
    if sig.len() < sig_len {
        return None;
    }

    // Message representative mu.
    let mu: [u8; 64] = if id.first() == Some(&HASH_ID_EXTMU[0]) {
        hv.try_into().ok()?
    } else {
        let hvk: &[u8; 64] = sign_key[sign_key.len() - 64..].try_into().unwrap();
        compute_mu(hvk, ctx, id, hv)?
    };

    let sig = &mut sig[..sig_len];
    let r = match logn {
        6 => sign_64(logn, sign_key, &mu, seed, sig),
        7 => sign_128(logn, sign_key, &mu, seed, sig),
        8 => sign_256(logn, sign_key, &mu, seed, sig),
        9 => sign_512(logn, sign_key, &mu, seed, sig),
        10 => sign_1024(logn, sign_key, &mu, seed, sig),
        _ => sign_32(logn, sign_key, &mu, seed, sig),
    };
    if r == 0 {
        None
    } else {
        Some(r)
    }
}

/// Sign (standard degrees 512 and 1024), using the system RNG. Returns the
/// signature length on success.
pub fn sign(sign_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8], sig: &mut [u8]) -> Option<usize> {
    sign_wrapper(false, sign_key, ctx, id, hv, None, sig)
}

/// Sign (standard degrees), with all randomness derived from `seed`. The
/// seed should have at least 320 bits of entropy (40 bytes).
pub fn sign_seeded(
    sign_key: &[u8],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
    seed: &[u8],
    sig: &mut [u8],
) -> Option<usize> {
    sign_wrapper(false, sign_key, ctx, id, hv, Some(seed), sig)
}

/// Sign with a weak-degree key (4 to 256); for tests and research only.
pub fn sign_weak(sign_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8], sig: &mut [u8]) -> Option<usize> {
    sign_wrapper(true, sign_key, ctx, id, hv, None, sig)
}

/// Seeded variant of [`sign_weak`].
pub fn sign_weak_seeded(
    sign_key: &[u8],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
    seed: &[u8],
    sig: &mut [u8],
) -> Option<usize> {
    sign_wrapper(true, sign_key, ctx, id, hv, Some(seed), sig)
}
