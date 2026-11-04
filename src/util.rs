//! Hashing helpers: hash-to-point, hashed verifying key, message
//! representative (mu).

use crate::sha3::Shake;

/// Hash a (pre-hashed) message into a polynomial with coefficients in
/// [0, q-1] (external representation).
pub fn hash_to_point(logn: u32, nonce: &[u8; 40], mu: &[u8; 64], c: &mut [u16]) {
    let mut sc = Shake::new(256);
    sc.inject(nonce);
    sc.inject(mu);
    sc.flip();

    let n = 1usize << logn;
    let c = &mut c[..n];
    let mut sbuf = [0u8; 136];
    let mut j = sbuf.len();
    let mut i = 0usize;
    while i < n {
        if j == 136 {
            sc.extract(&mut sbuf);
            j = 0;
        }
        let mut w = sbuf[j] as u32 | ((sbuf[j + 1] as u32) << 8);
        j += 2;
        if w < 61445 {
            while w >= 12289 {
                w -= 12289;
            }
            c[i] = w as u16;
            i += 1;
        }
    }
}

/// Compute the 64-byte hash of an encoded verifying key.
pub fn hashed_vrfykey_from_vrfykey(vrfy_key: &[u8]) -> [u8; 64] {
    let mut sc = Shake::new(256);
    sc.inject(vrfy_key);
    sc.flip();
    let mut out = [0u8; 64];
    sc.extract(&mut out);
    out
}

/// Compute the message representative mu from the hashed verifying key, the
/// context, the hash identifier and the (pre-hashed) message. Returns `None`
/// if the parameters are invalid (context longer than 255 bytes, unknown or
/// malformed hash identifier, including the "external mu" identifier).
///
/// An empty `id` means a raw message.
pub fn compute_mu(
    hashed_vkey: &[u8; 64],
    ctx: &[u8],
    id: &[u8],
    hv: &[u8],
) -> Option<[u8; 64]> {
    if ctx.len() > 255 {
        return None;
    }
    let mut sc = Shake::new(256);
    sc.inject(hashed_vkey);

    // An empty identifier is the C empty string, i.e. first byte 0x00.
    let first = id.first().copied().unwrap_or(0x00);
    let id_part: &[u8] = match first {
        0x00 => {
            sc.inject(&[0]);
            &[]
        }
        0x06 => {
            // Tag, length (at most 127), then the OID bytes.
            let len = *id.get(1)? as usize;
            if len > 127 || id.len() < len + 2 {
                return None;
            }
            sc.inject(&[1]);
            &id[..len + 2]
        }
        _ => return None,
    };

    sc.inject(&[ctx.len() as u8]);
    sc.inject(ctx);
    sc.inject(id_part);
    sc.inject(hv);
    sc.flip();
    let mut mu = [0u8; 64];
    sc.extract(&mut mu);
    Some(mu)
}

/// Obtain fresh randomness from the operating system. Returns false on
/// failure (or if no secure source is available on this platform).
pub fn sysrng(dst: &mut [u8]) -> bool {
    if dst.is_empty() {
        return true;
    }
    sysrng_impl(dst)
}

#[cfg(windows)]
fn sysrng_impl(dst: &mut [u8]) -> bool {
    #[link(name = "advapi32")]
    extern "system" {
        // RtlGenRandom
        fn SystemFunction036(buf: *mut u8, len: u32) -> u8;
    }
    for chunk in dst.chunks_mut(1 << 20) {
        // SAFETY: `chunk` is a valid writable buffer of the given length.
        if unsafe { SystemFunction036(chunk.as_mut_ptr(), chunk.len() as u32) } == 0 {
            return false;
        }
    }
    true
}

#[cfg(all(unix, not(windows)))]
fn sysrng_impl(dst: &mut [u8]) -> bool {
    use std::io::Read;
    match std::fs::File::open("/dev/urandom") {
        Ok(mut f) => f.read_exact(dst).is_ok(),
        Err(_) => false,
    }
}

#[cfg(not(any(windows, unix)))]
fn sysrng_impl(_dst: &mut [u8]) -> bool {
    false
}
