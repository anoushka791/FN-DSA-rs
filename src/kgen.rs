//! Key pair generation (port of kgen.c and kgen_gauss.c).

use crate::codec::{mqpoly_encode, trim_i8_encode};
use crate::mq;
use crate::mqsign as ms;
use crate::ntru::{check_ortho_norm, solve_ntru};
use crate::sha3::Shake;
use crate::util::sysrng;
use crate::{fg_nbits, sign_key_size, vrfy_key_size};

// ---------------------------------------------------------------------
// Gaussian generation of (f, g).

/// q = 12289, n = 512
static KGDIST_512: [u16; 17] =
    [29543, 23286, 17574, 12669, 8706, 5692, 3535, 2083, 1164, 615, 308, 146, 65, 28, 11, 4, 1];

/// q = 12289, n = 1024
static KGDIST_1024: [u16; 12] = [28207, 19623, 12472, 7198, 3753, 1761, 742, 280, 94, 28, 8, 2];

/// Sample f (or g) from the SHAKE-based PRNG; the polynomial has odd parity.
fn sample_f(logn: u32, pc: &mut Shake, f: &mut [i8]) {
    let n = 1usize << logn;

    // Lookup tables exist for degrees 512 and 1024; for smaller degrees,
    // the degree-512 table is used, adding 512/n sampled values together.
    let (tab, zz): (&[u16], u32) = if logn == 10 {
        (&KGDIST_1024, 1)
    } else {
        (&KGDIST_512, 1u32 << (9 - logn))
    };

    // Each value is sampled once, except the last one, which is sampled
    // repeatedly until the complete polynomial has odd parity.
    let mut parity = 0u32;
    let mut i = 0usize;
    while i < n {
        // Sample the value as an unsigned 32-bit word.
        let mut aa = 0u32;
        for _ in 0..zz {
            let u = pc.next_u16();
            let v = u >> 1;
            // Add 2^16 to a for each table element greater than v.
            let mut a = 0u32;
            for &t in tab.iter() {
                a = a.wrapping_sub(v.wrapping_sub(t as u32) & !0xFFFFu32);
            }
            // Apply the sign bit (lsb of u).
            a = a.wrapping_sub((a << 1) & (u & 1).wrapping_neg());
            aa = aa.wrapping_add(a);
        }

        // The sampled value is in the top 16 bits of aa.
        let ai = (aa as i32) >> 16;

        // At low (test) degrees the value may be outside [-127, +127],
        // which does not fit the storage format: resample.
        if !(-127..=127).contains(&ai) {
            continue;
        }
        f[i] = ai as i8;

        // For the last element, insist on an odd parity.
        if i == n - 1 {
            if ((parity ^ ai as u32) & 1) != 1 {
                continue;
            }
        } else {
            parity ^= ai as u32;
        }
        i += 1;
    }
}

/// Check whether f is invertible modulo q (and X^n+1). `tmp` needs n values.
fn is_invertible(logn: u32, f: &[i8], tmp: &mut [u16]) -> bool {
    let n = 1usize << logn;
    ms::small_to_int(logn, f, tmp);
    (mq::best().int_to_ntt)(logn, tmp);
    let mut r = 0xFFFF_FFFFu32;
    for &v in tmp[..n].iter() {
        r &= (v as u32).wrapping_sub(mq::Q);
    }
    (r >> 16) != 0
}

/// Core key generation. `tmp` is a 32-byte aligned area of 22*n bytes.
///
/// # Safety
/// `tmp` must point to `22 << logn` valid bytes, 32-byte aligned.
unsafe fn keygen_inner(
    logn: u32,
    seed: &[u8],
    mut sign_key: Option<&mut [u8]>,
    mut vrfy_key: Option<&mut [u8]>,
    tmp: *mut u8,
) {
    let n = 1usize << logn;
    // f and g at the end of the temporary area.
    let f: &mut [i8] = core::slice::from_raw_parts_mut(tmp.add(20 * n) as *mut i8, n);
    let g: &mut [i8] = core::slice::from_raw_parts_mut(tmp.add(21 * n) as *mut i8, n);
    unsafe fn t16<'a>(tmp: *mut u8, off: usize, len: usize) -> &'a mut [u16] {
        core::slice::from_raw_parts_mut(tmp.add(off) as *mut u16, len)
    }

    // PRNG with the provided seed.
    let mut pc = Shake::new(256);
    pc.inject(seed);
    pc.flip();
    loop {
        // Sample f, with odd parity; retry if not invertible mod q. (The
        // NTT of f is not kept, to save RAM; it is recomputed at the end.)
        sample_f(logn, &mut pc, f);
        if !is_invertible(logn, f, t16(tmp, 0, n)) {
            continue;
        }

        // Sample g, also with odd parity.
        sample_f(logn, &mut pc, g);

        // Ensure ||(g, -f)||^2 < (1.17^2)*q = 16822.4121.
        let mut sn = 0i32;
        for i in 0..n {
            let xf = f[i] as i32;
            let xg = g[i] as i32;
            sn += xf * xf + xg * xg;
        }
        if sn >= 16823 {
            continue;
        }

        // (f, g) must have an acceptable orthogonalized norm.
        if !check_ortho_norm(logn, f, g, tmp as *mut u64) {
            continue;
        }

        // Try to solve the NTRU equation.
        if !solve_ntru(logn, f, g, tmp as *mut u32) {
            continue;
        }

        // F and G are at the start of tmp. Encode the private key (header,
        // f, g, F); the hash of the public key comes after.
        let mut j = 1usize;
        if let Some(sk) = sign_key.as_deref_mut() {
            let big_f: &[i8] = core::slice::from_raw_parts(tmp as *const i8, n);
            sk[0] = 0x50 + logn as u8;
            let nbits = fg_nbits(logn);
            j += trim_i8_encode(logn, f, nbits, &mut sk[j..]);
            j += trim_i8_encode(logn, g, nbits, &mut sk[j..]);
            j += trim_i8_encode(logn, big_f, 8, &mut sk[j..]);
        }

        if sign_key.is_some() || vrfy_key.is_some() {
            // The signing key includes the hash of the verifying key.
            // F, G are already encoded, so they can be overwritten.
            let h = t16(tmp, 0, n);
            let ft = t16(tmp, 2 * n, n);
            let scratch = t16(tmp, 4 * n, n);
            ms::small_to_int(logn, f, ft);
            ms::small_to_int(logn, g, h);
            let k = mq::best();
            (k.int_to_ntt)(logn, ft);
            (k.int_to_ntt)(logn, h);
            ms::div_ntt(logn, h, ft, scratch);
            (k.int_to_ext)(logn, h);

            // Encode the public key in the temporary area (over ft).
            let buf = core::slice::from_raw_parts_mut(tmp.add(2 * n) as *mut u8, vrfy_key_size(logn));
            buf[0] = logn as u8;
            let vklen = 1 + mqpoly_encode(logn, t16(tmp, 0, n), &mut buf[1..]);

            if let Some(sk) = sign_key.as_deref_mut() {
                let mut sc = Shake::new(256);
                sc.inject(&buf[..vklen]);
                sc.flip();
                sc.extract(&mut sk[j..j + 64]);
            }
            if let Some(vk) = vrfy_key.as_deref_mut() {
                vk[..vklen].copy_from_slice(&buf[..vklen]);
            }
        }
        break;
    }
}

/// Temporary area, aligned on a 32-byte boundary.
#[repr(C, align(32))]
struct Tmp<const N: usize>([u8; N]);

// Wrappers allocating the temporary area on the stack, with several sizes so
// that stack use is not always the worst case.
macro_rules! keygen_wrap {
    ($name:ident, $sz:expr) => {
        #[inline(never)]
        fn $name(logn: u32, seed: &[u8], sign_key: Option<&mut [u8]>, vrfy_key: Option<&mut [u8]>) {
            let mut tmp = Tmp::<{ $sz * 22 }>([0u8; $sz * 22]);
            // SAFETY: the area is 32-byte aligned and 22*n bytes for n <= sz.
            unsafe { keygen_inner(logn, seed, sign_key, vrfy_key, tmp.0.as_mut_ptr()) }
        }
    };
}

keygen_wrap!(keygen_32, 32);
keygen_wrap!(keygen_64, 64);
keygen_wrap!(keygen_128, 128);
keygen_wrap!(keygen_256, 256);
keygen_wrap!(keygen_512, 512);
keygen_wrap!(keygen_1024, 1024);

fn keygen_common(
    logn: u32,
    seed: Option<&[u8]>,
    sign_key: Option<&mut [u8]>,
    vrfy_key: Option<&mut [u8]>,
) -> bool {
    if !(2..=10).contains(&logn) {
        return false;
    }
    // Output buffers must be large enough.
    if sign_key.as_ref().map_or(false, |k| k.len() < sign_key_size(logn))
        || vrfy_key.as_ref().map_or(false, |k| k.len() < vrfy_key_size(logn))
    {
        return false;
    }

    // Without a seed, use a 32-byte seed from the system RNG.
    let mut seedbuf = [0u8; 32];
    let seed: &[u8] = match seed {
        Some(s) => s,
        None => {
            if !sysrng(&mut seedbuf) {
                return false;
            }
            &seedbuf
        }
    };

    match logn {
        6 => keygen_64(logn, seed, sign_key, vrfy_key),
        7 => keygen_128(logn, seed, sign_key, vrfy_key),
        8 => keygen_256(logn, seed, sign_key, vrfy_key),
        9 => keygen_512(logn, seed, sign_key, vrfy_key),
        10 => keygen_1024(logn, seed, sign_key, vrfy_key),
        _ => keygen_32(logn, seed, sign_key, vrfy_key),
    }
    true
}

/// Generate a key pair using the system RNG. Either output may be `None`.
/// `logn` is 9 or 10 for the standard degrees (2 to 8 are "weak" test
/// degrees). Returns false on invalid parameters, short output buffers, or
/// RNG failure.
pub fn keygen(logn: u32, sign_key: Option<&mut [u8]>, vrfy_key: Option<&mut [u8]>) -> bool {
    keygen_common(logn, None, sign_key, vrfy_key)
}

/// Generate a key pair deterministically from `seed`.
pub fn keygen_seeded(
    logn: u32,
    seed: &[u8],
    sign_key: Option<&mut [u8]>,
    vrfy_key: Option<&mut [u8]>,
) -> bool {
    keygen_common(logn, Some(seed), sign_key, vrfy_key)
}
