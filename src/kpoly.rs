//! Polynomials with integer coefficients (key generation).
//!
//! Polynomials use an interleaved in-memory representation: n = 2^logn
//! coefficients of `len` words each, with stride n; the n coefficients start
//! on consecutive words. A polynomial may be in plain representation (signed
//! or unsigned) or in RNS (the words are the integer reduced modulo the first
//! `len` primes), possibly in NTT form. This code works on raw pointers since
//! the access patterns are strided and the buffers are carved out of one
//! shared temporary area.

use crate::fxp::Fxr;
use crate::mp31::*;
use crate::zint31::*;

/// q = x / 31 and r = x % 31 for x up to 63487 (constant-time).
#[inline(always)]
pub fn divrem31(x: u32) -> (u32, u32) {
    let q = x.wrapping_mul(67651) >> 21;
    (q, x - 31 * q)
}

pub fn poly_mp_set_small(logn: u32, d: &mut [u32], f: &[i8], p: u32) {
    let n = 1usize << logn;
    for i in 0..n {
        d[i] = mp_set(f[i] as i32, p);
    }
}

/// Convert one-word signed normal representation to RNS modulo p.
pub fn poly_mp_set(logn: u32, f: &mut [u32], p: u32) {
    let n = 1usize << logn;
    for v in f[..n].iter_mut() {
        let mut x = *v;
        x |= (x & 0x4000_0000) << 1;
        *v = mp_set(x as i32, p);
    }
}

/// Convert RNS modulo p to one-word signed normal representation.
pub fn poly_mp_norm(logn: u32, f: &mut [u32], p: u32) {
    let n = 1usize << logn;
    for v in f[..n].iter_mut() {
        *v = (mp_norm(*v, p) as u32) & 0x7FFF_FFFF;
    }
}

/// Convert to small integers; returns false if any coefficient exceeds
/// `lim` in absolute value.
pub fn poly_big_to_small(logn: u32, d: &mut [i8], s: &[u32], lim: i32) -> bool {
    let n = 1usize << logn;
    for i in 0..n {
        let mut x = s[i];
        x |= (x & 0x4000_0000) << 1;
        let z = x as i32;
        if z < -lim || z > lim {
            return false;
        }
        d[i] = z as i8;
    }
    true
}

/// Maximum bit length of all coefficients of a polynomial (constant-time).
///
/// # Safety
/// `f` must be valid for `flen << logn` reads.
pub unsafe fn poly_max_bitlength(logn: u32, f: *const u32, flen: usize) -> u32 {
    if flen == 0 {
        return 0;
    }
    let n = 1usize << logn;
    let mut t = 0u32;
    let mut tk = 0u32;
    for i in 0..n {
        let fp = f.add(i);
        // Extend the sign bit into a 31-bit mask.
        let m = (*fp.add((flen - 1) << logn) >> 30).wrapping_neg() & 0x7FFF_FFFF;

        // Top non-zero sign-adjusted word, with index.
        let mut c = 0u32;
        let mut ck = 0u32;
        for j in 0..flen {
            let w = *fp.add(j << logn) ^ m;
            let nz = (w.wrapping_sub(1) >> 31).wrapping_sub(1);
            c ^= nz & (c ^ w);
            ck ^= nz & (ck ^ j as u32);
        }

        // Replace (t, tk) with (c, ck) if ck > tk, or ck == tk and c > t.
        let rr = tbmask(tk.wrapping_sub(ck) | (((tk ^ ck).wrapping_sub(1)) & t.wrapping_sub(c)));
        t ^= rr & (t ^ c);
        tk ^= rr & (tk ^ ck);
    }
    31 * tk + 32 - t.leading_zeros()
}

/// Convert a polynomial to fixed-point approximations x/2^sc (constant-time).
///
/// # Safety
/// `f` must be valid for `len << logn` reads.
pub unsafe fn poly_big_to_fixed(logn: u32, d: &mut [Fxr], f: *const u32, len: usize, sc: u32) {
    let n = 1usize << logn;
    if len == 0 {
        d[..n].fill(0);
        return;
    }

    // sc = 31*sch + scl, with scl in 1..31 (sch may wrap around to -1).
    let (mut sch, mut scl) = divrem31(sc);
    let z = scl.wrapping_sub(1) >> 31;
    sch = sch.wrapping_sub(z);
    scl |= 31 & z.wrapping_neg();

    let t0 = sch.wrapping_sub(1) & 0xFF_FFFF;
    let t1 = sch & 0xFF_FFFF;
    let t2 = sch.wrapping_add(1) & 0xFF_FFFF;

    for i in 0..n {
        let fp = f.add(i);
        let mut w0 = 0u32;
        let mut w1 = 0u32;
        let mut w2 = 0u32;
        for j in 0..len {
            let w = *fp.add(j << logn);
            let t = (j as u32) & 0xFF_FFFF;
            w0 |= w & (((t ^ t0).wrapping_sub(1)) >> 31).wrapping_neg();
            w1 |= w & (((t ^ t1).wrapping_sub(1)) >> 31).wrapping_neg();
            w2 |= w & (((t ^ t2).wrapping_sub(1)) >> 31).wrapping_neg();
        }

        // Not enough words for the requested scaling: supply copies with
        // the proper sign.
        let ws = (*fp.add((len - 1) << logn) >> 30).wrapping_neg() >> 1;
        let l = len as u32;
        w0 |= ws & (l.wrapping_sub(sch) >> 31).wrapping_neg();
        w1 |= ws & (l.wrapping_sub(sch).wrapping_sub(1) >> 31).wrapping_neg();
        w2 |= ws & (l.wrapping_sub(sch).wrapping_sub(2) >> 31).wrapping_neg();

        // Assemble the 64-bit value (scl is in 1..31); sign-extend w2.
        w2 |= (w2 & 0x4000_0000) << 1;
        let xl = (w0 >> (scl - 1)) | (w1 << (32 - scl));
        let xh = (w1 >> scl) | (w2 << (31 - scl));
        d[i] = (xl as u64) | ((xh as u64) << 32);
    }
}

/// F <- F - (2^sc)*k*f, with F, f, k polynomials modulo X^n+1; quadratic
/// algorithm, no extra buffer.
///
/// # Safety
/// F valid for `Flen << logn` words, f for `flen << logn`, k for `n`.
pub unsafe fn poly_sub_scaled(
    logn: u32,
    f_big: *mut u32,
    flen_big: usize,
    f: *const u32,
    flen: usize,
    k: *const i32,
    sc: u32,
) {
    if flen == 0 {
        return;
    }
    let (sch, scl) = divrem31(sc);
    if sch as usize >= flen_big {
        return;
    }
    let big_f = f_big.add((sch as usize) << logn);
    let big_flen = flen_big - sch as usize;
    let n = 1usize << logn;

    if logn <= 3 {
        // Small degrees: fully unrolled in the C code; same arithmetic here.
        let mut t = [0u32; 8];
        let mut signf = [0u32; 8];
        let mut kk = [0i32; 8];
        let mut cc = [0i64; 8];
        for j in 0..n {
            signf[j] = (*f.add((flen << logn) - n + j) >> 30).wrapping_neg() >> 1;
            kk[j] = *k.add(j);
        }
        for i in 0..big_flen {
            let mut fs = [0u32; 8];
            for j in 0..n {
                let fw = if i < flen { *f.add((i << logn) + j) } else { signf[j] };
                fs[j] = ((fw << scl) & 0x7FFF_FFFF) | t[j];
                t[j] = fw >> (31 - scl);
            }
            let mut z = [0i64; 8];
            for j in 0..n {
                let mut zj = (*big_f.add((i << logn) + j) as i64).wrapping_add(cc[j]);
                for l in 0..n {
                    if l <= j {
                        zj = zj.wrapping_sub((fs[l] as i64).wrapping_mul(kk[j - l] as i64));
                    } else {
                        zj = zj.wrapping_add((fs[l] as i64).wrapping_mul(kk[n + j - l] as i64));
                    }
                }
                z[j] = zj;
            }
            for j in 0..n {
                *big_f.add((i << logn) + j) = (z[j] as u32) & 0x7FFF_FFFF;
                cc[j] = z[j] >> 31;
            }
        }
        return;
    }

    for i in 0..n {
        let mut kf = (*k.add(i)).wrapping_neg();
        let mut x = big_f.add(i);
        for j in 0..n {
            zint_add_scaled_mul_small(x, big_flen, f.add(j), flen, n, kf, 0, scl);
            if i + j == n - 1 {
                x = big_f;
                kf = kf.wrapping_neg();
            } else {
                x = x.add(1);
            }
        }
    }
}

/// F <- F - (2^sc)*k*f using the NTT. f MUST be in RNS+NTT over flen+1
/// words; F MUST be in plain representation.
///
/// # Safety
/// Pointers valid for the layouts documented in the C code; `tmp` provides
/// room for `2n + (flen+1)*n + n` words.
pub unsafe fn poly_sub_scaled_ntt(
    logn: u32,
    big_f: *mut u32,
    flen_big: usize,
    f: *const u32,
    flen: usize,
    k: *const i32,
    sc: u32,
    tmp: *mut u32,
) {
    let n = 1usize << logn;
    let tlen = flen + 1;
    let gm = tmp;
    let igm = gm.add(n);
    let fk = igm.add(n);
    let t1 = fk.add(tlen << logn);
    let (sch, scl) = divrem31(sc);

    // k*f in fk[], in RNS notation.
    for i in 0..tlen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        mp_mkgmigm(
            logn,
            core::slice::from_raw_parts_mut(gm, n),
            core::slice::from_raw_parts_mut(igm, n),
            PRIMES[i].g,
            PRIMES[i].ig,
            p,
            p0i,
        );
        for j in 0..n {
            *t1.add(j) = mp_set(*k.add(j), p);
        }
        mp_ntt(logn, core::slice::from_raw_parts_mut(t1, n), core::slice::from_raw_parts(gm, n), p, p0i);

        let fs = f.add(i << logn);
        let ff = fk.add(i << logn);
        for j in 0..n {
            *ff.add(j) = mp_mmul(mp_mmul(*t1.add(j), *fs.add(j), p, p0i), r2, p, p0i);
        }
        mp_intt(logn, core::slice::from_raw_parts_mut(ff, n), core::slice::from_raw_parts(igm, n), p, p0i);
    }

    // Rebuild k*f.
    zint_rebuild_crt(fk, tlen, n, 1, true, t1);

    // Subtract k*f, scaled, from F.
    for i in 0..n {
        zint_sub_scaled(big_f.add(i), flen_big, fk.add(i), tlen, n, sch, scl);
    }
}

/// Depth-1 specialization: F <- F - (2^sc)*k*ft, with ft = N(f). `k` is
/// consumed. See the C documentation. Uses 3*n words in `tmp`.
///
/// # Safety
/// Pointers valid for the layouts documented in the C code.
pub unsafe fn poly_sub_kf_scaled_depth1(
    logn_top: u32,
    big_f: *mut u32,
    fglen: usize,
    k: *mut u32,
    sc: u32,
    f: &[i8],
    tmp: *mut u32,
) {
    let logn = logn_top - 1;
    let n = 1usize << logn;
    let hn = n >> 1;
    let gm = tmp;
    let t1 = gm.add(n);
    let t2 = t1.add(n);
    let sl = |p: *mut u32| core::slice::from_raw_parts_mut(p, n);

    // Step 1: convert F to RNS (FGlen is 1 or 2).
    if fglen == 1 {
        let p = PRIMES[0].p;
        for i in 0..n {
            let mut xf = *big_f.add(i);
            xf |= (xf & 0x4000_0000) << 1;
            *big_f.add(i) = mp_set(xf as i32, p);
        }
    } else {
        let p0 = PRIMES[0].p;
        let p0_0i = PRIMES[0].p0i;
        let z0 = mp_half(PRIMES[0].r2, p0);
        let p1 = PRIMES[1].p;
        let p1_0i = PRIMES[1].p0i;
        let z1 = mp_half(PRIMES[1].r2, p1);
        for i in 0..n {
            let xl = *big_f.add(i);
            let fh = *big_f.add(i + n);
            let xh = fh | ((fh & 0x4000_0000) << 1);
            let yl0 = xl.wrapping_sub(p0 & !tbmask(xl.wrapping_sub(p0)));
            let yh0 = mp_set(xh as i32, p0);
            let r0 = mp_add(yl0, mp_mmul(yh0, z0, p0, p0_0i), p0);
            let yl1 = xl.wrapping_sub(p1 & !tbmask(xl.wrapping_sub(p1)));
            let yh1 = mp_set(xh as i32, p1);
            let r1 = mp_add(yl1, mp_mmul(yh1, z1, p1, p1_0i), p1);
            *big_f.add(i) = r0;
            *big_f.add(i + n) = r1;
        }
    }

    // Step 2: for FGlen small primes, convert F to RNS+NTT, and subtract
    // (2^sc)*ft, with ft computed in RNS+NTT dynamically.
    for i in 0..fglen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        let r3 = mp_mmul(r2, r2, p, p0i);
        mp_mkgm(logn, sl(gm), PRIMES[i].g, p, p0i);

        // k <- (2^sc)*k (and into NTT).
        let mut scv = mp_mmul(1u32 << (sc & 31), r2, p, p0i);
        let mut m = sc >> 5;
        while m > 0 {
            scv = mp_mmul(scv, r2, p, p0i);
            m -= 1;
        }
        for j in 0..n {
            let x = mp_set(*k.add(j) as i32, p);
            *k.add(j) = mp_mmul(scv, x, p, p0i);
        }
        mp_ntt(logn, sl(k), sl(gm), p, p0i);

        // Convert F to NTT.
        let fu = big_f.add(i << logn);
        mp_ntt(logn, sl(fu), sl(gm), p, p0i);

        // ft = N(f) = f_e^2 - X*f_o^2 with f = f_e(X^2) + X*f_o(X^2).
        for j in 0..n {
            *t1.add(j) = mp_set(f[(j << 1) + 0] as i32, p);
            *t2.add(j) = mp_set(f[(j << 1) + 1] as i32, p);
        }
        mp_ntt(logn, sl(t1), sl(gm), p, p0i);
        mp_ntt(logn, sl(t2), sl(gm), p, p0i);
        for j in 0..hn {
            let mut xe0 = *t1.add((j << 1) + 0);
            let mut xe1 = *t1.add((j << 1) + 1);
            let mut xo0 = *t2.add((j << 1) + 0);
            let mut xo1 = *t2.add((j << 1) + 1);
            let xv0 = *gm.add(hn + j);
            let xv1 = p.wrapping_sub(xv0);
            xe0 = mp_mmul(xe0, xe0, p, p0i);
            xe1 = mp_mmul(xe1, xe1, p, p0i);
            xo0 = mp_mmul(xo0, xo0, p, p0i);
            xo1 = mp_mmul(xo1, xo1, p, p0i);
            let xf0 = mp_sub(xe0, mp_mmul(xo0, xv0, p, p0i), p);
            let xf1 = mp_sub(xe1, mp_mmul(xo1, xv1, p, p0i), p);

            let xkf0 = mp_mmul(mp_mmul(xf0, *k.add((j << 1) + 0), p, p0i), r3, p, p0i);
            let xkf1 = mp_mmul(mp_mmul(xf1, *k.add((j << 1) + 1), p, p0i), r3, p, p0i);
            *fu.add((j << 1) + 0) = mp_sub(*fu.add((j << 1) + 0), xkf0, p);
            *fu.add((j << 1) + 1) = mp_sub(*fu.add((j << 1) + 1), xkf1, p);
        }

        // Convert F back to RNS.
        mp_mkigm(logn, sl(t1), PRIMES[i].ig, p, p0i);
        mp_intt(logn, sl(fu), sl(t1), p, p0i);

        // Put k back to its initial value if there is another iteration.
        if i + 1 < fglen {
            mp_intt(logn, sl(k), sl(t1), p, p0i);
            let mut scv = 1u32 << (sc.wrapping_neg() & 31);
            let mut m = sc >> 5;
            while m > 0 {
                scv = mp_mmul(scv, 1, p, p0i);
                m -= 1;
            }
            for j in 0..n {
                *k.add(j) = mp_norm(mp_mmul(scv, *k.add(j), p, p0i), p) as u32;
            }
        }
    }

    // Output F is in RNS (non-NTT); convert back to plain integers.
    if fglen == 1 {
        let p = PRIMES[0].p;
        for i in 0..n {
            *big_f.add(i) = (mp_norm(*big_f.add(i), p) as u32) & 0x7FFF_FFFF;
        }
    } else {
        let p0 = PRIMES[0].p;
        let p1 = PRIMES[1].p;
        let p1_0i = PRIMES[1].p0i;
        let s = PRIMES[1].s;
        let pp = (p0 as u64) * (p1 as u64);
        let hpp = pp >> 1;
        for i in 0..n {
            // CRT with two primes on the coefficient of F.
            let x0 = *big_f.add(i);
            let x1 = *big_f.add(i + n);
            let x0m1 = x0.wrapping_sub(p1 & !tbmask(x0.wrapping_sub(p1)));
            let y = mp_mmul(mp_sub(x1, x0m1, p1), s, p1, p1_0i);
            let mut z = (x0 as u64).wrapping_add((p0 as u64).wrapping_mul(y as u64));
            z = z.wrapping_sub(pp & (hpp.wrapping_sub(z) >> 63).wrapping_neg());
            *big_f.add(i) = (z as u32) & 0x7FFF_FFFF;
            *big_f.add(i + n) = ((z >> 31) as u32) & 0x7FFF_FFFF;
        }
    }
}
