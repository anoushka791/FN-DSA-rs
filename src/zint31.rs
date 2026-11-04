//! Big-integer support for key generation.
//!
//! Big integers are sequences of 31-bit limbs, each stored in a `u32` with
//! the top bit cleared, low limb first; signed integers use two's complement
//! (bit 30 of the last limb is the sign). Integers may be strided in memory
//! (`stride` words between limbs), so these functions use raw pointers.
//! RNS representation holds x mod p for the primes of `PRIMES`.

use crate::mp31::{mp_add, mp_half, mp_mmul, mp_sub, tbmask, PRIMES};
use core::ptr;

/// Multiply the big integer `m` (stride 1, unsigned) by a small value
/// `x < 2^31`. Returns the carry word.
///
/// # Safety
/// `m` must be valid for `len` reads and writes.
pub unsafe fn zint_mul_small(m: *mut u32, len: usize, x: u32) -> u32 {
    let mut cc = 0u32;
    for i in 0..len {
        let z = (*m.add(i) as u64) * (x as u64) + cc as u64;
        *m.add(i) = (z as u32) & 0x7FFF_FFFF;
        cc = (z >> 31) as u32;
    }
    cc
}

/// Reduce an unsigned big integer modulo the small prime p.
///
/// # Safety
/// `d` must be valid for reads at `d + i*stride`, `i < len`.
pub unsafe fn zint_mod_small_unsigned(
    d: *const u32,
    len: usize,
    stride: usize,
    p: u32,
    p0i: u32,
    r2: u32,
) -> u32 {
    let mut x = 0u32;
    let z = mp_half(r2, p);
    let mut d = d.add(len * stride);
    let mut i = len;
    while i > 0 {
        d = d.sub(stride);
        let mut w = (*d).wrapping_sub(p);
        w = w.wrapping_add(p & tbmask(w));
        x = mp_mmul(x, z, p, p0i);
        x = mp_add(x, w, p);
        i -= 1;
    }
    x
}

/// Same as `zint_mod_small_unsigned`, with signed convention; `rx` is
/// 2^(31*len) mod p.
///
/// # Safety
/// See `zint_mod_small_unsigned`.
pub unsafe fn zint_mod_small_signed(
    d: *const u32,
    len: usize,
    stride: usize,
    p: u32,
    p0i: u32,
    r2: u32,
    rx: u32,
) -> u32 {
    if len == 0 {
        return 0;
    }
    let mut z = zint_mod_small_unsigned(d, len, stride, p, p0i, r2);
    z = mp_sub(z, rx & (*d.add((len - 1) * stride) >> 30).wrapping_neg(), p);
    z
}

/// Add s*y to x. x and y initially have `len` words; the new x has `len+1`
/// words. `s` must fit on 31 bits. x uses `xstride`, y has stride 1.
///
/// # Safety
/// No overlap between x and y; x valid for `len+1` strided words.
pub unsafe fn zint_add_mul_small(x: *mut u32, len: usize, xstride: usize, y: *const u32, s: u32) {
    let mut x = x;
    let mut cc = 0u32;
    for i in 0..len {
        let xw = *x;
        let yw = *y.add(i);
        let z = (yw as u64) * (s as u64) + (xw as u64) + (cc as u64);
        *x = (z as u32) & 0x7FFF_FFFF;
        cc = (z >> 31) as u32;
        x = x.add(xstride);
    }
    *x = cc;
}

/// Normalize a modular integer around 0: if x > m/2 then x is replaced with
/// x - m. x is unsigned on input and signed on output.
///
/// # Safety
/// x valid for `len` strided words; m valid for `len` words.
pub unsafe fn zint_norm_zero(x: *mut u32, len: usize, xstride: usize, m: *const u32) {
    let mut r = 0u32;
    let mut bb = 0u32;
    let mut xp = x.add(len * xstride);
    let mut i = len;
    while i > 0 {
        i -= 1;
        xp = xp.sub(xstride);
        // The two words to compare, both over 31 bits exactly.
        let wx = *xp;
        let wp = (*m.add(i) >> 1) | (bb << 30);
        bb = *m.add(i) & 1;

        // cc = -1, 0 or 1 depending on whether wp is lower than, equal to,
        // or greater than wx.
        let cc = wp.wrapping_sub(wx);
        let cc = (cc.wrapping_neg() >> 31) | (cc >> 31).wrapping_neg();

        // If r != 0 then it is 1 or -1, and we keep it; else we take cc.
        r |= cc & ((r & 1).wrapping_sub(1));
    }

    // r = -1, 0 or 1 depending on (m-1)/2 vs x: subtract iff r = -1.
    let mut cc = 0u32;
    let s = tbmask(r);
    let mut xq = xp;
    for j in 0..len {
        let mut xw = *xq;
        let w = xw.wrapping_sub(*m.add(j)).wrapping_sub(cc);
        cc = w >> 31;
        xw ^= ((w & 0x7FFF_FFFF) ^ xw) & s;
        *xq = xw;
        xq = xq.add(xstride);
    }
}

/// Rebuild integers from their RNS representation. See the C documentation
/// of `zint_rebuild_CRT`.
///
/// # Safety
/// `xx` valid for `num_sets * n * xlen` words, `tmp` for `xlen` words.
pub unsafe fn zint_rebuild_crt(
    xx: *mut u32,
    xlen: usize,
    n: usize,
    num_sets: usize,
    normalize_signed: bool,
    tmp: *mut u32,
) {
    let mut uu = 0usize;
    *tmp = PRIMES[0].p;
    for i in 1..xlen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        let s = PRIMES[i].s;
        uu += n;
        let mut kk = 0usize;
        for _ in 0..num_sets {
            for j in 0..n {
                // xp = x mod p; xq = (x mod q) mod p
                let xp = *xx.add(kk + j + uu);
                let xq = zint_mod_small_unsigned(xx.add(kk + j), i, n, p, p0i, r2);
                // New value: (x mod q) + q*(s*(xp - xq) mod p)
                let xr = mp_mmul(s, mp_sub(xp, xq, p), p, p0i);
                zint_add_mul_small(xx.add(kk + j), i, n, tmp, xr);
            }
            kk += n * xlen;
        }
        // Update the product of primes in tmp.
        *tmp.add(i) = zint_mul_small(tmp, i, p);
    }

    if normalize_signed {
        let mut kk = 0usize;
        for _ in 0..num_sets {
            for j in 0..n {
                zint_norm_zero(xx.add(kk + j), xlen, n, tmp);
            }
            kk += n * xlen;
        }
    }
}

/// Negate a big integer (stride 1) iff `ctl == 0xFFFFFFFF`.
///
/// # Safety
/// `a` valid for `len` reads/writes.
pub unsafe fn zint_negate(a: *mut u32, len: usize, ctl: u32) {
    let mut cc = ctl.wrapping_neg();
    let m = ctl >> 1;
    for i in 0..len {
        let aw = (*a.add(i) ^ m).wrapping_add(cc);
        *a.add(i) = aw & 0x7FFF_FFFF;
        cc = aw >> 31;
    }
}

/// a <- (a*xa + b*xb)/2^31, b <- (a*ya + b*yb)/2^31 (low bits dropped);
/// negated if negative. Returns 0..3 (bit 0: a negated, bit 1: b negated).
unsafe fn zint_co_reduce(
    a: *mut u32,
    b: *mut u32,
    len: usize,
    xa: i64,
    xb: i64,
    ya: i64,
    yb: i64,
) -> u32 {
    let mut cca = 0i64;
    let mut ccb = 0i64;
    for i in 0..len {
        let wa = *a.add(i) as u64;
        let wb = *b.add(i) as u64;
        let za = wa
            .wrapping_mul(xa as u64)
            .wrapping_add(wb.wrapping_mul(xb as u64))
            .wrapping_add(cca as u64);
        let zb = wa
            .wrapping_mul(ya as u64)
            .wrapping_add(wb.wrapping_mul(yb as u64))
            .wrapping_add(ccb as u64);
        if i > 0 {
            *a.add(i - 1) = (za as u32) & 0x7FFF_FFFF;
            *b.add(i - 1) = (zb as u32) & 0x7FFF_FFFF;
        }
        cca = (za as i64) >> 31;
        ccb = (zb as i64) >> 31;
    }
    *a.add(len - 1) = (cca as u32) & 0x7FFF_FFFF;
    *b.add(len - 1) = (ccb as u32) & 0x7FFF_FFFF;

    let nega = (cca >> 63) as u32;
    let negb = (ccb >> 63) as u32;
    zint_negate(a, len, nega);
    zint_negate(b, len, negb);
    nega.wrapping_neg() | (negb.wrapping_neg() << 1)
}

/// Finish modular reduction: if neg = 1 then -m <= a < 0; if neg = 0 then
/// 0 <= a < 2*m (top word of a may use 32 bits). m odd; both of length len.
unsafe fn zint_finish_mod(a: *mut u32, m: *const u32, len: usize, neg: u32) {
    // First pass: compare a (assumed nonnegative) with m.
    let mut cc = 0u32;
    for i in 0..len {
        cc = (*a.add(i)).wrapping_sub(*m.add(i)).wrapping_sub(cc) >> 31;
    }

    // neg = 1: add m. neg = 0, cc = 0: subtract m. neg = 0, cc = 1: nothing.
    let xm = neg.wrapping_neg() >> 1;
    let ym = (neg | (1u32.wrapping_sub(cc))).wrapping_neg();
    cc = neg;
    for i in 0..len {
        let mw = (*m.add(i) ^ xm) & ym;
        let aw = (*a.add(i)).wrapping_sub(mw).wrapping_sub(cc);
        *a.add(i) = aw & 0x7FFF_FFFF;
        cc = aw >> 31;
    }
}

/// a <- (a*xa + b*xb)/2^31 mod m, b <- (a*ya + b*yb)/2^31 mod m, for odd m
/// with m0i = -1/m[0] mod 2^31. a, b, m have length len and stride 1.
unsafe fn zint_co_reduce_mod(
    a: *mut u32,
    b: *mut u32,
    m: *const u32,
    len: usize,
    m0i: u32,
    xa: i64,
    xb: i64,
    ya: i64,
    yb: i64,
) {
    // Four combined Montgomery multiplications.
    let mut cca = 0i64;
    let mut ccb = 0i64;
    let a0 = *a;
    let b0 = *b;
    let fa = (a0
        .wrapping_mul(xa as u32)
        .wrapping_add(b0.wrapping_mul(xb as u32))
        .wrapping_mul(m0i))
        & 0x7FFF_FFFF;
    let fb = (a0
        .wrapping_mul(ya as u32)
        .wrapping_add(b0.wrapping_mul(yb as u32))
        .wrapping_mul(m0i))
        & 0x7FFF_FFFF;
    for i in 0..len {
        let wa = *a.add(i) as u64;
        let wb = *b.add(i) as u64;
        let mi = *m.add(i) as u64;
        let za = wa
            .wrapping_mul(xa as u64)
            .wrapping_add(wb.wrapping_mul(xb as u64))
            .wrapping_add(mi.wrapping_mul(fa as u64))
            .wrapping_add(cca as u64);
        let zb = wa
            .wrapping_mul(ya as u64)
            .wrapping_add(wb.wrapping_mul(yb as u64))
            .wrapping_add(mi.wrapping_mul(fb as u64))
            .wrapping_add(ccb as u64);
        if i > 0 {
            *a.add(i - 1) = (za as u32) & 0x7FFF_FFFF;
            *b.add(i - 1) = (zb as u32) & 0x7FFF_FFFF;
        }
        cca = (za as i64) >> 31;
        ccb = (zb as i64) >> 31;
    }
    *a.add(len - 1) = cca as u32;
    *b.add(len - 1) = ccb as u32;

    // -m <= a, b < 2*m: add or subtract the modulus as required.
    zint_finish_mod(a, m, len, ((cca as u64) >> 63) as u32);
    zint_finish_mod(b, m, len, ((ccb as u64) >> 63) as u32);
}

/// For odd x, compute -1/x mod 2^31.
#[inline(always)]
fn mp_ninv31(x: u32) -> u32 {
    let mut y = 2u32.wrapping_sub(x);
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y = y.wrapping_mul(2u32.wrapping_sub(x.wrapping_mul(y)));
    y.wrapping_neg() & 0x7FFF_FFFF
}

/// Compute GCD(x, y) for odd x and y; returns true (and u, v with
/// x*u - y*v = 1, 0 <= u <= y, 0 <= v <= x) iff the GCD is 1. All values
/// have `len` words and stride 1; `tmp` needs `4*len` words; none may overlap.
///
/// # Safety
/// Pointers must be valid for the stated lengths and not overlap.
pub unsafe fn zint_bezout(
    u: *mut u32,
    v: *mut u32,
    x: *const u32,
    y: *const u32,
    len: usize,
    tmp: *mut u32,
) -> bool {
    if len == 0 {
        return false;
    }

    // Optimized binary GCD (https://eprint.iacr.org/2020/972), with four
    // coefficients to get a full Bezout relation. Invariants:
    //   a = x*u0 - y*v0,  b = x*u1 - y*v1
    //   0 <= a <= x, 0 <= b <= y, 0 <= u0 < y, 0 <= v0 < x,
    //   0 <= u1 <= y, 0 <= v1 < x
    let u0 = tmp;
    let v0 = u0.add(len);
    let u1 = u;
    let v1 = v;
    let a = v0.add(len);
    let b = a.add(len);

    let x0i = mp_ninv31(*x);
    let y0i = mp_ninv31(*y);

    // Initial values: a = x, u0 = 1, v0 = 0; b = y, u1 = y, v1 = x - 1.
    ptr::copy_nonoverlapping(x, a, len);
    ptr::copy_nonoverlapping(y, b, len);
    *u0 = 1;
    ptr::write_bytes(u0.add(1), 0, len - 1);
    ptr::write_bytes(v0, 0, len);
    ptr::copy_nonoverlapping(y, u1, len);
    ptr::copy_nonoverlapping(x, v1, len);
    *v1 = (*v1).wrapping_sub(1);

    let mut num = 62 * (len as u32) + 31;
    while num >= 31 {
        // Top 32 bits of a and b (see the C code for the case analysis).
        let mut c0 = 0xFFFF_FFFFu32;
        let mut c1 = 0xFFFF_FFFFu32;
        let mut cp = 0xFFFF_FFFFu32;
        let mut a0 = 0u32;
        let mut a1 = 0u32;
        let mut b0 = 0u32;
        let mut b1 = 0u32;
        let mut j = len;
        while j > 0 {
            j -= 1;
            let aw = *a.add(j);
            let bw = *b.add(j);
            a1 ^= c1 & (a1 ^ aw);
            a0 ^= c0 & (a0 ^ aw);
            b1 ^= c1 & (b1 ^ bw);
            b0 ^= c0 & (b0 ^ bw);
            cp = c0;
            c0 = c1;
            c1 &= (((aw | bw).wrapping_add(0x7FFF_FFFF)) >> 31).wrapping_sub(1);
        }

        let s = (a1 | b1 | ((cp & c0) >> 1)).leading_zeros();
        let mut ha = (a1 << s) | (a0 >> (31 - s));
        let mut hb = (b1 << s) | (b0 >> (31 - s));

        // If j <= 62, use the non-aligned bits.
        ha ^= cp & (ha ^ a1);
        hb ^= cp & (hb ^ b1);

        // If j <= 31, clear the upper bits.
        ha &= !c0;
        hb &= !c0;

        // Approximate values xa and xb (63 bits each).
        let mut xa = ((ha as u64) << 31) | (*a as u64);
        let mut xb = ((hb as u64) << 31) | (*b as u64);

        // Reduction factors: a' = a*pa + b*pb, b' = a*qa + b*qb, both
        // multiples of 2^31 and only marginally larger than a and b.
        let mut fg0 = 1u64;
        let mut fg1 = 1u64 << 32;
        for _ in 0..31 {
            let a_odd = (xa & 1).wrapping_neg();
            let dx = ((xa.wrapping_sub(xb)) as i64 >> 63) as u64;
            let swap = a_odd & dx;
            let t1 = swap & (xa ^ xb);
            xa ^= t1;
            xb ^= t1;
            let t2 = swap & (fg0 ^ fg1);
            fg0 ^= t2;
            fg1 ^= t2;
            xa = xa.wrapping_sub(a_odd & xb);
            fg0 = fg0.wrapping_sub(a_odd & fg1);
            xa >>= 1;
            fg1 <<= 1;
        }

        // Split update factors.
        fg0 = fg0.wrapping_add(0x7FFF_FFFF_7FFF_FFFF);
        fg1 = fg1.wrapping_add(0x7FFF_FFFF_7FFF_FFFF);
        let mut f0 = (fg0 & 0xFFFF_FFFF) as i64 - 0x7FFF_FFFF;
        let mut g0 = (fg0 >> 32) as i64 - 0x7FFF_FFFF;
        let mut f1 = (fg1 & 0xFFFF_FFFF) as i64 - 0x7FFF_FFFF;
        let mut g1 = (fg1 >> 32) as i64 - 0x7FFF_FFFF;

        // Apply the update factors.
        let negab = zint_co_reduce(a, b, len, f0, g0, f1, g1);
        f0 = f0.wrapping_sub((f0.wrapping_add(f0)) & -((negab & 1) as i64));
        g0 = g0.wrapping_sub((g0.wrapping_add(g0)) & -((negab & 1) as i64));
        f1 = f1.wrapping_sub((f1.wrapping_add(f1)) & -((negab >> 1) as i64));
        g1 = g1.wrapping_sub((g1.wrapping_add(g1)) & -((negab >> 1) as i64));
        zint_co_reduce_mod(u0, u1, y, len, y0i, f0, g0, f1, g1);
        zint_co_reduce_mod(v0, v1, x, len, x0i, f0, g0, f1, g1);

        num -= 31;
    }

    // b contains GCD(x, y), provided that x and y were indeed odd.
    let mut r = *b ^ 1;
    for j in 1..len {
        r |= *b.add(j);
    }
    r |= (*x & *y & 1) ^ 1;
    1u32.wrapping_sub((r | r.wrapping_neg()) >> 31) != 0
}

/// Add k*(2^sc)*y to x (result truncated to xlen words). sc = 31*sch + scl,
/// scl in [0, 30]. xlen >= ylen. Signed convention, same stride.
///
/// # Safety
/// x valid for `xlen` strided words, y for `ylen`.
pub unsafe fn zint_add_scaled_mul_small(
    x: *mut u32,
    xlen: usize,
    y: *const u32,
    ylen: usize,
    stride: usize,
    k: i32,
    sch: u32,
    scl: u32,
) {
    if ylen == 0 {
        return;
    }
    let ysign = (*y.add(stride * (ylen - 1)) >> 30).wrapping_neg() >> 1;
    let mut tw = 0u32;
    let mut cc = 0i32;
    let mut xp = x.add(sch as usize * stride);
    let mut yp = y;
    let mut yl = ylen;
    for _ in (sch as usize)..xlen {
        // Next word of (2^sc)*y.
        let wy;
        if yl > 0 {
            wy = *yp;
            yp = yp.add(stride);
            yl -= 1;
        } else {
            wy = ysign;
        }
        let wys = ((wy << scl) & 0x7FFF_FFFF) | tw;
        tw = wy >> (31 - scl);

        // The expression does not overflow.
        let z = ((wys as i64) * (k as i64) + (*xp as i64) + (cc as i64)) as u64;
        *xp = (z as u32) & 0x7FFF_FFFF;
        xp = xp.add(stride);

        // New carry word is a signed right shift of z.
        cc = ((z >> 31) as u32) as i32;
    }
}

/// Subtract y*2^sc from x (specialization of `zint_add_scaled_mul_small`
/// with k = -1).
///
/// # Safety
/// See `zint_add_scaled_mul_small`.
pub unsafe fn zint_sub_scaled(
    x: *mut u32,
    xlen: usize,
    y: *const u32,
    ylen: usize,
    stride: usize,
    sch: u32,
    scl: u32,
) {
    if ylen == 0 {
        return;
    }
    let ysign = (*y.add(stride * (ylen - 1)) >> 30).wrapping_neg() >> 1;
    let mut tw = 0u32;
    let mut cc = 0u32;
    let mut xp = x.add(sch as usize * stride);
    let mut yp = y;
    let mut yl = ylen;
    for _ in (sch as usize)..xlen {
        let wy;
        if yl > 0 {
            wy = *yp;
            yp = yp.add(stride);
            yl -= 1;
        } else {
            wy = ysign;
        }
        let wys = ((wy << scl) & 0x7FFF_FFFF) | tw;
        tw = wy >> (31 - scl);

        let w = (*xp).wrapping_sub(wys).wrapping_sub(cc);
        *xp = w & 0x7FFF_FFFF;
        cc = w >> 31;
        xp = xp.add(stride);
    }
}
