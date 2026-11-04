//! Computations modulo small 31-bit primes (key generation).
//!
//! Primes p satisfy (4/3)*2^30 < p < 2^31 and p = 1 mod 2048. Values are
//! held in `u32` in [0, p-1]. R = 2^32 mod p; the Montgomery representation
//! of x is x*R mod p. Polynomials in NTT representation use roots of X^n+1
//! (odd powers of g, a primitive 2048-th root of 1).
//!
//! Arrays are passed as raw pointers where the C code relies on strides or
//! in-place aliasing, and as slices elsewhere.

pub use crate::primes::PRIMES;

#[derive(Clone, Copy)]
pub struct SmallPrime {
    pub p: u32,
    /// -1/p mod 2^32
    pub p0i: u32,
    /// 2^64 mod p
    pub r2: u32,
    /// primitive 2048-th root of 1, Montgomery representation
    pub g: u32,
    /// 1/g mod p, Montgomery representation
    pub ig: u32,
    /// inverse mod p of the product of the previous primes (Montgomery)
    pub s: u32,
}

/// Expand the top bit of a 32-bit word into a full 32-bit mask.
#[inline(always)]
pub fn tbmask(x: u32) -> u32 {
    ((x as i32) >> 31) as u32
}

/// Given v in [-(p-1), +(p-1)], get v mod p in [0, p-1].
#[inline(always)]
pub fn mp_set(v: i32, p: u32) -> u32 {
    (v as u32).wrapping_add(p & tbmask(v as u32))
}

/// For x in [0, p-1], get its signed normalized value in [-p/2, +p/2].
#[inline(always)]
pub fn mp_norm(x: u32, p: u32) -> i32 {
    x.wrapping_sub(p & tbmask((p >> 1).wrapping_sub(x))) as i32
}

/// R = 2^32 mod p.
#[inline(always)]
pub fn mp_r(p: u32) -> u32 {
    (p << 1).wrapping_neg()
}

/// hR = 2^31 mod p.
#[inline(always)]
pub fn mp_hr(p: u32) -> u32 {
    (1u32 << 31).wrapping_sub(p)
}

#[inline(always)]
pub fn mp_add(a: u32, b: u32, p: u32) -> u32 {
    let d = a.wrapping_add(b).wrapping_sub(p);
    d.wrapping_add(p & tbmask(d))
}

#[inline(always)]
pub fn mp_sub(a: u32, b: u32, p: u32) -> u32 {
    let d = a.wrapping_sub(b);
    d.wrapping_add(p & tbmask(d))
}

#[inline(always)]
pub fn mp_half(a: u32, p: u32) -> u32 {
    a.wrapping_add(p & (a & 1).wrapping_neg()) >> 1
}

/// Montgomery multiplication: (a*b)/R mod p. Also works for out-of-range
/// inputs provided that a*b <= p + (p-1)*2^32.
#[inline(always)]
pub fn mp_mmul(a: u32, b: u32, p: u32, p0i: u32) -> u32 {
    let z = (a as u64).wrapping_mul(b as u64);
    let w = (z as u32).wrapping_mul(p0i);
    let d = ((z.wrapping_add((w as u64).wrapping_mul(p as u64)) >> 32) as u32).wrapping_sub(p);
    d.wrapping_add(p & tbmask(d))
}

/// Divide x by y modulo p; returns 0 if y is not invertible. x, y in
/// [0, p-1]; p odd.
pub fn mp_div(x: u32, y: u32, p: u32) -> u32 {
    let mut a = y;
    let mut b = p;
    let mut u = x;
    let mut v = 0u32;
    for _ in 0..62 {
        let a_odd = (a & 1).wrapping_neg();
        let swap = tbmask(a.wrapping_sub(b)) & a_odd;
        let t1 = swap & (a ^ b);
        a ^= t1;
        b ^= t1;
        let t2 = swap & (u ^ v);
        u ^= t2;
        v ^= t2;
        a = a.wrapping_sub(a_odd & b);
        u = mp_sub(u, a_odd & v, p);
        a >>= 1;
        u = mp_half(u, p);
    }
    // GCD is in b; it is 1 iff y was invertible.
    v & tbmask(b.wrapping_sub(2))
}

/// 2^(31*e) mod p (e is not secret).
pub fn mp_rx31(mut e: u32, p: u32, p0i: u32, r2: u32) -> u32 {
    let mut x = mp_half(r2, p);
    let mut d = 1u32;
    loop {
        if e & 1 != 0 {
            d = mp_mmul(d, x, p, p0i);
        }
        e >>= 1;
        if e == 0 {
            return d;
        }
        x = mp_mmul(x, x, p, p0i);
    }
}

/// Bit reversal over 10 bits.
#[inline(always)]
fn rev10(i: usize) -> usize {
    ((i as u32).reverse_bits() >> 22) as usize
}

/// gm[rev(i)] = g^i, igm[rev(i)] = (1/2)*(1/g)^i, for i < n (Montgomery).
pub fn mp_mkgmigm(logn: u32, gm: &mut [u32], igm: &mut [u32], mut g: u32, mut ig: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    for _ in logn..10 {
        g = mp_mmul(g, g, p, p0i);
        ig = mp_mmul(ig, ig, p, p0i);
    }
    let k = 10 - logn;
    let mut x1 = mp_r(p);
    let mut x2 = mp_hr(p);
    let mut i = 0usize;
    loop {
        let j = rev10(i << k);
        gm[j] = x1;
        igm[j] = x2;
        i += 1;
        if i >= n {
            break;
        }
        x1 = mp_mmul(x1, g, p, p0i);
        x2 = mp_mmul(x2, ig, p, p0i);
    }
}

/// Like `mp_mkgmigm`, computing only gm.
pub fn mp_mkgm(logn: u32, gm: &mut [u32], mut g: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    for _ in logn..10 {
        g = mp_mmul(g, g, p, p0i);
    }
    let k = 10 - logn;
    let mut x1 = mp_r(p);
    let mut i = 0usize;
    loop {
        gm[rev10(i << k)] = x1;
        i += 1;
        if i >= n {
            break;
        }
        x1 = mp_mmul(x1, g, p, p0i);
    }
}

/// Like `mp_mkgmigm`, computing only igm.
pub fn mp_mkigm(logn: u32, igm: &mut [u32], mut ig: u32, p: u32, p0i: u32) {
    let n = 1usize << logn;
    for _ in logn..10 {
        ig = mp_mmul(ig, ig, p, p0i);
    }
    let k = 10 - logn;
    let mut x2 = mp_hr(p);
    let mut i = 0usize;
    loop {
        igm[rev10(i << k)] = x2;
        i += 1;
        if i >= n {
            break;
        }
        x2 = mp_mmul(x2, ig, p, p0i);
    }
}

/// NTT of a polynomial, in place.
pub fn mp_ntt(logn: u32, a: &mut [u32], gm: &[u32], p: u32, p0i: u32) {
    if logn == 0 {
        return;
    }
    if logn == 1 {
        let s = gm[1];
        let x0 = a[0];
        let x1 = mp_mmul(a[1], s, p, p0i);
        a[0] = mp_add(x0, x1, p);
        a[1] = mp_sub(x0, x1, p);
        return;
    }
    let mut t: usize;
    let mut lm: u32;
    if logn & 1 != 0 {
        t = 1usize << (logn - 1);
        let s = gm[1];
        for i in 0..t {
            let x0 = a[i];
            let x1 = mp_mmul(a[i + t], s, p, p0i);
            a[i] = mp_add(x0, x1, p);
            a[i + t] = mp_sub(x0, x1, p);
        }
        lm = 1;
    } else {
        t = 1usize << logn;
        lm = 0;
    }

    // Two layers at a time.
    while lm + 3 < logn {
        let m = 1usize << lm;
        let ht = t >> 1;
        let qt = ht >> 1;
        let mut j0 = 0usize;
        for i in 0..m {
            let s = gm[i + m];
            let s0 = gm[((i + m) << 1) + 0];
            let s1 = gm[((i + m) << 1) + 1];
            for j in 0..qt {
                let k0 = j0 + j;
                let k1 = k0 + qt;
                let k2 = k0 + ht;
                let k3 = k0 + ht + qt;
                let (mut x0, mut x1, mut x2, mut x3) = (a[k0], a[k1], a[k2], a[k3]);

                let mut tt = mp_mmul(x2, s, p, p0i);
                x2 = mp_sub(x0, tt, p);
                x0 = mp_add(x0, tt, p);

                tt = mp_mmul(x3, s, p, p0i);
                x3 = mp_sub(x1, tt, p);
                x1 = mp_add(x1, tt, p);

                tt = mp_mmul(x1, s0, p, p0i);
                x1 = mp_sub(x0, tt, p);
                x0 = mp_add(x0, tt, p);

                tt = mp_mmul(x3, s1, p, p0i);
                x3 = mp_sub(x2, tt, p);
                x2 = mp_add(x2, tt, p);

                a[k0] = x0;
                a[k1] = x1;
                a[k2] = x2;
                a[k3] = x3;
            }
            j0 += t;
        }
        t = qt;
        lm += 2;
    }

    let m = 1usize << (logn - 2);
    for i in 0..m {
        let s = gm[i + m];
        let s0 = gm[((i + m) << 1) + 0];
        let s1 = gm[((i + m) << 1) + 1];
        let (mut x0, mut x1, mut x2, mut x3) =
            (a[(i << 2) + 0], a[(i << 2) + 1], a[(i << 2) + 2], a[(i << 2) + 3]);

        let mut tt = mp_mmul(x2, s, p, p0i);
        x2 = mp_sub(x0, tt, p);
        x0 = mp_add(x0, tt, p);

        tt = mp_mmul(x3, s, p, p0i);
        x3 = mp_sub(x1, tt, p);
        x1 = mp_add(x1, tt, p);

        tt = mp_mmul(x1, s0, p, p0i);
        x1 = mp_sub(x0, tt, p);
        x0 = mp_add(x0, tt, p);

        tt = mp_mmul(x3, s1, p, p0i);
        x3 = mp_sub(x2, tt, p);
        x2 = mp_add(x2, tt, p);

        a[(i << 2) + 0] = x0;
        a[(i << 2) + 1] = x1;
        a[(i << 2) + 2] = x2;
        a[(i << 2) + 3] = x3;
    }
}

/// Inverse NTT, in place.
pub fn mp_intt(logn: u32, a: &mut [u32], igm: &[u32], p: u32, p0i: u32) {
    if logn == 0 {
        return;
    }
    if logn == 1 {
        let s = igm[1];
        let x0 = a[0];
        let x1 = a[1];
        a[0] = mp_half(mp_add(x0, x1, p), p);
        a[1] = mp_mmul(mp_sub(x0, x1, p), s, p, p0i);
        return;
    }

    let qn = 1usize << (logn - 2);
    for i in 0..qn {
        let s0 = igm[((i + qn) << 1) + 0];
        let s1 = igm[((i + qn) << 1) + 1];
        let s = igm[i + qn];

        let (mut x0, mut x1, mut x2, mut x3) =
            (a[(i << 2) + 0], a[(i << 2) + 1], a[(i << 2) + 2], a[(i << 2) + 3]);

        let mut tt = mp_sub(x0, x1, p);
        x0 = mp_half(mp_add(x0, x1, p), p);
        x1 = mp_mmul(tt, s0, p, p0i);

        tt = mp_sub(x2, x3, p);
        x2 = mp_half(mp_add(x2, x3, p), p);
        x3 = mp_mmul(tt, s1, p, p0i);

        tt = mp_sub(x0, x2, p);
        x0 = mp_half(mp_add(x0, x2, p), p);
        x2 = mp_mmul(tt, s, p, p0i);

        tt = mp_sub(x1, x3, p);
        x1 = mp_half(mp_add(x1, x3, p), p);
        x3 = mp_mmul(tt, s, p, p0i);

        a[(i << 2) + 0] = x0;
        a[(i << 2) + 1] = x1;
        a[(i << 2) + 2] = x2;
        a[(i << 2) + 3] = x3;
    }

    let mut t = 4usize;
    let mut lm = 2u32;
    while lm + 1 < logn {
        let hm = 1usize << (logn - 1 - lm);
        let qm = hm >> 1;
        let dt = t << 1;
        let ft = t << 2;
        let mut j0 = 0usize;
        for i in 0..qm {
            let s0 = igm[((i + qm) << 1) + 0];
            let s1 = igm[((i + qm) << 1) + 1];
            let s = igm[i + qm];
            for j in 0..t {
                let k0 = j0 + j;
                let k1 = k0 + t;
                let k2 = k0 + dt;
                let k3 = k0 + t + dt;
                let (mut x0, mut x1, mut x2, mut x3) = (a[k0], a[k1], a[k2], a[k3]);

                let mut tt = mp_sub(x0, x1, p);
                x0 = mp_half(mp_add(x0, x1, p), p);
                x1 = mp_mmul(tt, s0, p, p0i);

                tt = mp_sub(x2, x3, p);
                x2 = mp_half(mp_add(x2, x3, p), p);
                x3 = mp_mmul(tt, s1, p, p0i);

                tt = mp_sub(x0, x2, p);
                x0 = mp_half(mp_add(x0, x2, p), p);
                x2 = mp_mmul(tt, s, p, p0i);

                tt = mp_sub(x1, x3, p);
                x1 = mp_half(mp_add(x1, x3, p), p);
                x3 = mp_mmul(tt, s, p, p0i);

                a[k0] = x0;
                a[k1] = x1;
                a[k2] = x2;
                a[k3] = x3;
            }
            j0 += ft;
        }
        t = ft;
        lm += 2;
    }

    if logn & 1 != 0 {
        let hn = 1usize << (logn - 1);
        let s = igm[1];
        for i in 0..hn {
            let x0 = a[i];
            let x1 = a[i + hn];
            a[i] = mp_half(mp_add(x0, x1, p), p);
            a[i + hn] = mp_mmul(mp_sub(x0, x1, p), s, p, p0i);
        }
    }
}
