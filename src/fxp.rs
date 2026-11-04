//! Fixed-point numbers for key generation.
//!
//! A value is a `u64` in two's complement: the top 32 bits are the integral
//! part and the low 32 bits the fractional part. All arithmetic wraps, and is
//! bit-exact with the C code (no floating point is involved).

use crate::fxp_tables::GM_TAB;

pub type Fxr = u64;

#[inline(always)]
pub fn fxr_of(j: i32) -> Fxr {
    ((j as i64) << 32) as u64
}

#[inline(always)]
pub fn fxr_add(x: Fxr, y: Fxr) -> Fxr {
    x.wrapping_add(y)
}

#[inline(always)]
pub fn fxr_sub(x: Fxr, y: Fxr) -> Fxr {
    x.wrapping_sub(y)
}

#[inline(always)]
pub fn fxr_neg(x: Fxr) -> Fxr {
    x.wrapping_neg()
}

#[inline(always)]
pub fn fxr_mul(x: Fxr, y: Fxr) -> Fxr {
    (((x as i64 as i128) * (y as i64 as i128)) >> 32) as u64
}

#[inline(always)]
pub fn fxr_sqr(x: Fxr) -> Fxr {
    let t = x as i64 as i128;
    ((t * t) >> 32) as u64
}

#[inline(always)]
pub fn fxr_round(x: Fxr) -> i32 {
    (x.wrapping_add(0x8000_0000) as i64 >> 32) as i32
}

/// Divide x by 2^e (with rounding).
#[inline(always)]
pub fn fxr_div2e(x: Fxr, e: u32) -> Fxr {
    let x = x.wrapping_add((1u64 << e) >> 1);
    ((x as i64) >> e) as u64
}

#[inline(always)]
pub fn fxr_mul2e(x: Fxr, e: u32) -> Fxr {
    x << e
}

#[inline(always)]
pub fn fxr_lt(x: Fxr, y: Fxr) -> bool {
    (x as i64) < (y as i64)
}

/// Bit-by-bit division on fixed-point values.
pub fn inner_fxr_div(x: u64, y: u64) -> u64 {
    // Absolute values and signs; from now on x and y fit on 63 bits.
    let mut sx = x >> 63;
    let x = (x ^ sx.wrapping_neg()).wrapping_add(sx);
    let sy = y >> 63;
    let y = (y ^ sy.wrapping_neg()).wrapping_add(sy);

    // Bit-by-bit division; the numerator starts at x*2^31.
    let mut q = 0u64;
    let mut num = x >> 31;
    for i in (33..=63).rev() {
        let b = 1 - (num.wrapping_sub(y) >> 63);
        q |= b << i;
        num = num.wrapping_sub(y & b.wrapping_neg());
        num <<= 1;
        num |= (x >> (i - 33)) & 1;
    }
    for i in (0..=32).rev() {
        let b = 1 - (num.wrapping_sub(y) >> 63);
        q |= b << i;
        num = num.wrapping_sub(y & b.wrapping_neg());
        num <<= 1;
    }

    // Rounding: if the remainder is at least y/2 (scaled), add 2^(-32).
    let b = 1 - (num.wrapping_sub(y) >> 63);
    q = q.wrapping_add(b);

    // Sign management.
    sx ^= sy;
    (q ^ sx.wrapping_neg()).wrapping_add(sx)
}

#[inline(always)]
pub fn fxr_div(x: Fxr, y: Fxr) -> Fxr {
    inner_fxr_div(x, y)
}

/// Complex value (re, im).
pub type Fxc = (Fxr, Fxr);

#[inline(always)]
fn fxc_add(x: Fxc, y: Fxc) -> Fxc {
    (fxr_add(x.0, y.0), fxr_add(x.1, y.1))
}

#[inline(always)]
fn fxc_sub(x: Fxc, y: Fxc) -> Fxc {
    (fxr_sub(x.0, y.0), fxr_sub(x.1, y.1))
}

#[inline(always)]
fn fxc_half(x: Fxc) -> Fxc {
    (fxr_div2e(x.0, 1), fxr_div2e(x.1, 1))
}

#[inline(always)]
fn fxc_conj(x: Fxc) -> Fxc {
    (x.0, fxr_neg(x.1))
}

/// Complex multiplication with three real multiplications (the imaginary
/// part may differ slightly from a*d + b*c, because of the approximation).
#[inline(always)]
fn fxc_mul(x: Fxc, y: Fxc) -> Fxc {
    let z0 = fxr_mul(x.0, y.0);
    let z1 = fxr_mul(x.1, y.1);
    let z2 = fxr_mul(fxr_add(x.0, x.1), fxr_add(y.0, y.1));
    (fxr_sub(z0, z1), fxr_sub(z2, fxr_add(z0, z1)))
}

/// Convert a real vector to FFT representation (real parts in f[i],
/// imaginary parts in f[i + n/2]).
pub fn vect_fft(logn: u32, f: &mut [Fxr]) {
    let hn = 1usize << (logn - 1);
    let mut t = hn;
    for lm in 1..logn {
        let m = 1usize << lm;
        let ht = t >> 1;
        let mut j0 = 0usize;
        let hm = m >> 1;
        for i in 0..hm {
            let s = GM_TAB[m + i];
            for j in j0..j0 + ht {
                let x = (f[j], f[j + hn]);
                let y = fxc_mul(s, (f[j + ht], f[j + ht + hn]));
                let z1 = fxc_add(x, y);
                f[j] = z1.0;
                f[j + hn] = z1.1;
                let z2 = fxc_sub(x, y);
                f[j + ht] = z2.0;
                f[j + ht + hn] = z2.1;
            }
            j0 += t;
        }
        t = ht;
    }
}

/// Convert back from FFT representation to a real vector.
pub fn vect_ifft(logn: u32, f: &mut [Fxr]) {
    let hn = 1usize << (logn - 1);
    let mut ht = 1usize;
    for lm in (1..logn).rev() {
        let m = 1usize << lm;
        let t = ht << 1;
        let mut j0 = 0usize;
        let hm = m >> 1;
        for i in 0..hm {
            let s = fxc_conj(GM_TAB[m + i]);
            for j in j0..j0 + ht {
                let x = (f[j], f[j + hn]);
                let y = (f[j + ht], f[j + ht + hn]);
                let z1 = fxc_half(fxc_add(x, y));
                f[j] = z1.0;
                f[j + hn] = z1.1;
                let z2 = fxc_mul(s, fxc_half(fxc_sub(x, y)));
                f[j + ht] = z2.0;
                f[j + ht + hn] = z2.1;
            }
            j0 += t;
        }
        ht = t;
    }
}

pub fn vect_set(logn: u32, d: &mut [Fxr], f: &[i8]) {
    let n = 1usize << logn;
    for i in 0..n {
        d[i] = fxr_of(f[i] as i32);
    }
}

pub fn vect_mul_realconst(logn: u32, a: &mut [Fxr], c: Fxr) {
    let n = 1usize << logn;
    for v in a[..n].iter_mut() {
        *v = fxr_mul(*v, c);
    }
}

/// a <- a*b (FFT representation; no overlap).
pub fn vect_mul_fft(logn: u32, a: &mut [Fxr], b: &[Fxr]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        let z = fxc_mul((a[i], a[i + hn]), (b[i], b[i + hn]));
        a[i] = z.0;
        a[i + hn] = z.1;
    }
}

pub fn vect_adj_fft(logn: u32, a: &mut [Fxr]) {
    let hn = 1usize << (logn - 1);
    let n = 1usize << logn;
    for v in a[hn..n].iter_mut() {
        *v = fxr_neg(*v);
    }
}

/// a <- a*b, with b self-adjoint (half-length, real).
pub fn vect_mul_selfadj_fft(logn: u32, a: &mut [Fxr], b: &[Fxr]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        a[i] = fxr_mul(a[i], b[i]);
        a[i + hn] = fxr_mul(a[i + hn], b[i]);
    }
}

/// a <- a/b, with b self-adjoint (half-length, real).
pub fn vect_div_selfadj_fft(logn: u32, a: &mut [Fxr], b: &[Fxr]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        a[i] = fxr_div(a[i], b[i]);
        a[i + hn] = fxr_div(a[i + hn], b[i]);
    }
}

/// d <- (2^e)/(a*adj(a) + b*adj(b)) (self-adjoint, half-length).
pub fn vect_invnorm_fft(logn: u32, d: &mut [Fxr], a: &[Fxr], b: &[Fxr], e: u32) {
    let hn = 1usize << (logn - 1);
    let fe = fxr_of(1i32 << e);
    for i in 0..hn {
        let z = fxr_add(
            fxr_add(fxr_sqr(a[i]), fxr_sqr(a[i + hn])),
            fxr_add(fxr_sqr(b[i]), fxr_sqr(b[i + hn])),
        );
        d[i] = fxr_div(fe, z);
    }
}

/// a <- (2^e)*adj(a)/(a*adj(a)).
pub fn vect_inv_mul2e_fft(logn: u32, a: &mut [Fxr], e: u32) {
    let hn = 1usize << (logn - 1);
    for u in 0..hn {
        let re = a[u];
        let im = fxr_neg(a[u + hn]);
        let z = fxr_add(fxr_sqr(re), fxr_sqr(im));
        a[u] = fxr_div(fxr_mul2e(re, e), z);
        a[u + hn] = fxr_div(fxr_mul2e(im, e), z);
    }
}
