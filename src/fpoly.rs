//! Floating-point polynomials (modulo X^n+1), in FFT representation.
//!
//! Values are IEEE-754 binary64 with round-to-nearest-even. Rust never fuses
//! a multiplication and an addition unless asked to, so the sequence of
//! operations below produces bit-for-bit the same results as the C code.
//!
//! In FFT representation a polynomial of degree n has n/2 complex
//! coefficients, with the real parts in `f[i]` and imaginary parts in
//! `f[i + n/2]`.
//!
//! Some functions convert integer polynomials *in place* inside the storage
//! of the floating-point result (the C implementation does this to save RAM);
//! those take raw pointers and are `unsafe`.

use crate::sign_tables::{
    IRNJ, RNJ, SIN, X_MINUS_1_OVER_X_PLUS_1_TIMES_I, X_PLUS_1_SQUARED_OVER_X,
};
use core::ptr;

pub const FPR_SQRT2: f64 = f64::from_bits(0x3FF6A09E667F3BCD);
/// q^2
pub const FPR_QQ: f64 = f64::from_bits(0x41A200C002000000);

#[inline(always)]
fn add_sub(x: f64, y: f64) -> (f64, f64) {
    (x + y, x - y)
}

/// Complex multiplication with four real multiplications.
#[inline(always)]
fn fpc_mul(a_re: f64, a_im: f64, b_re: f64, b_im: f64) -> (f64, f64) {
    (a_re * b_re - a_im * b_im, a_re * b_im + a_im * b_re)
}

/// Convert a self-adjoint integer polynomial (n/2 `i32` values at `a`) into
/// FFT representation (n/2 values at `d`). `a` may overlap `d`; it is most
/// efficient if `a` is at the end of the `d` area, i.e. at
/// `(d + n/2) as *i32 - n/2`.
///
/// # Safety
/// `d` must be valid for `n/2` `f64` writes, `a` for `n/2` `i32` reads.
pub unsafe fn fft_int32_selfadj(logn: u32, d: *mut f64, a: *const i32) {
    if logn == 1 {
        *d = *a as f64;
        return;
    }
    let n = 1usize << logn;
    let hn = n >> 1;
    let qn = hn >> 1;

    // Move the source to the end of the destination area.
    let aa = (d.add(hn) as *mut i32).sub(hn);
    if a != aa as *const i32 {
        ptr::copy(a, aa, hn);
    }

    // First pass: all even/odd splits, multiplying the odd half by X + 1.
    for lk in 0..(logn - 2) {
        let k = 1usize << lk;
        let dk = k << 1;
        let mut i = hn - dk;
        while i > 0 {
            for j in 0..k {
                let v = (*aa.add(i + j + k)).wrapping_add(*aa.add(i + j - k));
                *aa.add(i + j + k) = v;
            }
            i -= dk;
        }
        for j in 0..k {
            *aa.add(j + k) = (*aa.add(j + k)).wrapping_shl(1);
        }
    }

    // Convert to floating-point (ascending order: the writes never reach
    // unread source values).
    for j in 0..hn {
        let v = *aa.add(j) as f64;
        *d.add(j) = v;
    }

    // Second pass: final split and merge, folded in one multiplication by
    // sqrt(2).
    for j in 0..qn {
        let a0z = *d.add(j);
        let b1z = *d.add(j + qn);
        let xz = b1z * FPR_SQRT2;
        let (s, t) = add_sub(a0z, xz);
        *d.add(j) = s;
        *d.add(j + qn) = t;
    }

    // Third pass: keep merging back polynomials.
    for lm in 1..(logn - 1) {
        let lk = logn - 2 - lm;
        let k = 1usize << lk;
        let m = 1usize << lm;
        for j in 0..m {
            let mj = j << (lk + 1);
            let rz = RNJ[j + m];
            for i in 0..k {
                let a0z = *d.add(mj + i);
                let b1z = *d.add(mj + i + k);
                let xz = b1z * rz;
                let (s, t) = add_sub(a0z, xz);
                *d.add(mj + i) = s;
                *d.add(mj + i + k) = t;
            }
        }
    }
}

/// Same as `fft_int32_selfadj`, with `i16` source values.
///
/// # Safety
/// `d` must be valid for `n/2` `f64` writes, `a` for `n/2` `i16` reads.
pub unsafe fn fft_int16_selfadj(logn: u32, d: *mut f64, a: *const i16) {
    if logn == 1 {
        *d = *a as f64;
        return;
    }
    let hn = 1usize << (logn - 1);

    let aa = (d.add(hn) as *mut i16).sub(hn);
    if a != aa as *const i16 {
        ptr::copy(a, aa, hn);
    }

    // Widen to 32 bits, in place (ascending order is safe).
    let bb = (d.add(hn) as *mut i32).sub(hn);
    for i in 0..hn {
        let v = *aa.add(i) as i32;
        *bb.add(i) = v;
    }

    fft_int32_selfadj(logn, d, bb);
}

/// Convert an integer polynomial (n `i32` values at `a`) into FFT
/// representation (n values at `d`). `a` may overlap `d`; it is most
/// efficient if `a` is at the end of the `d` area.
///
/// # Safety
/// `d` must be valid for `n` `f64` writes, `a` for `n` `i32` reads.
pub unsafe fn fft_int32(logn: u32, d: *mut f64, a: *const i32) {
    let n = 1usize << logn;
    let hn = n >> 1;

    // Move the source to the end of the destination so that it can be
    // modified.
    let cc = (d.add(n) as *mut i32).sub(n);
    if cc as *const i32 != a {
        ptr::copy(a, cc, n);
    }

    // c = (a + b*X^(n/2))/2 with adj(a) = a and adj(b) = -b:
    //   a = c + adj(c),  b = (adj(c) - c)*X^(n/2)
    let aa = (d.add(hn) as *mut i32).sub(hn);
    *aa = (*cc).wrapping_shl(1);
    for j in 1..hn {
        *aa.add(j) = (*cc.add(j)).wrapping_sub(*cc.add(n - j));
    }
    fft_int32_selfadj(logn, d, aa);

    // bb[] is aliased on the second half of cc[].
    let bb = cc.add(hn);
    *bb = (*bb).wrapping_shl(1);
    for j in 1..hn {
        *bb.add(j) = (*bb.add(j)).wrapping_add(*cc.add(hn - j));
    }
    fft_int32_selfadj(logn, d.add(hn), bb);

    // Divisions by 2 were delayed so that the inputs stayed integral.
    for j in 0..n {
        *d.add(j) *= 0.5;
    }
}

pub fn add(logn: u32, a: &mut [f64], b: &[f64]) {
    let n = 1usize << logn;
    for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
        *x += *y;
    }
}

pub fn sub(logn: u32, a: &mut [f64], b: &[f64]) {
    let n = 1usize << logn;
    for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
        *x -= *y;
    }
}

/// a <- a*b (FFT representation).
pub fn mul_fft(logn: u32, a: &mut [f64], b: &[f64]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        let (re, im) = fpc_mul(a[i], a[i + hn], b[i], b[i + hn]);
        a[i] = re;
        a[i + hn] = im;
    }
}

/// a <- x/a (a self-adjoint, half-length).
pub fn invmul_selfadj_fft(logn: u32, a: &mut [f64], x: f64) {
    let hn = 1usize << (logn - 1);
    for v in a[..hn].iter_mut() {
        *v = x / *v;
    }
}

/// a <- a/b (b self-adjoint, half-length).
pub fn div_selfadj_fft(logn: u32, a: &mut [f64], b: &[f64]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        let xb = 1.0 / b[i];
        a[i] *= xb;
        a[i + hn] *= xb;
    }
}

/// d <- a*(X+1)*b, with b self-adjoint representing (X+1)*b. `d` and `a`
/// may be distinct slices only (no overlap).
pub fn mul_ixadj_fft(logn: u32, d: &mut [f64], a: &[f64], b: &[f64]) {
    let hn = 1usize << (logn - 1);
    for i in 0..hn {
        let a0 = a[i];
        let a1 = a[i + hn];
        let b0 = b[i];
        let theta = X_PLUS_1_SQUARED_OVER_X[hn + i];
        let mu = X_MINUS_1_OVER_X_PLUS_1_TIMES_I[hn + i];
        let bt = (b0 * theta) * 0.5;
        d[i] = bt * (a0 + a1 * mu);
        d[i + hn] = bt * (a1 - a0 * mu);
    }
}

/// LDL decomposition step; see the C implementation for the derivation.
pub fn ldl_fft(logn: u32, d: &mut [f64], e: &mut [f64], a0: &[f64], b1: &[f64]) {
    let hn = 1usize << (logn - 1);
    if logn == 1 {
        let a0z = a0[0];
        let b1z = b1[0];
        let ez = (b1z / a0z) * 0.5;
        e[0] = ez;
        d[0] = a0z - ez * b1z;
        return;
    }
    for i in 0..hn {
        let a0z = a0[i];
        let b1z = b1[i];
        let xz = X_PLUS_1_SQUARED_OVER_X[hn + i];
        let ez = b1z / (xz * a0z);
        let dz = a0z - ez * b1z;
        e[i] = ez;
        d[i] = dz;
    }
}

/// Split f into f0 (written over the first half of f) and f1.
pub fn split_fft(logn: u32, f: &mut [f64], f1: &mut [f64]) {
    let hn = 1usize << (logn - 1);
    let qn = hn >> 1;
    if logn == 1 {
        f1[0] = f[1];
        return;
    }
    // First pass; output f0 is written over f[2*i].
    for i in 0..qn {
        let j = i << 1;
        let a_re = f[j];
        let a_im = f[j + hn];
        let b_re = f[j + 1];
        let b_im = f[j + 1 + hn];

        let (t_re, u_re) = add_sub(a_re, b_re);
        let (t_im, u_im) = add_sub(a_im, b_im);
        f[j] = t_re * 0.5;
        f[j + hn] = t_im * 0.5;
        let (s_re, s_im) = if i & 1 == 0 {
            (SIN[hn - 1 - (i >> 1)], -SIN[qn + (i >> 1)])
        } else {
            (-SIN[qn + (i >> 1)], -SIN[hn - 1 - (i >> 1)])
        };
        let (m_re, m_im) = fpc_mul(u_re, u_im, s_re, s_im);
        f1[i] = m_re * 0.5;
        f1[i + qn] = m_im * 0.5;
    }
    // Second pass: squeeze holes out of f0.
    for i in 0..hn {
        f[i] = f[i << 1];
    }
}

/// Split a self-adjoint polynomial `a` (n/2 values) into a0 and b1 (n/4
/// values each).
pub fn split_selfadj_fft(logn: u32, a0: &mut [f64], b1: &mut [f64], a: &[f64]) {
    if logn == 1 {
        a0[0] = a[0];
        b1[0] = 0.0;
        return;
    }
    let hn = 1usize << (logn - 1);
    let qn = hn >> 1;
    for i in 0..qn {
        let a0z = a[i << 1];
        let b1z = a[(i << 1) + 1];
        let (t0, t1) = add_sub(a0z, b1z);
        a0[i] = t0 * 0.5;
        b1[i] = t1 * IRNJ[qn + i];
    }
}

/// Merge f0 and f1 into f = f0(x^2) + x*f1(x^2).
pub fn merge_fft(logn: u32, f: &mut [f64], f0: &[f64], f1: &[f64]) {
    let hn = 1usize << (logn - 1);
    let qn = hn >> 1;
    if logn == 1 {
        f[0] = f0[0];
        f[hn] = f1[0];
    }
    for i in 0..qn {
        let a_re = f0[i];
        let a_im = f0[i + qn];
        let b_re = f1[i];
        let b_im = f1[i + qn];
        let (s_re, s_im) = if i & 1 == 0 {
            (SIN[hn - 1 - (i >> 1)], SIN[qn + (i >> 1)])
        } else {
            (-SIN[qn + (i >> 1)], SIN[hn - 1 - (i >> 1)])
        };
        let (b_re, b_im) = fpc_mul(b_re, b_im, s_re, s_im);
        let (x, y) = add_sub(a_re, b_re);
        f[(i << 1) + 0] = x;
        f[(i << 1) + 1] = y;
        let (x, y) = add_sub(a_im, b_im);
        f[(i << 1) + 0 + hn] = x;
        f[(i << 1) + 1 + hn] = y;
    }
}
