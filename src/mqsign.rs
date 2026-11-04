//! Arithmetic used only by signing: extra operations modulo q = 12289, and
//! operations modulo q' = 18433 (used to compute f*adj(f) + g*adj(g) exactly
//! by combining the results modulo q and q').
//!
//! These are the plain scalar versions (the C code also has SIMD variants);
//! they are written so that the compiler can vectorize the simple loops.

use crate::mq::{mq_add, mq_mmul, mq_sub, B_INF, Q, R2};
use crate::sign_tables::{GM2, IGM2};

// ---------------------------------------------------------------------
// Modulo q.

pub fn small_to_int(logn: u32, f: &[i8], d: &mut [u16]) {
    let n = 1usize << logn;
    for (dv, &fv) in d[..n].iter_mut().zip(f[..n].iter()) {
        let x = (fv as i32 as u32).wrapping_neg();
        *dv = Q.wrapping_sub(x.wrapping_add(Q & (x >> 16))) as u16;
    }
}

pub fn int_to_signed(logn: u32, d: &mut [u16]) {
    let n = 1usize << logn;
    for v in d[..n].iter_mut() {
        let mut x = *v as u32;
        x = x.wrapping_sub(Q & (((Q - 1) >> 1).wrapping_sub(x) >> 16));
        *v = x as u16;
    }
}

pub fn neg(logn: u32, a: &mut [u16]) {
    let n = 1usize << logn;
    for v in a[..n].iter_mut() {
        *v = mq_sub(Q, *v as u32) as u16;
    }
}

/// a <- a*adj(a) + b*adj(b), all in NTT representation.
pub fn muladj_x2_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    let hn = n >> 1;
    for i in 0..hn {
        a[i] = mq_mmul(
            mq_add(
                mq_mmul(a[i] as u32, a[n - 1 - i] as u32),
                mq_mmul(b[i] as u32, b[n - 1 - i] as u32),
            ),
            R2,
        ) as u16;
    }
    for i in hn..n {
        a[i] = a[n - 1 - i];
    }
}

fn mq_div(x: u32, y: u32) -> u32 {
    let y = mq_mmul(y, R2);
    let y2 = mq_mmul(y, y);
    let y3 = mq_mmul(y2, y);
    let y5 = mq_mmul(y3, y2);
    let y10 = mq_mmul(y5, y5);
    let y20 = mq_mmul(y10, y10);
    let y40 = mq_mmul(y20, y20);
    let y80 = mq_mmul(y40, y40);
    let y160 = mq_mmul(y80, y80);
    let y163 = mq_mmul(y160, y3);
    let y323 = mq_mmul(y163, y160);
    let y646 = mq_mmul(y323, y323);
    let y1292 = mq_mmul(y646, y646);
    let y1455 = mq_mmul(y1292, y163);
    let y2910 = mq_mmul(y1455, y1455);
    let y5820 = mq_mmul(y2910, y2910);
    let y6143 = mq_mmul(y5820, y323);
    let y12286 = mq_mmul(y6143, y6143);
    let iy = mq_mmul(y12286, y);
    mq_mmul(x, iy)
}

/// Invert `a` in place (NTT representation), using `tmp` (n values) as
/// scratch. Returns true if `a` was invertible.
pub fn inv_ntt(logn: u32, a: &mut [u16], tmp: &mut [u16]) -> bool {
    let n = 1usize << logn;
    let mut ax = a[0] as u32;
    tmp[0] = ax as u16;
    for i in 1..n {
        ax = mq_mmul(ax, a[i] as u32);
        tmp[i] = ax as u16;
    }
    ax = mq_div(1, ax);
    let r = ax.wrapping_sub(Q);
    for i in (1..n).rev() {
        let cx = mq_mmul(ax, tmp[i - 1] as u32);
        ax = mq_mmul(ax, a[i] as u32);
        a[i] = cx as u16;
    }
    a[0] = ax as u16;
    (r >> 16) != 0
}

/// a <- a/b (NTT representation), using `tmp` (n values) as scratch.
pub fn div_ntt(logn: u32, a: &mut [u16], b: &[u16], tmp: &mut [u16]) -> bool {
    let n = 1usize << logn;
    let mut bx = b[0] as u32;
    tmp[0] = bx as u16;
    for i in 1..n {
        bx = mq_mmul(bx, b[i] as u32);
        tmp[i] = bx as u16;
    }
    bx = mq_div(1, bx);
    let r = bx.wrapping_sub(Q);
    for i in (1..n).rev() {
        let cx = mq_mmul(bx, tmp[i - 1] as u32);
        bx = mq_mmul(bx, b[i] as u32);
        a[i] = mq_mmul(mq_mmul(a[i] as u32, cx), R2) as u16;
    }
    a[0] = mq_mmul(mq_mmul(a[0] as u32, bx), R2) as u16;
    (r >> 16) != 0
}

/// a <- a*adj(b) + c*adj(d) (NTT representation; no overlap).
pub fn muladj_add_muladj(logn: u32, a: &mut [u16], b: &[u16], c: &[u16], d: &[u16]) {
    let n = 1usize << logn;
    for i in 0..n {
        let x = mq_mmul(a[i] as u32, b[n - 1 - i] as u32);
        let y = mq_mmul(c[i] as u32, d[n - 1 - i] as u32);
        a[i] = mq_mmul(mq_add(x, y), R2) as u16;
    }
}

/// Squared norm of a signed polynomial; saturates to 2^32-1 on overflow or
/// if any coefficient exceeds B_INF in absolute value.
pub fn sqnorm_binf_signed(logn: u32, a: &[i16]) -> u32 {
    let n = 1usize << logn;
    let mut s = 0u32;
    let mut sat = 0u32;
    let binf = B_INF as i32;
    for &v in a[..n].iter() {
        let y = v as i32;
        s = s.wrapping_add(y.wrapping_mul(y) as u32);
        sat |= s;
        sat |= binf.wrapping_sub(y) as u32;
        sat |= binf.wrapping_add(y) as u32;
    }
    s | (sat >> 31).wrapping_neg()
}

pub fn sqnorm_small(logn: u32, a: &[i8]) -> u32 {
    let n = 1usize << logn;
    let mut s = 0u32;
    for &v in a[..n].iter() {
        let y = v as i32;
        s = s.wrapping_add(y.wrapping_mul(y) as u32);
    }
    s
}

// ---------------------------------------------------------------------
// Modulo q' = 18433.

const Q2: u32 = 18433;
/// -1/q' mod 2^32
const Q2_1I: u32 = 3955247103;
/// 2^64 mod q'
const R2_2: u32 = 806;
/// 1/q mod q' (Montgomery representation)
const Q2_IQ_MR: u32 = 13692;
/// 1/q mod q'
pub const MQ2_IQ: u32 = 3;

#[inline(always)]
fn mq2_add(x: u32, y: u32) -> u32 {
    let mut x = Q2.wrapping_sub(x.wrapping_add(y));
    x = x.wrapping_add(Q2 & (x >> 16));
    Q2.wrapping_sub(x)
}

#[inline(always)]
fn mq2_sub(x: u32, y: u32) -> u32 {
    let mut y = y.wrapping_sub(x);
    y = y.wrapping_add(Q2 & (y >> 16));
    Q2.wrapping_sub(y)
}

#[inline(always)]
fn mq2_half(x: u32) -> u32 {
    x.wrapping_add(Q2 & (x & 1).wrapping_neg()) >> 1
}

#[inline(always)]
fn mq2_mred(x: u32) -> u32 {
    let x = x.wrapping_mul(Q2_1I);
    let x = (x >> 16) * Q2;
    (x >> 16) + 1
}

#[inline(always)]
fn mq2_mmul(x: u32, y: u32) -> u32 {
    mq2_mred(x.wrapping_mul(y))
}

pub fn mq2_small_to_mod(logn: u32, f: &[i8], d: &mut [u16]) {
    let n = 1usize << logn;
    for (dv, &fv) in d[..n].iter_mut().zip(f[..n].iter()) {
        let x = (fv as i32 as u32).wrapping_neg();
        *dv = Q2.wrapping_sub(x.wrapping_add(Q2 & (x >> 16))) as u16;
    }
}

pub fn mq2_signed_to_mod(logn: u32, d: &mut [u16]) {
    let n = 1usize << logn;
    for v in d[..n].iter_mut() {
        let x = (*v as i16 as i32 as u32).wrapping_neg();
        *v = Q2.wrapping_sub(x.wrapping_add(Q2 & (x >> 16))) as u16;
    }
}

pub fn mq2_mod_to_signed(logn: u32, d: &mut [u16]) {
    let n = 1usize << logn;
    for v in d[..n].iter_mut() {
        let mut x = *v as u32;
        x = x.wrapping_sub(Q2 & (((Q2 - 1) >> 1).wrapping_sub(x) >> 16));
        *v = x as u16;
    }
}

pub fn mq2_ntt(logn: u32, d: &mut [u16]) {
    if crate::mq::mq2_ntt_fast(logn, d) {
        return;
    }
    let n = 1usize << logn;
    let d = &mut d[..n];
    let mut t = n;
    for lm in 0..logn {
        let m = 1usize << lm;
        let ht = t >> 1;
        for i in 0..m {
            let s = GM2.0[i + m] as u32;
            let (lo, hi) = d[i * t..(i + 1) * t].split_at_mut(ht);
            for (a, b) in lo.iter_mut().zip(hi.iter_mut()) {
                let x1 = *a as u32;
                let x2 = mq2_mmul(*b as u32, s);
                *a = mq2_add(x1, x2) as u16;
                *b = mq2_sub(x1, x2) as u16;
            }
        }
        t = ht;
    }
}

pub fn mq2_intt(logn: u32, d: &mut [u16]) {
    if logn == 0 {
        return;
    }
    if crate::mq::mq2_intt_fast(logn, d) {
        return;
    }
    let n = 1usize << logn;
    let d = &mut d[..n];
    let mut t = 1usize;
    for lm in 0..logn {
        let hm = 1usize << (logn - 1 - lm);
        let dt = t << 1;
        for i in 0..hm {
            let s = IGM2.0[i + hm] as u32;
            let (lo, hi) = d[i * dt..(i + 1) * dt].split_at_mut(t);
            for (a, b) in lo.iter_mut().zip(hi.iter_mut()) {
                let x1 = *a as u32;
                let x2 = *b as u32;
                *a = mq2_half(mq2_add(x1, x2)) as u16;
                *b = mq2_mmul(mq2_sub(x1, x2), s) as u16;
            }
        }
        t = dt;
    }
}

pub fn mq2_mul_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
        *x = mq2_mmul(mq2_mmul(*x as u32, *y as u32), R2_2) as u16;
    }
}

/// a <- a*adj(b)
pub fn mq2_muladj_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    for i in 0..n {
        a[i] = mq2_mmul(mq2_mmul(a[i] as u32, b[n - 1 - i] as u32), R2_2) as u16;
    }
}

/// a <- a*adj(a) + b*adj(b)
pub fn mq2_muladj_x2_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    let hn = n >> 1;
    for i in 0..hn {
        a[i] = mq2_mmul(
            mq2_add(
                mq2_mmul(a[i] as u32, a[n - 1 - i] as u32),
                mq2_mmul(b[i] as u32, b[n - 1 - i] as u32),
            ),
            R2_2,
        ) as u16;
    }
    for i in hn..n {
        a[i] = a[n - 1 - i];
    }
}

pub fn mq2_sub_poly(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
        *x = mq2_sub(*x as u32, *y as u32) as u16;
    }
}

pub fn mq2_add_poly(logn: u32, a: &mut [u16], b: &[u16]) {
    let n = 1usize << logn;
    for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
        *x = mq2_add(*x as u32, *y as u32) as u16;
    }
}

pub fn mq2_mulconst(logn: u32, a: &mut [u16], m: u32) {
    let n = 1usize << logn;
    let m = mq2_mmul(m, R2_2);
    for v in a[..n].iter_mut() {
        *v = mq2_mmul(*v as u32, m) as u16;
    }
}

/// Given `a` (mod q', normalized) and `b` (mod q, normalized) representing
/// the same polynomial, adjust `a` so that it matches `b` as an integer
/// polynomial; adjustment bits are saved into `adj`.
pub fn mq2mq_adjust_save(logn: u32, a: &mut [i16], b: &[i16], adj: &mut [u8]) {
    let n = 1usize << logn;
    if logn >= 3 {
        adj[..n >> 3].fill(0);
    } else {
        adj[0] = 0;
    }
    for i in 0..n {
        let y1 = a[i] as i32;
        let y2 = b[i] as i32;
        let x0 = y1.wrapping_sub(y2) as u32;
        let xp = x0.wrapping_sub(Q);
        let xm = x0.wrapping_add(Q);
        let xa = (x0 | x0.wrapping_neg()) & (xp | xp.wrapping_neg()) & (xm | xm.wrapping_neg());
        adj[i >> 3] |= (((xa >> 16) & 1) << (i & 7)) as u8;
        let s = !((y1 >> 16) as u32);
        let d = (xa >> 16) & ((Q2 ^ s).wrapping_sub(s));
        a[i] = y1.wrapping_add(d as i32) as i16;
    }
}

/// Reduced variant of `mq2mq_adjust_save` (adjustment bits already known).
pub fn mq2mq_adjust(logn: u32, a: &mut [i16], adj: &[u8]) {
    let n = 1usize << logn;
    for i in 0..n {
        let y1 = a[i] as i32;
        let xa = ((adj[i >> 3] >> (i & 7)) & 1) as u32;
        let s = !((y1 >> 16) as u32);
        let d = xa.wrapping_neg() & ((Q2 ^ s).wrapping_sub(s));
        a[i] = y1.wrapping_add(d as i32) as i16;
    }
}

/// Given `aq` (mod q) and `aq2` (mod q'), both normalized, compute the
/// value modulo q*q' (normalized) into `d`.
pub fn mq2mq_crt(logn: u32, d: &mut [i32], aq: &[i16], aq2: &[i16]) {
    let n = 1usize << logn;
    for i in 0..n {
        let xq = aq[i] as i32;
        let xq2 = aq2[i] as i32;
        let mut u = xq2.wrapping_sub(xq) as u32;
        let w = u.wrapping_neg();
        u = Q2.wrapping_sub(w.wrapping_add(Q2 & (w >> 16)));
        u = mq2_mmul(u, Q2_IQ_MR);
        u = (xq as u32).wrapping_add(Q.wrapping_mul(u));
        let v = ((Q * Q2) >> 1).wrapping_sub(u);
        u = u.wrapping_sub((((v as i32) >> 30) as u32) & (Q * Q2));
        d[i] = u as i32;
    }
}
