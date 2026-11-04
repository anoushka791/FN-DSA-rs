//! Computations modulo q = 12289.
//!
//! Polynomials are slices of `u16`. Internal representation: zero is
//! represented by q (values are in [1, q]). Two kernel sets are provided: a
//! portable one (plain scalar code, left to the compiler's autovectorizer)
//! and an AVX2 one, selected at runtime.

use crate::tables::{GM, IGM};

pub const Q: u32 = 12289;
/// -1/q mod 2^32
const Q1I: u32 = 4143984639;
/// 2^64 mod q
pub(crate) const R2: u32 = 5664;
/// Maximum L-infinity norm for the signature elements.
pub const B_INF: u32 = 840;

// ---------------------------------------------------------------------
// Scalar primitives (all operands are u32 holding 16-bit-range values).

#[inline(always)]
pub(crate) fn mq_add(x: u32, y: u32) -> u32 {
    // -(x+y) in [-q, +q-2]
    let mut x = Q.wrapping_sub(x.wrapping_add(y));
    // normalize to [0, q-1]
    x = x.wrapping_add(Q & (x >> 16));
    // negate to get x+y in [1, q]
    Q.wrapping_sub(x)
}

#[inline(always)]
pub(crate) fn mq_sub(x: u32, y: u32) -> u32 {
    let mut y = y.wrapping_sub(x);
    y = y.wrapping_add(Q & (y >> 16));
    Q.wrapping_sub(y)
}

#[inline(always)]
pub(crate) fn mq_half(x: u32) -> u32 {
    x.wrapping_add(Q & (x & 1).wrapping_neg()) >> 1
}

#[inline(always)]
pub(crate) fn mq_mred(x: u32) -> u32 {
    let x = x.wrapping_mul(Q1I);
    let x = (x >> 16) * Q;
    (x >> 16) + 1
}

#[inline(always)]
pub(crate) fn mq_mmul(x: u32, y: u32) -> u32 {
    mq_mred(x.wrapping_mul(y))
}

// ---------------------------------------------------------------------
// Portable kernels.

mod portable {
    use super::*;

    pub fn signed_to_int(logn: u32, d: &mut [u16]) {
        let n = 1usize << logn;
        for v in d[..n].iter_mut() {
            let x = (*v as i16 as i32 as u32).wrapping_neg();
            *v = Q.wrapping_sub(x.wrapping_add(Q & (x >> 16))) as u16;
        }
    }

    pub fn ext_to_int(logn: u32, d: &mut [u16]) {
        let n = 1usize << logn;
        for v in d[..n].iter_mut() {
            let x = *v as u32;
            *v = x.wrapping_add(Q & (x.wrapping_sub(1) >> 16)) as u16;
        }
    }

    pub fn int_to_ext(logn: u32, d: &mut [u16]) {
        let n = 1usize << logn;
        for v in d[..n].iter_mut() {
            let mut x = (*v as u32).wrapping_sub(Q);
            x = x.wrapping_add(Q & (x >> 16));
            *v = x as u16;
        }
    }

    pub fn int_to_ntt(logn: u32, d: &mut [u16]) {
        if logn == 0 {
            return;
        }
        let n = 1usize << logn;
        let d = &mut d[..n];
        let mut t = n;
        for lm in 0..logn {
            let m = 1usize << lm;
            let ht = t >> 1;
            for i in 0..m {
                let s = GM.0[i + m] as u32;
                let (lo, hi) = d[i * t..(i + 1) * t].split_at_mut(ht);
                for (a, b) in lo.iter_mut().zip(hi.iter_mut()) {
                    let x1 = *a as u32;
                    let x2 = mq_mmul(*b as u32, s);
                    *a = mq_add(x1, x2) as u16;
                    *b = mq_sub(x1, x2) as u16;
                }
            }
            t = ht;
        }
    }

    pub fn ntt_to_int(logn: u32, d: &mut [u16]) {
        if logn == 0 {
            return;
        }
        let n = 1usize << logn;
        let d = &mut d[..n];
        let mut t = 1usize;
        for lm in 0..logn {
            let hm = 1usize << (logn - 1 - lm);
            let dt = t << 1;
            for i in 0..hm {
                let s = IGM.0[i + hm] as u32;
                let (lo, hi) = d[i * dt..(i + 1) * dt].split_at_mut(t);
                for (a, b) in lo.iter_mut().zip(hi.iter_mut()) {
                    let x1 = *a as u32;
                    let x2 = *b as u32;
                    *a = mq_half(mq_add(x1, x2)) as u16;
                    *b = mq_mmul(mq_sub(x1, x2), s) as u16;
                }
            }
            t = dt;
        }
    }

    pub fn mul_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
        let n = 1usize << logn;
        for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
            *x = mq_mmul(mq_mmul(*x as u32, *y as u32), R2) as u16;
        }
    }

    pub fn sub(logn: u32, a: &mut [u16], b: &[u16]) {
        let n = 1usize << logn;
        for (x, y) in a[..n].iter_mut().zip(b[..n].iter()) {
            *x = mq_sub(*x as u32, *y as u32) as u16;
        }
    }

    pub fn sqnorm_signed(logn: u32, a: &[u16]) -> u32 {
        let n = 1usize << logn;
        let mut s = 0u32;
        for &v in a[..n].iter() {
            let y = v as i16 as i32;
            s = s.wrapping_add(y.wrapping_mul(y) as u32);
        }
        s
    }

    /// Squared norm of a polynomial in external/internal representation
    /// (normalized to [-q/2, q/2]); saturates to 2^32-1 on overflow or if
    /// any coefficient exceeds B_INF.
    pub fn sqnorm_binf_ext(logn: u32, a: &[u16]) -> u32 {
        let n = 1usize << logn;
        let mut s = 0u32;
        let mut sat = 0u32;
        for &v in a[..n].iter() {
            let mut x = v as u32;
            x = x.wrapping_sub(Q & (((Q - 1) >> 1).wrapping_sub(x) >> 16));
            let y = x as i32;
            s = s.wrapping_add(y.wrapping_mul(y) as u32);
            sat |= s;
            sat |= (B_INF as i32).wrapping_sub(y) as u32;
            sat |= (B_INF as i32).wrapping_add(y) as u32;
        }
        s | (sat >> 31).wrapping_neg()
    }
}

// ---------------------------------------------------------------------
// AVX2 kernels.

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use super::*;
    use core::arch::x86_64::*;

    /// Parameters of a modulus for the vectorized Montgomery arithmetic.
    pub trait Mod {
        const Q: i16;
        /// Low and high 16-bit halves of -1/q mod 2^32, as signed values.
        const Q1ILO: i16;
        const Q1IHI: i16;
        fn gm() -> &'static [u16; 1024];
        fn igm() -> &'static [u16; 1024];
    }

    /// q = 12289
    pub struct M1;
    impl Mod for M1 {
        const Q: i16 = 12289;
        const Q1ILO: i16 = 12287;
        const Q1IHI: i16 = -2304;
        fn gm() -> &'static [u16; 1024] { &GM.0 }
        fn igm() -> &'static [u16; 1024] { &IGM.0 }
    }

    /// q' = 18433
    pub struct M2;
    impl Mod for M2 {
        const Q: i16 = 18433;
        const Q1ILO: i16 = 18431;
        const Q1IHI: i16 = -5184;
        fn gm() -> &'static [u16; 1024] { &crate::sign_tables::GM2.0 }
        fn igm() -> &'static [u16; 1024] { &crate::sign_tables::IGM2.0 }
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn add16<M: Mod>(x: __m256i, y: __m256i) -> __m256i {
        let qq = _mm256_set1_epi16(M::Q);
        let a = _mm256_sub_epi16(qq, _mm256_add_epi16(x, y));
        let b = _mm256_add_epi16(a, _mm256_and_si256(qq, _mm256_srai_epi16::<15>(a)));
        _mm256_sub_epi16(qq, b)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn sub16<M: Mod>(x: __m256i, y: __m256i) -> __m256i {
        let qq = _mm256_set1_epi16(M::Q);
        let a = _mm256_sub_epi16(y, x);
        let b = _mm256_add_epi16(a, _mm256_and_si256(qq, _mm256_srai_epi16::<15>(a)));
        _mm256_sub_epi16(qq, b)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn half16<M: Mod>(x: __m256i) -> __m256i {
        let qq = _mm256_set1_epi16(M::Q);
        let mut y = _mm256_and_si256(x, _mm256_set1_epi16(1));
        y = _mm256_sub_epi16(_mm256_setzero_si256(), y);
        y = _mm256_and_si256(y, qq);
        y = _mm256_add_epi16(x, y);
        _mm256_srli_epi16::<1>(y)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn mred16<M: Mod>(lo: __m256i, hi: __m256i) -> __m256i {
        let qx16 = _mm256_set1_epi16(M::Q);
        let q1ilo = _mm256_set1_epi16(M::Q1ILO);
        let q1ihi = _mm256_set1_epi16(M::Q1IHI);
        let x = _mm256_add_epi16(
            _mm256_add_epi16(
                _mm256_mulhi_epu16(lo, q1ilo),
                _mm256_mullo_epi16(lo, q1ihi),
            ),
            _mm256_mullo_epi16(hi, q1ilo),
        );
        let x = _mm256_mulhi_epu16(x, qx16);
        _mm256_add_epi16(x, _mm256_set1_epi16(1))
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn mmul16<M: Mod>(x: __m256i, y: __m256i) -> __m256i {
        mred16::<M>(_mm256_mullo_epi16(x, y), _mm256_mulhi_epu16(x, y))
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn ld(p: *const u16) -> __m256i {
        _mm256_loadu_si256(p as *const __m256i)
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn st(p: *mut u16, v: __m256i) {
        _mm256_storeu_si256(p as *mut __m256i, v)
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn signed_to_int(logn: u32, d: &mut [u16]) {
        if logn >= 4 {
            let n = 1usize << logn;
            let d = &mut d[..n];
            let qq = _mm256_set1_epi16(Q as i16);
            let p = d.as_mut_ptr();
            for i in 0..(n >> 4) {
                let mut y = ld(p.add(i << 4));
                y = _mm256_sub_epi16(_mm256_setzero_si256(), y);
                y = _mm256_add_epi16(y, _mm256_and_si256(qq, _mm256_srai_epi16::<15>(y)));
                y = _mm256_sub_epi16(qq, y);
                st(p.add(i << 4), y);
            }
        } else {
            portable::signed_to_int(logn, d);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn ext_to_int(logn: u32, d: &mut [u16]) {
        if logn >= 4 {
            let n = 1usize << logn;
            let d = &mut d[..n];
            let qq = _mm256_set1_epi16(Q as i16);
            let p = d.as_mut_ptr();
            for i in 0..(n >> 4) {
                let mut y = ld(p.add(i << 4));
                y = _mm256_add_epi16(
                    y,
                    _mm256_and_si256(qq, _mm256_cmpeq_epi16(y, _mm256_setzero_si256())),
                );
                st(p.add(i << 4), y);
            }
        } else {
            portable::ext_to_int(logn, d);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn int_to_ext(logn: u32, d: &mut [u16]) {
        if logn >= 4 {
            let n = 1usize << logn;
            let d = &mut d[..n];
            let qq = _mm256_set1_epi16(Q as i16);
            let p = d.as_mut_ptr();
            for i in 0..(n >> 4) {
                let y = ld(p.add(i << 4));
                st(p.add(i << 4), _mm256_andnot_si256(_mm256_cmpeq_epi16(y, qq), y));
            }
        } else {
            portable::int_to_ext(logn, d);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn sub(logn: u32, a: &mut [u16], b: &[u16]) {
        if logn >= 4 {
            let n = 1usize << logn;
            let a = &mut a[..n];
            let b = &b[..n];
            let ap = a.as_mut_ptr();
            let bp = b.as_ptr();
            for i in 0..(n >> 4) {
                let ya = ld(ap.add(i << 4));
                let yb = ld(bp.add(i << 4));
                st(ap.add(i << 4), sub16::<M1>(ya, yb));
            }
        } else {
            portable::sub(logn, a, b);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn mul_ntt(logn: u32, a: &mut [u16], b: &[u16]) {
        if logn >= 4 {
            let n = 1usize << logn;
            let a = &mut a[..n];
            let b = &b[..n];
            let ap = a.as_mut_ptr();
            let bp = b.as_ptr();
            let yr2 = _mm256_set1_epi16(R2 as i16);
            for i in 0..(n >> 4) {
                let ya = ld(ap.add(i << 4));
                let yb = ld(bp.add(i << 4));
                st(ap.add(i << 4), mmul16::<M1>(mmul16::<M1>(ya, yb), yr2));
            }
        } else {
            portable::mul_ntt(logn, a, b);
        }
    }

    // Forward NTT of 32 consecutive values (the last 5 levels).
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn ntt32<M: Mod>(a0: &mut __m256i, a1: &mut __m256i, k: usize) {
        let gm = M::gm().as_ptr();
        let mut ya0 = *a0;
        let mut ya1 = *a1;

        // t = 32, m = 1
        let mut yt1 = ya0;
        let mut yt2 = mmul16::<M>(ya1, _mm256_set1_epi16(*gm.add(k) as i16));
        ya0 = add16::<M>(yt1, yt2);
        ya1 = sub16::<M>(yt1, yt2);

        // t = 16, m = 2
        yt1 = _mm256_permute2x128_si256::<0x20>(ya0, ya1);
        yt2 = _mm256_permute2x128_si256::<0x31>(ya0, ya1);
        let g1_0 = *gm.add((k << 1) + 0) as i16;
        let g1_1 = *gm.add((k << 1) + 1) as i16;
        let yg1 = _mm256_setr_epi16(
            g1_0, g1_0, g1_0, g1_0, g1_0, g1_0, g1_0, g1_0, g1_1, g1_1, g1_1, g1_1, g1_1, g1_1,
            g1_1, g1_1,
        );
        yt2 = mmul16::<M>(yt2, yg1);
        ya0 = add16::<M>(yt1, yt2);
        ya1 = sub16::<M>(yt1, yt2);

        // t = 8, m = 4
        yt1 = _mm256_unpacklo_epi64(ya0, ya1);
        yt2 = _mm256_unpackhi_epi64(ya0, ya1);
        let mut yg2 = _mm256_setr_epi64x(
            *gm.add((k << 2) + 0) as i64,
            *gm.add((k << 2) + 1) as i64,
            *gm.add((k << 2) + 2) as i64,
            *gm.add((k << 2) + 3) as i64,
        );
        yg2 = _mm256_or_si256(yg2, _mm256_slli_epi64::<32>(yg2));
        yg2 = _mm256_or_si256(yg2, _mm256_slli_epi32::<16>(yg2));
        yt2 = mmul16::<M>(yt2, yg2);
        ya0 = add16::<M>(yt1, yt2);
        ya1 = sub16::<M>(yt1, yt2);

        // t = 4, m = 8
        let yt3 = _mm256_shuffle_epi32::<0xD8>(ya0);
        let yt4 = _mm256_shuffle_epi32::<0xD8>(ya1);
        yt1 = _mm256_unpacklo_epi32(yt3, yt4);
        yt2 = _mm256_unpackhi_epi32(yt3, yt4);
        let mut yg3 = _mm256_cvtepi16_epi32(_mm_loadu_si128(gm.add(k << 3) as *const __m128i));
        yg3 = _mm256_or_si256(yg3, _mm256_slli_epi32::<16>(yg3));
        yt2 = mmul16::<M>(yt2, yg3);
        ya0 = add16::<M>(yt1, yt2);
        ya1 = sub16::<M>(yt1, yt2);

        // t = 2, m = 16
        let ysk = _mm256_setr_epi8(
            0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3,
            6, 7, 10, 11, 14, 15,
        );
        let yt3 = _mm256_shuffle_epi8(ya0, ysk);
        let yt4 = _mm256_shuffle_epi8(ya1, ysk);
        yt1 = _mm256_unpacklo_epi16(yt3, yt4);
        yt2 = _mm256_unpackhi_epi16(yt3, yt4);
        yt2 = mmul16::<M>(yt2, ld(gm.add(k << 4)));
        ya0 = add16::<M>(yt1, yt2);
        ya1 = sub16::<M>(yt1, yt2);

        yt1 = _mm256_unpacklo_epi16(ya0, ya1);
        yt2 = _mm256_unpackhi_epi16(ya0, ya1);
        ya0 = _mm256_permute2x128_si256::<0x20>(yt1, yt2);
        ya1 = _mm256_permute2x128_si256::<0x31>(yt1, yt2);

        *a0 = ya0;
        *a1 = ya1;
    }

    // Inverse NTT of 32 consecutive values (the first 5 levels).
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn intt32<M: Mod>(a0: &mut __m256i, a1: &mut __m256i, k: usize) {
        let igm = M::igm().as_ptr();
        let mut ya0 = *a0;
        let mut ya1 = *a1;

        let mut yt1 = _mm256_permute2x128_si256::<0x20>(ya0, ya1);
        let mut yt2 = _mm256_permute2x128_si256::<0x31>(ya0, ya1);

        let ysk = _mm256_setr_epi8(
            0, 1, 4, 5, 8, 9, 12, 13, 2, 3, 6, 7, 10, 11, 14, 15, 0, 1, 4, 5, 8, 9, 12, 13, 2, 3,
            6, 7, 10, 11, 14, 15,
        );
        let yt3 = _mm256_shuffle_epi8(yt1, ysk);
        let yt4 = _mm256_shuffle_epi8(yt2, ysk);

        yt1 = _mm256_unpacklo_epi64(yt3, yt4);
        yt2 = _mm256_unpackhi_epi64(yt3, yt4);
        ya0 = half16::<M>(add16::<M>(yt1, yt2));
        ya1 = mmul16::<M>(sub16::<M>(yt1, yt2), ld(igm.add(k << 4)));

        yt1 = _mm256_blend_epi16::<0xAA>(ya0, _mm256_slli_epi32::<16>(ya1));
        yt2 = _mm256_blend_epi16::<0xAA>(_mm256_srli_epi32::<16>(ya0), ya1);
        let mut yig3 = _mm256_cvtepi16_epi32(_mm_loadu_si128(igm.add(k << 3) as *const __m128i));
        yig3 = _mm256_or_si256(yig3, _mm256_slli_epi32::<16>(yig3));
        ya0 = half16::<M>(add16::<M>(yt1, yt2));
        ya1 = mmul16::<M>(sub16::<M>(yt1, yt2), yig3);

        yt1 = _mm256_blend_epi16::<0xCC>(ya0, _mm256_slli_epi64::<32>(ya1));
        yt2 = _mm256_blend_epi16::<0xCC>(_mm256_srli_epi64::<32>(ya0), ya1);
        let mut yig2 = _mm256_setr_epi64x(
            *igm.add((k << 2) + 0) as i64,
            *igm.add((k << 2) + 1) as i64,
            *igm.add((k << 2) + 2) as i64,
            *igm.add((k << 2) + 3) as i64,
        );
        yig2 = _mm256_or_si256(yig2, _mm256_slli_epi64::<32>(yig2));
        yig2 = _mm256_or_si256(yig2, _mm256_slli_epi32::<16>(yig2));
        ya0 = half16::<M>(add16::<M>(yt1, yt2));
        ya1 = mmul16::<M>(sub16::<M>(yt1, yt2), yig2);

        yt1 = _mm256_unpacklo_epi64(ya0, ya1);
        yt2 = _mm256_unpackhi_epi64(ya0, ya1);
        let ig1_0 = *igm.add((k << 1) + 0) as i16;
        let ig1_1 = *igm.add((k << 1) + 1) as i16;
        let yig1 = _mm256_setr_epi16(
            ig1_0, ig1_0, ig1_0, ig1_0, ig1_0, ig1_0, ig1_0, ig1_0, ig1_1, ig1_1, ig1_1, ig1_1,
            ig1_1, ig1_1, ig1_1, ig1_1,
        );
        ya0 = half16::<M>(add16::<M>(yt1, yt2));
        ya1 = mmul16::<M>(sub16::<M>(yt1, yt2), yig1);

        yt1 = _mm256_permute2x128_si256::<0x20>(ya0, ya1);
        yt2 = _mm256_permute2x128_si256::<0x31>(ya0, ya1);
        let yig0 = _mm256_set1_epi16(*igm.add(k) as i16);
        ya0 = half16::<M>(add16::<M>(yt1, yt2));
        ya1 = mmul16::<M>(sub16::<M>(yt1, yt2), yig0);

        *a0 = ya0;
        *a1 = ya1;
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn int_to_ntt_g<M: Mod>(logn: u32, d: &mut [u16]) {
        debug_assert!(logn >= 5);
        if logn == 0 {
            return;
        }
        if logn >= 5 {
            let n = 1usize << logn;
            let d = &mut d[..n];
            let dp = d.as_mut_ptr(); // addressed in units of 16 values
            let mut t = n >> 4;
            for lm in 0..(logn - 5) {
                let m = 1usize << lm;
                let ht = t >> 1;
                let mut j0 = 0usize;
                for i in 0..m {
                    let ys = _mm256_set1_epi16(M::gm()[i + m] as i16);
                    for j in 0..ht {
                        let j1 = j0 + j;
                        let j2 = j1 + ht;
                        let y1 = ld(dp.add(j1 << 4));
                        let y2 = mmul16::<M>(ld(dp.add(j2 << 4)), ys);
                        st(dp.add(j1 << 4), add16::<M>(y1, y2));
                        st(dp.add(j2 << 4), sub16::<M>(y1, y2));
                    }
                    j0 += t;
                }
                t = ht;
            }
            let m = n >> 5;
            for i in 0..m {
                let mut ya0 = ld(dp.add(((i << 1) + 0) << 4));
                let mut ya1 = ld(dp.add(((i << 1) + 1) << 4));
                ntt32::<M>(&mut ya0, &mut ya1, i + m);
                st(dp.add(((i << 1) + 0) << 4), ya0);
                st(dp.add(((i << 1) + 1) << 4), ya1);
            }
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn ntt_to_int_g<M: Mod>(logn: u32, d: &mut [u16]) {
        debug_assert!(logn >= 5);
        if logn == 0 {
            return;
        }
        if logn >= 5 {
            let n = 1usize << logn;
            let d = &mut d[..n];
            let dp = d.as_mut_ptr();
            let m = n >> 5;
            for i in 0..m {
                let mut ya0 = ld(dp.add(((i << 1) + 0) << 4));
                let mut ya1 = ld(dp.add(((i << 1) + 1) << 4));
                intt32::<M>(&mut ya0, &mut ya1, i + m);
                st(dp.add(((i << 1) + 0) << 4), ya0);
                st(dp.add(((i << 1) + 1) << 4), ya1);
            }
            let mut t = 2usize;
            for lm in 5..logn {
                let hm = 1usize << (logn - 1 - lm);
                let dt = t << 1;
                let mut j0 = 0usize;
                for i in 0..hm {
                    let ys = _mm256_set1_epi16(M::igm()[i + hm] as i16);
                    for j in 0..t {
                        let j1 = j0 + j;
                        let j2 = j1 + t;
                        let y1 = ld(dp.add(j1 << 4));
                        let y2 = ld(dp.add(j2 << 4));
                        st(dp.add(j1 << 4), half16::<M>(add16::<M>(y1, y2)));
                        st(dp.add(j2 << 4), mmul16::<M>(ys, sub16::<M>(y1, y2)));
                    }
                    j0 += dt;
                }
                t = dt;
            }
        }
    }


    #[target_feature(enable = "avx2")]
    pub unsafe fn int_to_ntt(logn: u32, d: &mut [u16]) {
        if logn >= 5 {
            int_to_ntt_g::<M1>(logn, d);
        } else {
            portable::int_to_ntt(logn, d);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn ntt_to_int(logn: u32, d: &mut [u16]) {
        if logn >= 5 {
            ntt_to_int_g::<M1>(logn, d);
        } else {
            portable::ntt_to_int(logn, d);
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn sqnorm_signed(logn: u32, a: &[u16]) -> u32 {
        if logn >= 4 {
            let n = 1usize << logn;
            let a = &a[..n];
            let ap = a.as_ptr();
            let mut ys = _mm256_setzero_si256();
            for i in 0..(n >> 4) {
                let y = ld(ap.add(i << 4));
                let ylo = _mm256_mullo_epi16(y, y);
                let yhi = _mm256_mulhi_epi16(y, y);
                let y0 = _mm256_blend_epi16::<0xAA>(ylo, _mm256_bslli_epi128::<2>(yhi));
                let y1 = _mm256_blend_epi16::<0xAA>(_mm256_bsrli_epi128::<2>(ylo), yhi);
                ys = _mm256_add_epi32(ys, _mm256_add_epi32(y0, y1));
            }
            ys = _mm256_add_epi32(ys, _mm256_srli_epi64::<32>(ys));
            ys = _mm256_add_epi32(ys, _mm256_bsrli_epi128::<8>(ys));
            let xs = _mm_add_epi32(_mm256_castsi256_si128(ys), _mm256_extracti128_si256::<1>(ys));
            _mm_cvtsi128_si32(xs) as u32
        } else {
            portable::sqnorm_signed(logn, a)
        }
    }

    #[target_feature(enable = "avx2")]
    pub unsafe fn sqnorm_binf_ext(logn: u32, a: &[u16]) -> u32 {
        if logn >= 4 {
            let n = 1usize << logn;
            let a = &a[..n];
            let ap = a.as_ptr();
            let mut ys = _mm256_setzero_si256();
            let mut ysat = _mm256_setzero_si256();
            let qq = _mm256_set1_epi16(Q as i16);
            let hq = _mm256_set1_epi16(((Q - 1) >> 1) as i16);
            let bbp = _mm256_set1_epi16(B_INF as i16);
            let bbm = _mm256_set1_epi16(-(B_INF as i16));
            let mut ylif = _mm256_setzero_si256();
            for i in 0..(n >> 4) {
                let mut y = ld(ap.add(i << 4));

                // Normalize to [-q/2, +q/2].
                let ym = _mm256_cmpgt_epi16(y, hq);
                y = _mm256_sub_epi16(y, _mm256_and_si256(ym, qq));

                // Flag values outside of [-B_INF, +B_INF].
                ylif = _mm256_or_si256(
                    ylif,
                    _mm256_or_si256(_mm256_cmpgt_epi16(y, bbp), _mm256_cmpgt_epi16(bbm, y)),
                );

                // Compute and add coefficient squares.
                let ylo = _mm256_mullo_epi16(y, y);
                let yhi = _mm256_mulhi_epi16(y, y);
                let y0 = _mm256_blend_epi16::<0xAA>(ylo, _mm256_bslli_epi128::<2>(yhi));
                let y1 = _mm256_blend_epi16::<0xAA>(_mm256_bsrli_epi128::<2>(ylo), yhi);
                ys = _mm256_add_epi32(ys, _mm256_add_epi32(y0, y1));

                // Any overflow must first go through a value with its high
                // bit set.
                ysat = _mm256_or_si256(ysat, ys);
            }

            ylif = _mm256_or_si256(ylif, _mm256_slli_epi32::<16>(ylif));
            ysat = _mm256_or_si256(ysat, ylif);

            ys = _mm256_add_epi32(ys, _mm256_srli_epi64::<32>(ys));
            ysat = _mm256_or_si256(ysat, ys);
            ys = _mm256_add_epi32(ys, _mm256_bsrli_epi128::<8>(ys));
            ysat = _mm256_or_si256(ysat, ys);
            let xs = _mm_add_epi32(_mm256_castsi256_si128(ys), _mm256_extracti128_si256::<1>(ys));
            let r = _mm_cvtsi128_si32(xs) as u32;
            let mut sat = _mm256_movemask_epi8(ysat) as u32;
            sat = (sat & 0x8888_8888) | (r & 0x8000_0000);
            sat |= sat.wrapping_neg();
            r | (((sat as i32) >> 31) as u32)
        } else {
            portable::sqnorm_binf_ext(logn, a)
        }
    }
}

// ---------------------------------------------------------------------
// Kernel tables and runtime selection.

/// A set of polynomial kernels. All functions take the degree `logn` and
/// slices with at least `2^logn` elements.
pub struct Kernels {
    pub ext_to_int: fn(u32, &mut [u16]),
    pub int_to_ext: fn(u32, &mut [u16]),
    pub signed_to_int: fn(u32, &mut [u16]),
    pub int_to_ntt: fn(u32, &mut [u16]),
    pub ntt_to_int: fn(u32, &mut [u16]),
    pub mul_ntt: fn(u32, &mut [u16], &[u16]),
    pub sub: fn(u32, &mut [u16], &[u16]),
    pub sqnorm_signed: fn(u32, &[u16]) -> u32,
    pub sqnorm_binf_ext: fn(u32, &[u16]) -> u32,
}

pub static PORTABLE: Kernels = Kernels {
    ext_to_int: portable::ext_to_int,
    int_to_ext: portable::int_to_ext,
    signed_to_int: portable::signed_to_int,
    int_to_ntt: portable::int_to_ntt,
    ntt_to_int: portable::ntt_to_int,
    mul_ntt: portable::mul_ntt,
    sub: portable::sub,
    sqnorm_signed: portable::sqnorm_signed,
    sqnorm_binf_ext: portable::sqnorm_binf_ext,
};

#[cfg(target_arch = "x86_64")]
mod avx2_kernels {
    use super::*;

    // SAFETY (all wrappers): these functions are only reachable through
    // `AVX2`, which `avx2()` hands out only after a successful runtime
    // check for AVX2 support.
    pub static TABLE: Kernels = Kernels {
        ext_to_int: |l, d| unsafe { avx2::ext_to_int(l, d) },
        int_to_ext: |l, d| unsafe { avx2::int_to_ext(l, d) },
        signed_to_int: |l, d| unsafe { avx2::signed_to_int(l, d) },
        int_to_ntt: |l, d| unsafe { avx2::int_to_ntt(l, d) },
        ntt_to_int: |l, d| unsafe { avx2::ntt_to_int(l, d) },
        mul_ntt: |l, a, b| unsafe { avx2::mul_ntt(l, a, b) },
        sub: |l, a, b| unsafe { avx2::sub(l, a, b) },
        sqnorm_signed: |l, a| unsafe { avx2::sqnorm_signed(l, a) },
        sqnorm_binf_ext: |l, a| unsafe { avx2::sqnorm_binf_ext(l, a) },
    };
}

/// Forward NTT modulo q' = 18433 with AVX2, if available and `logn >= 5`.
/// Returns false (leaving `d` untouched) if the caller must use the scalar
/// code instead.
#[inline]
pub fn mq2_ntt_fast(logn: u32, d: &mut [u16]) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        if logn >= 5 && has_avx2() {
            assert!(d.len() >= 1usize << logn);
            // SAFETY: AVX2 support was just checked.
            unsafe { avx2::int_to_ntt_g::<avx2::M2>(logn, d) };
            return true;
        }
    }
    let _ = (logn, d);
    false
}

/// Inverse NTT modulo q' = 18433; see [`mq2_ntt_fast`].
#[inline]
pub fn mq2_intt_fast(logn: u32, d: &mut [u16]) -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        if logn >= 5 && has_avx2() {
            assert!(d.len() >= 1usize << logn);
            // SAFETY: AVX2 support was just checked.
            unsafe { avx2::ntt_to_int_g::<avx2::M2>(logn, d) };
            return true;
        }
    }
    let _ = (logn, d);
    false
}

/// Whether AVX2 is usable on this CPU (and enabled by the OS).
pub fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        std::is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// The AVX2 kernels, if this CPU supports them.
pub fn avx2() -> Option<&'static Kernels> {
    #[cfg(target_arch = "x86_64")]
    {
        if has_avx2() {
            return Some(&avx2_kernels::TABLE);
        }
    }
    None
}

/// The best kernel set for this machine.
#[inline]
pub fn best() -> &'static Kernels {
    avx2().unwrap_or(&PORTABLE)
}

/// Check whether a squared norm of a signature vector (s1, s2) is acceptable.
pub fn sqnorm_is_acceptable(logn: u32, norm: u32) -> bool {
    const SQBETA: [u32; 11] = [
        0, 101498, 208714, 428865, 892039, 1852696, 3842630, 7959734, 16468416, 34034726,
        70265242,
    ];
    norm <= SQBETA[logn as usize]
}
