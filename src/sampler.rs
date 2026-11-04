//! Gaussian sampling and Fast Fourier sampling.
//!
//! The sampler's centre and standard deviation are secret and dynamic. The
//! sampler state includes a PRNG (a SHAKE256 context) from which the random
//! bytes are drawn.
//!
//! This follows the C implementation's x86 (SSE2) path operation by
//! operation; Rust's `f64` arithmetic gives the same IEEE-754 results.

use crate::fpoly;
use crate::sha3::Shake;
use core::ptr;

pub const FLAG_INIT0: u32 = 1;
pub const FLAG_NOOUT: u32 = 2;

/// 1/q
const FPR_IQ: f64 = f64::from_bits(0x3F1554E39097A782);
/// log(2)
const LOG2: f64 = f64::from_bits(0x3FE62E42FEFA39EF);
/// 1/log(2)
const INV_LOG2: f64 = f64::from_bits(0x3FF71547652B82FE);
/// 1/(2*(1.8205^2))
const INV_2SQRSIGMA0: f64 = f64::from_bits(0x3FC34F8BC183BBC2);

/// 1/sigma, indexed by logn (1 to 10).
static INV_SIGMA: [f64; 11] = [
    0.0,
    f64::from_bits(0x3F7C48EB7E241699),
    f64::from_bits(0x3F7BE50A548CAED8),
    f64::from_bits(0x3F7B852EE09E762B),
    f64::from_bits(0x3F7AFC5ED3CADA35),
    f64::from_bits(0x3F7A7B3B0976B3EC),
    f64::from_bits(0x3F7A011282CA9C97),
    f64::from_bits(0x3F798D49CE5F2735),
    f64::from_bits(0x3F791F57C56ED9ED),
    f64::from_bits(0x3F78B6C2DE64C7C9),
    f64::from_bits(0x3F78531EF6311AE2),
];

/// sigma_min, indexed by logn (1 to 10).
static SIGMA_MIN: [f64; 11] = [
    0.0,
    f64::from_bits(0x3FF1DD380644568B),
    f64::from_bits(0x3FF21D2EDCAD8626),
    f64::from_bits(0x3FF25C46E1AA7C7A),
    f64::from_bits(0x3FF2B95C574AFB25),
    f64::from_bits(0x3FF314ABC7FE22B6),
    f64::from_bits(0x3FF36E4E3475D7C3),
    f64::from_bits(0x3FF3C65A66A1C224),
    f64::from_bits(0x3FF41CE5358CB3A0),
    f64::from_bits(0x3FF47201BF1F7A75),
    f64::from_bits(0x3FF4C5C19990C764),
];

/// Truncating float-to-int conversions with the hardware semantics of
/// `cvttsd2si` (out-of-range and NaN inputs give the "integer indefinite"
/// value), as in the C code; Rust's `as` is saturating and costs extra
/// instructions on this hot path.
#[inline(always)]
fn cvtt_i32(x: f64) -> i32 {
    #[cfg(target_arch = "x86_64")]
    {
        use core::arch::x86_64::{_mm_cvttsd_si32, _mm_set_sd};
        // SAFETY: SSE2 is part of the x86_64 baseline.
        unsafe { _mm_cvttsd_si32(_mm_set_sd(x)) }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        x as i32
    }
}

#[inline(always)]
fn cvtt_i64(x: f64) -> i64 {
    #[cfg(target_arch = "x86_64")]
    {
        use core::arch::x86_64::{_mm_cvttsd_si64, _mm_set_sd};
        // SAFETY: SSE2 is part of the x86_64 baseline.
        unsafe { _mm_cvttsd_si64(_mm_set_sd(x)) }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        x as i64
    }
}

pub struct SamplerState {
    pub pc: Shake,
    /// Extra input / destination for sampled integers.
    pub dz: *mut i16,
    /// Current offset (reversed) in `dz`.
    pub off: u32,
    /// Top-level degree.
    pub logn: u32,
}

const GAUSS0: [[u32; 3]; 18] = [
    [1375468055, 6936092, 9176346],
    [711562636, 1023934, 15582455],
    [289335016, 4826132, 14746371],
    [90749601, 12313548, 10843417],
    [21676598, 5732767, 9419414],
    [3908963, 11171514, 6206010],
    [529006, 11146683, 11891888],
    [53503, 15610582, 11661180],
    [4032, 7095613, 11671091],
    [225, 16701698, 407118],
    [9, 6787688, 1638204],
    [0, 4870085, 8687822],
    [0, 111406, 13946073],
    [0, 1887, 12017202],
    [0, 23, 11452285],
    [0, 0, 3689579],
    [0, 0, 25354],
    [0, 0, 129],
];

/// Sample the half-Gaussian, using `(lo >> 1) | (hi << 63)` as the 79-bit
/// look-up value (the low bit of `lo` is ignored).
#[inline(always)]
fn gaussian0_helper(lo: u64, hi: u32) -> i32 {
    let v0 = ((lo as u32) >> 1) & 0x00FF_FFFF;
    let v1 = ((lo >> 25) as u32) & 0x00FF_FFFF;
    let v2 = ((lo >> 49) as u32) | (hi << 15);
    let mut z = 0i32;
    for row in GAUSS0.iter() {
        let mut cc = v0.wrapping_sub(row[2]) >> 31;
        cc = v1.wrapping_sub(row[1]).wrapping_sub(cc) >> 31;
        cc = v2.wrapping_sub(row[0]).wrapping_sub(cc) >> 31;
        z += cc as i32;
    }
    z
}

const EXPM_COEFFS: [u64; 13] = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
];

#[inline(always)]
fn mulhi(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) >> 64) as u64
}

/// ccs*exp(-x)*2^63, rounded; requires 0 <= x < log(2) and 0 <= ccs <= 1.
#[inline(always)]
fn expm_p63(x: f64, ccs: f64) -> u64 {
    const TWO63: f64 = 9223372036854775808.0;
    let mut y = EXPM_COEFFS[0];
    let z = (cvtt_i64(x * TWO63) as u64) << 1;
    let w = (cvtt_i64(ccs * TWO63) as u64) << 1;
    for &c in &EXPM_COEFFS[1..] {
        y = c.wrapping_sub(mulhi(z, y));
    }
    mulhi(w, y)
}

/// Sample a bit with probability ccs*exp(-x) (for x >= 0).
#[inline(always)]
fn ber_exp(pc: &mut Shake, x: f64, ccs: f64) -> bool {
    // x = s*log(2) + r, with s an integer and 0 <= r < log(2).
    let si = cvtt_i32(x * INV_LOG2);
    let r = x - (si as f64) * LOG2;

    // Saturate s at 63.
    let mut s = si as u32;
    s |= 63u32.wrapping_sub(s) >> 26;
    s &= 63;

    let z = (expm_p63(r, ccs) << 1).wrapping_sub(1) >> s;

    // Lazy comparison of z with a uniform 64-bit integer.
    let mut i = 56i32;
    while i >= 0 {
        let w = pc.next_u8() as u32;
        let bz = ((z >> i) as u32) & 0xFF;
        if w != bz {
            return w < bz;
        }
        i -= 8;
    }
    false
}

/// Sample the next small integer, with centre `mu` and inverse standard
/// deviation `isigma`.
#[inline(always)]
pub fn sampler_next(ss: &mut SamplerState, mu: f64, isigma: f64) -> i32 {
    // Split centre mu into s + r, for an integer s, and 0 <= r < 1.
    let mut s = cvtt_i32(mu);
    if mu < (s as f64) {
        s -= 1;
    }
    let r = mu - (s as f64);

    // dss = 1/(2*sigma^2) = 0.5*(isigma^2)
    let dss = (isigma * isigma) * 0.5;

    // ccs = sigma_min / sigma = sigma_min * isigma
    let ccs = isigma * SIGMA_MIN[ss.logn as usize];

    loop {
        // Half-Gaussian sample, then a random bit to make it bimodal.
        let lo = ss.pc.next_u64();
        let hi = ss.pc.next_u16();
        let b = (lo as u32 & 1) as i32;
        let z0 = gaussian0_helper(lo, hi);
        let z = b + ((b << 1) - 1) * z0;

        // Rejection sampling.
        let d = (z as f64) - r;
        let mut x = (d * d) * dss;
        x -= ((z0 * z0) as f64) * INV_2SQRSIGMA0;
        if ber_exp(&mut ss.pc, x, ccs) {
            return s + z;
        }
    }
}

#[inline(never)]
unsafe fn ffsamp_fft_deepest(ss: &mut SamplerState, flags: u32, tmp: *mut f64) {
    // n = 2: t = t0 + t1*X; FFT is t0(i) = t0 + t1*i (already split).
    // a = a0 (since a is self-adjoint).
    let a0 = *tmp.add(2);
    let leaf = a0.sqrt() * INV_SIGMA[ss.logn as usize];

    // Elements of dz[] are explored by pairs; base offset is the
    // bit-reversal of a decrementing counter.
    ss.off -= 1;
    let mut j = ss.off;
    j = ((j & 0x1C0) >> 6) | (j & 0x038) | ((j & 0x007) << 6);
    j = ((j & 0x124) >> 2) | (j & 0x092) | ((j & 0x049) << 2);
    j >>= 10 - ss.logn;
    let j = j as usize;
    let j2 = j + (1usize << (ss.logn - 1));

    let mut t0 = (*ss.dz.add(j) as f64) * FPR_IQ;
    if flags & FLAG_INIT0 == 0 {
        t0 += *tmp;
    }
    let z0 = sampler_next(ss, t0, leaf);
    *ss.dz.add(j) = z0 as i16;
    if flags & FLAG_NOOUT == 0 {
        *tmp = t0 - (z0 as f64);
    }

    let mut t1 = (*ss.dz.add(j2) as f64) * FPR_IQ;
    if flags & FLAG_INIT0 == 0 {
        t1 += *tmp.add(1);
    }
    let z1 = sampler_next(ss, t1, leaf);
    *ss.dz.add(j2) = z1 as i16;
    if flags & FLAG_NOOUT == 0 {
        *tmp.add(1) = t1 - (z1 as f64);
    }
}

#[inline(always)]
unsafe fn sm<'a>(p: *mut f64, n: usize) -> &'a mut [f64] {
    core::slice::from_raw_parts_mut(p, n)
}

#[inline(always)]
unsafe fn sc<'a>(p: *const f64, n: usize) -> &'a [f64] {
    core::slice::from_raw_parts(p, n)
}

/// `tmp` contents on input: t (FFT, 8*n bytes; implicit if INIT0), a (FFT,
/// self-adjoint, 4*n bytes). On output: u = t - z (FFT, 8*n bytes) unless
/// NOOUT is set.
unsafe fn ffsamp_fft_inner(ss: &mut SamplerState, logn: u32, flags: u32, tmp: *mut f64) {
    if logn == 1 {
        ffsamp_fft_deepest(ss, flags, tmp);
        return;
    }

    let n = 1usize << logn;
    let hn = n >> 1;
    let qn = n >> 2;

    let tmp0 = tmp;
    let tmp4 = tmp0.add(hn);
    let tmp6 = tmp4.add(qn);
    let tmp8 = tmp4.add(hn);
    let tmp12 = tmp8.add(hn);
    let tmp14 = tmp12.add(qn);

    // (tmp12, tmp14) <- (a0, a1) = split(a)  (self-adj)
    fpoly::split_selfadj_fft(logn, sm(tmp12, qn), sm(tmp14, qn), sc(tmp8, hn));

    // (tmp0, tmp8) <- (t0, t1) = split(t)
    if flags & FLAG_INIT0 == 0 {
        fpoly::split_fft(logn, sm(tmp0, n), sm(tmp8, hn));
    }
    ptr::copy(tmp12, tmp4, qn);

    // tmp6 <- (X+1)*L10;  tmp12 <- a0 - a1*adj(a1)/a0
    fpoly::ldl_fft(logn - 1, sm(tmp12, qn), sm(tmp6, qn), sc(tmp4, qn), sc(tmp14, qn));

    // First recursive call.
    ffsamp_fft_inner(ss, logn - 1, flags & !FLAG_NOOUT, tmp8);

    // tmp12 <- a0
    ptr::copy(tmp4, tmp12, hn);
    // tmp4 <- u1
    ptr::copy(tmp8, tmp4, hn);
    // tmp8 <- t0
    if flags & FLAG_INIT0 == 0 {
        ptr::copy(tmp0, tmp8, hn);
    }

    // tmp0 <- t0' - t0 = u1*L10
    fpoly::mul_ixadj_fft(logn - 1, sm(tmp0, hn), sc(tmp4, hn), sc(tmp14, qn));

    // tmp8 <- t0'
    if flags & FLAG_INIT0 != 0 {
        ptr::copy(tmp0, tmp8, hn);
    } else {
        fpoly::add(logn - 1, sm(tmp8, hn), sc(tmp0, hn));
    }

    // Second recursive call.
    ffsamp_fft_inner(ss, logn - 1, flags & !FLAG_INIT0, tmp8);

    if flags & FLAG_NOOUT == 0 {
        // tmp8 <- u0 = v0 - (t0' - t0)
        fpoly::sub(logn - 1, sm(tmp8, hn), sc(tmp0, hn));
        // tmp12 <- u1
        ptr::copy(tmp4, tmp12, hn);
        // tmp0 <- u = merge(u0, u1)
        fpoly::merge_fft(logn, sm(tmp0, n), sc(tmp8, hn), sc(tmp12, hn));
    }
}

/// Apply Fast Fourier sampling. `tmp` has 16*n bytes: t (8*n bytes,
/// implicit if INIT0), a (4*n bytes, self-adjoint), free (4*n bytes).
///
/// # Safety
/// `tmp` must point to at least `2^(logn+1)` valid, 8-byte aligned `f64`
/// slots, and `ss.dz` to at least `2^logn` valid `i16` slots.
pub unsafe fn ffsamp_fft(ss: &mut SamplerState, flags: u32, tmp: *mut f64) {
    ss.off = 1u32 << (ss.logn - 1);
    let logn = ss.logn;
    ffsamp_fft_inner(ss, logn, flags, tmp);
}
