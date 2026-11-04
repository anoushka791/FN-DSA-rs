//! Solving the NTRU equation f*G - g*F = q (mod X^n+1) for key generation.
//!
//! Port of kgen_ntru.c. The recursion works in one shared temporary area of
//! `u32` words, split into overlapping views exactly as the C code does, so
//! the code is pointer-based.

use crate::fxp::*;
use crate::kpoly::*;
use crate::mp31::*;
use crate::zint31::*;
use core::ptr;
use core::slice::{from_raw_parts as sc, from_raw_parts_mut as sm};

const Q: u32 = 12289;

static MAX_BL_SMALL: [usize; 11] = [1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207];
static MAX_BL_LARGE: [usize; 10] = [1, 2, 3, 6, 11, 21, 40, 78, 155, 308];
static WORD_WIN: [usize; 10] = [1, 1, 2, 2, 2, 3, 3, 4, 5, 7];
static MIN_SAVE_FG: [usize; 11] = [0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 4];

/// Smallest degree for which the NTT variant of the F reduction is used.
const MIN_LOGN_FGNTT: u32 = 4;

/// Offset (in words) in tmp[] for saving intermediate (f,g) values.
#[inline(always)]
fn fg_save_offset(logn_top: u32) -> usize {
    5usize << logn_top
}

/// Convert source f and g into RNS+NTT at the start of tmp (one word per
/// coefficient).
unsafe fn make_fg_zero(logn: u32, f: &[i8], g: &[i8], tmp: *mut u32) {
    let n = 1usize << logn;
    let ft = tmp;
    let gt = ft.add(n);
    let gm = gt.add(n);
    let p = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    poly_mp_set_small(logn, sm(ft, n), f, p);
    poly_mp_set_small(logn, sm(gt, n), g, p);
    mp_mkgm(logn, sm(gm, n), PRIMES[0].g, p, p0i);
    mp_ntt(logn, sm(ft, n), sc(gm, n), p, p0i);
    mp_ntt(logn, sm(gt, n), sc(gm, n), p, p0i);
}

/// One step of computing (f,g) at a given depth: input (f,g) of degree
/// 2^(logn_top-depth), output (f',g') of degree 2^(logn_top-(depth+1)), both
/// at the start of tmp, in RNS+NTT notation.
unsafe fn make_fg_step(logn_top: u32, depth: u32, tmp: *mut u32) {
    let logn = logn_top - depth;
    let n = 1usize << logn;
    let hn = n >> 1;
    let slen = MAX_BL_SMALL[depth as usize];
    let tlen = MAX_BL_SMALL[depth as usize + 1];

    // Layout: fd (hn*tlen), gd (hn*tlen), fs (n*slen), gs (n*slen),
    // t1 (n), t2 (extra).
    let fd = tmp;
    let gd = fd.add(hn * tlen);
    let fs = gd.add(hn * tlen);
    let gs = fs.add(n * slen);
    let t1 = gs.add(n * slen);
    let t2 = t1.add(n);
    ptr::copy(tmp, fs, 2 * n * slen);

    // First slen words: use the input values directly, applying the inverse
    // NTT as we go to get the sources in RNS (non-NTT).
    let mut xf = fs;
    let mut xg = gs;
    let mut yf = fd;
    let mut yg = gd;
    for i in 0..slen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        for j in 0..hn {
            *yf.add(j) = mp_mmul(mp_mmul(*xf.add(2 * j), *xf.add(2 * j + 1), p, p0i), r2, p, p0i);
            *yg.add(j) = mp_mmul(mp_mmul(*xg.add(2 * j), *xg.add(2 * j + 1), p, p0i), r2, p, p0i);
        }
        mp_mkigm(logn, sm(t1, n), PRIMES[i].ig, p, p0i);
        mp_intt(logn, sm(xf, n), sc(t1, n), p, p0i);
        mp_intt(logn, sm(xg, n), sc(t1, n), p, p0i);
        xf = xf.add(n);
        xg = xg.add(n);
        yf = yf.add(hn);
        yg = yg.add(hn);
    }

    // Rebuild the plain integer coefficients of fs and gs.
    zint_rebuild_crt(fs, slen, n, 2, true, t1);

    // Remaining output words.
    for i in slen..tlen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        let rx = mp_rx31(slen as u32, p, p0i, r2);
        mp_mkgm(logn, sm(t1, n), PRIMES[i].g, p, p0i);
        for j in 0..n {
            *t2.add(j) = zint_mod_small_signed(fs.add(j), slen, n, p, p0i, r2, rx);
        }
        mp_ntt(logn, sm(t2, n), sc(t1, n), p, p0i);
        for j in 0..hn {
            *yf.add(j) = mp_mmul(mp_mmul(*t2.add(2 * j), *t2.add(2 * j + 1), p, p0i), r2, p, p0i);
        }
        yf = yf.add(hn);
        for j in 0..n {
            *t2.add(j) = zint_mod_small_signed(gs.add(j), slen, n, p, p0i, r2, rx);
        }
        mp_ntt(logn, sm(t2, n), sc(t1, n), p, p0i);
        for j in 0..hn {
            *yg.add(j) = mp_mmul(mp_mmul(*t2.add(2 * j), *t2.add(2 * j + 1), p, p0i), r2, p, p0i);
        }
        yg = yg.add(hn);
    }
}

/// Compute (f,g) at a depth below the top, in RNS+NTT, at the start of tmp.
unsafe fn make_fg_intermediate(logn_top: u32, f: &[i8], g: &[i8], depth: u32, tmp: *mut u32) {
    make_fg_zero(logn_top, f, g, tmp);
    for d in 0..depth {
        make_fg_step(logn_top, d, tmp);
    }
}

/// Compute (f,g) at the deepest level (the resultants). Intermediate values
/// below the save threshold are copied at the end of tmp.
unsafe fn make_fg_deepest(logn_top: u32, f: &[i8], g: &[i8], tmp: *mut u32, mut sav_off: usize) {
    make_fg_zero(logn_top, f, g, tmp);
    for d in 0..logn_top {
        make_fg_step(logn_top, d, tmp);

        // Save (f,g) for depth d+1 if it is at least at the save threshold
        // but not the deepest level.
        let d2 = d + 1;
        if d2 < logn_top && d2 as usize >= MIN_SAVE_FG[logn_top as usize] {
            let slen = MAX_BL_SMALL[d2 as usize];
            let fglen = slen << (logn_top + 1 - d2);
            sav_off -= fglen;
            ptr::copy(tmp, tmp.add(sav_off), fglen);
        }
    }
}

/// Solve the NTRU equation at the deepest level; F is written into tmp.
unsafe fn solve_ntru_deepest(logn_top: u32, f: &[i8], g: &[i8], tmp: *mut u32) -> bool {
    make_fg_deepest(logn_top, f, g, tmp, fg_save_offset(logn_top));

    // Reorganize memory: Fp (len), Gp (len), fp (len), gp (len), t1.
    let len = MAX_BL_SMALL[logn_top as usize];
    let fp_big = tmp;
    let gp_big = fp_big.add(len);
    let fp = gp_big.add(len);
    let gp = fp.add(len);
    let t1 = gp.add(len);
    ptr::copy(tmp, fp, 2 * len);

    // Convert the resultants back into plain integers.
    zint_rebuild_crt(fp, len, 1, 2, false, t1);

    // Binary GCD: f*G - g*F = 1.
    if !zint_bezout(gp_big, fp_big, fp, gp, len, t1) {
        return false;
    }

    // Multiply by q (only F; G is discarded).
    zint_mul_small(fp_big, len, Q) == 0
}

/// Solve the NTRU equation at an intermediate level: input (F,G) from one
/// level deeper (plain) at the start of tmp, output (F,G) for this level.
unsafe fn solve_ntru_intermediate(
    logn_top: u32,
    f: &[i8],
    g: &[i8],
    depth: u32,
    tmp: *mut u32,
) -> bool {
    let logn = logn_top - depth;
    let n = 1usize << logn;
    let hn = n >> 1;

    // slen: size for (f,g) at this level (also size of output F);
    // llen: size for unreduced F at this level; dlen: size for F from deeper.
    let slen = MAX_BL_SMALL[depth as usize];
    let llen = MAX_BL_LARGE[depth as usize];
    let dlen = MAX_BL_SMALL[depth as usize + 1];

    // Fd: F from deeper level (dlen*hn); fgt: f,g from this level (2*slen*n).
    let mut fd = tmp;
    let fgt = fd.add(dlen * hn);

    // Get (f,g) for this level (in RNS+NTT).
    if (depth as usize) < MIN_SAVE_FG[logn_top as usize] {
        make_fg_intermediate(logn_top, f, g, depth, fgt);
    } else {
        let mut sav_fg = tmp.add(fg_save_offset(logn_top));
        for d in MIN_SAVE_FG[logn_top as usize]..=(depth as usize) {
            sav_fg = sav_fg.sub(MAX_BL_SMALL[d] << (logn_top + 1 - d as u32));
        }
        ptr::copy(sav_fg, fgt, 2 * slen * n);
    }

    // Move buffers: Ft (llen*n), ft (slen*n), gt (slen*n), Fd (dlen*hn).
    let ft_big = tmp;
    let ft = ft_big.add(llen * n);
    let gt = ft.add(slen * n);
    fd = gt.add(slen * n);
    ptr::copy(fgt, ft, 2 * n * slen);
    ptr::copy(tmp, fd, hn * dlen);

    // Convert Fd to RNS, with output temporarily stored in Ft: the values
    // for each modulus p go in the last hn slots of the n-word line.
    for i in 0..llen {
        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;
        let rx = mp_rx31(dlen as u32, p, p0i, r2);
        let xt = ft_big.add(i * n + hn);
        for j in 0..hn {
            *xt.add(j) = zint_mod_small_signed(fd.add(j), dlen, hn, p, p0i, r2, rx);
        }
    }

    // Fd is no longer needed.
    let mut t1 = fd;

    // Compute F (unreduced) modulo sufficiently many small primes, un-NTT-ing
    // (f,g) as we go.
    for i in 0..llen {
        // Exactly slen primes processed: (f,g) are in RNS, rebuild them.
        if i == slen {
            zint_rebuild_crt(ft, slen, n, 2, true, t1);
        }

        let p = PRIMES[i].p;
        let p0i = PRIMES[i].p0i;
        let r2 = PRIMES[i].r2;

        // Memory: gm (n), igm (n), gx (n: temporary g mod p, NTT).
        let gm = t1;
        let igm = gm.add(n);
        let gx = igm.add(n);
        mp_mkgmigm(logn, sm(gm, n), sm(igm, n), PRIMES[i].g, PRIMES[i].ig, p, p0i);
        if i < slen {
            ptr::copy_nonoverlapping(gt.add(i * n), gx, n);
            mp_intt(logn, sm(ft.add(i * n), n), sc(igm, n), p, p0i);
            mp_intt(logn, sm(gt.add(i * n), n), sc(igm, n), p, p0i);
        } else {
            let rx = mp_rx31(slen as u32, p, p0i, r2);
            for j in 0..n {
                *gx.add(j) = zint_mod_small_signed(gt.add(j), slen, n, p, p0i, r2, rx);
            }
            mp_ntt(logn, sm(gx, n), sc(gm, n), p, p0i);
        }

        // F from the deeper level, in RNS, in Ft: apply the NTT modulo p.
        let fe = ft_big.add(i * n);
        mp_ntt(logn - 1, sm(fe.add(hn), hn), sc(gm, hn), p, p0i);

        // Compute F (unreduced) modulo p.
        for j in 0..hn {
            let ga = *gx.add((j << 1) + 0);
            let gb = *gx.add((j << 1) + 1);
            let m_fp = mp_mmul(*fe.add(j + hn), r2, p, p0i);
            *fe.add((j << 1) + 0) = mp_mmul(gb, m_fp, p, p0i);
            *fe.add((j << 1) + 1) = mp_mmul(ga, m_fp, p, p0i);
        }

        // The new F in RNS only (no NTT).
        mp_intt(logn, sm(fe, n), sc(igm, n), p, p0i);
    }

    // We no longer need g.
    t1 = gt;

    // Edge case: slen == llen, f not yet rebuilt into plain representation.
    if slen == llen {
        zint_rebuild_crt(ft, slen, n, 1, true, t1);
    }

    // Rebuild the plain representation of the unreduced F.
    zint_rebuild_crt(ft_big, llen, n, 1, true, t1);

    // Reduce F with Babai's nearest plane algorithm, using fixed-point
    // approximations k = round(F/f) (see the C comments for the derivation).

    // With depth >= 2 and the NTT used for k*f, f needs an extra word.
    let use_sub_ntt = depth > 1 && logn >= MIN_LOGN_FGNTT;
    if use_sub_ntt {
        t1 = t1.add(n);
    }

    // Layout: Ft (llen*n), ft (slen*n, +n if use_sub_ntt), rt3 (n fxr).
    let rt3 = t1 as *mut u64;

    // Only the top rlen words of f are considered.
    let mut rlen = WORD_WIN[depth as usize];
    if rlen > slen {
        rlen = slen;
    }
    let blen = slen - rlen;
    let ftb = ft.add(blen * n);
    let scale_fg = 31 * blen as u32;
    let mut scale_ff = 31 * llen as u32; // scale_FG

    // Convert f into fixed-point approximations (see the C comments on the
    // choice of scale_t, which keeps f*adj(f) from overflowing).
    let scale_x = poly_max_bitlength(logn, ftb, rlen);
    let mut scale_t = 15 - logn;
    scale_t ^= (scale_t ^ scale_x) & tbmask(scale_x.wrapping_sub(scale_t));
    let scdiff = scale_x.wrapping_sub(scale_t);

    poly_big_to_fixed(logn, sm(rt3, n), ftb, rlen, scdiff);

    // rt3 <- adj(f)/(f*adj(f))  (FFT)
    vect_fft(logn, sm(rt3, n));
    vect_inv_mul2e_fft(logn, sm(rt3, n), scale_t);

    // Layout: ..., rt3 (n fxr), rt1 (n fxr) | k (n) , t2 (3*n).
    let rt1 = rt3.add(n);
    let k = rt1 as *mut i32;
    let t2 = k.add(n) as *mut u32;

    // If poly_sub_scaled_ntt() will be used, convert f to the NTT
    // representation (over slen+1 words).
    if use_sub_ntt {
        let gm = t2;
        let mut tn = gm.add(n);
        for i in 0..=slen {
            let p = PRIMES[i].p;
            let p0i = PRIMES[i].p0i;
            let r2 = PRIMES[i].r2;
            let rx = mp_rx31(slen as u32, p, p0i, r2);
            mp_mkgm(logn, sm(gm, n), PRIMES[i].g, p, p0i);
            for j in 0..n {
                *tn.add(j) = zint_mod_small_signed(ft.add(j), slen, n, p, p0i, r2, rx);
            }
            mp_ntt(logn, sm(tn, n), sc(gm, n), p, p0i);
            tn = tn.add(n);
        }
        let tn0 = gm.add(n);
        ptr::copy(tn0, ft, (slen + 1) * n);
    }

    // Reduce F repeatedly.
    let mut fglen = llen;
    let reduce_bits: u32 = match logn_top {
        9 => 9,
        10 => 8,
        _ => 10,
    };
    loop {
        // Convert the current F into fixed-point (scaling scale_FG + scale_x).
        let (tlen, toff) = divrem31(scale_ff);
        poly_big_to_fixed(
            logn,
            sm(rt1, n),
            ft_big.add(tlen as usize * n),
            fglen - tlen as usize,
            scale_x.wrapping_add(toff),
        );

        // rt1 <- (F*adj(f)) / (f*adj(f))
        vect_fft(logn, sm(rt1, n));
        vect_mul_fft(logn, sm(rt1, n), sc(rt3, n));
        vect_ifft(logn, sm(rt1, n));

        // k <- round(rt1); k aliases rt1 (ascending order is safe).
        for i in 0..n {
            let v = fxr_round(*rt1.add(i));
            *k.add(i) = v;
        }

        // k is scaled by scale_FG - scale_fg, which is public.
        let scale_k = scale_ff.wrapping_sub(scale_fg);
        if depth == 1 {
            poly_sub_kf_scaled_depth1(logn_top, ft_big, fglen, k as *mut u32, scale_k, f, t2);
        } else if use_sub_ntt {
            poly_sub_scaled_ntt(logn, ft_big, fglen, ft, slen, k, scale_k, t2);
        } else {
            poly_sub_scaled(logn, ft_big, fglen, ft, slen, k, scale_k);
        }

        // Assume F has shrunk by at least reduce_bits; adjust FGlen.
        if scale_ff <= scale_fg {
            break;
        }
        if scale_ff <= scale_fg + reduce_bits {
            scale_ff = scale_fg;
        } else {
            scale_ff -= reduce_bits;
        }
        while fglen > slen && 31 * (fglen - slen) > (scale_ff - scale_fg + 30) as usize {
            // Decrement FGlen, checking that the removed words are
            // redundant with the remaining ones.
            fglen -= 1;
            let xp = ft_big.add((fglen - 1) << logn);
            for i in 0..n {
                let sw = (*xp.add(i) >> 30).wrapping_neg() >> 1;
                if *xp.add(i + n) != sw {
                    return false;
                }
            }
        }
    }

    true
}

/// Top recursion level (depth 0): a specialized variant with lower RAM and
/// better precision. F and G are returned at the start of tmp.
unsafe fn solve_ntru_depth0(logn: u32, f: &[i8], g: &[i8], tmp: *mut u32) -> bool {
    let n = 1usize << logn;
    let hn = n >> 1;

    // At depth 0 all values fit on 30 bits: a single modulus p.
    let p = PRIMES[0].p;
    let p0i = PRIMES[0].p0i;
    let r2 = PRIMES[0].r2;

    // On input, Fd from the upper level (hn words) is at the start of tmp.
    let t1 = tmp;
    let t2 = t1.add(n);
    let mut t3 = t2.add(n);
    let t4 = t3.add(n);
    let t5 = t4.add(n);

    // Convert Fd to RNS+NTT, into t3.
    let gm = t4;
    mp_mkgm(logn, sm(gm, n), PRIMES[0].g, p, p0i);
    poly_mp_set(logn - 1, sm(t1, hn), p);
    mp_ntt(logn - 1, sm(t1, hn), sc(gm, hn), p, p0i);
    ptr::copy_nonoverlapping(t1, t3, hn);

    // Compute F (unreduced, RNS+NTT) into t1.
    poly_mp_set_small(logn, sm(t2, n), g, p);
    mp_ntt(logn, sm(t2, n), sc(gm, n), p, p0i);
    for i in 0..hn {
        let ga = *t2.add((i << 1) + 0);
        let gb = *t2.add((i << 1) + 1);
        let m_f = mp_mmul(*t3.add(i), r2, p, p0i);
        *t1.add((i << 1) + 0) = mp_mmul(gb, m_f, p, p0i);
        *t1.add((i << 1) + 1) = mp_mmul(ga, m_f, p, p0i);
    }

    // Convert f to RNS+NTT (t3), checking that f is invertible modulo p.
    poly_mp_set_small(logn, sm(t3, n), f, p);
    mp_ntt(logn, sm(t3, n), sc(gm, n), p, p0i);
    for i in 0..n {
        if *t3.add(i) == 0 {
            return false;
        }
    }

    // Last reduction, with the full expression
    //   k = round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
    // with G = (q + g*F)/f computed modulo p.
    for i in 0..hn {
        let tf0 = *t3.add(i);
        let tf1 = *t3.add(n - 1 - i);
        let tg0 = *t2.add(i);
        let tg1 = *t2.add(n - 1 - i);
        let tf_big0 = *t1.add(i);
        let tf_big1 = *t1.add(n - 1 - i);
        let mf0 = mp_mmul(tf0, r2, p, p0i);
        let mf1 = mp_mmul(tf1, r2, p, p0i);
        let mg0 = mp_mmul(tg0, r2, p, p0i);
        let mg1 = mp_mmul(tg1, r2, p, p0i);
        let tg_big0 = mp_div(mp_add(Q, mp_mmul(mg0, tf_big0, p, p0i), p), tf0, p);
        let tg_big1 = mp_div(mp_add(Q, mp_mmul(mg1, tf_big1, p, p0i), p), tf1, p);
        let kn0 = mp_add(mp_mmul(mf1, tf_big0, p, p0i), mp_mmul(mg1, tg_big0, p, p0i), p);
        let kn1 = mp_add(mp_mmul(mf0, tf_big1, p, p0i), mp_mmul(mg0, tg_big1, p, p0i), p);
        let kd = mp_add(mp_mmul(mf0, tf1, p, p0i), mp_mmul(mg0, tg1, p, p0i), p);
        *t2.add(i) = kn0;
        *t2.add(n - 1 - i) = kn1;
        *t3.add(i) = kd;
        *t3.add(n - 1 - i) = kd;
    }

    // Convert back numerator and denominator to plain integers (no
    // truncation to 31 bits).
    mp_mkigm(logn, sm(t4, n), PRIMES[0].ig, p, p0i);
    mp_intt(logn, sm(t2, n), sc(t4, n), p, p0i);
    mp_intt(logn, sm(t3, n), sc(t4, n), p, p0i);
    for i in 0..n {
        *t2.add(i) = mp_norm(*t2.add(i), p) as u32;
        *t3.add(i) = mp_norm(*t3.add(i), p) as u32;
    }

    // Divide t2 by t3 and round, in FFT representation downscaled by 2^10.
    const DOWNSCALE: u32 = 10;
    let rt4 = t4 as *mut u64;
    for i in 0..n {
        let x = ((*t3.add(i) as i32 as i64) as u64) << (32 - DOWNSCALE);
        *rt4.add(i) = x;
    }
    vect_fft(logn, sm(rt4, n));
    ptr::copy(rt4, rt4.add(hn), hn);
    let rt5 = rt4.add(hn);
    let rt3 = t3 as *mut u64;
    for i in 0..n {
        let x = ((*t2.add(i) as i32 as i64) as u64) << (32 - DOWNSCALE);
        *rt3.add(i) = x;
    }
    vect_fft(logn, sm(rt3, n));

    // Divide F*adj(f) + G*adj(g) by f*adj(f) + g*adj(g), round into t2 (RNS).
    vect_div_selfadj_fft(logn, sm(rt3, n), sc(rt5, hn));
    vect_ifft(logn, sm(rt3, n));
    for i in 0..n {
        *t2.add(i) = mp_set(fxr_round(*rt3.add(i)), p);
    }

    // Get back f and g, and convert all polynomials to RNS+NTT.
    let gm = t5;
    mp_mkgm(logn, sm(gm, n), PRIMES[0].g, p, p0i);
    t3 = t4.sub(n);
    poly_mp_set_small(logn, sm(t3, n), f, p);
    poly_mp_set_small(logn, sm(t4, n), g, p);
    mp_ntt(logn, sm(t2, n), sc(gm, n), p, p0i);
    mp_ntt(logn, sm(t3, n), sc(gm, n), p, p0i);
    mp_ntt(logn, sm(t4, n), sc(gm, n), p, p0i);

    // Reduce F by subtracting k*F and recompute G = (q + g*F)/f.
    for i in 0..n {
        let mut t_f = *t1.add(i);
        let tk = *t2.add(i);
        let tf = *t3.add(i);
        let tg = *t4.add(i);
        let mf = mp_mmul(tf, r2, p, p0i);
        let mg = mp_mmul(tg, r2, p, p0i);
        t_f = mp_sub(t_f, mp_mmul(mf, tk, p, p0i), p);
        let t_g = mp_div(mp_add(Q, mp_mmul(mg, t_f, p, p0i), p), tf, p);
        *t1.add(i) = t_f;
        *t2.add(i) = t_g;
    }

    // Convert back F and G into normal representation.
    mp_mkigm(logn, sm(t3, n), PRIMES[0].ig, p, p0i);
    mp_intt(logn, sm(t1, n), sc(t3, n), p, p0i);
    mp_intt(logn, sm(t2, n), sc(t3, n), p, p0i);
    poly_mp_norm(logn, sm(t1, n), p);
    poly_mp_norm(logn, sm(t2, n), p);

    true
}

/// Solve the NTRU equation for (f,g). On success, (F,G) are returned at the
/// start of tmp as two consecutive `i8` arrays. tmp needs 6*n words.
///
/// # Safety
/// `tmp` must be valid for `6 << logn` words, 8-byte aligned.
pub unsafe fn solve_ntru(logn: u32, f: &[i8], g: &[i8], tmp: *mut u32) -> bool {
    let n = 1usize << logn;

    if !solve_ntru_deepest(logn, f, g, tmp) {
        return false;
    }
    let mut depth = logn;
    while depth > 1 {
        depth -= 1;
        if !solve_ntru_intermediate(logn, f, g, depth, tmp) {
            return false;
        }
    }
    if !solve_ntru_depth0(logn, f, g, tmp) {
        return false;
    }

    // F and G are at the start of tmp (plain, 31 bits per value): convert
    // to 8-bit representation and check the range.
    let big_f = tmp.add(2 * n) as *mut i8;
    let big_g = big_f.add(n);
    let lim = 127;
    if !poly_big_to_small(logn, sm(big_f, n), sc(tmp, n), lim) {
        return false;
    }
    if !poly_big_to_small(logn, sm(big_g, n), sc(tmp.add(n), n), lim) {
        return false;
    }
    ptr::copy(big_f as *const u8, tmp as *mut u8, 2 * n);
    true
}

/// Check that (f,g) has an acceptable orthogonalized norm. `tmp` has room
/// for 2.5*n fixed-point values.
///
/// # Safety
/// `tmp` must be valid for `5 << (logn - 1)` `u64` values.
pub unsafe fn check_ortho_norm(logn: u32, f: &[i8], g: &[i8], tmp: *mut u64) -> bool {
    let n = 1usize << logn;
    let rt1 = sm(tmp, n);
    let rt2 = sm(tmp.add(n), n);
    // rt3 is self-adjoint: only n/2 values (the f,g area starts right after).
    let rt3 = sm(tmp.add(2 * n), n >> 1);
    vect_set(logn, rt1, f);
    vect_set(logn, rt2, g);
    vect_fft(logn, rt1);
    vect_fft(logn, rt2);
    vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
    vect_adj_fft(logn, rt1);
    vect_adj_fft(logn, rt2);
    vect_mul_realconst(logn, rt1, fxr_of(12289));
    vect_mul_realconst(logn, rt2, fxr_of(12289));
    vect_mul_selfadj_fft(logn, rt1, rt3);
    vect_mul_selfadj_fft(logn, rt2, rt3);
    vect_ifft(logn, rt1);
    vect_ifft(logn, rt2);
    let mut sn: Fxr = 0;
    for i in 0..n {
        sn = fxr_add(sn, fxr_add(fxr_sqr(rt1[i]), fxr_sqr(rt2[i])));
    }
    // (0.999*1.17*sqrt(q))^2 scaled by 2^32.
    fxr_lt(sn, 72107278641426)
}
