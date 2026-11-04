//! Core signature generation (port of sign_core.c).
//!
//! All intermediate values live in a single temporary area of 20*n bytes,
//! which is carved into overlapping typed views exactly as the C code does.
//! The view names (`tmpK`) are the byte offset K*n inside that area. Each
//! view used in a given call is disjoint from the others used in that call.

use crate::codec::{comp_encode, trim_i8_decode};
use crate::fpoly;
use crate::mq::{self, sqnorm_is_acceptable};
use crate::mqsign as ms;
use crate::sampler::{ffsamp_fft, SamplerState, FLAG_INIT0, FLAG_NOOUT};
use crate::sha3::Shake;
use crate::util::hash_to_point;
use crate::{fg_nbits, signature_size};
use core::slice::{from_raw_parts as sl, from_raw_parts_mut as slm};

/// Sign the message representative `mu` with the encoded key material
/// `fg_f` (f, g, F, without the leading byte), using `seed` for all
/// randomness. Returns the signature length, or 0 on error.
///
/// # Safety
/// `tmp` must point to at least `20 << logn` bytes, aligned to 32 bytes.
/// `sig` must hold at least `signature_size(logn)` bytes, `fg_f` the whole
/// encoded (f, g, F).
pub unsafe fn sign_core(
    logn: u32,
    fg_f: &[u8],
    mu: &[u8; 64],
    seed: &[u8],
    sig: &mut [u8],
    tmp: *mut u8,
) -> usize {
    let k = mq::best();
    let n = 1usize << logn;
    let nbits = fg_nbits(logn);
    let flen = ((nbits as usize) << logn) >> 3;
    let fk = &fg_f[..flen];
    let gk = &fg_f[flen..2 * flen];
    // Storage format of F is already an array of int8_t.
    let big_f: &[i8] = sl(fg_f.as_ptr().add(2 * flen) as *const i8, n);

    // Typed views into the temporary area (byte offset = K*n).
    macro_rules! u16p { ($k:expr) => { tmp.add($k * n) as *mut u16 }; }
    macro_rules! u16v { ($k:expr) => { slm(u16p!($k), n) }; }
    macro_rules! i16v { ($k:expr) => { slm(tmp.add($k * n) as *mut i16, n) }; }
    macro_rules! i8v { ($k:expr, $off:expr) => { slm(tmp.add($k * n + $off) as *mut i8, n) }; }
    macro_rules! fp { ($k:expr) => { tmp.add($k * n) as *mut f64 }; }

    let mut ret = 0usize;
    let sig_len = signature_size(logn);

    // The loop usually runs once; we force an exit after 27 iterations.
    for counter in 0u8..27 {
        // PRNG for the nonce and the sampler.
        let mut ss = SamplerState { pc: Shake::new(256), dz: core::ptr::null_mut(), off: 0, logn };
        ss.pc.inject(seed);
        ss.pc.inject(&[counter]);
        ss.pc.flip();
        ss.pc.extract(&mut sig[1..41]);

        // Step 1: compute q^2/(f*adj(f) + g*adj(g)) (FFT) and c*f mod q.

        // tmp4 <- f (first half of tmp4, as int8)
        trim_i8_decode(logn, fk, i8v!(4, 0), nbits);
        let sqnf = ms::sqnorm_small(logn, i8v!(4, 0));

        // tmp0 <- f mod q' (NTT)
        ms::mq2_small_to_mod(logn, i8v!(4, 0), u16v!(0));
        ms::mq2_ntt(logn, u16v!(0));

        // tmp4 (second half) <- g
        trim_i8_decode(logn, gk, i8v!(4, n), nbits);
        let sqng = ms::sqnorm_small(logn, i8v!(4, n));

        // A key with ||(g,-f)|| > 1.17*sqrt(q) is invalid.
        if sqnf.wrapping_add(sqng) > 16822 {
            return 0;
        }

        // tmp2 <- g mod q' (NTT)
        ms::mq2_small_to_mod(logn, i8v!(4, n), u16v!(2));
        ms::mq2_ntt(logn, u16v!(2));

        // tmp0 <- f*adj(f) + g*adj(g) mod q' (normalized)
        ms::mq2_muladj_x2_ntt(logn, u16v!(0), u16v!(2));
        ms::mq2_intt(logn, u16v!(0));
        ms::mq2_mod_to_signed(logn, u16v!(0));

        // tmp18 <- f mod q (NTT)
        ms::small_to_int(logn, i8v!(4, 0), u16v!(18));
        (k.int_to_ntt)(logn, u16v!(18));

        // tmp2 <- g mod q (NTT)
        ms::small_to_int(logn, i8v!(4, n), u16v!(2));
        (k.int_to_ntt)(logn, u16v!(2));

        // tmp2 <- f*adj(f) + g*adj(g) mod q (normalized)
        ms::muladj_x2_ntt(logn, u16v!(2), u16v!(18));
        (k.ntt_to_int)(logn, u16v!(2));
        ms::int_to_signed(logn, u16v!(2));

        // Adjust the value in tmp0 to match the one in tmp2; save the
        // adjustment bits in sig + 41.
        ms::mq2mq_adjust_save(logn - 1, i16v!(0), i16v!(2), &mut sig[41..]);

        // tmp8 <- q^2/(f*adj(f) + g*adj(g)) (FFT, self-adjoint)
        fpoly::fft_int16_selfadj(logn, fp!(8), tmp as *const i16);
        fpoly::invmul_selfadj_fft(logn, slm(fp!(8), n >> 1), fpoly::FPR_QQ);

        // tmp16 <- c*f mod q (normalized)
        let nonce: [u8; 40] = sig[1..41].try_into().unwrap();
        hash_to_point(logn, &nonce, mu, u16v!(16));
        (k.ext_to_int)(logn, u16v!(16));
        (k.int_to_ntt)(logn, u16v!(16));
        (k.mul_ntt)(logn, u16v!(16), u16v!(18));
        (k.ntt_to_int)(logn, u16v!(16));
        ms::int_to_signed(logn, u16v!(16));

        // Step 2: ffSampling on t1.
        ss.dz = u16p!(16) as *mut i16;
        ffsamp_fft(&mut ss, FLAG_INIT0, tmp as *mut f64);

        // Step 3: we have t1 - z1 (FFT) in tmp0.

        // tmp12 <- f mod q (NTT); value already in tmp18.
        core::ptr::copy_nonoverlapping(u16p!(18), u16p!(12), n);

        // tmp18 <- 1/f mod q (NTT)
        ms::inv_ntt(logn, u16v!(18), u16v!(10));

        // tmp14 <- g mod q (NTT)
        trim_i8_decode(logn, gk, i8v!(10, 0), nbits);
        ms::small_to_int(logn, i8v!(10, 0), u16v!(14));
        (k.int_to_ntt)(logn, u16v!(14));

        // tmp8 <- F mod q (NTT)
        ms::small_to_int(logn, big_f, u16v!(8));
        (k.int_to_ntt)(logn, u16v!(8));

        // tmp10 <- G = (q + g*F)/f mod q (NTT)
        core::ptr::copy_nonoverlapping(u16p!(14), u16p!(10), n);
        (k.mul_ntt)(logn, u16v!(10), u16v!(8));
        (k.mul_ntt)(logn, u16v!(10), u16v!(18));

        // tmp8 <- F*adj(f) + G*adj(g) mod q (normalized)
        ms::muladj_add_muladj(logn, u16v!(8), u16v!(12), u16v!(10), u16v!(14));
        (k.ntt_to_int)(logn, u16v!(8));
        ms::int_to_signed(logn, u16v!(8));

        // tmp10 <- G (normalized)
        (k.ntt_to_int)(logn, u16v!(10));
        ms::int_to_signed(logn, u16v!(10));

        // tmp10 <- G mod q' (NTT)
        ms::mq2_signed_to_mod(logn, u16v!(10));
        ms::mq2_ntt(logn, u16v!(10));

        // tmp14 <- g mod q' (NTT)
        trim_i8_decode(logn, gk, i8v!(18, 0), nbits);
        ms::mq2_small_to_mod(logn, i8v!(18, 0), u16v!(14));
        ms::mq2_ntt(logn, u16v!(14));

        // tmp10 <- G*adj(g) mod q' (NTT)
        ms::mq2_muladj_ntt(logn, u16v!(10), u16v!(14));

        // tmp12 <- f mod q' (NTT)
        trim_i8_decode(logn, fk, i8v!(18, 0), nbits);
        ms::mq2_small_to_mod(logn, i8v!(18, 0), u16v!(12));
        ms::mq2_ntt(logn, u16v!(12));

        // tmp18 <- f*adj(f) + g*adj(g) mod q' (normalized)
        core::ptr::copy_nonoverlapping(u16p!(12), u16p!(18), n);
        ms::mq2_muladj_x2_ntt(logn, u16v!(18), u16v!(14));
        ms::mq2_intt(logn, u16v!(18));
        ms::mq2_mod_to_signed(logn, u16v!(18));

        // Adjust tmp18 to the true value of f*adj(f) + g*adj(g).
        ms::mq2mq_adjust(logn - 1, i16v!(18), &sig[41..]);

        // tmp14 <- F mod q' (NTT)
        ms::mq2_small_to_mod(logn, big_f, u16v!(14));
        ms::mq2_ntt(logn, u16v!(14));

        // tmp10 <- F*adj(f) + G*adj(g) mod q' (normalized)
        ms::mq2_muladj_ntt(logn, u16v!(14), u16v!(12));
        ms::mq2_add_poly(logn, u16v!(10), u16v!(14));
        ms::mq2_intt(logn, u16v!(10));
        ms::mq2_mod_to_signed(logn, u16v!(10));

        // tmp12 <- F*adj(f) + G*adj(g) (int32)
        ms::mq2mq_crt(logn, slm(tmp.add(12 * n) as *mut i32, n), i16v!(8), i16v!(10));

        // tmp8 <- F*adj(f) + G*adj(g) (FFT)
        fpoly::fft_int32(logn, fp!(8), tmp.add(12 * n) as *const i32);

        // tmp0 <- (t1 - z1)*(F*adj(f) + G*adj(g)) (FFT)
        fpoly::mul_fft(logn, slm(fp!(0), n), sl(fp!(8), n));

        // tmp8 <- f*adj(f) + g*adj(g) (FFT, self-adjoint)
        fpoly::fft_int16_selfadj(logn, fp!(8), tmp.add(18 * n) as *const i16);

        // tmp0 <- tmp0 / tmp8 (FFT)
        fpoly::div_selfadj_fft(logn, slm(fp!(0), n), sl(fp!(8), n >> 1));

        // tmp12 <- c (NTT mod q)
        hash_to_point(logn, &nonce, mu, u16v!(12));
        (k.ext_to_int)(logn, u16v!(12));
        (k.int_to_ntt)(logn, u16v!(12));

        // tmp18 <- -F (NTT mod q)
        ms::small_to_int(logn, big_f, u16v!(18));
        (k.int_to_ntt)(logn, u16v!(18));
        ms::neg(logn, u16v!(18));

        // tmp18 <- -c*F mod q (normalized)
        (k.mul_ntt)(logn, u16v!(18), u16v!(12));
        (k.ntt_to_int)(logn, u16v!(18));
        ms::int_to_signed(logn, u16v!(18));

        // Step 4: ffSampling on t0'.
        ss.dz = u16p!(18) as *mut i16;
        ffsamp_fft(&mut ss, FLAG_NOOUT, tmp as *mut f64);

        // Step 5: build the signature s1 = f*(z0 + (c*F mod q)/q)
        //                               + F*(z1 - (c*f mod q)/q),
        // computed modulo q'. z0 is in tmp18 and z1 in tmp16.

        // tmp2 <- c mod q (NTT)
        hash_to_point(logn, &nonce, mu, u16v!(2));
        (k.ext_to_int)(logn, u16v!(2));
        (k.int_to_ntt)(logn, u16v!(2));

        // tmp4 <- f mod q (NTT); f[] decoded into tmp12 first.
        trim_i8_decode(logn, fk, i8v!(12, 0), nbits);
        ms::small_to_int(logn, i8v!(12, 0), u16v!(4));
        (k.int_to_ntt)(logn, u16v!(4));

        // tmp6 <- F mod q (NTT)
        ms::small_to_int(logn, big_f, u16v!(6));
        (k.int_to_ntt)(logn, u16v!(6));

        // tmp8 <- c*f mod q (normalized)
        core::ptr::copy_nonoverlapping(u16p!(2), u16p!(8), n);
        (k.mul_ntt)(logn, u16v!(8), u16v!(4));
        (k.ntt_to_int)(logn, u16v!(8));
        ms::int_to_signed(logn, u16v!(8));

        // tmp6 <- c*F mod q (normalized)
        (k.mul_ntt)(logn, u16v!(6), u16v!(2));
        (k.ntt_to_int)(logn, u16v!(6));
        ms::int_to_signed(logn, u16v!(6));

        // tmp6 <- z0 + (c*F mod q)/q (NTT mod q')
        ms::mq2_signed_to_mod(logn, u16v!(6));
        ms::mq2_mulconst(logn, u16v!(6), ms::MQ2_IQ);
        ms::mq2_signed_to_mod(logn, u16v!(18));
        ms::mq2_add_poly(logn, u16v!(6), u16v!(18));
        ms::mq2_ntt(logn, u16v!(6));

        // tmp16 <- z1 - (c*f mod q)/q (NTT mod q')
        ms::mq2_signed_to_mod(logn, u16v!(8));
        ms::mq2_mulconst(logn, u16v!(8), ms::MQ2_IQ);
        ms::mq2_signed_to_mod(logn, u16v!(16));
        ms::mq2_sub_poly(logn, u16v!(16), u16v!(8));
        ms::mq2_ntt(logn, u16v!(16));

        // tmp0 <- f mod q' (NTT)
        ms::mq2_small_to_mod(logn, i8v!(12, 0), u16v!(0));
        ms::mq2_ntt(logn, u16v!(0));

        // tmp10 <- F mod q' (NTT)
        ms::mq2_small_to_mod(logn, big_f, u16v!(10));
        ms::mq2_ntt(logn, u16v!(10));

        // tmp6 <- s1 (normalized)
        ms::mq2_mul_ntt(logn, u16v!(6), u16v!(0));
        ms::mq2_mul_ntt(logn, u16v!(16), u16v!(10));
        ms::mq2_add_poly(logn, u16v!(6), u16v!(16));
        ms::mq2_intt(logn, u16v!(6));
        ms::mq2_mod_to_signed(logn, u16v!(6));

        // tmp0 <- g mod q (NTT); g[] decoded into tmp12 first.
        trim_i8_decode(logn, gk, i8v!(12, 0), nbits);
        ms::small_to_int(logn, i8v!(12, 0), u16v!(0));
        (k.int_to_ntt)(logn, u16v!(0));

        // tmp0 <- h = g/f mod q (NTT)
        ms::div_ntt(logn, u16v!(0), u16v!(4), u16v!(12));

        // tmp4 <- s1 mod q (NTT)
        core::ptr::copy_nonoverlapping(u16p!(6), u16p!(4), n);
        (k.signed_to_int)(logn, u16v!(4));
        (k.int_to_ntt)(logn, u16v!(4));

        // tmp2 <- s0 = c - s1*h mod q (normalized)
        (k.mul_ntt)(logn, u16v!(4), u16v!(0));
        (k.sub)(logn, u16v!(2), u16v!(4));
        (k.ntt_to_int)(logn, u16v!(2));
        ms::int_to_signed(logn, u16v!(2));

        // The putative signature (s0, s1) is in (tmp2, tmp6). Acceptable if
        // the L2 and infinity norms are low enough and s1 can be encoded.
        let sqn0 = ms::sqnorm_binf_signed(logn, i16v!(2));
        let sqn1 = ms::sqnorm_binf_signed(logn, i16v!(6));
        let sqn = sqn0.wrapping_add(sqn1) | ((sqn0 | sqn1) & 0x8000_0000);
        if !sqnorm_is_acceptable(logn, sqn) {
            continue;
        }

        if comp_encode(logn, i16v!(6), &mut sig[41..sig_len]) {
            sig[0] = 0x30 + logn as u8;
            ret = sig_len;
            break;
        }
    }
    ret
}
