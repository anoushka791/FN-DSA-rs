//! Decoding primitives needed for verification.

use crate::mq::B_INF;

/// Decode a polynomial with 14 bits per coefficient (values must be in
/// [0, q-1]). Returns the number of consumed bytes, or 0 on error.
/// `d` must contain at least `7 * n / 4` bytes and `h` must hold `n` values.
pub fn mqpoly_decode(logn: u32, d: &[u8], h: &mut [u16]) -> usize {
    let n = 1usize << logn;
    let h = &mut h[..n];
    let d = &d[..7 * (n >> 2)];
    let mut ov: u32 = 0xFFFF_FFFF;
    // logn >= 2, so n is a multiple of 4.
    for (hq, dq) in h.chunks_exact_mut(4).zip(d.chunks_exact(7)) {
        let d0 = dq[0] as u32;
        let d1 = dq[1] as u32;
        let d2 = dq[2] as u32;
        let d3 = dq[3] as u32;
        let d4 = dq[4] as u32;
        let d5 = dq[5] as u32;
        let d6 = dq[6] as u32;
        let h0 = (d0 | (d1 << 8)) & 0x3FFF;
        let h1 = ((d1 >> 6) | (d2 << 2) | (d3 << 10)) & 0x3FFF;
        let h2 = ((d3 >> 4) | (d4 << 4) | (d5 << 12)) & 0x3FFF;
        let h3 = (d5 >> 2) | (d6 << 6);
        hq[0] = h0 as u16;
        hq[1] = h1 as u16;
        hq[2] = h2 as u16;
        hq[3] = h3 as u16;
        ov &= h0.wrapping_sub(12289);
        ov &= h1.wrapping_sub(12289);
        ov &= h2.wrapping_sub(12289);
        ov &= h3.wrapping_sub(12289);
    }
    if (ov >> 16) == 0 {
        0
    } else {
        d.len()
    }
}

/// Decode a compressed (Golomb-Rice) signature polynomial. Returns false on
/// any error (truncated input, invalid encoding, non-canonical padding,
/// |value| > B_INF).
pub fn comp_decode(logn: u32, d: &[u8], s: &mut [u16]) -> bool {
    let n = 1usize << logn;
    let dlen = d.len();
    let mut acc: u32 = 0;
    let mut acc_len: u32 = 0;
    let mut j = 0usize;
    for si in s[..n].iter_mut() {
        // Invariant: acc_len <= 7.
        if j >= dlen {
            return false;
        }
        acc |= (d[j] as u32) << acc_len;
        j += 1;
        let t = acc & 1;
        let mut m = (acc >> 1) & 0x7F;
        acc >>= 8;

        if acc == 0 {
            if j >= dlen {
                return false;
            }
            acc |= (d[j] as u32) << acc_len;
            j += 1;
            acc_len += 8;
            if acc == 0 {
                return false;
            }
        }
        let tz = acc.trailing_zeros();
        m += tz << 7;
        if m > B_INF {
            return false;
        }
        acc >>= tz + 1;
        acc_len -= tz + 1;

        // Reject "-0".
        if m == 0 && t != 0 {
            return false;
        }
        m = (m ^ t.wrapping_neg()).wrapping_add(t);
        *si = m as u16;
    }
    if acc != 0 {
        return false;
    }
    d[j..].iter().all(|&b| b == 0)
}

/// Decode a small polynomial (`nbits` bits per element, sign-extended).
/// Returns the number of consumed bytes. An element equal to -2^(nbits-1)
/// is not reported as an error.
pub fn trim_i8_decode(logn: u32, d: &[u8], f: &mut [i8], nbits: u32) -> usize {
    let needed = ((nbits as usize) << logn) >> 3;
    let n = 1usize << logn;
    if nbits == 8 {
        for (fv, &dv) in f[..n].iter_mut().zip(d[..n].iter()) {
            *fv = dv as i8;
        }
        return needed;
    }
    let mask1 = (1u32 << nbits) - 1;
    let mask2 = 1u32 << (nbits - 1);
    let mut j = 0usize;
    let mut acc = 0u32;
    let mut acc_len = 0u32;
    for &b in d[..needed].iter() {
        acc |= (b as u32) << acc_len;
        acc_len += 8;
        while acc_len >= nbits {
            let mut w = acc & mask1;
            acc >>= nbits;
            acc_len -= nbits;
            w |= (w & mask2).wrapping_neg();
            f[j] = w as i32 as i8;
            j += 1;
        }
    }
    needed
}

/// Encode a signed polynomial with the compressed (Golomb-Rice) format into
/// `d` (the whole buffer is written, padded with zeros). Returns false if a
/// value is out of range or the buffer is too small.
pub fn comp_encode(logn: u32, s: &[i16], d: &mut [u8]) -> bool {
    let n = 1usize << logn;
    let dlen = d.len();
    let mut acc = 0u32;
    let mut acc_len = 0u32;
    let mut j = 0usize;
    for &sv in s[..n].iter() {
        // Invariant: acc_len <= 7
        let x = sv as i32;
        let sw = (x >> 16) as u32;
        let w = ((x as u32) ^ sw).wrapping_sub(sw);
        if w > B_INF {
            return false;
        }

        // Sign bit and low 7 bits of the absolute value.
        acc |= ((sw & 1) | ((w & 0x7F) << 1)) << acc_len;
        acc_len += 8;

        // High bits, in unary.
        let wh = w >> 7;
        acc |= 1u32 << (acc_len + wh);
        acc_len += wh + 1;

        while acc_len >= 8 {
            if j >= dlen {
                return false;
            }
            d[j] = acc as u8;
            j += 1;
            acc >>= 8;
            acc_len -= 8;
        }
    }
    if acc_len > 0 {
        if j >= dlen {
            return false;
        }
        d[j] = acc as u8;
        j += 1;
    }
    d[j..].fill(0);
    true
}

/// Encode a small polynomial with `nbits` bits per element; returns the
/// number of bytes written. `nbits * 2^logn` must be a multiple of 8.
pub fn trim_i8_encode(logn: u32, f: &[i8], nbits: u32, d: &mut [u8]) -> usize {
    let n = 1usize << logn;
    if nbits == 8 {
        for (dv, &fv) in d[..n].iter_mut().zip(f[..n].iter()) {
            *dv = fv as u8;
        }
        return n;
    }
    let mut j = 0usize;
    let mut acc = 0u32;
    let mut acc_len = 0u32;
    let mask = (1u32 << nbits) - 1;
    for &fv in f[..n].iter() {
        acc |= ((fv as i32 as u32) & mask) << acc_len;
        acc_len += nbits;
        if acc_len >= 8 {
            d[j] = acc as u8;
            j += 1;
            acc >>= 8;
            acc_len -= 8;
        }
    }
    j
}

/// Encode a polynomial (values in [0, q-1]) with 14 bits per value; returns
/// the number of bytes written.
pub fn mqpoly_encode(logn: u32, h: &[u16], d: &mut [u8]) -> usize {
    let n = 1usize << logn;
    let mut j = 0usize;
    for hq in h[..n].chunks_exact(4) {
        let h0 = hq[0] as u32;
        let h1 = hq[1] as u32;
        let h2 = hq[2] as u32;
        let h3 = hq[3] as u32;
        d[j + 0] = h0 as u8;
        d[j + 1] = ((h0 >> 8) | (h1 << 6)) as u8;
        d[j + 2] = (h1 >> 2) as u8;
        d[j + 3] = ((h1 >> 10) | (h2 << 4)) as u8;
        d[j + 4] = (h2 >> 4) as u8;
        d[j + 5] = ((h2 >> 12) | (h3 << 2)) as u8;
        d[j + 6] = (h3 >> 6) as u8;
        j += 7;
    }
    j
}
