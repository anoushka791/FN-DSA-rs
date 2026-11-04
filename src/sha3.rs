//! SHAKE on top of the (mechanically translated) Keccak-f[1600] permutation.

use crate::keccak_gen::process_block;

#[derive(Clone)]
pub struct Shake {
    a: [u64; 25],
    dptr: usize,
    rate: usize,
}

impl Shake {
    /// `size` = 128 or 256.
    pub fn new(size: u32) -> Self {
        Shake { a: [0; 25], dptr: 0, rate: 200 - (size as usize >> 2) }
    }

    pub fn inject(&mut self, mut data: &[u8]) {
        while !data.is_empty() {
            let clen = core::cmp::min(self.rate - self.dptr, data.len());
            let (chunk, rest) = data.split_at(clen);
            self.xor_in(chunk);
            self.dptr += clen;
            data = rest;
            if self.dptr == self.rate {
                process_block(&mut self.a);
                self.dptr = 0;
            }
        }
    }

    /// XOR `chunk` into the state at byte offset `dptr` (little-endian lanes).
    #[inline]
    fn xor_in(&mut self, mut chunk: &[u8]) {
        let mut pos = self.dptr;
        // Unaligned head.
        while pos & 7 != 0 && !chunk.is_empty() {
            self.a[pos >> 3] ^= (chunk[0] as u64) << ((pos & 7) << 3);
            chunk = &chunk[1..];
            pos += 1;
        }
        // Whole words.
        let mut it = chunk.chunks_exact(8);
        for w in &mut it {
            self.a[pos >> 3] ^= u64::from_le_bytes(w.try_into().unwrap());
            pos += 8;
        }
        // Tail.
        for &b in it.remainder() {
            self.a[pos >> 3] ^= (b as u64) << ((pos & 7) << 3);
            pos += 1;
        }
    }

    pub fn flip(&mut self) {
        let v = self.dptr;
        self.a[v >> 3] ^= 0x1Fu64 << ((v & 7) << 3);
        let v = self.rate - 1;
        self.a[v >> 3] ^= 0x80u64 << ((v & 7) << 3);
        self.dptr = self.rate;
    }

    pub fn extract(&mut self, out: &mut [u8]) {
        let mut off = 0;
        while off < out.len() {
            if self.dptr == self.rate {
                process_block(&mut self.a);
                self.dptr = 0;
            }
            let clen = core::cmp::min(self.rate - self.dptr, out.len() - off);
            for k in 0..clen {
                let p = self.dptr + k;
                out[off + k] = (self.a[p >> 3] >> ((p & 7) << 3)) as u8;
            }
            self.dptr += clen;
            off += clen;
        }
    }

    // Streaming reads (little-endian), used by the sampler PRNG. They are
    // equivalent to extracting the next 1, 2 or 8 bytes of output. The fast
    // paths read directly from the state; the (rare) block-boundary cases go
    // through `extract`.

    #[inline(always)]
    pub fn next_u8(&mut self) -> u8 {
        if self.dptr < self.rate {
            let v = self.byte_at(self.dptr);
            self.dptr += 1;
            v
        } else {
            self.next_slow::<1>() as u8
        }
    }

    #[inline(always)]
    pub fn next_u16(&mut self) -> u32 {
        if self.dptr + 1 < self.rate {
            let v = self.u16_at(self.dptr);
            self.dptr += 2;
            v
        } else {
            self.next_slow::<2>() as u32
        }
    }

    #[inline(always)]
    pub fn next_u64(&mut self) -> u64 {
        if self.dptr + 7 < self.rate {
            let v = self.u64_at(self.dptr);
            self.dptr += 8;
            v
        } else {
            self.next_slow::<8>()
        }
    }

    #[cold]
    #[inline(never)]
    fn next_slow<const N: usize>(&mut self) -> u64 {
        let mut x = [0u8; 8];
        self.extract(&mut x[..N]);
        u64::from_le_bytes(x)
    }

    // Raw access to the state as a little-endian byte string. Callers
    // guarantee `pos + width <= rate <= 168 <= 200`.

    #[cfg(target_endian = "little")]
    #[inline(always)]
    fn byte_at(&self, pos: usize) -> u8 {
        debug_assert!(pos < 200);
        // SAFETY: pos < rate <= 168, inside the 200-byte state.
        unsafe { *(self.a.as_ptr() as *const u8).add(pos) }
    }

    #[cfg(target_endian = "little")]
    #[inline(always)]
    fn u16_at(&self, pos: usize) -> u32 {
        debug_assert!(pos + 2 <= 200);
        // SAFETY: pos + 2 <= rate <= 168, inside the 200-byte state.
        unsafe { core::ptr::read_unaligned((self.a.as_ptr() as *const u8).add(pos) as *const u16) as u32 }
    }

    #[cfg(target_endian = "little")]
    #[inline(always)]
    fn u64_at(&self, pos: usize) -> u64 {
        debug_assert!(pos + 8 <= 200);
        // SAFETY: pos + 8 <= rate <= 168, inside the 200-byte state.
        unsafe { core::ptr::read_unaligned((self.a.as_ptr() as *const u8).add(pos) as *const u64) }
    }

    #[cfg(not(target_endian = "little"))]
    #[inline(always)]
    fn byte_at(&self, pos: usize) -> u8 {
        (self.a[pos >> 3] >> ((pos & 7) << 3)) as u8
    }

    #[cfg(not(target_endian = "little"))]
    #[inline(always)]
    fn u16_at(&self, pos: usize) -> u32 {
        self.byte_at(pos) as u32 | ((self.byte_at(pos + 1) as u32) << 8)
    }

    #[cfg(not(target_endian = "little"))]
    #[inline(always)]
    fn u64_at(&self, pos: usize) -> u64 {
        let lane = pos >> 3;
        let sh = (pos & 7) << 3;
        if sh == 0 {
            self.a[lane]
        } else {
            (self.a[lane] >> sh) | (self.a[lane + 1] << (64 - sh))
        }
    }
}
