# fndsa-alt (Rust)

A Rust port of [c-fn-dsa-alt](https://github.com/pornin/c-fn-dsa-alt), the RAM-efficient FN-DSA (Falcon, FIPS 206 draft) implementation. Key generation, signing and verification are all ported, with no dependencies beyond `std`.

Output is byte-identical to the C code, so keys and signatures work across both implementations.

## Status

| | Matches C output | Speed vs C (AVX2) |
|---|---|---|
| Verify | yes | 4 to 10% faster |
| Sign | yes | 7 to 10% slower |
| Keygen | yes | 50 to 75% slower |

Keygen is slower because the C code uses AVX2 for the prime-modulus NTT, big-integer routines and fixed-point FFT, and none of that is ported yet. Signing lacks the SSE2 FFT routines the C code has.

Not ported: Cortex-M4, RISC-V, NEON and SSE2 paths. Only tested on x86_64 Windows.

## Usage

```rust
pub fn keygen(logn: u32, sign_key: Option<&mut [u8]>, vrfy_key: Option<&mut [u8]>) -> bool;
pub fn sign(sign_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8], sig: &mut [u8]) -> Option<usize>;
pub fn verify(sig: &[u8], vrfy_key: &[u8], ctx: &[u8], id: &[u8], hv: &[u8]) -> bool;
```

There are also `_seeded` variants for deterministic output, and `_weak` variants for the test-only degrees 4 to 256. The `id` argument works as in `fndsa.h`: `HASH_ID_RAW` for a raw message, a DER-encoded OID for a pre-hashed one, or `HASH_ID_EXTMU` with a 64-byte `hv`.

| Degree | Verifying key | Signing key | Signature |
|---|---|---|---|
| 512 | 897 B | 1281 B | 666 B |
| 1024 | 1793 B | 2305 B | 1280 B |

Signing needs `20*n` bytes of stack and key generation `22*n` bytes, the same layout as the C code. Verification allocates nothing on the heap.

## How it differs from the C code

The algorithms and arithmetic are the same. The main differences:

- **Kernel tables instead of `#if`.** One `inner_verify` takes a table of function pointers, either `PORTABLE` or `AVX2`, picked at runtime. AVX2 is only used after `is_x86_feature_detected!` succeeds.
- **Generated constants.** The NTT tables and the unrolled Keccak code come from the C sources via scripts in `tools/`, not by hand.
- **Faster SHAKE input.** Input is XORed into the state a word at a time rather than a byte at a time, which probably explains most of the verify speedup.
- **Safer inputs.** Slices instead of pointer/length pairs, and a malformed OID is rejected where the C code reads past it.
- **Native `f64`.** Signing uses plain `f64` rather than emulated floats. Results are bit-identical to the C x86 path. Keygen is integer-only.
- **`unsafe`.** To keep the single shared buffer, signing and keygen use raw pointers (`sign_core.rs`, `fpoly.rs`, `sampler.rs`, `ntru.rs`, `kpoly.rs`, `zint31.rs`), and the AVX2 code in `mq.rs` is `unsafe` behind a runtime check.
- **Q' NTT in AVX2.** The C code has no AVX2 for signing. The Rust NTT kernels are generic over the modulus, so the q' = 18433 transform reuses them.

**Careful with `impl_id`.** Signing hashes the implementation identifier `20260830`, copied from C. That is only valid while this port does the same floating-point operations in the same order. If you change the sampler or FFT, change the identifier (see `src/sign.rs`).

## Testing

Run `cargo test --release` from this directory. Test vectors come from the C code (`tools/gen_*_vectors.c`).

- **Verify:** 102 C-made signatures are accepted. Bit flips and other bad inputs are rejected. AVX2 and portable kernels agree.
- **Sign:** 144 seeded signatures match the C output byte for byte.
- **Keygen:** 161 seeded key pairs match the C output byte for byte, and keygen, sign, verify round-trips for every degree.

Not covered: timing side channels (no review of the Rust code), peak stack use, and non-x86_64 targets.

## Benchmarks

AMD Ryzen 5 220, Windows 11, single thread. C is clang 22 `-O2`; Rust is rustc 1.99 `--release` with fat LTO. Median per call, best of three runs.

| | C | Rust |
|---|---|---|
| Verify n=512 | 17.7 µs | 15.9 µs |
| Verify n=1024 | 31.4 µs | 30.2 µs |
| Sign n=512 | 356 µs | 381 µs |
| Sign n=1024 | 725 µs | 777 µs |
| Keygen n=512 | 3.43 ms | 5.99 ms |
| Keygen n=1024 | 20.7 ms | 31.3 ms |

The sign numbers were taken with about 40% background CPU load, so trust the ratio, not the absolute values. Run-to-run noise is a few percent, so the 4% verify gap at n=1024 is within noise. The portable (non-AVX2) verify is about 1.9x slower than AVX2. Everything here is from one CPU and one OS.

To reproduce, build the matching C tool from `tools/` and run the Rust examples:

```sh
cargo run --release --example speed -- tests/vectors.txt
cargo run --release --example speed_sign -- tests/sign_vectors.txt
cargo run --release --example speed_keygen
```

## Next steps

1. Profile keygen by phase, then add AVX2 for the prime NTT, big-integer reduction and fixed-point FFT.
2. Add SSE2 FFT routines to signing and trim the sampling recursion.
3. Measure peak stack use in both C and Rust.
4. Review the secret-handling code for timing leaks.
5. Test and tune on aarch64.

## Layout

```
src/
  lib.rs                public API
  vrfy.rs, codec.rs     verification, decoding
  sign.rs, sign_core.rs signing
  sampler.rs, fpoly.rs  Gaussian sampler, FFT polynomials
  mqsign.rs             arithmetic mod q' = 18433
  kgen.rs, ntru.rs      key generation, NTRU solver
  kpoly.rs, zint31.rs, mp31.rs, fxp.rs   keygen arithmetic
  mq.rs                 arithmetic mod q, NTT, portable + AVX2 kernels
  util.rs, sha3.rs      hash_to_point, mu, SHAKE
  *_tables.rs, tables.rs, primes.rs, keccak_gen.rs   generated
tests/                  interop, signing and keygen tests with vectors
examples/               benchmarks
tools/                  table generators, C vector generators, C benchmarks
```

## License

Public domain (Unlicense), same as the parent repository.
