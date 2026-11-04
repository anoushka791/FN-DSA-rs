//! Verifies signatures produced by the C implementation (see
//! tools/gen_vectors.c), with both kernel sets, and checks that tampering is
//! rejected.

use fndsa::internals::{avx2_kernels, verify_with, Kernels, PORTABLE};

struct Vector {
    logn: u32,
    weak: bool,
    id: Vec<u8>,
    ctx: Vec<u8>,
    hv: Vec<u8>,
    vk: Vec<u8>,
    sig: Vec<u8>,
}

fn unhex(s: &str) -> Vec<u8> {
    if s == "-" {
        return Vec::new();
    }
    (0..s.len() / 2).map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap()).collect()
}

fn vectors() -> Vec<Vector> {
    let text = include_str!("vectors.txt");
    text.lines()
        .map(|l| {
            let f: Vec<&str> = l.split(' ').collect();
            assert_eq!(f.len(), 7);
            Vector {
                logn: f[0].parse().unwrap(),
                weak: f[1] == "1",
                id: unhex(f[2]),
                ctx: unhex(f[3]),
                hv: unhex(f[4]),
                vk: unhex(f[5]),
                sig: unhex(f[6]),
            }
        })
        .collect()
}

fn kernel_sets() -> Vec<(&'static str, &'static Kernels)> {
    let mut v = vec![("portable", &PORTABLE)];
    if let Some(k) = avx2_kernels() {
        v.push(("avx2", k));
    }
    v
}

fn check(k: &Kernels, v: &Vector, sig: &[u8], vk: &[u8], ctx: &[u8], hv: &[u8]) -> bool {
    verify_with(k, v.weak, sig, vk, ctx, &v.id, hv)
}

#[test]
fn accepts_c_signatures() {
    let vs = vectors();
    assert!(vs.len() > 90);
    for (name, k) in kernel_sets() {
        for v in &vs {
            assert_eq!(v.sig.len(), fndsa::signature_size(v.logn));
            assert_eq!(v.vk.len(), fndsa::vrfy_key_size(v.logn));
            assert!(check(k, v, &v.sig, &v.vk, &v.ctx, &v.hv), "{name}: logn={} rejected", v.logn);
        }
    }
}

#[test]
fn public_api_matches() {
    for v in vectors() {
        let ok = if v.weak {
            fndsa::verify_weak(&v.sig, &v.vk, &v.ctx, &v.id, &v.hv)
        } else {
            fndsa::verify(&v.sig, &v.vk, &v.ctx, &v.id, &v.hv)
        };
        assert!(ok, "logn={}", v.logn);
        // The wrong entry point must reject (degree range enforcement).
        let wrong = if v.weak {
            fndsa::verify(&v.sig, &v.vk, &v.ctx, &v.id, &v.hv)
        } else {
            fndsa::verify_weak(&v.sig, &v.vk, &v.ctx, &v.id, &v.hv)
        };
        assert!(!wrong);
    }
}

#[test]
fn rejects_tampering() {
    let vs = vectors();
    for (name, k) in kernel_sets() {
        for v in &vs {
            // Flip every bit of the signature (header, nonce, body) for a
            // sample of vectors: none may verify. (Bit flips in the padding
            // are rejected by the canonical-encoding check.)
            let step = if v.logn >= 9 { 61 } else { 1 };
            for bit in (0..v.sig.len() * 8).step_by(step) {
                let mut s = v.sig.clone();
                s[bit >> 3] ^= 1 << (bit & 7);
                assert!(
                    !check(k, v, &s, &v.vk, &v.ctx, &v.hv),
                    "{name}: logn={} accepted sig with bit {bit} flipped",
                    v.logn
                );
            }
            // Tampered message.
            if !v.hv.is_empty() {
                let mut hv = v.hv.clone();
                hv[0] ^= 1;
                assert!(!check(k, v, &v.sig, &v.vk, &v.ctx, &hv), "{name}: bad message accepted");
            }
            // Tampered context (only meaningful when a context is used).
            if !v.ctx.is_empty() {
                let mut ctx = v.ctx.clone();
                ctx[0] ^= 1;
                assert!(!check(k, v, &v.sig, &v.vk, &ctx, &v.hv), "{name}: bad ctx accepted");
            }
            // Tampered key.
            let mut vk = v.vk.clone();
            let last = vk.len() - 1;
            vk[last] ^= 1;
            assert!(!check(k, v, &v.sig, &vk, &v.ctx, &v.hv), "{name}: bad key accepted");
            // Truncated / extended inputs.
            assert!(!check(k, v, &v.sig[..v.sig.len() - 1], &v.vk, &v.ctx, &v.hv));
            let mut longer = v.sig.clone();
            longer.push(0);
            assert!(!check(k, v, &longer, &v.vk, &v.ctx, &v.hv));
            assert!(!check(k, v, &[], &v.vk, &v.ctx, &v.hv));
            assert!(!check(k, v, &v.sig, &[], &v.ctx, &v.hv));
        }
    }
}

#[test]
fn rejects_bad_parameters() {
    let v = vectors().into_iter().find(|v| v.logn == 9).unwrap();
    let k = &PORTABLE;
    // Context longer than 255 bytes.
    assert!(!check(k, &v, &v.sig, &v.vk, &[0u8; 256], &v.hv));
    // Unknown identifier, and a truncated OID.
    let bad_id = Vector { id: vec![0x07, 0x00], ..clone_vec(&v) };
    assert!(!verify_with(k, false, &v.sig, &v.vk, &[], &bad_id.id, &v.hv));
    assert!(!verify_with(k, false, &v.sig, &v.vk, &[], &[0x06, 0x09, 0x60], &v.hv));
    // External mu with the wrong length.
    assert!(!verify_with(k, false, &v.sig, &v.vk, &[], fndsa::HASH_ID_EXTMU, &[0u8; 63]));
}

fn clone_vec(v: &Vector) -> Vector {
    Vector {
        logn: v.logn,
        weak: v.weak,
        id: v.id.clone(),
        ctx: v.ctx.clone(),
        hv: v.hv.clone(),
        vk: v.vk.clone(),
        sig: v.sig.clone(),
    }
}

/// The AVX2 and portable kernels must agree exactly on random data.
#[test]
fn kernels_agree() {
    let Some(a) = avx2_kernels() else {
        eprintln!("AVX2 not available; skipping");
        return;
    };
    let p = &PORTABLE;
    let mut x: u64 = 0x1234_5678_9abc_def0;
    let mut rnd = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    for logn in 2..=10u32 {
        let n = 1usize << logn;
        for _ in 0..200 {
            let ext: Vec<u16> = (0..n).map(|_| (rnd() % 12289) as u16).collect();
            let sgn: Vec<u16> = (0..n).map(|_| ((rnd() % 4095) as i32 - 2047) as i16 as u16).collect();
            let other: Vec<u16> = (0..n).map(|_| (1 + rnd() % 12289) as u16).collect();
            // In-place unary kernels. Inputs are in the domain each expects.
            macro_rules! unary {
                ($f:ident, $inp:expr) => {{
                    let mut u = $inp.clone();
                    let mut w = $inp.clone();
                    (p.$f)(logn, &mut u);
                    (a.$f)(logn, &mut w);
                    assert_eq!(u, w, "{} logn={logn}", stringify!($f));
                }};
            }
            unary!(ext_to_int, ext);
            unary!(int_to_ext, other);
            unary!(signed_to_int, sgn);
            unary!(int_to_ntt, other);
            unary!(ntt_to_int, other);
            let (mut u, mut w) = (other.clone(), other.clone());
            (p.mul_ntt)(logn, &mut u, &ext.iter().map(|&e| e.max(1)).collect::<Vec<_>>());
            (a.mul_ntt)(logn, &mut w, &ext.iter().map(|&e| e.max(1)).collect::<Vec<_>>());
            assert_eq!(u, w, "mul_ntt logn={logn}");
            let (mut u, mut w) = (other.clone(), other.clone());
            (p.sub)(logn, &mut u, &other.iter().rev().cloned().collect::<Vec<_>>());
            (a.sub)(logn, &mut w, &other.iter().rev().cloned().collect::<Vec<_>>());
            assert_eq!(u, w, "sub logn={logn}");
            assert_eq!((p.sqnorm_signed)(logn, &sgn), (a.sqnorm_signed)(logn, &sgn));
            assert_eq!((p.sqnorm_binf_ext)(logn, &ext), (a.sqnorm_binf_ext)(logn, &ext));
            // Small coefficients (inside B_INF) exercise the non-saturated path.
            let small: Vec<u16> = (0..n).map(|_| ((rnd() % 1201) as i32 - 600).rem_euclid(12289) as u16).collect();
            assert_eq!((p.sqnorm_binf_ext)(logn, &small), (a.sqnorm_binf_ext)(logn, &small));
        }
    }
}

/// NTT followed by inverse NTT is the identity (internal representation).
#[test]
fn ntt_roundtrip() {
    let mut x: u64 = 99;
    for (name, k) in kernel_sets() {
        for logn in 2..=10u32 {
            let n = 1usize << logn;
            let orig: Vec<u16> = (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (1 + (x >> 33) % 12289) as u16
                })
                .collect();
            let mut d = orig.clone();
            (k.int_to_ntt)(logn, &mut d);
            (k.ntt_to_int)(logn, &mut d);
            assert_eq!(d, orig, "{name} logn={logn}");
        }
    }
}
