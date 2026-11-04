//! Key generation tests.
//!
//! Key generation uses only integer and fixed-point arithmetic, so a seeded
//! run must give exactly the same keys as the C code. `tests/keygen_vectors.txt`
//! was produced by the C implementation (tools/gen_keygen_vectors.c).

fn unhex(s: &str) -> Vec<u8> {
    (0..s.len() / 2).map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap()).collect()
}

struct V {
    logn: u32,
    seed: Vec<u8>,
    sk: Vec<u8>,
    vk: Vec<u8>,
}

fn vectors() -> Vec<V> {
    include_str!("keygen_vectors.txt")
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split(' ').collect();
            assert_eq!(f.len(), 4);
            V { logn: f[0].parse().unwrap(), seed: unhex(f[1]), sk: unhex(f[2]), vk: unhex(f[3]) }
        })
        .collect()
}

#[test]
fn byte_identical_to_c() {
    let vs = vectors();
    assert!(vs.len() > 100);
    for (i, v) in vs.iter().enumerate() {
        let mut sk = vec![0u8; fndsa::sign_key_size(v.logn)];
        let mut vk = vec![0u8; fndsa::vrfy_key_size(v.logn)];
        assert!(fndsa::keygen_seeded(v.logn, &v.seed, Some(&mut sk), Some(&mut vk)));
        assert_eq!(sk, v.sk, "vector {i} (logn={}): signing key differs from C", v.logn);
        assert_eq!(vk, v.vk, "vector {i} (logn={}): verifying key differs from C", v.logn);
    }
}

#[test]
fn only_one_key_requested() {
    let v = vectors().into_iter().find(|v| v.logn == 9).unwrap();
    let mut sk = vec![0u8; fndsa::sign_key_size(9)];
    let mut vk = vec![0u8; fndsa::vrfy_key_size(9)];
    assert!(fndsa::keygen_seeded(9, &v.seed, Some(&mut sk), None));
    assert_eq!(sk, v.sk);
    assert!(fndsa::keygen_seeded(9, &v.seed, None, Some(&mut vk)));
    assert_eq!(vk, v.vk);
}

#[test]
fn generated_keys_sign_and_verify() {
    // Rust keygen -> Rust sign -> Rust verify, for every standard and weak
    // degree, with the system RNG.
    for logn in 2..=10u32 {
        let mut sk = vec![0u8; fndsa::sign_key_size(logn)];
        let mut vk = vec![0u8; fndsa::vrfy_key_size(logn)];
        assert!(fndsa::keygen(logn, Some(&mut sk), Some(&mut vk)), "logn={logn}");
        let mut sig = vec![0u8; fndsa::signature_size(logn)];
        let msg = b"hello from the Rust port";
        let weak = logn <= 8;
        let r = if weak {
            fndsa::sign_weak(&sk, &[], fndsa::HASH_ID_RAW, msg, &mut sig)
        } else {
            fndsa::sign(&sk, &[], fndsa::HASH_ID_RAW, msg, &mut sig)
        };
        assert_eq!(r, Some(sig.len()), "logn={logn}");
        let ok = if weak {
            fndsa::verify_weak(&sig, &vk, &[], fndsa::HASH_ID_RAW, msg)
        } else {
            fndsa::verify(&sig, &vk, &[], fndsa::HASH_ID_RAW, msg)
        };
        assert!(ok, "logn={logn}");
    }
}

#[test]
fn rejects_bad_parameters() {
    let mut sk = vec![0u8; fndsa::sign_key_size(9)];
    assert!(!fndsa::keygen_seeded(1, b"x", Some(&mut sk), None));
    assert!(!fndsa::keygen_seeded(11, b"x", Some(&mut sk), None));
    // Output buffer too small.
    let mut small = vec![0u8; fndsa::sign_key_size(9) - 1];
    assert!(!fndsa::keygen_seeded(9, b"x", Some(&mut small), None));
}
