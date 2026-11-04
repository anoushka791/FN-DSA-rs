//! Key generation benchmark; same seeds and method as tools/bench_keygen.c.
//!
//!   cargo run --release --example speed_keygen

use std::hint::black_box;
use std::time::Instant;

const ROUNDS: usize = 9;

fn main() {
    for logn in 9u32..=10 {
        let keys = if logn == 9 { 40 } else { 10 };
        let mut sk = vec![0u8; fndsa::sign_key_size(logn)];
        let mut vk = vec![0u8; fndsa::vrfy_key_size(logn)];
        let mut r = Vec::with_capacity(ROUNDS);
        let mut acc = 0u8;
        for round in -1i32..ROUNDS as i32 {
            let t0 = Instant::now();
            for k in 0..keys {
                let mut seed = [0u8; 16];
                for (i, s) in seed.iter_mut().enumerate() {
                    *s = (k * 13 + i * 5 + logn as usize) as u8;
                }
                assert!(fndsa::keygen_seeded(
                    logn,
                    black_box(&seed),
                    Some(black_box(&mut sk)),
                    Some(black_box(&mut vk))
                ));
                acc ^= sk[7] ^ vk[3];
            }
            let dt = t0.elapsed().as_nanos() as f64 / keys as f64;
            if round >= 0 {
                r.push(dt);
            }
        }
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "Rust keygen n={:<5} median {:12.1} ns   min {:12.1} ns   ({:02x})",
            1u32 << logn,
            r[ROUNDS / 2],
            r[0],
            acc
        );
    }
}
