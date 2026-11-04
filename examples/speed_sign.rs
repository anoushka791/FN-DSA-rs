//! Signing benchmark; same method as tools/bench_sign.c.
//!
//!   cargo run --release --example speed_sign -- tests/sign_vectors.txt
//!
//! For each standard degree, takes the first raw-message vector (key,
//! message, seed) and times ROUNDS rounds of ITERS seeded signatures; the
//! seed changes at every call. Prints median and minimum per-call time.

use std::hint::black_box;
use std::time::Instant;

const ROUNDS: usize = 21;
const ITERS: usize = 200;

fn unhex(s: &str) -> Vec<u8> {
    if s == "-" {
        return Vec::new();
    }
    (0..s.len() / 2).map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap()).collect()
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "tests/sign_vectors.txt".into());
    let text = std::fs::read_to_string(path).expect("vectors file");
    let mut done = [false; 2];
    for l in text.lines() {
        let f: Vec<&str> = l.split(' ').collect();
        let logn: u32 = f[0].parse().unwrap();
        if !(9..=10).contains(&logn) || done[(logn - 9) as usize] || f[2] != "-" {
            continue;
        }
        done[(logn - 9) as usize] = true;
        let (ctx, hv, mut seed, sk) = (unhex(f[3]), unhex(f[4]), unhex(f[5]), unhex(f[6]));
        let mut sig = vec![0u8; fndsa::signature_size(logn)];

        let mut r = Vec::with_capacity(ROUNDS);
        let mut bad = 0usize;
        let mut ctr = 0u32;
        for round in -2i32..ROUNDS as i32 {
            let t0 = Instant::now();
            for _ in 0..ITERS {
                seed[0] = ctr as u8;
                seed[1] = (ctr >> 8) as u8;
                ctr += 1;
                let res = fndsa::sign_seeded(
                    black_box(&sk),
                    black_box(&ctx),
                    black_box(b""),
                    black_box(&hv),
                    black_box(&seed),
                    black_box(&mut sig),
                );
                bad += (res != Some(sig.len())) as usize;
            }
            let dt = t0.elapsed().as_nanos() as f64 / ITERS as f64;
            if round >= 0 {
                r.push(dt);
            }
        }
        r.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "Rust sign n={:<5} median {:10.1} ns   min {:10.1} ns{}",
            1u32 << logn,
            r[ROUNDS / 2],
            r[0],
            if bad != 0 { "  (FAILED!)" } else { "" }
        );
    }
}
