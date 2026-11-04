//! Verification benchmark; same method as tools/bench_verify.c.
//!
//!   cargo run --release --example speed -- tests/vectors.txt [portable]
//!
//! For each standard degree, takes the first raw-message vector and times
//! ROUNDS rounds of ITERS verifications; prints median and minimum per-call
//! time in nanoseconds.

use fndsa::internals::{avx2_kernels, verify_with, Kernels, PORTABLE};
use std::hint::black_box;
use std::time::Instant;

const ROUNDS: usize = 31;
const ITERS: usize = 20000;

fn unhex(s: &str) -> Vec<u8> {
    if s == "-" {
        return Vec::new();
    }
    (0..s.len() / 2).map(|i| u8::from_str_radix(&s[2 * i..2 * i + 2], 16).unwrap()).collect()
}

fn bench(label: &str, k: &Kernels, logn: u32, ctx: &[u8], hv: &[u8], vk: &[u8], sig: &[u8]) {
    let mut r = Vec::with_capacity(ROUNDS);
    let mut bad = 0usize;
    for round in -3i32..ROUNDS as i32 {
        let t0 = Instant::now();
        for _ in 0..ITERS {
            let ok = verify_with(
                k,
                false,
                black_box(sig),
                black_box(vk),
                black_box(ctx),
                black_box(b""),
                black_box(hv),
            );
            bad += !ok as usize;
        }
        let dt = t0.elapsed().as_nanos() as f64 / ITERS as f64;
        if round >= 0 {
            r.push(dt);
        }
    }
    r.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!(
        "{label:<5}verify n={:<5} median {:9.1} ns   min {:9.1} ns{}",
        1u32 << logn,
        r[ROUNDS / 2],
        r[0],
        if bad != 0 { "  (FAILED!)" } else { "" }
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).map(String::as_str).unwrap_or("tests/vectors.txt");
    let force_portable = args.get(2).map(String::as_str) == Some("portable");
    let k: &Kernels = if force_portable { &PORTABLE } else { avx2_kernels().unwrap_or(&PORTABLE) };
    let label = if std::ptr::eq(k, &PORTABLE) { "Rust-" } else { "Rust" };

    let text = std::fs::read_to_string(path).expect("vectors file");
    let mut done = [false; 2];
    for l in text.lines() {
        let f: Vec<&str> = l.split(' ').collect();
        let logn: u32 = f[0].parse().unwrap();
        if !(9..=10).contains(&logn) || done[(logn - 9) as usize] || f[2] != "-" {
            continue;
        }
        done[(logn - 9) as usize] = true;
        bench(label, k, logn, &unhex(f[3]), &unhex(f[4]), &unhex(f[5]), &unhex(f[6]));
    }
}
