/*
 * C-side key generation benchmark; same seeds and method as
 * rust/examples/speed_keygen.rs. Seeded key generation is deterministic, so
 * both programs do exactly the same work (including the same retries).
 *
 *   clang -O2 -I. rust/tools/bench_keygen.c <objs other than speed_fndsa.o> -o bench_keygen
 *   ./bench_keygen
 *
 * For each standard degree: ROUNDS rounds, each generating KEYS keys from
 * KEYS fixed seeds; reports the median and minimum per-key time.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "fndsa.h"

#define ROUNDS  9

static double
now_ns(void)
{
#ifdef _WIN32
	static LARGE_INTEGER freq;
	LARGE_INTEGER c;
	if (freq.QuadPart == 0) {
		QueryPerformanceFrequency(&freq);
	}
	QueryPerformanceCounter(&c);
	return (double)c.QuadPart * 1e9 / (double)freq.QuadPart;
#else
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
#endif
}

static int
cmp_d(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

int
main(void)
{
	static unsigned char sk[FNDSA_SIGN_KEY_SIZE(10)];
	static unsigned char vk[FNDSA_VRFY_KEY_SIZE(10)];
	for (unsigned logn = 9; logn <= 10; logn ++) {
		int keys = (logn == 9) ? 40 : 10;
		double r[ROUNDS];
		unsigned char acc = 0;
		for (int round = -1; round < ROUNDS; round ++) {
			double t0 = now_ns();
			for (int k = 0; k < keys; k ++) {
				unsigned char seed[16];
				for (int i = 0; i < 16; i ++) {
					seed[i] = (unsigned char)(k * 13 + i * 5 + logn);
				}
				fndsa_keygen_seeded(logn, seed, sizeof seed, sk, vk);
				acc ^= sk[7] ^ vk[3];
			}
			double t1 = now_ns();
			if (round >= 0) {
				r[round] = (t1 - t0) / keys;
			}
		}
		qsort(r, ROUNDS, sizeof r[0], cmp_d);
		printf("C    keygen n=%-5u median %12.1f ns   min %12.1f ns   (%02x)\n",
			1u << logn, r[ROUNDS / 2], r[0], acc);
	}
	return 0;
}
