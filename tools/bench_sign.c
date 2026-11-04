/*
 * C-side signing benchmark; same vectors and method as
 * rust/examples/speed_sign.rs.
 *
 *   clang -O2 -I. rust/tools/bench_sign.c <objs other than speed_fndsa.o> -o bench_sign
 *   ./bench_sign rust/tests/sign_vectors.txt
 *
 * For each standard degree, take the first raw-message vector (key, message
 * and seed). Time ROUNDS rounds of ITERS seeded signatures each (the seed is
 * changed at every call so that the rejection-sampling variance averages
 * out); report the median and minimum per-call time in nanoseconds.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "fndsa.h"

#define ROUNDS  21
#define ITERS   200

static size_t
unhex(const char *s, unsigned char *out)
{
	if (s[0] == '-' && s[1] == 0) {
		return 0;
	}
	size_t n = strlen(s) / 2;
	for (size_t i = 0; i < n; i ++) {
		unsigned v;
		sscanf(s + 2 * i, "%2x", &v);
		out[i] = (unsigned char)v;
	}
	return n;
}

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
main(int argc, char **argv)
{
	if (argc < 2) {
		return 1;
	}
	FILE *f = fopen(argv[1], "r");
	if (!f) {
		return 1;
	}
	static char line[1 << 17];
	int done9 = 0, done10 = 0;
	while (fgets(line, sizeof line, f)) {
		char *tok[9];
		char *p = strtok(line, " \n");
		int k = 0;
		while (p && k < 9) {
			tok[k ++] = p;
			p = strtok(NULL, " \n");
		}
		if (k != 9) {
			continue;
		}
		unsigned logn = (unsigned)atoi(tok[0]);
		if ((logn != 9 || done9) && (logn != 10 || done10)) {
			continue;
		}
		static unsigned char id[64], ctx[256], hv[256], seed[64];
		static unsigned char sk[4096], vk[4096], sig[4096];
		size_t id_len = unhex(tok[2], id);
		size_t ctx_len = unhex(tok[3], ctx);
		size_t hv_len = unhex(tok[4], hv);
		size_t seed_len = unhex(tok[5], seed);
		size_t sk_len = unhex(tok[6], sk);
		(void)unhex(tok[7], vk);
		if (id_len != 0) {
			continue;
		}
		id[0] = 0;
		if (logn == 9) done9 = 1; else done10 = 1;

		double r[ROUNDS];
		int bad = 0;
		unsigned ctr = 0;
		for (int round = -2; round < ROUNDS; round ++) {
			double t0 = now_ns();
			for (int i = 0; i < ITERS; i ++) {
				seed[0] = (unsigned char)ctr;
				seed[1] = (unsigned char)(ctr >> 8);
				ctr ++;
				size_t r2 = fndsa_sign_seeded(sk, sk_len,
					ctx_len ? ctx : NULL, ctx_len,
					(const char *)id, hv, hv_len,
					seed, seed_len,
					sig, FNDSA_SIGNATURE_SIZE(logn));
				bad += (r2 != FNDSA_SIGNATURE_SIZE(logn));
			}
			double t1 = now_ns();
			if (round >= 0) {
				r[round] = (t1 - t0) / ITERS;
			}
		}
		qsort(r, ROUNDS, sizeof r[0], cmp_d);
		printf("C    sign n=%-5u median %10.1f ns   min %10.1f ns%s\n",
			1u << logn, r[ROUNDS / 2], r[0], bad ? "  (FAILED!)" : "");
	}
	return 0;
}
