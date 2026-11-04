/*
 * C-side verify benchmark using the same vectors and the same method as
 * rust/examples/speed.rs, so the two can be compared directly.
 *
 *   clang -O2 -I. rust/tools/bench_verify.c <objs> -o bench_verify
 *   ./bench_verify rust/tests/vectors.txt
 *
 * Method: for each of the standard degrees (512, 1024) take the first
 * raw-message vector, then time ROUNDS rounds of ITERS verifications each,
 * and report the median and the minimum per-call time in nanoseconds.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "fndsa.h"

#define ROUNDS  31
#define ITERS   20000

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
	static char line[1 << 16];
	int done9 = 0, done10 = 0;
	while (fgets(line, sizeof line, f)) {
		char *tok[7];
		char *p = strtok(line, " \n");
		int k = 0;
		while (p && k < 7) {
			tok[k ++] = p;
			p = strtok(NULL, " \n");
		}
		if (k != 7) {
			continue;
		}
		unsigned logn = (unsigned)atoi(tok[0]);
		if ((logn != 9 || done9) && (logn != 10 || done10)) {
			continue;
		}
		static unsigned char id[64], ctx[256], hv[256], vk[2048], sig[2048];
		size_t id_len = unhex(tok[2], id);
		size_t ctx_len = unhex(tok[3], ctx);
		size_t hv_len = unhex(tok[4], hv);
		size_t vk_len = unhex(tok[5], vk);
		size_t sig_len = unhex(tok[6], sig);
		if (id_len != 0) {
			continue;  /* raw-message vectors only */
		}
		id[id_len] = 0;
		if (logn == 9) done9 = 1; else done10 = 1;

		double r[ROUNDS];
		int bad = 0;
		for (int round = -3; round < ROUNDS; round ++) {
			double t0 = now_ns();
			for (int i = 0; i < ITERS; i ++) {
				bad += !fndsa_verify(sig, sig_len, vk, vk_len,
					ctx_len ? ctx : NULL, ctx_len,
					(const char *)id, hv, hv_len);
			}
			double t1 = now_ns();
			if (round >= 0) {
				r[round] = (t1 - t0) / ITERS;
			}
		}
		qsort(r, ROUNDS, sizeof r[0], cmp_d);
		printf("C    verify n=%-5u median %9.1f ns   min %9.1f ns%s\n",
			1u << logn, r[ROUNDS / 2], r[0], bad ? "  (FAILED!)" : "");
	}
	return 0;
}
