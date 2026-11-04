/*
 * Generates deterministic key generation vectors using the C implementation,
 * for byte-for-byte comparison with the Rust port. Build from the repository
 * root, e.g.:
 *
 *   clang -O2 -I. rust/tools/gen_keygen_vectors.c <objects other than speed> \
 *       -o gen_keygen_vectors
 *   ./gen_keygen_vectors > rust/tests/keygen_vectors.txt
 *
 * Line format (hex): logn seed sk vk
 */
#include <stdio.h>
#include <string.h>
#include "fndsa.h"

static void
hex(const void *p, size_t n)
{
	const unsigned char *b = p;
	for (size_t i = 0; i < n; i ++) {
		printf("%02x", b[i]);
	}
}

int
main(void)
{
	static unsigned char sk[FNDSA_SIGN_KEY_SIZE(10)];
	static unsigned char vk[FNDSA_VRFY_KEY_SIZE(10)];
	for (unsigned logn = 2; logn <= 10; logn ++) {
		int reps = (logn >= 9) ? 6 : (logn >= 7 ? 12 : 25);
		for (int rep = 0; rep < reps; rep ++) {
			unsigned char seed[16];
			for (int i = 0; i < 16; i ++) {
				seed[i] = (unsigned char)(logn * 31 + rep * 7 + i);
			}
			fndsa_keygen_seeded(logn, seed, sizeof seed, sk, vk);
			printf("%u ", logn);
			hex(seed, sizeof seed);
			printf(" ");
			hex(sk, FNDSA_SIGN_KEY_SIZE(logn));
			printf(" ");
			hex(vk, FNDSA_VRFY_KEY_SIZE(logn));
			printf("\n");
		}
	}
	return 0;
}
