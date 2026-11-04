/*
 * Generates interoperability vectors for the Rust port, using the C
 * implementation as signer. Build from the repository root, e.g.:
 *
 *   clang -O2 -I. rust/tools/gen_vectors.c codec.o mq.o sha3.o sysrng.o \
 *       util.o vrfy.o kgen*.o sign*.o -o gen_vectors
 *   ./gen_vectors > rust/tests/vectors.txt
 *
 * Line format (all fields hex, "-" for empty):
 *   logn weak id ctx hv vk sig
 * "id" is the raw hash-identifier bytes (empty = raw message).
 */
#include <stdio.h>
#include <string.h>
#include "fndsa.h"

static void
hex(const void *p, size_t n)
{
	const unsigned char *b = p;
	if (n == 0) {
		printf("-");
		return;
	}
	for (size_t i = 0; i < n; i ++) {
		printf("%02x", b[i]);
	}
}

static void
emit(unsigned logn, int weak, const char *id, size_t id_len,
	const void *ctx, size_t ctx_len, const void *hv, size_t hv_len,
	const void *vk, const void *sig)
{
	printf("%u %d ", logn, weak);
	hex(id, id_len);
	printf(" ");
	hex(ctx, ctx_len);
	printf(" ");
	hex(hv, hv_len);
	printf(" ");
	hex(vk, FNDSA_VRFY_KEY_SIZE(logn));
	printf(" ");
	hex(sig, FNDSA_SIGNATURE_SIZE(logn));
	printf("\n");
}

int
main(void)
{
	static unsigned char sk[FNDSA_SIGN_KEY_SIZE(10)];
	static unsigned char vk[FNDSA_VRFY_KEY_SIZE(10)];
	static unsigned char sig[FNDSA_SIGNATURE_SIZE(10)];
	static const char sha256[] = FNDSA_HASH_ID_SHA256;
	static const char extmu[] = FNDSA_HASH_ID_EXTMU;
	unsigned char seed[40];

	for (unsigned logn = 2; logn <= 10; logn ++) {
		int weak = (logn <= 8);
		int reps = (logn >= 9) ? 3 : 4;
		for (int rep = 0; rep < reps; rep ++) {
			memset(seed, (int)(logn * 16 + rep), sizeof seed);
			fndsa_keygen_seeded(logn, seed, sizeof seed, sk, vk);
			size_t sklen = FNDSA_SIGN_KEY_SIZE(logn);
			unsigned char msg[48];
			for (int i = 0; i < 48; i ++) {
				msg[i] = (unsigned char)(i * 7 + rep);
			}
			for (int variant = 0; variant < 3; variant ++) {
				seed[0] ^= 0x5A;
				seed[1] += (unsigned char)variant;
				const char *id;
				size_t id_len;
				const unsigned char *ctx = (const unsigned char *)"ctx-ABC";
				size_t ctx_len;
				const unsigned char *hv = msg;
				size_t hv_len;
				unsigned char mu[64];
				switch (variant) {
				case 0: /* raw message, no context */
					id = ""; id_len = 0; ctx_len = 0;
					hv_len = 11;
					break;
				case 1: /* pre-hashed (SHA-256 OID), with context */
					id = sha256; id_len = 11; ctx_len = 7;
					hv_len = 32;
					break;
				default: { /* external mu */
					unsigned char hvk[64];
					fndsa_hashed_vrfykey_from_vrfykey(hvk, vk,
						FNDSA_VRFY_KEY_SIZE(logn));
					fndsa_compute_mu(mu, hvk, ctx, 7,
						FNDSA_HASH_ID_RAW, msg, 20);
					id = extmu; id_len = 1; ctx_len = 0;
					hv = mu; hv_len = 64;
					break;
				}
				}
				size_t r = weak
					? fndsa_sign_weak_seeded(sk, sklen,
						ctx_len ? ctx : NULL, ctx_len,
						id, hv, hv_len, seed, sizeof seed,
						sig, FNDSA_SIGNATURE_SIZE(logn))
					: fndsa_sign_seeded(sk, sklen,
						ctx_len ? ctx : NULL, ctx_len,
						id, hv, hv_len, seed, sizeof seed,
						sig, FNDSA_SIGNATURE_SIZE(logn));
				if (r != FNDSA_SIGNATURE_SIZE(logn)) {
					fprintf(stderr, "sign failed logn=%u\n", logn);
					return 1;
				}
				/* Sanity: the C verifier must accept it. */
				int ok = weak
					? fndsa_verify_weak(sig, r, vk,
						FNDSA_VRFY_KEY_SIZE(logn),
						ctx_len ? ctx : NULL, ctx_len,
						id, hv, hv_len)
					: fndsa_verify(sig, r, vk,
						FNDSA_VRFY_KEY_SIZE(logn),
						ctx_len ? ctx : NULL, ctx_len,
						id, hv, hv_len);
				if (!ok) {
					fprintf(stderr, "C verify failed logn=%u\n", logn);
					return 1;
				}
				emit(logn, weak, id, id_len, ctx, ctx_len, hv,
					hv_len, vk, sig);
			}
		}
	}
	return 0;
}
