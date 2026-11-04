/*
 * Generates deterministic signing vectors (key, seed, message, signature)
 * using the C implementation, for byte-for-byte comparison with the Rust
 * port. Build from the repository root, e.g.:
 *
 *   clang -O2 -I. rust/tools/gen_sign_vectors.c <objects other than speed> \
 *       -o gen_sign_vectors
 *   ./gen_sign_vectors > rust/tests/sign_vectors.txt
 *
 * Line format (all fields hex, "-" for empty):
 *   logn weak id ctx hv seed sk vk sig
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
		int reps = (logn >= 9) ? 3 : 6;
		for (int rep = 0; rep < reps; rep ++) {
			memset(seed, (int)(logn * 16 + rep + 1), sizeof seed);
			fndsa_keygen_seeded(logn, seed, sizeof seed, sk, vk);
			size_t sklen = FNDSA_SIGN_KEY_SIZE(logn);
			unsigned char msg[48];
			for (int i = 0; i < 48; i ++) {
				msg[i] = (unsigned char)(i * 11 + rep);
			}
			for (int variant = 0; variant < 3; variant ++) {
				seed[0] ^= 0xA5;
				seed[1] += (unsigned char)(variant + 1);
				const char *id;
				size_t id_len;
				const unsigned char *ctx = (const unsigned char *)"sign-ctx";
				size_t ctx_len;
				const unsigned char *hv = msg;
				size_t hv_len;
				unsigned char mu[64];
				switch (variant) {
				case 0:
					id = ""; id_len = 0; ctx_len = 0;
					hv_len = 13;
					break;
				case 1:
					id = sha256; id_len = 11; ctx_len = 8;
					hv_len = 32;
					break;
				default: {
					unsigned char hvk[64];
					fndsa_hashed_vrfykey_from_vrfykey(hvk, vk,
						FNDSA_VRFY_KEY_SIZE(logn));
					fndsa_compute_mu(mu, hvk, ctx, 8,
						FNDSA_HASH_ID_RAW, msg, 21);
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
				printf("%u %d ", logn, weak);
				hex(id, id_len); printf(" ");
				hex(ctx, ctx_len); printf(" ");
				hex(hv, hv_len); printf(" ");
				hex(seed, sizeof seed); printf(" ");
				hex(sk, sklen); printf(" ");
				hex(vk, FNDSA_VRFY_KEY_SIZE(logn)); printf(" ");
				hex(sig, r); printf("\n");
			}
		}
	}
	return 0;
}
