#include "aes256ctr.h"
#include "rng.h"
#include <stdint.h>
#include <stdio.h>

static inline uint16_t red_plantard_p1(uint32_t a) {
	uint32_t c;
	c = a * 226050911;
	return (((c >> 16) + 1) * 19) >> 16;
}

/* return (x < y) at the most significant bit */
static inline uint64_t ct_lt(const uint64_t x, const uint64_t y) {
	return x - y;
}

/* Sample ZZ_p, p=19 */
static const uint64_t sampler_zzp_check = 247;
void uni_sampler_zzp_1d_out64(uint64_t *out, const size_t len, aes256ctr_ctx *ctx) {
	size_t i;
	uint8_t r[AES256CTR_BLOCKBYTES];
	uint8_t *r_head = NULL;
	size_t rem = 0;
	uint8_t x;
	for (i = 0; i < len; i++) {
		do {
			if ((ctx) && (rem < 1)) {
				aes256ctr_squeezeblocks(r, 1, ctx);
				rem = AES256CTR_BLOCKBYTES;
				r_head = r;
			}
			x = *(r_head++);
			rem--;
		} while (1 ^ (ct_lt(x, sampler_zzp_check) >> 63));
		out[i] = red_plantard_p1(x);
	}
}

/* Sample ZZ_p, p=19 */
void uni_sampler_zzp_1d(uint16_t *out, const size_t len, aes256ctr_ctx *ctx) {
	size_t i;
	uint8_t r[AES256CTR_BLOCKBYTES];
	uint8_t *r_head = NULL;
	size_t rem = 0;
	uint8_t x;
	for (i = 0; i < len; i++) {
		do {
			if ((ctx) && (rem < 1)) {
				aes256ctr_squeezeblocks(r, 1, ctx);
				rem = AES256CTR_BLOCKBYTES;
				r_head = r;
			}
			x = *(r_head++);
			rem--;
		} while (1 ^ (ct_lt(x, sampler_zzp_check) >> 63));
		out[i] = red_plantard_p1(x);
	}
}

static inline uint64_t load_64(const unsigned char *x) {
	return ((uint64_t)(*x)) | (((uint64_t)(*(x + 1))) << 8) | (((uint64_t)(*(x + 2))) << 16) | (((uint64_t)(*(x + 3))) << 24) | (((uint64_t)(*(x + 4))) << 32) | (((uint64_t)(*(x + 5))) << 40) | (((uint64_t)(*(x + 6))) << 48) | (((uint64_t)(*(x + 7))) << 56);
}

static const uint64_t qbar = 478757428866084737ULL; // NTT friendly prime (Rq = Rp * Rqbar, qbar = 1 (mod 2d)) 
static const __uint128_t red_plantard_qbar_r = (((__uint128_t)2594488098315967366ULL)<<64)|9714240157448392833ULL; // r constant used in plantard modular reduction (mod qbar)
static inline uint64_t red_plantard_qbar(__uint128_t a, __uint128_t b) {
	__uint128_t c;
	c = a * b * red_plantard_qbar_r;
	return (((c >> 64) + 1) * qbar) >> 64;
}

/* Sample ZZ_qbar, qbar=478757428866084737 
 * The output is x*(-2**(-128)) mod qbar */
static const uint64_t sampler_zzqbar_check = 18192782296911220006ULL; // largest multiple of qbar < 2^64 - 1
void uni_sampler_zzqbar_1d(uint64_t *out, const size_t len, aes256ctr_ctx *ctx) {
	size_t i;
	uint8_t r[AES256CTR_BLOCKBYTES];
	uint8_t *r_head;
	size_t rem = 0;
	uint64_t x;
	for (i = 0; i < len; i++) {
		do {
			if (rem < 8) {
				aes256ctr_squeezeblocks(r, 1, ctx);
				rem = AES256CTR_BLOCKBYTES;
				r_head = r;
			}
			x = load_64(r_head);
			r_head += 8;
			rem -= 8;
		} while (x >= sampler_zzqbar_check);
		out[i] = red_plantard_qbar(x, 1);
	}
}

