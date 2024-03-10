//#include <immintrin.h>
#include "common.h"
#include "mathc.h"

#define COMPUTE_16X_SMALL_PROOTS_AVX2	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i tmp3;	\
				__m128i tmp4;	\
				__m128i zeros; \
				zeros = _mm_xor_si128(zeros, zeros); \
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					tmp1 = root1s; \
					tmp2 = root2s; \
					tmp3 = root1s; \
					tmp4 = root2s; \
					root1s = _mm_sub_epi16(root1s, ptrs); \
					root2s = _mm_sub_epi16(root2s, ptrs); \
					tmp1 = _mm_subs_epu16(tmp1, ptrs); \
					tmp2 = _mm_subs_epu16(tmp2, ptrs); \
					tmp3 = _mm_cmpeq_epi16(tmp3,ptrs);	\
					tmp4 = _mm_cmpeq_epi16(tmp4,ptrs);	\
					tmp1 = _mm_cmpeq_epi16(tmp1,zeros);	\
					tmp2 = _mm_cmpeq_epi16(tmp2,zeros);	\
					tmp3 = _mm_andnot_si128(tmp3, tmp1);	\
					tmp4 = _mm_andnot_si128(tmp4, tmp2);	\
					tmp3 = _mm_and_si128(tmp3, primes);	\
					tmp4 = _mm_and_si128(tmp4, primes);	\
					root1s = _mm_add_epi16(root1s, tmp3);	\
					root2s = _mm_add_epi16(root2s, tmp4);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epu16(tmp1, root1s); \
					root2s = _mm_min_epu16(root2s, root1s); \
					tmp2 = primes; \
					_mm_store_si128((__m128i *)(h.fbp1 + j), root2s);	\
					primes = _mm_sub_epi16(primes, root2s); \
					_mm_store_si128((__m128i *)(h.first_r1 + j), root2s);	\
					_mm_store_si128((__m128i *)(h.fbp2 + j), tmp1);	\
					tmp2 = _mm_sub_epi16(tmp2, tmp1); \
					_mm_store_si128((__m128i *)(h.first_r2 + j), tmp1);	\
					_mm_store_si128((__m128i *)(h.fbn1 + j), tmp2);	\
					_mm_store_si128((__m128i *)(h.fbn2 + j), primes);	\
				}	\
			} while(0);

#define COMPUTE_16X_SMALL_NROOTS_AVX2	\
		do {	\
				__m128i primes;	\
				__m128i root1s;	\
				__m128i root2s;	\
				__m128i ptrs;	\
				__m128i tmp1;	\
				__m128i tmp2;	\
				__m128i tmp3;	\
				__m128i tmp4;	\
				__m128i zeros; \
				zeros = _mm_xor_si128(zeros, zeros); \
				for (j = h.start; j < h.stop; j += 8) \
				{	\
					ptrs = _mm_load_si128((__m128i *)(h.updates + j)); \
					root1s = _mm_load_si128((__m128i *)(h.first_r1 + j)); \
					root2s = _mm_load_si128((__m128i *)(h.first_r2 + j)); \
					primes = _mm_load_si128((__m128i *)(h.primes + j)); \
					tmp3 = root1s; \
					tmp4 = root2s; \
					root1s = _mm_add_epi16(root1s, ptrs); \
					root2s = _mm_add_epi16(root2s, ptrs); \
					tmp3 = _mm_adds_epu16(tmp3, ptrs); \
					tmp4 = _mm_adds_epu16(tmp4, ptrs); \
					tmp1 = primes; \
					tmp2 = primes; \
					tmp1 = _mm_subs_epu16(tmp1, tmp3); \
					tmp2 = _mm_subs_epu16(tmp2, tmp4); \
					tmp1 = _mm_cmpeq_epi16(tmp1,zeros);	\
					tmp2 = _mm_cmpeq_epi16(tmp2,zeros);	\
					tmp1 = _mm_and_si128(tmp1, primes);	\
					tmp2 = _mm_and_si128(tmp2, primes);	\
					root1s = _mm_sub_epi16(root1s, tmp1);	\
					root2s = _mm_sub_epi16(root2s, tmp2);	\
					tmp1 = root2s;	\
					tmp1 = _mm_max_epu16(tmp1, root1s); \
					root2s = _mm_min_epu16(root2s, root1s); \
					tmp2 = primes; \
					_mm_store_si128((__m128i *)(h.fbp1 + j), root2s);	\
					primes = _mm_sub_epi16(primes, root2s); \
					_mm_store_si128((__m128i *)(h.first_r1 + j), root2s);	\
					_mm_store_si128((__m128i *)(h.fbp2 + j), tmp1);	\
					tmp2 = _mm_sub_epi16(tmp2, tmp1); \
					_mm_store_si128((__m128i *)(h.first_r2 + j), tmp1);	\
					_mm_store_si128((__m128i *)(h.fbn1 + j), tmp2);	\
					_mm_store_si128((__m128i *)(h.fbn2 + j), primes);	\
				}	\
			} while(0);




#ifdef USE_AVX512BW
#include <immintrin.h>

#define COMPUTE_32X_SMALL_PROOTS_AVX512	\
		do {	\
				__m512i primes;	\
				__m512i root1s;	\
				__m512i root2s;	\
				__m512i ptrs;	\
				__m512i tmp1;	\
				__m512i tmp2;	\
				__m512i tmp3;	\
				__m512i tmp4;	\
                __mmask32 m1, m2; \
				for (j = h.start; j < h.stop - 32; j += 32) \
				{	\
					ptrs = _mm512_loadu_si512((__m512i *)(h.updates + j)); \
					root1s = _mm512_loadu_si512((__m512i *)(h.first_r1 + j)); \
					root2s = _mm512_loadu_si512((__m512i *)(h.first_r2 + j)); \
					primes = _mm512_loadu_si512((__m512i *)(h.primes + j)); \
					m1 = _mm512_cmpgt_epu16_mask(ptrs, root1s);	\
					m2 = _mm512_cmpgt_epu16_mask(ptrs, root2s);	\
					root1s = _mm512_sub_epi16(root1s, ptrs); \
					root2s = _mm512_sub_epi16(root2s, ptrs); \
					root1s = _mm512_mask_add_epi16(root1s, m1, root1s, primes);	\
					root2s = _mm512_mask_add_epi16(root2s, m2, root2s, primes);	\
					tmp1 = _mm512_max_epu16(root1s, root2s); \
					tmp2 = _mm512_min_epu16(root1s, root2s); \
                    tmp3 = _mm512_sub_epi16(primes, tmp1); \
                    tmp4 = _mm512_sub_epi16(primes, tmp2); \
					_mm512_storeu_si512((__m512i *)(h.fbp1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.first_r1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.fbp2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.first_r2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.fbn1 + j), tmp3);	\
					_mm512_storeu_si512((__m512i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);

#define COMPUTE_32X_SMALL_NROOTS_AVX512	\
		do {	\
				__m512i primes;	\
				__m512i root1s;	\
				__m512i root2s;	\
				__m512i ptrs;	\
				__m512i tmp1;	\
				__m512i tmp2;	\
				__m512i tmp3;	\
				__m512i tmp4;	\
                __mmask32 m1, m2; \
				for (j = h.start; j < h.stop - 32; j += 32) \
				{	\
					ptrs = _mm512_loadu_si512((__m512i *)(h.updates + j)); \
					root1s = _mm512_loadu_si512((__m512i *)(h.first_r1 + j)); \
					root2s = _mm512_loadu_si512((__m512i *)(h.first_r2 + j)); \
					primes = _mm512_loadu_si512((__m512i *)(h.primes + j)); \
					tmp1 = _mm512_add_epi16(root1s, ptrs); \
					tmp2 = _mm512_add_epi16(root2s, ptrs); \
                    m1 = _mm512_cmplt_epu16_mask(tmp1, root1s);	\
					m2 = _mm512_cmplt_epu16_mask(tmp2, root2s);	\
                    m1 |= _mm512_cmpge_epu16_mask(tmp1, primes);	\
					m2 |= _mm512_cmpge_epu16_mask(tmp2, primes);	\
					root1s = _mm512_mask_sub_epi16(tmp1, m1, tmp1, primes);	\
					root2s = _mm512_mask_sub_epi16(tmp2, m2, tmp2, primes);	\
					tmp1 = _mm512_max_epu16(root1s, root2s); \
					tmp2 = _mm512_min_epu16(root1s, root2s); \
                    tmp3 = _mm512_sub_epi16(primes, tmp1); \
                    tmp4 = _mm512_sub_epi16(primes, tmp2); \
					_mm512_storeu_si512((__m512i *)(h.fbp1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.first_r1 + j), tmp2);	\
					_mm512_storeu_si512((__m512i *)(h.fbp2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.first_r2 + j), tmp1);	\
					_mm512_storeu_si512((__m512i *)(h.fbn1 + j), tmp3);	\
					_mm512_storeu_si512((__m512i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);

#define COMPUTE_16X_SMALL_PROOTS_AVX512	\
		do {	\
				__m256i primes;	\
				__m256i root1s;	\
				__m256i root2s;	\
				__m256i ptrs;	\
				__m256i tmp1;	\
				__m256i tmp2;	\
				__m256i tmp3;	\
				__m256i tmp4;	\
                __mmask16 m1, m2; \
				for (j = h.start; j < h.stop; j += 16) \
				{	\
					ptrs = _mm256_load_si256((__m256i *)(h.updates + j)); \
					root1s = _mm256_load_si256((__m256i *)(h.first_r1 + j)); \
					root2s = _mm256_load_si256((__m256i *)(h.first_r2 + j)); \
					primes = _mm256_load_si256((__m256i *)(h.primes + j)); \
					m1 = _mm256_cmpgt_epu16_mask(ptrs, root1s);	\
					m2 = _mm256_cmpgt_epu16_mask(ptrs, root2s);	\
					root1s = _mm256_sub_epi16(root1s, ptrs); \
					root2s = _mm256_sub_epi16(root2s, ptrs); \
					root1s = _mm256_mask_add_epi16(root1s, m1, root1s, primes);	\
					root2s = _mm256_mask_add_epi16(root2s, m2, root2s, primes);	\
					tmp1 = _mm256_max_epu16(root1s, root2s); \
					tmp2 = _mm256_min_epu16(root1s, root2s); \
                    tmp3 = _mm256_sub_epi16(primes, tmp1); \
                    tmp4 = _mm256_sub_epi16(primes, tmp2); \
					_mm256_store_si256((__m256i *)(h.fbp1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.first_r1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.fbp2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.first_r2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.fbn1 + j), tmp3);	\
					_mm256_store_si256((__m256i *)(h.fbn2 + j), tmp4);	\
				}	\
			} while(0);

#define COMPUTE_16X_SMALL_NROOTS_AVX512	\
		do {	\
				__m256i primes;	\
				__m256i root1s;	\
				__m256i root2s;	\
				__m256i ptrs;	\
				__m256i tmp1;	\
				__m256i tmp2;	\
				__m256i tmp3;	\
				__m256i tmp4;	\
                __mmask16 m1, m2; \
				for (j = h.start; j < h.stop; j += 16) \
				{	\
					ptrs = _mm256_load_si256((__m256i *)(h.updates + j)); \
					root1s = _mm256_load_si256((__m256i *)(h.first_r1 + j)); \
					root2s = _mm256_load_si256((__m256i *)(h.first_r2 + j)); \
					primes = _mm256_load_si256((__m256i *)(h.primes + j)); \
					tmp1 = _mm256_add_epi16(root1s, ptrs); \
					tmp2 = _mm256_add_epi16(root2s, ptrs); \
                    m1 = _mm256_cmplt_epu16_mask(tmp1, root1s);	\
					m2 = _mm256_cmplt_epu16_mask(tmp2, root2s);	\
                    m1 |= _mm256_cmpge_epu16_mask(tmp1, primes);	\
					m2 |= _mm256_cmpge_epu16_mask(tmp2, primes);	\
					root1s = _mm256_mask_sub_epi16(tmp1, m1, tmp1, primes);	\
					root2s = _mm256_mask_sub_epi16(tmp2, m2, tmp2, primes);	\
					tmp1 = _mm256_max_epu16(root1s, root2s); \
					tmp2 = _mm256_min_epu16(root1s, root2s); \
                    tmp3 = _mm256_sub_epi16(primes, tmp1); \
                    tmp4 = _mm256_sub_epi16(primes, tmp2); \
					_mm256_store_si256((__m256i *)(h.fbp1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.first_r1 + j), tmp2);	\
					_mm256_store_si256((__m256i *)(h.fbp2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.first_r2 + j), tmp1);	\
					_mm256_store_si256((__m256i *)(h.fbn1 + j), tmp3);	\
					_mm256_store_si256((__m256i *)(h.fbn2 + j), tmp4);	\
				}	\
                h.stop = j - 32; \
			} while(0);

#endif
