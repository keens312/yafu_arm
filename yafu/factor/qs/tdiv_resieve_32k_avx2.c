/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may
benefit from your work.

Some parts of the code (and also this header), included in this
distribution have been reused from other sources. In particular I
have benefitted greatly from the work of Jason Papadopoulos's msieve @
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St.
Denis Tom's Fast Math library.  Many thanks to their kind donation of
code to the public domain.
--bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "common.h"

#if defined( USE_AVX2 )

#include "tdiv_macros_common.h"
#include "qs_impl.h"
#include <immintrin.h>


#define INIT_CORRECTIONS \
	corrections[0] = 32768 - block_loc; \
	corrections[1] = 32768 - block_loc; \
	corrections[2] = 32768 - block_loc; \
	corrections[3] = 32768 - block_loc; \
	corrections[4] = 32768 - block_loc; \
	corrections[5] = 32768 - block_loc; \
	corrections[6] = 32768 - block_loc; \
	corrections[7] = 32768 - block_loc; \
	corrections[8] = 32768 - block_loc; \
	corrections[9] = 32768 - block_loc; \
	corrections[10] = 32768 - block_loc; \
	corrections[11] = 32768 - block_loc; \
	corrections[12] = 32768 - block_loc; \
	corrections[13] = 32768 - block_loc; \
	corrections[14] = 32768 - block_loc; \
	corrections[15] = 32768 - block_loc; 

#if defined(GCC_ASM64X)


#ifdef _WIN32
#define ASM_ ASM_M
#else
#define ASM_ ASM_G
#endif

#define TDIV_MED_CLEAN_AVX2 ASM_("vzeroupper   \n\t");

#define STEP_COMPARE_COMBINE \
		"vpsubw     %%xmm1, %%xmm2, %%xmm2 \n\t"		/* subtract primes from root1s */ \
		"vpsubw     %%xmm1, %%xmm3, %%xmm3 \n\t"		/* subtract primes from root2s */ \
		"vpcmpeqw   %%xmm2, %%xmm5, %%xmm5 \n\t"	/* root1s ?= 0 */ \
		"vpcmpeqw   %%xmm3, %%xmm6, %%xmm6 \n\t"	/* root2s ?= 0 (don't need to re-zero) */ \
		"vpor       %%xmm5, %%xmm7, %%xmm7 \n\t"		/* combine results (if == we are done, else) */ \
		"vpor       %%xmm6, %%xmm0, %%xmm0 \n\t"		/* combine results (xmm5/6 remain 0) */

#define STEP_COMPARE_COMBINE_AVX2 \
		"vpsubw     %%ymm1, %%ymm2, %%ymm2 \n\t"		/* subtract primes from root1s */ \
		"vpsubw     %%ymm1, %%ymm3, %%ymm3 \n\t"		/* subtract primes from root2s */ \
		"vpcmpeqw   %%ymm2, %%ymm5, %%ymm5 \n\t"	/* root1s ?= 0 */ \
		"vpcmpeqw   %%ymm3, %%ymm6, %%ymm6 \n\t"	/* root2s ?= 0 (don't need to re-zero) */ \
		"vpor       %%ymm5, %%ymm7, %%ymm7 \n\t"		/* combine results (if == we are done, else) */ \
		"vpor       %%ymm6, %%ymm0, %%ymm0 \n\t"		/* combine results (xmm5/6 remain 0) */

#define INIT_RESIEVE \
		"vmovdqa    (%4),   %%xmm4 \n\t"		/* bring in corrections to roots */				\
		"vpxor      %%xmm0, %%xmm0, %%xmm0 \n\t"		/* zero xmm8 */ \
		"vmovdqa    (%2),   %%xmm2 \n\t"		/* bring in 8 root1s */ \
		"vpaddw     %%xmm4, %%xmm2, %%xmm2 \n\t"		/* correct root1s */ \
		"vmovdqa    (%3),   %%xmm3 \n\t"		/* bring in 8 root2s */ \
		"vpaddw     %%xmm4, %%xmm3, %%xmm3 \n\t"		/* correct root2s */ \
		"vmovdqa    (%1),   %%xmm1 \n\t"		/* bring in 8 primes */ \
		"vpxor      %%xmm7, %%xmm7, %%xmm7 \n\t"		/* zero xmm7 */ \
		"vpxor      %%xmm5, %%xmm5, %%xmm5 \n\t"		/* zero xmm5 */ \
		"vpxor      %%xmm6, %%xmm6, %%xmm6 \n\t"		/* zero xmm6 */

#define INIT_RESIEVE_AVX2 \
		"vmovdqa    (%4),   %%ymm4 \n\t"		/* bring in corrections to roots */				\
		"vpxor      %%ymm0, %%ymm0, %%ymm0 \n\t"		/* zero xmm8 */ \
		"vmovdqa    (%2),   %%ymm2 \n\t"		/* bring in 8 root1s */ \
		"vpaddw     %%ymm4, %%ymm2, %%ymm2 \n\t"		/* correct root1s */ \
		"vmovdqa    (%3),   %%ymm3 \n\t"		/* bring in 8 root2s */ \
		"vpaddw     %%ymm4, %%ymm3, %%ymm3 \n\t"		/* correct root2s */ \
		"vmovdqa    (%1),   %%ymm1 \n\t"		/* bring in 8 primes */ \
		"vpxor      %%ymm7, %%ymm7, %%ymm7 \n\t"		/* zero xmm7 */ \
		"vpxor      %%ymm5, %%ymm5, %%ymm5 \n\t"		/* zero xmm5 */ \
		"vpxor      %%ymm6, %%ymm6, %%ymm6 \n\t"		/* zero xmm6 */

#define INIT_RESIEVE_2(offset) \
		"vmovdqa    " offset "(%4),   %%xmm4 \n\t"		/* bring in corrections to roots */				\
		"vpxor      %%xmm0, %%xmm0, %%xmm0 \n\t"		/* zero xmm8 */ \
		"vmovdqa    " offset "(%2),   %%xmm2 \n\t"		/* bring in 8 root1s */ \
		"vpaddw     %%xmm4, %%xmm2, %%xmm2 \n\t"		/* correct root1s */ \
		"vmovdqa    " offset "(%3),   %%xmm3 \n\t"		/* bring in 8 root2s */ \
		"vpaddw     %%xmm4, %%xmm3, %%xmm3 \n\t"		/* correct root2s */ \
		"vmovdqa    " offset "(%1),   %%xmm1 \n\t"		/* bring in 8 primes */ \
		"vpxor      %%xmm7, %%xmm7, %%xmm7 \n\t"		/* zero xmm7 */ \
		"vpxor      %%xmm5, %%xmm5, %%xmm5 \n\t"		/* zero xmm5 */ \
		"vpxor      %%xmm6, %%xmm6, %%xmm6 \n\t"		/* zero xmm6 */



#define GATHER_BIT_INDICES \
        "andl   $0xaaaaaaaa,%%r8d   \n\t"   /* mask the bits we don't care about */ \
        "movl	%0,%%r11d		\n\t"		/* initialize count of set bits */ \
        "xorq	%%r10,%%r10		\n\t"		/* initialize bit scan offset */ \
        "1:			\n\t"					/* top of bit scan loop */ \
        "bsfl	%%r8d,%%ecx		\n\t"		/* put least significant set bit index into rcx */ \
        "jz 2f	\n\t"						/* jump out if zero (no hits).  high percentage. */ \
        "addl	%%ecx,%%r10d	\n\t"			/* add in the offset of this index */ \
        "movl   %%r10d,%%r9d \n\t" \
        "shrl   $1,%%r9d \n\t"   \
        "addl   %6,%%r9d \n\t"   \
        "movw	%%r9w, (%5, %%r11, 2) \n\t"		/* put the bit index into the output buffer */ \
        "shrq	%%cl,%%r8	\n\t"			/* shift the bit scan register up to the bit we just processed */ \
        "incl	%%r11d		\n\t"			/* increment the count of set bits */ \
        "incq	%%r10		\n\t"			/* increment the index */ \
        "shrq	$1, %%r8 \n\t"				/* clear the bit */ \
        "jmp 1b		\n\t"					/* loop if so */ \
        "2:		\n\t"						/*  */ \
        "movl	%%r11d, %0 \n\t"			/* return the count of set bits */ \

#define RESIEVE_8X_14BIT_MAX \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_8X_15BIT_MAX \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_8X_16BIT_MAX \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %0 \n\t"		/* if one of these primes divides this location, this will be !0*/ \
			: "=r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections) \
			: "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_8X_14BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_8X_15BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_8X_16BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_14BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE_2("0") \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            INIT_RESIEVE_2("16") \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r9d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            "sall   $16,%%r9d    \n\t"   \
            "orl    %%r9d, %%r8d  \n\t"   /* now has 16 comparisons (taking 32 bits) */ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_15BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE_2("0") \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            INIT_RESIEVE_2("16") \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r9d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            "sall   $16,%%r9d    \n\t"   \
            "orl    %%r9d, %%r8d  \n\t"   /* now has 16 comparisons (taking 32 bits) */ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_16BIT_MAX_VEC \
		ASM_ ( \
			INIT_RESIEVE_2("0") \
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            INIT_RESIEVE_2("16") \
			STEP_COMPARE_COMBINE	\
			"vpor	%%xmm0, %%xmm7, %%xmm7 \n\t" \
			"vpmovmskb %%xmm7, %%r9d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            "sall   $16,%%r9d    \n\t"   \
            "orl    %%r9d, %%r8d  \n\t"   /* now has 16 comparisons (taking 32 bits) */ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_14BIT_MAX_VEC_AVX2 \
		ASM_ ( \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			"vpor	%%ymm0, %%ymm7, %%ymm7 \n\t" \
			"vpmovmskb %%ymm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_15BIT_MAX_VEC_AVX2 \
		ASM_ ( \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			"vpor	%%ymm0, %%ymm7, %%ymm7 \n\t" \
			"vpmovmskb %%ymm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);

#define RESIEVE_16X_16BIT_MAX_VEC_AVX2 \
		ASM_ ( \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			"vpor	%%ymm0, %%ymm7, %%ymm7 \n\t" \
			"vpmovmskb %%ymm7, %%r8d \n\t"		/* if one of these primes divides this location, this will be !0*/ \
            GATHER_BIT_INDICES \
			: "+r"(result) \
			: "r"(fbc->prime + i), "r"(fbc->root1 + i), "r"(fbc->root2 + i), "r"(corrections), "r"(buffer), "r"(i) \
			: "r8", "r9", "r10", "r11", "rcx", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7", "xmm0", "cc", "memory" \
			);


#else

#define TDIV_MED_CLEAN

#define STEP_COMPARE_COMBINE \
        v128_x2 = _mm_sub_epi16(v128_x2, v128_p); \
        v128_x3 = _mm_sub_epi16(v128_x3, v128_p); \
        v128_x5 = _mm_cmpeq_epi16(v128_x2, v128_x5); \
        v128_x6 = _mm_cmpeq_epi16(v128_x3, v128_x6); \
        v128_x7 = _mm_or_si128(v128_x7, v128_x5); \
        v128_x0 = _mm_or_si128(v128_x0, v128_x6);


#define STEP_COMPARE_COMBINE_AVX2 \
		v256_y2 = _mm256_sub_epi16(v256_y2, v256_p); \
        v256_y3 = _mm256_sub_epi16(v256_y3, v256_p); \
        v256_y5 = _mm256_cmpeq_epi16(v256_y2, v256_y5); \
        v256_y6 = _mm256_cmpeq_epi16(v256_y3, v256_y6); \
        v256_y7 = _mm256_or_si256(v256_y7, v256_y5); \
        v256_y0 = _mm256_or_si256(v256_y0, v256_y6);

#define INIT_RESIEVE \
__m128i v128_p, v128_x0, v128_x2, v128_x3, v128_x4, v128_x5, v128_x6, v128_x7; \
uint32_t msk32, pos; \
v128_x4 = _mm_load_si128(corrections);                  \
v128_x0 = _mm_xor_si128(v128_x0, v128_x0);              \
v128_x2 = _mm_load_si128(fbc->root1 + i);               \
v128_x3 = _mm_load_si128(fbc->root2 + i);               \
v128_x2 = _mm_add_epi16(v128_x4, v128_x2);              \
v128_x3 = _mm_add_epi16(v128_x4, v128_x3);              \
v128_p = _mm_load_si128(fbc->prime + i);                \
v128_x5 = _mm_xor_si128(v128_x5, v128_x5);              \
v128_x6 = _mm_xor_si128(v128_x6, v128_x6);              \
v128_x7 = _mm_xor_si128(v128_x7, v128_x7);              


#define INIT_RESIEVE_AVX2 \
__m256i v256_p, v256_y0, v256_y2, v256_y3, v256_y4, v256_y5, v256_y6, v256_y7; \
v256_p = v256_y0 = v256_y2 = v256_y3 = v256_y4 = v256_y5 = v256_y6 = v256_y7 = _mm256_setzero_si256(); \
uint32_t msk32, pos; \
v256_y4 = _mm256_load_si256(corrections);                   \
v256_y0 = _mm256_xor_si256(v256_y0, v256_y0);               \
v256_y2 = _mm256_load_si256(fbc->root1 + i);                \
v256_y3 = _mm256_load_si256(fbc->root2 + i);                \
v256_y2 = _mm256_add_epi16(v256_y4, v256_y2);               \
v256_y3 = _mm256_add_epi16(v256_y4, v256_y3);               \
v256_p = _mm256_load_si256(fbc->prime + i);                 \
v256_y5 = _mm256_xor_si256(v256_y5, v256_y5);               \
v256_y6 = _mm256_xor_si256(v256_y6, v256_y6);               \
v256_y7 = _mm256_xor_si256(v256_y7, v256_y7);    




#define GATHER_BIT_INDICES \
        msk32 &= 0xaaaaaaaa; \
        while (_BitScanForward(&pos, msk32)) { \
            buffer[result++] = (pos >> 1) + i; \
            _reset_lsb(msk32); \
        }


#define RESIEVE_8X_14BIT_MAX \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
            v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
            /* if one of these primes divides this location, this will be !0*/ \
            result = _mm_movemask_epi8(v128_x7);

#define RESIEVE_8X_15BIT_MAX \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
            /* if one of these primes divides this location, this will be !0*/ \
            result = _mm_movemask_epi8(v128_x7);

#define RESIEVE_8X_16BIT_MAX \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
            /* if one of these primes divides this location, this will be !0*/ \
            result = _mm_movemask_epi8(v128_x7);

#define RESIEVE_8X_14BIT_MAX_VEC \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
			msk32 = _mm_movemask_epi8(v128_x7); \
            GATHER_BIT_INDICES

#define RESIEVE_8X_15BIT_MAX_VEC \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			STEP_COMPARE_COMBINE	\
			v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
			msk32 = _mm_movemask_epi8(v128_x7); \
            GATHER_BIT_INDICES

#define RESIEVE_8X_16BIT_MAX_VEC \
			INIT_RESIEVE \
			STEP_COMPARE_COMBINE	\
			v128_x7 = _mm_xor_si128(v128_x0, v128_x7); \
			msk32 = _mm_movemask_epi8(v128_x7); \
            GATHER_BIT_INDICES

#define RESIEVE_16X_14BIT_MAX_VEC_AVX2 \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			v256_y7 = _mm256_xor_si256(v256_y0, v256_y7); \
			msk32 = _mm256_movemask_epi8(v256_y7); \
            GATHER_BIT_INDICES

#define RESIEVE_16X_15BIT_MAX_VEC_AVX2 \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			STEP_COMPARE_COMBINE_AVX2	\
			v256_y7 = _mm256_xor_si256(v256_y0, v256_y7); \
			msk32 = _mm256_movemask_epi8(v256_y7); \
            GATHER_BIT_INDICES

#define RESIEVE_16X_16BIT_MAX_VEC_AVX2 \
			INIT_RESIEVE_AVX2 \
			STEP_COMPARE_COMBINE_AVX2	\
			v256_y7 = _mm256_xor_si256(v256_y0, v256_y7); \
			msk32 = _mm256_movemask_epi8(v256_y7); \
            GATHER_BIT_INDICES


#endif



#define CHECK_8_RESULTS \
	if (result & 0x2) {				  \
		DIVIDE_RESIEVED_PRIME(0);	  \
    	}								  \
	if (result & 0x8){				  \
		DIVIDE_RESIEVED_PRIME(1);	  \
    	}								  \
	if (result & 0x20){				  \
		DIVIDE_RESIEVED_PRIME(2);	  \
    	}								  \
	if (result & 0x80){				  \
		DIVIDE_RESIEVED_PRIME(3);	  \
    	}								  \
	if (result & 0x200){			  \
		DIVIDE_RESIEVED_PRIME(4);	  \
    	}								  \
	if (result & 0x800){			  \
		DIVIDE_RESIEVED_PRIME(5);	  \
    	}								  \
	if (result & 0x2000){			  \
		DIVIDE_RESIEVED_PRIME(6);	  \
    	}								  \
	if (result & 0x8000){			  \
		DIVIDE_RESIEVED_PRIME(7);	  \
    	}


//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria
indicating they may factor completely over the factor base.

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14,
depending on the version (2^13 for 32k version, 2^14 for 64k).

4) Resieve primes between 2^{13|14} and 2^16, max.

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 4)


*/

void resieve_medprimes_32k_avx2(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    int i;
    uint32_t report_num;
    int smooth_num;
    uint32_t *fb_offsets;
    sieve_fb_compressed *fbc;
    uint32_t block_loc;
    ALIGNED_MEM uint16_t corrections[32]; // = dconf->corrections;
    uint16_t buffer[16];
    uint32_t result = 0;
    uint32_t bound14;
    uint32_t bound15;
    uint32_t bound16;

    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }

    // 16x resieving
    if ((sconf->factor_base->fb_14bit_B & 15) == 0)
    {
        bound14 = sconf->factor_base->fb_14bit_B;
    }
    else
    {
        bound14 = sconf->factor_base->fb_14bit_B - 8;
    }

    // determine the next bound
    if ((sconf->factor_base->fb_15bit_B & 15) == 0)
    {
        bound15 = sconf->factor_base->fb_15bit_B;
    }
    else
    {
        bound15 = sconf->factor_base->fb_15bit_B - 8;
    }

    // determine the next bound
    if ((sconf->factor_base->med_B & 15) == 0)
    {
        bound16 = sconf->factor_base->med_B;
    }
    else
    {
        bound16 = sconf->factor_base->med_B - 8;
    }

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {

        if (!dconf->valid_Qs[report_num])
            continue;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        // where tdiv_medprimes left off
        i = sconf->factor_base->fb_13bit_B;

        // the roots have already been advanced to the next block.
        // we need to correct them back to where they were before resieving.
        INIT_CORRECTIONS;

        result = 0;

        if ((i & 15) != 0)
        {
            RESIEVE_8X_14BIT_MAX_VEC;
            i += 8;
        }

        while ((uint32_t)i < bound14)
        {

            RESIEVE_16X_14BIT_MAX_VEC_AVX2;

            i += 16;
        }

        while ((uint32_t)i < bound15)
        {

            RESIEVE_16X_15BIT_MAX_VEC_AVX2;

            i += 16;
        }

        while ((uint32_t)i < bound16)
        {

            RESIEVE_16X_16BIT_MAX_VEC_AVX2;

            i += 16;
        }

        if (i != sconf->factor_base->med_B)
        {
            RESIEVE_8X_16BIT_MAX_VEC;
            i += 8;
        }


        CLEAN_AVX2;

        for (i = 0; i < result; i++)
        {
            if (buffer[i] < 2)
                continue;

            DIVIDE_RESIEVED_PRIME_2((buffer[i]));
        }

        CLEAN_AVX2;

        // after either resieving or standard trial division, record
        // how many factors we've found so far.
        dconf->smooth_num[report_num] = smooth_num;

    }

    TDIV_MED_CLEAN;

    return;
}

#endif // USE_AVX2


#ifdef USE_AVX512BW
void resieve_medprimes_32k_avx512bw(uint8_t parity, uint32_t poly_id, uint32_t bnum,
    static_conf_t* sconf, dynamic_conf_t* dconf)
{
    //we have flagged this sieve offset as likely to produce a relation
    //nothing left to do now but check and see.
    int i;
    uint32_t report_num;
    int smooth_num;
    uint32_t* fb_offsets;
    sieve_fb_compressed* fbc;
    uint32_t block_loc;
    ALIGNED_MEM uint16_t corrections[32]; // = dconf->corrections;
    uint16_t buffer[16];
    uint32_t result = 0;
    uint32_t bound14;
    uint32_t bound15;
    uint32_t bound16;

    if (parity)
    {
        fbc = dconf->comp_sieve_n;
    }
    else
    {
        fbc = dconf->comp_sieve_p;
    }
       
    // 32x resieving
    bound14 = sconf->factor_base->fb_14bit_B;
    bound15 = sconf->factor_base->fb_15bit_B;
    bound16 = sconf->factor_base->med_B;

    __m512i vzero = _mm512_setzero_epi32();

    for (report_num = 0; report_num < dconf->num_reports; report_num++)
    {

        if (!dconf->valid_Qs[report_num])
            continue;

        // pull the details of this report to get started.
        fb_offsets = &dconf->fb_offsets[report_num][0];
        smooth_num = dconf->smooth_num[report_num];
        block_loc = dconf->reports[report_num];

        // where tdiv_medprimes left off
        i = sconf->factor_base->fb_13bit_B;

        // the roots have already been advanced to the next block.
        // we need to correct them back to where they were before resieving.
        INIT_CORRECTIONS_32;

        result = 0;

        if ((i & 31) != 0)
        {
            RESIEVE_8X_14BIT_MAX_VEC_AVX512;
            i += 8;
        }

        while ((uint32_t)(i + 32) < bound14)
        {

            RESIEVE_32X_14BIT_MAX_VEC_AVX512;

            i += 32;
        }

        while ((uint32_t)(i + 32) < bound15)
        {

            RESIEVE_32X_15BIT_MAX_VEC_AVX512;

            i += 32;
        }

        while ((uint32_t)(i + 32) < bound16)
        {
            RESIEVE_32X_16BIT_MAX_VEC_AVX512;

            i += 32;
        }

        while (i < sconf->factor_base->med_B)
        {
            RESIEVE_8X_16BIT_MAX_VEC_AVX512;
            i += 8;
        }

        for (i = 0; i < result; i++)
        {
			if ((buffer[i] < 2) || (fbc->prime[buffer[i]] < 2))
			{
				// sometimes this happens for small inputs... TODO: figure
				// out why and fix it.  But for now detecting and ignoring
				// also works.
				// printf("attempting to divide by %u at fb index %u during resieve\n",
				//	fbc->prime[buffer[i]], buffer[i]);
				break;
			}

            DIVIDE_RESIEVED_PRIME_2((buffer[i]));
        }

        // after either resieving or standard trial division, record
        // how many factors we've found so far.
        dconf->smooth_num[report_num] = smooth_num;

    }



    TDIV_MED_CLEAN;

    return;
}
#endif

