#ifndef __MATHC__
#define __MATHC__
#include <inttypes.h>
#include <stdio.h>
typedef long long __m128i_u __attribute__ ((__vector_size__ (16), __may_alias__, __aligned__ (1)));
typedef long long __m128i __attribute ((__vector_size__ (16) , __may_alias__));
typedef unsigned long long __v2du __attribute__((__vector_size__ (16)));
typedef unsigned int __v4su __attribute__((__vector_size__ (16)));
typedef int __v4si __attribute__((__vector_size__ (16)));
typedef unsigned short __v8hu __attribute__ ((__vector_size__ (16)));
typedef short __v8hi __attribute__ ((__vector_size__ (16)));
typedef unsigned short __v8hu __attribute__ ((__vector_size__ (16)));
void _mm_storeu_si128(__m128i *p , __m128i a);
__m128i _mm_loadu_si128(__m128i const* mem_address);
__m128i _mm_set_epi32(int i3 , int i2 , int i1 , int i0);
__m128i _mm_load_si128(__m128i *p);
int _mm_extract_epi16(__m128i a , int imm);
__m128i _mm_add_epi16(__m128i a , __m128i b);
__m128i _mm_subs_epu16(__m128i a , __m128i b );
int _mm_movemask_epi8(__m128i a );
uint64_t _umul128(uint64_t multiplier , uint64_t multiplicand , uint64_t *highproduct);
uint8_t _addcarry_u64(uint64_t x ,uint8_t w , int64_t y , int64_t *sum);
uint8_t _addcarry_u64_window(uint8_t w,uint64_t x,int64_t y,int64_t *d);
__m128i _mm_and_si128(__m128i a ,__m128i b );
__m128i _mm_sub_epi16(__m128i a ,__m128i b );
__m128i _mm_mullo_epi16(__m128i a,__m128i b);
__m128i _mm_or_si128(__m128i a,__m128i b);
__m128i _mm_mulhi_epu16(__m128i a,__m128i b);
__m128i _mm_srli_epi16(__m128i a,int b);
uint8_t _subborrow_u64(uint8_t a , uint64_t b , uint64_t c , uint64_t *d);
uint8_t _subborrow_u32(uint8_t a,uint32_t b,uint32_t c,uint32_t *d);
__m128i _mm_sub_epi32(__m128i a,__m128i b );
__m128i _mm_cmpgt_epi32 (__m128i __A, __m128i __B);
__m128i _mm_cmpgt_epi16 (__m128i __A, __m128i __B);
__m128i _mm_xor_si128 (__m128i __A, __m128i __B);
__m128i _mm_add_epi32 (__m128i __A, __m128i __B);
__m128i _mm_cmpeq_epi16(__m128i a ,__m128i b);
__m128i _mm_max_epi16(__m128i a ,__m128i b);
__m128i _mm_min_epi16(__m128i a,__m128i b);
void _mm_store_si128(__m128i *p , __m128i a );
uint64_t _udiv128(uint64_t a, uint64_t b , uint64_t c, uint64_t* r);
void _mm256_zeroupper(void);
__m128i _mm_shuffle_epi32 (__m128i __A, const int __mask);
#endif 
