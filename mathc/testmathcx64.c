#include <stdio.h>
#include <inttypes.h>
#include <emmintrin.h>
#include <x86gprintrin.h>


uint64_t mulredc(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    if (n & 0x8000000000000000)
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "movq $0, %%r12 \n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "cmovae %4, %%r12 \n\t"
            "xorq %%rax, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %%r12, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "r12", "cc");
    }
    else
    {
        __asm__(
            "mulq %2	\n\t"
            "movq %%rax, %%r10		\n\t"
            "movq %%rdx, %%r11		\n\t"
            "mulq %3 \n\t"
            "mulq %4 \n\t"
            "addq %%r10, %%rax \n\t"
            "adcq %%r11, %%rdx \n\t"
            "movq $0, %%rax \n\t"
            "subq %4, %%rdx \n\t"
            "cmovc %4, %%rax \n\t"
            "addq %%rdx, %%rax \n\t"
            : "=&a"(x)
            : "0"(x), "r"(y), "r"(nhat), "r"(n)
            : "rdx", "r10", "r11", "cc");

    }
    return x;
}


uint64_t _umul128(uint64_t x, uint64_t y, uint64_t* hi)
{
    __asm__(
        "mulq %3	\n\t"
        : "=&a"(x), "=&d"(y)
        : "0"(x), "1"(y)
        : "cc"
    );

    *hi = y;
    return x;
}

uint64_t mulredc2(uint64_t x, uint64_t y, uint64_t n, uint64_t nhat)
{
    uint64_t th, tl, u, ah, al;
    tl = _umul128(x, y, &th);
    u = tl * nhat;
    al = _umul128(u, n, &ah);
    tl = _addcarry_u64(0, al, tl, &al);
    th = _addcarry_u64((uint8_t)tl, th, ah, &x);
    if (th || (x >= n)) x -= n;
    return x;
}



uint64_t u64div(uint64_t c, uint64_t n)
{
#if 1
    __asm__("divq %4"
        : "=a"(c), "=d"(n)
        : "1"(c), "0"(0ULL), "r"(n));
#else
// this should work if the above won't compile (e.g. on clang)
    uint64_t tmp = 0;
    __asm__("divq %4"
        : "=a"(tmp), "=d"(n)
        : "1"(c), "0"(tmp), "r"(n));
#endif
    return n;
}





uint64_t spDivide(uint64_t* q, uint64_t* r, uint64_t u[2], uint64_t v)
{
    *r = u[1];
    *q = u[0];
    __asm__("divq %4"
        : "=a"(*q), "=d"(*r)
        : "1"(*r), "0"(*q), "r"(v));

    return 0;
}


void spMultiply(uint64_t u, uint64_t v, uint64_t* product, uint64_t* carry)
{
    *product = v;
    *carry = u;

    __asm__("movq %2, %%rax	\n\t"
        "mulq %3	\n\t"
        "movq %%rax, %0		\n\t"
        "movq %%rdx, %1		\n\t"
        : "=r"(*product), "=r"(*carry)
        : "1"(*carry), "0"(*product)
        : "rax", "rdx", "cc");

    return;
}



void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t* w)
{
    uint64_t p[2];
    uint64_t q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
}

uint8_t _addcarry_u642(uint64_t x, uint8_t w, uint64_t y, uint64_t *sum)
{
    uint64_t s, c; 
    s = y;
    c = 0;

    __asm__("movq %2, %%rax		\n\t"
        "addq %3, %%rax		\n\t"
        "adcq $0, %5		\n\t"
        "addq %%rax, %4		\n\t"
        "adcq $0, %5		\n\t"
        : "=r"(s), "=r"(c)
        : "r"(x), "r"((uint64_t)w), "0"(s), "1"(c)
        : "rax", "memory", "cc");

    *sum = s;
    return c;
}



int main(){
	__m128i c1=_mm_set_epi32(0xBA987654,0x3210FEDC,0xFEDCBA98,0x76543210);
	
	__m128i c2=_mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	
	__m128i c3=_mm_set_epi32(0x76543210,0xFEDCBA98,0x3210FEDC,0xBA987654);
	
		uint64_t* o2 = (uint64_t*)(&c2);
		uint64_t* o1 = (uint64_t*)(&c1);
		uint64_t* o3 = (uint64_t*)(&c3);	
		
		printf("orgin-----------------%"PRIx64"\n",*o2);
		printf("orgin-----------------%"PRIx64"\n",*(o2+1));
		printf("orgin-----------------%"PRIx64"\n",*o1);
		printf("orgin-----------------%"PRIx64"\n",*(o1+1));
		printf("orgin-----------------%"PRIx64"\n",*o3);
		printf("orgin-----------------%"PRIx64"\n",*(o3+1));
	
#ifdef test_mm_storeu_si128  // test ok 
		 
		_mm_storeu_si128(&c2 ,c1);
		
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		
#endif 

#ifdef test_mm_loadu_si128 // test ok 
		c2 = _mm_loadu_si128(&c1);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 
 

#ifdef test_mm_load_si128 // test ok 
		c2 = _mm_load_si128(&c1);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 

#ifdef test_mm_extract_epi16 // test ok
		int result = _mm_extract_epi16(c1 , 0x5);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		
		printf("-----------------%"PRIx32"\n",result); // test ok
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 

#ifdef test_mm_add_epi16 // test ok
				
		c2=	_mm_add_epi16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));

#endif 


#ifdef test_mm_subs_epu16 // test ok
		c2=	_mm_subs_epu16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 


#ifdef test_mm_movemask_epi8 // test ok 
		int result = _mm_movemask_epi8(c1);
		uint64_t* m1 = (uint64_t*)(&c1);
		printf("-----------------%"PRIx32"\n",result); // test ok
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 


#ifdef test_umul128   // test ok 
	uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
	uint64_t* m1 = (uint64_t*)(&c1);
	uint64_t result = _umul128(*m1,*(m1+1),t1);
	printf("-----------------%"PRIx64"\n",result); 
	printf("-----------------%"PRIx64"\n",*t1); 
			
#endif 


#ifdef test_addcarry_u64 // test not ok
	uint64_t* m1 = (uint64_t*)(&c1);
	uint64_t* m3 = (uint64_t*)(&c3);
	uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
	
	uint8_t result =_addcarry_u64(*(m1+1),*m1,*m3,t1);
	printf("-----------------%"PRIx8"\n",result);
	printf("-----------------%"PRIx64"\n",*t1);
	
	result =_addcarry_u642(*m1,*(m1+1),*m3,t1);
	printf("-----------------%"PRIx8"\n",result);
	printf("-----------------%"PRIx64"\n",*t1);
	
	
#endif 

#ifdef test_mm_and_si128 // test ok 

		c2=	_mm_and_si128(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_sub_epi16 // test ok 
		c2=	_mm_sub_epi16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 
 
#ifdef test_subborrow_u64  // test not ok 
			uint64_t* m1 = (uint64_t*)(&c1);
			uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
			uint8_t t2 = 0x15 ;
			uint8_t result = 	_subborrow_u64(t2 ,*m1 , *(m1+1) ,t1 );
			printf("-----------------%"PRIx8"\n",result); 
			printf("-----------------%"PRIx64"\n",*m1);
			printf("-----------------%"PRIx64"\n",*(m1+1));
			printf("-----------------%"PRIx64"\n",*t1);
			printf("-----------------%"PRIx8"\n",t2);
			
#endif 

#ifdef test_subborrow_u32 // test ok 
			uint32_t* m1 = (uint32_t*)(&c1);
			uint32_t* t1 = (uint32_t*)malloc(sizeof(uint32_t));
			uint8_t t2 = 0x15 ;
			uint8_t result =_subborrow_u32(t2 ,*m1 , *(m1+1) ,t1);
			printf("-----------------%"PRIx8"\n",result); 
			printf("-----------------%"PRIx32"\n",*m1);
			printf("-----------------%"PRIx32"\n",*(m1+1));
			printf("-----------------%"PRIx32"\n",*t1);
			printf("-----------------%"PRIx8"\n",t2);
#endif 

#ifdef test_mm_sub_epi32  // test ok
		c2=	_mm_sub_epi32(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_cmpgt_epi32  // test ok

		c2=	_mm_cmpgt_epi32(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 
 
#ifdef test_mm_xor_si128  // test ok

		c2=	_mm_xor_si128(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_add_epi32  // test ok 

			c2=	_mm_add_epi32(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_cmpeq_epi16  // test ok

		c2=	_mm_cmpeq_epi16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_max_epi16 // test ok 

		c2=	_mm_max_epi16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_min_epi16 // test ok 
	
		c2=	_mm_min_epi16(c1,c3);
		uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		uint64_t* m3 = (uint64_t*)(&c3);	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif 

#ifdef test_mm_store_si128 // test ok 
				_mm_min_epi16(&c2,c1); 
				uint64_t* m2 = (uint64_t*)(&c2);
		uint64_t* m1 = (uint64_t*)(&c1);
		 	
		printf("-----------------%"PRIx64"\n",*m2);
		printf("-----------------%"PRIx64"\n",*(m2+1));
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 

#ifdef test_udiv128  // test ok 
	
	uint64_t* m2 = (uint64_t*)(&c2);
	uint64_t* m1 = (uint64_t*)(&c1);
	uint64_t* m3 = (uint64_t*)(&c3);
	uint64_t result = u64div(*m3,*(m1));
	printf("-----------------%"PRIx64"\n",result); 
	result = spDivide(m2,(m2+1) ,m3 , *m1);	
	printf("-----------------%"PRIx64"\n",result); 
	printf("-----------------%"PRIx64"\n",*m2); 
	printf("-----------------%"PRIx64"\n",*(m2+1)); 
#endif 

#ifdef test_mm256_zeroupper

			printf("-----------------%"PRIx32"\n",result); // test ok
#endif 

#ifdef test_mm_shuffle_epi32  // test ok 
		__m128i result =_mm_shuffle_epi32(c1,0xe4);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif

#ifdef test_mm_srli_epi16  //  
		__m128i result =_mm_srli_epi16(c1,0x4);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif

#ifdef test_mm_mulhi_epu16  //  
		__m128i result =_mm_mulhi_epu16(c1,c3);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif

#ifdef test_mm_or_si128  //  
		__m128i result =_mm_or_si128(c1,c3);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif


#ifdef test_mm_mullo_epi16  //  
		__m128i result =_mm_mullo_epi16(c1,c3);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif


#ifdef test_spMulMod
	uint64_t a = 46380 ;
	uint64_t b = 46381 ;
	spMulMod(a, a, b, &a);
	printf("-----------------%"PRIx64"\n",a);
	
#endif 
#ifdef test_cmpara
	char c = 0xff ; 
	if(c>0){
		printf("%d > 0 \n",c);
	}else{
		printf("%d <0 \n" ,c);
	}
#endif
	
	return 0 ;
}
