#include <stdio.h>
#include <inttypes.h>
#include "mathc.h"
#include <stdlib.h>

uint64_t spDivide(uint64_t * q, uint64_t * r, uint64_t u[2], uint64_t v)
{
    *q = _udiv128(u[1], u[0], v, r);
    return 0;
}


void spMultiply(uint64_t u, uint64_t v, uint64_t* product, uint64_t* carry)
{
    *product = _umul128(u, v, carry);
    return;
}


uint64_t u64div(uint64_t c, uint64_t n)
{
    uint64_t r;
    _udiv128(c, 0, n, &r);
    return r;
}



void spMulMod(uint64_t u, uint64_t v, uint64_t m, uint64_t* w)
{
    uint64_t p[2];
    uint64_t q;

    spMultiply(u, v, &p[0], &p[1]);
    spDivide(&q, w, p, m);

    return;
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
	
#ifdef test_mm_storeu_si128 // test ok 
		 
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
 

#ifdef test_mm_load_si128	// test ok 
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
		
		printf("-----------------%"PRIx32"\n",result);  
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
		printf("-----------------%"PRIx32"\n",result);
		printf("-----------------%"PRIx64"\n",*m1);
		printf("-----------------%"PRIx64"\n",*(m1+1));
#endif 


#ifdef test_umul128 // test ok 
	uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
	uint64_t* m1 = (uint64_t*)(&c1);
	uint64_t result = _umul128(*m1,*(m1+1),t1);
	printf("-----------------%"PRIx64"\n",result); 
	printf("-----------------%"PRIx64"\n",*t1); 
#endif 


#ifdef test_addcarry_u64 // test ok
	uint64_t* m1 = (uint64_t*)(&c1);
	uint64_t* m3 = (uint64_t*)(&c3);
	uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
	
		uint8_t result =_addcarry_u64_window(*(m1+1),*m1,*m3,t1);
	printf("-----------------%"PRIx8"\n",result);
	printf("-----------------%"PRIx64"\n",*t1);
	result =_addcarry_u64(*m1,*(m1+1),*m3,t1);
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

#ifdef test_subborrow_u64 // wait to test

			uint64_t* m1 = (uint64_t*)(&c1);
			uint64_t* t1 = (uint64_t*)malloc(sizeof(uint64_t));
			uint8_t t2 = 0x15 ;
			uint8_t result = _subborrow_u64(t2 ,*m1 , *(m1+1) ,t1);
			printf("-----------------%"PRIx8"\n",result); 
			printf("-----------------%"PRIx64"\n",*m1);
			printf("-----------------%"PRIx64"\n",*(m1+1));
			printf("-----------------%"PRIx64"\n",*t1);
			printf("-----------------%"PRIx8"\n",t2);
#endif 

#ifdef test_subborrow_u32 // wait to test

			uint32_t* m1 = (uint32_t*)(&c1);
			uint32_t* t1 = (uint32_t*)malloc(sizeof(uint32_t));
			uint8_t t2 = 0x15 ;
			uint8_t result =_subborrow_u32(t2 ,*m1 , *(m1+1) ,t1  );
			printf("-----------------%"PRIx8"\n",result); 
			printf("-----------------%"PRIx32"\n",*m1);
			printf("-----------------%"PRIx32"\n",*(m1+1));
			printf("-----------------%"PRIx32"\n",*t1);
			printf("-----------------%"PRIx8"\n",t2);
#endif 

#ifdef test_mm_sub_epi32								 // test ok
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

#ifdef test_mm_cmpgt_epi32 		// test ok

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

#ifdef test_mm_xor_si128 		// test ok

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

#ifdef test_mm_add_epi32 		// test ok 

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

#ifdef test_mm_max_epi16 	// test ok 

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

#ifdef test_mm_min_epi16  // test ok 
	
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

#ifdef test_mm_store_si128		// test ok 
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

#ifdef test_mm_shuffle_epi32
	__m128i result =_mm_shuffle_epi32(c1,0xe4);
	uint64_t* m3 = (uint64_t*)(&result);
	printf("-----------------%"PRIx64"\n",*m3);
	printf("-----------------%"PRIx64"\n",*(m3+1));	
#endif



#ifdef test_mm_srli_epi16  //   test ok
		__m128i result =_mm_srli_epi16(c1,0x4);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif

#ifdef test_mm_mulhi_epu16  //  test ok 
		__m128i result =_mm_mulhi_epu16(c1,c3);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif

#ifdef test_mm_or_si128  //  test ok
		__m128i result =_mm_or_si128(c1,c3);
		uint64_t* m3 = (uint64_t*)(&result);
		printf("-----------------%"PRIx64"\n",*m3);
		printf("-----------------%"PRIx64"\n",*(m3+1));
#endif


#ifdef test_mm_mullo_epi16  //  test ok 
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
