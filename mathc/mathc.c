#include "mathc.h"





void _mm_storeu_si128(__m128i_u *p,__m128i a){
        *p =a;	
}

__m128i _mm_loadu_si128(__m128i const* mem_address){
    return *mem_address;
}

__m128i _mm_set_epi32(int i3,int i2,int i1,int i0){
    __m128i m ;
    int* x= (int*)(&m);
    x[0] = i0 ; 
    x[1] = i1 ;
    x[2] =i2 ;
    x[3] = i3 ;
    return m;
}

__m128i _mm_load_si128(__m128i *p){
    return *p;
}


int _mm_extract_epi16(__m128i a,int imm){
    uint16_t* x = (uint16_t*)(&a);
    return (int)x[imm];
}

__m128i _mm_add_epi16(__m128i __A,__m128i __B){
   /*__m128i r  ;
    int16_t* i = (int16_t*)(&r);
    int16_t* j = (int16_t*)(&a) ;
    int16_t* z = (int16_t*)(&b);
    i[0] = j[0] + z[0] ;
    i[1] = j[1] + z[1] ;
    i[2] = j[2] + z[2] ;
    i[3] = j[3] + z[3] ;
    i[4] = j[4] + z[4] ;
    i[5] = j[5] + z[5] ;
    i[6] = j[6] + z[6] ;
    i[7] = j[7] + z[7] ;
    return r;*/
    return (__m128i) ((__v8hu)__A + (__v8hu)__B);
}



// has error 
__m128i _mm_subs_epu16(__m128i a,__m128i b ){
     __m128i r ;
    uint16_t* i = (uint16_t*)(&r);
    uint16_t* j = (uint16_t*)(&a) ;
    uint16_t* z = (uint16_t*)(&b);
    i[0] =(j[0] > z[0])?(j[0] - z[0]):0;
    i[1] =(j[1] > z[1])?(j[1] - z[1]):0;
    i[2] =(j[2] > z[2])?(j[2] - z[2]):0 ;
    i[3] =(j[3] > z[3])?(j[3] - z[3]):0 ;
    i[4] =(j[4] > z[4])?(j[4] - z[4]):0 ;
    i[5] =(j[5] > z[5])?(j[5] - z[5]):0 ;
    i[6] =(j[6] > z[6])?(j[6] - z[6]):0 ;
    i[7] =(j[7] > z[7])?(j[7] - z[7]):0 ;
    return r;
}

int _mm_movemask_epi8(__m128i a ){
    uint8_t* pa = (uint8_t*)(&a);
    
    return ((pa[15]>>7) <<15) | (pa[14]>>7) <<14| (pa[13]>>7)<<13| (pa[12]>>7)<<12| (pa[11]>>7)<<11| (pa[10]>>7)<<10|(pa[9]>>7)<<9 |(pa[8]>>7)<<8 |(pa[7]>>7) <<7|(pa[6]>>7) <<6|(pa[5]>>7)<<5 |(pa[4]>>7) <<4|(pa[3]>>7)<<3 |(pa[2]>>7)<<2 |(pa[1]>>7) <<1| pa[0] >>7;
}

uint64_t _umul128(uint64_t a,uint64_t b,uint64_t *highproduct){

    uint64_t high = 0 ;
    uint64_t low = 0;
    __asm__(
    "sub sp,sp,16 \n\t"
    "mov x1,%2 \n\t"
    "mov x2,%3 \n\t"
    "umulh x3,x1,x2 \n\t"
    "umull x4,w1,w2 \n\t"
    "str x3,[sp] \n\t"
    "str x4,[sp,8] \n\t"
    "lsr x3,x1,#32 \n\t"
    "lsr x4,x2,#32 \n\t"
    "umull x4,w4,w1 \n\t"
    "umull x3,w3,w2 \n\t "
    "LSL x4,x4,32 \n\t"
    "LSL x3,x3,32 \n\t"
    "add x1,x4,x3 \n\t"
    "ldr x2,[sp,8] \n\t"
    "add %0,x1,x2 \n\t"
    "ldr %1,[sp] \n\t "
    "add sp,sp,16 \n\t"
    :"=r"(low),"=r"(high)
    :"r"(a),"r"(b)
    : "x1","x2","x3","x4","cc");

    *highproduct = high ;
    return low;
}

uint64_t _udiv128(uint64_t high,uint64_t low,uint64_t divisor,uint64_t* r){
    
    uint64_t b1 = 0x0 ;
	uint64_t b2 = divisor ;
	uint64_t a0 = 64 ;
	uint64_t result = 0 ;
	while(a0>0){
		result <<=1;
		a0-=1;
		
		b1 = divisor<<a0;
		if(a0==0){
			b2=0;
		}else{
			b2 = divisor>>(64-a0);
		}
		
		if(high > b2){
			result+=1 ;
		}else if(high < b2){
			continue;
		}else if(low < b1){
			continue;
		}else{
			result+=1;
		}
		
		if(low < b1){
			high-=b2;
			high-=1;
		}else{
			high-=b2;
		}
		low-=b1;
	}
	*r = low;
	return result;
}

uint8_t _subborrow_u64(uint8_t a,uint64_t b,uint64_t c,uint64_t *d){
    //uint64_t s ;
   // uint8_t dc =0 ;
    uint8_t bv  = a!=0?1:0;
   /* __asm__(
		"subs %2,%2,%3 \n\t"
        "subs %2,%2,%4 \n\t"
        "adc %1,%1,xzr \n\t"
        "mov %0,%2  \n\t"
        :"=&r"(s),"+&r"(dc)
        :"r"(b),"r"(c),"r"(bv)
        :"cc");
    *d =s ;
    return dc;*/
	
	uint64_t Diff = b - c - bv;
	uint64_t CarryVector = (Diff & c) ^ ((Diff ^ c) & ~b);
	*d = Diff;
 return CarryVector >> 63;
}

__m128i _mm_and_si128(__m128i a,__m128i b ){
    return (__m128i)((__v2du)a&(__v2du)b);
}


__m128i _mm_add_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4su)__A + (__v4su)__B);
}

__m128i _mm_xor_si128 (__m128i __A, __m128i __B){
   return (__m128i) ((__v2du)__A ^ (__v2du)__B);
}

 __m128i _mm_cmpgt_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4si)__A > (__v4si)__B);
}

 __m128i _mm_cmpgt_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hi)__A > (__v8hi)__B);
}

__m128i _mm_sub_epi16(__m128i __A,__m128i __B ){
    /*__m128i r  ;
    int16_t* i = (int16_t*)(&r);
    int16_t* j = (int16_t*)(&a) ;
    int16_t* z = (int16_t*)(&b);
    i[0] = j[0] - z[0] ;
    i[1] = j[1] - z[1] ;
    i[2] = j[2] - z[2] ;
    i[3] = j[3] - z[3] ;
    i[4] = j[4] - z[4] ;
    i[5] = j[5] - z[5] ;
    i[6] = j[6] - z[6] ;
    i[7] = j[7] - z[7] ;
    return r;*/
    return (__m128i) ((__v8hu)__A - (__v8hu)__B);
}

// wait to test
__m128i _mm_sub_epi32(__m128i __A,__m128i __B ){
   /* __m128i r  ;
    int32_t* i = (int32_t*)(&r);
    int32_t* j = (int32_t*)(&a) ;
    int32_t* z = (int32_t*)(&b);
    i[0] = j[0] - z[0] ;
    i[1] = j[1] - z[1] ;
    i[2] = j[2] - z[2] ;
    i[3] = j[3] - z[3] ;
    return r;*/
    return (__m128i) ((__v4su)__A - (__v4su)__B);
}

uint8_t _subborrow_u32(uint8_t a,uint32_t b,uint32_t c,uint32_t *d){
	 uint8_t bv  = a!=0?1:0;
    /*uint32_t s ;
    uint8_t dc =0 ;
   
    __asm__(
        "subs %w2,%w2,%w3 \n\t"
        "subs %w2,%w2,%w4 \n\t"
        "adc %w1,%w1,wzr \n\t"
        "mov %w0,w1  \n\t"
        :"=&r"(s),"=&r"(dc)
        :"r"(b),"r"(c),"r"(bv)
        :"x1","cc");
    *d =s ;
    return dc;*/
	
  uint32_t Diff = b - c - bv;
  uint32_t CarryVector = (Diff & c) ^ ((Diff ^ c) & ~b);
   *d = Diff;
   return CarryVector >> 31;
}



uint8_t _addcarry_u64(uint64_t x,uint8_t w,int64_t y,int64_t *sum){
    uint64_t s =y;
    uint8_t c =0;

    __asm__(" mov x2,%2  \n\t"
        "adds x2,x2,%3 \n\t"
        "adc %1,%1,xzr \n\t"
        "adds %0,%0,x2 \n\t" 
        "adc %1,%1,xzr \n\t" 
    :"+&r"(s),"+&r"(c)
    :"r"(x),"r"((uint64_t)w)
    :"x2","memory","cc");
    *sum =s ;
    return c;
}

uint8_t _addcarry_u64_window(uint8_t w,uint64_t x,int64_t y,int64_t *d){
    uint8_t bv  = w!=0?1:0;
	uint64_t Sum = bv + x + y;
	uint64_t CarryVector = (x & y) ^ ((x ^ y) & ~Sum);
	*d = Sum;
	return CarryVector >> 63;
}


__m128i _mm_cmpeq_epi16(__m128i __A,__m128i __B){
    /*
    __m128i m ;
    int16_t* im= (int16_t*)(&m);
    int16_t* ia =(int16_t*)(&a);
    int16_t* ib =(int16_t*)(&b);
    im[0] = (ia[0] ==ib[0])?0xffff:0;
    im[1] =(ia[1] ==ib[1])?0xffff:0;
    im[2] =(ia[2] ==ib[2])?0xffff:0;
    im[3] =(ia[3] ==ib[3])?0xffff:0;
    im[4] =(ia[4] ==ib[4])?0xffff:0;
    im[5] =(ia[5] ==ib[5])?0xffff:0;
    im[6] = (ia[6] ==ib[6])?0xffff:0;
    im[7] = (ia[7] ==ib[7])?0xffff:0;
    return m;*/
    return (__m128i) ((__v8hi)__A == (__v8hi)__B);
}


__m128i _mm_max_epi16(__m128i a,__m128i b){
       __m128i m ;
    int16_t* im= (int16_t*)(&m);
    int16_t* ia =(int16_t*)(&a);
    int16_t* ib =(int16_t*)(&b);
    im[0] =(ia[0]>ib[0])?ia[0]:ib[0];
    im[1] =(ia[1]>ib[1])?ia[1]:ib[1];
    im[2] =(ia[2]>ib[2])?ia[2]:ib[2];
    im[3] =(ia[3]>ib[3])?ia[3]:ib[3];
    im[4] =(ia[4]>ib[4])?ia[4]:ib[4];
    im[5] =(ia[5]>ib[5])?ia[5]:ib[5];
    im[6] =(ia[6]>ib[6])?ia[6]:ib[6];
    im[7] =(ia[7]>ib[7])?ia[7]:ib[7];
    return m;
}
__m128i _mm_min_epi16(__m128i a,__m128i b){
      __m128i m ;
    int16_t* im= (int16_t*)(&m);
    int16_t* ia =(int16_t*)(&a);
    int16_t* ib =(int16_t*)(&b);
    im[0] =(ia[0]<ib[0])?ia[0]:ib[0];
    im[1] =(ia[1]<ib[1])?ia[1]:ib[1];
    im[2] =(ia[2]<ib[2])?ia[2]:ib[2];
    im[3] =(ia[3]<ib[3])?ia[3]:ib[3];
    im[4] =(ia[4]<ib[4])?ia[4]:ib[4];
    im[5] =(ia[5]<ib[5])?ia[5]:ib[5];
    im[6] =(ia[6]<ib[6])?ia[6]:ib[6];
    im[7] =(ia[7]<ib[7])?ia[7]:ib[7];
    return m;
}
void _mm_store_si128(__m128i *p,__m128i a ){
    *p =a ;
}

void _mm256_zeroupper(void){
	
}

__m128i _mm_shuffle_epi32 (__m128i __A, const int __mask){
	__m128i m ;
	uint32_t* m32 = (uint32_t*)&m ;
	uint32_t* t = (uint32_t*)(&__A) ;
	m32[0] = t[(__mask&0x3)] ;
	m32[1] = t[(__mask&0xC)>>2];
	m32[2] = t[(__mask&0x30)>>4];
	m32[3] = t[(__mask&0xC0)>>6];
	return  m;
}

__m128i _mm_mullo_epi16(__m128i a,__m128i b){
	return (__m128i)((__v8hu)a *(__v8hu)b);
}

__m128i _mm_or_si128(__m128i a,__m128i b){
	return (__m128i)((__v2du)a | (__v2du)b);
}

__m128i _mm_mulhi_epu16(__m128i a,__m128i b){
	__m128i m ; 
	__asm__(
	"ld1 {v0.2d} , [%1] \n  " 
	"ld1 {v1.2d} , [%2] \n  "
	"mov v2.D[0],v0.D[1] \n "
	"mov v3.D[0],v1.D[1] \n "
	"umull v4.4S,v0.4H,v1.4H \n "
	"umull v5.4S,v2.4H,v3.4H \n "
	" mov v0.H[0],v4.H[1] \n"
	" mov v0.H[1],v4.H[3] \n"
	" mov v0.H[2],v4.H[5] \n"
	" mov v0.H[3],v4.H[7] \n"
	" mov v0.H[4],v5.H[1] \n"
	" mov v0.H[5],v5.H[3] \n"
	" mov v0.H[6],v5.H[5] \n"
	" mov v0.H[7],v5.H[7] \n"
	" st1 {v0.2d} , [%0] \n"
	:
	:"r"(&m), "r"(&a), "r"(&b)
	:"v0" ,"v1" ,"v2","v3" ,"v4" ,"v5");
	return m;
}


__m128i _mm_srli_epi16(__m128i a,int b){
	 uint16_t* t = (uint16_t*)(&a);
	 __m128i m ;
	uint16_t* m32 = (uint16_t*)&m ;
	m32[0] = t[0] >> b ;
	m32[1] = t[1] >> b ;
	m32[2] = t[2] >> b ;
	m32[3] = t[3] >> b ;
	m32[4] = t[4] >> b ;
	m32[5] = t[5] >> b ;
	m32[6] = t[6] >> b ;
	m32[7] = t[7] >> b ;
	return m ;
}
