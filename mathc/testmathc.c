#include "mathc.h"
#include <stdlib.h>

int main(){
	
	uint64_t r= 0;
	uint64_t p=_udiv128(0, 566582809 ,46573, &r);
	printf("-----------------%"PRIx64"\n",p); //test ok
	printf("-----------------%"PRIx64"\n",r); //test ok
	return 0 ;
}