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

#include <stdint.h>
#include "qs_impl.h"
#include "ytools.h"
#include "common.h"
#include "poly_macros_32k.h"
#include "poly_macros_common.h"

#define COMPUTE_NEXT_ROOTS_BATCH(i) \
    if (gray[numB + i] > 0) { \
        root1 = (int)root1 - rootupdates[(nu[numB + i] - 1) * bound + j]; \
        root2 = (int)root2 - rootupdates[(nu[numB + i] - 1) * bound + j]; \
        root1 += ((root1 >> 31) * prime); \
        root2 += ((root2 >> 31) * prime); \
            } else { \
        root1 = (int)root1 + rootupdates[(nu[numB + i] - 1) * bound + j]; \
        root2 = (int)root2 + rootupdates[(nu[numB + i] - 1) * bound + j]; \
        root1 -= ((root1 >= prime) * prime); \
        root2 -= ((root2 >= prime) * prime); \
                }

#define COMPUTE_NEXT_ROOTS_BATCH_P(i) \
    root1 = (int)root1 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
    root2 = (int)root2 - rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
    root1 += ((root1 >> 31) * prime); \
    root2 += ((root2 >> 31) * prime);

#define COMPUTE_NEXT_ROOTS_BATCH_N(i) \
    root1 = (int)root1 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
    root2 = (int)root2 + rootupdates[(nu[numB + i] - 1) * bound + j + k]; \
    root1 -= ((root1 >= prime) * prime); \
    root2 -= ((root2 >= prime) * prime); \

//this is in the poly library, even though the bulk of the time is spent
//bucketizing large primes, because it's where the roots of a poly are updated
void nextRoots_32k(static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//update the roots 
	sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
	sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
	int *rootupdates = dconf->rootupdates;

	update_t update_data = dconf->update_data;

	uint32_t startprime = 2;
	uint32_t bound = sconf->factor_base->B;

	char v = dconf->curr_poly->nu[dconf->numB];
	char sign = dconf->curr_poly->gray[dconf->numB];
	int *ptr;
	uint16_t *sm_ptr;

	lp_bucket *lp_bucket_p = dconf->buckets;
	uint32_t med_B = sconf->factor_base->med_B;
	uint32_t large_B = sconf->factor_base->large_B;

	uint32_t j, interval; //, fb_offset;
	int k = 0,numblocks;
    uint32_t root1, root2, prime;

	int bound_index=0;
	int check_bound = BUCKET_ALLOC/2 - 1;
	uint32_t bound_val = med_B;
	uint32_t *numptr_p, *numptr_n, *sliceptr_p,*sliceptr_n;
    uint8_t* slicelogp_ptr = NULL;
    uint32_t* slicebound_ptr = NULL;
	
#if 1 //!defined(USE_POLY_SSE2_ASM) || defined(PROFILING)
	uint32_t *bptr;
	int bnum, room;
#endif

	uint8_t logp=0;
	polysieve_t helperstruct;

	numblocks = sconf->num_blocks;
	interval = numblocks << 15;
	
	if (lp_bucket_p->alloc_slices != 0) // != NULL)
	{
		lp_bucket_p->fb_bounds[0] = med_B;

		sliceptr_p = lp_bucket_p->list;
		sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

		numptr_p = lp_bucket_p->num;
		numptr_n = lp_bucket_p->num + numblocks;
		
        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
            numptr_p[j] = 0;
	
		lp_bucket_p->num_slices = 0;

        slicelogp_ptr = lp_bucket_p->logp;
        slicebound_ptr = lp_bucket_p->fb_bounds;

	}
	else
	{
		sliceptr_p = NULL;
		sliceptr_n = NULL;
		numptr_p = NULL;
		numptr_n = NULL;
	}

	ptr = &rootupdates[(v-1) * bound + startprime];	

	if (sign > 0)
	{
		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16_t)root1;
			fb_p->root2[j] = (uint16_t)root2;
			fb_n->root1[j] = (uint16_t)(prime - root2);
			fb_n->root2[j] = (uint16_t)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).		
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint32_t)update_data.sm_firstroots1[j];
			root2 = (uint32_t)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}
		
#if defined(D_HAS_SSE2) && (defined(GCC_ASM64X) || defined(_MSC_VER)) //NOTDEF //GCC_ASM64X
		
		// update 8 at a time using SSE2 and no branching
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = sconf->factor_base->fb_10bit_B;		// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_PROOTS;
			
			j = h.stop;
		}	

#else


		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}

#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16_t)update_data.sm_firstroots1[j];
			root2 = (uint16_t)update_data.sm_firstroots2[j];

			if ((prime > 32768) && ((j&7) == 0))
				break;

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}	

		// continue one at a time once we exceed 15 bits, because the 8x SSE2
		// code has a hard time with unsigned 16 bit comparisons
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		

		
#if 0 //defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		
#else

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
            uint32_t nroot1, nroot2;

			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_P;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

			FILL_ONE_PRIME_LOOP_P(j);
			FILL_ONE_PRIME_LOOP_N(j);
		}

#endif

		
#if 0 //defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
	 

#else

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_P;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

	}
	else
	{
		for (j=startprime;j<sconf->sieve_small_fb_start;j++,ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			//we don't sieve these, so ordering doesn't matter
			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;

			fb_p->root1[j] = (uint16_t)root1;
			fb_p->root2[j] = (uint16_t)root2;
			fb_n->root1[j] = (uint16_t)(prime - root2);
			fb_n->root2[j] = (uint16_t)(prime - root1);
			if (fb_n->root1[j] == prime)
				fb_n->root1[j] = 0;
			if (fb_n->root2[j] == prime)
				fb_n->root2[j] = 0;

		}

		// do one at a time up to the 10bit boundary, where
		// we can start doing things 8 at a time and be
		// sure we can use aligned moves (static_data_init).	
		for (j=sconf->sieve_small_fb_start; 
			j < sconf->factor_base->fb_10bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint32_t)update_data.sm_firstroots1[j];
			root2 = (uint32_t)update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}
		
#if defined(D_HAS_SSE2) && (defined(GCC_ASM64X) || defined(_MSC_VER)) //NOTDEF //GCC_ASM64X
		// update 8 at a time using SSE2 and no branching		
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B];
		{
			small_update_t h;
			
			h.first_r1 = update_data.sm_firstroots1;		// 0
			h.first_r2 = update_data.sm_firstroots2;		// 8
			h.fbp1 = fb_p->root1;							// 16
			h.fbp2 = fb_p->root2;							// 24
			h.fbn1 = fb_n->root1;							// 32
			h.fbn2 = fb_n->root2;							// 40
			h.primes = fb_p->prime;							// 48
			h.updates = sm_ptr;								// 56
			h.start = sconf->factor_base->fb_10bit_B;		// 64
			h.stop = sconf->factor_base->fb_15bit_B;		// 68
			if ((h.stop - 8) > h.start)
				h.stop -= 8;

			COMPUTE_8X_SMALL_NROOTS;
			
			j = h.stop;
		}
		sm_ptr = &dconf->sm_rootupdates[(v-1) * med_B + j];
		

#else

		ptr = &dconf->rootupdates[(v-1) * bound + sconf->factor_base->fb_10bit_B];
		for (j=sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}
#endif		

		// assembly code may not get all the way to 15 bits since we 
		// do things in blocks of 8 there.  Make sure we are at the 15 bit
		// boundary before we switch to using update_data.firstroots1/2.
		// this should only run a few iterations, if any.
		ptr = &dconf->rootupdates[(v-1) * bound + j];
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = (uint16_t)update_data.sm_firstroots1[j];
			root2 = (uint16_t)update_data.sm_firstroots2[j];

			if ((prime > 32768) && ((j&7) == 0))
				break;

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}	

		// continue one at a time once we exceed 15 bits, because the 8x SSE2
		// code has a hard time with unsigned 16 bit comparisons
		for ( ; j < med_B; j++, ptr++)
		{
			prime = update_data.prime[j];
			root1 = update_data.sm_firstroots1[j];
			root2 = update_data.sm_firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			if (root2 < root1)
			{
				update_data.sm_firstroots1[j] = (uint16_t)root2;
				update_data.sm_firstroots2[j] = (uint16_t)root1;

				fb_p->root1[j] = (uint16_t)root2;
				fb_p->root2[j] = (uint16_t)root1;
				fb_n->root1[j] = (uint16_t)(prime - root1);
				fb_n->root2[j] = (uint16_t)(prime - root2);
			}
			else
			{
				update_data.sm_firstroots1[j] = (uint16_t)root1;
				update_data.sm_firstroots2[j] = (uint16_t)root2;

				fb_p->root1[j] = (uint16_t)root1;
				fb_p->root2[j] = (uint16_t)root2;
				fb_n->root1[j] = (uint16_t)(prime - root2);
				fb_n->root2[j] = (uint16_t)(prime - root1);
			}
		}	

		bound_index = 0;
		bound_val = med_B;
		check_bound = med_B + BUCKET_ALLOC/2;
		
		
#if 0 // defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		 


#else

		logp = update_data.logp[j-1];
		for (j=med_B;j<large_B;j++,ptr++)
		{
            uint32_t nroot1, nroot2;

			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];
			root1 = update_data.firstroots1[j];
			root2 = update_data.firstroots2[j];

			COMPUTE_NEXT_ROOTS_N;

			update_data.firstroots1[j] = root1;
			update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
		}

#endif
	
#if 0// defined(USE_POLY_SSE2_ASM) && defined(GCC_ASM64X) && !defined(PROFILING)
		 
#else

		logp = update_data.logp[j-1];
		for (j=large_B;j<bound;j++,ptr++)				
		{				
			CHECK_NEW_SLICE(j);

			prime = update_data.prime[j];			
			root1 = update_data.firstroots1[j];	
			root2 = update_data.firstroots2[j];	

			COMPUTE_NEXT_ROOTS_N;		

			update_data.firstroots1[j] = root1;	
			update_data.firstroots2[j] = root2;	

			FILL_ONE_PRIME_P(j);	

			root1 = (prime - root1);		
			root2 = (prime - root2);	
			
			FILL_ONE_PRIME_N(j);
		}

#endif

	}

	if (lp_bucket_p->list != NULL)
	{
		lp_bucket_p->num_slices = bound_index + 1;
		lp_bucket_p->logp[bound_index] = logp;
	}

	return;
}

void nextRoots_32k_generic(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
    int *rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;
    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int *ptr;
    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t large_B = sconf->factor_base->large_B;
    uint32_t j, interval; //, fb_offset;
    int k, numblocks;
    uint32_t root1, root2, nroot1, nroot2, prime;
    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32_t bound_val = med_B;
    uint32_t *numptr_p, *numptr_n, *sliceptr_p, *sliceptr_n;
    uint8_t* slicelogp_ptr = NULL;
    uint32_t* slicebound_ptr = NULL;
    uint32_t *bptr;
    int bnum, room;
    uint8_t logp = 0;

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {
        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
            numptr_p[j] = 0;

        lp_bucket_p->num_slices = 0;

        slicelogp_ptr = lp_bucket_p->logp;
        slicebound_ptr = lp_bucket_p->fb_bounds;

    }
    else
    {
        sliceptr_p = NULL;
        sliceptr_n = NULL;
        numptr_p = NULL;
        numptr_n = NULL;
    }

    k = 0;
    ptr = &rootupdates[(v - 1) * bound + startprime];

    if (sign > 0)
    {
        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16_t)root1;
            fb_p->root2[j] = (uint16_t)root2;
            fb_n->root1[j] = (uint16_t)(prime - root2);
            fb_n->root2[j] = (uint16_t)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        // do one at a time up to the 10bit boundary, where
        // we can start doing things 8 at a time and be
        // sure we can use aligned moves (static_data_init).		
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32_t)update_data.sm_firstroots1[j];
            root2 = (uint32_t)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        // assembly code may not get all the way to 15 bits since we 
        // do things in blocks of 8 there.  Make sure we are at the 15 bit
        // boundary before we switch to using update_data.firstroots1/2.
        // this should only run a few iterations, if any.
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16_t)update_data.sm_firstroots1[j];
            root2 = (uint16_t)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
                break;

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        // continue one at a time once we exceed 15 bits, because the 8x SSE2
        // code has a hard time with unsigned 16 bit comparisons
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;



#if defined(D_HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = med_B; j<large_B;)
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_PROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;
        }


#else

        ptr = &rootupdates[(v - 1) * bound + med_B];
        logp = update_data.logp[med_B - 1];
        for (j = med_B; j<large_B; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }

#endif

#if defined(D_HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = large_B; j<bound;)
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_PROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;
        }


#else

        ptr = &rootupdates[(v - 1) * bound + large_B];
        logp = update_data.logp[large_B - 1];
        for (j = large_B; j<bound; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);
        }

#endif

    }
    else
    {
        for (j = startprime; j<sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16_t)root1;
            fb_p->root2[j] = (uint16_t)root2;
            fb_n->root1[j] = (uint16_t)(prime - root2);
            fb_n->root2[j] = (uint16_t)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        // do one at a time up to the 10bit boundary, where
        // we can start doing things 8 at a time and be
        // sure we can use aligned moves (static_data_init).	
        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->fb_10bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32_t)update_data.sm_firstroots1[j];
            root2 = (uint32_t)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        ptr = &dconf->rootupdates[(v - 1) * bound + sconf->factor_base->fb_10bit_B];
        for (j = sconf->factor_base->fb_10bit_B; j < sconf->factor_base->fb_15bit_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        // assembly code may not get all the way to 15 bits since we 
        // do things in blocks of 8 there.  Make sure we are at the 15 bit
        // boundary before we switch to using update_data.firstroots1/2.
        // this should only run a few iterations, if any.
        ptr = &dconf->rootupdates[(v - 1) * bound + j];
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint16_t)update_data.sm_firstroots1[j];
            root2 = (uint16_t)update_data.sm_firstroots2[j];

            if ((prime > 32768) && ((j & 7) == 0))
                break;

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        // continue one at a time once we exceed 15 bits, because the 8x SSE2
        // code has a hard time with unsigned 16 bit comparisons
        for (; j < med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.sm_firstroots1[j];
            root2 = update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

        bound_index = 0;
        bound_val = med_B;
        check_bound = med_B + BUCKET_ALLOC / 2;


#if defined(D_HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = med_B; j<large_B;)
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_NROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];
            nroot1 = (prime - root1);
            nroot2 = (prime - root2);
            

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);

            j++;
        }


#else

        ptr = &rootupdates[(v - 1) * bound + med_B];
        logp = update_data.logp[med_B - 1];
        for (j = med_B; j<large_B; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;
            nroot1 = prime - root1;
            nroot2 = prime - root2;

            FILL_ONE_PRIME_LOOP_P(j);
            FILL_ONE_PRIME_LOOP_N(j);
        }

#endif

#if defined(D_HAS_SSE2)

        logp = update_data.logp[j - 1];
        for (j = large_B; j<bound;)
        {
            CHECK_NEW_SLICE(j);

            COMPUTE_4_NROOTS(j);

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;

            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];
            prime = update_data.prime[j];

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);

            j++;
        }


#else

        ptr = &rootupdates[(v - 1) * bound + large_B];
        logp = update_data.logp[large_B - 1];
        for (j = large_B; j<bound; j++, ptr++)
        {
            CHECK_NEW_SLICE(j);

            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            FILL_ONE_PRIME_P(j);

            root1 = (prime - root1);
            root2 = (prime - root2);

            FILL_ONE_PRIME_N(j);
        }

#endif

    }

    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
    }

    return;
}

void nextRoots_32k_generic_small(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    //update the roots 
    sieve_fb_compressed *fb_p = dconf->comp_sieve_p;
    sieve_fb_compressed *fb_n = dconf->comp_sieve_n;
    int *rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;
    uint32_t startprime = 2;
    uint32_t bound = sconf->factor_base->B;
    char v = dconf->curr_poly->nu[dconf->numB];
    char sign = dconf->curr_poly->gray[dconf->numB];
    int *ptr;
    uint32_t j;
    uint32_t root1, root2, prime;

    ptr = &rootupdates[(v - 1) * bound + startprime];
    //ptr = rootupdates;

    if (sign > 0)
    {
        for (j = startprime; j < sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16_t)root1;
            fb_p->root2[j] = (uint16_t)root2;
            fb_n->root1[j] = (uint16_t)(prime - root2);
            fb_n->root2[j] = (uint16_t)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32_t)update_data.sm_firstroots1[j];
            root2 = (uint32_t)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_P;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }
    }
    else
    {
        for (j = startprime; j < sconf->sieve_small_fb_start; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = update_data.firstroots1[j];
            root2 = update_data.firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            //we don't sieve these, so ordering doesn't matter
            update_data.firstroots1[j] = root1;
            update_data.firstroots2[j] = root2;

            fb_p->root1[j] = (uint16_t)root1;
            fb_p->root2[j] = (uint16_t)root2;
            fb_n->root1[j] = (uint16_t)(prime - root2);
            fb_n->root2[j] = (uint16_t)(prime - root1);
            if (fb_n->root1[j] == prime)
                fb_n->root1[j] = 0;
            if (fb_n->root2[j] == prime)
                fb_n->root2[j] = 0;

        }

        for (j = sconf->sieve_small_fb_start;
            j < sconf->factor_base->med_B; j++, ptr++)
        {
            prime = update_data.prime[j];
            root1 = (uint32_t)update_data.sm_firstroots1[j];
            root2 = (uint32_t)update_data.sm_firstroots2[j];

            COMPUTE_NEXT_ROOTS_N;

            if (root2 < root1)
            {
                update_data.sm_firstroots1[j] = (uint16_t)root2;
                update_data.sm_firstroots2[j] = (uint16_t)root1;

                fb_p->root1[j] = (uint16_t)root2;
                fb_p->root2[j] = (uint16_t)root1;
                fb_n->root1[j] = (uint16_t)(prime - root1);
                fb_n->root2[j] = (uint16_t)(prime - root2);
            }
            else
            {
                update_data.sm_firstroots1[j] = (uint16_t)root1;
                update_data.sm_firstroots2[j] = (uint16_t)root2;

                fb_p->root1[j] = (uint16_t)root1;
                fb_p->root2[j] = (uint16_t)root2;
                fb_n->root1[j] = (uint16_t)(prime - root2);
                fb_n->root2[j] = (uint16_t)(prime - root1);
            }
        }

    }


    return;
}

void nextRoots_32k_generic_polybatch(static_conf_t *sconf, dynamic_conf_t *dconf)
{
    int *rootupdates = dconf->rootupdates;
    update_t update_data = dconf->update_data;

    uint32_t bound = sconf->factor_base->B;
    char *nu = dconf->curr_poly->nu;
    char *gray = dconf->curr_poly->gray;
    int numB = dconf->numB;
    uint32_t poly_offset = 2 * sconf->num_blocks * dconf->buckets->alloc_slices;

    lp_bucket *lp_bucket_p = dconf->buckets;
    uint32_t med_B = sconf->factor_base->med_B;
    uint32_t large_B = sconf->factor_base->large_B;

    uint32_t j, interval;
    int k, numblocks;
    uint32_t root1, root2, nroot1, nroot2, prime;

    int bound_index = 0;
    int check_bound = BUCKET_ALLOC / 2 - 1;
    uint32_t bound_val = med_B;
    uint32_t *numptr_p, *numptr_n, *sliceptr_p, *sliceptr_n;
    uint8_t* slicelogp_ptr;
    uint32_t* slicebound_ptr;

    uint32_t *bptr;
    int bnum, room;
    uint8_t logp = 0;

    numblocks = sconf->num_blocks;
    interval = numblocks << 15;

    if (lp_bucket_p->alloc_slices != 0) // != NULL)
    {
        lp_bucket_p->fb_bounds[0] = med_B;

        sliceptr_p = lp_bucket_p->list;
        sliceptr_n = lp_bucket_p->list + (numblocks << BUCKET_BITS);

        numptr_p = lp_bucket_p->num;
        numptr_n = lp_bucket_p->num + numblocks;

        // reset bucket counts
        for (j = 0; j < lp_bucket_p->list_size; j++)
            numptr_p[j] = 0;

        lp_bucket_p->num_slices = 0;

        slicelogp_ptr = lp_bucket_p->logp;
        slicebound_ptr = lp_bucket_p->fb_bounds;

    }
    else
    {
        sliceptr_p = NULL;
        sliceptr_n = NULL;
        numptr_p = NULL;
        numptr_n = NULL;
    }
   
    bound_index = 0;
    bound_val = med_B;
    check_bound = med_B + BUCKET_ALLOC / 2;

    logp = update_data.logp[med_B];
    for (j = med_B; j < large_B; j += 16)
    {
        int p;

        CHECK_NEW_SLICE_BATCH(j);

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            if (gray[numB + p] > 0)
            {
                for (k = 0; k < 16; k++)
                {
                    prime = update_data.prime[j + k];
                    root1 = update_data.firstroots1[j + k];
                    root2 = update_data.firstroots2[j + k];

                    COMPUTE_NEXT_ROOTS_BATCH_P(p);

                    update_data.firstroots1[j + k] = root1;
                    update_data.firstroots2[j + k] = root2;
                    nroot1 = prime - root1;
                    nroot2 = prime - root2;
                    FILL_ONE_PRIME_LOOP_P(j + k);
                    FILL_ONE_PRIME_LOOP_N(j + k);
                }
            }
            else
            {
                for (k = 0; k < 16; k++)
                {
                    prime = update_data.prime[j + k];
                    root1 = update_data.firstroots1[j + k];
                    root2 = update_data.firstroots2[j + k];

                    COMPUTE_NEXT_ROOTS_BATCH_N(p);

                    update_data.firstroots1[j + k] = root1;
                    update_data.firstroots2[j + k] = root2;
                    nroot1 = prime - root1;
                    nroot2 = prime - root2;
                    FILL_ONE_PRIME_LOOP_P(j + k);
                    FILL_ONE_PRIME_LOOP_N(j + k);
                }
            }

            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;
        }

        // reset pointers
        sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
        sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
        numptr_p -= p * poly_offset;
        numptr_n -= p * poly_offset;

    }

    logp = update_data.logp[j - 1];
    for (j = large_B; j < bound; j += 16)
    {
        int p;

        CHECK_NEW_SLICE_BATCH(j);

        for (p = 0; (p < dconf->poly_batchsize) && ((numB + p) < dconf->maxB); p++)
        {
            if (gray[numB + p] > 0)
            {
                for (k = 0; k < 16; k++)
                {
                    prime = update_data.prime[j + k];
                    root1 = update_data.firstroots1[j + k];
                    root2 = update_data.firstroots2[j + k];

                    COMPUTE_NEXT_ROOTS_BATCH_P(p);

                    update_data.firstroots1[j + k] = root1;
                    update_data.firstroots2[j + k] = root2;

                    FILL_ONE_PRIME_P(j + k);
                    root1 = prime - root1;
                    root2 = prime - root2;
                    FILL_ONE_PRIME_N(j + k);
                }
            }
            else
            {
                for (k = 0; k < 16; k++)
                {
                    prime = update_data.prime[j + k];
                    root1 = update_data.firstroots1[j + k];
                    root2 = update_data.firstroots2[j + k];

                    COMPUTE_NEXT_ROOTS_BATCH_N(p);

                    update_data.firstroots1[j + k] = root1;
                    update_data.firstroots2[j + k] = root2;

                    FILL_ONE_PRIME_P(j + k);
                    root1 = prime - root1;
                    root2 = prime - root2;
                    FILL_ONE_PRIME_N(j + k);
                }
            }

            // advance pointers
            sliceptr_p += poly_offset * BUCKET_ALLOC;
            sliceptr_n += poly_offset * BUCKET_ALLOC;
            numptr_p += poly_offset;
            numptr_n += poly_offset;

        }

        // reset pointers
        sliceptr_p -= p * poly_offset * BUCKET_ALLOC;
        sliceptr_n -= p * poly_offset * BUCKET_ALLOC;
        numptr_p -= p * poly_offset;
        numptr_n -= p * poly_offset;

    }


    if (lp_bucket_p->list != NULL)
    {
        lp_bucket_p->num_slices = bound_index + 1;
        lp_bucket_p->logp[bound_index] = logp;
    }

    return;
}
