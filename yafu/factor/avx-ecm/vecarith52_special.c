/*
Copyright (c) 2014, Ben Buhrow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.


Copyright (c) 2018 by The Mayo Clinic, though its Special Purpose
 Processor Development Group (SPPDG). All Rights Reserved Worldwide.
 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain
 a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.
Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied,
           including conditions of title, non-infringement, merchantability,
          or fitness for a particular purpose
 See the License for the specific language governing permissions and
 limitations under the License.
This file is a snapshot of a work in progress, originated by Mayo
 Clinic SPPDG.
*/

#include "avx_ecm.h"

#ifdef _MSC_VER
//#define USE_AVX512F
#endif

#ifdef USE_AVX512F
#include <immintrin.h>

//#define PRINT_DEBUG 2
//#define NPRINT_DEBUG 0

// test cases:
// is stage 2 broken for mersennes? - yes looks that way

// ./yafu "ecm(2^823-1,8)" -v -v -B1ecm 3000000 -sigma 4248094344 finds factor 1460915248436556406607
// ./yafu "ecm(2^823-1,8)" -v -v -B1ecm 3000000 -sigma 2102985689 finds both factors 1460915248436556406607 and 1534086200463688788034864584049
// ./yafu "ecm((2^1063+1)/3,8)" -v -v -B1ecm 3000000 -sigma 3299028348 finds both factors 114584129081 and 26210488518118323164267329859
// ./yafu "ecm((2^827+1)/3,8)" -v -v -B1ecm 28000000 -B2ecm 760000000 -sigma 4029008539 fails to find factor!
// ./yafu "ecm((2^773-1)/6864241/9461521,8)" -v -v -B1ecm 2500000 -B2ecm 1800000000 -sigma 8170945836124664 fails to find factor!

// need to fix
// ./yafu "ecm((2^843+1)/9,8)" -v -v -B1ecm 1000000
// Mersenne input 2^843 - 1 determined to be faster by REDC
// conflicting parameters:
// Choosing MAXBITS = 1040, NWORDS = 12, NBLOCKS = 3 based on input size 561
// then crashes.

void vecmulmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecmulmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half mul
    //for (i = 0; i < NBLOCKS; i++)
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);


        // now, do a carry-propagating column summation and store the results.
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c00 = a3;
        c01 = a2;
        c02 = a1;
        c03 = a0;
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);


        // now, do a carry-propagating column summation and store the results.
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c04 = a3;
        c05 = a2;
        c06 = a1;
        c07 = a0;
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half mul
    //for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    i = 2;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);

        // final column accumulation
        {
            ACCUM_SHIFT(a3, te0, te1);
            ACCUM_SHIFT(a2, te2, te3);
            ACCUM_SHIFT(a1, te4, te5);
            ACCUM_SHIFT(a0, te6, te7);

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

            c08 = a3;
            c09 = a2;
            c10 = a1;
            c11 = a0;

        }

    }

    i = 3;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }

        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);

        // final column accumulation
        {
            ACCUM_SHIFT(a3, te0, te1);
            ACCUM_SHIFT(a2, te2, te3);
            ACCUM_SHIFT(a1, te4, te5);
            ACCUM_SHIFT(a0, te6, te7);

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            // i = 3, a3..0 = c[0..3]
            // i = 4, a3..0 = c[4..7]
            // i = 5, a3..0 = c[8..11]
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            //_mm512_store_epi64(c->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

            c12 = a3;
            c13 = a2;
            c14 = a1;
            c15 = a0;

        }

    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        // for different lanes computing with different N.
        // right now that is only for avxppm1
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        _mm512_store_epi64(s->data + 0 * VECLEN, c00);
        _mm512_store_epi64(s->data + 1 * VECLEN, c01);
        _mm512_store_epi64(s->data + 2 * VECLEN, c02);
        _mm512_store_epi64(s->data + 3 * VECLEN, c03);
        _mm512_store_epi64(s->data + 4 * VECLEN, c04);
        _mm512_store_epi64(s->data + 5 * VECLEN, c05);
        _mm512_store_epi64(s->data + 6 * VECLEN, c06);
        _mm512_store_epi64(s->data + 7 * VECLEN, c07);
        _mm512_store_epi64(s->data + 8 * VECLEN, c08);
        _mm512_store_epi64(s->data + 9 * VECLEN, c09);
        _mm512_store_epi64(s->data + 10 * VECLEN, c10);
        _mm512_store_epi64(s->data + 11 * VECLEN, c11);
        _mm512_store_epi64(s->data + 12 * VECLEN, c12);
        _mm512_store_epi64(s->data + 13 * VECLEN, c13);
        _mm512_store_epi64(s->data + 14 * VECLEN, c14);
        _mm512_store_epi64(s->data + 15 * VECLEN, c15);

        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);


        //int wshift = mdata->nbits / 52;
        //int bshift = mdata->nbits % 52;
        //__m512i lowmask = _mm512_set1_epi64((1ULL << (uint64_t)bshift) - 1ULL);
        //
        //// split off the high half into the top words.  
        //c15 = _mm512_or_epi64(_mm512_slli_epi64(c15, (DIGITBITS - bshift)), _mm512_srli_epi64(c14, bshift));
        //c14 = _mm512_or_epi64(_mm512_slli_epi64(c14, (DIGITBITS - bshift)), _mm512_srli_epi64(c13, bshift));
        //c13 = _mm512_or_epi64(_mm512_slli_epi64(c13, (DIGITBITS - bshift)), _mm512_srli_epi64(c12, bshift));
        //c12 = _mm512_or_epi64(_mm512_slli_epi64(c12, (DIGITBITS - bshift)), _mm512_srli_epi64(c11, bshift));
        //c11 = _mm512_or_epi64(_mm512_slli_epi64(c11, (DIGITBITS - bshift)), _mm512_srli_epi64(c10, bshift));
        //c10 = _mm512_or_epi64(_mm512_slli_epi64(c10, (DIGITBITS - bshift)), _mm512_srli_epi64(c09, bshift));
        //c09 = _mm512_or_epi64(_mm512_slli_epi64(c09, (DIGITBITS - bshift)), _mm512_srli_epi64(c08, bshift));
        //c08 = _mm512_or_epi64(_mm512_slli_epi64(c08, (DIGITBITS - bshift)), _mm512_srli_epi64(c07, bshift));
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        //c00 = _mm512_adc_epi52(c00, scarry, c08, &scarry);
        //c01 = _mm512_adc_epi52(c01, scarry, c09, &scarry);
        //c02 = _mm512_adc_epi52(c02, scarry, c10, &scarry);
        //c03 = _mm512_adc_epi52(c03, scarry, c11, &scarry);
        //c04 = _mm512_adc_epi52(c04, scarry, c12, &scarry);
        //c05 = _mm512_adc_epi52(c05, scarry, c13, &scarry);
        //c06 = _mm512_adc_epi52(c06, scarry, c14, &scarry);
        //c07 = _mm512_adc_epi52(c07, scarry, c15, &scarry);

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);

        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        //scarry = _mm512_test_epi64_mask(c07, vbpshift);
        //c00 = _mm512_addcarry_epi52(c00, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c01 = _mm512_addcarry_epi52(c01, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c02 = _mm512_addcarry_epi52(c02, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c03 = _mm512_addcarry_epi52(c03, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c04 = _mm512_addcarry_epi52(c04, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c05 = _mm512_addcarry_epi52(c05, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c06 = _mm512_addcarry_epi52(c06, scarry, &scarry); if (scarry == 0) goto carryadddone;
        //c07 = _mm512_addcarry_epi52(c07, scarry, &scarry); if (scarry == 0) goto carryadddone;

    carryadddone:


        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        //c07 = _mm512_and_epi64(vbshift, c07);

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        //_mm512_store_epi64(c->data + 0 * VECLEN, c00);
        //_mm512_store_epi64(c->data + 1 * VECLEN, c01);
        //_mm512_store_epi64(c->data + 2 * VECLEN, c02);
        //_mm512_store_epi64(c->data + 3 * VECLEN, c03);
        //_mm512_store_epi64(c->data + 4 * VECLEN, c04);
        //_mm512_store_epi64(c->data + 5 * VECLEN, c05);
        //_mm512_store_epi64(c->data + 6 * VECLEN, c06);
        //_mm512_store_epi64(c->data + 7 * VECLEN, c07);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mersenne_416(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    vec_bignum_t* b = a;

    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19
    __m512i c00, c01, c02, c03, c04, c05, c06, c07,
        c08, c09, c10, c11, c12, c13, c14, c15;

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecsqrmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a0, a2, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, a3, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a1, a2, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a0, a1);    // te2
            //prod2_e = _mm512_mul_epu32(a0, a2);    // te4
            //prod3_e = _mm512_mul_epu32(a0, a3);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);
                prod3_ld = _mm512_cvtepu64_pd(a2);
                prod4_ld = _mm512_cvtepu64_pd(a3);

                prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a3 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }

            //prod1_e = _mm512_mul_epu32(a1, a2);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a2 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
            SUB_BIAS_LO(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a0, a0);
            //prod2_e = _mm512_mul_epu32(a1, a1);
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te0/1
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            }
#endif

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c00 = a3;
        c01 = a2;
        c02 = a1;
        c03 = a0;
    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        {
            // i odd
            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
            // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
            // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a2, b2);   // te0
            //prod2_e = _mm512_mul_epu32(a1, b2);   // te2
            //prod3_e = _mm512_mul_epu32(a1, b3);   // te4
            //prod4_e = _mm512_mul_epu32(a0, b3);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a0);
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);
                prod3_ld = _mm512_cvtepu64_pd(b2);
                prod4_ld = _mm512_cvtepu64_pd(b3);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1 
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b3, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a2, b4, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b3);   // te0
            //prod2_e = _mm512_mul_epu32(a2, b3);   // te2
            //prod3_e = _mm512_mul_epu32(a2, b4);   // te4
            //prod4_e = _mm512_mul_epu32(a1, b4);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a1);
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b3);
                prod4_ld = _mm512_cvtepu64_pd(b4);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b3 -> to te0/1
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b4 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b4, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a3, b5, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a2, b5, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a3, b6, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b4);  // te2
            //prod2_e = _mm512_mul_epu32(a3, b5);  // te4
            //prod3_e = _mm512_mul_epu32(a2, b5);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b4);
                prod4_ld = _mm512_cvtepu64_pd(b5);

                prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b4 -> to te2/3
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b5 -> to te4/5
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b5 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            //prod1_e = _mm512_mul_epu32(a3, b6);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a3);
                prod2_ld = _mm512_cvtepu64_pd(b6);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
            SUB_BIAS_LO(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
            //prod2_e = _mm512_mul_epu32(a0, a0);    // te4
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te0/1
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
            }
#endif

        }


#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c04 = a3;
        c05 = a2;
        c06 = a1;
        c07 = a0;
    }
#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half sqr
    i = 0;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is even
        {
            // i even, block shape 1.
            {

                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                //prod1_e = _mm512_mul_epu32(a0, b0);    // te0
                //prod2_e = _mm512_mul_epu32(a0, b1);    // te2
                //prod3_e = _mm512_mul_epu32(a0, b2);    // te4
                //prod4_e = _mm512_mul_epu32(a0, b3);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);

                //prod1_e = _mm512_mul_epu32(a1, b1);    // te0
                //prod2_e = _mm512_mul_epu32(a1, b2);    // te2
                //prod3_e = _mm512_mul_epu32(a1, b3);    // te4
                //prod4_e = _mm512_mul_epu32(a1, b4);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);    // te0
                //prod2_e = _mm512_mul_epu32(a2, b3);    // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
                SUB_BIAS_LO(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);   // te0
                //prod2_e = _mm512_mul_epu32(a2, a2);   // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c08 = a3;
        c09 = a2;
        c10 = a1;
        c11 = a0;


    }

    i = 1;
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        // NBLOCKS is even
        {
            // i odd, block shape 1.
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);   // te0
                //prod2_e = _mm512_mul_epu32(a0, a1);   // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
                SUB_BIAS_LO(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);       // te0
                //prod2_e = _mm512_mul_epu32(a0, a0);       // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        ACCUM_SHIFT(a3, te0, te1);
        ACCUM_SHIFT(a2, te2, te3);
        ACCUM_SHIFT(a1, te4, te5);
        ACCUM_SHIFT(a0, te6, te7);

        a3 = _mm512_and_epi64(vlmask, a3);
        a2 = _mm512_and_epi64(vlmask, a2);
        a1 = _mm512_and_epi64(vlmask, a1);
        a0 = _mm512_and_epi64(vlmask, a0);

        // i = 3, a3..0 = c[0..3]
        // i = 4, a3..0 = c[4..7]
        // i = 5, a3..0 = c[8..11]
        c12 = a3;
        c13 = a2;
        c14 = a1;
        c15 = a0;


    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        _mm512_store_epi64(s->data + 0 * VECLEN, c00);
        _mm512_store_epi64(s->data + 1 * VECLEN, c01);
        _mm512_store_epi64(s->data + 2 * VECLEN, c02);
        _mm512_store_epi64(s->data + 3 * VECLEN, c03);
        _mm512_store_epi64(s->data + 4 * VECLEN, c04);
        _mm512_store_epi64(s->data + 5 * VECLEN, c05);
        _mm512_store_epi64(s->data + 6 * VECLEN, c06);
        _mm512_store_epi64(s->data + 7 * VECLEN, c07);
        _mm512_store_epi64(s->data + 8 * VECLEN, c08);
        _mm512_store_epi64(s->data + 9 * VECLEN, c09);
        _mm512_store_epi64(s->data + 10 * VECLEN, c10);
        _mm512_store_epi64(s->data + 11 * VECLEN, c11);
        _mm512_store_epi64(s->data + 12 * VECLEN, c12);
        _mm512_store_epi64(s->data + 13 * VECLEN, c13);
        _mm512_store_epi64(s->data + 14 * VECLEN, c14);
        _mm512_store_epi64(s->data + 15 * VECLEN, c15);

        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul.  It's possible these hi bits
        // are located in two words
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        b0 = _mm512_set1_epi64(mdata->isMersenne);
        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        if (scarry > 0)
        {
            // it is unlikely, but possible that we still have a final carry to propagate
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        // zero any remaining lo words.
        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after sub:");
#endif

        // if there was a carry, add 1 to the lo bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after bit clear:");
        exit(1);
#endif

    }

    c->size = NWORDS;
    return;
}

void vecmulmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    int i, j;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;

    //if (NWORDS == 8)
    //{
    //    vecmulmod52_mersenne_416(a, b, c, n, s, mdata);
    //    return;
    //}


    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = a->signmask ^ b->signmask;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecmulmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half mul
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i; j > 0; j--)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // finish each triangluar shaped column sum
        a0 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 0) * VECLEN);
        a1 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + (i * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + 3 * VECLEN);
        b1 = _mm512_load_epi64(b->data + 2 * VECLEN);
        b2 = _mm512_load_epi64(b->data + 1 * VECLEN);
        b3 = _mm512_load_epi64(b->data + 0 * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a0, b3, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b3, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b2, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b3, te4, te5);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod1_hd));  // a1 * b3 -> to te2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));  // a2 * b3 -> to te4/5 
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a0 * b2 -> to te2/3
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a0 * b3 -> to te0/1

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod1_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b2, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b1, te4, te5);
        //VEC_MUL_ACCUM_LOHI_PD(a1, b1, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a0, b0, te6, te7);
        {
            prod5_ld = _mm512_cvtepu64_pd(a0);
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(b0);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b2);

            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b3, te6, te7);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b2, te6, te7);
        {
            prod1_ld = _mm512_cvtepu64_pd(a2);
            prod2_ld = _mm512_cvtepu64_pd(a3);
            prod3_ld = _mm512_cvtepu64_pd(b2);
            prod4_ld = _mm512_cvtepu64_pd(b3);

            prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
            te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod2_hd));

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod2_ld));
        }

        // subtract out all of the bias at once.  these are
        // counts of how many times we have mul_accumulated into each column
        // during the block a*b and final triangular a*b.
        SUB_BIAS_HI(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
        SUB_BIAS_LO(
            i * 4 + 1,
            i * 4 + 2,
            i * 4 + 3,
            i * 4 + 4);
#endif



        // now, do a carry-propagating column summation and store the results.
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half mul
    for (i = NBLOCKS; i < 2 * NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (j = i - NBLOCKS + 1; j < NBLOCKS; j++)
        {
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }


        // finish each triangluar shaped column sum (a * b)
        a1 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 1) * VECLEN);
        a2 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 2) * VECLEN);
        a3 = _mm512_load_epi64(a->data + ((i - NBLOCKS) * BLOCKWORDS + 3) * VECLEN);

        b0 = _mm512_load_epi64(b->data + (NWORDS - 1) * VECLEN);
        b1 = _mm512_load_epi64(b->data + (NWORDS - 2) * VECLEN);
        b2 = _mm512_load_epi64(b->data + (NWORDS - 3) * VECLEN);

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
#else
        // ======
        //VEC_MUL_ACCUM_LOHI_PD(a1, b0, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b1, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a2, b0, te2, te3);
        {
            prod1_ld = _mm512_cvtepu64_pd(a1);
            prod2_ld = _mm512_cvtepu64_pd(a2);
            prod3_ld = _mm512_cvtepu64_pd(b0);
            prod4_ld = _mm512_cvtepu64_pd(b1);

            prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));  // a1*b0 -> t0/1
            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod4_hd));  // a2*b1 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a2*b0 -> t2/3

            prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod4_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
        }
#endif

#ifdef IFMA
        VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
#else
        //VEC_MUL_ACCUM_LOHI_PD(a3, b2, te0, te1);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b1, te2, te3);
        //VEC_MUL_ACCUM_LOHI_PD(a3, b0, te4, te5);
        {
            prod1_ld = _mm512_cvtepu64_pd(a3);
            prod2_ld = _mm512_cvtepu64_pd(b2);
            prod3_ld = _mm512_cvtepu64_pd(b1);
            prod4_ld = _mm512_cvtepu64_pd(b0);

            prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));  // a3*b2 -> t0/1
            te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));  // a3*b1 -> t2/3
            te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod4_hd));  // a3*b0 -> t4/5

            prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
            prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
            prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

            prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
            prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

            te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
            te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
            te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod4_ld));
        }


        // subtract out all of the bias at once.
        SUB_BIAS_HI(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
        SUB_BIAS_LO(
            (2 * NBLOCKS - i - 1) * 4 + 3,
            (2 * NBLOCKS - i - 1) * 4 + 2,
            (2 * NBLOCKS - i - 1) * 4 + 1,
            (2 * NBLOCKS - i - 1) * 4 + 0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + 3) * VECLEN, a0);

        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);
                
                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);
            
            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        __mmask8 carryout = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo
            
            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            
            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = NWORDS;
    return;
}

void vecsqrmod52_mersenne(vec_bignum_t* a, vec_bignum_t* c, vec_bignum_t* n, vec_bignum_t* s, vec_monty_t* mdata)
{
    // 8x sqr:
    // input 8 bignums in the even lanes of a.
    // output 8 squaremod bignums in the even lanes of c.
    int i, j, k;
    uint32_t NWORDS = mdata->NWORDS;
    uint32_t NBLOCKS = mdata->NBLOCKS;
    vec_bignum_t* b = a;

    //if (NWORDS == 8)
    //{
    //    vecsqrmod52_mersenne_416(a, c, n, s, mdata);
    //    return;
    //}

    // needed in loops
    __m512i a0, a1, a2, a3;                                     // 4
    __m512i b0, b1, b2, b3, b4, b5, b6;                         // 11
    __m512i te0, te1, te2, te3, te4, te5, te6, te7;             // 19

#ifndef IFMA
    __m512d prod1_hd, prod2_hd, prod3_hd, prod4_hd;                 // 23
    __m512d prod1_ld, prod2_ld, prod3_ld, prod4_ld, prod5_ld;        // 28
    __m512d dbias = _mm512_castsi512_pd(_mm512_set1_epi64(0x4670000000000000ULL));
    __m512i vbias1 = _mm512_set1_epi64(0x4670000000000000ULL);  // 31
    __m512i vbias2 = _mm512_set1_epi64(0x4670000000000001ULL);  // 31
    __m512i vbias3 = _mm512_set1_epi64(0x4330000000000000ULL);  // 31
#endif

    // needed after loops
    __m512i vlmask = _mm512_set1_epi64(0x000fffffffffffffULL);
    __m512i acc_e0, acc_e1, acc_e2;
    __m512i zero = _mm512_set1_epi64(0);
    __mmask8 scarry;

    // deal with the sign
    c->size = NWORDS;
    c->signmask = 0;

    // zero the accumulator
    acc_e0 = zero;
    acc_e1 = zero;
    acc_e2 = zero;

#ifdef DEBUG_MERSENNE
    printf("commencing vecsqrmod52_mersenne\n");
    print_vechexbignum(a, "input a: ");
    print_vechexbignum(b, "input b: ");
    print_vechexbignum(s, "input s: ");
#endif

    // first half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = i; j > (i + 1) / 2; j--, k++)
        {
            // when i = 0, j = 0, j > 0 --> 0 iterations
            // when i = 1, j = 1, j > 1 --> 0 iterations
            // when i = 2, j = 2, j > 1 --> 1 iteration @ a[3..0], b[5..11]
            // when i = 3, j = 3, j > 2 --> 1 iteration @ a[3..0], b[9..15]
            // when i = 4, j = 4, j > 2 --> 2 iteration @ a[3..0], b[13..19] and a[7..4], b[9..15]
            // when i = 5, j = 5, j > 3 --> 2 iteration @ a[3..0], b[17..23] and a[7..4], b[13..19]
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        if (i & 1)
        {
            // i odd
            // for 384-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}
            // for 512-bit inputs when i == 3, j = 1, a = {7,6,5,4} and b = {6,7,8,9,a,b}
            // for 512-bit inputs when i == 1, j = 0, a = {3,2,1,0} and b = {2,3,4,5,6,7}

            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);

            b1 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j - 1) * BLOCKWORDS + 7) * VECLEN);


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a2, b2);   // te0
            //prod2_e = _mm512_mul_epu32(a1, b2);   // te2
            //prod3_e = _mm512_mul_epu32(a1, b3);   // te4
            //prod4_e = _mm512_mul_epu32(a0, b3);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a0);
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);
                prod3_ld = _mm512_cvtepu64_pd(b2);
                prod4_ld = _mm512_cvtepu64_pd(b3);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1 
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b3, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a2, b4, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b3);   // te0
            //prod2_e = _mm512_mul_epu32(a2, b3);   // te2
            //prod3_e = _mm512_mul_epu32(a2, b4);   // te4
            //prod4_e = _mm512_mul_epu32(a1, b4);   // te6
            //ACCUM_4X_DOUBLED_PROD;
            {
                prod5_ld = _mm512_cvtepu64_pd(a1);
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b3);
                prod4_ld = _mm512_cvtepu64_pd(b4);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b3 -> to te0/1
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b4 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod5_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }
#endif


#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a3, b4, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a3, b5, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a2, b5, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a3, b6, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a3, b4);  // te2
            //prod2_e = _mm512_mul_epu32(a3, b5);  // te4
            //prod3_e = _mm512_mul_epu32(a2, b5);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a2);
                prod2_ld = _mm512_cvtepu64_pd(a3);
                prod3_ld = _mm512_cvtepu64_pd(b4);
                prod4_ld = _mm512_cvtepu64_pd(b5);

                prod3_hd = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b4 -> to te2/3
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b5 -> to te4/5
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b5 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                prod3_ld = _mm512_fmadd_round_pd(prod2_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod4_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            //prod1_e = _mm512_mul_epu32(a3, b6);  // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a3);
                prod2_ld = _mm512_cvtepu64_pd(b6);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * b6 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
            SUB_BIAS_LO(
                k * 4 + 2,
                k * 4 + 3,
                k * 4 + 3,
                k * 4 + 4);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
            //prod2_e = _mm512_mul_epu32(a0, a0);    // te4
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te0/1
                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
            }
#endif

        }
        else
        {
            // i even
            a0 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 0) * VECLEN);
            a1 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 1) * VECLEN);
            a2 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 2) * VECLEN);
            a3 = _mm512_load_epi64(a->data + ((i - j) * BLOCKWORDS + 3) * VECLEN);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
            VEC_MUL_ACCUM_LOHI_PD(a0, a2, te4, te5);
            VEC_MUL_ACCUM_LOHI_PD(a0, a3, te6, te7);
            VEC_MUL_ACCUM_LOHI_PD(a1, a2, te6, te7);
#else
            //prod1_e = _mm512_mul_epu32(a0, a1);    // te2
            //prod2_e = _mm512_mul_epu32(a0, a2);    // te4
            //prod3_e = _mm512_mul_epu32(a0, a3);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);
                prod3_ld = _mm512_cvtepu64_pd(a2);
                prod4_ld = _mm512_cvtepu64_pd(a3);

                prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3
                prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te4/5
                prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a3 -> to te6/7

                te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod3_hd));
                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));

                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod3_ld));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
            }

            //prod1_e = _mm512_mul_epu32(a1, a2);    // te6
            {
                prod1_ld = _mm512_cvtepu64_pd(a1);
                prod2_ld = _mm512_cvtepu64_pd(a2);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a2 -> to te6/7

                te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod1_hd));
                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod1_ld));
            }

            // all terms so far need to be doubled.  
            // but to do that we first need to remove all bias
            // that has been accumulated so far.
            SUB_BIAS_HI(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
            SUB_BIAS_LO(
                k * 4 + 0,
                k * 4 + 1,
                k * 4 + 1,
                k * 4 + 2);
#endif

            // now double
            te1 = _mm512_slli_epi64(te1, 1);
            te3 = _mm512_slli_epi64(te3, 1);
            te5 = _mm512_slli_epi64(te5, 1);
            te7 = _mm512_slli_epi64(te7, 1);
            te0 = _mm512_slli_epi64(te0, 1);
            te2 = _mm512_slli_epi64(te2, 1);
            te4 = _mm512_slli_epi64(te4, 1);
            te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
            VEC_MUL_ACCUM_LOHI_PD(a0, a0, te0, te1);
            VEC_MUL_ACCUM_LOHI_PD(a1, a1, te4, te5);
#else
            // finally, accumulate the two non-doubled terms.
            //prod1_e = _mm512_mul_epu32(a0, a0);
            //prod2_e = _mm512_mul_epu32(a1, a1);
            {
                prod1_ld = _mm512_cvtepu64_pd(a0);
                prod2_ld = _mm512_cvtepu64_pd(a1);

                prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a0 -> to te0/1
                prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                    (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * a1 -> to te4/5

                te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
            }
#endif

        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum:
            // carry propagate low to high.
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);

            // store the lo word
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + j) * VECLEN, acc_e0);

            // now shift.
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;
        }
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after lo half:");
#endif

    // second half sqr
    for (i = 0; i < NBLOCKS; i++)
    {
        te0 = te1 = te2 = te3 = te4 = te5 = te6 = te7 = zero;

        for (k = 0, j = 0; j < (NBLOCKS - i - 1) / 2; j++, k++)
        {
            // Compute a solid block (all matching terms are in the lower
            // half triangle of the expansion).
            a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
            a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
            a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
            a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

            b0 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 1) * VECLEN);
            b1 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 2) * VECLEN);
            b2 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 3) * VECLEN);
            b3 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 4) * VECLEN);
            b4 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 5) * VECLEN);
            b5 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 6) * VECLEN);
            b6 = _mm512_load_epi64(b->data + ((j + i) * BLOCKWORDS + 7) * VECLEN);

            VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);
            VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);
            VEC_MUL4_ACCUM(a2, b2, b3, b4, b5);
            VEC_MUL4_ACCUM(a3, b3, b4, b5, b6);
        }

        // The final block shape depends on the parity of i and NBLOCKS
        // Block shape 1 (small upper triangle) if 'i' is odd and NBLOCKS is even,
        // or if 'i' is even and NBLOCKS is odd.
        // Block shape 2 if 'i' is odd and NBLOCKS is odd or if 'i' is even
        // and NBLOCKS is even.
        if (NBLOCKS & 1)		// NBLOCKS is odd
        {
            // i odd, block shape 2.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN);
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b0, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a1, b2, te2, te3);
                VEC_MUL_ACCUM_LOHI_PD(a1, b1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, b1, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a0, b0);        // te0
                //prod1_e = _mm512_mul_epu32(a1, b1);        // te0
                //prod1_e = _mm512_mul_epu32(a0, b1);        // te2
                //prod1_e = _mm512_mul_epu32(a1, b2);        // te2
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b0);
                    prod3_ld = _mm512_cvtepu64_pd(b1);
                    prod4_ld = _mm512_cvtepu64_pd(b2);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te0/1
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b2 -> to te2/3
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod4_hd));
                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod4_ld));
                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }
#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);        // te0
                //prod1_e = _mm512_mul_epu32(a2, b3);        // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }
#endif

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, b2, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a1, b4, te6, te7);
                VEC_MUL_ACCUM_LOHI_PD(a1, b3, te4, te5);
                VEC_MUL_ACCUM_LOHI_PD(a0, b3, te6, te7);
#else
                // prod1_e = _mm512_mul_epu32(a0, b2);       // te4
                // prod1_e = _mm512_mul_epu32(a1, b3);       // te4
                // prod1_e = _mm512_mul_epu32(a0, b3);       // te6
                // prod1_e = _mm512_mul_epu32(a1, b4);       // te6
                {
                    prod5_ld = _mm512_cvtepu64_pd(a0);
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);
                    prod4_ld = _mm512_cvtepu64_pd(b4);

                    prod2_hd = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b2 -> to te4/5
                    prod4_hd = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b4 -> to te6/7
                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b3 -> to te4/5
                    prod3_hd = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b3 -> to te6/7

                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod4_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod1_hd));
                    te7 = _mm512_add_epi64(te7, _mm512_castpd_si512(prod3_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);
                    prod4_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod4_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod5_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod4_ld = _mm512_fmadd_round_pd(prod1_ld, prod4_ld, prod4_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod5_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod4_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod1_ld));
                    te6 = _mm512_add_epi64(te6, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);
                SUB_BIAS_LO(
                    k * 4 + 3,
                    k * 4 + 3,
                    k * 4 + 2,
                    k * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);    // te0
                //prod1_e = _mm512_mul_epu32(a2, a2);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif

            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);      // te0
                //prod1_e = _mm512_mul_epu32(a0, a1);      // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);
                SUB_BIAS_LO(
                    k * 4 + 1,
                    k * 4 + 1,
                    k * 4 + 0,
                    k * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);    // te0
                //prod1_e = _mm512_mul_epu32(a0, a0);    // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif

            }

        }
        else				// NBLOCKS is even
        {
            // i odd, block shape 1.
            if (i & 1)
            {
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);  // {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);  // {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);  // {d, 9}

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a0, a2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a1, te2, te3);
#else
                //k == 0;
                //prod1_e = _mm512_mul_epu32(a0, a2);   // te0
                //prod2_e = _mm512_mul_epu32(a0, a1);   // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a0);
                    prod2_ld = _mm512_cvtepu64_pd(a1);
                    prod3_ld = _mm512_cvtepu64_pd(a2);

                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a2 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * a1 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod3_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod2_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod3_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod2_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
                SUB_BIAS_LO(
                    j * 4 + 1,
                    j * 4 + 1,
                    j * 4 + 0,
                    j * 4 + 0);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a1, a1, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a0, a0, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a1, a1);       // te0
                //prod2_e = _mm512_mul_epu32(a0, a0);       // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a1);
                    prod2_ld = _mm512_cvtepu64_pd(a0);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a1 * b1 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a0 * b0 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
            else
            {
                // i even, block shape 1.
                // always a continuation of the full-block loop, so use the same 
                // loading pattern.  Only now we don't need as many b-terms.
                a0 = _mm512_load_epi64(a->data + (NWORDS - 1 - j * BLOCKWORDS) * VECLEN);		// {f, b}
                a1 = _mm512_load_epi64(a->data + (NWORDS - 2 - j * BLOCKWORDS) * VECLEN);		// {e, a}
                a2 = _mm512_load_epi64(a->data + (NWORDS - 3 - j * BLOCKWORDS) * VECLEN);		// {d, 9}
                a3 = _mm512_load_epi64(a->data + (NWORDS - 4 - j * BLOCKWORDS) * VECLEN);		// {c, 8}

                b0 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 1) * VECLEN); // {9, 5}
                b1 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 2) * VECLEN);	// {a, 6}
                b2 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 3) * VECLEN);	// {b, 7}
                b3 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 4) * VECLEN);	// {c, 8}
                b4 = _mm512_load_epi64(b->data + (j * BLOCKWORDS + i * BLOCKWORDS + 5) * VECLEN); // {d, 9}

                //prod1_e = _mm512_mul_epu32(a0, b0);    // te0
                //prod2_e = _mm512_mul_epu32(a0, b1);    // te2
                //prod3_e = _mm512_mul_epu32(a0, b2);    // te4
                //prod4_e = _mm512_mul_epu32(a0, b3);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a0, b0, b1, b2, b3);

                //prod1_e = _mm512_mul_epu32(a1, b1);    // te0
                //prod2_e = _mm512_mul_epu32(a1, b2);    // te2
                //prod3_e = _mm512_mul_epu32(a1, b3);    // te4
                //prod4_e = _mm512_mul_epu32(a1, b4);    // te6
                //ACCUM_4X_DOUBLED_PROD;
                VEC_MUL4_ACCUM(a1, b1, b2, b3, b4);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a2, b2, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, b3, te2, te3);
#else
                //prod1_e = _mm512_mul_epu32(a2, b2);    // te0
                //prod2_e = _mm512_mul_epu32(a2, b3);    // te2
                {
                    prod1_ld = _mm512_cvtepu64_pd(a2);
                    prod2_ld = _mm512_cvtepu64_pd(b2);
                    prod3_ld = _mm512_cvtepu64_pd(b3);

                    prod2_hd = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b2 -> to te0/1
                    prod3_hd = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * b3 -> to te2/3

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod2_hd));
                    te3 = _mm512_add_epi64(te3, _mm512_castpd_si512(prod3_hd));

                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);
                    prod3_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod3_hd);

                    prod2_ld = _mm512_fmadd_round_pd(prod1_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod3_ld = _mm512_fmadd_round_pd(prod1_ld, prod3_ld, prod3_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod2_ld));
                    te2 = _mm512_add_epi64(te2, _mm512_castpd_si512(prod3_ld));
                }

                // all terms so far need to be doubled.  
                // but to do that we first need to remove all bias
                // that has been accumulated so far.
                SUB_BIAS_HI(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
                SUB_BIAS_LO(
                    j * 4 + 3,
                    j * 4 + 3,
                    j * 4 + 2,
                    j * 4 + 2);
#endif

                // now double
                te1 = _mm512_slli_epi64(te1, 1);
                te3 = _mm512_slli_epi64(te3, 1);
                te5 = _mm512_slli_epi64(te5, 1);
                te7 = _mm512_slli_epi64(te7, 1);
                te0 = _mm512_slli_epi64(te0, 1);
                te2 = _mm512_slli_epi64(te2, 1);
                te4 = _mm512_slli_epi64(te4, 1);
                te6 = _mm512_slli_epi64(te6, 1);

#ifdef IFMA
                VEC_MUL_ACCUM_LOHI_PD(a3, a3, te0, te1);
                VEC_MUL_ACCUM_LOHI_PD(a2, a2, te4, te5);
#else
                // finally, accumulate the two non-doubled terms.
                //prod1_e = _mm512_mul_epu32(a3, a3);   // te0
                //prod2_e = _mm512_mul_epu32(a2, a2);   // te4
                {
                    prod1_ld = _mm512_cvtepu64_pd(a3);
                    prod2_ld = _mm512_cvtepu64_pd(a2);

                    prod1_hd = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a3 * a3 -> to te0/1
                    prod2_hd = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, dbias,
                        (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC)); // a2 * a2 -> to te4/5

                    te1 = _mm512_add_epi64(te1, _mm512_castpd_si512(prod1_hd));
                    te5 = _mm512_add_epi64(te5, _mm512_castpd_si512(prod2_hd));

                    prod1_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod1_hd);
                    prod2_hd = _mm512_sub_pd(_mm512_castsi512_pd(vbias2), prod2_hd);

                    prod1_ld = _mm512_fmadd_round_pd(prod1_ld, prod1_ld, prod1_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));
                    prod2_ld = _mm512_fmadd_round_pd(prod2_ld, prod2_ld, prod2_hd, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC));

                    te0 = _mm512_add_epi64(te0, _mm512_castpd_si512(prod1_ld));
                    te4 = _mm512_add_epi64(te4, _mm512_castpd_si512(prod2_ld));
                }
#endif
            }
        }

#ifndef IFMA
        // need to remove bias from the two non-doubled terms of the a*a loop.
        SUB_BIAS_HI(
            1,
            0,
            1,
            0);
        SUB_BIAS_LO(
            1,
            0,
            1,
            0);
#endif

        // final column accumulation
        {
            j = 0;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te0);
            acc_e1 = _mm512_add_epi64(acc_e1, te1);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a3 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN, 
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 1;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te2);
            acc_e1 = _mm512_add_epi64(acc_e1, te3);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a2 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 2;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te4);
            acc_e1 = _mm512_add_epi64(acc_e1, te5);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a1 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            j = 3;
            // accumulate this column-sum
            acc_e0 = _mm512_add_epi64(acc_e0, te6);
            acc_e1 = _mm512_add_epi64(acc_e1, te7);
            VEC_CARRYPROP_LOHI(acc_e0, acc_e1);
            acc_e2 = _mm512_srli_epi64(acc_e1, 52);
            acc_e1 = _mm512_and_epi64(acc_e1, vlmask);
            a0 = acc_e0;

            // store the low-word final result and shift
            //_mm512_store_pd(s->data + (i * BLOCKWORDS + j) * VECLEN,
            //    _mm512_cvtepu64_pd(_mm512_and_epi64(vlmask, acc_e0)));
            acc_e0 = acc_e1;
            acc_e1 = acc_e2;
            acc_e2 = zero;

            a3 = _mm512_and_epi64(vlmask, a3);
            a2 = _mm512_and_epi64(vlmask, a2);
            a1 = _mm512_and_epi64(vlmask, a1);
            a0 = _mm512_and_epi64(vlmask, a0);

            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 0) * VECLEN, a3);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 1) * VECLEN, a2);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 2) * VECLEN, a1);
            _mm512_store_epi64(s->data + (i * BLOCKWORDS + NWORDS + 3) * VECLEN, a0);

        }


    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(s, "after hi half:");
#endif

    // reduce by adding hi to lo.  first right shift hi into output.
    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));

        vec_bignum52_mask_rshift_vn(s, c, mdata->vnbits, 0xff);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
        vec_bignum52_mask_rshift_n(s, c, mdata->nbits, 0xff);
    }

#ifdef DEBUG_MERSENNE
    print_vechexbignum(c, "hi part:");
    print_vechexbignum(s, "lo part:");
#endif

    if (mdata->isMersenne > 1)
    {
        // this number has been identified as a "pseudo"-Mersenne: 2^n-c, with
        // 'c' small.  Fast reduction is still possible.

        // multiply c * hi
        b0 = _mm512_set1_epi64(mdata->isMersenne);
        vecmul52_1(c, b0, c, n, s, mdata);

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after mul_1 with c:");
#endif

        // clear any hi-bits in the lo result
        a1 = _mm512_load_epi64(s->data + wshift * VECLEN);
        _mm512_store_epi64(s->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(s, "cleared low part:");
#endif

        // add c * hi + lo
        scarry = 0;
        for (i = 0; i <= NWORDS; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // now we need to do it again, but the hi part is just a single word so
        // we only have a single-precision vecmul
        //vec_bignum52_mask_rshift_n(c, s, mdata->nbits, 0xff);
        //a0 = _mm512_load_epi64(s->data);

        // It's possible these hi bits are located in two words
        a0 = _mm512_load_epi64(c->data + wshift * VECLEN);
        a0 = _mm512_srli_epi64(a0, bshift);
        a1 = _mm512_load_epi64(c->data + (wshift + 1) * VECLEN);
        a0 = _mm512_or_epi64(a0, _mm512_slli_epi64(a1, 52 - bshift));

        //print_regvechex(a0, 0, "hi part:");

        b0 = _mm512_set1_epi64(mdata->isMersenne);

        //print_regvechex(b0, 0, "c:");

        VEC_MUL_LOHI_PD(a0, b0, acc_e0, acc_e1);

        // clear any hi-bits now that we have the multiplier
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

        // now add that into the lo part again.  Here we add in 
        // the two words of the hi-mul result.
        b0 = _mm512_load_epi64(c->data + 0 * VECLEN);
        b1 = _mm512_load_epi64(c->data + 1 * VECLEN);
        b2 = zero;
        b3 = zero;
        acc_e0 = _mm512_add_epi64(acc_e0, b0);
        VEC_CARRYPROP_LOHI(acc_e0, b2);         // lo word and carry for 1st add.
        acc_e1 = _mm512_add_epi64(acc_e1, b1);
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // lo word and carry for 2nd add.
        acc_e1 = _mm512_add_epi64(acc_e1, b2);  // add 1st carry to 2nd lo word.
        VEC_CARRYPROP_LOHI(acc_e1, b3);         // that could have a carry, add it to b3.

        // save first two result words 
        _mm512_store_epi64(c->data + 0 * VECLEN, acc_e0);
        _mm512_store_epi64(c->data + 1 * VECLEN, acc_e1);

        // propagate the carry through the rest of the lo result.
        b2 = _mm512_load_epi64(c->data + 2 * VECLEN);
        scarry = 0;
        a0 = _mm512_adc_epi52(b2, scarry, b3, &scarry);
        _mm512_store_epi64(c->data + 2 * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i = 3; (i < wshift) && (scarry > 0); i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "lo + hi * c:");
#endif

        // it is unlikely, but possible that we still have a final carry to propagate
        if (scarry > 0)
        {
            printf("unlikely add\n");
            b0 = _mm512_set1_epi64(mdata->isMersenne);
            a0 = _mm512_load_epi64(c->data + 0 * VECLEN);
            a0 = _mm512_adc_epi52(a0, 0, b0, &scarry);
            _mm512_store_epi64(c->data + 0 * VECLEN, a0);

            i = 1;
            while (scarry > 0)
            {
                a1 = _mm512_load_epi64(c->data + i * VECLEN);

                //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
                __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
                scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
                a0 = _mm512_and_epi64(t, vlmask);

                _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
                i++;
            }
        }

        // finally clear any hi-bits in the result
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

    }
    else if (mdata->isMersenne > 0)
    {
        // now add the low part into the high.
        scarry = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);
            b0 = _mm512_load_epi64(s->data + i * VECLEN);

            //a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
            __m512i t = _mm512_add_epi64(a1, b0);
            t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpgt_epu64_mask(t, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_adc_epi52(a1, scarry, b0, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add it back in.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // clear the potential hi-bit
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));

        for (i = wshift + 1; i <= NWORDS; i++)
        {
            _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif
    }
    else
    {
        // now subtract the hi part from the lo
        scarry = 0;
        __mmask8 carryout = 0;
        for (i = 0; i < wshift; i++)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);   // hi
            b0 = _mm512_load_epi64(s->data + i * VECLEN);   // lo

            //a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
            __m512i t = _mm512_sub_epi64(b0, a1);
            carryout = _mm512_cmpgt_epu64_mask(a1, b0);
            __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
            a0 = _mm512_and_epi64(t2, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
        }

        a1 = _mm512_load_epi64(c->data + i * VECLEN);
        b0 = _mm512_load_epi64(s->data + i * VECLEN);
        b0 = _mm512_and_epi64(vbshift, b0);
        a0 = _mm512_sbb_epi52(b0, scarry, a1, &scarry);
        _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));

        for (i++; i < NWORDS; i++)
        {
            _mm512_store_epi64(s->data + i * VECLEN, _mm512_set1_epi64(0));
        }

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after add:");
#endif

        // if there was a carry, add 1.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        scarry = _mm512_test_epi64_mask(a1, vbpshift);
        i = 0;
        while (scarry > 0)
        {
            a1 = _mm512_load_epi64(c->data + i * VECLEN);

            //a0 = _mm512_addcarry_epi52(a1, scarry, &scarry);
            __m512i t = _mm512_add_epi64(a1, _mm512_maskz_set1_epi64(scarry, 1));
            scarry = _mm512_cmpeq_epu64_mask(a1, vlmask);
            a0 = _mm512_and_epi64(t, vlmask);

            _mm512_store_epi64(c->data + i * VECLEN, _mm512_and_epi64(vlmask, a0));
            i++;
        }

        // and add one to the hi-bit.  This should resolve the borrow bit.
        a1 = _mm512_load_epi64(c->data + wshift * VECLEN);
        //a1 = _mm512_add_epi64(_mm512_set1_epi64((1ULL << (uint64_t)(bshift))), a1);
        //_mm512_store_epi64(c->data + wshift * VECLEN,
        //    _mm512_and_epi64(_mm512_set1_epi64(0xfffffffffffffULL), a1));

#ifdef DEBUG_MERSENNE
        print_vechexbignum(c, "after carry add:");
        exit(1);
#endif

        _mm512_store_epi64(c->data + wshift * VECLEN,
            _mm512_and_epi64(vbshift, a1));
    }

    c->size = NWORDS;
    return;
}

void vecaddmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }

    // add
    for (i = 0; i < NWORDS; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);

        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_add_epi64(avec, bvec);
        t = _mm512_add_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_cmpgt_epu64_mask(t, lomask);
        cvec = _mm512_and_epi64(t, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    // check for a carry.
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, vbpshift);

    if (mdata->isMersenne > 0)
    {
        // the modulo is just the low part plus 1 (the carry, if present).
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);

        bvec = _mm512_set1_epi64(mdata->isMersenne);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // the modulo is just the low part minus 1 (the carry, if present).
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_set1_epi64(-mdata->isMersenne);
        bvec = _mm512_mask_sbb_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    _mm512_store_epi64(c->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vecsubmod52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* c, vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // s1 is of length VECLEN
    // a, b, c, n, and s1 are aligned
    // a and b are both positive
    // a >= b
    // n is the montgomery base
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 mask = 0;
    __mmask8 mask2 = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
    }

    // subtract
    carry = 0;
    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        
        //cvec = _mm512_sbb_epi52(avec, carry, bvec, &carry);
        __m512i t = _mm512_sub_epi64(avec, bvec);
        carryout = _mm512_cmpgt_epu64_mask(bvec, avec);
        __m512i t2 = _mm512_sub_epi64(t, _mm512_maskz_set1_epi64(carry, 1));
        carry = _mm512_kor(carryout, _mm512_cmpgt_epu64_mask(t2, t));
        cvec = _mm512_and_epi64(t2, lomask);

        _mm512_store_epi64(c->data + i * VECLEN, cvec);
    }

    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        nvec = lomask;
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, mask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        mask = carry;
        carry = 0;
        nvec = _mm512_set1_epi64(1);
        cvec = _mm512_load_epi64(c->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, mask, carry, nvec, &carry);
        _mm512_store_epi64(c->data + 0 * VECLEN, bvec);
        for (i = 1; (i <= wshift) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(c->data + i * VECLEN);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(c->data + i * VECLEN, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(c->data + i * VECLEN, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(c->data + wshift * VECLEN);
    _mm512_store_epi64(c->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    return;
}

void vec_simul_addsub52_mersenne(vec_bignum_t* a, vec_bignum_t* b, vec_bignum_t* sum, vec_bignum_t* diff,
    vec_monty_t* mdata)
{
    // assumptions:
    // a, b, c are of length VECLEN * NWORDS
    // a, b, c, and n are aligned
    // a and b are both positive
    // n is the montgomery base
    // produce sum = a + b and diff = a - b at the same time which
    // saves 3N loads (only have to load a,b, and n once)
    int i;
    uint32_t NWORDS = a->WORDS_ALLOC;
    __mmask8 carry = 0;
    __mmask8 carryout = 0;
    __mmask8 borrow = 0;
    __mmask8 bmask = 0;
    __m512i avec;
    __m512i bvec;
    __m512i cvec;
    __m512i nvec;
    __m512i lomask = _mm512_set1_epi64(0xfffffffffffffULL);
    __m512i zero = _mm512_set1_epi64(0);

    __m512i vbshift, vbpshift;
    int bshift = mdata->nbits % 52;
    int wshift = mdata->nbits / 52;

    if (mdata->use_vnbits)
    {
        int bshift_a[VECLEN];
        for (i = 0; i < VECLEN; i++)
        {
            bshift_a[i] = mdata->vnbits[i] % 52;
        }

        //vbshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
        //    (1ULL << (uint64_t)(bshift_a[0])) - 1ULL);
        //
        //vbpshift = _mm512_set_epi64(
        //    (1ULL << (uint64_t)(bshift_a[7])),
        //    (1ULL << (uint64_t)(bshift_a[6])),
        //    (1ULL << (uint64_t)(bshift_a[5])),
        //    (1ULL << (uint64_t)(bshift_a[4])),
        //    (1ULL << (uint64_t)(bshift_a[3])),
        //    (1ULL << (uint64_t)(bshift_a[2])),
        //    (1ULL << (uint64_t)(bshift_a[1])),
        //    (1ULL << (uint64_t)(bshift_a[0])));

        vbshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[1])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[2])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[3])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[4])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[5])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[6])) - 1ULL,
            (1ULL << (uint64_t)(bshift_a[7])) - 1ULL);

        vbpshift = _mm512_set_epi64(
            (1ULL << (uint64_t)(bshift_a[0])),
            (1ULL << (uint64_t)(bshift_a[1])),
            (1ULL << (uint64_t)(bshift_a[2])),
            (1ULL << (uint64_t)(bshift_a[3])),
            (1ULL << (uint64_t)(bshift_a[4])),
            (1ULL << (uint64_t)(bshift_a[5])),
            (1ULL << (uint64_t)(bshift_a[6])),
            (1ULL << (uint64_t)(bshift_a[7])));
    }
    else
    {
        vbshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)) - 1ULL);
        vbpshift = _mm512_set1_epi64((1ULL << (uint64_t)(bshift)));
    }


    //for (i = 0; i < NWORDS; i++)
    for (i = 0; i <= wshift; i++)
    {
        // add/sub
        avec = _mm512_load_epi64(a->data + i * VECLEN);
        bvec = _mm512_load_epi64(b->data + i * VECLEN);
        
        //cvec = _mm512_adc_epi52(avec, carry, bvec, &carry);
        //cvec = _mm512_sbb_epi52(avec, borrow, bvec, &borrow);
        _mm512_adcsbb_epi52(avec, carry, borrow, bvec, lomask, avec, bvec);
        
        _mm512_store_epi64(sum->data + i * VECLEN, avec);
        _mm512_store_epi64(diff->data + i * VECLEN, bvec);
    }

    bmask = borrow;     // result too small, need to add n

    // check for a carry.
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    carry = _mm512_test_epi64_mask(avec, vbpshift);

    // deal with the sum:
    if (mdata->isMersenne > 0)
    {
        // sum produced a carry, subtract n, which is all FF's, so
        // equivalent to adding 1.
        // the modulo is just the low part plus 1 (the carry, if present).
        cvec = _mm512_load_epi64(sum->data + 0 * VECLEN);
        //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
        bvec = _mm512_set1_epi64(mdata->isMersenne);
        bvec = _mm512_mask_adc_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(sum->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(sum->data + i * VECLEN);
            
            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // the modulo is just the low part minus 1 (the carry, if present).
        cvec = _mm512_load_epi64(sum->data + 0 * VECLEN);
        //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
        //mm512_mask_sbb_epi52(__m512i a, __mmask8 m, __mmask8 c, __m512i b, __mmask8* cout)
        bvec = _mm512_set1_epi64(-mdata->isMersenne);
        bvec = _mm512_mask_sbb_epi52(cvec, carry, 0, bvec, &carry);
        _mm512_store_epi64(sum->data + 0 * VECLEN, bvec);

        for (i = 1; (i < NWORDS) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(sum->data + i * VECLEN);
            
            //bvec = _mm512_subborrow_epi52(cvec, carry, &carry);
            __m512i t = _mm512_sub_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, zero);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(sum->data + i * VECLEN, bvec);
        }
    }

    // clear the potential hi-bit
    avec = _mm512_load_epi64(sum->data + wshift * VECLEN);
    _mm512_store_epi64(sum->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));

    // deal with the diff:
    if (mdata->isMersenne > 0)
    {
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        nvec = _mm512_set1_epi64(0x10000000000000ULL);
        nvec = _mm512_sub_epi64(nvec, _mm512_set1_epi64(mdata->isMersenne));
        cvec = _mm512_load_epi64(diff->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, 0, nvec, &carry);
        _mm512_store_epi64(diff->data + 0 * VECLEN, bvec);
        nvec = _mm512_set1_epi64(0xfffffffffffffULL);
        for (i = 1; i <= wshift; i++)
        {
            cvec = _mm512_load_epi64(diff->data + i * VECLEN);

            //bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
            __m512i t = _mm512_add_epi64(cvec, nvec);
            t = _mm512_mask_add_epi64(cvec, bmask, t, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpgt_epu64_mask(t, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN, bvec);
        }
    }
    else
    {
        // diff produced a carry, add n, which is all 1000...0001, so
        // equivalent to adding 1.
        // if we had a final carry, then b was bigger than a so we need to re-add n.
        carry = 0;
        nvec = _mm512_set1_epi64(1);
        cvec = _mm512_load_epi64(diff->data + 0 * VECLEN);
        bvec = _mm512_mask_adc_epi52(cvec, bmask, carry, nvec, &carry);
        _mm512_store_epi64(diff->data + 0 * VECLEN, bvec);
        for (i = 1; (i <= wshift) && (carry > 0); i++)
        {
            cvec = _mm512_load_epi64(diff->data + i * VECLEN);

            //bvec = _mm512_addcarry_epi52(cvec, carry, &carry);
            __m512i t = _mm512_add_epi64(cvec, _mm512_maskz_set1_epi64(carry, 1));
            carry = _mm512_cmpeq_epu64_mask(cvec, lomask);
            bvec = _mm512_and_epi64(t, lomask);

            _mm512_store_epi64(diff->data + i * VECLEN, bvec);
        }
    }

    //for (; i < NWORDS; i++)
    //{
    //    _mm512_store_epi64(diff->data + i * VECLEN, _mm512_set1_epi64(0));
    //}

    // clear the potential hi-bit
    avec = _mm512_load_epi64(diff->data + wshift * VECLEN);
    _mm512_store_epi64(diff->data + wshift * VECLEN,
        _mm512_and_epi64(vbshift, avec));


    return;
}


#endif