/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Quantized Local Frequency Coding functions                */
/*-----------------------------------------------------------*/

/*--

This file is a part of bsc and/or libbsc, a program and a library for
lossless, block-sorting data compression.

   Copyright (c) 2009-2021 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information and file AUTHORS
for full list of contributors.

See also the bsc and libbsc web site:
  http://libbsc.com/ for more information.

--*/

#include "fc.h"
#include "fc_model.h"
#include "predictor.h"
#include "tables.h"
#include "rangecoder.h"
#include <memory>
#include <stddef.h>
#include <stdio.h>

// Squeeze 2 bits from contextRank4 for avgRank bucket
// Map avgRank to 4 coarse-grained buckets (0..3)
static INLINE int fc_avgRank_bucket(int avgRank)
{
    // avgRank ~ exponential moving average of MTF rank, range approximately 0..maxRank
    if (avgRank < 2)  return 0;  // extremely local: almost always at tid=0/1
    if (avgRank < 4)  return 1;  // somewhat local: mainly small ranks
    if (avgRank < 8)  return 2;  // medium
    return 3;                    // rank often large
}

// #define FC_DEBUG

#if FC_CPU_TYPES >= FC_CPU_TYPES_SSE41
static const __m128i ALIGNED(64) tid_reorder[16] =
{
    _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 0, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 0, 9, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 11, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 13, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 14, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 15),
    _mm_setr_epi8(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0),
};
#elif FC_CPU_TYPES == FC_CPU_TYPES_A64
static const unsigned char ALIGNED(64) tid_reorder[16][16] =
{
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 0, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 0, 8, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 0, 9, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 11, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 12, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 13, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 14, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 15},
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0},
};
#endif

#if FC_CPU_TYPES >= FC_CPU_TYPES_SSE2

INLINE ptrdiff_t transform_ (const unsigned char * RESTRICT input, ptrdiff_t i, unsigned char currentChar)
{
#if FC_CPU_TYPES >= FC_CPU_TYPES_AVX2
    __m256i v = _mm256_set1_epi8(currentChar);

    while (i >= 32)
    {
        i -= 32; int m = _mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i *)(input + i)), v));
        if (m != (int)0xffffffff) { return i + fc_bit_scan_reverse(((unsigned int)(~m))); }
    }
#elif FC_CPU_TYPES >= FC_CPU_TYPES_SSE2
    __m128i v = _mm_set1_epi8(currentChar);

    while (i >= 16)
    {
        i -= 16; int m = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_loadu_si128((const __m128i *)(input + i)), v));
        if (m != 0xffff) { return i + fc_bit_scan_reverse((unsigned int)(m ^ 0xffff)); }
    }
#endif

    do {} while ((--i >= 0) && (input[i] == currentChar)); return i;
}

unsigned char * fctransform (const unsigned char * RESTRICT input, unsigned char * RESTRICT buffer, int n, unsigned char * RESTRICT MTFTable)
{
    signed char ALIGNED(64) tids[CHAR_SIZE];
    signed char ALIGNED(64) flags[CHAR_SIZE];

    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { tids[i] = (signed char)(i - 128); }
    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { flags[i] = 0; }

    ptrdiff_t i = (ptrdiff_t)n - 1, j = n; signed char nSymbols = 0;

    for (; i >= 0; )
    {
        unsigned char currentChar1 = input[i]; i = transform_(input, i, currentChar1); if (i < 0) { i = 0; break; }
        unsigned char currentChar2 = input[i]; i = transform_(input, i, currentChar2);

        signed char tid1 = tids[currentChar1], tid2 = tids[currentChar2]; tid2 += tid1 > tid2;

        buffer[--j] = tid1 + 128; if (flags[currentChar1] == 0) { flags[currentChar1] = 1; buffer[j] = nSymbols++; }
        buffer[--j] = tid2 + 128; if (flags[currentChar2] == 0) { flags[currentChar2] = 1; buffer[j] = nSymbols++; }

        for (int t = 0 * 32; t < 1 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 1 * 32; t < 2 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 2 * 32; t < 3 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 3 * 32; t < 4 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 4 * 32; t < 5 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 5 * 32; t < 6 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 6 * 32; t < 7 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }
        for (int t = 7 * 32; t < 8 * 32; ++t) { tids[t] -= (tid1 > tids[t] ? (signed char)-1 : (signed char)0) + (tid2 > tids[t] ? (signed char)-1 : (signed char)0); }

        tids[currentChar1] = -127; tids[currentChar2] = -128;
    }

    if (i >= 0)
    {
        unsigned char currentChar = input[0]; signed char tid = tids[currentChar];

        buffer[--j] = tid + 128; if (flags[currentChar] == 0) { flags[currentChar] = 1; buffer[j] = nSymbols++; }

        for (int t = 0 * 32; t < 1 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 1 * 32; t < 2 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 2 * 32; t < 3 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 3 * 32; t < 4 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 4 * 32; t < 5 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 5 * 32; t < 6 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 6 * 32; t < 7 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        for (int t = 7 * 32; t < 8 * 32; ++t) { tids[t] -= (tids[t] < tid ? -1 : 0); }
        tids[currentChar] = -128;
    }

    buffer[n - 1] = 1;

    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { MTFTable[tids[i] + 128] = (unsigned char)i; }
    for (ptrdiff_t i = 1; i < CHAR_SIZE; ++i) { if (flags[MTFTable[i]] == 0) { MTFTable[i] = MTFTable[i - 1]; break; } }

    return buffer + j;
}

#elif FC_CPU_TYPES == FC_CPU_TYPES_A64

INLINE ptrdiff_t transform_ (const unsigned char * RESTRICT input, ptrdiff_t i, unsigned long long currentChar)
{
    unsigned long long v = currentChar; v |= (v << 8); v |= (v << 16); v |= (v << 32);

    while (i >= 8)
    {
        i -= 8; unsigned long long m = (*(unsigned long long const *)(input + i)) ^ v;
        if (m != 0) { return i + (fc_bit_scan_reverse64(m) / 8); }
    }

    do {} while ((--i >= 0) && (input[i] == currentChar)); return i;
}

unsigned char * fctransform (const unsigned char * RESTRICT input, unsigned char * RESTRICT buffer, int n, unsigned char * RESTRICT MTFTable)
{
    signed char ALIGNED(64) tids[CHAR_SIZE];
    signed char ALIGNED(64) flags[CHAR_SIZE];

    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { tids[i] = (signed char)(i - 128); }
    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { flags[i] = 0; }

    ptrdiff_t i = (ptrdiff_t)n - 1, j = n; signed char nSymbols = 0;

    for (; i >= 0;)
    {
        unsigned char currentChar1 = input[i]; i = transform_(input, i, currentChar1); if (i < 0) { i = 0; break; }
        unsigned char currentChar2 = input[i]; i = transform_(input, i, currentChar2);

        signed char tid1 = tids[currentChar1], tid2 = tids[currentChar2]; tid2 += tid1 > tid2;

        buffer[--j] = tid1 + 128; if (flags[currentChar1] == 0) { flags[currentChar1] = 1; buffer[j] = nSymbols++; }
        buffer[--j] = tid2 + 128; if (flags[currentChar2] == 0) { flags[currentChar2] = 1; buffer[j] = nSymbols++; }

        int8x16_t r1 = vdupq_n_s8(tid1), r2 = vdupq_n_s8(tid2), x, y;

        x = vld1q_s8((int8_t const *)(tids + 16 * 0)); y = vld1q_s8((int8_t const *)(tids + 16 * 1));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 0), x); vst1q_s8((int8_t *)(tids + 16 * 1), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 2)); y = vld1q_s8((int8_t const *)(tids + 16 * 3));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 2), x); vst1q_s8((int8_t *)(tids + 16 * 3), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 4)); y = vld1q_s8((int8_t const *)(tids + 16 * 5));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 4), x); vst1q_s8((int8_t *)(tids + 16 * 5), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 6)); y = vld1q_s8((int8_t const *)(tids + 16 * 7));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 6), x); vst1q_s8((int8_t *)(tids + 16 * 7), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 8)); y = vld1q_s8((int8_t const *)(tids + 16 * 9));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 8), x); vst1q_s8((int8_t *)(tids + 16 * 9), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 10)); y = vld1q_s8((int8_t const *)(tids + 16 * 11));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 10), x); vst1q_s8((int8_t *)(tids + 16 * 11), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 12)); y = vld1q_s8((int8_t const *)(tids + 16 * 13));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 12), x); vst1q_s8((int8_t *)(tids + 16 * 13), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 14)); y = vld1q_s8((int8_t const *)(tids + 16 * 15));
        x = vsubq_s8(vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r1, x))), vreinterpretq_s8_u8(vcgtq_s8(r2, x)));
        y = vsubq_s8(vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r1, y))), vreinterpretq_s8_u8(vcgtq_s8(r2, y)));
        vst1q_s8((int8_t *)(tids + 16 * 14), x); vst1q_s8((int8_t *)(tids + 16 * 15), y);

        tids[currentChar1] = -127; tids[currentChar2] = -128;
    }

    if (i >= 0)
    {
        unsigned char currentChar = input[0]; signed char tid = tids[currentChar];

        buffer[--j] = tid + 128; if (flags[currentChar] == 0) { flags[currentChar] = 1; buffer[j] = nSymbols++; }

        int8x16_t r = vdupq_n_s8(tid), x, y;

        x = vld1q_s8((int8_t const *)(tids + 16 * 0)); y = vld1q_s8((int8_t const *)(tids + 16 * 1));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 0), x); vst1q_s8((int8_t *)(tids + 16 * 1), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 2)); y = vld1q_s8((int8_t const *)(tids + 16 * 3));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 2), x); vst1q_s8((int8_t *)(tids + 16 * 3), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 4)); y = vld1q_s8((int8_t const *)(tids + 16 * 5));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 4), x); vst1q_s8((int8_t *)(tids + 16 * 5), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 6)); y = vld1q_s8((int8_t const *)(tids + 16 * 7));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 6), x); vst1q_s8((int8_t *)(tids + 16 * 7), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 8)); y = vld1q_s8((int8_t const *)(tids + 16 * 9));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 8), x); vst1q_s8((int8_t *)(tids + 16 * 9), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 10)); y = vld1q_s8((int8_t const *)(tids + 16 * 11));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 10), x); vst1q_s8((int8_t *)(tids + 16 * 11), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 12)); y = vld1q_s8((int8_t const *)(tids + 16 * 13));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 12), x); vst1q_s8((int8_t *)(tids + 16 * 13), y);

        x = vld1q_s8((int8_t const *)(tids + 16 * 14)); y = vld1q_s8((int8_t const *)(tids + 16 * 15));
        x = vsubq_s8(x, vreinterpretq_s8_u8(vcgtq_s8(r, x)));
        y = vsubq_s8(y, vreinterpretq_s8_u8(vcgtq_s8(r, y)));
        vst1q_s8((int8_t *)(tids + 16 * 14), x); vst1q_s8((int8_t *)(tids + 16 * 15), y);

        tids[currentChar] = -128;
    }

    buffer[n - 1] = 1;

    for (ptrdiff_t i = 0; i < CHAR_SIZE; ++i) { MTFTable[tids[i] + 128] = (unsigned char)i; }
    for (ptrdiff_t i = 1; i < CHAR_SIZE; ++i) { if (flags[MTFTable[i]] == 0) { MTFTable[i] = MTFTable[i - 1]; break; } }

    return buffer + j;
}

#else

unsigned char * fctransform (const unsigned char * RESTRICT input, unsigned char * RESTRICT buffer, int n, unsigned char * RESTRICT MTFTable)
{
    unsigned char Flag[CHAR_SIZE];

    for (int i = 0; i < CHAR_SIZE; ++i) Flag[i] = 0;
    for (int i = 0; i < CHAR_SIZE; ++i) MTFTable[i] = i;

    if (input[n - 1] == 0)
    {
        MTFTable[0] = 1; MTFTable[1] = 0;
    }

    int index = n, nSymbols = 0;
    for (int i = n - 1; i >= 0;)
    {
        unsigned char currentChar = input[i--];
        for (; (i >= 0) && (input[i] == currentChar); --i) ;

        unsigned char previousChar = MTFTable[0], tid = 1; MTFTable[0] = currentChar;
        while (true)
        {
            unsigned char temporaryChar0 = MTFTable[tid + 0]; MTFTable[tid + 0] = previousChar;
            if (temporaryChar0 == currentChar) { tid += 0; break; }

            unsigned char temporaryChar1 = MTFTable[tid + 1]; MTFTable[tid + 1] = temporaryChar0;
            if (temporaryChar1 == currentChar) { tid += 1; break; }

            unsigned char temporaryChar2 = MTFTable[tid + 2]; MTFTable[tid + 2] = temporaryChar1;
            if (temporaryChar2 == currentChar) { tid += 2; break; }

            unsigned char temporaryChar3 = MTFTable[tid + 3]; MTFTable[tid + 3] = temporaryChar2;
            if (temporaryChar3 == currentChar) { tid += 3; break; }

            tid += 4; previousChar = temporaryChar3;
        }

        if (Flag[currentChar] == 0)
        {
            Flag[currentChar] = 1;
            tid = nSymbols++;
        }

        buffer[--index] = tid;
    }

    buffer[n - 1] = 1;

    for (int tid = 1; tid < CHAR_SIZE; ++tid)
    {
        if (Flag[MTFTable[tid]] == 0)
        {
            MTFTable[tid] = MTFTable[tid - 1];
            break;
        }
    }

    return buffer + index;
}

#endif


int fc_encode_do (const unsigned char * input, unsigned char * output, unsigned char * buffer, int inputSize, int outputSize, fcModel1 * model)
{
    unsigned char MTFTable[CHAR_SIZE];

    fcinit_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char tidHistory[CHAR_SIZE], runHistory[CHAR_SIZE];
    for (int i = 0; i < CHAR_SIZE; ++i)
    {
        tidHistory[i] = runHistory[i] = 0;
    }

    // Global frequency sorting optimizes the symbol distribution of MTF output
    unsigned char * tidArray = fctransform(input, buffer, inputSize, MTFTable);

#ifdef FC_DEBUG
    fprintf(stderr, "\ninput after bwt: \n");
    for (int i = 0; i < 20 && i < inputSize; i++)fprintf(stderr, "%u|", input[i]);
    fprintf(stderr, "\n");
#endif

    RangeCoder coder;

    coder.InitEncoder(output, outputSize);
    coder.EncodeWord((unsigned int)inputSize);

    unsigned char usedChar[CHAR_SIZE];
    for (int i = 0; i < CHAR_SIZE; ++i) usedChar[i] = 0;

#ifdef FC_DEBUG
    fprintf(stderr, "\nMTFTable: \n");
#endif

    int prevChar = -1;
    for (int tid = 0; tid < CHAR_SIZE; ++tid)
    {
        int currentChar = MTFTable[tid];
#ifdef FC_DEBUG
        fprintf(stderr, "%u|", currentChar); if (((tid + 1) & 0xF) == 0)fprintf(stderr, "\n");
#endif

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < CHAR_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if ((currentChar >> (bit + 1)) == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                coder.EncodeBit(currentChar & (1 << bit));
            }
        }

        if (currentChar == prevChar)
        {
            maxRank = fc_bit_scan_reverse(tid - 1); // Get the number of valid bits for the last tid, i.e., the distance between 0 and the leftmost 1, e.g., 5 for 101000
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }
#ifdef FC_DEBUG
    fprintf(stderr, " @ maxRank %d\n", maxRank);
#endif

    const unsigned char * inputEnd      = input  + inputSize;
    const unsigned char * tidArrayEnd  = buffer + inputSize;

    for (; tidArray < tidArrayEnd; ) // Process each tid in loop, i.e., the content of transform
    {
        if (coder.is_end())
        {
            return FC_FAILED_ZIP;
        }

        // Calculate how many times the current character repeats, recorded as runSize, and input advances to the next non-repeating character
        int currentChar = *input, runSize;
        {
            const unsigned char * inputStart = input++;

            if (tidArray >= tidArrayEnd - 16)
            {
                while ((input < inputEnd) && (*input == currentChar)) { input++; }
            }
            else
            {
#if FC_CPU_TYPES >= FC_CPU_TYPES_SSE2
                __m128i v = _mm_set1_epi8(currentChar);

                while (true)
                {
                   int m = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_loadu_si128((const __m128i *)input), v));
                   if (m != 0xffff)
                   {
                      input += fc_bit_scan_forward((unsigned int)(~m));
                      break;
                   }

                   input += 16;
                }
#elif FC_CPU_TYPES == FC_CPU_TYPES_A64
                unsigned long long v = currentChar; v |= (v << 8); v |= (v << 16); v |= (v << 32);

                while (true)
                {
                    unsigned long long m = (*(unsigned long long const *)input) ^ v;
                    if (m != 0)
                    {
                        input += fc_bit_scan_forward64(m) / 8;
                        break;
                    }

                    input += 8;
                }
#else
                while (*input == currentChar) { input++; }
#endif
            }

            runSize = (int)(input - inputStart);
        }

        // Modeling approach: First describe the variables involved in the current process: tid, current character (MTFTable), runSize (character repetition count)
        int tid = *tidArray++;
        int history = tidHistory[currentChar]; // This history value is the number of valid bits for the tid corresponding to the current character minus 1, tidHistory should be less than 8, i.e., only 3 bits
        int state = model_tid_state(contextRank4, contextRun, history);
        // State transition table, contextRanke :3 bits, contextRun : 4 bits, rank : clipped to [0,7], 3 bits, runsizeHistory : clipped to [0,7], 3 bits
        // Total: 3+4+3+3=13 bits 8192 entries
        // contextRank4 is the history record of the lowest 2 bits of tid (i.e., save the last 4 records at most), when tid > 3 take 3, otherwise take the entire tid
        // contextRun is the history record of whether runSize is less than 3 (i.e., save the last 4 records at most)

        short *            RESTRICT statePredictor  = & model->tid_t.state_model[state]; //
        short *            RESTRICT charPredictor   = & model->tid_t.char_model[currentChar];
        short *            RESTRICT staticPredictor = & model->tid_t.static_model;
        pMixer * RESTRICT mixer           = & model->tid_mixer[currentChar];

        // @@@ --- In summary: encode condition flag bits for decompression restoration; encode bit width of tid (corresponding to model bit_width); encode real bits of tid (corresponding to model bits_value)
        if (avgRank < 32)
        {
            if (tid == 1) // @@@ --- step : tid of 1 does not need encoding, only need to encode 0, similar to encoding a condition, this encoding is completely for decompression restoration, the meaning of tid being 1 is not fully understood, but the approach is understandable
            {
                tidHistory[currentChar] = 0;

                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                fc_count::UpdateBit0(*statePredictor,  FC_RANK_TS_TH0, FC_RANK_TS_AR0);
                fc_count::UpdateBit0(*charPredictor,   FC_RANK_TC_TH0, FC_RANK_TC_AR0);
                fc_count::UpdateBit0(*staticPredictor, FC_RANK_TP_TH0, FC_RANK_TP_AR0);

                coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RANK_TM_LR0, FC_RANK_TM_LR1, FC_RANK_TM_LR2, FC_RANK_TM_TH0, FC_RANK_TM_AR0));
            }
            else
            {
                { // @@@ --- step : current tid needs encoding, first encode 1, which also represents a condition, this encoding is completely for decompression restoration
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    fc_count::UpdateBit1(*statePredictor,  FC_RANK_TS_TH1, FC_RANK_TS_AR1);
                    fc_count::UpdateBit1(*charPredictor,   FC_RANK_TC_TH1, FC_RANK_TC_AR1);
                    fc_count::UpdateBit1(*staticPredictor, FC_RANK_TP_TH1, FC_RANK_TP_AR1);

                    coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RANK_TM_LR0, FC_RANK_TM_LR1, FC_RANK_TM_LR2, FC_RANK_TM_TH1, FC_RANK_TM_AR1));
                }

                // __builtin_clz (x) ^ 31: __builtin_clz (x) represents the number of zeros before the first 1 on the left of x, e.g., when x is 40, it is 26, 31 - 26 = 5
                int bitRankSize = fc_bit_scan_reverse(tid); tidHistory[currentChar] = bitRankSize;
                // printf("tid %d, bitRankSize %d\n", tid, bitRankSize);

                statePredictor  = & model->tid_t.bit_width.state_model[state][0];
                charPredictor   = & model->tid_t.bit_width.char_model[currentChar][0];
                staticPredictor = & model->tid_t.bit_width.static_model[0];
                mixer           = & model->tid_mixerbit_width[history < 1 ? 1 : history][1];//history is the number of valid bits for the current character's last time minus 1

                // @@@ --- step : encode the valid bits of tid value, the highest bit is not encoded during encoding because the highest bit is definitely 1
                // For example, if tid is 40, then bitRankSize is 5, note that tid's bits are recorded in Rank.bit_width
                // bit range is [1, 4], why not 5 is because the highest bit is definitely 1
                for (int bit = 1; bit < bitRankSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
                {
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    fc_count::UpdateBit1(*statePredictor,  FC_RANK_ES_TH1, FC_RANK_ES_AR1);
                    fc_count::UpdateBit1(*charPredictor,   FC_RANK_EC_TH1, FC_RANK_EC_AR1);
                    fc_count::UpdateBit1(*staticPredictor, FC_RANK_EP_TH1, FC_RANK_EP_AR1);

                    coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RANK_EM_LR0, FC_RANK_EM_LR1, FC_RANK_EM_LR2, FC_RANK_EM_TH1, FC_RANK_EM_AR1));

                    mixer = & model->tid_mixerbit_width[history <= bit ? bit + 1 : history][bit + 1];
                }
                if (bitRankSize < maxRank) // Encode a flag bit when less than maxRank
                {
                    int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                    fc_count::UpdateBit0(*statePredictor,  FC_RANK_ES_TH0, FC_RANK_ES_AR0);
                    fc_count::UpdateBit0(*charPredictor,   FC_RANK_EC_TH0, FC_RANK_EC_AR0);
                    fc_count::UpdateBit0(*staticPredictor, FC_RANK_EP_TH0, FC_RANK_EP_AR0);

                    coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RANK_EM_LR0, FC_RANK_EM_LR1, FC_RANK_EM_LR2, FC_RANK_EM_TH0, FC_RANK_EM_AR0));
                }

                // @@@ --- step : actually encode the tid value, the highest bit is not encoded during encoding because the highest bit is definitely 1
                // Encode the mantissa, for example if bank is 40, i.e., 0010 1000, then bitRankSize==31-26=5, bit range is [4, 0]
                statePredictor  = & model->tid_t.bits_value[bitRankSize].state_model[state][0];
                charPredictor   = & model->tid_t.bits_value[bitRankSize].char_model[currentChar][0];
                staticPredictor = & model->tid_t.bits_value[bitRankSize].static_model[0];
                mixer           = & model->tid_mixerbits_value[bitRankSize];

                // bit range is [4, 0], why not 5 is because the highest bit is definitely 1
                for (int context = 1, bit = bitRankSize - 1; bit >= 0; --bit)  // context is actually the value of tid accurate to bit
                {
                    if (tid & (1 << bit))
                    {
                        int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                        fc_count::UpdateBit1(statePredictor[context],  FC_RANK_MS_TH1, FC_RANK_MS_AR1);
                        fc_count::UpdateBit1(charPredictor[context],   FC_RANK_MC_TH1, FC_RANK_MC_AR1);
                        fc_count::UpdateBit1(staticPredictor[context], FC_RANK_MP_TH1, FC_RANK_MP_AR1);

                        coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RANK_MM_LR0, FC_RANK_MM_LR1, FC_RANK_MM_LR2, FC_RANK_MM_TH1, FC_RANK_MM_AR1));

                        context += context + 1;
                    }
                    else
                    {
                        int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                        fc_count::UpdateBit0(statePredictor[context],  FC_RANK_MS_TH0, FC_RANK_MS_AR0);
                        fc_count::UpdateBit0(charPredictor[context],   FC_RANK_MC_TH0, FC_RANK_MC_AR0);
                        fc_count::UpdateBit0(staticPredictor[context], FC_RANK_MP_TH0, FC_RANK_MP_AR0);

                        coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RANK_MM_LR0, FC_RANK_MM_LR1, FC_RANK_MM_LR2, FC_RANK_MM_TH0, FC_RANK_MM_AR0));

                        context += context;
                    }
                }
            }
        }
        else // tid is too large, use out_of_range corresponding model for encoding
        {
            // __builtin_clz (x) ^ 31: __builtin_clz (x) represents the number of zeros before the first 1 on the left of x, e.g., when x is 40, it is 26, 31 - 26 = 5
            // Record the bit width of the current character minus 1 into tidHistory
            tidHistory[currentChar] = (unsigned char)fc_bit_scan_reverse(tid);

            statePredictor  = & model->tid_t.out_of_range.state_model[state][0];
            charPredictor   = & model->tid_t.out_of_range.char_model[currentChar][0];
            staticPredictor = & model->tid_t.out_of_range.static_model[0];

            for (int context = 1, bit = maxRank; bit >= 0; --bit) // Align to the maximum bank bits, each tid encodes [0, maxRank] bits
            {
                mixer = & model->tid_mixerout_of_range[context];

                if (tid & (1 << bit))
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    fc_count::UpdateBit1(statePredictor[context],  FC_RANK_PS_TH1, FC_RANK_PS_AR1);
                    fc_count::UpdateBit1(charPredictor[context],   FC_RANK_PC_TH1, FC_RANK_PC_AR1);
                    fc_count::UpdateBit1(staticPredictor[context], FC_RANK_PP_TH1, FC_RANK_PP_AR1);

                    coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RANK_PM_LR0, FC_RANK_PM_LR1, FC_RANK_PM_LR2, FC_RANK_PM_TH1, FC_RANK_PM_AR1));

                    context += context + 1;
                }
                else
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    fc_count::UpdateBit0(statePredictor[context],  FC_RANK_PS_TH0, FC_RANK_PS_AR0);
                    fc_count::UpdateBit0(charPredictor[context],   FC_RANK_PC_TH0, FC_RANK_PC_AR0);
                    fc_count::UpdateBit0(staticPredictor[context], FC_RANK_PP_TH0, FC_RANK_PP_AR0);

                    coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RANK_PM_LR0, FC_RANK_PM_LR1, FC_RANK_PM_LR2, FC_RANK_PM_TH0, FC_RANK_PM_AR0));

                    context += context;
                }
            }
        }

        avgRank         =   (avgRank * 125 + tid * 3) >> 7; // Current tid encoding is completed, put into avgRank, this is actually also a context,
        tid            =   tid - 1;
        history         =   runHistory[currentChar];
        int avgBucket   =   fc_avgRank_bucket(avgRank); // Add avgRank bucket (0..3)
        state           =   lstate(contextRank0, contextRun, tid, history, avgBucket); // lstate internally compresses tid into 0,1,2,>=3 as 2 bits
        statePredictor  = & model->run_t.state_model[state];
        charPredictor   = & model->run_t.char_model[currentChar];
        staticPredictor = & model->run_t.static_model;
        mixer           = & model->run_mixer[currentChar];

        if (runSize == 1)
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;

            int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

            fc_count::UpdateBit0(*statePredictor,  FC_RUN_TS_TH0, FC_RUN_TS_AR0);
            fc_count::UpdateBit0(*charPredictor,   FC_RUN_TC_TH0, FC_RUN_TC_AR0);
            fc_count::UpdateBit0(*staticPredictor, FC_RUN_TP_TH0, FC_RUN_TP_AR0);

            coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RUN_TM_LR0, FC_RUN_TM_LR1, FC_RUN_TM_LR2, FC_RUN_TM_TH0, FC_RUN_TM_AR0));
        }
        else
        {
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                fc_count::UpdateBit1(*statePredictor,  FC_RUN_TS_TH1, FC_RUN_TS_AR1);
                fc_count::UpdateBit1(*charPredictor,   FC_RUN_TC_TH1, FC_RUN_TC_AR1);
                fc_count::UpdateBit1(*staticPredictor, FC_RUN_TP_TH1, FC_RUN_TP_AR1);

                coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RUN_TM_LR0, FC_RUN_TM_LR1, FC_RUN_TM_LR2, FC_RUN_TM_TH1, FC_RUN_TM_AR1));
            }

            int bitRunSize = fc_bit_scan_reverse(runSize); runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->run_t.bit_width.state_model[state][0];
            charPredictor   = & model->run_t.bit_width.char_model[currentChar][0];
            staticPredictor = & model->run_t.bit_width.static_model[0];
            mixer           = & model->run_mixerbit_width[history < 1 ? 1 : history][1];

            for (int bit = 1; bit < bitRunSize; ++bit, ++statePredictor, ++charPredictor, ++staticPredictor)
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                fc_count::UpdateBit1(*statePredictor,  FC_RUN_ES_TH1, FC_RUN_ES_AR1);
                fc_count::UpdateBit1(*charPredictor,   FC_RUN_EC_TH1, FC_RUN_EC_AR1);
                fc_count::UpdateBit1(*staticPredictor, FC_RUN_EP_TH1, FC_RUN_EP_AR1);

                coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RUN_EM_LR0, FC_RUN_EM_LR1, FC_RUN_EM_LR2, FC_RUN_EM_TH1, FC_RUN_EM_AR1));

                mixer = & model->run_mixerbit_width[history <= bit ? bit + 1 : history][bit + 1];
            }
            {
                int probability0 = *charPredictor, probability1 = *statePredictor, probability2 = *staticPredictor;

                fc_count::UpdateBit0(*statePredictor,  FC_RUN_ES_TH0, FC_RUN_ES_AR0);
                fc_count::UpdateBit0(*charPredictor,   FC_RUN_EC_TH0, FC_RUN_EC_AR0);
                fc_count::UpdateBit0(*staticPredictor, FC_RUN_EP_TH0, FC_RUN_EP_AR0);

                coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RUN_EM_LR0, FC_RUN_EM_LR1, FC_RUN_EM_LR2, FC_RUN_EM_TH0, FC_RUN_EM_AR0));
            }

            statePredictor  = & model->run_t.bits_value[bitRunSize].state_model[state][0];
            charPredictor   = & model->run_t.bits_value[bitRunSize].char_model[currentChar][0];
            staticPredictor = & model->run_t.bits_value[bitRunSize].static_model[0];
            mixer           = & model->run_mixerbits_value[bitRunSize];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                if (runSize & (1 << bit))
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    fc_count::UpdateBit1(statePredictor[context],  FC_RUN_MS_TH1, FC_RUN_MS_AR1);
                    fc_count::UpdateBit1(charPredictor[context],   FC_RUN_MC_TH1, FC_RUN_MC_AR1);
                    fc_count::UpdateBit1(staticPredictor[context], FC_RUN_MP_TH1, FC_RUN_MP_AR1);

                    coder.EncodeBit1(mixer->MixupAndUpdateBit1(probability0, probability1, probability2, FC_RUN_MM_LR0, FC_RUN_MM_LR1, FC_RUN_MM_LR2, FC_RUN_MM_TH1, FC_RUN_MM_AR1));

                    if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    int probability0 = charPredictor[context], probability1 = statePredictor[context], probability2 = staticPredictor[context];

                    fc_count::UpdateBit0(statePredictor[context],  FC_RUN_MS_TH0, FC_RUN_MS_AR0);
                    fc_count::UpdateBit0(charPredictor[context],   FC_RUN_MC_TH0, FC_RUN_MC_AR0);
                    fc_count::UpdateBit0(staticPredictor[context], FC_RUN_MP_TH0, FC_RUN_MP_AR0);

                    coder.EncodeBit0(mixer->MixupAndUpdateBit0(probability0, probability1, probability2, FC_RUN_MM_LR0, FC_RUN_MM_LR1, FC_RUN_MM_LR2, FC_RUN_MM_TH0, FC_RUN_MM_AR0));

                    if (bitRunSize <= 5) context += context + 0; else context++;
                }
            }
        }

        contextRank0 = ((contextRank0 << 1) | (tid == 0 ? 1 : 0)) & 0x7;

        // 1) Compress tid history into 3 events * 2bit = 6bit, placed in low 6 bits
        int tidBucket = (tid < 3 ? tid : 3);        // 0,1,2,3
        avgBucket = fc_avgRank_bucket(avgRank); // 0..3
        int rank_hist = ((contextRank4 & 0x3f) << 2) | tidBucket; // Only low 6 bits shifted left and new bucket inserted
        rank_hist &= 0x3f;                                        // Prevent overflow

        // 2) High 2 bits store avgRank bucket
        contextRank4 = rank_hist | (avgBucket << 6);

        // run history remains unchanged
        contextRun = ((contextRun << 1) | (runSize < 3 ? 1 : 0)) & 0xf;
    }

    return coder.FinishEncoder();
}


int fc_decode_do (const unsigned char * input, unsigned char * output, fcModel1 * model)
{
    RangeCoder coder;

    unsigned char ALIGNED(64) MTFTable[CHAR_SIZE];

    fcinit_model(model);

    int contextRank0 = 0;
    int contextRank4 = 0;
    int contextRun   = 0;
    int maxRank      = 7;
    int avgRank      = 0;

    unsigned char tidHistory[CHAR_SIZE], runHistory[CHAR_SIZE];
    for (int i = 0; i < CHAR_SIZE; ++i)
    {
        tidHistory[i] = runHistory[i] = 0;
    }

    coder.InitDecoder(input);
    int n = (int)coder.DecodeWord();

    unsigned char usedChar[CHAR_SIZE];
    for (int i = 0; i < CHAR_SIZE; ++i) usedChar[i] = 0;

    int prevChar = -1;
    for (int tid = 0; tid < CHAR_SIZE; ++tid)
    {
        int currentChar = 0;

        for (int bit = 7; bit >= 0; --bit)
        {
            bool bit0 = false, bit1 = false;

            for (int c = 0; c < CHAR_SIZE; ++c)
            {
                if (c == prevChar || usedChar[c] == 0)
                {
                    if (currentChar == (c >> (bit + 1)))
                    {
                        if (c & (1 << bit)) bit1 = true; else bit0 = true;
                        if (bit0 && bit1) break;
                    }
                }
            }

            if (bit0 && bit1)
            {
                currentChar += currentChar + coder.DecodeBit();
            }
            else
            {
                if (bit0) currentChar += currentChar + 0;
                if (bit1) currentChar += currentChar + 1;
            }
        }

        MTFTable[tid] =  currentChar;

        if (currentChar == prevChar)
        {
            maxRank = fc_bit_scan_reverse(tid - 1);
            break;
        }

        prevChar = currentChar; usedChar[currentChar] = 1;
    }

    for (int i = 0; i < n;)
    {
        // current = (i* 100) / n;
        // fprintf(stderr, "%-s -[%d%%]-\r", "progress", current);
        int                 currentChar     =   MTFTable[0];
        int                 history         =   tidHistory[currentChar];
        int                 state           =   model_tid_state(contextRank4, contextRun, history);

        short *            RESTRICT statePredictor  = & model->tid_t.state_model[state];
        short *            RESTRICT charPredictor   = & model->tid_t.char_model[currentChar];
        short *            RESTRICT staticPredictor = & model->tid_t.static_model;
        pMixer * RESTRICT mixer           = & model->tid_mixer[currentChar];

        int tid = 1;
        if (avgRank < 32)
        {
            if (coder.DecodeBit(mixer->Mixup(*charPredictor, *statePredictor, *staticPredictor)))
            {
                fc_count::UpdateBit1(*statePredictor,  FC_RANK_TS_TH1, FC_RANK_TS_AR1);
                fc_count::UpdateBit1(*charPredictor,   FC_RANK_TC_TH1, FC_RANK_TC_AR1);
                fc_count::UpdateBit1(*staticPredictor, FC_RANK_TP_TH1, FC_RANK_TP_AR1);
                mixer->UpdateBit1(FC_RANK_TM_LR0, FC_RANK_TM_LR1, FC_RANK_TM_LR2, FC_RANK_TM_TH1, FC_RANK_TM_AR1);

                statePredictor  = & model->tid_t.bit_width.state_model[state][0];
                charPredictor   = & model->tid_t.bit_width.char_model[currentChar][0];
                staticPredictor = & model->tid_t.bit_width.static_model[0];
                mixer           = & model->tid_mixerbit_width[history < 1 ? 1 : history][1];

                int bitRankSize = 1;
                while (true)
                {
                    if (bitRankSize == maxRank) break;
                    if (coder.DecodeBit(mixer->Mixup(*charPredictor, *statePredictor, *staticPredictor)))
                    {
                        fc_count::UpdateBit1(*statePredictor,  FC_RANK_ES_TH1, FC_RANK_ES_AR1); statePredictor++;
                        fc_count::UpdateBit1(*charPredictor,   FC_RANK_EC_TH1, FC_RANK_EC_AR1); charPredictor++;
                        fc_count::UpdateBit1(*staticPredictor, FC_RANK_EP_TH1, FC_RANK_EP_AR1); staticPredictor++;
                        mixer->UpdateBit1(FC_RANK_EM_LR0, FC_RANK_EM_LR1, FC_RANK_EM_LR2, FC_RANK_EM_TH1, FC_RANK_EM_AR1);
                        bitRankSize++;
                        mixer = & model->tid_mixerbit_width[history < bitRankSize ? bitRankSize : history][bitRankSize];
                    }
                    else
                    {
                        fc_count::UpdateBit0(*statePredictor,  FC_RANK_ES_TH0, FC_RANK_ES_AR0);
                        fc_count::UpdateBit0(*charPredictor,   FC_RANK_EC_TH0, FC_RANK_EC_AR0);
                        fc_count::UpdateBit0(*staticPredictor, FC_RANK_EP_TH0, FC_RANK_EP_AR0);
                        mixer->UpdateBit0(FC_RANK_EM_LR0, FC_RANK_EM_LR1, FC_RANK_EM_LR2, FC_RANK_EM_TH0, FC_RANK_EM_AR0);
                        break;
                    }
                }

                tidHistory[currentChar] = bitRankSize;

                statePredictor  = & model->tid_t.bits_value[bitRankSize].state_model[state][0];
                charPredictor   = & model->tid_t.bits_value[bitRankSize].char_model[currentChar][0];
                staticPredictor = & model->tid_t.bits_value[bitRankSize].static_model[0];
                mixer           = & model->tid_mixerbits_value[bitRankSize];

                for (int bit = bitRankSize - 1; bit >= 0; --bit)
                {
                    if (coder.DecodeBit(mixer->Mixup(charPredictor[tid], statePredictor[tid], staticPredictor[tid])))
                    {
                        fc_count::UpdateBit1(statePredictor[tid],  FC_RANK_MS_TH1, FC_RANK_MS_AR1);
                        fc_count::UpdateBit1(charPredictor[tid],   FC_RANK_MC_TH1, FC_RANK_MC_AR1);
                        fc_count::UpdateBit1(staticPredictor[tid], FC_RANK_MP_TH1, FC_RANK_MP_AR1);
                        mixer->UpdateBit1(FC_RANK_MM_LR0, FC_RANK_MM_LR1, FC_RANK_MM_LR2, FC_RANK_MM_TH1, FC_RANK_MM_AR1);
                        tid += tid + 1;
                    }
                    else
                    {
                        fc_count::UpdateBit0(statePredictor[tid],  FC_RANK_MS_TH0, FC_RANK_MS_AR0);
                        fc_count::UpdateBit0(charPredictor[tid],   FC_RANK_MC_TH0, FC_RANK_MC_AR0);
                        fc_count::UpdateBit0(staticPredictor[tid], FC_RANK_MP_TH0, FC_RANK_MP_AR0);
                        mixer->UpdateBit0(FC_RANK_MM_LR0, FC_RANK_MM_LR1, FC_RANK_MM_LR2, FC_RANK_MM_TH0, FC_RANK_MM_AR0);
                        tid += tid;
                    }
                }
            }
            else
            {
                tidHistory[currentChar] = 0;
                fc_count::UpdateBit0(*statePredictor, FC_RANK_TS_TH0,  FC_RANK_TS_AR0);
                fc_count::UpdateBit0(*charPredictor, FC_RANK_TC_TH0,   FC_RANK_TC_AR0);
                fc_count::UpdateBit0(*staticPredictor, FC_RANK_TP_TH0, FC_RANK_TP_AR0);
                mixer->UpdateBit0(FC_RANK_TM_LR0, FC_RANK_TM_LR1, FC_RANK_TM_LR2, FC_RANK_TM_TH0, FC_RANK_TM_AR0);
            }
        }
        else
        {
            statePredictor  = & model->tid_t.out_of_range.state_model[state][0];
            charPredictor   = & model->tid_t.out_of_range.char_model[currentChar][0];
            staticPredictor = & model->tid_t.out_of_range.static_model[0];

            tid = 0;
            for (int context = 1, bit = maxRank; bit >= 0; --bit)
            {
                mixer = & model->tid_mixerout_of_range[context];

                if (coder.DecodeBit(mixer->Mixup(charPredictor[context], statePredictor[context], staticPredictor[context])))
                {
                    fc_count::UpdateBit1(statePredictor[context],  FC_RANK_PS_TH1, FC_RANK_PS_AR1);
                    fc_count::UpdateBit1(charPredictor[context],   FC_RANK_PC_TH1, FC_RANK_PC_AR1);
                    fc_count::UpdateBit1(staticPredictor[context], FC_RANK_PP_TH1, FC_RANK_PP_AR1);
                    mixer->UpdateBit1(FC_RANK_PM_LR0, FC_RANK_PM_LR1, FC_RANK_PM_LR2, FC_RANK_PM_TH1, FC_RANK_PM_AR1);
                    context += context + 1; tid += tid + 1;
                }
                else
                {
                    fc_count::UpdateBit0(statePredictor[context],  FC_RANK_PS_TH0, FC_RANK_PS_AR0);
                    fc_count::UpdateBit0(charPredictor[context],   FC_RANK_PC_TH0, FC_RANK_PC_AR0);
                    fc_count::UpdateBit0(staticPredictor[context], FC_RANK_PP_TH0, FC_RANK_PP_AR0);
                    mixer->UpdateBit0(FC_RANK_PM_LR0, FC_RANK_PM_LR1, FC_RANK_PM_LR2, FC_RANK_PM_TH0, FC_RANK_PM_AR0);
                    context += context; tid += tid;
                }
            }

            tidHistory[currentChar] = (unsigned char)fc_bit_scan_reverse(tid);
        }

        {
#if FC_CPU_TYPES >= FC_CPU_TYPES_SSE41
            __m128i * MTFTable_p = (__m128i *)&MTFTable[tid & (-16)];
            __m128i r = _mm_load_si128(MTFTable_p); _mm_store_si128(MTFTable_p, _mm_shuffle_epi8(_mm_insert_epi8(r, currentChar, 0), tid_reorder[tid & 15]));

            while ((--MTFTable_p) >= (__m128i *)MTFTable)
            {
                __m128i t = _mm_load_si128(MTFTable_p); _mm_store_si128(MTFTable_p, _mm_alignr_epi8(r, t, 1)); r = t;
            }
#elif FC_CPU_TYPES == FC_CPU_TYPES_A64
            uint8x16_t * MTFTable_p = (uint8x16_t *)&MTFTable[tid & (-16)];
            uint8x16_t r = vld1q_u8((const unsigned char *)MTFTable_p); vst1q_u8((unsigned char *)MTFTable_p, vqtbl1q_u8(vsetq_lane_u8((unsigned char)currentChar, r, 0), vld1q_u8((const unsigned char *)&tid_reorder[tid & 15][0])));

            while ((--MTFTable_p) >= (uint8x16_t *)MTFTable)
            {
                uint8x16_t t = vld1q_u8((const unsigned char *)MTFTable_p); vst1q_u8((unsigned char *)MTFTable_p, vextq_u8(t, r, 1)); r = t;
            }
#else
            for (int r = 0; r < tid; ++r)
            {
                MTFTable[r] = MTFTable[r + 1];
            }
            MTFTable[tid] = currentChar;
#endif
        }

        avgRank         =   (avgRank * 125 + tid * 3) >> 7;
        tid             =   tid - 1;
        history         =   runHistory[currentChar];
        int avgBucket   =   fc_avgRank_bucket(avgRank);
        state           =   lstate(contextRank0, contextRun, tid, history, avgBucket);
        statePredictor  = & model->run_t.state_model[state];
        charPredictor   = & model->run_t.char_model[currentChar];
        staticPredictor = & model->run_t.static_model;
        mixer           = & model->run_mixer[currentChar];

        int runSize = 1;
        if (coder.DecodeBit(mixer->Mixup(*charPredictor, *statePredictor, *staticPredictor)))
        {
            fc_count::UpdateBit1(*statePredictor,  FC_RUN_TS_TH1, FC_RUN_TS_AR1);
            fc_count::UpdateBit1(*charPredictor,   FC_RUN_TC_TH1, FC_RUN_TC_AR1);
            fc_count::UpdateBit1(*staticPredictor, FC_RUN_TP_TH1, FC_RUN_TP_AR1);
            mixer->UpdateBit1(FC_RUN_TM_LR0, FC_RUN_TM_LR1, FC_RUN_TM_LR2, FC_RUN_TM_TH1, FC_RUN_TM_AR1);

            statePredictor  = & model->run_t.bit_width.state_model[state][0];
            charPredictor   = & model->run_t.bit_width.char_model[currentChar][0];
            staticPredictor = & model->run_t.bit_width.static_model[0];
            mixer           = & model->run_mixerbit_width[history < 1 ? 1 : history][1];

            int bitRunSize = 1;
            while (true)
            {
                if (coder.DecodeBit(mixer->Mixup(*charPredictor, *statePredictor, *staticPredictor)))
                {
                    fc_count::UpdateBit1(*statePredictor,  FC_RUN_ES_TH1, FC_RUN_ES_AR1); statePredictor++;
                    fc_count::UpdateBit1(*charPredictor,   FC_RUN_EC_TH1, FC_RUN_EC_AR1); charPredictor++;
                    fc_count::UpdateBit1(*staticPredictor, FC_RUN_EP_TH1, FC_RUN_EP_AR1); staticPredictor++;
                    mixer->UpdateBit1(FC_RUN_EM_LR0, FC_RUN_EM_LR1, FC_RUN_EM_LR2, FC_RUN_EM_TH1, FC_RUN_EM_AR1);
                    bitRunSize++; mixer = & model->run_mixerbit_width[history < bitRunSize ? bitRunSize : history][bitRunSize];
                }
                else
                {
                    fc_count::UpdateBit0(*statePredictor,  FC_RUN_ES_TH0, FC_RUN_ES_AR0);
                    fc_count::UpdateBit0(*charPredictor,   FC_RUN_EC_TH0, FC_RUN_EC_AR0);
                    fc_count::UpdateBit0(*staticPredictor, FC_RUN_EP_TH0, FC_RUN_EP_AR0);
                    mixer->UpdateBit0(FC_RUN_EM_LR0, FC_RUN_EM_LR1, FC_RUN_EM_LR2, FC_RUN_EM_TH0, FC_RUN_EM_AR0);
                    break;
                }
            }

            runHistory[currentChar] = (runHistory[currentChar] + 3 * bitRunSize + 3) >> 2;

            statePredictor  = & model->run_t.bits_value[bitRunSize].state_model[state][0];
            charPredictor   = & model->run_t.bits_value[bitRunSize].char_model[currentChar][0];
            staticPredictor = & model->run_t.bits_value[bitRunSize].static_model[0];
            mixer           = & model->run_mixerbits_value[bitRunSize];

            for (int context = 1, bit = bitRunSize - 1; bit >= 0; --bit)
            {
                if (coder.DecodeBit(mixer->Mixup(charPredictor[context], statePredictor[context], staticPredictor[context])))
                {
                    fc_count::UpdateBit1(statePredictor[context],  FC_RUN_MS_TH1, FC_RUN_MS_AR1);
                    fc_count::UpdateBit1(charPredictor[context],   FC_RUN_MC_TH1, FC_RUN_MC_AR1);
                    fc_count::UpdateBit1(staticPredictor[context], FC_RUN_MP_TH1, FC_RUN_MP_AR1);
                    mixer->UpdateBit1(FC_RUN_MM_LR0, FC_RUN_MM_LR1, FC_RUN_MM_LR2, FC_RUN_MM_TH1, FC_RUN_MM_AR1);
                    runSize += runSize + 1; if (bitRunSize <= 5) context += context + 1; else context++;
                }
                else
                {
                    fc_count::UpdateBit0(statePredictor[context],  FC_RUN_MS_TH0, FC_RUN_MS_AR0);
                    fc_count::UpdateBit0(charPredictor[context],   FC_RUN_MC_TH0, FC_RUN_MC_AR0);
                    fc_count::UpdateBit0(staticPredictor[context], FC_RUN_MP_TH0, FC_RUN_MP_AR0);
                    mixer->UpdateBit0(FC_RUN_MM_LR0, FC_RUN_MM_LR1, FC_RUN_MM_LR2, FC_RUN_MM_TH0, FC_RUN_MM_AR0);
                    runSize += runSize; if (bitRunSize <= 5) context += context; else context++;
                }
            }

        }
        else
        {
            runHistory[currentChar] = (runHistory[currentChar] + 2) >> 2;
            fc_count::UpdateBit0(*statePredictor,  FC_RUN_TS_TH0, FC_RUN_TS_AR0);
            fc_count::UpdateBit0(*charPredictor,   FC_RUN_TC_TH0, FC_RUN_TC_AR0);
            fc_count::UpdateBit0(*staticPredictor, FC_RUN_TP_TH0, FC_RUN_TP_AR0);
            mixer->UpdateBit0(FC_RUN_TM_LR0, FC_RUN_TM_LR1, FC_RUN_TM_LR2, FC_RUN_TM_TH0, FC_RUN_TM_AR0);
        }

        contextRank0 = ((contextRank0 << 1) | (tid == 0 ? 1 : 0)) & 0x7;
        int tidBucket = (tid < 3 ? tid : 3);
        avgBucket = fc_avgRank_bucket(avgRank);
        int rank_hist = ((contextRank4 & 0x3f) << 2) | tidBucket;
        rank_hist &= 0x3f;
        contextRank4 = rank_hist | (avgBucket << 6);
        contextRun = ((contextRun << 1) | (runSize < 3 ? 1 : 0)) & 0xf;

        for (; runSize > 0; --runSize) output[i++] = currentChar;
    }

    return n;
}

int fcinit()
{
    return fcinit_model();
}

int fc_encode(const unsigned char * input, unsigned char * output, int inputSize, int outputSize)
{
    if (fcModel1 * model = (fcModel1 *)malloc(sizeof(fcModel1)))
    {
        if (unsigned char * buffer = (unsigned char *)malloc(inputSize * sizeof(unsigned char)))
        {
            int result = fc_encode_do(input, output, buffer, inputSize, outputSize, model);

            free(buffer); free(model);

            return result;
        };
        free(model);
    }
    return FC_LACK_OF_MEMORY;
}


int fc_decode(const unsigned char * input, unsigned char * output)
{
    if (fcModel1 * model = (fcModel1 *)malloc(sizeof(fcModel1)))
    {
        int result = fc_decode_do(input, output, model);

        free(model);

        return result;
    }
    return FC_LACK_OF_MEMORY;
}