/*
 * base_pack.h - map A/C/G/T bases to 2-bit codes and pack four per byte.
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <cstdint>
#include <cstring>

#if defined(__SSSE3__) && defined(__SSE4_1__)
#include <tmmintrin.h>
#include <smmintrin.h>
#endif
#if defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#endif

namespace pbgz_base_pack {

static inline uint8_t baseCode(uint8_t c)
{
    return (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : (c == 'T') ? 3 : 0xFF;
}

/* Scalar fallback / tail. Packs ceil(n/4) bytes, returns false on a non-ACGT byte. */
static inline bool scalar(const uint8_t *in, uint32_t n, uint8_t *out)
{
    uint32_t i = 0, o = 0;
    for (; i + 4 <= n; i += 4) {
        const uint8_t a = baseCode(in[i + 0]);
        const uint8_t b = baseCode(in[i + 1]);
        const uint8_t c = baseCode(in[i + 2]);
        const uint8_t d = baseCode(in[i + 3]);
        if ((a | b | c | d) > 3) {
            return false;
        }
        out[o++] = (uint8_t)(a | (b << 2) | (c << 4) | (d << 6));
    }
    if (i < n) {
        uint8_t pk = 0;
        uint32_t shift = 0;
        for (; i < n; ++i, shift += 2) {
            const uint8_t a = baseCode(in[i]);
            if (a > 3) {
                return false;
            }
            pk = (uint8_t)(pk | (a << shift));
        }
        out[o++] = pk;
    }
    return true;
}

#if defined(__SSSE3__) && defined(__SSE4_1__)

static inline bool x86(const uint8_t *in, uint32_t n, uint8_t *out)
{
    const __m128i vA = _mm_set1_epi8('A');
    const __m128i vC = _mm_set1_epi8('C');
    const __m128i vG = _mm_set1_epi8('G');
    const __m128i vT = _mm_set1_epi8('T');
    const __m128i one = _mm_set1_epi8(1);
    const __m128i two = _mm_set1_epi8(2);
    const __m128i three = _mm_set1_epi8(3);
    const __m128i w1 = _mm_set1_epi16(0x0401);     /* pair: 1*lo + 4*hi        */
    const __m128i w2 = _mm_set1_epi32(0x00100001); /* nibble: 1*lo + 16*hi     */

    uint32_t i = 0, o = 0;
    for (; i + 16 <= n; i += 16) {
        const __m128i x = _mm_loadu_si128(reinterpret_cast<const __m128i *>(in + i));
        const __m128i cA = _mm_cmpeq_epi8(x, vA);
        const __m128i cC = _mm_cmpeq_epi8(x, vC);
        const __m128i cG = _mm_cmpeq_epi8(x, vG);
        const __m128i cT = _mm_cmpeq_epi8(x, vT);
        const __m128i valid = _mm_or_si128(_mm_or_si128(cA, cC), _mm_or_si128(cG, cT));
        if (_mm_movemask_epi8(valid) != 0xFFFF) {
            /* some byte in this chunk is not A/C/G/T: fall back to scalar for it
               and for the rest, so the whole stream is handled consistently. */
            return scalar(in + i, n - i, out + o);
        }
        const __m128i code = _mm_or_si128(
            _mm_and_si128(cC, one),
            _mm_or_si128(_mm_and_si128(cG, two), _mm_and_si128(cT, three)));
        const __m128i nib = _mm_maddubs_epi16(code, w1);
        const __m128i byt = _mm_madd_epi16(nib, w2);
        __m128i pk = _mm_packus_epi32(byt, byt);
        pk = _mm_packus_epi16(pk, pk);
        memcpy(out + o, &pk, 4);
        o += 4;
    }
    if (i < n) {
        return scalar(in + i, n - i, out + o);
    }
    return true;
}

#endif /* SSSE3 && SSE4.1 */

#if defined(__aarch64__) && defined(__ARM_NEON)

static inline bool arm(const uint8_t *in, uint32_t n, uint8_t *out)
{
    const uint8x16_t vA = vdupq_n_u8('A');
    const uint8x16_t vC = vdupq_n_u8('C');
    const uint8x16_t vG = vdupq_n_u8('G');
    const uint8x16_t vT = vdupq_n_u8('T');
    const uint8x16_t one = vdupq_n_u8(1);
    const uint8x16_t two = vdupq_n_u8(2);
    const uint8x16_t three = vdupq_n_u8(3);

    uint32_t i = 0, o = 0;
    for (; i + 16 <= n; i += 16) {
        const uint8x16_t x = vld1q_u8(in + i);
        const uint8x16_t cA = vceqq_u8(x, vA);
        const uint8x16_t cC = vceqq_u8(x, vC);
        const uint8x16_t cG = vceqq_u8(x, vG);
        const uint8x16_t cT = vceqq_u8(x, vT);
        const uint8x16_t valid = vorrq_u8(vorrq_u8(cA, cC), vorrq_u8(cG, cT));
        if (vminvq_u8(valid) != 0xFF) {
            return scalar(in + i, n - i, out + o);
        }
        const uint8x16_t code = vorrq_u8(vandq_u8(cC, one),
                                         vorrq_u8(vandq_u8(cG, two), vandq_u8(cT, three)));
        /* code is 16 bytes, each 0..3. Pack two codes into a 4-bit nibble, then
           two nibbles into a byte, entirely in the 16/32-bit lanes. */
        const uint16x8_t c16 = vreinterpretq_u16_u8(code); /* lo=code[2j], hi=code[2j+1] */
        const uint16x8_t nib = vorrq_u16(vandq_u16(c16, vdupq_n_u16(0x0003)),
                                         vandq_u16(vshrq_n_u16(c16, 6), vdupq_n_u16(0x000C)));
        const uint32x4_t n32 = vreinterpretq_u32_u16(nib);
        const uint32x4_t byt = vorrq_u32(vandq_u32(n32, vdupq_n_u32(0x0000000Fu)),
                                         vandq_u32(vshrq_n_u32(n32, 12), vdupq_n_u32(0x000000F0u)));
        const uint16x4_t t16 = vmovn_u32(byt);
        const uint8x8_t t8 = vmovn_u16(vcombine_u16(t16, t16));
        vst1_lane_u32(reinterpret_cast<uint32_t *>(out + o), vreinterpret_u32_u8(t8), 0);
        o += 4;
    }
    if (i < n) {
        return scalar(in + i, n - i, out + o);
    }
    return true;
}

#endif /* aarch64 && NEON */

static inline bool packBases(const uint8_t *in, uint32_t n, uint8_t *out)
{
#if defined(__SSSE3__) && defined(__SSE4_1__)
    return x86(in, n, out);
#elif defined(__aarch64__) && defined(__ARM_NEON)
    return arm(in, n, out);
#else
    return scalar(in, n, out);
#endif
}

} /* namespace pbgz_base_pack */
