/*
 * md5_util.h - Header file for pbgz project
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

#include "md5_util.h"
#include <string.h>

/* Basic MD5 function definitions */
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) ((y) ^ ((x) | (~z)))

/* Rotate left */
#define ROTATE_LEFT(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

/* One step in four-round transformation */
#define FF(a, b, c, d, x, s, ac) { \
    (a) += F((b), (c), (d)) + (x) + (uint32_t)(ac); \
    (a) = ROTATE_LEFT((a), (s)); \
    (a) += (b); \
}

#define GG(a, b, c, d, x, s, ac) { \
    (a) += G((b), (c), (d)) + (x) + (uint32_t)(ac); \
    (a) = ROTATE_LEFT((a), (s)); \
    (a) += (b); \
}

#define HH(a, b, c, d, x, s, ac) { \
    (a) += H((b), (c), (d)) + (x) + (uint32_t)(ac); \
    (a) = ROTATE_LEFT((a), (s)); \
    (a) += (b); \
}

#define II(a, b, c, d, x, s, ac) { \
    (a) += I((b), (c), (d)) + (x) + (uint32_t)(ac); \
    (a) = ROTATE_LEFT((a), (s)); \
    (a) += (b); \
}

/* Byte order conversion (little-endian system) */
#ifndef WORDS_BIGENDIAN
#define byte_reverse(words, byte_count) /* Little-endian system doesn't need reversal */
#else
static void byte_reverse(uint8_t *buf, size_t words) {
    uint32_t t;
    do {
        t = (uint32_t)((uint16_t)(buf[3] << 8 | buf[2])) << 16 |
            ((uint16_t)(buf[1] << 8 | buf[0]));
        *(uint32_t *)buf = t;
        buf += 4;
    } while (--words);
}
#endif

/* S table - shift amounts for each round */
static const uint8_t S[4][4] = {
    { 7, 12, 17, 22 },
    { 5,  9, 14, 20 },
    { 4, 11, 16, 23 },
    { 6, 10, 15, 21 }
};

/* T table - constants */
static const uint32_t T[64] = {
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee,
    0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
    0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
    0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
    0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
    0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
    0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
    0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
    0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
    0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
    0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
    0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039,
    0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
    0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
    0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391
};

/* MD5 core transformation, processes one 64-byte block */
static void md5_transform(uint32_t state[4], const uint8_t block[64]) {
    uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
    uint32_t x[16];
    
    /* Convert block to 16 32-bit words */
#ifdef WORDS_BIGENDIAN
    memcpy(x, block, 64);
    byte_reverse((uint8_t *)x, 16);
#else
    memcpy(x, block, 64);
#endif
    
    /* Round 1 */
    FF(a, b, c, d, x[ 0], S[0][0], T[ 0]);
    FF(d, a, b, c, x[ 1], S[0][1], T[ 1]);
    FF(c, d, a, b, x[ 2], S[0][2], T[ 2]);
    FF(b, c, d, a, x[ 3], S[0][3], T[ 3]);
    FF(a, b, c, d, x[ 4], S[0][0], T[ 4]);
    FF(d, a, b, c, x[ 5], S[0][1], T[ 5]);
    FF(c, d, a, b, x[ 6], S[0][2], T[ 6]);
    FF(b, c, d, a, x[ 7], S[0][3], T[ 7]);
    FF(a, b, c, d, x[ 8], S[0][0], T[ 8]);
    FF(d, a, b, c, x[ 9], S[0][1], T[ 9]);
    FF(c, d, a, b, x[10], S[0][2], T[10]);
    FF(b, c, d, a, x[11], S[0][3], T[11]);
    FF(a, b, c, d, x[12], S[0][0], T[12]);
    FF(d, a, b, c, x[13], S[0][1], T[13]);
    FF(c, d, a, b, x[14], S[0][2], T[14]);
    FF(b, c, d, a, x[15], S[0][3], T[15]);
    
    /* Round 2 */
    GG(a, b, c, d, x[ 1], S[1][0], T[16]);
    GG(d, a, b, c, x[ 6], S[1][1], T[17]);
    GG(c, d, a, b, x[11], S[1][2], T[18]);
    GG(b, c, d, a, x[ 0], S[1][3], T[19]);
    GG(a, b, c, d, x[ 5], S[1][0], T[20]);
    GG(d, a, b, c, x[10], S[1][1], T[21]);
    GG(c, d, a, b, x[15], S[1][2], T[22]);
    GG(b, c, d, a, x[ 4], S[1][3], T[23]);
    GG(a, b, c, d, x[ 9], S[1][0], T[24]);
    GG(d, a, b, c, x[14], S[1][1], T[25]);
    GG(c, d, a, b, x[ 3], S[1][2], T[26]);
    GG(b, c, d, a, x[ 8], S[1][3], T[27]);
    GG(a, b, c, d, x[13], S[1][0], T[28]);
    GG(d, a, b, c, x[ 2], S[1][1], T[29]);
    GG(c, d, a, b, x[ 7], S[1][2], T[30]);
    GG(b, c, d, a, x[12], S[1][3], T[31]);
    
    /* Round 3 */
    HH(a, b, c, d, x[ 5], S[2][0], T[32]);
    HH(d, a, b, c, x[ 8], S[2][1], T[33]);
    HH(c, d, a, b, x[11], S[2][2], T[34]);
    HH(b, c, d, a, x[14], S[2][3], T[35]);
    HH(a, b, c, d, x[ 1], S[2][0], T[36]);
    HH(d, a, b, c, x[ 4], S[2][1], T[37]);
    HH(c, d, a, b, x[ 7], S[2][2], T[38]);
    HH(b, c, d, a, x[10], S[2][3], T[39]);
    HH(a, b, c, d, x[13], S[2][0], T[40]);
    HH(d, a, b, c, x[ 0], S[2][1], T[41]);
    HH(c, d, a, b, x[ 3], S[2][2], T[42]);
    HH(b, c, d, a, x[ 6], S[2][3], T[43]);
    HH(a, b, c, d, x[ 9], S[2][0], T[44]);
    HH(d, a, b, c, x[12], S[2][1], T[45]);
    HH(c, d, a, b, x[15], S[2][2], T[46]);
    HH(b, c, d, a, x[ 2], S[2][3], T[47]);
    
    /* Round 4 */
    II(a, b, c, d, x[ 0], S[3][0], T[48]);
    II(d, a, b, c, x[ 7], S[3][1], T[49]);
    II(c, d, a, b, x[14], S[3][2], T[50]);
    II(b, c, d, a, x[ 5], S[3][3], T[51]);
    II(a, b, c, d, x[12], S[3][0], T[52]);
    II(d, a, b, c, x[ 3], S[3][1], T[53]);
    II(c, d, a, b, x[10], S[3][2], T[54]);
    II(b, c, d, a, x[ 1], S[3][3], T[55]);
    II(a, b, c, d, x[ 8], S[3][0], T[56]);
    II(d, a, b, c, x[15], S[3][1], T[57]);
    II(c, d, a, b, x[ 6], S[3][2], T[58]);
    II(b, c, d, a, x[13], S[3][3], T[59]);
    II(a, b, c, d, x[ 4], S[3][0], T[60]);
    II(d, a, b, c, x[11], S[3][1], T[61]);
    II(c, d, a, b, x[ 2], S[3][2], T[62]);
    II(b, c, d, a, x[ 9], S[3][3], T[63]);
    
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
}

/* Initialize MD5 context */
void md5_init(md5_ctx *context) {
    context->count[0] = context->count[1] = 0;
    
    /* Load magic numbers */
    context->state[0] = 0x67452301;
    context->state[1] = 0xefcdab89;
    context->state[2] = 0x98badcfe;
    context->state[3] = 0x10325476;
}

/* Update MD5 calculation */
void md5_update(md5_ctx *context, const uint8_t *input, size_t input_len) {
    size_t i, index, part_len;
    
    /* Calculate number of bytes in current buffer */
    index = (size_t)((context->count[0] >> 3) & 0x3F);
    
    /* Update bit count */
    if ((context->count[0] += ((uint32_t)input_len << 3)) < ((uint32_t)input_len << 3))
        context->count[1]++;
    context->count[1] += ((uint32_t)input_len >> 29);
    
    part_len = 64 - index;
    
    /* Fill buffer as much as possible and transform */
    if (input_len >= part_len) {
        memcpy(&context->buffer[index], input, part_len);
        md5_transform(context->state, context->buffer);
        
        for (i = part_len; i + 63 < input_len; i += 64)
            md5_transform(context->state, &input[i]);
        
        index = 0;
    } else {
        i = 0;
    }
    
    /* Remaining input buffer */
    memcpy(&context->buffer[index], &input[i], input_len - i);
}

/* Finalize MD5 calculation */
void md5_final(uint8_t digest[16], md5_ctx *context) {
    uint8_t bits[8];
    size_t index, pad_len;
    
    /* Save bit count */
    for (int i = 0; i < 8; i++)
        bits[i] = (uint8_t)((context->count[i >> 2] >> ((i & 3) * 8)) & 0xFF);
    
    /* Padding: add 1 bit followed by enough 0 bits to make length ≡ 448 mod 512 */
    index = (size_t)((context->count[0] >> 3) & 0x3F);
    pad_len = (index < 56) ? (56 - index) : (120 - index);
    md5_update(context, (const uint8_t *)"\x80", 1);
    while (pad_len-- > 1)
        md5_update(context, (const uint8_t *)"\0", 1);
    
    /* Append length (in bits) */
    md5_update(context, bits, 8);
    
    /* Store state to digest */
    for (int i = 0; i < 4; i++) {
        digest[i * 4] = (uint8_t)(context->state[i] & 0xFF);
        digest[i * 4 + 1] = (uint8_t)((context->state[i] >> 8) & 0xFF);
        digest[i * 4 + 2] = (uint8_t)((context->state[i] >> 16) & 0xFF);
        digest[i * 4 + 3] = (uint8_t)((context->state[i] >> 24) & 0xFF);
    }
    
    /* Clear sensitive information */
    memset(context, 0, sizeof(*context));
}

/* Convenient wrapper function */
void md5_compute(const uint8_t *input, size_t length, uint8_t digest[16]) {
    md5_ctx context;
    md5_init(&context);
    md5_update(&context, input, length);
    md5_final(digest, &context);
}

/* Convert to hexadecimal string */
void md5_to_hex(const uint8_t digest[16], char output[33]) {
    static const char hex[] = "0123456789abcdef";
    
    for (int i = 0; i < 16; i++) {
        output[i * 2] = hex[digest[i] >> 4];
        output[i * 2 + 1] = hex[digest[i] & 0x0F];
    }
    output[32] = '\0';
}

std::string md5_to_string(const uint8_t digest[16]) {
    static const char hex[] = "0123456789abcdef";
    
    char output[33];
    for (int i = 0; i < 16; i++) {
        output[i * 2] = hex[digest[i] >> 4];
        output[i * 2 + 1] = hex[digest[i] & 0x0F];
    }
    output[32] = '\0';
    return std::string(output);
}
