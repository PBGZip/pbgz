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

#pragma once

#include <string>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stddef.h>

/* MD5 context structure */
typedef struct {
    uint32_t state[4];    /* state (ABCD) */
    uint32_t count[2];    /* bit count (low, high) */
    uint8_t buffer[64];   /* input buffer */
} md5_ctx;

/* Function declarations */

/**
 * Initialize MD5 context
 * @param context The context to initialize
 */
void md5_init(md5_ctx *context);

/**
 * Update MD5 calculation (can be called multiple times)
 * @param context MD5 context
 * @param input Input data
 * @param input_len Length of input data
 */
void md5_update(md5_ctx *context, const uint8_t *input, size_t input_len);

/**
 * Finalize MD5 calculation and output digest
 * @param digest Output 16-byte digest (binary format)
 * @param context MD5 context
 */
void md5_final(uint8_t digest[16], md5_ctx *context);

/**
 * Convenience wrapper function: directly compute MD5 of string
 * @param input Input string
 * @param length Input length (0 means automatically calculate null-terminated string length)
 * @param digest Output 16-byte digest
 */
void md5_compute(const uint8_t *input, size_t length, uint8_t digest[16]);

/**
 * Convert 16-byte binary digest to 33-byte hexadecimal string (including trailing null)
 * @param digest 16-byte binary digest
 * @param output 33-byte output buffer
 */
void md5_to_hex(const uint8_t digest[16], char output[33]);

/**
 * Convert 16-byte binary digest to string object
 * @param digest 16-byte binary digest
 * @return std::string
 */
std::string md5_to_string(const uint8_t digest[16]);

#ifdef __cplusplus
}
#endif

inline void calcMd5sum(std::string& md5, const uint8_t* data, uint32_t dataLen) {
    uint8_t digest[16];
    char hex[33];
    md5_compute((const uint8_t *)data, dataLen, digest);
    md5_to_hex(digest, hex);
    md5 = std::string(hex);
}
