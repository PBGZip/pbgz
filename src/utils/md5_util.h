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

/* MD5上下文结构 */
typedef struct {
    uint32_t state[4];    /* 状态 (ABCD) */
    uint32_t count[2];    /* 比特数 (低位, 高位) */
    uint8_t buffer[64];   /* 输入缓冲区 */
} md5_ctx;

/* 函数声明 */

/**
 * 初始化MD5上下文
 * @param context 要初始化的上下文
 */
void md5_init(md5_ctx *context);

/**
 * 更新MD5计算（可多次调用）
 * @param context MD5上下文
 * @param input 输入数据
 * @param input_len 输入数据的长度
 */
void md5_update(md5_ctx *context, const uint8_t *input, size_t input_len);

/**
 * 完成MD5计算，输出摘要
 * @param digest 输出的16字节摘要（二进制格式）
 * @param context MD5上下文
 */
void md5_final(uint8_t digest[16], md5_ctx *context);

/**
 * 方便的包装函数：直接计算字符串的MD5
 * @param input 输入字符串
 * @param length 输入长度（0表示自动计算以null结尾的字符串长度）
 * @param digest 输出的16字节摘要
 */
void md5_compute(const uint8_t *input, size_t length, uint8_t digest[16]);

/**
 * 将16字节二进制摘要转换为33字节的十六进制字符串（包含结尾的null）
 * @param digest 16字节的二进制摘要
 * @param output 33字节的输出缓冲区
 */
void md5_to_hex(const uint8_t digest[16], char output[33]);

/**
 * 将16字节二进制摘要转换为字符串对象
 * @param digest 16字节的二进制摘要
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