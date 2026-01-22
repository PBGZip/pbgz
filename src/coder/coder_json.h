/*
 * coder_json.h - Header file for pbgz project
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

#ifndef _CODER_JSON_H_
#define _CODER_JSON_H_

#include <json/json.h>
#include <zstd.h>

#include "coder.h"

/* JSON encoder using zstd streaming compression algorithm */
class coder_json
{
public:
    /* level [1, 7] */
    coder_json(int32_t level=6);
    virtual ~coder_json();

    /* Compress JSON data, return compressed length, negative value indicates insufficient out space, e.g., return -28 means out_len needs 28 more bytes
     * This function assumes out_len is sufficient, otherwise it will exit with error if out_len is insufficient
     */
    virtual int64_t encoder(const Json::Value &in, uint8_t *out, const int64_t out_len);

    /* Compress JSON data, allocate space internally, compressed data stored in out */
    virtual void encoder(const Json::Value &in, std::string &out);

    /* Decompress JSON format data */
    virtual void decoder(const uint8_t *in, const int64_t in_len, Json::Value &out);

private:
    /* Corresponding compression level */
    int32_t cLevel;
};

#endif
