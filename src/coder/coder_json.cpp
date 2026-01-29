/*
 * coder_json.cpp - Source file for pbgz project
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

#include "coder_json.h"

coder_json::coder_json(int32_t level) : cLevel(level)
{
}

coder_json::~coder_json()
{
}

/* Compress JSON data, return compressed length. If negative, it means out space is insufficient, 
 * e.g., return -28 means out_len needs to be increased by 28.
 * This function assumes out_len is sufficient, otherwise it will report an error and exit if out_len is insufficient.
 */
int64_t coder_json::encoder(const Json::Value &in, uint8_t *out, const int64_t out_len)
{
    Json::StreamWriterBuilder builder;
    if (in.empty())
        return 0;
    const std::string encjson = Json::writeString(builder, in);
    ZSTD_CCtx *const cctx = ZSTD_createCCtx();
    if (cctx ==  nullptr) {
        return -1;
    }
    
    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, cLevel), coder_ns::CODER_ERR_BAD_ARGS, "ZSTD_CCtx_setParameter cLevel failed!");
    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1), coder_ns::CODER_ERR_INNER, "ZSTD_CCtx_setParameter checksum flag failed!");

    ZSTD_EndDirective const mode = ZSTD_e_end;
    ZSTD_inBuffer input = {encjson.c_str(), encjson.length(), 0};
    ZSTD_outBuffer output = {out, (size_t)(out_len), 0};
    size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
    if (remaining != 0)
    {
        ZSTD_freeCCtx(cctx);
        coder_exit(coder_ns::CODER_ERR_INNER, "ZSTD_compressStream2 error : out buffer lack of %ld bytes\n", remaining);
    }
    if (input.pos != input.size)
    {
        coder_exit(coder_ns::CODER_ERR_INNER, "ZSTD_compressStream2 error : input.pos %ld != input.size %ld\n", input.pos, input.size);
    }

    ZSTD_freeCCtx(cctx);
    return output.pos;
}

/* Compress JSON data, allocate space internally, compressed data is stored in out */
void coder_json::encoder(const Json::Value &in, std::string &out)
{
    out.clear();
    Json::StreamWriterBuilder builder;
    if (in.empty())
        return;
    const std::string encjson = Json::writeString(builder, in);
    ZSTD_CCtx *const cctx = ZSTD_createCCtx();
    check_exit(cctx != NULL, coder_ns::CODER_ERR_INNER, "ZSTD_createCCtx() failed!");

    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, cLevel), coder_ns::CODER_ERR_INNER, "ZSTD_CCtx_setParameter cLevel failed!");
    check_exit(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1), coder_ns::CODER_ERR_INNER, "ZSTD_CCtx_setParameter checksum flag failed!");

    ZSTD_EndDirective const mode = ZSTD_e_end;
    ZSTD_inBuffer input = {encjson.c_str(), encjson.length(), 0};
    uint8_t *buffOut;
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    buffOut = static_cast<uint8_t*>(safe_alloc(buffOutSize));
    check_exit(buffOut, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory alloc failed");

    int32_t finished;
    do
    {
        /* Compress into the output buffer and write all of the output to
             * the file so we can reuse the buffer next iteration.
             */
        ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
        size_t const remaining = ZSTD_compressStream2(cctx, &output, &input, mode);
        if (ZSTD_isError(remaining))
        {
            safe_free((void**)&buffOut);
            coder_exit(coder_ns::CODER_ERR_INNER, "ZSTD_compressStream2 error : %s \n", ZSTD_getErrorName(remaining));
        }
        out += std::string((char *)buffOut, output.pos);
        /* If we're on the last chunk we're finished when zstd returns 0,
             * which means its consumed all the input AND finished the frame.
             * Otherwise, we're finished when we've consumed all the input.
             */
        finished = (remaining == 0);
    } while (!finished);

    if (input.pos != input.size)
    {
        safe_free((void**)&buffOut);
        coder_exit(coder_ns::CODER_ERR_INNER, "ZSTD_compressStream2 error : input.pos %ld != input.size %ld\n", input.pos, input.size);
    }

    ZSTD_freeCCtx(cctx);
    safe_free((void**)&buffOut);
}

/* Decompress JSON format data */
void coder_json::decoder(const uint8_t *in, const int64_t in_len, Json::Value &out)
{
    JSONCPP_STRING err;
    Json::CharReaderBuilder builder;
    if (in_len <= 0)
        return;

    ZSTD_DCtx *const dctx = ZSTD_createDCtx();
    check_exit(dctx != NULL, coder_ns::CODER_ERR_INNER, "ZSTD_createDCtx() failed!");

    size_t lastRet = 0;
    std::string outdata = "";
    size_t const buffOutSize = ZSTD_DStreamOutSize();
    uint8_t *buffOut = static_cast<uint8_t*>(safe_alloc(buffOutSize));
    ZSTD_inBuffer input = {in, (size_t)in_len, 0};
    while (input.pos < input.size)
    {
        ZSTD_outBuffer output = {buffOut, buffOutSize, 0};
        size_t const ret = ZSTD_decompressStream(dctx, &output, &input);

        if (ZSTD_isError(ret))
        {
            safe_free((void**)&buffOut);
            coder_exit(coder_ns::CODER_ERR_INNER, "ZSTD_decompress error : %s ", ZSTD_getErrorName(ret));
        }

        outdata += std::string((char*)buffOut, output.pos);
        lastRet = ret;
    }

    if (lastRet != 0)
    {
        /* The last return value from ZSTD_decompressStream did not end on a
         * frame, but we reached the end of the file! We assume this is an
         * error, and the input was truncated.
         */
        ZSTD_freeDCtx(dctx);
        safe_free((void**)&buffOut);
        coder_exit(coder_ns::CODER_ERR_INNER, "EOF before end of stream: %zu\n", lastRet);
    }

    const std::unique_ptr<Json::CharReader> reader(builder.newCharReader());
    if (!reader->parse((char *)(outdata.c_str()), (char *)(outdata.c_str()) + outdata.size(), &out,
                       &err))
    {
        ZSTD_freeDCtx(dctx);
        safe_free((void**)&buffOut);
        check_exit(false, coder_ns::CODER_ERR_INNER, "json decoder failed");
    }

    ZSTD_freeDCtx(dctx);
    safe_free((void**)&buffOut);
}
