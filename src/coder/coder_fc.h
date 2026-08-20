/*
 * coder_fc.h - Header file for pbgz project
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

#include "coder_io.h"
#include "coder.h"
#include "fc/pp.h"
#include "fc/transform.h"
#include "fc/fc.h"

#define FC_MIN_LEN 32
// #define FC_MIN_LEN 2048 // Set a conservative value first, because if too small, the encoded length may be greater than the original length, which causes fc compression errors and makes the above process difficult to handle
#define FC_MAX_LEN 2146435072

class coder_fc : public coder
{
public:
    coder_fc(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_fc");
        flushed = false;
        lzp_valid = true;

        // if (this->io->get_level() == 1)
        if (true)
        {
            this->ppHashSize = 17;
            this->ppMinLen = 28;
            // 29 -> 123616168, better than 28
            // 27 -> 75368634
            // 26 -> 75365398
            // 20 -> 75350908
            // 18 -> 75347944
            // 17 -> 75361248
            // 16 -> 75540925
        }
        else
            check_exit(false,  coder_ns::CODER_ERR_INNER, "Currently level %d is not enabled", this->io->get_level());

        check_exit(!(ppMinLen < 4 || ppMinLen > 255),  coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, ppMinLen);
        check_exit(!(ppHashSize < 10 || ppHashSize > 28),  coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, ppMinLen);
    }

    virtual ~coder_fc()
    {
        if (io->m == coder_io::MENC && !flushed)
            encode_flush();
    }

    /*
     * coder_fc 只支持整块压缩：encode_line 内部一次性做完 LZP、BWT 和熵编码，
     * 首行的 check_exit 会拦住第二次调用，输入长度也必须大于 FC_MIN_LEN。
     */
    bool supportsLineMode() const override { return false; }

    void encode_line(const uint8_t *in, const uint32_t in_len, bool need2hold __attribute__ ((unused)) = false) override
    {
        check_exit(io->m != coder_io::MENC, coder_ns::CODER_ERR_INNER, "only support block compress, not support line method"); // Not yet extended to line method
        check_exit(in_len > FC_MIN_LEN && in_len < FC_MAX_LEN,  coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, in_len);
        io->m = coder_io::MENC;

        // First do pp prehandler
        uint8_t *lout;
        lout = static_cast<uint8_t*>(safe_alloc(in_len));
        check_exit(lout, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");
        int lzSize;
        int result = fc_preprocess(in, in + in_len, lout + 1, lout + in_len - 1, ppHashSize, ppMinLen);
        if (result >= 0)
        {
            result = (lout[0] = 1, result + 1);
            lzSize = result;
        }
        else {
            memcpy(lout, in, in_len);
            lzSize = in_len;
            lzp_valid = false;
        }
        check_exit(lzSize > 0, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d, in_len %d", __LINE__, lzSize, in_len);

        // Then do bwt
        index = fc_transform(lout, lzSize, &num_indexes, indexes);
        if (in_len < 64 * 1024)
            num_indexes = 0;
        check_exit(index >= FC_OK, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d, in_len %d", __LINE__, index, in_len);

        // Then compress
        uint8_t *fc_buff;
        bool need_alloc = (static_cast<uint32_t>(lzSize) + 4096) > in_len;
        if (need_alloc)
        {
            fc_buff = static_cast<uint8_t*>(safe_alloc((lzSize << 1) + 4096));
            check_exit(fc_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");
        }
        else{
            fc_buff = this->io->data;
        }

        // result = fcinit();
        // check_exit(FC_OK == result, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d, in_len %d", __LINE__, result, in_len);
        result = fc_encode(lout, fc_buff + 1, lzSize, (lzSize << 1) - 1);
        check_exit(result >= 0, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d, in_len %d", __LINE__, result, in_len);
        result = (fc_buff[0] = 1, result + 1);

        if (need_alloc)
        {
            memcpy(this->io->data, fc_buff, result);
            free(fc_buff);
        }
        this->io->data_len += result;
        free(lout);
        // fprintf(stderr, "done\n");
    }

    void encode_flush() override
    {
        if (io->m != coder_io::MENC || flushed)
            return;

        io->meta["magic"] = io->get_magic();
        io->meta["bi"] = (Json::Value::Int)(index);
        io->meta["bn"] = (Json::Value::Int)(num_indexes);
        io->meta["lv"] = (Json::Value::Int)((lzp_valid ? 1 : 0));

        memcpy(io->data + this->io->data_len, indexes, 4 * num_indexes);
        this->io->data_len += 4 * num_indexes;

        flushed = true;
    }

    /* External decompression interface, returns actual decompressed length, exits when encountering split_ch during decompression */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold __attribute__ ((unused)) = false)
    {
        check_exit(io->m != coder_io::MDEC, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "only support block decompress, not support line method"); // Not yet extended to line method
        check_exit(split_ch == UINT8_MAX, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, split_ch);

        // fprintf(stderr, "decompress...");
        // int result = fcinit();
        // check_exit(FC_OK == result, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, result);

        uint8_t *lout;
        lout = static_cast<uint8_t*>(safe_alloc(out_len));
        check_exit(lout, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");

        int result = fc_decode(io->data + 1, lout);
        int lzSize = result;
        check_exit(lzSize > 0 && static_cast<uint32_t>(lzSize) <= out_len, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, lzSize);

        index = io->meta["coder"]["bi"].asInt();
        num_indexes = io->meta["coder"]["bn"].asInt();
        lzp_valid = io->meta["coder"]["lv"].asInt();
        memcpy(indexes, io->data + io->meta["dstlen"].asInt() - (num_indexes * 4), (num_indexes * 4));

        // fprintf(stderr, "done\n");
        // fprintf(stderr, "untransform...");
        result = fc_untransform(lout, lzSize, index, num_indexes, indexes);
        check_exit(result >= FC_OK, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, result);

        if (lzp_valid)
        {
            result = fc_unpreprocess(lout + 1, lout + lzSize, out, ppHashSize, ppMinLen);
            check_exit(result >= 0, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, result);
        } else {
            memcpy(out, lout, lzSize);
            result = lzSize;
        }

        if (lout)
            free(lout);
        io->m = coder_io::MDEC;
        // fprintf(stderr, "done\n");

        return result;
    }

private:
    int ppHashSize;
    int ppMinLen;
    bool flushed;

    // used lzp
    bool lzp_valid;

    // bwt
    int index = FC_BAD_ARGS;
    uint8_t num_indexes;
    int indexes[256];
};
