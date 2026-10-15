/*
 * coder_arith_k.h - lightweight entropy coders for the FC/base streams.
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
#include "base_pack.h"

#include <cstdlib>
#include <cstring>
#include <vector>

extern "C" {
unsigned char *arith_compress(unsigned char *in, unsigned int in_size,
                              unsigned int *out_size, int order);
unsigned char *arith_uncompress_to(unsigned char *in, unsigned int in_size,
                                   unsigned char *out, unsigned int *out_sz);
unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
                             unsigned int *out_size, int order);
unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
                               unsigned int *out_size);
}

class coder_arith_k : public coder
{
public:
    coder_arith_k(coder_io *io, int order = 1, bool mapBases = true, bool useRans = false,
                  bool ransPack = true)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_arith_k");
        this->order = order;
        this->mapBases = mapBases;
        this->useRans = useRans;
        this->ransPack = ransPack;
        heldSrc = nullptr;
        heldLen = 0;
        nCalls = 0;
        packedThis = false;
        flushed = false;
    }

    virtual ~coder_arith_k()
    {
        if (io->m == coder_io::MENC && !flushed)
            encode_flush();
    }

    bool supportsLineMode() const override { return false; }

    void encode_line(const uint8_t *in, const uint32_t in_len,
                     bool need2hold __attribute__((unused)) = false) override
    {
        if (io->m != coder_io::MENC) {
            io->m = coder_io::MENC;
        }
        if (in_len == 0) {
            return;
        }
        if (nCalls == 1) {
            /* promote the held single-call buffer before prepare() reuses it */
            buff.insert(buff.end(), heldSrc, heldSrc + heldLen);
        }
        const uint8_t *src = in;
        uint32_t n = in_len;
        prepare(in, in_len, src, n);
        if (nCalls == 0) {
            /* single-call (the SEQ match stream): hold it, no copy into buff */
            heldSrc = src;
            heldLen = n;
        } else {
            buff.insert(buff.end(), src, src + n);
        }
        ++nCalls;
    }

    void encode_flush() override
    {
        if (io->m != coder_io::MENC || flushed) {
            return;
        }
        const uint8_t *data = (nCalls == 1) ? heldSrc : buff.data();
        const uint32_t len = (nCalls == 1) ? heldLen : (uint32_t)buff.size();
        if (len > 0) {
            unsigned int outsz = 0;
            unsigned char *cout = useRans
                ? rans_compress(const_cast<uint8_t *>(data), len, &outsz, order)
                : arith_compress(const_cast<uint8_t *>(data), len, &outsz, order);
            check_exit(cout != nullptr, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_arith_k: compress failed");
            if ((int32_t)outsz > (io->data_capacity - io->data_len)) {
                free(cout);
                check_exit(false, coder_ns::CODER_ERR_INNER,
                           "coder_arith_k: output overflow (need %u, have %d)", outsz,
                           io->data_capacity - io->data_len);
            }
            memcpy(io->data + io->data_len, cout, outsz);
            io->data_len += (int32_t)outsz;
            free(cout);
            io->meta["akorder"] = (Json::Value::Int)order;
            io->meta["akmap"] = (Json::Value::Int)(packedThis ? 1 : 0);
            io->meta["akrans"] = (Json::Value::Int)(useRans ? 1 : 0);
            io->meta["akpack"] = (Json::Value::Int)((useRans && packedThis && ransPack) ? 1 : 0);
        }
        flushed = true;
    }

    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX,
                        bool need2hold __attribute__((unused)) = false) override
    {
        check_exit(split_ch == UINT8_MAX, coder_ns::CODER_ERR_INNER,
                   "coder_arith_k: split not supported");
        const unsigned int in_size =
            (io->meta["dstlen"].isUInt()) ? io->meta["dstlen"].asUInt() : (unsigned int)io->data_len;
        const bool map = (io->meta["coder"]["akmap"].asInt() == 1);
        const bool rans = (io->meta["coder"]["akrans"].asInt() == 1);

        uint8_t *dst = out;
        unsigned int want = (unsigned int)out_len;
        if (rans) {
            /* rANS decodes either the packed bytes (ACGT) or the raw bytes. */
            unsigned int osz = 0;
            unsigned char *r = rans_uncompress(io->data, in_size, &osz);
            check_exit(r != nullptr, coder_ns::CODER_ERR_INNER, "coder_arith_k: rans decode failed");
            if (map && io->meta["coder"]["akpack"].asInt() == 1) {
                for (unsigned int i = 0; i < out_len; ++i) {
                    const uint8_t b = r[i / 4];
                    out[i] = (uint8_t)((b >> (2 * (i % 4))) & 3u);
                }
            } else {
                memcpy(out, r, out_len);
            }
            free(r);
        } else {
            unsigned int osz = want;
            unsigned char *r = arith_uncompress_to(io->data, in_size, dst, &osz);
            check_exit(r == dst && (int32_t)osz == (int32_t)out_len, coder_ns::CODER_ERR_INNER,
                       "coder_arith_k: decompress failed (%d vs %u)", (int)osz, out_len);
        }
        if (map) {
            static const uint8_t kBase[4] = {'A', 'C', 'G', 'T'};
            for (uint32_t i = 0; i < (uint32_t)out_len; ++i) {
                out[i] = kBase[out[i] & 3u];
            }
        }
        io->m = coder_io::MDEC;
        return (int32_t)out_len;
    }

    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) override { return 0; }

private:
    /* Map/pack the input for the entropy stage and return the buffer to code. */
    void prepare(const uint8_t *in, uint32_t in_len, const uint8_t *&src, uint32_t &n)
    {
        src = in;
        n = in_len;
        packedThis = false;
        if (mapBases && useRans && ransPack) {
            /* One fused pass: ACGT check + map + pack four bases per byte. */
            const uint32_t pn = (in_len + 3) / 4;
            bufPacked.resize(pn);
            if (pbgz_base_pack::packBases(in, in_len, bufPacked.data())) {
                src = bufPacked.data();
                n = pn;
                packedThis = true;
            }
        } else if (mapBases) {
            /* arith path: map to 0..3 only (its own X_PACK packs the bits). */
            mapped.resize(in_len);
            bool acgt = true;
            for (uint32_t i = 0; i < in_len; ++i) {
                const uint8_t c = pbgz_base_pack::baseCode(in[i]);
                if (c > 3) { acgt = false; }
                mapped[i] = c;
            }
            if (acgt) {
                src = mapped.data();
                n = in_len;
                packedThis = true;
            }
        }
    }

    int order;
    bool mapBases;
    bool useRans;
    bool ransPack;
    std::vector<uint8_t> buff;      /* multi-call accumulation (unused by the SEQ path) */
    std::vector<uint8_t> mapped;    /* arith path scratch */
    std::vector<uint8_t> bufPacked; /* fused pack output */
    const uint8_t *heldSrc;         /* single-call prepared buffer (no copy into buff) */
    uint32_t heldLen;
    int nCalls;
    bool packedThis;
    bool flushed;
};
