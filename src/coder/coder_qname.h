/*
 * coder_qname.h - Codec dedicated to QNAME.
 *
 * Cross-line dedup: a direct-mapped recent-line hash table lets a QNAME that
 * occurred within RING_SIZE lines be coded as a 1-bit marker plus distance.
 * Misses are coded per position (line length + one adaptive model per position).
 * Input is the whole QNAME including the trailing '\t'.
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
#ifndef _CODER_QNAME_H_
#define _CODER_QNAME_H_

#include <stdint.h>
#include <cstring>

#include "simple_model.h"
#include "coder_io.h"
#include "coder.h"

class coder_qname : public coder
{
public:
    static const uint32_t DM_BITS = 18;
    static const uint32_t DM_SIZE = 1u << DM_BITS;
    static const uint32_t RING_SIZE = 4096;
    /* Longest supported QNAME (including '\t'): SAM spec QNAME<=254, +tab=255. */
    static const uint32_t MAX_POS = 256;

    coder_qname(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qname");
        this->flushed = false;
        this->line_idx = 0;
        this->last_len = 0;

        void* ptr = safe_alloc(MAX_POS * sizeof(SIMPLE_MODEL<256>));
        model_pos = static_cast<SIMPLE_MODEL<256>*>(ptr);
        check_exit(model_pos != nullptr, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
            "coder_qname: model_pos alloc failed");
        for (uint32_t n = 0; n < MAX_POS; n++) {
            model_pos[n].reset();
        }
        for (uint32_t n = 0; n < 256; n++) {
            model_len[n].reset();
        }
        model_copy.reset();
        model_dist.reset();
        memset(dm, 0xff, sizeof(dm)); /* line=0xffffffff means empty */
        memset(rlen, 0, sizeof(rlen));
    }

    virtual ~coder_qname()
    {
        if (model_pos) {
            safe_free((void**)&model_pos);
        }
        if (!flushed && io->m == coder_io::MENC) {
            encode_flush();
        }
    }

    /* Encode one line of QNAME (including the trailing '\t'). */
    void encode_line(const uint8_t *in, const uint32_t in_len, [[maybe_unused]] bool need2hold = false)
    {
        check_exit(in_len <= MAX_POS, coder_ns::CODER_ERR_BAD_ARGS,
            "coder_qname: QNAME too long (%u)", in_len);

        if (io->m != coder_io::MENC) {
            rc.output((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartEncode();
            io->m = coder_io::MENC;
        }

        uint32_t h = hash(in, in_len);
        const uint32_t ridx = line_idx % RING_SIZE;
        const uint32_t slot = h & (DM_SIZE - 1);

        /* Probe for the most recent occurrence: direct-mapped slot hit + line within window + identical content -> COPY. */
        int32_t dist = -1;
        if (dm[slot].hash == h && dm[slot].line != UINT32_MAX) {
            uint32_t lastLine = dm[slot].line;
            if (lastLine < line_idx && line_idx - lastLine <= RING_SIZE - 1 &&
                rline[lastLine % RING_SIZE] == lastLine &&
                rlen[lastLine % RING_SIZE] == in_len &&
                memcmp(recent[lastLine % RING_SIZE], in, in_len) == 0) {
                dist = (int32_t)(line_idx - lastLine);
            }
        }

        if (dist >= 0) {
            model_copy.encodeSymbol(&rc, 1);
            model_dist.encodeSymbol(&rc, (uint16_t)dist);
        } else {
            model_copy.encodeSymbol(&rc, 0);
            model_len[last_len].encodeSymbolOrder(&rc, in_len);
            for (uint32_t p = 0; p < in_len; p++) {
                model_pos[p].encodeSymbol(&rc, in[p]);
            }
        }

        last_len = (int32_t)in_len;

        /* Update the ring buffer and the direct-mapping table. */
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)in_len;
        memcpy(recent[ridx], in, in_len);

        dm[slot].hash = h;
        dm[slot].line = line_idx;

        line_idx++;
    }

    /* Decode one line of QNAME (including the trailing '\t'), returning the number of bytes written to out. */
    int32_t decode_line(uint8_t *out, uint32_t out_len,
                        [[maybe_unused]] uint8_t split_ch = UINT8_MAX,
                        [[maybe_unused]] bool need2hold = false)
    {
        if (io->m != coder_io::MDEC) {
            rc.input((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartDecode();
            io->m = coder_io::MDEC;
        }

        uint16_t copy = model_copy.decodeSymbol(&rc);
        uint32_t len;
        if (copy) {
            uint32_t dist = model_dist.decodeSymbol(&rc);
            /* dist can be line_idx (copying line 0), but must be >= 1 and <= line_idx. */
            if (dist == 0 || dist > line_idx) {
                return coder_ns::CODER_ERR_BAD_ARGS;
            }
            uint32_t ridx = (line_idx - dist) % RING_SIZE;
            len = rlen[ridx];
            if (len > out_len) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            memcpy(out, recent[ridx], len);
        } else {
            len = model_len[last_len].decodeSymbolOrder(&rc);
            if (len > out_len || len > MAX_POS) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            for (uint32_t p = 0; p < len; p++) {
                out[p] = (uint8_t)model_pos[p].decodeSymbol(&rc);
            }
        }
        last_len = (int32_t)len;

        /* Update the ring buffer, as the encoder does. */
        uint32_t ridx = line_idx % RING_SIZE;
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)len;
        memcpy(recent[ridx], out, len);

        line_idx++;
        return (int32_t)len;
    }

    /* fake */
    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) { return 0; }

    void encode_flush()
    {
        if (flushed) {
            return;
        }
        if (io->m == coder_io::MENC) {
            if (rc.FinishEncode() < 0) {
                io->set_err(coder_io::IO_BUF_FULL);
            }
            io->data_len += rc.size_out();
        }
        flushed = true;
    }

private:
    static uint32_t hash(const uint8_t *s, uint32_t n)
    {
        uint32_t h = 2166136261u;
        for (uint32_t i = 0; i < n; i++) {
            h ^= s[i];
            h *= 16777619u;
        }
        return h;
    }

    RangeCoder rc;
    SIMPLE_MODEL<2> model_copy;
    SIMPLE_MODEL<4096> model_dist;
    SIMPLE_MODEL<256> model_len[256];
    SIMPLE_MODEL<256> *model_pos;

    struct HashSlot {
        uint32_t hash;
        uint32_t line;
    } dm[DM_SIZE];
    uint8_t rlen[RING_SIZE];
    uint32_t rline[RING_SIZE];
    uint8_t recent[RING_SIZE][MAX_POS];

    uint32_t line_idx;
    int32_t last_len;
    bool flushed;
};

#endif
