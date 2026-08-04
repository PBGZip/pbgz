/*
 * coder_qual.h - Header file for pbgz project
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

#ifndef _CODER_QUAL_H_
#define _CODER_QUAL_H_

#include <math.h>
#include "qual_model.h"
#include "coder_io.h"

#define ROUND_UP(x) (x == 0 ? 0 : 1 << ((x)-1))
/* Quantisize quality scores q' = ((q + ROUND_UP(QUANT_Q)) >> QUANT_Q) */
#define QUANT_Q 1
/* Keep as a power of 2 */
/*QMAX = size of the quality score's alphabet*/
#define Q_LOG (7 - QUANT_Q)

/*AMAX = size of the base call's alphabet*/
#define AMAX 5 /* use N */
#define A_LOG 3

class coder_qual
{
public:
    coder_qual(coder_io *io, bool is_gen2, const std::vector<std::pair<uint16_t, uint16_t>> &qual_freq_table)
    {
        flushed = false;
        this->is_gen2 = is_gen2;

        if (is_gen2) {
            s_ctx_len = 1;
            q_ctx_len = 2;
        } else {
            s_ctx_len = 3;
            q_ctx_len = 2;
        }

        Q_MASK = (q_ctx_len > 0) ? ((1 << (Q_LOG * q_ctx_len)) - 1) : 0;
        S_CTX = pow(AMAX, s_ctx_len);

        if (is_gen2) {
            SSE_CTX = 1 << 4;
            ctx_cnt = S_CTX * (1 << (Q_LOG * q_ctx_len)) * SSE_CTX;
        } else
            ctx_cnt = S_CTX * (1 << (Q_LOG * q_ctx_len));

        model_qual = new QUAL_MODEL_ENGINE(qual_freq_table, ctx_cnt);
        check_exit(model_qual, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "initialize coder for quality failed: no enough memory.");

        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qual");

        L[static_cast<unsigned char>('A')] = L[static_cast<unsigned char>('a')] = 0;
        L[static_cast<unsigned char>('C')] = L[static_cast<unsigned char>('c')] = 1;
        L[static_cast<unsigned char>('G')] = L[static_cast<unsigned char>('g')] = 2;
        L[static_cast<unsigned char>('T')] = L[static_cast<unsigned char>('t')] = 3;
        L[static_cast<unsigned char>('N')] = L[static_cast<unsigned char>('n')] = 4;
    }

    ~coder_qual()
    {
        if (!flushed)
            if (io->m == coder_io::MENC)
                encode_flush();
        if (model_qual)
            delete model_qual;
    }

    uint32_t get_context_gen2(uint8_t s, uint8_t q, uint32_t &s_prev_ctx, uint32_t &q_prev_ctx, uint32_t sse_ctx)
    {

        /* Get quantized value of q, quantization method
        * QUANT_Q == 1: (q + 1) / 2 , this quantization parameter is used here
        * QUANT_Q == 2: (q + 2) / 4
        * QUANT_Q == 3: (q + 4) / 8
        * QUANT_Q == 4: (q + 8) / 16
        * */
        //q = (q + ROUND_UP(QUANT_Q)) >> QUANT_Q;
        q = Q_QUANT_GEN2[q];

        q_prev_ctx = ((q_prev_ctx << Q_LOG) + q) & Q_MASK;
        s_prev_ctx = ((s_prev_ctx * AMAX) + L[s]) % S_CTX;
        return (((s_prev_ctx << (q_ctx_len * Q_LOG)) + q_prev_ctx ) << 4) | (sse_ctx);
    }

    uint32_t get_context_gen3(uint8_t s, uint8_t q, uint32_t &s_prev_ctx, uint32_t &q_prev_ctx, uint32_t __attribute__ ((unused)) sse_ctx)
    {
        /* Get quantized value of q, quantization method
        * QUANT_Q == 1: (q + 1) / 2 , this quantization parameter is used here
        * QUANT_Q == 2: (q + 2) / 4
        * QUANT_Q == 3: (q + 4) / 8
        * QUANT_Q == 4: (q + 8) / 16
        * */
        q = (q + ROUND_UP(QUANT_Q)) >> QUANT_Q;

        q_prev_ctx = ((q_prev_ctx << Q_LOG) + q) & Q_MASK;
        s_prev_ctx = ((s_prev_ctx * AMAX) + L[s]) % S_CTX;
        return (s_prev_ctx << (q_ctx_len * Q_LOG)) + q_prev_ctx;
    }

    void encode_qual_gen2(uint8_t *seq, uint8_t *qual, uint32_t len) {
        uint32_t i, next_s = 1 + s_ctx_len/2;
        uint32_t s_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;
        uint64_t tot_qual = 0; /* Sum of quality scores in current row */
        uint32_t sse_ctx_cur = 0;

        if (io->m != coder_io::MENC) {
            rc.output((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartEncode();
            io->m = coder_io::MENC;
        }

        /* get first context */
        for (i = 0; i < next_s; i++) {
            if (i < len)
                ctx = get_context_gen2(seq[i], 0, s_prev_ctx, q_prev_ctx, 0);
            else
                ctx = get_context_gen2('A', 0, s_prev_ctx, q_prev_ctx, 0);
        }

        for (i = 0; i < len; i++, next_s++) {
            uint8_t q = (qual[i] - '!');
            model_qual->encodeSymbol(&rc, q, ctx);

            // 1) sse ctx 1, initializing context at start has improvement, but used for all contexts
            // sse_ctx_cur = ((i+1) < qual_len) ? (this->base_match_vec[cur_len][i+1]) : 0;

            // 2) sse ctx 2, take average of all quality rows before current symbol to be compressed, then quantize, obvious improvement
            tot_qual += qual[i];
            sse_ctx_cur = tot_qual / (i + 1);

            // 3) Based on 2, shifting sse_ctr_cur << 4 can slightly improve by 0.1% (34.906%->34.798%)
            // sse_ctx_cur = sse_ctx_cur & 0x3FFF; // quality : 33928344/97500000 . ratio : 34.798%

            if (next_s < len)
                ctx = get_context_gen2(seq[next_s], q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
            else
                ctx = get_context_gen2('A', q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
        }
    }

    void encode_qual_gen3(uint8_t *seq, uint8_t *qual, uint32_t len) {
        // int32_t delta = 5;
        uint32_t i, next_s = 1 + s_ctx_len / 2;
        uint32_t s_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;
        // uint64_t tot_qual = 0; // Sum of quality scores in current row
        uint32_t sse_ctx_cur = 0;

        // Get first context
        for (i = 0; i < next_s; i++) {
            if (i < len)
                ctx = get_context_gen3(seq[i], 0, s_prev_ctx, q_prev_ctx, 0);
            else
                ctx = get_context_gen3('A', 0, s_prev_ctx, q_prev_ctx, 0);
        }

        for (i = 0; i < len; i++, next_s++) {
            uint8_t q = (qual[i] - '!');
            model_qual->encodeSymbol(&rc, q, ctx);

            // 1) sse ctx 1, initializing context at start has improvement, but used for all contexts
            // sse_ctx_cur = ((i+1) < qual_len) ? (this->base_match_vec[cur_len][i+1]) : 0;

            // 2) sse ctx 2, take average of all quality rows before current symbol to be compressed, then quantize, obvious improvement
            // tot_qual += qual[i];
            // sse_ctx_cur = tot_qual / (i + 1);

            // 3) Based on 2, shifting sse_ctr_cur << 4 can slightly improve by 0.1% (34.906%->34.798%)
            // sse_ctx_cur = sse_ctx_cur & 0x3FFF; // quality : 33928344/97500000 . ratio : 34.798%

            if (next_s < len)
                ctx = get_context_gen3(seq[next_s], q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
            else
                ctx = get_context_gen3('A', q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
        }
    }

    void decode_qual_gen2(uint8_t *seq, uint8_t *qual, uint32_t len) {
        // uint32_t last = 0;
        uint32_t i;
        // int32_t delta = 5;
        //int32_t q1 = 0, q2 = 0;
        uint32_t next_s = 1 + s_ctx_len/2;
        uint32_t s_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;
        uint64_t tot_qual = 0;

        if (io->m != coder_io::MDEC)
        {
            rc.input((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartDecode();
            io->m = coder_io::MDEC;
        }

        for (i = 0; i < next_s; i++) {
            if (i < len)
                ctx = get_context_gen2(seq[i], 0, s_prev_ctx, q_prev_ctx, 0);
            else
                ctx = get_context_gen2('A', 0, s_prev_ctx, q_prev_ctx, 0);
        }

        for (i = 0; i < len; i++, next_s++) {
            uint8_t q = (uint8_t) (model_qual->decodeSymbol(&rc, ctx));

            uint32_t sse_ctx_cur = 0;

            // 1) sse ctx 1, initializing context at start has improvement, but used for all contexts
            // sse_ctx_cur = ((i+1) < qual_len) ? (this->base_match_vec[cur_len][i+1]) : 0;

            // 2) sse ctx 2, take average of all quality rows before current symbol to be compressed, then quantize, obvious improvement
            tot_qual += q + '!';
            sse_ctx_cur = tot_qual / (i + 1);

            // sse_ctx_cur = sse_ctx_cur & 0x3FFF;

            if (next_s < len)
                ctx = get_context_gen2(seq[next_s], q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
            else
                ctx = get_context_gen2('A', q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);

            qual[i] = q + '!';

            // q2 = q1;
            // q1 = q;

        }

        if (rc.FinishDecode() < 0) {
            io->err = coder_io::IO_READ_EMPTY;
        }
    }

    void decode_qual_gen3(uint8_t *seq, uint8_t *qual, uint32_t len) {
        // uint32_t last = 0;
        uint32_t i;
        //  int32_t delta = 5;
        // int32_t q1 = 0, q2 = 0;
        uint32_t next_s = 1 + s_ctx_len / 2;
        uint32_t s_prev_ctx = 0, q_prev_ctx = 0, ctx = 0;
        // uint64_t tot_qual = 0;
        uint32_t sse_ctx_cur = 0;

        for (i = 0; i < next_s; i++) {
            if (i < len)
                ctx = get_context_gen3(seq[i], 0, s_prev_ctx, q_prev_ctx, 0);
            else
                ctx = get_context_gen3('A', 0, s_prev_ctx, q_prev_ctx, 0);
        }

        for (i = 0; i < len; i++, next_s++) {
            uint8_t q = (uint8_t) (model_qual->decodeSymbol(&rc, ctx));

            // 1) sse ctx 1, initializing context at start has improvement, but used for all contexts
            // sse_ctx_cur = ((i+1) < qual_len) ? (this->base_match_vec[cur_len][i+1]) : 0;

            // 2) sse ctx 2, take average of all quality rows before current symbol to be compressed, then quantize, obvious improvement
            // tot_qual += q + '!';
            // sse_ctx_cur = tot_qual / (i + 1);

            // sse_ctx_cur = sse_ctx_cur & 0x3FFF;

            if (next_s < len)
                ctx = get_context_gen3(seq[next_s], q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);
            else
                ctx = get_context_gen3('A', q, s_prev_ctx, q_prev_ctx, sse_ctx_cur);

            qual[i] = q + '!';

            // q2 = q1;
            // q1 = q;

        }

        if (rc.FinishDecode() < 0) {
            io->err = coder_io::IO_READ_EMPTY;
        }
    }

    /* Encode unprocessed data */
    void encode_flush()
    {
        if (flushed)
            return;
        if (io->m == coder_io::MENC)
        {
            if (rc.FinishEncode() < 0) {
                io->err = coder_io::IO_BUF_FULL;
            }
            io->data_len += rc.size_out();
        }
        flushed = true;
    }

private:
    QUAL_MODEL_ENGINE *model_qual;
    bool is_gen2;
    int32_t q_ctx_len, s_ctx_len, ctx_cnt, Q_MASK, S_CTX;
    int32_t SSE_CTX; /* Second generation quality score compression extended context */
    coder_io *io;
    RangeCoder rc;
    bool flushed;

private:

     /* sequence table lookups ACGTN->0..4 */
    int32_t L[256];

    /*  generate Q_QUANT_GEN2 */
    // uint8_t Q_QUANT[128] = {0};
    // for (uint8_t i = 33; i < 128; i++) {

    //     uint8_t qq = i - '!';

    //     int32_t high = qq >> 1;
    //     int32_t low = qq % 10;
    //     if(low <= 3)qq = high + 1;
    //     else if (low > 3 && low <= 6) qq = high + 2;
    //     else qq = high + 3;

    //     Q_QUANT[ i - '!'] = qq;
    // }
    // fprintf(stderr, "\n");
    // for (uint8_t i = 0; i < 128; i++) {
    //     fprintf(stderr, "%3u, ", Q_QUANT[i]);
    //     if (( i + 1) % 16 == 0)
    //     fprintf(stderr, "\n");
    // }
    const uint8_t Q_QUANT_GEN2[128] = {
        1, 1, 2, 2, 4, 4, 5, 6, 7, 7, 6, 6, 7, 7, 9, 9,
        10, 11, 12, 12, 11, 11, 12, 12, 14, 14, 15, 16, 17, 17, 16, 16,
        17, 17, 19, 19, 20, 21, 22, 22, 21, 21, 22, 22, 24, 24, 25, 26,
        27, 27, 26, 26, 27, 27, 29, 29, 30, 31, 32, 32, 31, 31, 32, 32,
        34, 34, 35, 36, 37, 37, 36, 36, 37, 37, 39, 39, 40, 41, 42, 42,
        41, 41, 42, 42, 44, 44, 45, 46, 47, 47, 46, 46, 47, 47, 49, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
};

#endif
