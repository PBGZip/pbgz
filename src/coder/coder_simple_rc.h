/*
 * coder_simple_rc.h - Header file for pbgz project
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

#ifndef _CODER_SIMPLE_RC_H_
#define _CODER_SIMPLE_RC_H_

/*
 * coder_simple_rc: BWT + MTF + adaptive bit-level range coder
 *
 * Pipeline: BWT → MTF → simple adaptive RC
 * Optimized for quality score data where MTF output is highly skewed
 * (40%+ symbol 0). Outperforms the context-based RC in coder_fc for
 * such distributions due to lower model overhead.
 *
 * No LZP preprocessing (ineffective for quality data).
 */

#include <stdint.h>
#include <string.h>
#include "coder_io.h"
#include "coder.h"
#include "fc/transform.h"
#include "fc/rangecoder.h"

#define SIMPLE_RC_MIN_LEN 32
#define SIMPLE_RC_MAX_LEN 2146435072

class coder_simple_rc : public coder
{
public:
    coder_simple_rc(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_simple_rc");
        flushed = false;
        index = 0;
        num_indexes = 0;
    }

    virtual ~coder_simple_rc()
    {
        if (io->m == coder_io::MENC && !flushed)
            encode_flush();
    }

    void encode_line(const uint8_t *in, const uint32_t in_len, bool need2hold __attribute__((unused)) = false)
    {
        check_exit(io->m != coder_io::MENC, coder_ns::CODER_ERR_INNER, "only support block compress, not support line method");
        check_exit(in_len > SIMPLE_RC_MIN_LEN && in_len < SIMPLE_RC_MAX_LEN, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, in_len);
        io->m = coder_io::MENC;

        // Step 1: BWT transform (in-place)
        uint8_t *bwt_buf = static_cast<uint8_t*>(safe_alloc(in_len));
        check_exit(bwt_buf, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");
        memcpy(bwt_buf, in, in_len);

        unsigned char num_idx_byte;
        index = fc_transform(bwt_buf, (int)in_len, &num_idx_byte, indexes);
        num_indexes = num_idx_byte;
        if (in_len < 64 * 1024)
            num_indexes = 0;
        check_exit(index >= FC_OK, coder_ns::CODER_ERR_INNER, "BWT failed (%d) : %d", __LINE__, index);

        // Step 2: MTF transform (standard, easily invertible)
        uint8_t *mtf_buf = static_cast<uint8_t*>(safe_alloc(in_len));
        check_exit(mtf_buf, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");
        mtf_encode(bwt_buf, mtf_buf, (int)in_len);

        // Step 3: Simple adaptive RC encode
        int out_cap = (int)in_len * 2 + 4096;
        uint8_t *rc_buf = static_cast<uint8_t*>(safe_alloc(out_cap));
        check_exit(rc_buf, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");

        int rc_size = adaptive_rc_encode(mtf_buf, (int)in_len, rc_buf, out_cap);
        check_exit(rc_size > 0 && rc_size < out_cap, coder_ns::CODER_ERR_INNER, "RC encode failed (%d) : %d", __LINE__, rc_size);

        // Write to output
        check_exit((uint32_t)rc_size <= (uint32_t)(io->data_capacity - io->data_len), coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Output buffer too small");
        memcpy(io->data + io->data_len, rc_buf, rc_size);
        io->data_len += rc_size;

        free(bwt_buf);
        free(mtf_buf);
        free(rc_buf);
    }

    void encode_flush()
    {
        if (io->m != coder_io::MENC || flushed)
            return;

        io->meta["magic"] = io->get_magic();
        io->meta["bi"] = (Json::Value::Int)(index);
        io->meta["bn"] = (Json::Value::Int)(num_indexes);

        // Append BWT indexes after compressed data
        if (num_indexes > 0) {
            memcpy(io->data + io->data_len, indexes, 4 * num_indexes);
            io->data_len += 4 * num_indexes;
        }

        flushed = true;
    }

    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold __attribute__((unused)) = false)
    {
        check_exit(io->m != coder_io::MDEC, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "only support block decompress, not support line method");
        check_exit(split_ch == UINT8_MAX, coder_ns::CODER_ERR_INNER, "check failed (%d) : %d", __LINE__, split_ch);

        // Read metadata
        index = io->meta["coder"]["bi"].asInt();
        num_indexes = io->meta["coder"]["bn"].asInt();
        if (num_indexes > 0) {
            memcpy(indexes, io->data + io->meta["dstlen"].asInt() - (num_indexes * 4), (num_indexes * 4));
        }

        // Step 1: Simple adaptive RC decode
        uint8_t *mtf_buf = static_cast<uint8_t*>(safe_alloc(out_len));
        check_exit(mtf_buf, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");

        int decoded_len = adaptive_rc_decode(io->data, mtf_buf, (int)out_len);
        check_exit(decoded_len == (int)out_len, coder_ns::CODER_ERR_INNER, "RC decode size mismatch (%d) : got %d expect %d", __LINE__, decoded_len, out_len);

        // Step 2: Inverse MTF
        uint8_t *bwt_buf = static_cast<uint8_t*>(safe_alloc(out_len));
        check_exit(bwt_buf, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Memory not enough");
        mtf_decode(mtf_buf, bwt_buf, (int)out_len);

        // Step 3: Inverse BWT
        int result = fc_untransform(bwt_buf, (int)out_len, index, num_indexes, indexes);
        check_exit(result >= FC_OK, coder_ns::CODER_ERR_INNER, "Inverse BWT failed (%d) : %d", __LINE__, result);

        memcpy(out, bwt_buf, out_len);

        free(mtf_buf);
        free(bwt_buf);

        return (int32_t)out_len;
    }

private:
    bool flushed;
    int index;
    int num_indexes;
    int indexes[65536];

    // Standard MTF encode: for each byte, output its position in table, move to front
    static void mtf_encode(const uint8_t *in, uint8_t *out, int n)
    {
        uint8_t table[256];
        for (int i = 0; i < 256; i++) table[i] = (uint8_t)i;

        for (int i = 0; i < n; i++) {
            uint8_t sym = in[i];
            // Find position
            int pos = 0;
            while (table[pos] != sym) pos++;
            out[i] = (uint8_t)pos;
            // Move to front
            if (pos > 0) {
                memmove(table + 1, table, pos);
                table[0] = sym;
            }
        }
    }

    // Standard MTF decode: for each position, output table[pos], move to front
    static void mtf_decode(const uint8_t *in, uint8_t *out, int n)
    {
        uint8_t table[256];
        for (int i = 0; i < 256; i++) table[i] = (uint8_t)i;

        for (int i = 0; i < n; i++) {
            int pos = in[i];
            uint8_t sym = table[pos];
            out[i] = sym;
            // Move to front
            if (pos > 0) {
                memmove(table + 1, table, pos);
                table[0] = sym;
            }
        }
    }

    // Adaptive bit-level range coder encode
    // Each bit position has an independent adaptive probability
    static int adaptive_rc_encode(const uint8_t *data, int n, uint8_t *out, int out_cap)
    {
        int prob[8];
        for (int i = 0; i < 8; i++) prob[i] = 2048; // P(bit=0) in [0, 4096]

        RangeCoder coder;
        coder.InitEncoder(out, out_cap);
        coder.EncodeWord((unsigned int)n);

        for (int i = 0; i < n; i++) {
            uint8_t sym = data[i];
            for (int bit = 7; bit >= 0; bit--) {
                int b = (sym >> bit) & 1;
                int p = prob[bit];
                if (b == 0) {
                    coder.EncodeBit0<12>(p);
                    prob[bit] += (4096 - p) >> 5;
                } else {
                    coder.EncodeBit1<12>(p);
                    prob[bit] -= p >> 5;
                }
            }
        }

        return coder.FinishEncoder();
    }

    // Adaptive bit-level range coder decode (mirror of encode)
    static int adaptive_rc_decode(const uint8_t *comp, uint8_t *out, int out_cap)
    {
        int prob[8];
        for (int i = 0; i < 8; i++) prob[i] = 2048;

        RangeCoder coder;
        coder.InitDecoder(comp);
        int n = (int)coder.DecodeWord();
        if (n > out_cap) return -1;

        for (int i = 0; i < n; i++) {
            uint8_t sym = 0;
            for (int bit = 7; bit >= 0; bit--) {
                int p = prob[bit];
                int b = coder.DecodeBit<12>(p);
                if (b == 0) {
                    prob[bit] += (4096 - p) >> 5;
                } else {
                    sym |= (1 << bit);
                    prob[bit] -= p >> 5;
                }
            }
            out[i] = sym;
        }

        return n;
    }
};

#endif // _CODER_SIMPLE_RC_H_
