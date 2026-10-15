/*
 * coder_bwt_cm.h - Header file for pbgz project
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

#ifndef _CODER_BWT_CM_H_
#define _CODER_BWT_CM_H_

/* Context model encoder based on BWT transformation */
#include <stdint.h>
#include <algorithm>
#include <cinttypes>

#include "coder_io.h"
#include "coder.h"
#include "libsais.h"


class coder_bwt_cm: public coder
{
private:
    template <int32_t RATE> /* Smaller RATE = more sensitive probability adaptation */
    class _counter
    {
    private:
        uint16_t c;

    public:
        _counter()
        {
            c = 1 << 15; /* 0.5, max value */
        }

        inline void update_high()
        {
            c += (c ^ 0xFFFF) >> RATE;
        }

        inline void update_low()
        {
            c -= c >> RATE;
        }

        /* Update towards 1 when the bit is set, towards 0 otherwise. */
        inline void update(int32_t bit)
        {
            const uint16_t up = (uint16_t)((c ^ 0xFFFF) >> RATE);
            const uint16_t down = (uint16_t)(c >> RATE);
            c = (uint16_t)(c + (bit ? up : (uint16_t)(0 - down)));
        }

        inline void set(uint16_t count)
        {
            c = count;
        }

        inline uint16_t get() const
        {
            return c;
        }
    };

    class _coder
    {
    public:
        _coder(coder_io *io)
        {
            low = 0;
            high = UINT32_MAX;
            code = 0;
            this->io = io;
            flushed = false;
        }

        virtual ~_coder()
        {
            if (!flushed && io->m == coder_io::MENC)
                encode_flush();
        }

        void encode_flush()
        {
            if (flushed)
                return;
            for (int32_t i = 0; i < 4; ++i)
            {
                /* On overflow, latch IO_BUF_FULL and abandon the rest. */
                if (io->data_len >= io->data_capacity) { io->set_err(coder_io::IO_BUF_FULL); break; }
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
            }
            flushed = true;
        }

        template <int32_t P_LOG>
        inline void encode_bit(int32_t bit, uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG); /* p is the probability of the high-probability symbol */

            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* Normalization */
            while ((low ^ high) < (1 << 24))
            {
                /* Boundary check inside the byte pump; on overflow latch IO_BUF_FULL. */
                if (io->data_len >= io->data_capacity) { io->set_err(coder_io::IO_BUF_FULL); return; }
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
                high = (high << 8) + 255;
            }
        }

        template <int32_t P_LOG>
        inline int32_t decode_bit(uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG);

            const int32_t bit = (code <= mid);
            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* Normalization */
            while ((low ^ high) < (1 << 24))
            {
                /* Read-side exhaustion latches IO_READ_EMPTY; the caller reports it. */
                if (io->data_len >= io->data_capacity) { io->set_err(coder_io::IO_READ_EMPTY); return bit; }
                low <<= 8;
                high = (high << 8) + 255;
                code = (code << 8) + *(io->data + io->data_len++);
            }
            return bit;
        }

        void decode_init()
        {
            int32_t i;
            for (i = 0; i < 4; ++i)
            {
                if (io->data_len >= io->data_capacity) { io->set_err(coder_io::IO_READ_EMPTY); return; }
                code = (code << 8) + *(io->data + io->data_len++);
            }
        }

    private:
        uint32_t low;
        uint32_t high;
        uint32_t code;
        bool flushed;
        coder_io *io;
    };

public:
    coder_bwt_cm(coder_io *io)
    {
        int32_t i, j, k;
        coder_buff = nullptr;
        buf_bwt = nullptr;
        buf_arr = nullptr;
        cap = 0;
        flushed = false;
        run = 0;
        c1 = 0;
        c2 = 0;
        curr_incache = 0;
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_bwt_cm");

        for (i = 0; i < 2; ++i)
        {
            for (j = 0; j < 256; ++j)
            {
                for (k = 0; k <= 16; ++k)
                    counter2[i][j][k].set((k << 12) - (k == 16));
            }
        }

        coder = new _coder(io);
        check_exit(coder, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "coder initialize failed: no enough memory");
    }

    virtual ~coder_bwt_cm()
    {
        if (!flushed && io->m == coder_io::MENC)
            encode_flush();
        if (coder)
            delete coder;
        if (coder_buff)
            safe_free((void**)&coder_buff);
    }

    /* Size the work buffers for at least n bytes, never above bsize. Cached data is preserved. */
    void ensureCapacity(int32_t need)
    {
        if (need <= cap)
            return;
        int32_t want = (cap > 0) ? (cap << 1) : 4096;
        if (want < need)
            want = need;
        want = (want + 7) & ~7; /* keep buf_arr 4-byte aligned */
        if (want > bsize)
            want = bsize; /* the one size the alignment rounding must not push past */

        uint8_t* grown = static_cast<uint8_t*>(safe_alloc((uint64_t)want * 6));
        check_exit(grown, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                   "Error: Insufficient memory: need %" PRIu64 " MB\n",
                   ((uint64_t)want * 6) >> 20);
        if (coder_buff != nullptr)
        {
            if (curr_incache > 0)
                memcpy(grown, coder_buff, curr_incache);
            safe_free((void**)&coder_buff);
        }

        coder_buff = grown;
        buf_bwt = coder_buff + want;
        buf_arr = (uint32_t *)(buf_bwt + want); /* 4 bytes per position */
        cap = want;
    }

    /* External compression interface */
    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        if (io->m != coder_io::MENC)
        { /* Initialize when not initialized */
            const int32_t tab[10] =
                {
                    0,
                    1 << 20,    // -1 - 1 MB
                    1 << 22,    // -2 - 4 MB
                    1 << 23,    // -3 - 8 MB
                    0x00FFFFFF, // -4 - ~16 MB (Default)
                    1 << 25,    // -5 - 32 MB
                    1 << 26,    // -6 - 64 MB
                    1 << 27,    // -7 - 128 MB
                    1 << 28,    // -8 - 256 MB
                    0x7FFFFFFF, // -9 - ~2 GB
                };

            io->m = coder_io::MENC;
            int32_t level = (io->meta["level"].isInt()) ? (io->meta["level"].asInt()) : 4;
            io->meta["level"] = (Json::Value::Int)level;

            check_exit(level <= 9 && level >= 0,  coder_ns::CODER_ERR_BAD_ARGS, "coder level should in [0, 9], current is %d", level);
            const int32_t BLOCK_SIZE = (268435456);

            /* Block size only; the work buffers follow the data (see ensureCapacity). */
            bsize = std::min(tab[level], BLOCK_SIZE);
        }

        /* Feed the whole input in bsize chunks without writing past coder_buff. */
        const auto emitBlock = [&](int32_t n) {
            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, n, 0, nullptr);
            check_exit(idx >= 1,  coder_ns::CODER_ERR_INNER, "coder transform failed: idx %d", idx);

            put32(n); /* Encode current length */
            put32(idx);   /* Encode BWT index */

            for (int32_t k = 0; k < n; k++) /* Encode transformed data */
                put(buf_bwt[k]);
        };

        int32_t left = (int32_t)in_len;
        const uint8_t* src = in;

        /* 1. Top up the partially filled buffer; flush once it reaches bsize. */
        while (curr_incache > 0 && left > 0)
        {
            int32_t take = bsize - curr_incache;
            if (left < take)
                take = left;
            ensureCapacity(curr_incache + take);
            memcpy(coder_buff + curr_incache, src, take);
            curr_incache += take;
            src += take;
            left -= take;
            if (curr_incache == bsize)
            {
                emitBlock(bsize);
                curr_incache = 0;
            }
        }

        /* 2. Compress the bulk of the input in whole bsize blocks. */
        while (left >= bsize)
        {
            ensureCapacity(bsize);
            memcpy(coder_buff, src, bsize);
            emitBlock(bsize);
            src += bsize;
            left -= bsize;
        }

        /* 3. Keep the tail (< bsize) cached for encode_flush(). */
        if (left > 0)
        {
            ensureCapacity(left);
            memcpy(coder_buff, src, left);
            curr_incache = left;
        }
    }

    /* Call this after encoding to compress cached data */
    void encode_flush() override
    {
        int32_t i;
        if (io->m != coder_io::MENC || flushed) {
            return;
        }
        /* First compress all cached data */
        if (curr_incache > 0)
        {
            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, curr_incache, 0, nullptr);
            check_exit(idx >= 1,  coder_ns::CODER_ERR_INNER, "coder transform failed: idx %d", idx);
            put32(curr_incache); /* Encode current length */
            put32(idx);          /* Encode BWT index */

            for (i = 0; i < curr_incache; i++) /* Encode transformed data */
                put(buf_bwt[i]);
        }
        /* Then set EOF marker */
        put32(0);
        coder->encode_flush();
        flushed = true;
    }

    /* Decompress; returns actual length, stops at split_ch. */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold __attribute__ ((unused)) = false) override
    {
        uint8_t ch;
        int32_t len = 0, n;

        /* Error convention: <0 error, ==0 clean EOF, >0 decoded length. */
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;

        if (io->m != coder_io::MDEC)
        {
            decode_init();
            if (io->err != coder_io::IO_OK)
                return coder_ns::CODER_ERR_STREAM_END;
            /* A negative return means a corrupt stream or OOM; propagate it. */
            int32_t ret = decode_one_block();
            if (ret < 0)
                return ret;
            if (ret == 0)
                return 0;
            io->m = coder_io::MDEC;
        }

        if (split_ch != UINT8_MAX)
        {
            for (;;)
            {
                if (curr_out_offset == bsize)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < bsize; n++) {
                    ch = coder_buff[n];
                    *(out + len++) = ch;
                    if (ch == split_ch) {
                        curr_out_offset = ++n;
                        return len;
                    }
                    if (len >= (int32_t)out_len) {
                        /* Output buffer too small. */
                        curr_out_offset = ++n;
                        io->set_err(coder_io::IO_BUF_FULL);
                        return coder_ns::CODER_ERR_BUF_SMALL;
                    }
                }
                curr_out_offset = n;
            }
        } else {
            for (;;)
            {
                if (curr_out_offset == bsize)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < bsize; n++) {
                    *(out + len++) = coder_buff[n];
                    if (len >= (int32_t)out_len) {
                        /* Fixed-length mode: filling out_len is a normal end; the rest is left for the next call. */
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        }

stream_end:
        /* An early stream end means corruption, not a normal end. */
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool) override
    {
        return 0;
    }

    /* Decompression initialization */
    void decode_init()
    {
        bsize = 0;
        coder->decode_init();
    }

    /* Length of currently consumed data. */
    int32_t decode_inlen()
    {
        return this->io->data_len;
    }

private:
    inline void put32(uint32_t x)
    {
        uint32_t i;
        for (i = 1 << 31; i > 0; i >>= 1)
            coder->encode_bit<1>(x & i, 1);
    }

    inline uint32_t get32()
    {
        int32_t i;
        uint32_t x = 0;
        for (i = 0; i < 32; ++i)
            x += x + coder->decode_bit<1>(1);
        return x;
    }

    inline void put(int32_t c)
    {
        const int32_t f = (run > 2); /* run = consecutive repeats of the previous character */

        /* c1/c2 rows and this run's counter2 plane are fixed for the whole byte. */
        _counter<4>* const row1_c1 = counter1[c1];
        _counter<4>* const row1_c2 = counter1[c2];
        _counter<6>* const plane2 = &counter2[f][0][0];

        int32_t ctx = 1;
        for (int32_t i = 128; i > 0; i >>= 1) /* Get all bits of current byte c, starting from high bit */
        {
            const int32_t p0 = counter0[ctx].get();
            const int32_t p1 = row1_c1[ctx].get();
            const int32_t p2 = row1_c2[ctx].get();
            const int32_t p = (((p0 + p1) * 7) + p2 + p2) >> 4;

            /* Get probability of auxiliary symbols through linear interpolation */
            const int32_t j = p >> 12; /* SSE index: highest 4 bits of mixed probability */
            _counter<6>* const row2 = plane2 + (size_t)ctx * 17;
            uint32_t x12;
            memcpy(&x12, (const void*)(row2 + j), sizeof(uint32_t));
            const int32_t x1 = (int32_t)(uint16_t)x12;
            const int32_t x2 = (int32_t)(uint16_t)(x12 >> 16);
            const int32_t ssep = x1 + (((x2 - x1) * (p & 4095)) >> 12);

            /* p + 3*ssep is the probability that the next bit is 1. */
            const int32_t bit = (c & i) ? 1 : 0;
            coder->encode_bit<18>(bit, p + ssep + ssep + ssep);

            counter0[ctx].update(bit);
            row1_c1[ctx].update(bit);
            row2[j].update(bit);
            row2[j + 1].update(bit);

            ctx += ctx + bit; /* ctx left shift by one bit (+1 when the bit was 1) */
        }

        c2 = c1;        /* c2 = previous byte */
        c1 = ctx - 256; /* c1 = current byte */

        if (c1 == c2)   /* Increment run if the byte repeats. */
            ++run;
        else
            run = 0;
    }

    inline int32_t get()
    {
        const int32_t f = (run > 2);

        /* Same factoring as put(). */
        _counter<4>* const row1_c1 = counter1[c1];
        _counter<4>* const row1_c2 = counter1[c2];
        _counter<6>* const plane2 = &counter2[f][0][0];

        int32_t ctx = 1;
        while (ctx < 256)
        {
            const int32_t p0 = counter0[ctx].get();
            const int32_t p1 = row1_c1[ctx].get();
            const int32_t p2 = row1_c2[ctx].get();
            const int32_t p = (((p0 + p1) * 7) + p2 + p2) >> 4;

            /* Get probability of auxiliary symbols through linear interpolation */
            const int32_t j = p >> 12;
            _counter<6>* const row2 = plane2 + (size_t)ctx * 17;
            uint32_t x12;
            memcpy(&x12, (const void*)(row2 + j), sizeof(uint32_t));
            const int32_t x1 = (int32_t)(uint16_t)x12;
            const int32_t x2 = (int32_t)(uint16_t)(x12 >> 16);
            const int32_t ssep = x1 + (((x2 - x1) * (p & 0xFFF)) >> 12);

            /* Same counter updates as put(). */
            const int32_t bit = coder->decode_bit<18>(p + ssep + ssep + ssep);
            counter0[ctx].update(bit);
            row1_c1[ctx].update(bit);
            row2[j].update(bit);
            row2[j + 1].update(bit);

            ctx += ctx + bit; /* ctx left shift by one bit (+1 when the bit was 1) */
        }

        c2 = c1;
        c1 = ctx - 256;
        if (c1 == c2)
            ++run;
        else
            run = 0;
        return c1;
    }

    int32_t decode_one_block()
    {
        int32_t i, p;
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        bsize = get32();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;  // stream truncated while reading bsize
        if (bsize == 0)
            return 0;
        // check_exit(bsize > 0,  coder_ns::CODER_ERR_INNER, "block size %d is invalid", bsize);
        if (!(bsize > 0)) {
            coder_logger(coder_ns::ERROR, "block size %d is invalid", bsize);
            return coder_ns::CODER_ERR_INNER;
        }


        bidx = get32();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        // check_exit(bidx >= 1 && bidx <= bsize,  coder_ns::CODER_ERR_INNER, "corrupt input");
        if (!(bidx >= 1 && bidx <= bsize)) {
            coder_logger(coder_ns::ERROR, "corrupt input");
            return coder_ns::CODER_ERR_INNER;
        }

        if (!buf_arr)
        {   /* Space not allocated, allocate space first */
            if (bsize >= (1 << 24)) /* 5*N */
            {
                i = bsize + (bsize << 2) + bsize;
                coder_buff = static_cast<uint8_t*>(safe_alloc(i));
                // check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Error: Insufficient memory: need %" PRIu64 " MB\n",
                //         (static_cast<uint64_t>(sizeof(uint8_t)) * (i)) >> 20);
                if (!coder_buff) {
                    coder_logger(coder_ns::ERROR, "Error: Insufficient memory: need %" PRIu64 " MB\n",
                        (static_cast<uint64_t>(sizeof(uint8_t)) * (i)) >> 20);
                    return coder_ns::CODER_ERR_MEM_ALLOC_FAIL;
                }
                buf_arr = (uint32_t *)(coder_buff + bsize);
                buf_bwt = coder_buff + bsize + (bsize << 2);
            }
            else
            { /* 4*N */
                i = bsize + (bsize << 2);
                coder_buff = static_cast<uint8_t*>(safe_alloc(i));
                //check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Error: Insufficient memory: need %" PRIu64 " MB\n",
                //    (static_cast<uint64_t>(sizeof(uint8_t)) * (i)) >> 20);
                if (!coder_buff) {
                    coder_logger(coder_ns::ERROR, "Error: Insufficient memory: need %" PRIu64 " MB\n",
                        (static_cast<uint64_t>(sizeof(uint8_t)) * (i)) >> 20);
                    return coder_ns::CODER_ERR_MEM_ALLOC_FAIL;
                }
                buf_arr = (uint32_t *)(coder_buff + bsize);
                buf_bwt = nullptr;
            }
        }

        if (bsize >= (1 << 24)) /* 5*N */
        {
            /* BW inverse transformation */
            memset(cnt, 0, sizeof(cnt));
            for (i = 0; i < bsize; ++i)
                ++cnt[(buf_bwt[i] = get()) + 1];
            for (i = 1; i < 256; ++i)
                cnt[i] += cnt[i - 1];

            for (i = 0; i < bidx; ++i)
                buf_arr[cnt[buf_bwt[i]]++] = i;
            for (i = bidx + 1; i <= bsize; ++i)
                buf_arr[cnt[buf_bwt[i - 1]]++] = i;

            p = bidx;
            for (i = 0; i < bsize; ++i)
            {
                p = buf_arr[p - 1];
                coder_buff[i] = buf_bwt[p - (p >= bidx)];
            }
        }
        else
        {   /* 4*N */
            /* BW inverse transformation */
            memset(cnt, 0, sizeof(cnt));
            for (i = 0; i < bsize; ++i)
                ++cnt[(buf_arr[i] = get()) + 1];
            for (i = 1; i < 256; ++i)
                cnt[i] += cnt[i - 1];

            for (i = 0; i < bidx; ++i)
                buf_arr[cnt[buf_arr[i] & 255]++] |= i << 8;
            for (i = bidx + 1; i <= bsize; ++i)
                buf_arr[cnt[buf_arr[i - 1] & 255]++] |= i << 8;

            p = bidx;
            for (i = 0; i < bsize; ++i)
            {
                p = buf_arr[p - 1] >> 8;
                coder_buff[i] = buf_arr[p - (p >= bidx)];
            }
        }
        curr_out_offset = 0;
        /* Stream exhausted during the inverse transform: report corruption. */
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        return bsize;
    }

private:
    int32_t run;
    int32_t c1;
    int32_t c2;
    _counter<2> counter0[256]; /* counter0 has the greatest impact */
    _counter<4> counter1[256][256];
    _counter<6> counter2[2][256][17];

    uint8_t *coder_buff;  /* Cached original data, transformed once a full block is ready */
    int32_t curr_incache; /* Length of currently cached data */
    uint8_t *buf_bwt;
    uint32_t *buf_arr;
    int32_t cap;          /* Bytes each of coder_buff/buf_bwt currently holds (buf_arr holds 4x) */
    int32_t bsize;
    int32_t bidx;
    int32_t curr_out_offset; /* Offset position during decompression */
    int32_t cnt[257];
    bool flushed;

    _coder *coder;
};
#endif
