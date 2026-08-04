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
#include "bcm_libsais.h"


class coder_bwt_cm: public coder
{
private:
    template <int32_t RATE> /* Use count to represent probability, smaller RATE value means more sensitive and greater impact, i.e., high probability when each symbol arrives */
    class _counter
    {
    private:
        uint16_t c;

    public:
        _counter()
        {
            c = 1 << 15; /* 0.5, 32768, left shift by one (multiply by 2) is the maximum value */
        }

        inline void update_high()
        {
            c += (c ^ 0xFFFF) >> RATE;
        }

        inline void update_low()
        {
            c -= c >> RATE;
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
                /* 见 encode_bit 的边界说明：写满即 latch，剩余字节放弃，由出口统一上报。 */
                if (io->data_len >= io->data_capacity) { io->err = coder_io::IO_BUF_FULL; break; }
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
            }
            flushed = true;
        }

        template <int32_t P_LOG>
        inline void encode_bit(int32_t bit, uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG); /* High probability symbols are placed first, low probability symbols last, P represents the probability of high probability symbols */

            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* Normalization */
            while ((low ^ high) < (1 << 24))
            {
                /*
                 * 边界检查沉在字节泵里，fqz/arith_dynamic 同款（RC_ShiftLowCheck）：
                 * 只在真正吐字节的归一化分支查一次，单次寄存器比较，符号层零开销。
                 * 写满则 latch io->err = IO_BUF_FULL 并放弃本次写出；此后编码结果已
                 * 无意义，由 encode_flush 出口统一查 err 报错，热循环不做逐符号判断。
                 */
                if (io->data_len >= io->data_capacity) { io->err = coder_io::IO_BUF_FULL; return; }
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
                /*
                 * 读端耗尽置 IO_READ_EMPTY 并返回：原来会越过码流末尾读到相邻字段的
                 * 数据，产出错误却"成功"的内容。latch 后后续 decode_bit 空转，
                 * 出口 decode_line 查 err 统一报错（同 RC_Decode 的 in_end 检查）。
                 */
                if (io->data_len >= io->data_capacity) { io->err = coder_io::IO_READ_EMPTY; return bit; }
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
                if (io->data_len >= io->data_capacity) { io->err = coder_io::IO_READ_EMPTY; return; }
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

    /* External compression interface */
    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        int32_t i;
        if (io->m != coder_io::MENC)
        { /* Initialize when not initialized */
            int32_t len_buf, len_bwt, len_arr;
            // int32_t len2enc;
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

            bsize = std::min(tab[level], BLOCK_SIZE); /* Block size */
            len_buf = bsize;
            len_bwt = bsize;
            len_arr = bsize << 2; /* int32_t type */
            coder_buff = static_cast<uint8_t*>(safe_alloc(len_buf + len_bwt + len_arr));
            check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "Error: Insufficient memory: need %" PRIu64 " MB\n",
                       (static_cast<uint64_t>(sizeof(uint8_t)) *
                        ((uint64_t)len_buf + len_bwt + len_arr)) >> 20);

            buf_bwt = coder_buff + len_buf;
            buf_arr = (uint32_t *)(buf_bwt + len_bwt);
        }

        int32_t len2add = bsize - curr_incache;
        if ((int32_t)in_len < len2add)
        { /* Cache when not enough for one block */
            memcpy(coder_buff + curr_incache, in, in_len);
            curr_incache += in_len;
        }
        else
        { /* Compress when enough for one block and cache uncompressed part */
            memcpy(coder_buff + curr_incache, in, len2add);

            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, bsize);
            check_exit(idx >= 1,  coder_ns::CODER_ERR_INNER, "coder transform failed: idx %d", idx);
            // if (!(idx >= 1)) {
            //     coder_logger(coder_ns::ERROR, "coder transform failed: idx %d", idx);
            //     return coder_ns::CODER_ERR_INNER;
            // }

            put32(bsize); /* Encode current length */
            put32(idx);   /* Encode BWT index */

            for (i = 0; i < bsize; i++) /* Encode transformed data */
                put(buf_bwt[i]);

            curr_incache = in_len - len2add;
            if (curr_incache > 0)
                memcpy(coder_buff, in + len2add, curr_incache);
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
            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, curr_incache);
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

    /* External decompression interface, returns actual decompressed length, exits when encountering split_ch during decompression */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold __attribute__ ((unused)) = false)
    {
        uint8_t ch;
        int32_t len = 0, n;

        /*
         * 错误约定（全 coder 统一）：返回值 < 0 是错误（码流损坏/缓冲不足），
         * == 0 是干净 EOF，> 0 是实际解码长度。decode_one_block 的负值原样上抛，
         * 不能与 EOF 混为一谈——否则损坏流会被当成正常结尾，产出截断却"成功"的文件。
         */
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;

        if (io->m != coder_io::MDEC)
        {
            decode_init();
            if (io->err != coder_io::IO_OK)
                return coder_ns::CODER_ERR_STREAM_END;
            /* decode_one_block 返回负值代表码流损坏或内存不足，此时 coder_buff 可能
               仍是空指针，必须与正常的 EOF（返回 0）区分并原样上抛。 */
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
                        /* 没找到分隔符就写满：输出缓冲不足。原实现静默返回 len，
                           上层会把截断行当完整行写进文件。 */
                        curr_out_offset = ++n;
                        io->err = coder_io::IO_BUF_FULL;
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
                        /* 定长模式：写满 out_len 即正常结束，剩余留待下次调用。 */
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        }

stream_end:
        /* 流提前耗尽（_coder 读端 latch 了 IO_READ_EMPTY）是码流损坏，不是正常结束。 */
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool)
    {
        return 0;
    }

    /* Decompression initialization */
    void decode_init()
    {
        bsize = 0;
        coder->decode_init();
    }

    /* External interface to get the length of currently consumed data */
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
        const int32_t f = (run > 2); /* run represents the number of consecutive repetitions of previous characters */

        int32_t ctx = 1;
        for (int32_t i = 128; i > 0; i >>= 1) /* Get all bits of current byte c, starting from high bit */
        {
            const int32_t p0 = counter0[ctx].get();
            const int32_t p1 = counter1[c1][ctx].get();
            const int32_t p2 = counter1[c2][ctx].get();
            const int32_t p = (((p0 + p1) * 7) + p2 + p2) >> 4;

            /* Get probability of auxiliary symbols through linear interpolation */
            const int32_t j = p >> 12; /* j takes the highest 4 bits of mixed probability, uses SSE combined with ctx (i.e., current character c to be processed) to correct p_mixed, making another prediction on the probability itself using ctx */
            const int32_t x1 = counter2[f][ctx][j].get();
            const int32_t x2 = counter2[f][ctx][j + 1].get();
            const int32_t ssep = x1 + (((x2 - x1) * (p & 4095)) >> 12);

            if (c & i) /* Current bit is 1 */
            {
                coder->encode_bit<18>(1, p + ssep + ssep + ssep); /* p + ssep + ssep + ssep represents the probability of high probability symbols, i.e., the probability that next bit is 1 */

                counter0[ctx].update_high();
                counter1[c1][ctx].update_high();
                counter2[f][ctx][j].update_high();
                counter2[f][ctx][j + 1].update_high();

                ctx += ctx + 1; /* ctx left shift by one bit + 1, here ctx is the bit context of current byte */
            }
            else /* Current bit is 0 */
            {
                coder->encode_bit<18>(0, p + ssep + ssep + ssep);

                counter0[ctx].update_low();
                counter1[c1][ctx].update_low();
                counter2[f][ctx][j].update_low();
                counter2[f][ctx][j + 1].update_low();

                ctx += ctx; /* ctx left shift by one bit, here ctx is the bit context of current byte */
            }
        }

        c2 = c1;        /* c2 updates to ctx of previous byte, here ctx is the actual value of processed character */
        c1 = ctx - 256; /* c1 updates to ctx of current processed byte. Other: after processing current byte, ctx is at least 256, because even in extreme case where current byte is all 0, ctx == 1 << 8 == 256 */

        if (c1 == c2)   /* If current byte repeats with previous byte, then increment run */
            ++run;
        else
            run = 0;
    }

    inline int32_t get()
    {
        const int32_t f = (run > 2);

        int32_t ctx = 1;
        while (ctx < 256)
        {
            const int32_t p0 = counter0[ctx].get();
            const int32_t p1 = counter1[c1][ctx].get();
            const int32_t p2 = counter1[c2][ctx].get();
            const int32_t p = (((p0 + p1) * 7) + p2 + p2) >> 4;

            /* Get probability of auxiliary symbols through linear interpolation */
            const int32_t j = p >> 12;
            const int32_t x1 = counter2[f][ctx][j].get();
            const int32_t x2 = counter2[f][ctx][j + 1].get();
            const int32_t ssep = x1 + (((x2 - x1) * (p & 0xFFF)) >> 12);

            if (coder->decode_bit<18>(p + ssep + ssep + ssep))
            {
                counter0[ctx].update_high();
                counter1[c1][ctx].update_high();
                counter2[f][ctx][j].update_high();
                counter2[f][ctx][j + 1].update_high();

                ctx += ctx + 1; /* ctx left shift by one bit, here ctx is the bit context of current byte */
            }
            else
            {
                counter0[ctx].update_low();
                counter1[c1][ctx].update_low();
                counter2[f][ctx][j].update_low();
                counter2[f][ctx][j + 1].update_low();

                ctx += ctx; /* ctx left shift by one bit, here ctx is the bit context of current byte */
            }
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
            return coder_ns::CODER_ERR_STREAM_END;  // 读 bsize 时流已耗尽：码流被截断
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
        /* 反变换期间若码流耗尽（_coder latch 了 IO_READ_EMPTY），结果已是垃圾，按损坏上报。 */
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

    uint8_t *coder_buff;  /* Buffer for caching original data, transform and process only when enough for one block */
    int32_t curr_incache; /* Length of currently cached data */
    uint8_t *buf_bwt;
    uint32_t *buf_arr;
    int32_t bsize;
    int32_t bidx;
    int32_t curr_out_offset; /* Offset position during decompression */
    int32_t cnt[257];
    bool flushed;

    _coder *coder;
};
#endif
