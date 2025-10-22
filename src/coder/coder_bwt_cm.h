#ifndef _CODER_BWT_CM_H_
#define _CODER_BWT_CM_H_

/* 基于bwt变换后的context model编码器 */
#include <stdint.h>
#include <algorithm>
#include <cinttypes> 

#include "coder_io.h"
#include "coder.h"
#include "bcm_libsais.h"


class coder_bwt_cm: public coder
{
private:
    template <int32_t RATE> /* 用count表示概率，RATE速率值越小，表示越灵敏，影响越大，即每次来一个符号时概率很大 */
    class _counter
    {
    private:
        uint16_t c;

    public:
        _counter()
        {
            c = 1 << 15; /* 0.5 ， 32768，再左移一位（乘以2）就是最大值 */
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
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
            }
            flushed = true;
        }

        template <int32_t P_LOG>
        inline void encode_bit(int32_t bit, uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG); /* 这里是把高概率符号放前面了，低概率放后面，P表示高概率符号的概率 */

            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* 归一化 */
            while ((low ^ high) < (1 << 24))
            {
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

            /* 归一化 */
            while ((low ^ high) < (1 << 24))
            {
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
                code = (code << 8) + *(io->data + io->data_len++);
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

    /* 对外压缩接口 */
    void encode_line(const uint8_t *in, const int32_t in_len)
    {
        int32_t i;
        if (io->m != coder_io::MENC)
        { /* 未初始化时先初始化 */
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
            // if (!(level <= 9 && level >= 0)) {
            //     coder_logger(coder_ns::ERROR, "coder level should in [0, 9], current is %d", level);
            //     return coder_ns::CODER_ERR_BAD_ARGS;
            // }
            const int32_t BLOCK_SIZE = (268435456);

            bsize = std::min(tab[level], BLOCK_SIZE); /* 块大小 */
            len_buf = bsize;
            len_bwt = bsize;
            len_arr = bsize << 2; /* int32_t类型 */
            coder_buff = static_cast<uint8_t*>(safe_alloc(len_buf + len_bwt + len_arr));

            buf_bwt = coder_buff + len_buf;
            buf_arr = (uint32_t *)(buf_bwt + len_bwt);
        }

        int32_t len2add = bsize - curr_incache;
        if (in_len < len2add)
        { /* 不够一个block时缓存 */
            memcpy(coder_buff + curr_incache, in, in_len);
            curr_incache += in_len;
        }
        else
        { /* 够一个block时压缩再缓存未压缩部分 */
            memcpy(coder_buff + curr_incache, in, len2add);

            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, bsize);
            check_exit(idx >= 1,  coder_ns::CODER_ERR_INNER, "coder transform failed: idx %d", idx);
            // if (!(idx >= 1)) {
            //     coder_logger(coder_ns::ERROR, "coder transform failed: idx %d", idx);
            //     return coder_ns::CODER_ERR_INNER;
            // }

            put32(bsize); /* 编码当前长度 */
            put32(idx);   /* 编码bwt的index */

            for (i = 0; i < bsize; i++) /* 编码变换后的数据 */
                put(buf_bwt[i]);

            curr_incache = in_len - len2add;
            if (curr_incache > 0)
                memcpy(coder_buff, in + len2add, curr_incache);
        }
    }

    /* 编码完成后需要调用它将缓存中的数据压缩 */
    void encode_flush()
    {
        int32_t i;
        if (io->m != coder_io::MENC || flushed) {
            return;
        }
        /* 先将缓存中的数据压缩完 */
        if (curr_incache > 0)
        {
            const int32_t idx = libsais_bwt(coder_buff, buf_bwt, (int32_t *)buf_arr, curr_incache);
            check_exit(idx >= 1,  coder_ns::CODER_ERR_INNER, "coder transform failed: idx %d", idx);
            // if (!(idx >= 1)) {
            //     coder_logger(coder_ns::ERROR, "coder transform failed: idx %d", idx);
            //     return coder_ns::CODER_ERR_INNER;
            // }

            put32(curr_incache); /* 编码当前长度 */
            put32(idx);          /* 编码bwt的index */

            for (i = 0; i < curr_incache; i++) /* 编码变换后的数据 */
                put(buf_bwt[i]);
        }
        /* 再置结束标识 EOF */
        put32(0);
        coder->encode_flush();
        flushed = true;
    }

    /* 对外解压接口，返回实际解压出来的长度，当解压遇到split_ch时退出 */
    int32_t decode_line(uint8_t *out, int32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold = false)
    {
        uint8_t ch;
        int32_t len = 0, n;

        if (io->m != coder_io::MDEC)
        {
            decode_init();
            if (decode_one_block() == 0)
                return 0;
            io->m = coder_io::MDEC;
        }

        if (split_ch != UINT8_MAX)
        {
            for (;;)
            {
                if (curr_out_offset == bsize)
                { /* 当前block已经解压完成了，需要重新解压 */
                    if (decode_one_block() == 0)
                        return len;
                }
                for (n = curr_out_offset; n < bsize; n++) {
                    ch = coder_buff[n];
                    *(out + len++) = ch;
                    if (ch == split_ch) {
                        curr_out_offset = ++n;
                        return len;
                    }
                    if (len >= out_len) {
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        } else {
            for (;;)
            {
                if (curr_out_offset == bsize)
                { /* 当前block已经解压完成了，需要重新解压 */
                    if (decode_one_block() == 0)
                        return len;
                }
                for (n = curr_out_offset; n < bsize; n++) {
                    *(out + len++) = coder_buff[n];
                    if (len >= out_len) {
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        }
        return 0;
    }

    /* fake */
    int32_t decode_line(uint8_t *out, int32_t out_len, uint8_t *rely = nullptr, uint8_t split_ch = UINT8_MAX, bool need2hold = false)
    {
        return 0;
    }

    /* 解压初始化 */
    void decode_init()
    {
        bsize = 0;
        coder->decode_init();
    }

    /* 对外接口，获取当前已经消耗数据的长度 */
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
        const int32_t f = (run > 2); /* run表示前面字符连续重复的次数 */

        int32_t ctx = 1;
        for (int32_t i = 128; i > 0; i >>= 1) /* 取当前byte c的所有bit，从高位开始取 */
        {
            const int32_t p0 = counter0[ctx].get();
            const int32_t p1 = counter1[c1][ctx].get();
            const int32_t p2 = counter1[c2][ctx].get();
            const int32_t p = (((p0 + p1) * 7) + p2 + p2) >> 4;

            /* 通过线性插值得到辅助符号的概率 */
            const int32_t j = p >> 12; /* j取混合概率的最高4位 , 借助sse结合ctx(即当前待处理字符c)对p_mixed做一个修正，将概率本身做ctx再做一次预测 */
            const int32_t x1 = counter2[f][ctx][j].get();
            const int32_t x2 = counter2[f][ctx][j + 1].get();
            const int32_t ssep = x1 + (((x2 - x1) * (p & 4095)) >> 12);

            if (c & i) /*当前bit为1 */
            {
                coder->encode_bit<18>(1, p + ssep + ssep + ssep); /* p + ssep + ssep + ssep表示高概率符号的概率，即下个bit是1的概率 */

                counter0[ctx].update_high();
                counter1[c1][ctx].update_high();
                counter2[f][ctx][j].update_high();
                counter2[f][ctx][j + 1].update_high();

                ctx += ctx + 1; /* ctx左移一位 + 1，这里ctx为当前byte的bit的ctx */
            }
            else /*当前bit为0 */
            {
                coder->encode_bit<18>(0, p + ssep + ssep + ssep);

                counter0[ctx].update_low();
                counter1[c1][ctx].update_low();
                counter2[f][ctx][j].update_low();
                counter2[f][ctx][j + 1].update_low();

                ctx += ctx; /* ctx左移一位，这里ctx为当前byte的bit的ctx */
            }
        }

        c2 = c1;        /* c2更新为上个字节的ctx，这里ctx即为处理字符的实际值 */
        c1 = ctx - 256; /* c1更新为当前处理字节的ctx。其他：处理完当前字节后，ctx至少为256，因为即使当前字节全为0这种极端情况ctx==1 << 8==256 */

        if (c1 == c2)   /* 如果上当前字节与上个字节重复了，那么设置run自加1 */
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

            /* 通过线性插值得到辅助符号的概率 */
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

                ctx += ctx + 1; /* ctx左移一位，这里ctx为当前byte的bit的ctx */
            }
            else
            {
                counter0[ctx].update_low();
                counter1[c1][ctx].update_low();
                counter2[f][ctx][j].update_low();
                counter2[f][ctx][j + 1].update_low();

                ctx += ctx; /* ctx左移一位，这里ctx为当前byte的bit的ctx */
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
        bsize = get32();
        if (bsize == 0)
            return 0;
        // check_exit(bsize > 0,  coder_ns::CODER_ERR_INNER, "block size %d is invalid", bsize);
        if (!(bsize > 0)) {
            coder_logger(coder_ns::ERROR, "block size %d is invalid", bsize);
            return coder_ns::CODER_ERR_INNER;
        }


        bidx = get32();
        // check_exit(bidx >= 1 && bidx <= bsize,  coder_ns::CODER_ERR_INNER, "corrupt input");
        if (!(bidx >= 1 && bidx <= bsize)) {
            coder_logger(coder_ns::ERROR, "corrupt input");
            return coder_ns::CODER_ERR_INNER;
        }

        if (!buf_arr)
        {                           /* 没有分配空间，先分配空间 */
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
            /* BW逆变换 */
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
            /* BW逆变换 */
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
        return bsize;
    }
    
private:
    int32_t run;
    int32_t c1;
    int32_t c2;
    _counter<2> counter0[256]; /* counter0影响最大 */
    _counter<4> counter1[256][256];
    _counter<6> counter2[2][256][17];

    uint8_t *coder_buff;  /* 缓存原始数据的buffer，累够一个块时才变换然后处理 */
    int32_t curr_incache; /* 当前缓存数据的长度 */
    uint8_t *buf_bwt;
    uint32_t *buf_arr;
    int32_t bsize;
    int32_t bidx;
    int32_t curr_out_offset; /* 解压时offset位置 */
    int32_t cnt[257];
    bool flushed;

    _coder *coder;
};
#endif
