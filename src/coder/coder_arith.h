/*
 * coder_arith.h - order-0 adaptive arithmetic byte-stream coder.
 *
 * Two backends: flags 0x01 legacy cumulative model (decode-only) and flags
 * 0x02 current htscodecs-style carry-less range coder (encode+decode).
 * Input smaller than one block is cached until a whole block can be committed.
 * Block format: [u32 LE length][u8 flags][payload].
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

#ifndef _CODER_ARITH_H_
#define _CODER_ARITH_H_

#include <stdint.h>
#include <algorithm>
#include <vector>

#include "coder_io.h"
#include "coder.h"

class coder_arith : public coder
{
public:
    /* Prior: 256 little-endian uint16 symbol weights. Encode and decode must use the same prior. */
    static constexpr uint32_t kPriorWeights = 256;
    static constexpr uint32_t kPriorBytes = kPriorWeights * 2;

    coder_arith(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_arith");

        coder_buff = nullptr;
        buff_capacity = 0;
        curr_incache = 0;
        curr_out_offset = 0;
        curr_blk_len = 0;
        bsize = 0;
        flushed = false;
    }

    /* Supply the file-level order-0 byte-frequency prior before first use. */
    void set_prior(const uint8_t *prior, uint32_t len)
    {
        prior_.assign(prior, prior + len);
    }

    virtual ~coder_arith()
    {
        if (!flushed && io->m == coder_io::MENC)
            encode_flush();
        if (coder_buff)
            safe_free((void**)&coder_buff);
    }

    /* External compression interface. */
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
            const int32_t BLOCK_SIZE = (268435456);

            io->m = coder_io::MENC;
            int32_t level = (io->meta["level"].isInt()) ? (io->meta["level"].asInt()) : 4;
            io->meta["level"] = (Json::Value::Int)level;
            /* Clamp level to 1..9 (0 would break buffering). */
            if (level < 1) level = 1;
            if (level > 9) level = 9;

            bsize = std::min(tab[level], BLOCK_SIZE); /* Block size (split threshold) */
            /* Grow the buffer on demand instead of pre-allocating bsize. */
            const int32_t initCap = std::min(bsize, ARITH_INIT_BUFFER);
            coder_buff = static_cast<uint8_t*>(safe_alloc(initCap));
            check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_arith: insufficient memory for %d-byte block buffer", initCap);
            buff_capacity = initCap;
        }

        uint32_t off = 0;
        while (off < in_len)
        {
            if (curr_incache == buff_capacity)
            {
                if (buff_capacity < bsize)
                {
                    /* Grow instead of splitting so the field stays one coded block. */
                    const int32_t ncap = (int32_t)std::min((int64_t)buff_capacity * 2, (int64_t)bsize);
                    uint8_t *nb = static_cast<uint8_t *>(safe_alloc(ncap));
                    if (nb != nullptr) {
                        memcpy(nb, coder_buff, (size_t)curr_incache);
                        safe_free((void **)&coder_buff);
                        coder_buff = nb;
                        buff_capacity = ncap;
                    } else {
                        /* Out of memory: fall back to splitting here. */
                        encode_one_block(curr_incache);
                        curr_incache = 0;
                    }
                }
                else
                {
                    encode_one_block(bsize);
                    curr_incache = 0;
                }
            }
            /* Commit whole blocks as they fill; keep only the tail cached. */
            const int32_t len2add = buff_capacity - curr_incache;
            const uint32_t take = std::min((uint32_t)(in_len - off), (uint32_t)len2add);
            memcpy(coder_buff + curr_incache, in + off, take);
            curr_incache += (int32_t)take;
            off += take;
        }
    }

    /* Call this after encoding to compress cached data */
    void encode_flush() override
    {
        if (io->m != coder_io::MENC || flushed)
            return;
        /* First compress all cached data */
        if (curr_incache > 0)
            encode_one_block(curr_incache);
        /* Then set EOF marker */
        put32(0);
        flushed = true;
    }

    /* Decompress; returns actual length, stops at split_ch. */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX,
                        bool need2hold [[maybe_unused]] = false) override
    {
        uint8_t ch;
        int32_t len = 0, n;

        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;

        if (io->m != coder_io::MDEC)
        {
            io->m = coder_io::MDEC;
            int32_t ret = decode_one_block();
            if (ret < 0)
                return ret;
            if (ret == 0)
                return 0;
        }

        if (split_ch != UINT8_MAX)
        {
            for (;;)
            {
                if (curr_out_offset == curr_blk_len)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < curr_blk_len; n++)
                {
                    ch = coder_buff[n];
                    *(out + len++) = ch;
                    if (ch == split_ch)
                    {
                        curr_out_offset = ++n;
                        return len;
                    }
                    if (len >= (int32_t)out_len)
                    {
                        /* Output buffer too small. */
                        curr_out_offset = ++n;
                        io->set_err(coder_io::IO_BUF_FULL);
                        return coder_ns::CODER_ERR_BUF_SMALL;
                    }
                }
                curr_out_offset = n;
            }
        }
        else
        {
            for (;;)
            {
                if (curr_out_offset == curr_blk_len)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < curr_blk_len; n++)
                {
                    *(out + len++) = coder_buff[n];
                    if (len >= (int32_t)out_len)
                    {
                        /* Fixed-length mode: filling out_len is a normal end; the rest is left for the next call. */
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        }

stream_end:
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool) override
    {
        return 0;
    }

    /* External interface to get the length of currently consumed data */
    int32_t decode_inlen()
    {
        return this->io->data_len;
    }

private:
    /* Order-0 adaptive model with the project's Ppmd8 carryless range coder.
       Total kept <= 2^14; rescale keeps the table strictly increasing. */

    static constexpr uint32_t kTop = 1u << 24;
    static constexpr uint32_t kBot = 1u << 15;

    struct ArithEnc {
        std::vector<uint8_t> *out;
        uint32_t low;
        uint32_t range;

        explicit ArithEnc(std::vector<uint8_t> *o) : out(o), low(0), range(0xFFFFFFFFu) {}

        void normalize()
        {
            while ((low ^ (low + range)) < kTop ||
                   (range < kBot && ((range = (0 - low) & (kBot - 1)), 1)))
            {
                out->push_back((uint8_t)(low >> 24));
                range <<= 8;
                low <<= 8;
            }
        }

        void encode(uint32_t start, uint32_t size, uint32_t total)
        {
            low += start * (range /= total);
            range *= size;
            normalize();
        }

        void flush()
        {
            for (int i = 0; i < 4; ++i)
            {
                out->push_back((uint8_t)(low >> 24));
                low <<= 8;
            }
        }
    };

    struct ArithDec {
        coder_io *io;
        uint32_t low;
        uint32_t range;
        uint32_t code;

        explicit ArithDec(coder_io *i) : io(i), low(0), range(0xFFFFFFFFu), code(0) {}

        void init()
        {
            low = 0;
            range = 0xFFFFFFFFu;
            code = 0;
            for (int i = 0; i < 4; ++i)
                code = (code << 8) | (uint32_t)io->getc();
        }

        uint32_t getThreshold(uint32_t total)
        {
            return code / (range /= total);
        }

        void decode(uint32_t start, uint32_t size, uint32_t total)
        {
            (void)total;
            start *= range;
            low += start;
            code -= start;
            range *= size;
            while ((low ^ (low + range)) < kTop ||
                   (range < kBot && ((range = (0 - low) & (kBot - 1)), 1)))
            {
                code = (code << 8) | (uint32_t)io->getc();
                range <<= 8;
                low <<= 8;
            }
        }
    };

    struct Order0Model {
        uint32_t freq[257]; /* cumulative: freq[0]=0, freq[256]=total */

        Order0Model()
        {
            for (int i = 0; i <= 256; ++i)
                freq[i] = (uint32_t)i; /* uniform */
        }

        /* Rebuild from a 256 x uint16 LE weight blob; zero blob keeps uniform. */
        void init_from_weights(const uint8_t *w)
        {
            uint32_t wt[256];
            uint64_t total = 0;
            for (int i = 0; i < 256; ++i)
            {
                wt[i] = (uint32_t)w[2 * i] | ((uint32_t)w[2 * i + 1] << 8);
                total += wt[i];
            }
            if (total == 0)
                return;
            const uint64_t target = 1u << 14;
            uint64_t acc = 0, last = 0;
            freq[0] = 0;
            for (int i = 0; i < 256; ++i)
            {
                acc += wt[i];
                uint64_t v = (acc * target) / total;
                if (v <= last)
                    v = last + 1;
                freq[i + 1] = (uint32_t)v;
                last = v;
            }
        }

        void rescale()
        {
            /* Halve values, keeping the table strictly increasing. */
            freq[0] = 0;
            for (int i = 1; i <= 256; ++i)
            {
                const uint32_t h = (freq[i] + 1) >> 1;
                freq[i] = (h > freq[i - 1] + 1) ? h : freq[i - 1] + 1;
            }
            if (freq[256] > (1u << 16)) /* safety: still too large, halve again */
                rescale();
        }

        void adapt(uint8_t sym)
        {
            for (int i = sym + 1; i <= 256; ++i)
                ++freq[i];
            if (freq[256] > (1u << 14))
                rescale();
        }

        void encode(ArithEnc &rc, uint8_t sym)
        {
            const uint32_t total = freq[256];
            const uint32_t lo = freq[sym];
            const uint32_t hi = freq[sym + 1];
            rc.encode(lo, hi - lo, total);
            adapt(sym);
        }

        uint8_t decode(ArithDec &rc)
        {
            const uint32_t total = freq[256];
            const uint32_t target = rc.getThreshold(total);
            int lo = 0, hi = 255;
            while (lo < hi)
            {
                const int mid = (lo + hi + 1) >> 1;
                if (freq[mid] <= target)
                    lo = mid;
                else
                    hi = mid - 1;
            }
            const uint8_t sym = (uint8_t)lo;
            const uint32_t loC = freq[sym];
            const uint32_t hiC = freq[sym + 1];
            rc.decode(loC, hiC - loC, total);
            adapt(sym);
            return sym;
        }
    };

    /* V2 backend: carry-less range coder (2^24 floor) + 12-bit order-0 model. */
    static constexpr uint32_t kTopV2 = 1u << 24;
    static constexpr uint32_t kThresV2 = 255u * kTopV2;

    /* Initial encode buffer size; grows on demand up to bsize. */
    static constexpr int32_t ARITH_INIT_BUFFER = 1 << 20; /* 1 MB */

    struct RcEncV2 {
        coder_io *io;
        uint32_t low;
        uint32_t range;
        uint32_t ffNum; /* pending 0xFF bytes */
        uint32_t cache; /* pending top byte, not yet emitted */
        uint32_t carry; /* carry into cache */

        explicit RcEncV2(coder_io *i)
            : io(i), low(0), range(0xFFFFFFFFu), ffNum(0), cache(0), carry(0) {}

        void shift_low()
        {
            if (low < kThresV2 || carry) {
                io->putc((uint8_t)(cache + carry));
                while (ffNum) {
                    io->putc((uint8_t)(carry - 1)); /* 0xFF + carry */
                    --ffNum;
                }
                cache = low >> 24;
                carry = 0;
            } else {
                /* 0xFFxxxxxx: a carry can still be absorbed into this 0xFF */
                ++ffNum;
            }
            low <<= 8;
        }

        void encode(uint32_t cum, uint32_t freq, uint32_t total)
        {
            const uint32_t prev = low;
            low += cum * (range /= total);
            range *= freq;
            /* 32-bit overflow of low is the pending carry, folded in by shift_low. */
            carry += (low < prev) ? 1u : 0u;
            while (range < kTopV2) {
                range <<= 8;
                shift_low();
            }
        }

        void flush()
        {
            for (int i = 0; i < 5; ++i)
                shift_low();
        }
    };

    struct RcDecV2 {
        coder_io *io;
        uint32_t range;
        uint32_t code;

        explicit RcDecV2(coder_io *i) : io(i), range(0xFFFFFFFFu), code(0) {}

        void init()
        {
            range = 0xFFFFFFFFu;
            code = 0;
            for (int i = 0; i < 5; ++i)
                code = (code << 8) | (uint32_t)io->getc();
        }

        uint32_t get_freq(uint32_t total)
        {
            return (total != 0 && range >= total) ? (code / (range /= total)) : 0;
        }

        void decode(uint32_t cum, uint32_t freq, [[maybe_unused]] uint32_t total)
        {
            code -= cum * range;
            range *= freq;
            while (range < kTopV2) {
                code = (code << 8) | (uint32_t)io->getc();
                range <<= 8;
            }
        }
    };

    struct Order0ModelV2 {
        static constexpr uint32_t kMaxTotal = 1u << 16;
        static constexpr uint32_t kStep = 16;

        uint32_t freq[257]; /* cumulative: freq[0]=0, freq[256]=total */

        Order0ModelV2()
        {
            for (int i = 0; i <= 256; ++i)
                freq[i] = (uint32_t)i; /* uniform start, total = 256 */
        }

        /* Uniform weight 1 over [0, max_sym), zero above: impossible symbols get no probability mass. */
        void init_default(uint32_t max_sym)
        {
            if (max_sym == 0 || max_sym > 256)
                max_sym = 256;
            for (uint32_t i = 0; i <= 256; ++i)
                freq[i] = (i <= max_sym) ? i : max_sym;
        }

        /* Rebuild from a 256 x uint16 LE weight blob; zero blob keeps uniform. */
        void init_from_weights(const uint8_t *w)
        {
            uint32_t wt[256];
            uint64_t total = 0;
            for (int i = 0; i < 256; ++i) {
                wt[i] = (uint32_t)w[2 * i] | ((uint32_t)w[2 * i + 1] << 8);
                total += wt[i];
            }
            if (total == 0)
                return;
            const uint64_t target = kMaxTotal;
            uint64_t acc = 0, last = 0;
            freq[0] = 0;
            for (int i = 0; i < 256; ++i) {
                acc += wt[i];
                uint64_t v = (acc * target) / total;
                /* Keep every weight >= 1 so no symbol gets a zero-width interval. */
                if (v <= last)
                    v = last + 1;
                freq[i + 1] = (uint32_t)v;
                last = v;
            }
        }

        void rescale()
        {
            /* Ceil-halve each weight; never-seen symbols stay at 0 (no probability mass). */
            uint32_t acc = 0;
            for (int i = 0; i < 256; ++i) {
                uint32_t w = freq[i + 1] - freq[i];
                w = (w + 1) >> 1;
                freq[i] = acc;
                acc += w;
            }
            freq[256] = acc;
            if (acc > kMaxTotal)
                rescale();
        }

        void adapt(uint8_t sym)
        {
            for (int i = (int)sym + 1; i <= 256; ++i)
                freq[i] += kStep;
            if (freq[256] >= kMaxTotal)
                rescale();
        }

        void encode(RcEncV2 &rc, uint8_t sym)
        {
            const uint32_t total = freq[256];
            rc.encode(freq[sym], freq[sym + 1] - freq[sym], total);
            adapt(sym);
        }

        uint8_t decode(RcDecV2 &rc)
        {
            const uint32_t total = freq[256];
            const uint32_t target = rc.get_freq(total);
            /* First symbol whose upper bound exceeds target, skipping empty intervals. */
            int lo = 0, hi = 255;
            while (lo < hi) {
                const int mid = (lo + hi) >> 1;
                if (freq[mid + 1] > target)
                    hi = mid;
                else
                    lo = mid + 1;
            }
            const uint8_t sym = (uint8_t)lo;
            rc.decode(freq[sym], freq[sym + 1] - freq[sym], total);
            adapt(sym);
            return sym;
        }
    };

    /* Write one block plus header; falls back to raw if not shorter. */
    void encode_one_block(int32_t n)
    {
        put32((uint32_t)n);

        /* Temporary buffer sized worst case plus flush; bounds every write. */
        std::vector<uint8_t> coded((size_t)n + 64);
        coder_io tmp(coded.data(), (int32_t)coded.size());
        const bool use_prior = (prior_.size() >= kPriorBytes);
        {
            RcEncV2 enc(&tmp);
            Order0ModelV2 model;
            if (use_prior)
                model.init_from_weights(prior_.data());
            else {
                uint32_t maxv = 0;
                for (int32_t i = 0; i < n; i++)
                    if ((uint32_t)coder_buff[i] > maxv)
                        maxv = coder_buff[i];
                model.init_default(maxv + 1);
            }
            for (int32_t i = 0; i < n; i++)
                model.encode(enc, coder_buff[i]);
            enc.flush();
        }
        const int32_t clen = tmp.data_len;

        if (tmp.err != coder_io::IO_OK || clen >= n)
        {
            io->putc(0x00); /* flags: raw */
            for (int32_t i = 0; i < n; i++)
                io->putc(coder_buff[i]);
            return;
        }

        if (use_prior)
        {
            /* 0x03: model rebuilt from the file-level prior, no alphabet byte */
            io->putc(0x03);
        }
        else
        {
            /* 0x02: alphabet bound m = max byte + 1 in one byte; 0 means 256 */
            uint32_t maxv = 0;
            for (int32_t i = 0; i < n; i++)
                if ((uint32_t)coder_buff[i] > maxv)
                    maxv = coder_buff[i];
            const uint32_t m = maxv + 1;
            io->putc(0x02);
            io->putc((uint8_t)(m & 0xff));
        }
        for (int32_t i = 0; i < clen; i++)
            io->putc(coded[i]);
    }

    /* Decode one block into coder_buff; returns length, 0 at end, <0 on error. */
    int32_t decode_one_block()
    {
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        int32_t n = (int32_t)get32();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END; /* stream truncated while reading length */
        if (n == 0)
            return 0;
        if (n < 0 || n > (1 << 30)) /* sanity bound for a corrupt header */
        {
            coder_logger(coder_ns::ERROR, "coder_arith: block size %d is invalid", n);
            return coder_ns::CODER_ERR_INNER;
        }

        ensure_buffer((uint32_t)n);

        const int32_t flags = io->getc();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        if (flags == 0x00)
        {
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = io->getc();
        }
        else if (flags == 0x01)
        { /* legacy backend, kept so files written before the switch still decode */
            ArithDec dec(io);
            dec.init();
            Order0Model model;
            if (prior_.size() >= kPriorBytes)
                model.init_from_weights(prior_.data());
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = model.decode(dec);
        }
        else if (flags == 0x02)
        {
            const uint32_t mb = io->getc();
            const uint32_t max_sym = mb ? mb : 256; /* 0 encodes 256 */
            RcDecV2 dec(io);
            dec.init();
            Order0ModelV2 model;
            model.init_default(max_sym);
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = model.decode(dec);
        }
        else if (flags == 0x03)
        { /* model rebuilt from the same file-level prior the encoder used */
            RcDecV2 dec(io);
            dec.init();
            Order0ModelV2 model;
            if (prior_.size() >= kPriorBytes)
                model.init_from_weights(prior_.data());
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = model.decode(dec);
        }
        else
        {
            coder_logger(coder_ns::ERROR, "coder_arith: invalid block flags %d", flags);
            return coder_ns::CODER_ERR_INNER;
        }

        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        curr_blk_len = n;
        curr_out_offset = 0;
        return n;
    }

    /* Grow the internal buffer (capacity only ever increases). */
    void ensure_buffer(uint32_t need)
    {
        if (need <= (uint32_t)buff_capacity)
            return;
        uint8_t *nb = static_cast<uint8_t*>(safe_alloc(need));
        if (!nb)
        {
            coder_logger(coder_ns::ERROR, "coder_arith: insufficient memory for %u-byte block buffer", need);
            io->set_err(coder_io::IO_BUF_FULL);
            return;
        }
        if (coder_buff)
        {
            memcpy(nb, coder_buff, buff_capacity);
            safe_free((void**)&coder_buff);
        }
        coder_buff = nb;
        buff_capacity = (int32_t)need;
    }

    inline void put32(uint32_t x)
    {
        io->putc((uint8_t)(x & 0xff));
        io->putc((uint8_t)((x >> 8) & 0xff));
        io->putc((uint8_t)((x >> 16) & 0xff));
        io->putc((uint8_t)((x >> 24) & 0xff));
    }

    inline uint32_t get32()
    {
        uint32_t x = 0;
        x |= (uint32_t)io->getc();
        x |= (uint32_t)io->getc() << 8;
        x |= (uint32_t)io->getc() << 16;
        x |= (uint32_t)io->getc() << 24;
        return x;
    }

private:
    uint8_t *coder_buff;  /* Block cache for encoding / decoded block for decoding */
    int32_t buff_capacity; /* Capacity of coder_buff */
    int32_t curr_incache; /* Length of currently cached encode data */
    int32_t curr_out_offset; /* Offset position during decompression */
    int32_t curr_blk_len;    /* Source length of the currently decoded block */
    int32_t bsize;           /* Encode block size selected by level */
    bool flushed;
    std::vector<uint8_t> prior_; /* File-level frequency prior (empty = none) */
};
#endif
