/*
 * coder_qcm.h - order-2 context model for quality-value streams
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _CODER_QCM_H_
#define _CODER_QCM_H_

#include <stdint.h>
#include <string.h>

#include <vector>

#include "coder_io.h"
#include "coder.h"

/*
 * Tuning constants, macros so a measurement build can sweep them.
 *   QCM_THETA  weight of the lower-order estimate in the interpolation.
 *   QCM_CAP    per-node count cap; a node reaching it is halved.
 */
#ifndef QCM_THETA
#define QCM_THETA 64
#endif
#ifndef QCM_CAP
#define QCM_CAP 16384
#endif

/*
 * coder_qcm - an order-2 context model for the quality-value column.
 *
 * Same external contract as the other byte-stream coders: encode_line
 * accumulates a column and encode_flush commits one independent block. It needs
 * no per-record input (no record boundaries or strand).
 *
 * A bit tree over the alphabet index, with three nested orders (global,
 * previous value, previous two values) combined by interpolated PPM smoothing.
 * Order 2 is enabled only for small alphabets (<=128); the block carries the
 * sorted list of distinct values and codes symbols as indices into it.
 *
 * Block format: [u8 version][u16 LE alphabet size K][K bytes alphabet]
 * [u32 LE value count][range-coded symbols]. Error convention matches the other
 * coders.
 */
class coder_qcm : public coder
{
private:
    /*
     * Range coder: 32-bit range, 64-bit low with a delayed carry byte. The split
     * is P(bit == 1) scaled to 1 << 12.
     */
    class _range_enc
    {
    public:
        void init(coder_io *io_ptr)
        {
            io = io_ptr;
            low = 0;
            range = 0xFFFFFFFFu;
            cache = 0;
            cache_size = 1;
        }

        /* p is P(bit == 1), so the lower interval (bit 0) is the complement. */
        inline void encode_bit(int bit, uint32_t p)
        {
            const uint32_t bound = (range >> 12) * (4096u - p);
            if (bit) {
                low += bound;
                range -= bound;
            } else {
                range = bound;
            }
            /* Renormalize once: range stays in [2^24, 2^32), so one shift suffices. */
            if (range < (1u << 24)) {
                range <<= 8;
                shift_low();
            }
        }

        void flush()
        {
            for (int32_t i = 0; i < 5; ++i) {
                shift_low();
            }
        }

    private:
        inline void shift_low()
        {
            if ((uint32_t)(low >> 32) != 0 || low < 0xFF000000ull) {
                const uint8_t carry = (uint8_t)(low >> 32);
                for (;;) {
                    io->putc((uint8_t)(cache + carry));
                    cache = 0xFF;
                    if (--cache_size == 0) {
                        break;
                    }
                }
                cache = (uint8_t)((low >> 24) & 0xFF);
            }
            ++cache_size;
            low = (uint32_t)low << 8;
        }

        coder_io *io;
        uint64_t low;
        uint32_t range;
        uint8_t cache;
        uint64_t cache_size;
    };

    class _range_dec
    {
    public:
        void init(coder_io *io_ptr)
        {
            io = io_ptr;
            range = 0xFFFFFFFFu;
            code = 0;
            for (int32_t i = 0; i < 5; ++i) {
                code = (code << 8) | (uint32_t)io->getc();
            }
        }

        inline int decode_bit(uint32_t p)
        {
            const uint32_t bound = (range >> 12) * (4096u - p);
            int bit;
            if (code < bound) {
                bit = 0;
                range = bound;
            } else {
                bit = 1;
                code -= bound;
                range -= bound;
            }
            /* Renorm, unrolled to the single shift the encoder also performs. */
            if (range < (1u << 24)) {
                range <<= 8;
                code = (code << 8) | (uint32_t)io->getc();
            }
            return bit;
        }

    private:
        coder_io *io;
        uint32_t range;
        uint32_t code;
    };

    /* One binary node: a = observations of bit 0, b = observations of bit 1. */
    struct _pair
    {
        uint16_t a;
        uint16_t b;
    };

    /*
     * Three nested orders over one bit tree: t0 (order 0),
     * t1[ctx1*nodes+node] (order 1), t2[ctx2*nodes+node] (order 2). ctx1/ctx2
     * each have one extra slot for "not enough history yet".
     */
    struct _model
    {
        int32_t bits = 0;      /* bits per symbol */
        int32_t nodes = 0;     /* 1 << bits, index 0 unused */
        int32_t k = 0;         /* alphabet size */
        bool order2 = false;   /* false: previous value only */
        std::vector<_pair> t0;
        std::vector<_pair> t1;
        std::vector<_pair> t2;
        const uint32_t *recip = nullptr; /* reciprocal table, set once in init() */

        void init(int32_t alphabet)
        {
            k = alphabet;
            bits = 0;
            while ((1 << bits) < k) {
                ++bits;
            }
            nodes = 1 << bits;
            order2 = (k <= 128);
            recip = recip_table();
            const _pair zero = {0, 0};
            t0.assign((size_t)nodes, zero);
            t1.assign((size_t)(k + 1) * (size_t)nodes, zero);
            if (order2) {
                t2.assign(((size_t)k * (size_t)k + 1) * (size_t)nodes, zero);
            } else {
                t2.assign((size_t)nodes, zero);
            }
        }

        inline void bump(_pair &p, int bit) const
        {
            if ((uint32_t)p.a + (uint32_t)p.b >= (uint32_t)QCM_CAP) {
                p.a = (uint16_t)((p.a + 1) >> 1);
                p.b = (uint16_t)((p.b + 1) >> 1);
            }
            if (bit) {
                ++p.b;
            } else {
                ++p.a;
            }
        }

        /*
         * Division by a small denominator (bounded by QCM_CAP + QCM_THETA + 2) is
         * the hot operation; denominators become multiplies by a reciprocal.
         */
        inline uint32_t div_fixed(uint64_t num, uint32_t den) const
        {
            return (uint32_t)((num * (uint64_t)recip[den]) >> 32);
        }

        /* P(bit == 1) in 1 << 16, interpolating order 2 towards order 1 towards
           order 0; an empty context reduces to the estimate below it. */
        inline uint32_t prob(const _pair &p0, const _pair &p1, const _pair &p2) const
        {
            const uint32_t n0 = (uint32_t)p0.a + (uint32_t)p0.b;
            uint32_t pr = div_fixed(((uint64_t)p0.b + 1u) << 16, n0 + 2u);
            const uint32_t n1 = (uint32_t)p1.a + (uint32_t)p1.b;
            pr = div_fixed(((uint64_t)p1.b << 16) + (uint64_t)QCM_THETA * (uint64_t)pr,
                           n1 + (uint32_t)QCM_THETA);
            if (order2) {
                const uint32_t n2 = (uint32_t)p2.a + (uint32_t)p2.b;
                pr = div_fixed(((uint64_t)p2.b << 16) + (uint64_t)QCM_THETA * (uint64_t)pr,
                               n2 + (uint32_t)QCM_THETA);
            }
            return pr;
        }

        inline void encode_symbol(_range_enc &rc, uint32_t ctx1, uint32_t ctx2, int32_t sym)
        {
            _pair *r0 = t0.data();
            _pair *r1 = t1.data() + (size_t)ctx1 * (size_t)nodes;
            _pair *r2 = order2 ? (t2.data() + (size_t)ctx2 * (size_t)nodes) : t0.data();
            uint32_t node = 1;
            for (int32_t b = bits - 1; b >= 0; --b) {
                const int bit = (sym >> b) & 1;
                uint32_t p = prob(r0[node], r1[node], r2[node]) >> 4;
                if (p < 1) {
                    p = 1;
                } else if (p > 4095) {
                    p = 4095;
                }
                rc.encode_bit(bit, p);
                bump(r0[node], bit);
                bump(r1[node], bit);
                if (order2) {
                    /* r2 already holds t2.data() + ctx2*nodes, so the per-bit bump reuses
                       it instead of recomputing the ctx2*nodes+node index. */
                    bump(r2[node], bit);
                }
                node = node * 2 + (uint32_t)bit;
            }
        }

        /* Returns the symbol index, or -1 if it is outside the alphabet (only
           reachable from a corrupt stream). */
        inline int32_t decode_symbol(_range_dec &rc, uint32_t ctx1, uint32_t ctx2)
        {
            _pair *r0 = t0.data();
            _pair *r1 = t1.data() + (size_t)ctx1 * (size_t)nodes;
            _pair *r2 = order2 ? (t2.data() + (size_t)ctx2 * (size_t)nodes) : t0.data();
            uint32_t node = 1;
            int32_t sym = 0;
            for (int32_t b = 0; b < bits; ++b) {
                uint32_t p = prob(r0[node], r1[node], r2[node]) >> 4;
                if (p < 1) {
                    p = 1;
                } else if (p > 4095) {
                    p = 4095;
                }
                const int bit = rc.decode_bit(p);
                bump(r0[node], bit);
                bump(r1[node], bit);
                if (order2) {
                    /* r2 already holds t2.data() + ctx2*nodes (see encode_symbol). */
                    bump(r2[node], bit);
                }
                node = node * 2 + (uint32_t)bit;
                sym = (sym << 1) | bit;
            }
            return (sym < k) ? sym : -1;
        }

    private:
        static const size_t k_recip_size = (size_t)QCM_CAP + (size_t)QCM_THETA + 8;

        static const uint32_t *recip_table()
        {
            static const std::vector<uint32_t> tab = []() {
                std::vector<uint32_t> t(k_recip_size, 0);
                for (size_t i = 1; i < k_recip_size; ++i) {
                    t[i] = (uint32_t)(0xFFFFFFFFull / (uint64_t)i);
                }
                return t;
            }();
            return tab.data();
        }
    };

public:
    explicit coder_qcm(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qcm");

        coder_buff = nullptr;
        buff_capacity = 0;
        buff_len = 0;
        flushed = false;
        inited = false;
        k = 0;
        total = 0;
        produced = 0;
        prev1 = -1;
        prev2 = -1;
    }

    ~coder_qcm() override
    {
        if (!flushed && io->m == coder_io::MENC) {
            encode_flush();
        }
        if (coder_buff != nullptr) {
            safe_free((void**)&coder_buff);
        }
    }

    /* Accumulate the column; the block is produced on encode_flush. */
    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        if (io->m != coder_io::MENC) {
            io->m = coder_io::MENC;
            const int32_t initCap = 1 << 20;
            coder_buff = static_cast<uint8_t*>(safe_alloc(initCap));
            check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_qcm: insufficient memory for %d-byte buffer", initCap);
            buff_capacity = initCap;
        }
        if (in_len == 0) {
            return;
        }
        if ((int32_t)(buff_len + in_len) > buff_capacity) {
            const int32_t ncap = (int32_t)((int64_t)buff_capacity * 2 > (int64_t)buff_len + in_len
                                               ? (int64_t)buff_capacity * 2
                                               : (int64_t)buff_len + in_len);
            uint8_t *nb = static_cast<uint8_t*>(safe_alloc(ncap));
            check_exit(nb, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_qcm: insufficient memory for %d-byte buffer", ncap);
            if (buff_len > 0) {
                memcpy(nb, coder_buff, buff_len);
            }
            safe_free((void**)&coder_buff);
            coder_buff = nb;
            buff_capacity = ncap;
        }
        memcpy(coder_buff + buff_len, in, in_len);
        buff_len += (int32_t)in_len;
    }

    /* Commit everything accumulated into one independent coded block. */
    void encode_flush() override
    {
        if (io->m != coder_io::MENC || flushed) {
            return;
        }
        if (buff_len > 0) {
            encode_block();
        }
        flushed = true;
    }

    /* Serve the decoded bytes; the block is decoded as it is asked for. */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX,
                        [[maybe_unused]] bool need2hold = false) override
    {
        if (io->err != coder_io::IO_OK) {
            return coder_ns::CODER_ERR_STREAM_END;
        }
        if (!inited) {
            if (!read_block_header()) {
                return coder_ns::CODER_ERR_INNER;
            }
        }

        int32_t len = 0;
        while (len < (int32_t)out_len && produced < total) {
            const int32_t sym = model.decode_symbol(rcd, ctx1(), ctx2());
            if (sym < 0) {
                coder_logger(coder_ns::ERROR, "coder_qcm: symbol outside the alphabet");
                return coder_ns::CODER_ERR_INNER;
            }
            const uint8_t ch = alpha[sym];
            advance(sym);
            out[len++] = ch;
            ++produced;
            if (split_ch != UINT8_MAX && ch == split_ch) {
                break;
            }
        }
        if (io->err != coder_io::IO_OK) {
            /* The range decoder read past the end of the block: truncated or
               corrupt input, reported rather than silently returning zeros. */
            coder_logger(coder_ns::ERROR, "coder_qcm: stream ended inside the block");
            return coder_ns::CODER_ERR_STREAM_END;
        }
        return len;
    }

    /* Second overload of the base class; not used by this coder. */
    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) override
    {
        return 0;
    }

private:
    void encode_block()
    {
        uint8_t present[256];
        memset(present, 0, sizeof(present));
        for (int32_t i = 0; i < buff_len; ++i) {
            present[coder_buff[i]] = 1;
        }
        int32_t n_alpha = 0;
        for (int32_t v = 0; v < 256; ++v) {
            if (present[v]) {
                alpha[n_alpha++] = (uint8_t)v;
            }
        }

        int32_t index[256];
        for (int32_t v = 0; v < n_alpha; ++v) {
            index[alpha[v]] = v;
        }

        model.init(n_alpha);

        io->putc(kBlockVersion);
        io->putc((uint8_t)(n_alpha & 0xff));
        io->putc((uint8_t)((n_alpha >> 8) & 0xff));
        for (int32_t v = 0; v < n_alpha; ++v) {
            io->putc(alpha[v]);
        }
        put32((uint32_t)buff_len);

        rc.init(io);
        prev1 = -1;
        prev2 = -1;
        for (int32_t i = 0; i < buff_len; ++i) {
            const int32_t sym = index[coder_buff[i]];
            model.encode_symbol(rc, ctx1(), ctx2(), sym);
            advance(sym);
        }
        rc.flush();

        if (io->err != coder_io::IO_OK) {
            coder_logger(coder_ns::ERROR, "coder_qcm: output buffer full during encode");
        }
    }

    bool read_block_header()
    {
        const uint8_t ver = io->getc();
        if (io->err != coder_io::IO_OK) {
            return false;
        }
        if (ver != kBlockVersion) {
            coder_logger(coder_ns::ERROR, "coder_qcm: unsupported block version %u", (unsigned)ver);
            return false;
        }
        const uint32_t n_alpha = (uint32_t)io->getc() | ((uint32_t)io->getc() << 8);
        if (io->err != coder_io::IO_OK) {
            return false;
        }
        if (n_alpha < 1 || n_alpha > 256) {
            coder_logger(coder_ns::ERROR, "coder_qcm: alphabet size %u is invalid", n_alpha);
            return false;
        }
        k = (int32_t)n_alpha;
        for (int32_t v = 0; v < k; ++v) {
            alpha[v] = io->getc();
        }
        total = get32();
        if (io->err != coder_io::IO_OK) {
            return false;
        }
        if (total > (1u << 30)) {
            coder_logger(coder_ns::ERROR, "coder_qcm: block holds %u values, which is invalid", total);
            return false;
        }
        model.init(k);
        rcd.init(io);
        produced = 0;
        prev1 = -1;
        prev2 = -1;
        inited = true;
        return true;
    }

    /* The context indices for the next symbol, in the model's slot layout. */
    inline uint32_t ctx1() const
    {
        return (prev1 >= 0) ? (uint32_t)prev1 : (uint32_t)model.k;
    }

    inline uint32_t ctx2() const
    {
        return (prev1 >= 0 && prev2 >= 0) ? (uint32_t)(prev2 * model.k + prev1)
                                          : (uint32_t)(model.k * model.k);
    }

    inline void advance(int32_t sym)
    {
        prev2 = prev1;
        prev1 = sym;
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

    static const uint8_t kBlockVersion = 1;

    uint8_t *coder_buff;    /* encode accumulation buffer */
    int32_t buff_capacity;
    int32_t buff_len;

    _range_enc rc;          /* encode direction */
    _range_dec rcd;         /* decode direction */
    _model model;

    uint8_t alpha[256];     /* block alphabet, ascending */
    int32_t k;              /* alphabet size */
    uint32_t total;         /* values in the block (decode: to be produced) */
    uint32_t produced;      /* decode: values produced so far */
    int32_t prev1;          /* previous value, -1 before the first */
    int32_t prev2;          /* value before that, -1 when unavailable */
    bool flushed;
    bool inited;
};

#endif
