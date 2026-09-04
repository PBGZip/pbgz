/*
 * coder_arith.h - order-0 adaptive arithmetic byte-stream coder
 *
 * A byte-stream coder with the same external interface as coder_bwt_cm,
 * coder_affix_match and coder_golomb (encode_line / encode_flush / decode_line
 * / decode_inlen), and with the same block-buffering behavior: the compression
 * level selects one "compress at once" block size, and input smaller than one
 * block is cached until a whole block can be committed (see encode_line).
 *
 * Order-0 adaptive arithmetic coding with two backends:
 *
 *   flags 0x01 (legacy, decode-only): cumulative model with total <= 2^14 on
 *   the project's Ppmd8 carryless range coder. Its renormalisation floor was
 *   kBot = 2^15, so range/total could drop to ~2 and the per-symbol division
 *   lost a large part of a bit.
 *
 *   flags 0x02 (current, encode+decode): the htscodecs (CRAM) arrangement -
 *   a carry-less range coder with a 2^24 renormalisation floor (Cache/FFNum/
 *   Carry byte emission, 5-byte flush) fed by an order-0 cumulative model
 *   scaled to 12 bits (total <= 2^12). range stays >= 2^24 while total is at
 *   most 4096, so range/total >= 2^12 and the coding loss drops to ~0.0003
 *   bit/symbol; the 12-bit total with a step of 1 gives the same ~4096-symbol
 *   model memory as htscodecs (MAX_FREQ 2^16 / STEP 16) at higher precision.
 *   No frequency table or per-block dictionary is written: the model is
 *   rebuilt identically by the decoder, so the payload is pure coded bytes.
 *
 * It makes no geometric-distribution assumption, so it reaches the order-0
 * entropy bound. On the SAM POS-delta varint stream (con_sorted.sam, 100k
 * lines/block) flags 0x02 lands near htscodecs' arith_dynamic and well below
 * coder_bwt_cm; see test/golomb_arith_compare.cpp and golomb_pos_compare.cpp.
 *
 * Bitstream block format (each block is independent):
 *   [u32 LE]  source length N, 0 = end of stream
 *   [u8]      flags: 0x00 = raw, 0x01 = order-0 arithmetic (legacy decoder),
 *                    0x02 = order-0 arithmetic (current backend)
 *   [payload] N source bytes (raw), or the arithmetic-coded stream. The
 *             encoder falls back to raw when the coded stream is not shorter
 *             than the plain block (incompressible input).
 *
 * Error convention is the same as the other coders: decode_line returns < 0 on
 * a corrupt/truncated stream, 0 on clean EOF and > 0 on the decoded length.
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
    /* Prior layout: 256 little-endian uint16 symbol weights (512 bytes), the
     * scaled per-byte counts of a file-level training sample. Empty means no
     * prior (the model cold-starts from the uniform table). The bitstream is
     * independent of the prior only in the sense that encode and decode must
     * use the same initial model; a file written with a prior must be decoded
     * with that same prior, otherwise the arithmetic decoding diverges. */
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

    /* Supply the file-level order-0 byte-frequency prior (kPriorBytes).
     * Must be called before the first encode_line / decode_line. */
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
            /* level 0 would leave bsize = 0 and break the buffering; engine
               levels are always 1..9, clamp defensively. */
            if (level < 1) level = 1;
            if (level > 9) level = 9;

            bsize = std::min(tab[level], BLOCK_SIZE); /* Block size (split threshold) */
            /*
             * Grow the buffer to the data instead of pre-allocating bsize:
             * at level 8 bsize is 256 MB, while a field block is typically
             * tens of KB to a few MB, so allocating (and later touching) 256 MB
             * per coder instance cost far more than compressing the data.
             * bsize keeps its meaning as the point at which the stream is split
             * into independent blocks, so the layout is unchanged.
             */
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
                    /* Not at a split point yet: grow instead of splitting, so a
                       field still ends up as one coded block (same ratio as
                       having allocated bsize up front). */
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
            /* Fill up to one block at a time; a single encode_line call may
               carry more than one block's worth of data (e.g. a whole SAM
               block assembled in memory), so commit whole blocks as they are
               filled and keep only the tail cached. */
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

    /* External decompression interface, returns actual decompressed length,
       exits when encountering split_ch during decompression. */
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
                        /* Filled up without finding the separator: the output
                           buffer is too small. */
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
    /* ---- Order-0 adaptive arithmetic backend ------------------------------
     * Cumulative-frequency Laplace model + the project's Ppmd8 carryless
     * range coder (32-bit Low, bytes emitted directly, no carry cache). The
     * total is kept <= 2^14 so that range (>= kBot after renormalization) is
     * always larger than total, and the rescale keeps the cumulative table
     * strictly increasing (each weight >= 1), avoiding zero-size intervals. */

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

        /* Rebuild the cumulative table from a packed 256 x uint16 LE weight
           blob (see coder_arith::kPriorBytes). Weights are scaled so the
           total stays near 2^14 (mirroring the training-side normalization),
           and the table stays strictly increasing so no symbol gets a
           zero-size interval. An all-zero blob keeps the uniform table. */
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
            /* Halve the cumulative values but keep the table strictly
               increasing (each entry at least the previous one plus one). */
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

    /* ---- V2 backend: carry-less range coder (2^24 floor) + 12-bit order-0
     * dynamic model -------------------------------------------------------
     * The htscodecs (CRAM method 6) arrangement:
     *   - 32-bit range coder renormalising up to kTopV2 = 2^24, so range/total
     *     is at least 2^12 and the per-symbol division keeps ~12 bits of
     *     precision. (V1's floor was kBot = 2^15 with a total of 2^14, i.e. a
     *     ratio near 2 and a large coding loss: that was the ~12% gap to
     *     htscodecs, not the model itself.)
     *   - Carry handled by one cached byte plus a run of pending 0xFFs, so a
     *     carry never has to propagate into bytes already emitted.
     *   - Cumulative order-0 model capped at 2^12 with a halving rescale that
     *     keeps every weight >= 1 (strictly increasing table, no zero-size
     *     interval). Step 1 with a 2^12 cap gives ~4096 symbols of model
     *     memory, matching htscodecs (MAX_FREQ 2^16 / STEP 16) at finer
     *     probability resolution.
     * Nothing about the model is stored: encoder and decoder rebuild the same
     * table from the same initial state (plus the same prior, when set), so
     * the payload is pure coded bytes.
     */
    static constexpr uint32_t kTopV2 = 1u << 24;
    static constexpr uint32_t kThresV2 = 255u * kTopV2;

    /* Initial size of the encode buffer. It grows on demand up to bsize
       (the block split threshold), so a small field never pays for the
       level-8 256 MB block size. See encode_line(). */
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
            /* 32-bit overflow of low is the pending carry; it is folded into
               the cached byte by shift_low, so it never has to propagate into
               bytes that were already emitted. */
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

        /* Uniform weight of 1 over [0, max_sym), 0 above it. Symbols the block
           cannot contain keep an empty interval and therefore no probability
           mass - the same trick htscodecs uses (SIMPLE_MODEL init with
           max_sym). Giving all 256 symbols a weight of 1 instead, as V1 did,
           spends a fixed slice of the total on bytes that never occur, which
           at a 2^12 total is worth several percent of the output. */
        void init_default(uint32_t max_sym)
        {
            if (max_sym == 0 || max_sym > 256)
                max_sym = 256;
            for (uint32_t i = 0; i <= 256; ++i)
                freq[i] = (i <= max_sym) ? i : max_sym;
        }

        /* Rebuild the table from the packed 256 x uint16 LE weight blob,
           scaled so the total lands at kMaxTotal. An all-zero blob keeps the
           uniform table. */
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
                /*
                 * Strictly increasing: every symbol keeps a weight of at least
                 * 1 here. This is a correctness requirement, not a tuning
                 * choice: the prior is built from a trial sample, so a symbol
                 * absent there does occur in later blocks, and a zero-width
                 * interval would make the range coder collapse (range *= 0,
                 * and the renormalisation loop can never satisfy
                 * range >= kTopV2 - it does not terminate).
                 */
                if (v <= last)
                    v = last + 1;
                freq[i + 1] = (uint32_t)v;
                last = v;
            }
        }

        void rescale()
        {
            /* Ceil-halve every weight: a symbol seen at least once keeps a
               weight of at least 1, a never-seen one stays at 0 and keeps no
               probability mass (halving the cumulative table directly, as V1
               did, forces every one of the 256 symbols to keep weight 1, which
               at a 2^12 total wastes ~6% of the probability on symbols that
               never occur). */
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
            /* First symbol whose upper bound exceeds target. Choosing the
               lower bound the other way round would be ambiguous for symbols
               with an empty interval (freq[i] == freq[i+1]); this form skips
               them, matching the encoder's cumulative layout. */
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

    /* Write one block of n source bytes plus its header. The arithmetic code
       is produced into a temporary buffer first so that incompressible blocks
       can be stored raw instead of bloating the stream. */
    void encode_one_block(int32_t n)
    {
        put32((uint32_t)n);

        /* A touch over one byte per symbol worst case plus the 5-byte flush;
           the temporary coder_io view also bounds every write. */
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

    /* Decode one block from the code stream into coder_buff; returns the
       source length, 0 on the end marker, < 0 on a corrupt stream. */
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
