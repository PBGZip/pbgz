/*
 * coder_rans.h - static order-0 entropy byte-stream coder for the fast path
 *
 * Same external interface and buffering behaviour as coder_bwt_cm /
 * coder_arith (encode_line / encode_flush / decode_line / decode_inlen):
 * the caller caches the whole column with encode_line and commits it with a
 * single encode_flush, producing one independent coded block.
 *
 * Why this exists: the fast preset fixes one coder per column and skips
 * preselection, so per-symbol speed is what matters. coder_arith rebuilds its
 * order-0 adaptive model per symbol by adding a step to every suffix of the
 * cumulative table, i.e. O(256) work per byte - fine for archive ratio, not a
 * fast coder. This coder instead fits one static order-0 frequency model per
 * block (a single O(n) counting pass plus an O(256) normalisation) and then
 * codes every symbol in O(1): the entropy backend is the carry-less range
 * coder used by coder_arith's current (htscodecs-style) backend, with all
 * per-symbol model updates removed. The model is fully described by the
 * 256 x uint16 weight table written in the block header, so the decoder needs
 * no prior state, no side information and no alphabet inference: decode is
 * deterministic from the block alone.
 *
 * Bitstream block format (each block is independent):
 *   [u32 LE] source length N, 0 = end of stream
 *   [u8]     flags: 0x00 = raw, 0x01 = static order-0 coded
 *   [raw]    N source bytes                                  (flags 0x00)
 *   or
 *   [512B]   packed 256 x uint16 LE symbol weights           (flags 0x01)
 *   [payload] the entropy-coded stream (carry-less range coder, 5-byte flush)
 *
 * The encoder falls back to raw when the coded stream is not shorter than the
 * plain block (incompressible input).
 *
 * Error convention is the same as the other coders: decode_line returns < 0 on
 * a corrupt/truncated stream, 0 on clean EOF and > 0 on the decoded length.
 */

#ifndef _CODER_RANS_H_
#define _CODER_RANS_H_

#include <stdint.h>
#include <algorithm>
#include <cstring>
#include <vector>

#include "coder_io.h"
#include "coder.h"

class coder_rans : public coder
{
public:
    static constexpr uint32_t kModelTotal = 1u << 16; /* cumulative total the weights are scaled to */
    static constexpr uint32_t kTableBytes = 256 * 2;  /* 256 x uint16 LE weights */

    coder_rans(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_rans");

        coder_buff = nullptr;
        buff_capacity = 0;
        buff_len = 0;
        flushed = false;
    }

    virtual ~coder_rans()
    {
        if (!flushed && io->m == coder_io::MENC)
            encode_flush();
        if (coder_buff)
            safe_free((void**)&coder_buff);
    }

    /* External compression interface: accumulate the whole field. */
    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        if (io->m != coder_io::MENC)
        {
            io->m = coder_io::MENC;
            const int32_t initCap = 1 << 20; /* grow on demand */
            coder_buff = static_cast<uint8_t*>(safe_alloc(initCap));
            check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_rans: insufficient memory for %d-byte buffer", initCap);
            buff_capacity = initCap;
        }
        if (in_len == 0)
            return;
        if ((int32_t)(buff_len + in_len) > buff_capacity)
        {
            const int32_t ncap = (int32_t)std::max((int64_t)buff_capacity * 2,
                                                   (int64_t)buff_len + in_len);
            uint8_t *nb = static_cast<uint8_t*>(safe_alloc(ncap));
            check_exit(nb, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_rans: insufficient memory for %d-byte buffer", ncap);
            if (buff_len > 0)
                memcpy(nb, coder_buff, buff_len);
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
        if (io->m != coder_io::MENC || flushed)
            return;
        encode_one_block(buff_len);
        put32(0); /* end-of-stream marker */
        flushed = true;
    }

    /* External decompression interface; returns decoded length. */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX,
                        bool need2hold [[maybe_unused]] = false) override
    {
        int32_t len = 0;

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
                if (out_off == blk_len)
                {
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        return len;
                }
                for (int32_t n = out_off; n < blk_len; n++)
                {
                    uint8_t ch = coder_buff[n];
                    *(out + len++) = ch;
                    if (ch == split_ch)
                    {
                        out_off = ++n;
                        return len;
                    }
                    if (len >= (int32_t)out_len)
                    {
                        out_off = ++n;
                        io->set_err(coder_io::IO_BUF_FULL);
                        return coder_ns::CODER_ERR_BUF_SMALL;
                    }
                }
                out_off = blk_len;
            }
        }

        for (;;)
        {
            if (out_off == blk_len)
            {
                int32_t ret = decode_one_block();
                if (ret < 0)
                    return ret;
                if (ret == 0)
                    return len;
            }
            for (int32_t n = out_off; n < blk_len; n++)
            {
                *(out + len++) = coder_buff[n];
                if (len >= (int32_t)out_len)
                {
                    out_off = ++n;
                    return len;
                }
            }
            out_off = blk_len;
        }
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool) override
    {
        return 0;
    }

    int32_t decode_inlen()
    {
        return this->io->data_len;
    }

private:
    /* ---- Carry-less range coder (htscodecs/CRAM-style arrangement) -----
     * 32-bit range with a 2^24 renormalisation floor, a cached byte plus a
     * run of pending 0xFFs for carries, and a 5-byte flush. The model total
     * stays at 2^16, so range/total >= 2^8 at every symbol: the per-symbol
     * 32x32 divide keeps enough precision and there is no coding-loss trap. */
    static constexpr uint32_t kTop = 1u << 24;
    static constexpr uint32_t kThres = 255u * kTop;

    struct RcEnc {
        coder_io *io;
        uint32_t low;
        uint32_t range;
        uint32_t ffNum; /* pending 0xFF bytes */
        uint32_t cache; /* pending top byte, not yet emitted */
        uint32_t carry; /* carry into cache */

        explicit RcEnc(coder_io *i)
            : io(i), low(0), range(0xFFFFFFFFu), ffNum(0), cache(0), carry(0) {}

        void shift_low()
        {
            if (low < kThres || carry) {
                io->putc((uint8_t)(cache + carry));
                while (ffNum) {
                    io->putc((uint8_t)(carry - 1));
                    --ffNum;
                }
                cache = low >> 24;
                carry = 0;
            } else {
                ++ffNum;
            }
            low <<= 8;
        }

        void encode(uint32_t cum, uint32_t freq, uint32_t total)
        {
            const uint32_t prev = low;
            low += cum * (range /= total);
            range *= freq;
            carry += (low < prev) ? 1u : 0u;
            while (range < kTop) {
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

    struct RcDec {
        coder_io *io;
        uint32_t range;
        uint32_t code;

        explicit RcDec(coder_io *i) : io(i), range(0xFFFFFFFFu), code(0) {}

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

        void decode(uint32_t cum, uint32_t freq)
        {
            code -= cum * range;
            range *= freq;
            while (range < kTop) {
                code = (code << 8) | (uint32_t)io->getc();
                range <<= 8;
            }
        }
    };

    /* Static order-0 cumulative model: freq[256] is the total, and
       init_from_weights rebuilds the exact same table on both sides. */
    struct StaticModel {
        uint32_t freq[257]; /* cumulative: freq[0]=0, freq[256]=total */

        StaticModel() { reset_uniform(); }

        void reset_uniform()
        {
            for (int i = 0; i <= 256; ++i)
                freq[i] = (uint32_t)i;
        }

        /* Rebuild from the packed 256 x uint16 LE weight blob, scaled so the
           total lands at kModelTotal and every weight is >= 1 (strictly
           increasing table, no zero-size interval). An all-zero blob keeps
           the uniform table. */
        void init_from_weights(const uint8_t *w)
        {
            uint32_t wt[256];
            uint64_t total = 0;
            for (int i = 0; i < 256; ++i) {
                wt[i] = (uint32_t)w[2 * i] | ((uint32_t)w[2 * i + 1] << 8);
                total += wt[i];
            }
            if (total == 0) {
                reset_uniform();
                return;
            }
            const uint64_t target = kModelTotal;
            uint64_t acc = 0, last = 0;
            freq[0] = 0;
            for (int i = 0; i < 256; ++i) {
                acc += wt[i];
                uint64_t v = (acc * target) / total;
                if (v <= last)
                    v = last + 1;
                freq[i + 1] = (uint32_t)v;
                last = v;
            }
        }

        void encode(RcEnc &rc, uint8_t sym)
        {
            const uint32_t total = freq[256];
            rc.encode(freq[sym], freq[sym + 1] - freq[sym], total);
        }

        uint8_t decode(RcDec &rc)
        {
            const uint32_t total = freq[256];
            const uint32_t target = rc.get_freq(total);
            int lo = 0, hi = 255;
            while (lo < hi) {
                const int mid = (lo + hi) >> 1;
                if (freq[mid + 1] > target)
                    hi = mid;
                else
                    lo = mid + 1;
            }
            const uint8_t sym = (uint8_t)lo;
            rc.decode(freq[sym], freq[sym + 1] - freq[sym]);
            return sym;
        }
    };


    /* Sparse order-1 header encoding (shared with the range order-1 path and
       the rANS order-1 path): context = previous byte. */
    static size_t makeO1Header(const uint32_t cnt[256][256], uint8_t* hdr, size_t cap)
    {
        size_t p = 0;
        for (int ctx = 0; ctx < 256; ++ctx) {
            uint32_t n = 0;
            for (int sym = 0; sym < 256; ++sym)
                if (cnt[ctx][sym] > 0)
                    ++n;
            if (p + 1 + (size_t)n * 3 > cap)
                break;
            hdr[p++] = (uint8_t)n;
            for (int sym = 0; sym < 256; ++sym) {
                if (cnt[ctx][sym] == 0)
                    continue;
                uint32_t w = cnt[ctx][sym];
                if (w > 0xFFFF)
                    w = 0xFFFF;
                hdr[p++] = (uint8_t)sym;
                hdr[p++] = (uint8_t)(w & 0xff);
                hdr[p++] = (uint8_t)(w >> 8);
            }
        }
        return p;
    }

    /* Legacy range-coder order-1 model, kept to decode files whose flags are
       0x02 (produced before the rANS backend replaced the range coder). */
    struct Order1Model {
        std::vector<uint32_t> freq;   /* [ctx*257+sym] cumulative, total ~2^16 */

        bool init_from_header(const uint8_t* h, size_t len)
        {
            freq.assign(256 * 257, 0);
            size_t p = 0;
            for (int ctx = 0; ctx < 256; ++ctx) {
                if (p >= len)
                    return false;
                const uint32_t n = h[p++];
                if (n == 0)
                    continue;
                if (n > 256 || p + (size_t)n * 3 > len)
                    return false;
                uint32_t w[256];
                for (uint32_t i = 0; i < 256; ++i)
                    w[i] = 0;
                uint64_t total = 0;
                for (uint32_t i = 0; i < n; ++i) {
                    const uint32_t sym = h[p];
                    const uint32_t wt = h[p + 1] | ((uint32_t)h[p + 2] << 8);
                    p += 3;
                    if (sym >= 256)
                        return false;
                    w[sym] = wt;
                    total += wt;
                }
                if (total == 0)
                    continue;
                uint32_t* f = &freq[(size_t)ctx * 257];
                f[0] = 0;
                uint64_t acc = 0, last = 0;
                for (int sym = 0; sym < 256; ++sym) {
                    acc += w[sym];
                    uint64_t v = (acc * kModelTotal) / total;
                    if (v <= last)
                        v = last + 1;
                    f[sym + 1] = (uint32_t)v;
                    last = v;
                }
            }
            return true;
        }

        uint8_t decode(RcDec& rc, uint8_t ctx)
        {
            const uint32_t* f = &freq[(size_t)ctx * 257];
            const uint32_t total = f[256];
            const uint32_t target = rc.get_freq(total);
            int lo = 0, hi = 255;
            while (lo < hi) {
                const int mid = (lo + hi) >> 1;
                if (f[mid + 1] > target)
                    hi = mid;
                else
                    lo = mid + 1;
            }
            const uint8_t sym = (uint8_t)lo;
            rc.decode(f[sym], f[sym + 1] - f[sym]);
            return sym;
        }
    };

    /* ---- ryg-style rANS (order-0 / order-1) ------------------------------
     * 32-bit state, L = 2^23, 8-bit renormalisation, frequency scale 2^12.
     * Encoding runs symbols in reverse order (last to first) writing bytes to
     * the end of a temp buffer; the final state is flushed as 4 LE bytes in
     * front of the renorm bytes. Decoding is forward.
     *   put : renorm while x >= x_max; x = (x/freq)<<12 + x%freq + cum
     *   get : cum = x & 4095
     *   adv : x = freq*(x>>12) + (x&4095) - cum; renorm by reading bytes
     * Reference: Fabian 'ryg' Giesen, rans_byte.h (public domain).
     */
    static constexpr uint32_t kRansL = 1u << 23;
    static constexpr uint32_t kRansScaleBits = 12;
    static constexpr uint32_t kRansMask = (1u << kRansScaleBits) - 1;

    /* Build a 4096-total cumulative table from raw weights: symbols with a
       zero weight keep an empty interval (0 width), present symbols are scaled
       so the total is exactly 4096 and every interval is >= 1. Both encoder
       and decoder call this on the same stored weights, so it is the single
       source of truth for the model. */
    static void cum4096(const uint32_t *w, uint32_t f[257])
    {
        /* Every present symbol keeps an interval of width >= 1 and the widths
           sum to exactly 2^12 (rANS requires total == scale). Naive last+1
           scaling can overshoot 4096 when many rare symbols are forced to
           width 1, which breaks the x&mask algebra; give each present symbol a
           base width of 1 first and distribute the remaining budget by
           weighted largest remainders. */
        uint64_t tot = 0;
        uint32_t used[256];
        uint32_t base[256] = {0};
        uint32_t usedCount = 0;
        for (int i = 0; i < 256; ++i) {
            if (w[i] == 0)
                continue;
            tot += w[i];
            used[usedCount++] = (uint32_t)i;
        }
        if (tot == 0) {
            for (int i = 1; i <= 256; ++i)
                f[i] = (uint32_t)i;   /* uniform fallback */
            return;
        }
        const uint32_t extra = (1u << kRansScaleBits) - usedCount;
        uint32_t extra_i[256] = {0};
        uint64_t frac[256] = {0};
        {
            uint64_t sum = 0;
            for (uint32_t k = 0; k < usedCount; ++k) {
                const uint32_t i = used[k];
                extra_i[k] = (uint32_t)(((uint64_t)w[i] * extra) / tot);
                sum += extra_i[k];
                frac[k] = ((uint64_t)w[i] * extra) % tot;
            }
            uint32_t deficit = extra - (uint32_t)sum;
            while (deficit-- > 0) {
                uint32_t best = 0;
                for (uint32_t k = 1; k < usedCount; ++k)
                    if (frac[k] > frac[best])
                        best = k;
                extra_i[best]++;
                frac[best] = 0;
            }
        }
        for (uint32_t k = 0; k < usedCount; ++k)
            base[used[k]] = 1 + extra_i[k];
        f[0] = 0;
        for (int i = 0; i < 256; ++i)
            f[i + 1] = f[i] + base[i];
    }

    /* order-0 model: f is the 4096-total cumulative table. */
    struct RansModel0 {
        uint32_t f[257];

        void init(const uint8_t *wtBlob)   /* 256 x u16 LE */
        {
            uint32_t w[256];
            for (int i = 0; i < 256; ++i)
                w[i] = (uint32_t)wtBlob[2 * i] | ((uint32_t)wtBlob[2 * i + 1] << 8);
            cum4096(w, f);
        }
    };

    /* order-1 model: sparse per-context headers (same layout as the range
       order-1 header), rebuilt into 4096-total cumulative per context. */
    struct RansModel1 {
        std::vector<uint32_t> f;   /* [ctx*257 + sym] cumulative */

        bool init(const uint8_t *h, size_t len)
        {
            f.assign((size_t)256 * 257, 0);
            size_t p = 0;
            for (int ctx = 0; ctx < 256; ++ctx) {
                if (p >= len)
                    return false;
                const uint32_t n = h[p++];
                if (n == 0)
                    continue;
                if (n > 256 || p + (size_t)n * 3 > len)
                    return false;
                uint32_t w[256];
                for (uint32_t i = 0; i < 256; ++i)
                    w[i] = 0;
                for (uint32_t i = 0; i < n; ++i) {
                    const uint32_t sym = h[p];
                    const uint32_t wt = h[p + 1] | ((uint32_t)h[p + 2] << 8);
                    p += 3;
                    if (sym >= 256)
                        return false;
                    w[sym] = wt;
                }
                uint32_t *fctx = &f[(size_t)ctx * 257];
                cum4096(w, fctx);
            }
            return true;
        }
    };

    /* Four-way interleaved rANS: symbol i goes to lane i&3; symbols of a lane
       keep their original relative order, each lane is an independent rANS32
       stream. Payload layout (as written/read by the caller):
         [u32 LE lane0 len][lane1][lane2][lane3] then the four lane payloads
       concatenated lane-major, each starting with its 4-byte flush state.
       order1 uses the per-context models (ctx = the byte preceding the symbol
       in the ORIGINAL stream); order0 uses the single global table. */
    void ransEncode4(bool order1, int32_t n, const RansModel0 &m0,
                     const RansModel1 &m1, uint32_t lens[4],
                     std::vector<uint8_t> &payload)
    {
        std::vector<uint8_t> lsym[4], lctx[4];
        for (int32_t i = 0; i < n; ++i) {
            const uint8_t s = coder_buff[i];
            const uint8_t c = (i > 0) ? coder_buff[i - 1] : 0;
            lsym[i & 3].push_back(s);
            lctx[i & 3].push_back(c);
        }
        payload.clear();
        payload.reserve((size_t)n * 2 + 32);
        for (int j = 0; j < 4; ++j) {
            const size_t m = lsym[j].size();
            if (m == 0) {
                lens[j] = 0;
                continue;
            }
            std::vector<uint8_t> cap(m * 2 + 64);
            uint8_t *base = cap.data();
            uint8_t *ptr = base + cap.size();
            uint32_t x = kRansL;
            for (size_t k = m; k-- > 0;) {
                const uint8_t sym = lsym[j][k];
                uint32_t st, fr;
                if (order1) {
                    const uint32_t ctx = lctx[j][k];
                    const uint32_t *f = &m1.f[(size_t)ctx * 257];
                    st = f[sym];
                    fr = f[sym + 1] - st;
                } else {
                    st = m0.f[sym];
                    fr = m0.f[sym + 1] - st;
                }
                const uint32_t xmax = ((kRansL >> kRansScaleBits) << 8) * fr;
                while (x >= xmax) {
                    *--ptr = (uint8_t)(x & 0xff);
                    x >>= 8;
                }
                x = ((x / fr) << kRansScaleBits) + (x % fr) + st;
            }
            ptr -= 4;
            ptr[0] = (uint8_t)(x & 0xff);
            ptr[1] = (uint8_t)((x >> 8) & 0xff);
            ptr[2] = (uint8_t)((x >> 16) & 0xff);
            ptr[3] = (uint8_t)((x >> 24) & 0xff);
            lens[j] = (uint32_t)(base + cap.size() - ptr);
            payload.insert(payload.end(), ptr, base + cap.size());
        }
    }

    /* Encode one order-0/order-1 run into `out` (prefix bytes). Returns the
       number of coded bytes; the caller then compares with raw size. */
    size_t ransEncode(const uint32_t *cum /* [s+1]=upper, cum[0]=0 */, size_t n,
                      std::vector<uint8_t> &out) const
    {
        const size_t cap = n * 2 + 64;
        out.resize(cap);
        uint8_t *base = out.data();
        uint8_t *ptr = base + cap;
        uint32_t x = kRansL;
        for (size_t i = n; i-- > 0;) {
            const uint8_t sym = coder_buff[i];
            const uint32_t st = cum[sym];
            const uint32_t fr = cum[sym + 1] - st;
            const uint32_t xmax = ((kRansL >> kRansScaleBits) << 8) * fr;
            while (x >= xmax) {
                *--ptr = (uint8_t)(x & 0xff);
                x >>= 8;
            }
            x = ((x / fr) << kRansScaleBits) + (x % fr) + st;
        }
        /* flush final state: 4 LE bytes before the renorm bytes */
        ptr -= 4;
        ptr[0] = (uint8_t)(x & 0xff);
        ptr[1] = (uint8_t)((x >> 8) & 0xff);
        ptr[2] = (uint8_t)((x >> 16) & 0xff);
        ptr[3] = (uint8_t)((x >> 24) & 0xff);
        const size_t len = (size_t)(base + cap - ptr);
        memmove(base, ptr, len);
        out.resize(len);
        return len;
    }

    /* One block of n source bytes plus its header. Candidates: order-0 rANS,
       order-1 rANS, raw - the shortest is written (old range-encoded files
       keep decoding through their own flags). */
    void encode_one_block(int32_t n)
    {
        put32((uint32_t)n);
        if (n == 0)
            return;

        uint32_t cnt0[256] = {0};
        for (int32_t i = 0; i < n; i++)
            cnt0[coder_buff[i]]++;
        uint64_t cnt0Total = 0;
        for (int i = 0; i < 256; i++)
            cnt0Total += cnt0[i];

        /* order-0 weight blob (512 B, self-describing) */
        uint8_t wt0[kTableBytes];
        if (cnt0Total > 0) {
            for (int i = 0; i < 256; i++) {
                uint32_t w = 0;
                if (cnt0[i] > 0) {
                    uint64_t s = ((uint64_t)cnt0[i] * kModelTotal) / cnt0Total;
                    w = (uint32_t)(s > 0 ? s : 1);
                    if (w > 0xFFFF)
                        w = 0xFFFF;
                }
                wt0[2 * i] = (uint8_t)(w & 0xff);
                wt0[2 * i + 1] = (uint8_t)(w >> 8);
            }
        } else {
            memset(wt0, 0, sizeof(wt0));
        }
        RansModel0 m0;
        m0.init(wt0);

        /* order-1 sparse header + model (context = previous byte) */
        uint32_t cnt1[256][256] = {{0}};
        {
            uint8_t prev = 0;
            for (int32_t i = 0; i < n; i++) {
                cnt1[prev][coder_buff[i]]++;
                prev = coder_buff[i];
            }
        }
        std::vector<uint8_t> hdr1((size_t)256 * 769);
        const size_t hdr1Len = makeO1Header(cnt1, hdr1.data(), hdr1.size());
        RansModel1 m1;
        const bool o1ok = m1.init(hdr1.data(), hdr1Len);

        /* order choice and coding: encode BOTH 4-way candidates and take the
           shorter (4-way == single-stream in size; it only parallelises the
           state machine, so single-stream candidates are not produced). */
        std::vector<uint8_t> p0, p1;
        uint32_t len4o0[4], len4o1[4];
        ransEncode4(false, n, m0, m1, len4o0, p0);
        const int64_t lenO0_4T = (int64_t)p0.size() + 16 + (int64_t)kTableBytes;
        int64_t lenO1_4T = INT64_MAX;
        if (o1ok && hdr1Len > 0) {
            ransEncode4(true, n, m0, m1, len4o1, p1);
            lenO1_4T = (int64_t)p1.size() + 16 + (int64_t)hdr1Len;
        }
        const int64_t lenRaw = n;

        if (lenO1_4T < lenO0_4T && lenO1_4T < lenRaw) {
            io->putc(0x06);
            for (size_t i = 0; i < hdr1Len; i++)
                io->putc(hdr1[i]);
            for (int j = 0; j < 4; ++j) {
                io->putc((uint8_t)(len4o1[j] & 0xff));
                io->putc((uint8_t)((len4o1[j] >> 8) & 0xff));
                io->putc((uint8_t)((len4o1[j] >> 16) & 0xff));
                io->putc((uint8_t)((len4o1[j] >> 24) & 0xff));
            }
            for (size_t i = 0; i < p1.size(); i++)
                io->putc(p1[i]);
            return;
        }
        if (lenO0_4T < lenRaw) {
            io->putc(0x05);
            for (uint32_t i = 0; i < kTableBytes; i++)
                io->putc(wt0[i]);
            for (int j = 0; j < 4; ++j) {
                io->putc((uint8_t)(len4o0[j] & 0xff));
                io->putc((uint8_t)((len4o0[j] >> 8) & 0xff));
                io->putc((uint8_t)((len4o0[j] >> 16) & 0xff));
                io->putc((uint8_t)((len4o0[j] >> 24) & 0xff));
            }
            for (size_t i = 0; i < p0.size(); i++)
                io->putc(p0[i]);
            return;
        }
        io->putc(0x00);
        for (int32_t i = 0; i < n; i++)
            io->putc(coder_buff[i]);
    }

    /* Decode one block from the code stream into coder_buff; returns the
       source length, 0 on the end marker, < 0 on a corrupt stream. */
    int32_t decode_one_block()
    {
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        int32_t n = (int32_t)get32();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        if (n == 0)
            return 0;
        if (n < 0 || n > (1 << 30))
        {
            coder_logger(coder_ns::ERROR, "coder_rans: block size %d is invalid", n);
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
        {
            uint8_t wt[kTableBytes];
            for (uint32_t i = 0; i < kTableBytes; i++)
                wt[i] = io->getc();
            if (io->err != coder_io::IO_OK)
                return coder_ns::CODER_ERR_STREAM_END;
            RcDec dec(io);
            dec.init();
            StaticModel model;
            model.init_from_weights(wt);
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = model.decode(dec);
        }
        else if (flags == 0x02)
        {
            std::vector<uint8_t> hdr1((size_t)256 * 769);
            size_t hp = 0;
            for (int ctx = 0; ctx < 256; ++ctx) {
                const int ch = io->getc();
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
                const uint32_t cnt = (uint32_t)ch;
                hdr1[hp++] = (uint8_t)cnt;
                for (uint32_t k = 0; k < cnt; ++k) {
                    for (int b = 0; b < 3; ++b) {
                        const int cb = io->getc();
                        if (io->err != coder_io::IO_OK)
                            return coder_ns::CODER_ERR_STREAM_END;
                        hdr1[hp++] = (uint8_t)cb;
                    }
                }
            }
            RcDec dec(io);
            dec.init();
            Order1Model model;
            if (!model.init_from_header(hdr1.data(), hp)) {
                coder_logger(coder_ns::ERROR, "coder_rans: malformed order-1 header");
                return coder_ns::CODER_ERR_INNER;
            }
            uint8_t prev = 0;
            for (int32_t i = 0; i < n; i++) {
                coder_buff[i] = model.decode(dec, prev);
                prev = coder_buff[i];
            }
        }
        else if (flags == 0x03)
        {
            uint8_t wt[kTableBytes];
            for (uint32_t i = 0; i < kTableBytes; i++)
                wt[i] = io->getc();
            if (io->err != coder_io::IO_OK)
                return coder_ns::CODER_ERR_STREAM_END;
            RansModel0 m0;
            m0.init(wt);
            uint32_t x = 0;
            for (int i = 0; i < 4; ++i)
                x |= (uint32_t)io->getc() << (8 * i);
            for (int32_t i = 0; i < n; i++) {
                const uint32_t get = x & kRansMask;
                int lo = 0, hi = 255;
                while (lo < hi) {
                    const int mid = (lo + hi) >> 1;
                    if (m0.f[mid + 1] > get)
                        hi = mid;
                    else
                        lo = mid + 1;
                }
                const uint8_t sym = (uint8_t)lo;
                const uint32_t st = m0.f[sym];
                const uint32_t fr = m0.f[sym + 1] - st;
                x = fr * (x >> kRansScaleBits) + (x & kRansMask) - st;
                while (x < kRansL) {
                    x = (x << 8) | (uint32_t)io->getc();
                    if (io->err != coder_io::IO_OK)
                        return coder_ns::CODER_ERR_STREAM_END;
                }
                coder_buff[i] = sym;
            }
        }
        else if (flags == 0x04)
        {
            std::vector<uint8_t> hdr1((size_t)256 * 769);
            size_t hp = 0;
            for (int ctx = 0; ctx < 256; ++ctx) {
                const int ch = io->getc();
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
                const uint32_t cnt = (uint32_t)ch;
                hdr1[hp++] = (uint8_t)cnt;
                for (uint32_t k = 0; k < cnt; ++k) {
                    for (int b2 = 0; b2 < 3; ++b2) {
                        const int cb = io->getc();
                        if (io->err != coder_io::IO_OK)
                            return coder_ns::CODER_ERR_STREAM_END;
                        hdr1[hp++] = (uint8_t)cb;
                    }
                }
            }
            RansModel1 m1;
            if (!m1.init(hdr1.data(), hp)) {
                coder_logger(coder_ns::ERROR, "coder_rans: malformed order-1 header");
                return coder_ns::CODER_ERR_INNER;
            }
            uint32_t x = 0;
            for (int i = 0; i < 4; ++i)
                x |= (uint32_t)io->getc() << (8 * i);
            uint8_t prev = 0;
            for (int32_t i = 0; i < n; i++) {
                const uint32_t *fctx = &m1.f[(size_t)prev * 257];
                const uint32_t get = x & kRansMask;
                int lo = 0, hi = 255;
                while (lo < hi) {
                    const int mid = (lo + hi) >> 1;
                    if (fctx[mid + 1] > get)
                        hi = mid;
                    else
                        lo = mid + 1;
                }
                const uint8_t sym = (uint8_t)lo;
                const uint32_t st = fctx[sym];
                const uint32_t fr = fctx[sym + 1] - st;
                x = fr * (x >> kRansScaleBits) + (x & kRansMask) - st;
                while (x < kRansL)
                    x = (x << 8) | (uint32_t)io->getc();
                coder_buff[i] = sym;
                prev = sym;
            }
        }
        else if (flags == 0x05 || flags == 0x06)
        {
            /* four-way interleaved rANS. */
            RansModel0 m0;
            RansModel1 m1;
            if (flags == 0x06) {
                std::vector<uint8_t> hdr1((size_t)256 * 769);
                size_t hp = 0;
                for (int ctx = 0; ctx < 256; ++ctx) {
                    const int ch = io->getc();
                    if (io->err != coder_io::IO_OK)
                        return coder_ns::CODER_ERR_STREAM_END;
                    const uint32_t cnt = (uint32_t)ch;
                    hdr1[hp++] = (uint8_t)cnt;
                    for (uint32_t k = 0; k < cnt; ++k) {
                        for (int b2 = 0; b2 < 3; ++b2) {
                            const int cb = io->getc();
                            if (io->err != coder_io::IO_OK)
                                return coder_ns::CODER_ERR_STREAM_END;
                            hdr1[hp++] = (uint8_t)cb;
                        }
                    }
                }
                if (!m1.init(hdr1.data(), hp)) {
                    coder_logger(coder_ns::ERROR, "coder_rans: malformed order-1 header");
                    return coder_ns::CODER_ERR_INNER;
                }
            } else {
                uint8_t wt[kTableBytes];
                for (uint32_t i = 0; i < kTableBytes; i++)
                    wt[i] = io->getc();
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
                m0.init(wt);
            }
            uint32_t len4[4];
            for (int j = 0; j < 4; ++j) {
                uint32_t v = 0;
                for (int b2 = 0; b2 < 4; ++b2)
                    v |= (uint32_t)io->getc() << (8 * b2);
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
                len4[j] = v;
            }
            const uint64_t total = (uint64_t)len4[0] + len4[1] + len4[2] + len4[3];
            if (total > (uint64_t)(1u << 30))
                return coder_ns::CODER_ERR_INNER;
            std::vector<uint8_t> pay((size_t)total);
            for (size_t i = 0; i < (size_t)total; ++i) {
                pay[i] = io->getc();
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
            }
            uint32_t x4[4];
            size_t pos4[4], end4[4];
            size_t off = 0;
            for (int j = 0; j < 4; ++j) {
                if (len4[j] < 4)
                    return coder_ns::CODER_ERR_INNER;
                x4[j] = (uint32_t)pay[off] | ((uint32_t)pay[off + 1] << 8) |
                        ((uint32_t)pay[off + 2] << 16) | ((uint32_t)pay[off + 3] << 24);
                pos4[j] = off + 4;
                end4[j] = off + len4[j];
                off += len4[j];
            }
            uint8_t prev = 0;
            for (int32_t i = 0; i < n; ++i) {
                const int lane = i & 3;
                const uint32_t *f = (flags == 0x06) ? &m1.f[(size_t)prev * 257] : m0.f;
                uint32_t x = x4[lane];
                const uint32_t g = x & kRansMask;
                int lo = 0, hi = 255;
                while (lo < hi) {
                    const int mid = (lo + hi) >> 1;
                    if (f[mid + 1] > g)
                        hi = mid;
                    else
                        lo = mid + 1;
                }
                const uint8_t sym = (uint8_t)lo;
                const uint32_t st = f[sym];
                const uint32_t fr = f[sym + 1] - st;
                x = fr * (x >> kRansScaleBits) + (x & kRansMask) - st;
                while (x < kRansL) {
                    if (pos4[lane] >= end4[lane])
                        return coder_ns::CODER_ERR_STREAM_END;
                    x = (x << 8) | pay[pos4[lane]++];
                }
                x4[lane] = x;
                coder_buff[i] = sym;
                prev = sym;
            }
        }
        else
        {
            coder_logger(coder_ns::ERROR, "coder_rans: invalid block flags %d", flags);
            return coder_ns::CODER_ERR_INNER;
        }

        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        blk_len = n;
        out_off = 0;
        return n;
    }

    void ensure_buffer(uint32_t need)
    {
        if (need <= (uint32_t)buff_capacity)
            return;
        uint8_t *nb = static_cast<uint8_t*>(safe_alloc(need));
        if (!nb)
        {
            coder_logger(coder_ns::ERROR, "coder_rans: insufficient memory for %u-byte block buffer", need);
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
    uint8_t *coder_buff;  /* Accumulated encode input / decoded block output */
    int32_t buff_capacity; /* Capacity of coder_buff */
    int32_t buff_len;      /* Bytes currently buffered (encoding) */
    int32_t blk_len;       /* Source length of the current decoded block */
    int32_t out_off;       /* Read offset inside the current decoded block */
    bool flushed;
};
#endif
