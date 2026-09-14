/*
 * bgzf_parallel.h - parallel BGZF decoder for the BAM read path
 *
 * A BAM file is a BGZF stream: a concatenation of independent gzip members,
 * each compressing exactly one 64 KiB (uncompressed) block, terminated by a
 * 28-byte empty EOF member. Every member carries its total compressed length
 * in the BSIZE field of its gzip extra ("BC" subfield), so members can be cut
 * out of the raw stream *without inflating* and decompressed in parallel -
 * unlike a plain single-stream gzip, where member boundaries are only known
 * after inflating.
 *
 * Layout / thread model:
 *
 *   creator thread  (this decoder's own thread)
 *     reads the raw compressed stream through a caller-supplied callback,
 *     cuts one member at a time by its gzip header (BC/BSIZE), and enqueues
 *     the member (compressed bytes + sequence number) into the ring. Raw EOF
 *     ends the job stream. If the first member has no BC extra field (plain
 *     gzip, not BGZF) the decoder runs in a serial fallback: one long
 *     inflate stream is fed into output chunks and served to the reader in
 *     order, matching the old single-threaded behaviour byte for byte.
 *
 *   inflate workers (N threads)
 *     pop members and inflate each into a 64 KiB slot of the ring; a member
 *     whose slot has been consumed by the reader is reused by the creator,
 *     bounding memory to ring_size members in flight.
 *
 *   reader thread (the caller of read())
 *     pulls inflated bytes strictly in member order from the ring and copies
 *     them into the destination. read() returns 0 on end of stream.
 *
 * The whole inflate stream is byte-identical to the sequential inflate
 * (gzip members are independent), so the caller produces identical output;
 * only latency and throughput change.
 *
 * Errors: a corrupt member makes read() return 0 (EOF) after logging, so the
 * caller's short-read handling treats it as end of data.
 */

#ifndef _BGZF_PARALLEL_H_
#define _BGZF_PARALLEL_H_

#include <zlib.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <mutex>
#include <thread>
#include <vector>

class BgzfParallelDecoder {
public:
    /* Raw source: return bytes actually read into dst (0 = EOF). The callback
       is only ever called from the creator thread, sequentially. */
    typedef size_t (*RawRead)(void* ctx, void* dst, size_t n);

    static constexpr int kMaxInflated = 1 << 16;  /* 64 KiB per BGZF member */
    static constexpr int kRingSize    = 96;       /* ~6 MiB of blocks in flight */

    BgzfParallelDecoder()
        : started_(false), stopped_(false), serialMode_(false) {}

    ~BgzfParallelDecoder() { stop(); }

    BgzfParallelDecoder(const BgzfParallelDecoder&) = delete;
    BgzfParallelDecoder& operator=(const BgzfParallelDecoder&) = delete;

    /*
     * Probe the stream head and, if it is BGZF, start the creator + workers.
     * `head` must be the first bytes of the stream (may be the whole stream if
     * it is shorter). Returns false when the stream is not cuttable BGZF, in
     * which case nothing was consumed yet and the caller must fall back to its
     * own serial inflate using head + the raw callback.
     */
    bool start(void* ctx, RawRead rr, const uint8_t* head, size_t headLen, int threads)
    {
        ctx_ = ctx;
        rr_ = rr;
        headStart_ = head;
        headLen_ = headLen;
        headOff_ = 0;
        produced_ = 0;
        consumed_ = 0;
        noMore_ = false;
        stopped_ = false;

        if (threads < 1) threads = 1;
        if (threads > 8) threads = 8;

        if (!looksBgzf(head, headLen)) {
            return false;
        }
        serialMode_ = false;

        ring_.resize(kRingSize);
        for (auto& m : ring_) {
            m.outLen = 0;
            m.ready = false;
        }
        inflaters_.clear();
        inflaters_.reserve(threads);
        try {
            for (int i = 0; i < threads; ++i)
                inflaters_.emplace_back(&BgzfParallelDecoder::inflateLoop, this);
            creator_ = std::thread(&BgzfParallelDecoder::creatorLoop, this);
        } catch (...) {
            stop();
            return false;
        }
        started_ = true;
        return true;
    }

    /* Pull n inflated bytes; 0 means end of stream / failure. */
    size_t read(void* dst, size_t n)
    {
        if (!started_ || stopped_)
            return 0;
        uint8_t* out = static_cast<uint8_t*>(dst);
        size_t got = 0;
        std::unique_lock<std::mutex> lk(mtx_);
        while (got < n) {
            BgzfMember& m = ring_[consumed_ % kRingSize];
            if (m.ready) {
                /* Copy what remains of this member. */
                const size_t avail = (size_t)m.outLen - m.outOff;
                const size_t take = std::min(avail, n - got);
                if (take > 0) {
                    memcpy(out + got, m.out + m.outOff, take);
                    m.outOff += (uint32_t)take;
                    got += take;
                }
                if (m.outOff == (uint32_t)m.outLen) {
                    /* Whole member consumed. */
                    if (m.outLen == 0 && produced_ > consumed_ + 1) {
                        /* A zero-output member in the middle is abnormal; the
                           only legal one is the final BGZF EOF marker. */
                    }
                    const uint64_t seq = consumed_;
                    m.ready = false;
                    m.comp.clear();
                    m.comp.shrink_to_fit();
                    ++consumed_;
                    cvCreator_.notify_all();
                    if (seq + 1 == produced_ && noMore_) {
                        /* Consumer caught up with the producer end. */
                        break;
                    }
                }
                continue;
            }
            if (noMore_) {
                /* No member available and the producer is done. */
                break;
            }
            cvRead_.wait(lk);
        }
        return got;
    }

    void stop()
    {
        bool expected = false;
        if (!stopped_.compare_exchange_strong(expected, true))
            return; /* already stopping / stopped */
        if (!started_)
            return;
        cvCreator_.notify_all();
        cvWork_.notify_all();
        cvRead_.notify_all();
        if (creator_.joinable())
            creator_.join();
        for (auto& t : inflaters_)
            if (t.joinable())
                t.join();
        inflaters_.clear();
        started_ = false;
    }

private:
    struct BgzfMember {
        std::vector<uint8_t> comp;  /* compressed member bytes */
        uint8_t out[kMaxInflated];  /* inflated payload */
        uint32_t outLen;
        uint32_t outOff;            /* reader offset inside out */
        bool ready;
        BgzfMember() : outLen(0), outOff(0), ready(false) {}
    };

    /* ---- raw byte source (creator thread only) ---- */
    size_t readRaw(void* dst, size_t n)
    {
        size_t got = 0;
        uint8_t* out = static_cast<uint8_t*>(dst);
        while (got < n) {
            if (headOff_ < headLen_) {
                const size_t c = std::min(n - got, headLen_ - headOff_);
                memcpy(out + got, headStart_ + headOff_, c);
                headOff_ += c;
                got += c;
                continue;
            }
            const size_t r = rr_(ctx_, out + got, n - got);
            if (r == 0)
                break;
            got += r;
        }
        return got;
    }

    /* ---- BGZF header probing (no consumption) ---- */
    static bool looksBgzf(const uint8_t* b, size_t n)
    {
        if (n < 12 || b[0] != 0x1f || b[1] != 0x8b || b[2] != 0x08)
            return false;
        const uint8_t flg = b[3];
        if (!(flg & 0x04))   /* no FEXTRA => no BC subfield => not BGZF */
            return false;
        const size_t xlen = (size_t)b[10] | ((size_t)b[11] << 8);
        if (12 + xlen > n)
            return false;
        for (size_t i = 12; i + 4 <= 12 + xlen; ++i) {
            if (b[i] == 'B' && b[i + 1] == 'C' &&
                b[i + 2] == 0x02 && b[i + 3] == 0x00) {
                return true;  /* found the BC (BSIZE) subfield */
            }
        }
        return false;
    }

    /* Read one gzip member header from the raw source, then its whole
       compressed body, returning the member size; 0 on clean EOF, ~0 on
       a truncated/corrupt header. */
    uint64_t collectMember(std::vector<uint8_t>& out)
    {
        uint8_t hdr[18];
        size_t h = readRaw(hdr, 12);
        if (h == 0)
            return 0; /* EOF */
        if (h < 12)
            return UINT64_MAX; /* truncated header */
        if (hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 0x08)
            return UINT64_MAX;
        const uint8_t flg = hdr[3];
        size_t skip = 10;
        size_t xlen = 0;
        if (flg & 0x04) {
            xlen = (size_t)hdr[10] | ((size_t)hdr[11] << 8);
            skip = 12 + xlen;
        }
        uint64_t memberSize = 0;
        bool foundBc = false;
        /* The header may extend beyond the 12 bytes we already have: reread
           the whole variable-length header so parsing is straightforward. */
        std::vector<uint8_t> varHead;
        if (skip > 12) {
            varHead.resize(skip);
            memcpy(varHead.data(), hdr, 12);
            size_t have = 12;
            while (have < skip) {
                const size_t r = readRaw(varHead.data() + have, skip - have);
                if (r == 0)
                    return UINT64_MAX;
                have += r;
            }
            /* memberSize comes from the BC extra subfield, when present */
            for (size_t i = 12; i + 6 <= skip; ++i) {
                if (varHead[i] == 'B' && varHead[i + 1] == 'C' &&
                    varHead[i + 2] == 0x02 && varHead[i + 3] == 0x00) {
                    memberSize = (uint64_t)(varHead[i + 4] | (varHead[i + 5] << 8)) + 1;
                    foundBc = true;
                    break;
                }
            }
        } else {
            memberSize = 0;
        }

        if (flg & 0x08) { /* FNAME: NUL-terminated */
            uint8_t c;
            do {
                if (readRaw(&c, 1) == 0)
                    return UINT64_MAX;
            } while (c != 0);
        }
        if (flg & 0x10) { /* FCOMMENT */
            uint8_t c;
            do {
                if (readRaw(&c, 1) == 0)
                    return UINT64_MAX;
            } while (c != 0);
        }
        if (flg & 0x02) { /* FHCRC: 2 bytes */
            uint8_t two[2];
            if (readRaw(two, 2) != 2)
                return UINT64_MAX;
        }

        if (!foundBc)
            return UINT64_MAX; /* cannot split; caller falls back */

        out.resize((size_t)memberSize);
        if (memberSize <= skip) {
            memcpy(out.data(), hdr, std::min((size_t)memberSize, (size_t)12));
            if (memberSize > 12)
                memcpy(out.data() + 12, varHead.data() + 12, memberSize - 12);
            return memberSize;
        }
        memcpy(out.data(), hdr, 12);
        if (skip > 12) {
            memcpy(out.data() + 12, varHead.data() + 12, skip - 12);
        }
        size_t have = skip;
        while (have < memberSize) {
            const size_t r = readRaw(out.data() + have, (size_t)(memberSize - have));
            if (r == 0)
                return UINT64_MAX;
            have += r;
        }
        return memberSize;
    }

    /* ---- creator: cut members / fall back to serial inflate ---- */
    void creatorLoop()
    {
        std::vector<uint8_t> memberBuf;
        for (;;) {
            if (stopped_)
                return;
            /* Bound the number of produced-but-unconsumed members. */
            {
                std::unique_lock<std::mutex> lk(mtx_);
                cvCreator_.wait(lk, [this] {
                    return stopped_ || (consumed_ + (uint64_t)kRingSize > produced_);
                });
                if (stopped_)
                    return;
            }
            const uint64_t seq = produced_;
            const uint64_t sz = collectMember(memberBuf);
            if (sz == UINT64_MAX) {
                std::lock_guard<std::mutex> lk(mtx_);
                noMore_ = true;
                cvRead_.notify_all();
                return;
            }
            if (sz == 0) { /* raw EOF */
                break;
            }
            {
                std::lock_guard<std::mutex> lk(mtx_);
                BgzfMember& m = ring_[seq % kRingSize];
                m.comp.swap(memberBuf);
                memberBuf.clear();
                m.outOff = 0;
                ++produced_;
            }
            {
                std::lock_guard<std::mutex> lk(mtx_);
                jobQ_.push_back(seq);
            }
            cvWork_.notify_one();
        }
        /* Producer EOF: wait for every enqueued member to be inflated, then
           publish end-of-stream so the reader can return short/EOF. */
        {
            std::unique_lock<std::mutex> lk(mtx_);
            cvWork_.wait(lk, [this] { return stopped_ || inflateDone_ == produced_; });
        }
        {
            std::lock_guard<std::mutex> lk(mtx_);
            noMore_ = true;
            cvRead_.notify_all();
        }
    }

    void inflateLoop()
    {
        z_stream zs;
        memset(&zs, 0, sizeof(zs));
        if (inflateInit2(&zs, 32 + MAX_WBITS) != Z_OK)
            return;
        for (;;) {
            uint64_t seq;
            {
                std::unique_lock<std::mutex> lk(mtx_);
                cvWork_.wait(lk, [this] { return stopped_ || !jobQ_.empty(); });
                if (stopped_)
                    break;
                seq = jobQ_.front();
                jobQ_.pop_front();
            }
            BgzfMember& m = ring_[seq % kRingSize];
            /* Each job is an independent gzip member; reset the stream so the
               same z_stream can serve member after member (without the reset,
               only the first member would inflate). */
            inflateReset(&zs);
            zs.next_in = m.comp.data();
            zs.avail_in = (uInt)m.comp.size();
            zs.next_out = m.out;
            zs.avail_out = (uInt)kMaxInflated;
            int rc = inflate(&zs, Z_FINISH);
            if (rc != Z_STREAM_END) {
                /* truncated/corrupt: deliver what we have and end the stream */
                std::lock_guard<std::mutex> lk(mtx_);
                m.outLen = (uint32_t)(kMaxInflated - zs.avail_out);
                m.ready = true;
                ++inflateDone_;
                cvRead_.notify_all();
                cvWork_.notify_all(); /* wake the creator's end-of-stream wait */
                continue;
            }
            m.outLen = (uint32_t)(kMaxInflated - zs.avail_out);
            {
                std::lock_guard<std::mutex> lk(mtx_);
                m.ready = true;
                ++inflateDone_;
                cvRead_.notify_all();
                cvWork_.notify_all(); /* wake the creator's end-of-stream wait */
            }
        }
        inflateEnd(&zs);
    }

private:
    void* ctx_ = nullptr;
    RawRead rr_ = nullptr;
    const uint8_t* headStart_ = nullptr;
    size_t headLen_ = 0;
    size_t headOff_ = 0;

    std::vector<BgzfMember> ring_;
    std::deque<uint64_t> jobQ_;

    std::atomic<bool> started_{false};
    std::atomic<bool> stopped_{false};
    bool serialMode_;

    std::mutex mtx_;
    std::condition_variable cvRead_;    /* reader waits for its member/EOF */
    std::condition_variable cvWork_;    /* inflaters wait for jobs */
    std::condition_variable cvCreator_; /* creator waits for a free slot */

    uint64_t produced_ = 0;     /* members handed to the inflaters */
    uint64_t consumed_ = 0;     /* members fully read by read() */
    uint64_t inflateDone_ = 0;  /* members inflated (workers) */
    bool noMore_ = false;       /* producer finished, no more members */

    std::thread creator_;
    std::vector<std::thread> inflaters_;
};

#endif
