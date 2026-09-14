/*
 * profile_stats.h - lightweight pipeline phase profiler
 *
 * Header-only timer used to see where wall time goes in the
 * "read -> parallel compress -> write" pipeline. Every category is a pair of
 * atomic counters, so instrumentation is thread safe and cheap: only a couple
 * of clock reads per block, never per line.
 *
 * Enabled at run time with PBGZ_PROF=1; when disabled the scope objects do not
 * read the clock at all, so a profiled build behaves like a normal one.
 *
 * Categories:
 *   READ_*   reader thread      (raw read/split, codec pre-selection, QUAL
 *                                prior accumulation and training)
 *   CODER_*  worker threads     (queue waits, actuator create/pre-analysis,
 *                                compression)
 *   WRITE_*  writer thread      (queue wait, reorder sort, write)
 *   FIELD_BASE + fieldIdx       per-column cost inside the SAM actuator
 */

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <vector>

namespace pbgzprof {

enum Cat {
    READ_TOTAL = 0,
    READ_BLOCK,
    READ_DECISION,
    READ_PRETRAIN,
    READ_PRETRAIN_TRAIN,
    READ_QUEUE_WAIT,
    READ_MD5,
    READ_INFLATE,
    READ_PARSE,

    /* SEQ / QUAL internals: preparation vs coder construction vs encoding */
    SEQ_PREP = 20,
    SEQ_CTOR,
    SEQ_CODEC,
    QUAL_PREP = 24,
    QUAL_CTOR,
    QUAL_CODEC,

    CODER_TOTAL,
    CODER_BARRIER,
    CODER_WAIT_INPUT,
    CODER_WAIT_OUTPUT,
    CODER_CREATE,
    CODER_COMPRESS,
    CODER_PARSE,      /* worker-side BAM record stream -> BamColumns parsing */

    WRITE_TOTAL,
    WRITE_WAIT,
    WRITE_SORT,
    WRITE_BLOCK,

    FIELD_BASE = 64,
    FIELD_HDR = 90,
    FIELD_META = 91,
    FIELD_INDEX = 92,
    FIELD_QNAME_TRIAL = 93,

    /* codec pre-selection: one slot per field trial, READ_TRIAL_BASE + fieldIdx */
    READ_TRIAL_BASE = 100,
    READ_TRIAL_QUAL_EXTRACT = 116,
    READ_TRIAL_QUAL_SELECT = 117,

    CAT_MAX = 128
};

inline const int kFieldNameCount = 12;

inline const char* const* trialNames()
{
    static const char* kNames[] = {
        "TRIAL_F0_QNAME", "TRIAL_F1_FLAG", "TRIAL_F2_RNAME", "TRIAL_F3_POS", "TRIAL_F4_MAPQ",
        "TRIAL_F5_CIGAR", "TRIAL_F6_RNEXT", "TRIAL_F7_PNEXT", "TRIAL_F8_TLEN", "TRIAL_F9_SEQ",
        "TRIAL_F10_QUAL", "TRIAL_F11_OPT"
    };
    return kNames;
}

inline const char* const* fieldNames()
{
    static const char* kNames[] = {
        "F0_QNAME", "F1_FLAG", "F2_RNAME", "F3_POS", "F4_MAPQ", "F5_CIGAR",
        "F6_RNEXT", "F7_PNEXT", "F8_TLEN", "F9_SEQ", "F10_QUAL", "F11_OPT"
    };
    return kNames;
}

struct Slot {
    std::atomic<uint64_t> us;
    std::atomic<uint64_t> cnt;
};

inline Slot* slots()
{
    static Slot s[CAT_MAX];
    return s;
}

inline bool enabled()
{
    static const bool on = []() -> bool {
        const char* env = getenv("PBGZ_PROF");
        return (env != nullptr && (strcmp(env, "1") == 0 || strcmp(env, "on") == 0));
    }();
    return on;
}

inline std::chrono::steady_clock::time_point nowOrZero()
{
    return enabled() ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
}

inline void add(int cat, uint64_t us)
{
    if (cat < 0 || cat >= CAT_MAX) {
        return;
    }
    slots()[cat].us.fetch_add(us, std::memory_order_relaxed);
    slots()[cat].cnt.fetch_add(1, std::memory_order_relaxed);
}

/* Adds the elapsed time since t0 to cat; used where a scope object would
   hide declarations from the code that follows. A timer that was never started reads as the clock's
   epoch, and "since the epoch" is not a duration anyone wants in a profile - it is reported as not
   measured instead, so a path that forgets to start its timer under-reports rather than swamping
   the table with the machine's uptime. */
inline void addSince(int cat, const std::chrono::steady_clock::time_point& t0)
{
    if (!enabled() || t0 == std::chrono::steady_clock::time_point()) {
        return;
    }
    const auto d = std::chrono::steady_clock::now() - t0;
    add(cat, (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(d).count());
}


class Scope {
public:
    explicit Scope(int cat) : cat_(cat), on_(enabled())
    {
        if (on_) {
            t0_ = std::chrono::steady_clock::now();
        }
    }

    ~Scope()
    {
        if (on_) {
            const auto d = std::chrono::steady_clock::now() - t0_;
            add(cat_, (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(d).count());
        }
    }

private:
    int cat_;
    bool on_;
    std::chrono::steady_clock::time_point t0_;
};

inline const char* catName(int cat)
{
    switch (cat) {
        case READ_TOTAL:          return "READ_total";
        case READ_BLOCK:          return "READ_readBlock";
        case READ_DECISION:       return "READ_fileDecision(codec trial)";
        case READ_PRETRAIN:       return "READ_pretrainAccum";
        case READ_PRETRAIN_TRAIN: return "READ_pretrainTrain(incl above)";
        case READ_QUEUE_WAIT:     return "READ_waitFreeBlock";
        case READ_MD5:            return "READ_md5(reader side)";
        case READ_INFLATE:        return "READ_bgzf(inflate+io)";
        case READ_PARSE:          return "READ_parseRec2Cols";
        case SEQ_PREP:            return "SEQ_prep(buf+collect)";
        case SEQ_CTOR:            return "SEQ_ctor";
        case SEQ_CODEC:           return "SEQ_encode+flush";
        case QUAL_PREP:           return "QUAL_prep";
        case QUAL_CTOR:           return "QUAL_ctor";
        case QUAL_CODEC:          return "QUAL_encode+flush";
        case READ_TRIAL_QUAL_EXTRACT: return "TRIAL_qual_extract";
        case READ_TRIAL_QUAL_SELECT:  return "TRIAL_qual_select";
        case CODER_TOTAL:         return "CODER_total";
        case CODER_BARRIER:       return "CODER_startBarrier";
        case CODER_WAIT_INPUT:    return "CODER_waitInput";
        case CODER_WAIT_OUTPUT:   return "CODER_waitOutput";
        case CODER_CREATE:        return "CODER_createActuator(preAnalysis)";
        case CODER_COMPRESS:      return "CODER_compress";
        case CODER_PARSE:         return "CODER_parseRec2Cols(worker)";
        case WRITE_TOTAL:         return "WRITE_total";
        case WRITE_WAIT:          return "WRITE_waitOutput";
        case WRITE_SORT:          return "WRITE_sort";
        case WRITE_BLOCK:         return "WRITE_writeBlock";
        case FIELD_HDR:           return "F_samHeader";
        case FIELD_META:          return "F_meta(md5+json)";
        case FIELD_INDEX:         return "F_buildIndex";
        case FIELD_QNAME_TRIAL:     return "F0_QNAME_trial";
        default:
            if (cat >= FIELD_BASE && cat < FIELD_BASE + kFieldNameCount) {
                return fieldNames()[cat - FIELD_BASE];
            }
            if (cat >= READ_TRIAL_BASE && cat < READ_TRIAL_BASE + kFieldNameCount) {
                return trialNames()[cat - READ_TRIAL_BASE];
            }
            return "unknown";
    }
}

/*
 * Print every non-empty category, hottest first.
 * wallSec is the pipeline wall time used as the 100 percent reference; because
 * the stages run concurrently, per-thread numbers add up to more than it.
 */
inline void dump(FILE* f, double wallSec)
{
    if (!enabled()) {
        return;
    }

    struct Row {
        const char* name;
        uint64_t us;
        uint64_t cnt;
    };

    std::vector<Row> rows;
    for (int i = 0; i < CAT_MAX; ++i) {
        const uint64_t us = slots()[i].us.load(std::memory_order_relaxed);
        if (us == 0 && slots()[i].cnt.load(std::memory_order_relaxed) == 0) {
            continue;
        }
        Row r;
        r.name = catName(i);
        r.us = us;
        r.cnt = slots()[i].cnt.load(std::memory_order_relaxed);
        rows.push_back(r);
    }

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.us > b.us; });

    fprintf(f, "\n===== PBGZ pipeline profile (PBGZ_PROF) =====\n");
    fprintf(f, "wall = %.3f s\n", wallSec);
    fprintf(f, "%-34s %12s %10s %8s %12s\n", "category", "total_ms", "count", "%wall", "avg_us");
    fprintf(f, "--------------------------------------------------------------------------------\n");
    for (size_t i = 0; i < rows.size(); ++i) {
        const double ms = (double)rows[i].us / 1000.0;
        const double pct = (wallSec > 0.0) ? (ms / 10.0) / wallSec : 0.0;
        const double avg = (rows[i].cnt > 0) ? (double)rows[i].us / (double)rows[i].cnt : 0.0;
        fprintf(f, "%-34s %12.1f %10llu %7.1f%% %12.1f\n", rows[i].name, ms,
                (unsigned long long)rows[i].cnt, pct, avg);
    }
    fprintf(f, "===========================================\n\n");
    fflush(f);
}

inline void reset()
{
    for (int i = 0; i < CAT_MAX; ++i) {
        slots()[i].us.store(0, std::memory_order_relaxed);
        slots()[i].cnt.store(0, std::memory_order_relaxed);
    }
}

}  // namespace pbgzprof

#define PBGZ_PROF_SCOPE(cat) pbgzprof::Scope profScope##__LINE__(cat)
