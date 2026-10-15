/*
 * qual_selector.cpp - encoder evaluation for the quality-value column. Split from
 * codec_selector.cpp so coder_qual.h and coder_fc.h never share a translation unit.
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

#include "qual_selector.h"

#include <string.h>
#include <algorithm>
#include <chrono>
#include <atomic>
#include <cmath>
#include <functional>
#include <thread>

#include "profile_stats.h"

#include "coder/coder_io.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_qcm.h"
#include "field_coder_config.h"
#include "log/logger.h"

/* Mean record length at/above which fcv2 is not trialled; the fcv2 model count is fixed. */
const uint32_t QUAL_LONG_READ_LEN = 500;

/* Whether a candidate list a caller passed in holds this coder. */
static bool contains(const std::vector<CoderType>& list, CoderType type)
{
    return std::find(list.begin(), list.end(), type) != list.end();
}

namespace {

/* Below this sample size, skip evaluation and keep the default encoder. */
const uint32_t MIN_QUAL_SAMPLE = 64u << 10;

/* Slack for the trial-compression output buffer. Quality values are almost
 * never incompressible, so a 2x margin is more than safe. */
inline size_t trialCapacity(size_t srcLen)
{
    return (srcLen << 1) + (1u << 16);
}

/* Lead the winner must hold to settle at the first rung (1 MB). */
const double SETTLE_MARGIN = 0.03;

/* Smaller lead accepted at rungs >= 2 MB, where the ranking has held across a doubling. */
const double LATER_SETTLE_MARGIN = 0.003;

/* The margin used from the 2 MB rung up. */
inline double resolveLaterSettleMargin()
{
    return LATER_SETTLE_MARGIN;
}

/* Return the largest n such that the cumulative QUAL bytes of the first n
 * records do not exceed budget (only complete records are kept). */
size_t recordsForBudget(const std::vector<QualSampleRecord>& records, size_t budget)
{
    size_t used = 0;
    size_t n = 0;
    for (; n < records.size(); n++) {
        if (used + records[n].qual.size() > budget) {
            break;
        }
        used += records[n].qual.size();
    }
    return n;
}

/* Internal block sizes of coder_bwt_cm, kept consistent with BWT_LEVEL_SIZE in
 * codec_selector.cpp. */
const uint32_t kBwtLevelSize[10] = {
    0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
    1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
};

/* Pick the smallest tier whose internal block can hold this round's sample. */
int bwtLevelFor(uint32_t sampleLen)
{
    for (int lv = 1; lv <= 9; ++lv) {
        if (kBwtLevelSize[lv] >= sampleLen) {
            return lv;
        }
    }
    return 9;
}

/* Trial-compress with coder_qual; each record is fed with its SEQ as context. */
bool trialQual(const std::vector<QualSampleRecord>& records, size_t recordCount,
               const std::vector<uint32_t>& freqByByte,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<std::pair<uint16_t, uint16_t>> freqTable;
    for (uint32_t b = 0; b < 256; b++) {
        if (b < freqByByte.size() && freqByByte[b] > 0) {
            freqTable.push_back(std::make_pair((uint16_t)(b - '!'), (uint16_t)1));
        }
    }
    if (freqTable.empty()) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_qual coder(&io, true, freqTable);
        for (size_t i = 0; i < recordCount; i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_qual_gen2((uint8_t*)r.seq.data(),
                                   (uint8_t*)r.qual.data(),
                                   (uint32_t)r.qual.size());
        }
        coder.encode_flush();
        outLen = (uint32_t)io.data_len;
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/* Trial-compress with fcv2; each record is fed with its length and strand. */
bool trialFcv2(const std::vector<QualSampleRecord>& records, size_t recordCount,
               const std::vector<uint32_t>& freqByByte, const Fcv2Cfg& cfg,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_fcv2 coder(&io, freqByByte, cfg);
        for (size_t i = 0; i < recordCount; i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_record((const uint8_t*)r.qual.data(),
                                (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        outLen = (uint32_t)coder.encode_flush();
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/* Trial-compress with coder_bwt_cm, one encode_line per record, at the round's block size. */
bool trialBwtCm(const std::vector<QualSampleRecord>& records, size_t recordCount,
                int bwtLevel, uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        io.set_level(bwtLevel);
        coder_bwt_cm coder(&io);
        for (size_t i = 0; i < recordCount; i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_line((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size());
        }
        coder.encode_flush();
        outLen = (io.data_len > 0) ? (uint32_t)io.data_len : 0;
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/* Trial-compress with coder_qcm, an ordinary byte-stream coder fed record by record. */
bool trialQcm(const std::vector<QualSampleRecord>& records, size_t recordCount,
              uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_qcm coder(&io);
        for (size_t i = 0; i < recordCount; i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_line((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size());
        }
        coder.encode_flush();
        outLen = (io.data_len > 0) ? (uint32_t)io.data_len : 0;
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/* Translate the coder-layer tier into the parameters carried by PreprocessInfo
 * (fields correspond one-to-one). */
QualFcv2Params toQualParams(const Fcv2Cfg& cfg)
{
    QualFcv2Params p;
    p.cycleMax = cfg.cycleMax;
    p.cycleBucket = cfg.cycleBucket;
    p.prevShift = cfg.prevShift;
    p.useDedup = cfg.useDedup;
    return p;
}

/* Map compression level to the SAM block read-count tier: 1-5 -> ~10000 reads/block,
   6-7 -> ~25000, 8-9 -> ~100000. */
static int levelVolumeTier(uint8_t compressLevel)
{
    if (compressLevel >= 8) return 3;
    if (compressLevel >= 6) return 2;
    return 1;
}

/* fcv2 candidate tiers for the given compression level, filtered by block-volume tier. */
std::vector<Fcv2Cfg> candidateFcv2Cfgs(uint8_t compressLevel)
{
    std::vector<Fcv2Cfg> cfgs;
    const int volumeTier = levelVolumeTier(compressLevel);
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        if (volumeTier < FCV2_PRESETS[i].minVolumeTier) {
            continue;
        }
        cfgs.push_back(FCV2_PRESETS[i].cfg);
    }
    return cfgs;
}

/* One round that compresses every candidate under a byte budget. Returns each
 * candidate's compressed size and whether it succeeded. */
struct QualRoundResult {
    bool     qualOk = false;
    uint32_t qualLen = 0;
    bool     fcv2Ok = false;
    uint32_t fcv2Len = 0;
    Fcv2Cfg  fcv2Cfg;          /* the winning one among the several fcv2 tiers */
    bool     cmOk = false;
    uint32_t cmLen = 0;
    bool     qcmOk = false;
    uint32_t qcmLen = 0;
    uint32_t bestLen = 0;
    CoderType bestCoder = CoderType::QUAL;
    bool     anyOk = false;
    uint32_t fcv2Us = 0;
    uint32_t qualUs = 0;
    uint32_t cmUs = 0;
    uint32_t qcmUs = 0;
};

/* Run one rung's trial-compressions across the idle cores. */
void runInParallel(const std::vector<std::function<void()>>& tasks)
{
    if (tasks.empty()) {
        return;
    }
    const unsigned hw = std::thread::hardware_concurrency();
    size_t workers = (hw == 0) ? 1u : (size_t)hw;
    if (workers > tasks.size()) {
        workers = tasks.size();
    }
    if (workers <= 1) {
        for (size_t i = 0; i < tasks.size(); ++i) {
            tasks[i]();
        }
        return;
    }

    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (size_t i = 1; i < workers; ++i) {
        pool.emplace_back([&tasks, &next]() {
            while (true) {
                const size_t k = next.fetch_add(1, std::memory_order_relaxed);
                if (k >= tasks.size()) {
                    break;
                }
                tasks[k]();
            }
        });
    }
    while (true) {
        const size_t k = next.fetch_add(1, std::memory_order_relaxed);
        if (k >= tasks.size()) {
            break;
        }
        tasks[k]();
    }
    for (size_t i = 0; i < pool.size(); ++i) {
        pool[i].join();
    }
}

/* One rung's raw measurements, written by the parallel tasks. */
struct RoundSlots {
    std::vector<bool>     fcv2Ok;
    std::vector<uint32_t> fcv2Len;
    std::vector<uint32_t> fcv2Us;
    bool     qualOk = false;
    uint32_t qualLen = 0;
    uint32_t qualUs = 0;
    bool     cmOk = false;
    uint32_t cmLen = 0;
    uint32_t cmUs = 0;
    bool     qcmOk = false;
    uint32_t qcmLen = 0;
    uint32_t qcmUs = 0;
};

/* Reduce one rung to a verdict - the same comparison the serial code did. */
QualRoundResult assembleRound(const RoundSlots& s, const std::vector<Fcv2Cfg>& cfgs)
{
    QualRoundResult r;
    if (s.qualOk) {
        r.qualOk = true;
        r.qualLen = s.qualLen;
        r.qualUs = s.qualUs;
        r.anyOk = true;
        r.bestLen = r.qualLen;
        r.bestCoder = CoderType::QUAL;
    }

    bool picked = false;
    for (size_t i = 0; i < cfgs.size(); i++) {
        if (!s.fcv2Ok[i]) {
            continue;
        }
        if (!picked || s.fcv2Len[i] < r.fcv2Len) {
            r.fcv2Len = s.fcv2Len[i];
            r.fcv2Us = s.fcv2Us[i];
            r.fcv2Cfg = cfgs[i];
            picked = true;
        }
    }
    if (picked) {
        r.fcv2Ok = true;
        r.anyOk = true;
        if (r.fcv2Len < r.bestLen || !r.qualOk) {
            r.bestLen = r.fcv2Len;
            r.bestCoder = CoderType::FCV2;
        }
    }

    if (s.cmOk) {
        r.cmOk = true;
        r.cmUs = s.cmUs;
        r.cmLen = s.cmLen;
        r.anyOk = true;
        if (r.cmLen < r.bestLen || (!r.qualOk && !r.fcv2Ok)) {
            r.bestLen = r.cmLen;
            r.bestCoder = CoderType::BWT_CM;
        }
    }

    /* qcm is compared last with a strict `<`, so a tie keeps bwt_cm. */
    if (s.qcmOk) {
        r.qcmOk = true;
        r.qcmUs = s.qcmUs;
        r.qcmLen = s.qcmLen;
        r.anyOk = true;
        if (r.qcmLen < r.bestLen || (!r.qualOk && !r.fcv2Ok && !r.cmOk)) {
            r.bestLen = r.qcmLen;
            r.bestCoder = CoderType::QCM;
        }
    }
    return r;
}

} /* namespace */

FieldCodecSelection QualSelector::select(const std::vector<QualSampleRecord>& records,
                                         const std::vector<uint32_t>& freqByByte,
                                         uint8_t compressLevel)
{
    return select(records, freqByByte, compressLevel,
                  samFieldCandidates(SAM_QUAL, PBGZ_MODE_ARCHIVE));
}

FieldCodecSelection QualSelector::select(const std::vector<QualSampleRecord>& records,
                                         const std::vector<uint32_t>& freqByByte,
                                         uint8_t compressLevel,
                                         const std::vector<CoderType>& candidates)
{
    FieldCodecSelection sel;

    uint32_t sampleLen = 0;
    for (size_t i = 0; i < records.size(); i++) {
        sampleLen += (uint32_t)records[i].qual.size();
    }
    sel.sampleLen = sampleLen;
    sel.decidedLen = sampleLen;
    sel.rounds = 1;

    if (records.empty() || sampleLen < MIN_QUAL_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    /* The sample's mean record length prunes the fcv2 candidate set. */
    const uint64_t meanLen = (uint64_t)sampleLen / (uint64_t)records.size();

    /* Drop fcv2 when its Huffman tree cannot be built (alphabet > FCV2_TREE_CAP) or on long
       reads, where coder_qual / coder_qcm win. */
    uint32_t qualityAlphabet = 0;
    for (size_t b = 0; b < freqByByte.size(); ++b) {
        if (freqByByte[b] > 0) {
            ++qualityAlphabet;
        }
    }
    std::vector<CoderType> activeCandidates = candidates;
    if ((qualityAlphabet > (uint32_t)FCV2_TREE_CAP || meanLen >= (uint64_t)QUAL_LONG_READ_LEN) &&
        contains(activeCandidates, CoderType::FCV2)) {
        activeCandidates.erase(
            std::remove(activeCandidates.begin(), activeCandidates.end(), CoderType::FCV2),
            activeCandidates.end());
    }

    /* Only the coders the caller's candidate row lists are trialled. */
    const bool tryQual  = contains(activeCandidates, CoderType::QUAL);
    const bool tryFcv2  = contains(activeCandidates, CoderType::FCV2);
    const bool tryBwtCm = contains(activeCandidates, CoderType::BWT_CM);
    const bool tryQcm   = contains(activeCandidates, CoderType::QCM);

    /* fcv2 tiers are built only when fcv2 is in play. */
    const std::vector<Fcv2Cfg> cfgs =
        tryFcv2 ? candidateFcv2Cfgs(compressLevel) : std::vector<Fcv2Cfg>();

    /* Multi-round convergence: double the sample each round, finalizing once the leader opens
       enough of a gap, and never below MIN_SETTLE_PROBE where the coders have not converged. */
    const uint32_t MIN_SETTLE_PROBE = 1u << 20;   /* 1 MB */

    /* Enumerate the probe ladder up front, then trial one rung at a time in parallel. */
    std::vector<uint32_t> probes;
    {
        /* Start the ladder at the oldest rung any settle rule reads; earlier rungs are unused. */
        const uint32_t firstUsefulProbe = MIN_SETTLE_PROBE / 2;
        uint32_t p = (MIN_QUAL_SAMPLE < sampleLen) ? MIN_QUAL_SAMPLE : sampleLen;
        while (p < firstUsefulProbe && p < sampleLen) {
            p = (p > sampleLen / 2) ? sampleLen : (p * 2);
        }
        while (true) {
            probes.push_back(p);
            if (p >= sampleLen) {
                break;
            }
            p = (p > sampleLen / 2) ? sampleLen : (p * 2);
        }
    }

    /* Slots for the rung being trialled, fresh per rung. */

    QualRoundResult final;
    bool finalSet = false;
    uint32_t settleProbe = sampleLen;

    /* Ratio-plateau early exit: stop once the top two ratios stop moving across doublings. */
    const uint32_t MIN_PLATEAU_PROBE = 2u << 20;   /* 2 MB */
    const double plPlateauBand = 0.002;            /* 0.2% absolute ratio movement */
    const bool plateauEnabled = true;
    const double laterSettleMargin = resolveLaterSettleMargin();
    QualRoundResult prevRound;
    QualRoundResult prevPrevRound;
    uint32_t prevProbe = 0;
    bool havePrev = false;
    bool havePrevPrev = false;

    /* Walk the ladder lazily: a rung is trialled only if the one before it did not settle. */
    for (size_t i = 0; i < probes.size(); i++) {
        const uint32_t probe = probes[i];
        const size_t count = recordsForBudget(records, probe);
        RoundSlots slot;
        slot.fcv2Ok.assign(cfgs.size(), false);
        slot.fcv2Len.assign(cfgs.size(), 0);
        slot.fcv2Us.assign(cfgs.size(), 0);
        std::vector<std::function<void()>> tasks;
        tasks.reserve(cfgs.size() + 3);
        if (tryQual) {
            tasks.push_back([&records, &freqByByte, &slot, count]() {
                slot.qualOk = trialQual(records, count, freqByByte, slot.qualLen, slot.qualUs);
            });
        }
        for (size_t c = 0; c < cfgs.size(); c++) {
            tasks.push_back([&records, &freqByByte, &slot, &cfgs, c, count]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_QUAL_FCV2_BASE + (int)c);
                uint32_t len = 0;
                uint32_t us = 0;
                if (trialFcv2(records, count, freqByByte, cfgs[c], len, us)) {
                    slot.fcv2Ok[c] = true;
                    slot.fcv2Len[c] = len;
                    slot.fcv2Us[c] = us;
                }
            });
        }
        if (tryBwtCm) {
            const int bwtLevel = bwtLevelFor(probe);
            tasks.push_back([&records, &slot, count, bwtLevel]() {
                slot.cmOk = trialBwtCm(records, count, bwtLevel, slot.cmLen, slot.cmUs);
            });
        }
        if (tryQcm) {
            tasks.push_back([&records, &slot, count]() {
                slot.qcmOk = trialQcm(records, count, slot.qcmLen, slot.qcmUs);
            });
        }
        runInParallel(tasks);

        QualRoundResult r = assembleRound(slot, cfgs);
        final = r;
        finalSet = true;
        settleProbe = probe;
        sel.rounds++;

        /* Top-two (by compressed length) of the candidates that actually encoded. */
        struct TopTwo {
            int  n;
            uint32_t t[2];
            double rat[2];
        };
        auto topTwo = [](const QualRoundResult& rr, uint32_t pb) {
            struct Cand { uint32_t len; uint32_t type; double rat; };
            Cand c[4];
            int m = 0;
            const auto add = [&](bool ok, uint32_t len, uint32_t type) {
                if (ok) {
                    c[m++] = { len, type, (double)len / (double)pb };
                }
            };
            add(rr.qualOk, rr.qualLen, 0);
            add(rr.fcv2Ok, rr.fcv2Len, 1);
            add(rr.cmOk, rr.cmLen, 2);
            add(rr.qcmOk, rr.qcmLen, 3);
            for (int a = 1; a < m; ++a) {
                Cand t = c[a];
                int b = a - 1;
                while (b >= 0 && c[b].len > t.len) {
                    c[b + 1] = c[b];
                    --b;
                }
                c[b + 1] = t;
            }
            TopTwo tt;
            tt.n = m;
            tt.t[0] = (m > 0) ? c[0].type : 0;
            tt.t[1] = (m > 1) ? c[1].type : 0;
            tt.rat[0] = (m > 0) ? c[0].rat : 0.0;
            tt.rat[1] = (m > 1) ? c[1].rat : 0.0;
            return tt;
        };

        if (plateauEnabled && havePrev && havePrevPrev &&
            probe >= MIN_PLATEAU_PROBE &&
            probe == prevProbe * 2) {
            const TopTwo cur = topTwo(r, probe);
            const TopTwo prev = topTwo(prevRound, prevProbe);
            if (cur.n >= 2 && prev.n >= 2 &&
                r.bestCoder == prevRound.bestCoder &&
                prevRound.bestCoder == prevPrevRound.bestCoder &&
                cur.t[0] == prev.t[0] && cur.t[1] == prev.t[1] &&
                fabs(cur.rat[0] - prev.rat[0]) <= plPlateauBand &&
                fabs(cur.rat[1] - prev.rat[1]) <= plPlateauBand) {
                break;
            }
        }

        prevPrevRound = prevRound;
        prevRound = r;
        prevProbe = probe;
        havePrevPrev = havePrev;
        havePrev = true;

        if (!r.anyOk) {
            break;
        }
        if (probe >= sampleLen) {
            break;
        }
        if (probe < MIN_SETTLE_PROBE) {
            continue;   /* too early to trust the ranking, just keep doubling */
        }
        /* Only one candidate encoded; more data cannot change the comparison. */
        {
            /* The runner-up is the smallest candidate that is not the leader. */
            uint32_t runnerUp = UINT32_MAX;
            if (r.qualOk && r.bestCoder != CoderType::QUAL && r.qualLen < runnerUp) {
                runnerUp = r.qualLen;
            }
            if (r.fcv2Ok && r.bestCoder != CoderType::FCV2 && r.fcv2Len < runnerUp) {
                runnerUp = r.fcv2Len;
            }
            if (r.cmOk && r.bestCoder != CoderType::BWT_CM && r.cmLen < runnerUp) {
                runnerUp = r.cmLen;
            }
            if (r.qcmOk && r.bestCoder != CoderType::QCM && r.qcmLen < runnerUp) {
                runnerUp = r.qcmLen;
            }
            if (runnerUp == UINT32_MAX) {
                break;
            }
            /* The leader is far enough ahead that more data will not flip the ranking. */
            const double settleMargin =
                (probe <= (1u << 20)) ? SETTLE_MARGIN : laterSettleMargin;
            if ((double)(runnerUp - r.bestLen) / (double)runnerUp >= settleMargin) {
                break;
            }
        }
    }

    sel.decidedLen = finalSet ? settleProbe : sampleLen;
    sel.trialCount = 0;
    /* Only the coders put in play are reported. */
    if (tryQual) {
        sel.addTrial(CoderType::QUAL, final.qualOk ? final.qualLen : 0, final.qualUs);
    }
    if (tryFcv2) {
        sel.addTrial(CoderType::FCV2, final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us);
    }
    if (tryBwtCm) {
        sel.addTrial(CoderType::BWT_CM, final.cmOk ? final.cmLen : 0, final.cmUs);
    }
    if (tryQcm) {
        sel.addTrial(CoderType::QCM, final.qcmOk ? final.qcmLen : 0, final.qcmUs);
    }

    if (!final.anyOk) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    sel.selectedCoder = final.bestCoder;
    sel.bestCompLen = final.bestLen;
    if (final.bestCoder == CoderType::FCV2) {
        /* Hand the winning fcv2 tier to the caller. */
        sel.fcv2Params = toQualParams(final.fcv2Cfg);
    }
    sel.status = FieldStatus::SELECTED;

    /* Name the winning preset, not just the coder. */
    const char* preset = (final.bestCoder == CoderType::FCV2) ? fcv2PresetName(final.fcv2Cfg) : nullptr;
    LOG_DEBUG("Qual codec trial: coder_qual=%u (%u us), fcv2=%u (%u us), bwt_cm=%u (%u us), "
              "qcm=%u (%u us), mean_read_len=%llu, picked=%s%s%s",
              final.qualOk ? final.qualLen : 0, final.qualUs,
              final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us,
              final.cmOk ? final.cmLen : 0, final.cmUs,
              final.qcmOk ? final.qcmLen : 0, final.qcmUs,
              (unsigned long long)meanLen,
              coderTypeToMagic(sel.selectedCoder),
              (preset != nullptr) ? "/" : "",
              (preset != nullptr) ? preset : "");
    return sel;
}