/*
 * codec_selector.cpp - Codec pre-selection for the file preprocessing module
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

#include "codec_selector.h"

#include <chrono>

#include <climits>
#include <string>
#include <vector>

#include "io_block.h"
#include "block_wrapper.h"
#include "log/logger.h"
#include "safe_line_reader.h"
#include "reference.h"
#include "sam_info.h"
#include "sam_qname_column.h"
#include "sam_seq_payload.h"
#include "seq_stream_util.h"   /* SeqRleSplit, used by the SEQ trials below */
#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_arith.h"
#include "coder/coder_io.h"
#include "coder_factory.h"
#include "field_coder_config.h"

#include <atomic>
#include <functional>
#include <thread>
#include "profile_stats.h"

namespace {

/* Raw line bytes scanned from the start of the first block. */
const uint32_t SAMPLE_TARGET = 4u << 20;   /* 4 MB */

/* Each candidate is measured on PROBE_START of the SAMPLE_TARGET sample (see trialProbeStart). */
static uint32_t trialSampleTarget()
{
    return SAMPLE_TARGET;
}

/* Minimum per-field sample before a codec comparison is trusted; smaller fields are SKIPPED. */
const uint32_t MIN_SELECT_SAMPLE = 64u << 10;   /* 64 KB */

/* SEQ bases taken for the reference trial. */
const uint32_t kSeqRefSampleBytes = 4u << 20;   /* 4 MB of SEQ */

/* Cap on the QUAL sample for the codec trial: bounds the ladder's trial cost while
   still letting the l8 fcv2 preset compete. Separate from QUAL_PRIOR_TRAIN_MAX. */
const uint32_t QUAL_TRIAL_MAX = 32u << 20;   /* 32 MB of QUAL */

/* Lines of the QNAME column encoded with each candidate in the layout trial. */
const uint32_t QNAME_TRIAL_LINES = 20000;

/* coder_bwt_cm internal block sizes per level, mirroring coder_bwt_cm. */
const uint32_t BWT_LEVEL_SIZE[10] = {
    0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
    1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
};

/* Minimum POS delta-varint bytes before trusting the bwt_cm vs coder_arith comparison. */
const uint32_t MIN_POS_DELTA_SAMPLE = 8u << 10;   /* 8 KB of varint bytes */

/* Upper bound on POS delta trial lines: the trial must run near block granularity, but
   its cost must stay bounded. */
const uint32_t POS_TRIAL_MAX_LINES = 32u << 10;   /* 32768 lines */

/* Minimum estimated whole-file POS varint bytes before writing the ~512-byte prior. */
const uint64_t POS_PRIOR_MIN_ESTIMATED = 100ull << 10;   /* 100 KB varint */

/* Pack the varint histogram into the 512-byte prior layout (256 LE uint16 weights,
   total ~2^14). */
std::vector<uint8_t> packPosPriorBlob(const std::vector<uint64_t>& counts)
{
    std::vector<uint8_t> blob(coder_arith::kPriorBytes, 0);
    uint64_t total = 0;
    for (int i = 0; i < 256; ++i)
        total += counts[i];
    if (total == 0)
        return {}; /* no data: no prior */
    for (int i = 0; i < 256; ++i) {
        uint64_t w = (counts[i] * (1u << 14)) / total;
        if (w == 0 && counts[i] != 0)
            w = 1;
        if (w > 0xFFFFu)
            w = 0xFFFFu;
        blob[2 * i] = (uint8_t)(w & 0xff);
        blob[2 * i + 1] = (uint8_t)(w >> 8);
    }
    return blob;
}

/* Trial-compress a stream with one coder type and report the compressed size. The
   output buffer is generously sized (2x input + margin). */
template <typename CoderT>
bool trialEncode(const uint8_t* data, uint32_t len, int bwtLevel, uint32_t& outLen, uint32_t& usec)
{
    const auto t0 = std::chrono::steady_clock::now();
    uint32_t cap = (len << 1) + 65536;
    std::vector<uint8_t> outBuf(cap, 0);
    coder_io io(outBuf.data(), (int32_t)cap);
    if (bwtLevel > 0) {
        io.set_level(bwtLevel);   /* only coder_bwt_cm consumes this */
    }
    {
        CoderT coder(&io);
        coder.encode_line(data, len);
        coder.encode_flush();
    }
    outLen = (uint32_t)io.data_len;
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    /* A coder that overran its trial buffer is not a usable candidate. */
    return (outLen > 0 && outLen < cap && io.err == coder_io::IO_OK);
}

/* Trial-compress a stream with a coder chosen at runtime through the registry. Output
   buffer sized 2x input + margin. */
bool trialStream(CoderType type, const uint8_t* data, uint32_t len, int level,
                 uint32_t& outLen, uint32_t& usec)
{
    const auto t0 = std::chrono::steady_clock::now();
    uint32_t cap = (len << 1) + 65536;
    std::vector<uint8_t> outBuf(cap, 0);
    coder_io io(outBuf.data(), (int32_t)cap);
    if (level > 0) {
        io.set_level(level);
    }
    {
        std::shared_ptr<coder> c = CoderFactory::makeEncoder(type, &io);
        if (c == nullptr) {
            return false;
        }
        c->encode_line(data, len);
        c->encode_flush();
    }
    outLen = (uint32_t)io.data_len;
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    /* A coder that overran its trial buffer is not a usable candidate. */
    return (outLen > 0 && outLen < cap && io.err == coder_io::IO_OK);
}

/* Line-by-line trial for line-based coders (descriptor lineMode); preserving line
   boundaries is required to measure them faithfully. */
bool trialLines(CoderType type, const std::vector<LineSample>& lines,
                uint32_t& outLen, uint32_t& usec)
{
    if (lines.empty()) {
        return false;
    }

    /* Every sampled line, not just the probe's: a line-based coder is measured on the
       sampled lines it is handed either way. */
    const auto t0 = std::chrono::steady_clock::now();
    uint64_t total = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        total += lines[i].len;
    }
    uint32_t cap = (uint32_t)((total << 1) + 65536);
    std::vector<uint8_t> outBuf(cap, 0);
    coder_io io(outBuf.data(), (int32_t)cap);
    {
        std::shared_ptr<coder> c = CoderFactory::makeEncoder(type, &io);
        if (c == nullptr) {
            return false;
        }
        for (size_t i = 0; i < lines.size(); ++i) {
            c->encode_line(lines[i].data, lines[i].len);
        }
        c->encode_flush();
    }
    outLen = (uint32_t)io.data_len;
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    /* A coder that overran its trial buffer is not a usable candidate. */
    return (outLen > 0 && outLen < cap && io.err == coder_io::IO_OK);
}

/* Whether a field's candidates include a line-based coder that needs per-line views. */
bool fieldNeedsLineSample(uint32_t fieldIdx, uint8_t mode)
{
    const std::vector<CoderType>& candidates = samFieldCandidates(fieldIdx, mode);
    for (size_t i = 0; i < candidates.size(); ++i) {
        const CoderDescriptor* d = CoderFactory::descriptor(candidates[i]);
        if (d != nullptr && d->trialLineBased) {
            return true;
        }
    }
    return false;
}

} // namespace

int CodecSelector::pickBwtLevel(uint32_t sampleLen)
{
    /* Smallest level whose block size fits the whole sample in one block. */
    for (int lv = 1; lv <= 9; ++lv) {
        if (BWT_LEVEL_SIZE[lv] >= sampleLen) {
            return lv;
        }
    }
    return 9;
}

/* Sample bytes each candidate is measured on. A small probe is biased: coder_fc's
   long-range redundancy is invisible while bwt_cm's block is sized to the probe. */
const uint32_t PROBE_START = 256u << 10;

/* Probe size, overridable from the environment for measurement (PBGZ_TRIAL_PROBE). */
static uint32_t trialProbeStart()
{
    return PROBE_START;
}

/* How much size penalty the verdict will take for a much faster coder; the trial time
   only evidences order-of-magnitude differences. */
struct SpeedTrade {
    bool   enabled;
    double maxPenaltyPct;   /* extra output tolerated on the faster candidate, in percent */
    double minFactor;       /* how many times faster it has to be */
};

static const SpeedTrade& speedTrade()
{
    static const SpeedTrade t = { true, 1.0, 2.0 };
    return t;
}

/* Fixed verdict for a field whose fast path has exactly one coder: no trial is run, so no
   archive-profile coder is exercised during preprocessing. */
static FieldCodecSelection fixedSelection(CoderType type, uint32_t sampleLen)
{
    FieldCodecSelection sel;
    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = type;
    sel.sampleLen = sampleLen;
    sel.decidedLen = sampleLen;
    sel.bestCompLen = 0;
    sel.rounds = 1;
    return sel;
}

/* Winner among the trialled candidates: smallest output unless a clearly faster coder
   is close enough. outLen receives the chosen coder's measured length. */
static CoderType pickTrialWinner(const FieldCodecSelection& sel, uint32_t& outLen)
{
    CoderType best = CoderType::FC;
    uint32_t bestLen = UINT32_MAX;
    uint32_t bestUs = 0;
    for (uint32_t i = 0; i < sel.trialCount; ++i) {
        if (sel.trialLen[i] == 0) {
            continue;
        }
        if (sel.trialLen[i] < bestLen) {
            bestLen = sel.trialLen[i];
            best = sel.trialCoder[i];
            bestUs = sel.trialUs[i];
        }
    }
    outLen = (bestLen == UINT32_MAX) ? 0 : bestLen;

    const SpeedTrade& trade = speedTrade();
    if (!trade.enabled || bestLen == UINT32_MAX || bestUs == 0) {
        return best;
    }

    const double limit = (double)bestLen * (1.0 + trade.maxPenaltyPct / 100.0);
    CoderType fastest = best;
    uint32_t fastestLen = bestLen;
    uint32_t fastestUs = bestUs;
    for (uint32_t i = 0; i < sel.trialCount; ++i) {
        if (sel.trialLen[i] == 0 || sel.trialUs[i] == 0) {
            continue;
        }
        if ((double)sel.trialLen[i] > limit) {
            continue;
        }
        if (sel.trialUs[i] < fastestUs) {
            fastestUs = sel.trialUs[i];
            fastestLen = sel.trialLen[i];
            fastest = sel.trialCoder[i];
        }
    }
    if (fastest != best && (double)fastestUs * trade.minFactor <= (double)bestUs) {
        outLen = fastestLen;
        return fastest;
    }
    return best;
}

/* Sort candidates into registry trialPriority order so the verdict is order-independent
   (plain insertion sort; lists are tiny). */
static void orderCandidates(std::vector<CoderType>& candidates)
{
    for (size_t i = 1; i < candidates.size(); ++i) {
        for (size_t j = i; j > 0; --j) {
            const CoderDescriptor* a = CoderFactory::descriptor(candidates[j - 1]);
            const CoderDescriptor* b = CoderFactory::descriptor(candidates[j]);
            const uint8_t pa = (a != nullptr) ? a->trialPriority : 0xFF;
            const uint8_t pb = (b != nullptr) ? b->trialPriority : 0xFF;
            if (pa <= pb) {
                break;
            }
            CoderType tmp = candidates[j - 1];
            candidates[j - 1] = candidates[j];
            candidates[j] = tmp;
        }
    }
}

FieldCodecSelection CodecSelector::selectCoder(const uint8_t* data, uint32_t len,
                                               const std::vector<LineSample>* lines,
                                               const std::vector<CoderType>& candidates)
{
    FieldCodecSelection sel;
    sel.sampleLen = len;
    if (data == nullptr || len < MIN_SELECT_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    /* Trial each candidate on [0, probe) with a fresh instance; the candidate set comes
       from field_coder_config.h and CoderFactory. */
    std::vector<CoderType> ordered = candidates;
    orderCandidates(ordered);

    const uint32_t probeTarget = trialProbeStart();
    const uint32_t probe = (probeTarget < len) ? probeTarget : len;
    sel.rounds = 1;
    sel.decidedLen = probe;

    /* Coders that can actually be trialled, in priority order: creatable, a trial
       candidate, and line-based coders only when sample lines are present. */
    std::vector<CoderType> eligible;
    eligible.reserve(ordered.size());
    for (size_t i = 0; i < ordered.size(); ++i) {
        const CoderType type = ordered[i];
        const CoderDescriptor* d = CoderFactory::descriptor(type);
        if (d == nullptr || !d->trialCandidate || !CoderFactory::canMake(type)) {
            continue;
        }
        if (d->trialLineBased && (lines == nullptr || lines->empty())) {
            continue;
        }
        eligible.push_back(type);
    }

    /* One eligible coder means nothing to compare, so skip the trial and take it
       directly (bestCompLen 0, trialCount 0). */
    if (eligible.size() == 1) {
        sel.status = FieldStatus::SELECTED;
        sel.selectedCoder = eligible[0];
        return sel;
    }

    uint32_t bestLen = UINT32_MAX;
    bool haveBest = false;

    for (size_t i = 0; i < eligible.size(); ++i) {
        const CoderType type = eligible[i];
        const CoderDescriptor* d = CoderFactory::descriptor(type);
        bool ok = false;
        uint32_t outLen = 0;
        uint32_t outUs = 0;
        if (d->trialLineBased) {
            ok = trialLines(type, *lines, outLen, outUs);
        } else {
            /* trialUsesLevel: block level chosen to hold this sample (coder_bwt_cm). */
            const int level = d->trialUsesLevel ? pickBwtLevel(probe) : 0;
            ok = trialStream(type, data, probe, level, outLen, outUs);
        }

        sel.addTrial(type, ok ? outLen : 0u, outUs);
        if (ok && (!haveBest || outLen < bestLen)) {
            bestLen = outLen;
            haveBest = true;
        }
    }

    if (!haveBest) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    /* Pick the winner from recorded sizes and times together (see pickTrialWinner). */
    uint32_t chosenLen = 0;
    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = pickTrialWinner(sel, chosenLen);
    sel.bestCompLen = chosenLen;
    return sel;
}

FieldCodecSelection CodecSelector::selectCoder(const uint8_t* data, uint32_t len, bool trialAffix,
                                               const std::vector<LineSample>* lines, bool allowBwtCm)
{
    /* Historical candidate set: coder_fc; plus affix_match when trialAffix; plus bwt_cm when
       allowBwtCm. */
    std::vector<CoderType> candidates;
    candidates.push_back(CoderType::FC);
    if (trialAffix) {
        candidates.push_back(CoderType::AFFIX_MATCH);
    }
    if (allowBwtCm) {
        candidates.push_back(CoderType::BWT_CM);
    }
    return selectCoder(data, len, lines, candidates);
}

/* Rebuild the LEB128 delta-varint stream that compressPosFieldDelta feeds the POS coder
   (baseline reset on RNAME change), then trial bwt_cm and coder_arith on it. */
FieldCodecSelection CodecSelector::selectPosDeltaCoder(const std::vector<LineSample>& posLines,
                                                       const std::vector<LineSample>& chrLines,
                                                       std::vector<uint64_t>* varintCounts)
{
    FieldCodecSelection sel;

    if (varintCounts != nullptr) {
        varintCounts->assign(256, 0);
    }

    /* Rebuild the varint delta stream. */
    std::vector<uint8_t> varint;
    int64_t prevPos = 0;
    std::string prevChr;
    for (size_t i = 0; i < posLines.size(); ++i) {
        std::string chr((const char*)chrLines[i].data, chrLines[i].len - 1); /* drop trailing tab */
        if (i > 0 && !prevChr.empty() && chr != prevChr) {
            prevPos = 0;
        }
        prevChr = chr;

        int64_t pos = 0;
        bool valid = false;
        const uint8_t* p = posLines[i].data;
        uint32_t n = posLines[i].len - 1; /* drop trailing tab */
        if (n > 0) {
            int64_t v = 0;
            uint32_t k = 0;
            while (k < n && p[k] >= '0' && p[k] <= '9') {
                v = v * 10 + (p[k] - '0');
                ++k;
            }
            if (k == n) {
                pos = v;
                valid = true;
            }
        }

        int64_t delta = 0;
        if (valid) {
            delta = pos - prevPos;
            prevPos = pos;
        }
        uint64_t u = (uint64_t)delta;
        do {
            uint8_t b = (uint8_t)(u & 0x7f);
            u >>= 7;
            if (u) b |= 0x80;
            varint.push_back(b);
            if (varintCounts != nullptr)
                (*varintCounts)[b]++;
        } while (u);
    }

    sel.sampleLen = (uint32_t)varint.size();
    if (varint.size() < MIN_POS_DELTA_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }
    sel.decidedLen = (uint32_t)varint.size();
    sel.rounds = 1;

    uint32_t bwtLen = 0, bwtUs = 0, arithLen = 0, arithUs = 0, outLen = 0, outUs = 0;
    bool okBwt = trialEncode<coder_bwt_cm>(varint.data(), (uint32_t)varint.size(),
                                           pickBwtLevel((uint32_t)varint.size()), outLen, outUs);
    if (okBwt) {
        bwtLen = outLen;
        bwtUs = outUs;
    }
    bool okArith = trialEncode<coder_arith>(varint.data(), (uint32_t)varint.size(), 0, outLen, outUs);
    if (okArith) {
        arithLen = outLen;
        arithUs = outUs;
    }

    if (!okBwt && !okArith) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    sel.selectedCoder = CoderType::BWT_CM;
    sel.bestCompLen = bwtLen;
    if (okArith && arithLen < bwtLen) {
        sel.selectedCoder = CoderType::ARITH;
        sel.bestCompLen = arithLen;
    }
    sel.status = FieldStatus::SELECTED;
    sel.addTrial(CoderType::BWT_CM, bwtLen, bwtUs);
    sel.addTrial(CoderType::ARITH, arithLen, arithUs);
    return sel;
}

/* Append QUAL records from a block; supports single-shot and cross-block
   accumulation. qualBudget is the total cap, not the remaining budget. */
namespace {
void collectQualSamplesInto(RoughIOBlock* block,
                            std::vector<QualSampleRecord>& records,
                            std::vector<uint32_t>& freqByByte,
                            uint32_t sampleBudget,
                            uint64_t qualBudget,
                            uint64_t& collectedQualBytes)
{
    SafeLineReader reader(block);

    const uint8_t* line = nullptr;
    uint32_t lineLen = 0;
    while (reader.nextLine(line, lineLen)) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }
        if (lineLen == 0 || line[0] == '@') {
            continue;
        }

        /* Split on tabs to locate the first SAM_FIELD_COUNT fields. */
        uint32_t beg[SAM_FIELD_COUNT] = {0};
        uint32_t end[SAM_FIELD_COUNT] = {0};
        uint32_t fieldIdx = 0;
        uint32_t pos = 0;
        while (pos <= lineLen && fieldIdx < SAM_FIELD_COUNT) {
            uint32_t tabPos = pos;
            while (tabPos < lineLen && line[tabPos] != '\t') {
                ++tabPos;
            }
            beg[fieldIdx] = pos;
            end[fieldIdx] = tabPos;
            ++fieldIdx;
            pos = tabPos + 1;
        }
        if (fieldIdx < SAM_FIELD_COUNT) {
            continue;   /* incomplete fields: skip the line */
        }

        uint32_t qBeg = beg[SAM_QUAL], qEnd = end[SAM_QUAL];
        if (qEnd <= qBeg) {
            continue;
        }
        /* QUAL '*' means missing quality. */
        if (qEnd - qBeg == 1 && line[qBeg] == '*') {
            continue;
        }

        uint32_t qualLen = qEnd - qBeg;
        /* Only complete records (no truncation); over-cap records are skipped but the
           scan continues. */
        if (collectedQualBytes + qualLen > qualBudget) {
            continue;
        }

        for (uint32_t p = qBeg; p < qEnd; ++p) {
            freqByByte[line[p]]++;
        }

        QualSampleRecord rec;
        rec.qual.assign((const char*)(line + qBeg), qualLen);
        rec.seq.assign((const char*)(line + beg[SAM_SEQ]), end[SAM_SEQ] - beg[SAM_SEQ]);

        /* Parse FLAG manually for the 0x10 (reverse) bit; the field is not NUL-terminated. */
        long flagVal = 0;
        for (uint32_t p = beg[SAM_FLAG]; p < end[SAM_FLAG]; ++p) {
            uint8_t ch = line[p];
            if (ch < '0' || ch > '9') {
                break;
            }
            flagVal = flagVal * 10 + (ch - '0');
        }
        rec.rev = (flagVal & 16) != 0;

        records.push_back(rec);
        collectedQualBytes += qualLen;
    }
}

/* Append a FASTQ record's QUAL with its base line as context and no strand (FASTQ has
   no FLAG). Only complete records; over-cap records are dropped without stopping. */
void appendFastqQualRecord(const uint8_t* qualLine, uint32_t qualLen,
                           const uint8_t* baseLine, uint32_t baseLen,
                           std::vector<QualSampleRecord>& records,
                           std::vector<uint32_t>& freqByByte,
                           uint64_t& collectedQualBytes,
                           uint64_t qualBudget)
{
    if (qualLen == 0) {
        return;
    }
    if (collectedQualBytes + qualLen > qualBudget) {
        return;
    }
    for (uint32_t p = 0; p < qualLen; ++p) {
        freqByByte[qualLine[p]]++;
    }
    QualSampleRecord rec;
    rec.qual.assign((const char*)qualLine, qualLen);
    rec.seq.assign((const char*)baseLine, baseLen);
    rec.rev = false;
    records.push_back(rec);
    collectedQualBytes += qualLen;
}

/* FASTQ counterpart of collectQualSamplesInto: consumes records four lines at a time. */
void collectFastqQualSamplesInto(RoughIOBlock* block,
                                 std::vector<QualSampleRecord>& records,
                                 std::vector<uint32_t>& freqByByte,
                                 uint32_t sampleBudget,
                                 uint64_t qualBudget,
                                 uint64_t& collectedQualBytes)
{
    SafeLineReader reader(block);

    const uint8_t* lines[4];
    uint32_t lens[4];
    uint32_t slot = 0;
    while (reader.nextLine(lines[slot], lens[slot])) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }
        ++slot;
        if (slot < 4) {
            continue;
        }
        slot = 0;
        if (lens[0] == 0 || lines[0][0] != '@') {
            continue;
        }
        appendFastqQualRecord(lines[3], lens[3], lines[1], lens[1], records, freqByByte,
                              collectedQualBytes, qualBudget);
    }
}
} /* namespace */

void CodecSelector::extractQualSamples(RoughIOBlock* block,
                                       std::vector<QualSampleRecord>& records,
                                       std::vector<uint32_t>& freqByByte,
                                       uint32_t sampleBudget,
                                       uint64_t qualBudget)
{
    records.clear();
    freqByByte.assign(256, 0);
    uint64_t collected = 0;
    collectQualSamplesInto(block, records, freqByByte, sampleBudget, qualBudget, collected);
}

/* Append a block's QUAL records to the cross-block accumulator, capped at
   QUAL_PRIOR_TRAIN_MAX. */
void CodecSelector::accumulateQualPrior(RoughIOBlock* block, QualPriorAccum& acc)
{
    if (block == nullptr || block->getDataLen() <= 0 || acc.full()) {
        return;
    }
    if (acc.freqByByte.empty()) {
        acc.freqByByte.assign(256, 0);
    }
    /* The budget is the total cap, not a remaining quota (see collectQualSamplesInto). */
    if (BlockUtil::isFastqBlock(block->getBlockType())) {
        collectFastqQualSamplesInto(block, acc.records, acc.freqByByte, UINT32_MAX,
                                    QUAL_PRIOR_TRAIN_MAX, acc.collectedBytes);
        return;
    }
    collectQualSamplesInto(block, acc.records, acc.freqByByte, UINT32_MAX, QUAL_PRIOR_TRAIN_MAX, acc.collectedBytes);
}

/* Train fcv2 on the accumulated QUAL and export a model snapshot. Only produces the
   snapshot, so coder/ stays unaware of blocks and file offsets. */
std::vector<uint8_t> CodecSelector::trainQualPriorModel(const QualPriorAccum& acc,
                                                        const QualFcv2Params& params,
                                                        uint64_t* trainedBytes)
{
    uint64_t trained = 0;
    for (const QualSampleRecord& r : acc.records) {
        trained += r.qual.size();
    }
    if (trainedBytes != nullptr) {
        *trainedBytes = trained;
    }
    if (trained == 0) {
        return {};
    }

    try {
        std::vector<uint8_t> scratch((size_t)trained * 2 + (1u << 16), 0);
        coder_io io(scratch.data(), (int32_t)scratch.size());
        /* Translate the QualSelector tier into coder-layer parameters. */
        Fcv2Cfg cfg;
        cfg.cycleMax = params.cycleMax;
        cfg.cycleBucket = params.cycleBucket;
        cfg.prevShift = params.prevShift;
        cfg.useDedup = params.useDedup;
        coder_fcv2 coder(&io, acc.freqByByte, cfg);
        for (const QualSampleRecord& r : acc.records) {
            coder.encode_record((const uint8_t*)r.qual.data(),
                                (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        if (coder.encode_flush() <= 0) {
            return {};
        }

        std::vector<uint8_t> snapshot;
        if (!coder.export_model(snapshot) || snapshot.empty()) {
            return {};
        }
        return snapshot;
    } catch (...) {
        /* Optimization only; any failure falls back to the prior-less path. */
        return {};
    }
}

uint32_t CodecSelector::extractSamFieldSamples(RoughIOBlock* block,
                                           std::vector<std::string>& fieldBufs,
                                           std::vector<std::vector<LineSample>>& fieldLines,
                                           uint32_t sampleBudget,
                                           uint8_t mode)
{
    fieldBufs.assign(SAM_FIELD_COUNT_SELECT, std::string());
    fieldLines.assign(SAM_FIELD_COUNT_SELECT, std::vector<LineSample>());
    SafeLineReader reader(block);

    const uint8_t* line = nullptr;
    uint32_t lineLen = 0;
    while (reader.nextLine(line, lineLen)) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }
        if (lineLen == 0 || line[0] == '@') {
            continue;
        }

        uint32_t fieldIdx = 0;
        uint32_t pos = 0;
        while (pos <= lineLen && fieldIdx < SAM_FIELD_COUNT_SELECT) {
            uint32_t tabPos = pos;
            while (tabPos < lineLen && line[tabPos] != '\t') {
                ++tabPos;
            }
            fieldBufs[fieldIdx].append((const char*)(line + pos), (size_t)(tabPos - pos));
            /* Record line views for line-based trials and for the POS delta trial
               (POS + RNAME). Only pointers are stored; the trailing tab is included. */
            if (fieldNeedsLineSample(fieldIdx, mode) ||
                fieldIdx == (uint32_t)SAM_POS || fieldIdx == (uint32_t)SAM_RNAME) {
                LineSample ls;
                ls.data = line + pos;
                ls.len = (uint32_t)(tabPos - pos) + 1;
                fieldLines[fieldIdx].push_back(ls);
            }
            ++fieldIdx;
            pos = tabPos + 1;
        }
    }
    return reader.scannedBytes();
}

/* Collect RNAME+POS line views for the POS delta trial, bounded by line count so the
   trial sees block-like volume; the two vectors stay index-aligned. */
uint32_t CodecSelector::extractPosDeltaSamples(RoughIOBlock* block,
                                               std::vector<LineSample>& posLines,
                                               std::vector<LineSample>& chrLines,
                                               uint32_t maxLines)
{
    posLines.clear();
    chrLines.clear();
    if (maxLines == 0) {
        return 0;
    }

    SafeLineReader reader(block);

    const uint8_t* line = nullptr;
    uint32_t lineLen = 0;
    while (reader.nextLine(line, lineLen)) {
        if (posLines.size() >= maxLines) {
            break;
        }
        if (lineLen == 0 || line[0] == '@') {
            continue;
        }

        uint32_t pos = 0;
        uint32_t fieldIdx = 0;
        bool haveChr = false;
        LineSample chr;
        while (pos <= lineLen && fieldIdx <= (uint32_t)SAM_POS) {
            uint32_t tabPos = pos;
            while (tabPos < lineLen && line[tabPos] != '\t') {
                ++tabPos;
            }
            if (fieldIdx == (uint32_t)SAM_RNAME) {
                chr.data = line + pos;
                chr.len = (uint32_t)(tabPos - pos) + 1; /* trailing tab, as compression feeds it */
                haveChr = true;
            } else if (fieldIdx == (uint32_t)SAM_POS) {
                /* Malformed/truncated line: drop it so the two views stay aligned. */
                if (!haveChr) {
                    break;
                }
                LineSample p;
                p.data = line + pos;
                p.len = (uint32_t)(tabPos - pos) + 1;
                chrLines.push_back(chr);
                posLines.push_back(p);
                break;
            }
            ++fieldIdx;
            pos = tabPos + 1;
        }
    }
    return reader.scannedBytes();
}

uint32_t CodecSelector::extractFastqFieldSamples(RoughIOBlock* block,
                                             std::vector<std::string>& fieldBufs,
                                             std::vector<std::vector<LineSample>>& fieldLines,
                                             std::vector<QualSampleRecord>& qualRecords,
                                             std::vector<uint32_t>& qualFreq,
                                             uint32_t sampleBudget,
                                             uint64_t qualBudget)
{
    fieldBufs.assign(FQ_FIELD_COUNT, std::string());
    fieldLines.assign(FQ_FIELD_COUNT, std::vector<LineSample>());
    qualRecords.clear();
    qualFreq.assign(256, 0);
    SafeLineReader reader(block);

    const uint8_t* lines[4];
    uint32_t lens[4];
    uint32_t slot = 0;
    uint64_t collectedQual = 0;

    while (reader.nextLine(lines[slot], lens[slot])) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }

        ++slot;
        if (slot < 4) {
            continue;
        }

        /* A complete 4-line FASTQ record is now buffered in lines[0..3]. */
        const uint8_t* id = lines[0];
        const uint8_t* comment = lines[2];
        if (lens[0] > 0 && lens[2] > 0 && id[0] == '@' && comment[0] == '+') {
            /* Line views are the block's own bytes including the newline. */
            static const uint32_t kFieldOfLine[4] = {FQ_ID, FQ_SEQ, FQ_COMMENT, FQ_QUAL};
            for (uint32_t l = 0; l < 4; ++l) {
                LineSample s;
                s.data = lines[l];
                s.len = lens[l] + 1;   /* the trailing newline; see the header */
                fieldLines[kFieldOfLine[l]].push_back(s);
            }

            if (lens[0] > 1) fieldBufs[FQ_ID].append((const char*)(id + 1), (size_t)(lens[0] - 1));
            if (lens[1] > 0) fieldBufs[FQ_SEQ].append((const char*)lines[1], (size_t)lens[1]);
            if (lens[3] > 0) fieldBufs[FQ_QUAL].append((const char*)lines[3], (size_t)lens[3]);
            if (lens[2] > 1) fieldBufs[FQ_COMMENT].append((const char*)(comment + 1), (size_t)(lens[2] - 1));

            appendFastqQualRecord(lines[3], lens[3], lines[1], lens[1], qualRecords, qualFreq,
                                  collectedQual, qualBudget);
        }

        slot = 0;
    }
    return reader.scannedBytes();
}

/* Whether the QUAL prior is worth writing: a fixed cost versus a benefit that grows
   with total QUAL volume. Unknown input size writes; the threshold is the break-even. */
static bool qualPriorPaysOff(uint64_t qualSampleBytes, uint64_t scannedBytes, uint64_t inputTotalBytes)
{
    /* Break-even total QUAL volume, taken at the high end of the measured range so the
       prior is never written at a loss. Size trade-off only. */
    const uint64_t QUAL_PRIOR_MIN_TOTAL = 150ull * 1024ull * 1024ull;

    if (inputTotalBytes == 0 || scannedBytes == 0) {
        LOG_INFO("Input size unknown, keep QUAL prior.");
        return true;
    }

    const double qualShare = (double)qualSampleBytes / (double)scannedBytes;
    const uint64_t estimatedQual = (uint64_t)(qualShare * (double)inputTotalBytes);

    LOG_INFO("Estimated total QUAL %llu bytes (share %.4f of %llu), prior threshold %llu.",
             (unsigned long long)estimatedQual, qualShare,
             (unsigned long long)inputTotalBytes, (unsigned long long)QUAL_PRIOR_MIN_TOTAL);

    return estimatedQual >= QUAL_PRIOR_MIN_TOTAL;
}

/* Run the independent per-field trials in parallel on the idle cores. The verdict
   depends only on measured sizes, never on order or timing. */
void runTrialsInParallel(const std::vector<std::function<void()>>& tasks)
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
    /* The calling thread takes part too, otherwise one core stays unused. */
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

/* Decimal value of a field view, stopping at the first non-digit (FLAG, POS). */
static uint64_t parseDecimalField(const LineSample& field)
{
    uint64_t value = 0;
    for (uint32_t i = 0; i < field.len; ++i) {
        if (field.data[i] < '0' || field.data[i] > '9') {
            break;
        }
        value = value * 10 + (uint64_t)(field.data[i] - '0');
    }
    return value;
}

/* Candidate coders for a SAM field: its mode-specific list intersected with the
   registry's profile/level gate. Empty means a fixed strategy or none available. */
static std::vector<CoderType> eligibleCandidates(uint32_t fieldIdx, uint8_t mode, uint8_t level)
{
    /* The mode's configured list, plus whatever the extreme-ratio path (-l9) adds for this field. */
    std::vector<CoderType> configured = samFieldCandidates(fieldIdx, mode);
    const std::vector<CoderType> extra = samFieldExtremeCandidates(fieldIdx, mode, level);
    configured.insert(configured.end(), extra.begin(), extra.end());
    std::vector<CoderType> candidates;
    for (size_t c = 0; c < configured.size(); ++c) {
        const CoderType type = configured[c];
        const CoderDescriptor* desc = CoderFactory::descriptor(type);
        if (desc == nullptr || !desc->trialCandidate || !CoderFactory::canMake(type) ||
            !CoderFactory::eligible(type, mode, level)) {
            continue;
        }
        candidates.push_back(type);
    }
    return candidates;
}

/* The same intersection for a FASTQ field, read off kFastqFieldCoderConfig. */
static std::vector<CoderType> eligibleFastqCandidates(uint32_t fieldIdx, uint8_t mode, uint8_t level)
{
    std::vector<CoderType> configured = fastqFieldCandidates(fieldIdx, mode);
    const std::vector<CoderType> extra = fastqFieldExtremeCandidates(fieldIdx, mode, level);
    configured.insert(configured.end(), extra.begin(), extra.end());
    std::vector<CoderType> candidates;
    for (size_t c = 0; c < configured.size(); ++c) {
        const CoderType type = configured[c];
        const CoderDescriptor* desc = CoderFactory::descriptor(type);
        if (desc == nullptr || !desc->trialCandidate || !CoderFactory::canMake(type) ||
            !CoderFactory::eligible(type, mode, level)) {
            continue;
        }
        candidates.push_back(type);
    }
    return candidates;
}

/* Chromosome index of an RNAME view, or 0xFFFF for '*', '=', or an undeclared name. */
static uint16_t chrIndexOf(const LineSample& field)
{
    if (field.len == 0) {
        return 0xFFFF;
    }
    if (field.len == 1 && (field.data[0] == '*' || field.data[0] == '=')) {
        return 0xFFFF;
    }
    const uint16_t idx =
        SamInfo::getInstance().getChrNameIndex(std::string((const char*)field.data, field.len));
    return (idx != 65535) ? idx : 0xFFFF;
}

/* The SEQ match payload coder: coder_fc by default (single-entry list, no trial). At the
   extreme-ratio level (-l9) coder_bwt_cm joins it and the smaller of the two is taken.
   Returned as a FieldCodecSelection so the per-field table reports the choice. */
static FieldCodecSelection seqMatchCoderChoice(const uint8_t* sample, uint32_t segmentLen,
                                               uint8_t level)
{
    FieldCodecSelection sel;
    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = CoderType::FC;
    sel.sampleLen = segmentLen;
    sel.decidedLen = segmentLen;
    sel.rounds = 1;

    if (level < 9) {
        return sel;
    }

    const CoderType kCandidates[2] = {CoderType::BWT_CM, CoderType::FC};
    for (int i = 0; i < 2; ++i) {
        const CoderType type = kCandidates[i];
        /* coder_fc refuses inputs of FC_MIN_LEN bytes or less (and above its own maximum). */
        if (type == CoderType::FC &&
            (segmentLen <= FC_MIN_LEN || segmentLen >= (uint32_t)FC_MAX_LEN)) {
            continue;
        }
        uint32_t outLen = 0, usec = 0;
        if (!trialStream(type, sample, segmentLen, level, outLen, usec)) {
            continue;   /* this coder cannot carry the segment */
        }
        sel.addTrial(type, outLen, usec);
    }

    uint32_t chosenLen = 0;
    sel.selectedCoder = pickTrialWinner(sel, chosenLen);
    sel.bestCompLen = chosenLen;
    return sel;
}

uint32_t CodecSelector::extractSeqReferenceSamples(RoughIOBlock* block,
                                                   std::vector<SeqRefRecord>& records,
                                                   std::vector<uint8_t>& seq,
                                                   std::vector<uint32_t>& nPositions,
                                                   uint32_t sampleBudget)
{
    records.clear();
    seq.clear();
    nPositions.clear();
    if (block == nullptr || block->getBuffer() == nullptr) {
        return 0;
    }
    const uint8_t* buffer = block->getBuffer();
    const std::vector<size_t>& npos = block->getNpos();
    seq.reserve(sampleBudget);

    for (size_t lineIdx = 0; lineIdx < npos.size(); ++lineIdx) {
        if (seq.size() >= sampleBudget) {
            break;   /* sample full; the rest of the block is not read */
        }
        const size_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        const size_t lineEnd = npos[lineIdx];
        if (lineStart >= lineEnd || buffer[lineStart] == '@') {
            continue;   /* header line */
        }

        size_t colBeg[SAM_FIELD_COUNT];
        size_t colEnd[SAM_FIELD_COUNT];
        uint32_t col = 0;
        size_t p = lineStart;
        while (col < SAM_FIELD_COUNT && p < lineEnd) {
            colBeg[col] = p;
            while (p < lineEnd && buffer[p] != '\t') {
                ++p;
            }
            colEnd[col] = p;
            ++col;
            ++p;    /* step over the tab, or past the end of the line */
        }
        if (col < 10) {
            continue;   /* a record line carries at least ten columns */
        }

        const size_t seqBeg = colBeg[SAM_SEQ];
        const size_t seqEnd = colEnd[SAM_SEQ];
        if (seqEnd <= seqBeg) {
            continue;
        }
        const size_t avail = seqEnd - seqBeg;
        const size_t left = sampleBudget - seq.size();
        const size_t take = (avail < left) ? avail : left;

        SeqRefRecord rec;
        rec.seqBeg = (uint32_t)seq.size();
        rec.seqEnd = (uint32_t)(seq.size() + take);
        rec.flag.data = buffer + colBeg[SAM_FLAG];
        rec.flag.len = (uint32_t)(colEnd[SAM_FLAG] - colBeg[SAM_FLAG]);
        rec.rname.data = buffer + colBeg[SAM_RNAME];
        rec.rname.len = (uint32_t)(colEnd[SAM_RNAME] - colBeg[SAM_RNAME]);
        rec.pos.data = buffer + colBeg[SAM_POS];
        rec.pos.len = (uint32_t)(colEnd[SAM_POS] - colBeg[SAM_POS]);
        rec.cigar.data = buffer + colBeg[SAM_CIGAR];
        rec.cigar.len = (uint32_t)(colEnd[SAM_CIGAR] - colBeg[SAM_CIGAR]);

        const size_t base = seq.size();
        for (size_t i = 0; i < take; ++i) {
            const uint8_t ch = buffer[seqBeg + i];
            if (ch == 'N') {
                nPositions.push_back((uint32_t)(base + i));
            }
            seq.push_back(ch);
        }
        records.push_back(rec);
    }
    return (uint32_t)seq.size();
}

/* Compressed size of one form of a class's N positions (0 if uncodable); only the size
   is wanted, so the bytes go to a scratch buffer. */
static uint32_t nposFormSize(const uint8_t* data, uint32_t len, uint8_t level, uint8_t mode)
{
    if (len == 0) {
        return 0;
    }
    std::vector<uint8_t> scratch((size_t)len * 2 + 65536);
    coder_io trialIo(scratch.data(), (int32_t)scratch.size());
    const CoderType type = CoderFactory::fieldCoderForMode(SAM_SEQ, mode, CoderType::BWT_CM);
    CoderFactory::applyLevel(&trialIo, type, level);
    std::shared_ptr<coder> trialCoder = CoderFactory::makeEncoder(type, &trialIo);
    trialCoder->encode_line(data, len);
    trialCoder->encode_flush();
    return (trialIo.err == coder_io::IO_OK && trialIo.data_len > 0) ? (uint32_t)trialIo.data_len : 0;
}

/* Best of the three N-position forms (0 absolute, 1 deltas, 2 runs); ties keep the
   earlier form. `exc` must have keepForms set. */
static int32_t pickSeqExcForm(SeqExceptionClass& exc, uint8_t level, uint8_t mode)
{
    exc.closeRun();   /* close the open run so it is measurable */

    const uint32_t absSize = nposFormSize((const uint8_t*)exc.abs.data(),
                                          (uint32_t)(exc.abs.size() * sizeof(uint32_t)), level, mode);
    const uint32_t varSize = nposFormSize(exc.varint.data(), (uint32_t)exc.varint.size(), level, mode);
    std::vector<uint8_t> runBuf;
    runBuf.reserve((exc.runGaps.size() << 1) + 16);
    buildRunForm(exc, runBuf);
    const uint32_t runSize = nposFormSize(runBuf.data(), (uint32_t)runBuf.size(), level, mode);

    int32_t best = 0;
    uint32_t bestSize = absSize;
    if (varSize > 0 && (bestSize == 0 || varSize < bestSize)) {
        best = 1;
        bestSize = varSize;
    }
    if (runSize > 0 && (bestSize == 0 || runSize < bestSize)) {
        best = 2;
    }
    return best;
}

bool CodecSelector::selectSeqReferenceCoder(RoughIOBlock* block, Reference* reference, uint8_t level,
                                            uint8_t mode, PreprocessInfo& info)
{
    if (block == nullptr || reference == nullptr) {
        return false;
    }

    /* The @SQ lines RNAME resolves against; the reader thread usually already registered
       them, and re-parsing is a no-op. */
    SamUtil::prefetchHeaderChromosomes(block->getBuffer(), block->getDataLen(), block->getNpos());

    std::vector<SeqRefRecord> records;
    std::vector<uint8_t> seq;
    std::vector<uint32_t> nPositions;
    extractSeqReferenceSamples(block, records, seq, nPositions, kSeqRefSampleBytes);
    if (seq.empty()) {
        return false;
    }

    std::vector<uint8_t> match;
    match.reserve(seq.size());
    std::vector<uint8_t> ref2bit;
    std::vector<uint8_t> coded;
    std::vector<CigarOp> ops;
    uint32_t refCoded = 0;
    for (size_t r = 0; r < records.size(); ++r) {
        const SeqRefRecord& rec = records[r];
        const uint32_t seqLength = rec.seqEnd - rec.seqBeg;
        if (seqLength == 0) {
            continue;
        }
        const uint8_t* seqStart = seq.data() + rec.seqBeg;

        const uint16_t flag = (uint16_t)parseDecimalField(rec.flag);
        const uint64_t startPos = parseDecimalField(rec.pos);
        const uint16_t chrId = chrIndexOf(rec.rname);
        const uint32_t refConsumed = parseCigarOps(rec.cigar.data, rec.cigar.len, ops);

        /* Same payload branch as the block pass (buildSeqRecordPayload); N positions
           were counted over the whole sample. */
        coded.assign(seqLength, 0);
        ref2bit.assign(seqLength + 8, 0);
        const SeqRecordPayload payload =
            buildSeqRecordPayload(chrId, flag, startPos, refConsumed, ops, seqStart, seqLength,
                                  reference, coded.data(), ref2bit.data());
        if (payload.usedReference) {
            ++refCoded;
        }
        match.insert(match.end(), payload.bytes, payload.bytes + payload.length);
    }
    if (match.empty()) {
        return false;
    }

    /* Same RLE decision as the block pass, so the reported segment is the one encoded. */
    const SeqRleSplit rle = splitSeqMatchStream(match.data(), (uint32_t)match.size());
    const uint8_t* payload = rle.useRle ? rle.run.get() : match.data();
    const uint32_t payloadLen = rle.useRle ? rle.runLength : (uint32_t)match.size();
    if (payloadLen == 0) {
        return false;
    }

    const FieldCodecSelection verdict =
        (mode == PBGZ_MODE_FAST)
            ? fixedSelection(CoderType::RANS, payloadLen)
            : seqMatchCoderChoice(payload, payloadLen, level);
    info.fields[SAM_SEQ] = verdict;

    /* A caller-published verdict (>= 0) wins over this trial; -1 is undecided. */
    if (info.seqMatchCoder.load(std::memory_order_relaxed) < 0) {
        info.setSeqMatchCoder((int32_t)verdict.selectedCoder);
    }

    /* The smallest N-position form is a file property, so measure it once here. */
    if (!nPositions.empty() && info.nposForm.load(std::memory_order_relaxed) < 0) {
        SeqExceptionClass exc;
        exc.keepForms = true;   /* all three columns must be kept for comparison */
        exc.abs.reserve(nPositions.size() + 16);
        exc.varint.reserve(nPositions.size() + 64);
        for (size_t i = 0; i < nPositions.size(); ++i) {
            exc.add(nPositions[i]);
        }
        info.setNposForm(pickSeqExcForm(exc, level, mode));
    }

    LOG_DEBUG("SEQ payload: %u of %u records reference-coded, payload %u -> %s (match %u), "
              "npos form %d",
              refCoded, (unsigned)records.size(), payloadLen,
              coderTypeToMagic(verdict.selectedCoder), (uint32_t)match.size(),
              info.nposForm.load(std::memory_order_relaxed));
    return true;
}

bool CodecSelector::selectQnameLayout(RoughIOBlock* block, PreprocessInfo& info)
{
    if (block == nullptr || block->getBuffer() == nullptr) {
        return false;
    }
    const uint8_t* buffer = block->getBuffer();
    const std::vector<size_t>& npos = block->getNpos();
    if (npos.empty()) {
        return false;
    }

    /* Drive the ID analysis and tab collection the same way the block's own
       pre-analysis does (SamCodecActuator::scanDataLine). */
    IdSplitAnalysis analysis;
    std::vector<LineTabs> fieldTabs;
    size_t headEndLine = 0;
    while (headEndLine < npos.size()) {
        const size_t lineStart = (headEndLine == 0) ? 0 : npos[headEndLine - 1] + 1;
        if (lineStart >= npos[headEndLine] || buffer[lineStart] != '@') {
            break;
        }
        ++headEndLine;   /* the header is the block's leading lines */
    }

    bool sawFirstRecord = false;
    for (size_t lineIdx = headEndLine; lineIdx < npos.size(); ++lineIdx) {
        const size_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        const size_t lineEnd = npos[lineIdx];
        if (lineStart >= lineEnd || buffer[lineStart] == '@') {
            fieldTabs.emplace_back();
            continue;
        }
        LineTabs tabs;
        for (size_t i = lineStart; i < lineEnd; ++i) {
            if (buffer[i] == '\t') {
                tabs.push_back((int64_t)(i - lineStart));
            }
        }
        if (tabs.empty()) {
            return false;   /* not a record line: the block is not the SAM text this expects */
        }
        /* The QNAME field runs up to and including its tab. */
        const uint32_t qnameLen = (uint32_t)tabs[0] + 1;
        if (!sawFirstRecord) {
            analyzeQnameFirstLine(buffer + lineStart, qnameLen, analysis);
            sawFirstRecord = true;
        } else {
            analyzeQnameLine(buffer + lineStart, qnameLen, analysis);
        }
        fieldTabs.push_back(tabs);
    }

    if (!sawFirstRecord || analysis.symbols.empty() || analysis.posLength == UINT32_MAX) {
        /* Nothing to choose (no record, or the split is unavailable): the column is
           written whole. */
        return false;
    }

    /* Run both encoders as the block pass does, into a throwaway scratch block. */
    RoughIOBlock scratch(block->getBufferSize());
    scratch.setBlockType(block->getBlockType());
    coder_err_sink sink;

    QnameColumnInput in;
    in.buffer = buffer;
    in.npos = &npos;
    in.headEndLine = (int64_t)headEndLine;
    in.fieldTabs = &fieldTabs;

    Json::Value metaAffix;
    Json::Value metaQname;
    uint32_t srcAffix = 0;
    uint32_t srcQname = 0;
    const auto t0 = std::chrono::steady_clock::now();
    const int32_t lenAffix = encodeQnameSplit(in, analysis, &scratch, &sink, QNAME_TRIAL_LINES,
                                              metaAffix, srcAffix);
    const uint32_t affixUs = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
                                 std::chrono::steady_clock::now() - t0).count();
    scratch.setDataLen(0);
    const auto t1 = std::chrono::steady_clock::now();
    const int32_t lenQname = encodeQnameByQname(in, &scratch, &sink, QNAME_TRIAL_LINES, metaQname,
                                                srcQname);
    const uint32_t qnameUs = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
                                 std::chrono::steady_clock::now() - t1).count();
    scratch.setDataLen(0);

    /* Affix layout unless qname is smaller; a failure of one side leaves the other. */
    const bool useAffix = (lenQname < 0) || (lenAffix >= 0 && lenAffix <= lenQname);

    FieldCodecSelection sel;
    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = useAffix ? CoderType::AFFIX_MATCH : CoderType::QNAME;
    sel.sampleLen = srcAffix;      /* both layouts are handed the same lines */
    sel.decidedLen = srcAffix;
    sel.rounds = 1;
    if (lenAffix >= 0) {
        sel.addTrial(CoderType::AFFIX_MATCH, (uint32_t)lenAffix, affixUs);
    }
    if (lenQname >= 0) {
        sel.addTrial(CoderType::QNAME, (uint32_t)lenQname, qnameUs);
    }
    sel.bestCompLen = (uint32_t)(useAffix ? lenAffix : lenQname);
    info.fields[SAM_QNAME] = sel;
    info.setQnameUseAffix(useAffix ? 1 : 0);

    LOG_DEBUG("QNAME trial: affix=%d qname=%d (up to %u of the block's %u lines) -> %s", lenAffix,
              lenQname, QNAME_TRIAL_LINES, (unsigned)fieldTabs.size(),
              useAffix ? "affix" : "qname");
    return true;
}

int32_t CodecSelector::analyzeSam(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                                  uint8_t compressLevel, uint8_t mode, Reference* reference)
{
    std::vector<std::string> fieldBufs;
    std::vector<std::vector<LineSample>> fieldLines;
    info.scannedBytes = extractSamFieldSamples(block, fieldBufs, fieldLines, trialSampleTarget(), mode);

    info.fields.resize(SAM_FIELD_COUNT_SELECT);
    uint64_t totalSample = 0;

    /* Three phases: build one trial per field, run them in parallel, then apply the
       file-level consequences (QUAL/POS prior) serially. */
    std::vector<std::function<void()>> trials;
    std::vector<uint64_t> posCounts;
    uint64_t qualSampleBytes = 0;
    uint64_t posSampleBytes = 0;
    bool qualSelectedFcv2 = false;

    for (uint32_t f = 0; f < SAM_FIELD_COUNT_SELECT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();

        /* With a reference, SEQ is trialled on the match payload rather than its raw
           text; this trial also settles the N-position form. Falls back to the raw-text
           trial when it can take no verdict. */
        if (f == (uint32_t)SAM_SEQ && reference != nullptr) {
            const uint8_t* sampleData = (const uint8_t*)buf.data();
            const uint32_t sampleLen = (uint32_t)buf.size();
            const std::vector<LineSample>* sampleLines = &fieldLines[f];
            std::vector<CoderType> seqCandidates = eligibleCandidates(f, mode, compressLevel);
            trials.push_back([block, reference, compressLevel, mode, &info, f, sampleData,
                              sampleLen, sampleLines, seqCandidates]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + f);
                if (selectSeqReferenceCoder(block, reference, compressLevel, mode, info)) {
                    return;
                }
                if (seqCandidates.empty()) {
                    info.fields[f].status = FieldStatus::SKIPPED;
                    info.fields[f].sampleLen = sampleLen;
                } else {
                    info.fields[f] = selectCoder(sampleData, sampleLen, sampleLines, seqCandidates);
                }
            });
            continue;
        }

        /* QNAME is never written as raw text but as one of two layouts; the layout
           trial replaces the raw-text one and falls back to it on no verdict. */
        if (f == (uint32_t)SAM_QNAME) {
            const uint8_t* sampleData = (const uint8_t*)buf.data();
            const uint32_t sampleLen = (uint32_t)buf.size();
            const std::vector<LineSample>* sampleLines = &fieldLines[f];
            std::vector<CoderType> qnameCandidates = eligibleCandidates(f, mode, compressLevel);
            trials.push_back([block, mode, &info, f, sampleData, sampleLen, sampleLines,
                              qnameCandidates]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + f);
                /* Fast mode has one coder for QNAME: skip the affix/qname layout trial. */
                if (mode == PBGZ_MODE_FAST) {
                    info.fields[f] = fixedSelection(CoderType::RANS, sampleLen);
                    info.setQnameUseAffix(0);
                    return;
                }
                if (selectQnameLayout(block, info)) {
                    return;
                }
                if (qnameCandidates.empty()) {
                    info.fields[f].status = FieldStatus::SKIPPED;
                    info.fields[f].sampleLen = sampleLen;
                } else {
                    info.fields[f] = selectCoder(sampleData, sampleLen, sampleLines, qnameCandidates);
                }
            });
            continue;
        }

        /* POS is exempt from the raw-text threshold (its varint stream is ~1/4 of the
           column); the POS trial re-checks on MIN_POS_DELTA_SAMPLE. */
        if (buf.size() < MIN_SELECT_SAMPLE && f != (uint32_t)SAM_POS) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }

        if (f == (uint32_t)SAM_QUAL) {
            /* QUAL uses its dedicated record-level path (coder_qual/fcv2) on the whole
               block up to QUAL_TRIAL_MAX. */
            const uint32_t qualBudget = (uint32_t)block->getDataLen();
            const uint64_t bufSize = buf.size();
            trials.push_back([block, qualBudget, compressLevel, mode, &info, &qualSelectedFcv2,
                              &qualSampleBytes, bufSize]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + SAM_QUAL);
                /* Fast mode uses coder_rans for QUAL; the dedicated QUAL trial only
                   chooses among archive-profile coders (qual/fcv2/qcm). */
                if (mode == PBGZ_MODE_FAST) {
                    info.fields[SAM_QUAL] =
                        fixedSelection(CoderType::RANS, (uint32_t)bufSize);
                    qualSampleBytes = bufSize;
                    qualSelectedFcv2 = false;
                    return;
                }
                std::vector<QualSampleRecord> qualRecs;
                std::vector<uint32_t> qualFreq;
                {
                    PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_QUAL_EXTRACT);
                    extractQualSamples(block, qualRecs, qualFreq, qualBudget, QUAL_TRIAL_MAX);
                }
                {
                    PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_QUAL_SELECT);
                    info.fields[SAM_QUAL] = QualSelector::select(qualRecs, qualFreq, compressLevel);
                }
                qualSampleBytes = bufSize;
                qualSelectedFcv2 = (info.fields[SAM_QUAL].status == FieldStatus::SELECTED &&
                                    info.fields[SAM_QUAL].selectedCoder == CoderType::FCV2);
            });
            continue;
        }

        /* POS has its own trial on the rebuilt delta-varint stream
           (selectPosDeltaCoder), not on raw text. */
        if (f == (uint32_t)SAM_POS) {
            const uint32_t posSampleLen = (uint32_t)fieldBufs[SAM_POS].size();
            trials.push_back([block, compressLevel, mode, posSampleLen, &info, &fieldLines,
                              &posCounts, &posSampleBytes]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + SAM_POS);
                /* Fast mode has no POS candidate trial: the delta stream is coded with rans. */
                if (mode == PBGZ_MODE_FAST) {
                    info.fields[SAM_POS] =
                        fixedSelection(CoderType::RANS, posSampleLen);
                    posSampleBytes = posSampleLen;
                    return;
                }
                std::vector<LineSample> posLines, chrLines;
                uint32_t posTrialLines = BlockFactory::samReadsPerBlockOfLevel(compressLevel);
                /* Trial at block granularity, but cap the trial cost (see POS_TRIAL_MAX_LINES). */
                if (posTrialLines > POS_TRIAL_MAX_LINES) {
                    posTrialLines = POS_TRIAL_MAX_LINES;
                }
                extractPosDeltaSamples(block, posLines, chrLines, posTrialLines);
                if (posLines.size() > fieldLines[SAM_POS].size()) {
                    info.fields[SAM_POS] = selectPosDeltaCoder(posLines, chrLines, &posCounts);
                } else {
                    /* Small file: the generic sampler's views already cover the whole block. */
                    info.fields[SAM_POS] = selectPosDeltaCoder(fieldLines[SAM_POS],
                                                               fieldLines[SAM_RNAME], &posCounts);
                }
                posSampleBytes = info.fields[SAM_POS].sampleLen;
            });
            continue;
        }

        /* Candidate list captured by value; the trials run after this loop, in parallel. */
        const std::vector<CoderType> trialCandidates = eligibleCandidates(f, mode, compressLevel);

        if (trialCandidates.empty()) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }

        const uint8_t* sampleData = (const uint8_t*)buf.data();
        const uint32_t sampleLen = (uint32_t)buf.size();
        const std::vector<LineSample>* sampleLines = &fieldLines[f];
        trials.push_back([&info, f, sampleData, sampleLen, sampleLines, trialCandidates]() {
            PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + f);
            info.fields[f] = selectCoder(sampleData, sampleLen, sampleLines, trialCandidates);
        });
    }

    runTrialsInParallel(trials);

    /* Record the prior request only; training is deferred to CompressEngine. */
    if (qualSelectedFcv2 && qualPriorPaysOff(qualSampleBytes, info.scannedBytes, inputTotalBytes)) {
        info.setQualPriorRequested(true);
    }

    /* The prior only helps the arith backend; write it only when arith won POS and
       the estimated whole-file varint volume covers the ~512-byte prior. */
    if (info.fields[SAM_POS].status == FieldStatus::SELECTED &&
        info.fields[SAM_POS].selectedCoder == CoderType::ARITH && !posCounts.empty()) {
        uint64_t estimated = 0;
        if (inputTotalBytes > 0 && info.scannedBytes > 0) {
            estimated = (posSampleBytes * inputTotalBytes) / info.scannedBytes;
        } else {
            /* Unknown input size: conservatively assume the sample covers a small share. */
            estimated = posSampleBytes * 16;
        }
        if (estimated >= POS_PRIOR_MIN_ESTIMATED) {
            info.setPosPrior(packPosPriorBlob(posCounts));
            LOG_DEBUG("Preprocess SAM field %u: write POS prior (sample %llu varint bytes, est. total %llu)",
                      (uint32_t)SAM_POS, (unsigned long long)posSampleBytes,
                      (unsigned long long)estimated);
        }
    }

    for (uint32_t f = 0; f < SAM_FIELD_COUNT_SELECT; ++f) {
        LOG_DEBUG("Preprocess SAM field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }

    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyzeFastq(RoughIOBlock* block, uint64_t inputTotalBytes,
                                    PreprocessInfo& info, uint8_t compressLevel, uint8_t mode)
{
    std::vector<std::string> fieldBufs;
    std::vector<std::vector<LineSample>> fieldLines;
    std::vector<QualSampleRecord> qualRecords;
    std::vector<uint32_t> qualFreq;
    info.scannedBytes = extractFastqFieldSamples(block, fieldBufs, fieldLines, qualRecords,
                                                 qualFreq, trialSampleTarget(), QUAL_PRIOR_TRAIN_MAX);

    info.fields.resize(FQ_FIELD_COUNT);
    uint64_t totalSample = 0;
    bool qualSelectedFcv2 = false;
    uint64_t qualSampleBytes = 0;

    std::vector<std::function<void()>> trials;

    for (uint32_t f = 0; f < FQ_FIELD_COUNT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();

        /* QUAL uses its dedicated record-level path, as on the SAM side; it needs the
           read bases and record length, so it is not trialled as a byte stream. */
        if (f == (uint32_t)FQ_QUAL) {
            const uint64_t bufSize = buf.size();
            if (bufSize < MIN_SELECT_SAMPLE) {
                info.fields[f].status = FieldStatus::SKIPPED;
                info.fields[f].sampleLen = (uint32_t)buf.size();
                continue;
            }
            std::vector<CoderType> qualCandidates = fastqFieldCandidates(FQ_QUAL, mode);
            const std::vector<CoderType> qualExtra =
                fastqFieldExtremeCandidates(FQ_QUAL, mode, compressLevel);
            qualCandidates.insert(qualCandidates.end(), qualExtra.begin(), qualExtra.end());
            if (qualCandidates.empty()) {
                /* The mode's QUAL row lists no coder: nothing to compare, so mark it
                   skipped and let the row default apply. */
                info.fields[f].status = FieldStatus::SKIPPED;
                info.fields[f].sampleLen = (uint32_t)buf.size();
                continue;
            }
            trials.push_back([&qualRecords, &qualFreq, compressLevel, &info, &qualSelectedFcv2,
                              &qualSampleBytes, bufSize, qualCandidates]() {
                PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + FQ_QUAL);
                info.fields[FQ_QUAL] = QualSelector::select(qualRecords, qualFreq, compressLevel,
                                                            qualCandidates);
                qualSampleBytes = bufSize;
                qualSelectedFcv2 = (info.fields[FQ_QUAL].status == FieldStatus::SELECTED &&
                                    info.fields[FQ_QUAL].selectedCoder == CoderType::FCV2);
            });
            continue;
        }

        if (buf.size() < MIN_SELECT_SAMPLE) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }

        const std::vector<CoderType> trialCandidates = eligibleFastqCandidates(f, mode, compressLevel);
        if (trialCandidates.empty()) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }

        const uint8_t* sampleData = (const uint8_t*)buf.data();
        const uint32_t sampleLen = (uint32_t)buf.size();
        const std::vector<LineSample>* sampleLines = &fieldLines[f];
        trials.push_back([&info, f, sampleData, sampleLen, sampleLines, trialCandidates]() {
            PBGZ_PROF_SCOPE(pbgzprof::READ_TRIAL_BASE + f);
            info.fields[f] = selectCoder(sampleData, sampleLen, sampleLines, trialCandidates);
        });
    }

    runTrialsInParallel(trials);

    /* Same QUAL prior trade-off as the SAM side. */
    if (qualSelectedFcv2 && qualPriorPaysOff(qualSampleBytes, info.scannedBytes, inputTotalBytes)) {
        info.setQualPriorRequested(true);
    }

    for (uint32_t f = 0; f < FQ_FIELD_COUNT; ++f) {
        LOG_DEBUG("Preprocess FASTQ field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }

    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyze(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                               uint8_t compressLevel, uint8_t mode, Reference* reference)
{
    if (block == nullptr) {
        return -1;
    }
    BlockType type = block->getBlockType();
    info.reset(type);

    if (BlockUtil::isSAMBlock(type)) {
        return analyzeSam(block, inputTotalBytes, info, compressLevel, mode, reference);
    }
    if (BlockUtil::isFastqBlock(type)) {
        return analyzeFastq(block, inputTotalBytes, info, compressLevel, mode);
    }

    /* Unsupported type: leave info empty; actuators use their defaults. */
    return 0;
}
