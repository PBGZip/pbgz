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
#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_affix_match.h"
#include "field_coder_config.h"

namespace {

/*
 * Amount of raw line data scanned from the start of the first block.
 * 4MB is enough that the dominant fields (SEQ/QUAL ~ 37% each for typical
 * short reads) each yield a ~1.5MB sample, which is comfortably large enough
 * to separate the candidate coders.
 */
const uint32_t SAMPLE_TARGET = 4u << 20;   /* 4 MB */

/*
 * Minimum per-field sample required before we trust a codec comparison.
 * Each coder carries a small fixed overhead (BWT index, range-coder flush,
 * EOF marker, on the order of tens of bytes). At 64KB a real 1.7pp coder
 * difference is ~1.1KB, i.e. an order of magnitude above that overhead, so
 * the comparison is meaningful. Fields with a smaller sample are SKIPPED and
 * keep their default coder.
 */
const uint32_t MIN_SELECT_SAMPLE = 64u << 10;   /* 64 KB */

/* coder_bwt_cm internal block sizes per level, mirroring coder_bwt_cm. */
const uint32_t BWT_LEVEL_SIZE[10] = {
    0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
    1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
};

/*
 * Trial-compress a byte stream with one coder type and report the compressed
 * size. All coders follow the same encode_line()/encode_flush() contract and
 * write into coder_io::data, so a single template covers them.
 *
 * The output buffer is sized generously (2x input + margin): arithmetic coding
 * never expands much, but coder_fc writes up to 2x input internally, and the
 * coder_io sink is not bounds-checked.
 */
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
    return (outLen > 0 && outLen < cap);
}

/*
 * Trial-compress coder_affix_match line by line. affix's prefix/suffix matching
 * happens between adjacent lines, so it must be fed line by line; feeding the
 * whole column in one encode_line call degrades it to an ordinary context model
 * and fails to measure its real performance.
 */
bool trialAffixLines(const std::vector<LineSample>& lines, uint32_t& outLen, uint32_t& usec)
{
    if (lines.empty()) {
        return false;
    }
    const auto t0 = std::chrono::steady_clock::now();
    uint64_t total = 0;
    for (size_t i = 0; i < lines.size(); ++i) {
        total += lines[i].len;
    }
    uint32_t cap = (uint32_t)((total << 1) + 65536);
    std::vector<uint8_t> outBuf(cap, 0);
    coder_io io(outBuf.data(), (int32_t)cap);
    {
        coder_affix_match coder(&io);
        for (size_t i = 0; i < lines.size(); ++i) {
            coder.encode_line(lines[i].data, lines[i].len);
        }
        coder.encode_flush();
    }
    outLen = (uint32_t)io.data_len;
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return (outLen > 0 && outLen < cap);
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

/*
 * How far the leader must be ahead for the winner to be settled.
 *
 * Uses the relative difference in compressed size: if the leader is more than
 * 3% smaller than the runner-up, finalize. If the gap is below this, re-testing
 * with more samples could flip the ranking, so another round is worthwhile;
 * above this the ranking is unlikely to flip.
 */
const double SETTLE_MARGIN = 0.03;

/* Starting sample size for the evaluation. Too small a sample is easily biased
 * by local data, so it matches the minimum trustworthy sample size. */
const uint32_t PROBE_START = MIN_SELECT_SAMPLE;

FieldCodecSelection CodecSelector::selectCoder(const uint8_t* data, uint32_t len, bool trialAffix,
                                               const std::vector<LineSample>* lines)
{
    FieldCodecSelection sel;
    sel.sampleLen = len;
    if (data == nullptr || len < MIN_SELECT_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    /*
     * Start from a small sample and double it each round, stopping once the
     * leader opens up enough of a gap.
     *
     * Each round trial-compresses the fresh range [0, probe) with a new
     * instance, rather than incrementally feeding on top of the previous round.
     * There are two reasons: coder_fc only supports whole-block compression, so
     * a second encode_line call fails outright; and re-compressing mimics how a
     * block is actually compressed in production ("one data block from start to
     * finish"), so the measured numbers are closer to real behavior.
     *
     * The extra cost of re-compression is bounded: the sample doubles each
     * round, so all rounds together cost no more than twice the last round.
     * And most fields settle in the first round or two, so in practice this
     * costs less than compressing the entire sample in one shot.
     */
    uint32_t bwtCmLen = 0, fcLen = 0, affixLen = 0;
    uint32_t bwtCmUs = 0, fcUs = 0, affixUs = 0, outUs = 0, outLen = 0;
    uint32_t bestLen = UINT32_MAX;
    CoderType bestCoder = CoderType::BWT_CM;
    bool anyOk = false;
    uint32_t probe = (PROBE_START < len) ? PROBE_START : len;

    while (true) {
        bwtCmLen = fcLen = 0;
        bestLen = UINT32_MAX;
        bestCoder = CoderType::BWT_CM;
        anyOk = false;
        uint32_t runnerUp = UINT32_MAX;
        sel.rounds++;

        if (trialEncode<coder_bwt_cm>(data, probe, pickBwtLevel(probe), outLen, outUs)) {
            bwtCmLen = outLen; bwtCmUs = outUs;
            bestLen = outLen; bestCoder = CoderType::BWT_CM;
            anyOk = true;
        }
        if (trialEncode<coder_fc>(data, probe, 0, outLen, outUs)) {
            fcLen = outLen; fcUs = outUs;
            if (outLen < bestLen) { runnerUp = bestLen; bestLen = outLen; bestCoder = CoderType::FC; }
            else if (outLen < runnerUp) { runnerUp = outLen; }
            anyOk = true;
        }

        if (!anyOk) {
            break;
        }
        /* The sample is exhausted; there is no more data to add, so settle. */
        if (probe >= len) {
            break;
        }
        /* Only one candidate compresses successfully; adding more data gives
         * nothing to compare against. */
        if (runnerUp == UINT32_MAX) {
            break;
        }
        /* The leader is far enough ahead that more data will not flip the
         * ranking. */
        if ((double)(runnerUp - bestLen) / (double)runnerUp >= SETTLE_MARGIN) {
            break;
        }
        probe = (probe > len / 2) ? len : (probe * 2);
    }

    /*
     * affix's trial compression is done separately: it is a line-based encoder
     * that needs line boundaries, so it cannot be merged into the whole-stream
     * trial above. It compresses all line samples once (without per-round
     * doubling) and is compared against the settled bwt/fc result by size. If
     * lines is empty, the caller supplied no line samples, so affix simply
     * abstains.
     */
    if (trialAffix && lines != nullptr && !lines->empty() &&
        trialAffixLines(*lines, outLen, outUs)) {
        affixLen = outLen; affixUs = outUs;
        if (affixLen < bestLen) {
            bestLen = affixLen;
            bestCoder = CoderType::AFFIX_MATCH;
        }
    }

    sel.decidedLen = probe;
    sel.trialCount = 0;
    sel.addTrial(CoderType::BWT_CM, bwtCmLen, bwtCmUs);
    sel.addTrial(CoderType::FC, fcLen, fcUs);
    if (trialAffix && affixLen > 0) {
        sel.addTrial(CoderType::AFFIX_MATCH, affixLen, affixUs);
    }

    if (!anyOk && bestCoder != CoderType::AFFIX_MATCH) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = bestCoder;
    sel.bestCompLen = bestLen;
    return sel;
}

/*
 * Append QUAL records collected from a block. records/freqByByte are not
 * cleared, and collectedQualBytes is a running total, so this supports both
 * "clear then collect once" (extractQualSamples) and cross-block accumulation
 * (accumulateQualPrior). qualBudget is the **total cap** (not the remaining
 * budget for this call): the internal bounds check is
 * collectedQualBytes + qualLen > qualBudget, and since collectedQualBytes is
 * the accumulated total when spanning blocks, passing a remaining budget would
 * wrongly skip later blocks.
 */
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

        /* Split on each tab to locate the start and end of the first 11
         * required fields. */
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
            continue;   /* skip lines with incomplete fields; they do not take part in evaluation */
        }

        uint32_t qBeg = beg[SAM_QUAL], qEnd = end[SAM_QUAL];
        if (qEnd <= qBeg) {
            continue;
        }
        /* A QUAL of a single '*' means the value is missing; such records have
         * no quality values to compress. */
        if (qEnd - qBeg == 1 && line[qBeg] == '*') {
            continue;
        }

        uint32_t qualLen = qEnd - qBeg;
        /*
         * Only accept complete records; do not truncate QUAL to fill the cap:
         * fcv2 relies on the real read length to recover the cycle position,
         * and truncation would silently train the wrong context. Records beyond
         * the cap are not kept, but scanning continues to the end of the block.
         */
        if (collectedQualBytes + qualLen > qualBudget) {
            continue;
        }

        QualSampleRecord rec;
        rec.qual.assign((const char*)(line + qBeg), qualLen);
        rec.seq.assign((const char*)(line + beg[SAM_SEQ]), end[SAM_SEQ] - beg[SAM_SEQ]);

        /* Manually parse FLAG and extract the 0x10 bit. The field is not
         * NUL-terminated, so strtol cannot be used directly. */
        long flagVal = 0;
        for (uint32_t p = beg[SAM_FLAG]; p < end[SAM_FLAG]; ++p) {
            uint8_t ch = line[p];
            if (ch < '0' || ch > '9') {
                break;
            }
            flagVal = flagVal * 10 + (ch - '0');
        }
        rec.rev = (flagVal & 16) != 0;

        for (uint32_t p = qBeg; p < qEnd; ++p) {
            freqByByte[line[p]]++;
        }
        records.push_back(rec);
        collectedQualBytes += qualLen;
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

/*
 * Append a block's QUAL records to the cross-block training accumulator. Once
 * QUAL_PRIOR_TRAIN_MAX is reached no more are collected; within a single block
 * only complete records are kept, and over-cap records are dropped without
 * stopping the scan early.
 */
void CodecSelector::accumulateQualPrior(RoughIOBlock* block, QualPriorAccum& acc)
{
    if (block == nullptr || block->getDataLen() <= 0 || acc.full()) {
        return;
    }
    if (acc.freqByByte.empty()) {
        acc.freqByByte.assign(256, 0);
    }
    /*
     * The budget is passed as the total cap QUAL_PRIOR_TRAIN_MAX: the bounds
     * check inside collectQualSamplesInto is
     * collectedQualBytes + qualLen > qualBudget, and collectedQualBytes is the
     * cross-block accumulated value, so passing a "remaining quota" would let
     * the accumulated amount exceed the remaining one and skip subsequent
     * blocks entirely.
     */
    collectQualSamplesInto(block, acc.records, acc.freqByByte, UINT32_MAX, QUAL_PRIOR_TRAIN_MAX, acc.collectedBytes);
}

/*
 * Train the fcv2 prior on the accumulated QUAL and export a model snapshot.
 * Codec selection still uses only a small sample: a small sample suffices to
 * compare candidates stably, while model learning needs many more quality
 * values. This deliberately only produces a snapshot and touches neither blocks
 * nor file offsets, keeping the coder/ layer unaware of block and file offsets
 * and preserving the existing layering. Returns an empty vector on failure or
 * when the samples are empty.
 */
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
        /* Translate the tier selected by QualSelector into coder-layer
         * parameters so the prior stays consistent with the compression side. */
        Fcv2Cfg cfg;
        cfg.cycleMax = params.cycleMax;
        cfg.cycleBucket = params.cycleBucket;
        cfg.deltaMax = params.deltaMax;
        cfg.deltaBucket = params.deltaBucket;
        cfg.prevShift = params.prevShift;
        cfg.useDelta = params.useDelta;
        cfg.useDedup = params.useDedup;
        cfg.useQa = params.useQa;
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
        /* Training is only an optional optimization; any allocation or
         * unknown-data failure must fall back to the existing prior-less
         * compression path. */
        return {};
    }
}

uint32_t CodecSelector::extractSamFieldSamples(RoughIOBlock* block,
                                           std::vector<std::string>& fieldBufs,
                                           std::vector<std::vector<LineSample>>& fieldLines,
                                           uint32_t sampleBudget)
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
            /*
             * affix's line-by-line trial compression needs line boundaries,
             * which are lost when the whole column is concatenated. Only the
             * pointers are recorded, no copy: SafeLineReader returns views into
             * the block buffer, and the block contents stay unchanged during
             * analyze. Including the trailing tab matches how compression feeds
             * the data (see compressRegularField).
             */
            if (samFieldCandidate(fieldIdx, CoderType::AFFIX_MATCH)) {
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

uint32_t CodecSelector::extractFastqFieldSamples(RoughIOBlock* block,
                                             std::vector<std::string>& fieldBufs,
                                             uint32_t sampleBudget)
{
    fieldBufs.assign(FQ_FIELD_COUNT, std::string());
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

        /* A complete 4-line FASTQ record is now buffered in lines[0..3]. */
        const uint8_t* id = lines[0];
        const uint8_t* comment = lines[2];
        if (lens[0] > 0 && lens[2] > 0 && id[0] == '@' && comment[0] == '+') {
            if (lens[0] > 1) fieldBufs[FQ_ID].append((const char*)(id + 1), (size_t)(lens[0] - 1));
            if (lens[1] > 0) fieldBufs[FQ_SEQ].append((const char*)lines[1], (size_t)lens[1]);
            if (lens[3] > 0) fieldBufs[FQ_QUAL].append((const char*)lines[3], (size_t)lens[3]);
            if (lens[2] > 1) fieldBufs[FQ_COMMENT].append((const char*)(comment + 1), (size_t)(lens[2] - 1));
        }

        slot = 0;
    }
    return reader.scannedBytes();
}

/*
 * Whether the prior is worth writing.
 *
 * The cost of the prior is fixed: one auxiliary block, measured in this file at
 * about 0.6 MB packed, independent of input size. The benefit is that every
 * block's QUAL saves a roughly fixed fraction, growing linearly with the total
 * QUAL volume. A constant cost line and a benefit line through the origin must
 * intersect; the intersection is the break-even point, and below it writing the
 * prior is always a loss.
 *
 * The threshold comes from the measurement in benchmark/HANDSOFF.md section
 * 20.14: QUAL of about 150 MB (corresponding to SAM of about 540 MB). Before
 * this check, the code effectively assumed every file sits to the right of the
 * break-even point.
 *
 * The total is extrapolated from the share of QUAL in the scanned raw bytes of
 * the first block's sample. That share is determined by the record format (read
 * length, field composition) and is fairly stable within a file; the measured
 * extrapolation error is about 1%, while the decision compares against a
 * threshold an order of magnitude away, so this precision is sufficient.
 *
 * When the input length is unknown (piped input), treat it as "write": the cost
 * is bounded and small, while the loss from omitting it grows unboundedly with
 * file size, so the two kinds of misjudgment are not symmetric.
 */
static bool qualPriorPaysOff(uint64_t qualSampleBytes, uint64_t scannedBytes, uint64_t inputTotalBytes)
{
    /* Break-even point for the total QUAL volume. */
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

int32_t CodecSelector::analyzeSam(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info)
{
    std::vector<std::string> fieldBufs;
    std::vector<std::vector<LineSample>> fieldLines;
    info.scannedBytes = extractSamFieldSamples(block, fieldBufs, fieldLines, SAMPLE_TARGET);

    info.fields.resize(SAM_FIELD_COUNT_SELECT);
    uint64_t totalSample = 0;
    for (uint32_t f = 0; f < SAM_FIELD_COUNT_SELECT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();
        if (buf.size() < MIN_SELECT_SAMPLE) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }
        if (f == (uint32_t)SAM_QUAL) {
            /*
             * The quality-value column goes through dedicated evaluation: the
             * candidates are coder_qual and fcv2, not generic byte-stream
             * compressors. Previously this compared coder_bwt_cm and coder_fc,
             * while compression used the other two, so the selection result was
             * never actually used.
             */
            std::vector<QualSampleRecord> qualRecs;
            std::vector<uint32_t> qualFreq;
            /*
             * Quality-value selection scans the whole block instead of 4 MB:
             * fcv2 is an adaptive context mixer and is underestimated on small
             * samples before it converges (measured 48.27% on a 970 KB sample;
             * at 7.7 MB, 48.33%, already beating bwt_cm). Only by letting it
             * compete at real data volumes does the chosen coder match actual
             * compression.
             */
            extractQualSamples(block, qualRecs, qualFreq,
                               (uint32_t)block->getDataLen(), QUAL_PRIOR_TRAIN_MAX);
            info.fields[f] = QualSelector::select(qualRecs, qualFreq);
            /*
             * If the prior is worth training, record the request; the actual
             * training is deferred until the read thread finishes the
             * pre-training blocks and the cross-block accumulation is complete
             * (see CompressEngine). This only produces the decision; training is
             * not done here.
             */
            if (info.fields[f].status == FieldStatus::SELECTED &&
                info.fields[f].selectedCoder == CoderType::FCV2 &&
                qualPriorPaysOff(buf.size(), info.scannedBytes, inputTotalBytes)) {
                info.setQualPriorRequested(true);
            }
            continue;
        }
        /*
         * Which candidate coders are tried for a field is decided by the config
         * table (field_coder_config.h). An empty candidate list means the field
         * uses a fixed strategy (PNEXT/TLEN differencing/inference) and does not
         * take part in the generic selection.
         */
        const FieldCoderConfig* cfg = samFieldCoderConfig(f);
        if (cfg == nullptr || cfg->candidates.empty()) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }
        const bool trialAffix = samFieldCandidate(f, CoderType::AFFIX_MATCH);
        info.fields[f] = selectCoder((const uint8_t*)buf.data(), (uint32_t)buf.size(), trialAffix,
                                     trialAffix ? &fieldLines[f] : nullptr);
        LOG_DEBUG("Preprocess SAM field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }
    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyzeFastq(RoughIOBlock* block, PreprocessInfo& info)
{
    std::vector<std::string> fieldBufs;
    info.scannedBytes = extractFastqFieldSamples(block, fieldBufs, SAMPLE_TARGET);

    info.fields.resize(FQ_FIELD_COUNT);
    uint64_t totalSample = 0;
    for (uint32_t f = 0; f < FQ_FIELD_COUNT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();
        if (buf.size() < MIN_SELECT_SAMPLE) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }
        info.fields[f] = selectCoder((const uint8_t*)buf.data(), (uint32_t)buf.size());
        LOG_DEBUG("Preprocess FASTQ field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }
    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyze(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info)
{
    if (block == nullptr) {
        return -1;
    }
    BlockType type = block->getBlockType();
    info.reset(type);

    if (BlockUtil::isSAMBlock(type)) {
        return analyzeSam(block, inputTotalBytes, info);
    }
    if (BlockUtil::isFastqBlock(type)) {
        return analyzeFastq(block, info);
    }

    /* Unsupported type: leave info empty; actuators use their defaults. */
    return 0;
}
