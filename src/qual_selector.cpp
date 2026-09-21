/*
 * qual_selector.cpp - encoder evaluation for the quality-value column
 *
 * This file deliberately does not include coder_fc.h or any header that pulls in
 * fc/rangecoder.h: its RangeCoder conflicts with the same-named class in clr.h,
 * which coder_qual.h brings in indirectly. That is one reason the quality-value
 * evaluation is split out of codec_selector.cpp; see qual_selector.h for
 * details.
 *
 * Mirroring selectCoder for generic fields, the quality-value column also uses
 * multi-round convergence (strategy 7): start at 64 KB, double the sample each
 * round, and finalize once the leader opens a gap of more than 3%, avoiding
 * being misled on small samples by adaptive encoders that have not yet
 * converged. In addition, context parameter tiers are picked for fcv2 based on
 * data characteristics (strategy 2, following the fqz_pick_parameters idea from
 * fqzcomp), and the tiers participate as candidates in the trial compression;
 * the winning tier is handed to the compression side and prior training via
 * FieldCodecSelection.fcv2Params.
 */

#include "qual_selector.h"

#include <string.h>
#include <chrono>
#include <memory>
#include <atomic>
#include <functional>
#include <thread>

#include "coder/coder_io.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_qcm.h"
#include "field_coder_config.h"
#include "log/logger.h"

/* How many of m0..m6 the long-read mix keeps: the count is close to a pure speed
 * lever there (4 models is ~28% faster than 6 for ~0.1% more bytes, and 5 is not
 * smaller than 4), and fcv2's per-bit mix-and-update accounts for ~99% of the
 * value selection's wall time on long-read files. */
const int LONG_READ_MODEL_COUNT = 4;

/*
 * The model count is the one fcv2 knob that follows the data regime rather than
 * the block volume, so it is decided before the trial instead of inside it: the
 * count has to be fixed for every candidate of one trial, and the short- and
 * long-read regimes are far apart (~90 vs ~14000 bytes per record), so the
 * threshold only has to land inside that gap rather than nail a boundary. Short
 * reads keep the default mix because they are the regime fcv2 wins on, and
 * trimming m5 there costs 0.08%-0.46% of size.
 */
const uint32_t QUAL_LONG_READ_LEN = 500;

int qualModelCountForMeanLen(uint64_t meanLen)
{
    return (meanLen >= (uint64_t)QUAL_LONG_READ_LEN) ? LONG_READ_MODEL_COUNT : FCV2_MODEL_COUNT;
}

namespace {

/* If the sample is smaller than this many bytes, skip evaluation and keep the
 * default encoder. Too little sample data makes the decision unreliable. */
const uint32_t MIN_QUAL_SAMPLE = 64u << 10;

/* Slack for the trial-compression output buffer. Quality values are almost
 * never incompressible, so a 2x margin is more than safe. */
inline size_t trialCapacity(size_t srcLen)
{
    return (srcLen << 1) + (1u << 16);
}

/* How far the leader must be ahead for the winner to be settled; the same
 * threshold as the generic-field selectCoder. */
const double SETTLE_MARGIN = 0.03;

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

/*
 * Trial-compress with coder_qual.
 *
 * It uses the corresponding base sequence as context, so records must be fed
 * one at a time with both seq and qual. The frequency-table format follows the
 * existing convention in sam_actuator: an alphabet sorted by descending
 * frequency of occurrence, with the second element fixed at 1 (the actual count
 * is not preserved).
 */
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

/*
 * Trial-compress with fcv2.
 *
 * It does not need the base sequence, but it needs each record's length and
 * strand direction. The strand direction is written into the bitstream by the
 * encoder itself, so once passed in here the decoder does not need to supply it
 * again. cfg selects the context parameter tier; see strategy 2.
 */
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

/*
 * Trial-compress with coder_bwt_cm.
 *
 * It is in the QUAL candidate set (see field_coder_config.h) because it is the
 * best of the three on long-read quality values, and because it is what the
 * column falls back to where fcv2 is not applicable: fcv2 needs each record's
 * length and strand direction, which only the QUAL column of an aligned SAM can
 * provide (see CoderFactory::coderSupports).
 *
 * Records are fed one encode_line at a time, matching how sam_actuator actually
 * calls it. The internal block size is chosen per this round's sample size
 * (same rationale as pickBwtLevel in the generic-field selectCoder); small
 * samples do not need large blocks.
 */
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

/*
 * Trial-compress with coder_qcm.
 *
 * coder_qcm is an ordinary byte-stream coder: it takes the quality values as
 * they come and needs neither the record boundaries nor the strand, so the
 * sample is simply fed record by record (which is also how the QUAL stream
 * itself is written) and committed with one flush. The records are not
 * delimited anywhere in the stream, and the coder does not depend on where the
 * boundaries fall - the same records are fetched back one by one on the decode
 * side, where their lengths come from the already-decoded SEQ.
 */
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
    p.deltaMax = cfg.deltaMax;
    p.deltaBucket = cfg.deltaBucket;
    p.prevShift = cfg.prevShift;
    p.useDelta = cfg.useDelta;
    p.useDedup = cfg.useDedup;
    p.useQa = cfg.useQa;
    p.modelCount = cfg.modelCount;
    return p;
}

/*
 * Pick candidate tiers for fcv2 based on data characteristics (strategy 2).
 *
 * The idea comes from fqz_pick_parameters in fqzcomp: when the quality-value
 * alphabet is small (NovaSeq/HiSeqX-like, few symbols), the context is not
 * sparse, so finer positional bucketing can be used and the predecessor quality
 * values do not need quantization; when the alphabet is large or the sample is
 * small, the context is sparse, so bucketing is coarser and the predecessor is
 * shifted more.
 *
 * That statistical rule is NOT implemented as a rule. It is the rationale for
 * which tiers exist, but the winner is decided by measurement: several tiers are
 * put into the candidate set and select() trial-compresses them all, keeping the
 * smallest. The tie between the two is not reliable - the winner is not monotone
 * in the alphabet, and the preferred cycle span depends on the read length - so
 * the candidates are chosen to span the plausible winners rather than to encode a
 * prediction.
 *
 * The only data-characteristic decision that prunes rather than spans is the
 * block-volume tier: the fine tier needs enough sample to converge, and with a
 * small sample including it as a candidate only adds trial-compression time and
 * risks being wrongly chosen by a not-yet-converged illusion, so with small
 * samples only the always-kept tiers compete.
 */
/* compress level -> SAM data-block read-count tier (see
 * BlockFactory::createBlockReader): 1-5 -> ~10000 reads/block, 6-7 -> ~25000,
 * 8-9 -> ~100000. The per-block QUAL volume (reads x read length) decides how
 * fine the context buckets may be before they go sparse. */
static int levelVolumeTier(uint8_t compressLevel)
{
    if (compressLevel >= 8) return 3;
    if (compressLevel >= 6) return 2;
    return 1;
}

/*
 * fcv2 candidate tiers for the given compression level and model count.
 *
 * The model count is not a tier: every candidate in one trial carries the same
 * count (see qualModelCountForMeanLen), because it is a property of the read length
 * rather than of the block volume, and letting the candidates differ in it would
 * make one trial compare two different questions at once.
 *
 * The presets themselves are the FCV2_PRESETS table in coder_fcv2.h, so they can
 * be printed, measured and compared as data; this function only decides which of
 * them take part. That pruning is by the block-volume tier and by the mean read
 * length (FCV2_PRESETS[i].minMeanReadLen / maxMeanReadLen), so trial time is not
 * spent on presets that cannot win on this kind of data:
 *
 *   - default / ultra / ultra-80 / coarse: always kept. ultra is the coarse-binned
 *     variant (cycleBucket=4, prevShift=0), which leads at every measured block
 *     volume by 0.4%-1.2% over coarse; ultra-80 is the same over a shorter cycle
 *     span, coarse is the coarse tier with quantized predecessors, and default is
 *     the conservative baseline.
 *   - fine (24/12/ps0): needs a per-block QUAL volume large enough for its fine
 *     buckets to stay populated (it only wins at ~27 MB and up), so it is tried
 *     only at -l 8/9 (100k reads/block).
 *   - short-trained: bounded to its own read-length range (maxMeanReadLen), so
 *     trial time is not spent on 55-90 bp settings for a long-read block.
 *
 * Which of them actually wins is not predicted here; see the table's comment for
 * why a feature rule is not used.
 *
 * The per-read average-quality tier (qa) and duplicate-read dedup are not
 * presets: qa loses 1%-1.7% of ratio, dedup is neutral, and neither depends on
 * the block volume, so leaving them out cannot be recovered by a tier choice at
 * another volume.
 */
std::vector<Fcv2Cfg> candidateFcv2Cfgs(uint8_t compressLevel, int modelCount, uint64_t meanLen)
{
    std::vector<Fcv2Cfg> cfgs;
    const int volumeTier = levelVolumeTier(compressLevel);
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        if (volumeTier < FCV2_PRESETS[i].minVolumeTier) {
            continue;
        }
        if (FCV2_PRESETS[i].minMeanReadLen > 0 &&
            meanLen < (uint64_t)FCV2_PRESETS[i].minMeanReadLen) {
            continue;
        }
        if (FCV2_PRESETS[i].maxMeanReadLen > 0 &&
            meanLen >= (uint64_t)FCV2_PRESETS[i].maxMeanReadLen) {
            continue;
        }
        Fcv2Cfg cfg = FCV2_PRESETS[i].cfg;
        /*
         * A trained row carries the model count it was trained with (ownModelCount): the count and
         * the context parameters interact, so the pair only means what was measured when both
         * travel together - that is why it is not overwritten by the read-length rule here.
         */
        if (!FCV2_PRESETS[i].ownModelCount) {
            cfg.modelCount = modelCount;
        }
        cfgs.push_back(cfg);
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

/*
 * Run the trial-compressions on every core that is idle at this point.
 *
 * Codec pre-selection happens before any data block is dispatched, so all
 * worker threads are parked in workStartBarrier and the machine is empty. The
 * quality-value selection dominates the decision's wall time, and it is a ladder
 * of independent measurements: each rung (sample size) x each candidate coder
 * writes only its own slot.
 *
 * Every rung is executed even when an earlier one would already have settled.
 * That spends more CPU but no more wall time - the rungs run concurrently - and
 * the verdict is taken from the rung the ladder logic stops at, so running them
 * together cannot change the decision.
 */
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

    /*
     * coder_qcm is compared last, and the comparison is strict, so a tie keeps
     * coder_bwt_cm: it is the coder that has been in production longest, and a
     * tie means the newer model bought nothing measurable on this block.
     */
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

    /* The mean record length of the sample picks the model count (see
     * qualModelCountForMeanLen); the winning candidate - count included - is
     * handed to the compression side and to prior training through fcv2Params
     * below, and the decoder reads it back from the stream header. The same mean
     * length also prunes the candidate set (see the preset table). */
    const uint64_t meanLen = (uint64_t)sampleLen / (uint64_t)records.size();
    const int modelCount = qualModelCountForMeanLen(meanLen);

    /*
     * Which quality coders are in play is the config table's decision, not this
     * file's: only what the QUAL row of kSamFieldCoderConfig lists is trialled
     * (see qualCoderCandidate). An empty row therefore means "no selection" and
     * the column falls through to the field's fallback coder.
     */
    const bool tryQual = qualCoderCandidate(CoderType::QUAL);
    const bool tryFcv2 = qualCoderCandidate(CoderType::FCV2);
    const bool tryBwtCm = qualCoderCandidate(CoderType::BWT_CM);
    const bool tryQcm = qualCoderCandidate(CoderType::QCM);

    /* The candidate tiers are only built when fcv2 is in play; the slots below
     * keep the same shape either way, with every tier marked as not-tried. */
    const std::vector<Fcv2Cfg> cfgs =
        tryFcv2 ? candidateFcv2Cfgs(compressLevel, modelCount, meanLen) : std::vector<Fcv2Cfg>();

    /*
     * Multi-round convergence (strategy 7): start from a small sample and double
     * it each round, stopping once the leader opens enough of a gap. Each round
     * re-compresses all candidates at the new sample size rather than feeding
     * incrementally — re-compression mimics the "one data block from start to
     * finish" case of real compression, so the measured numbers are closer to
     * real behavior. The rationale is the same as the generic-field selectCoder;
     * see selectCoder in codec_selector.cpp.
     *
     * Early finalization has a precondition: the sample must be large enough
     * that the adaptive encoders have converged. On small samples bwt_cm
     * compresses unusually well and leads fcv2, a lead it loses only once the
     * sample grows, so finalizing at the lead threshold before that point would
     * wrongly pick bwt_cm. Hence a minimum finalization sample size: below it,
     * only double, never finalize.
     */
    const uint32_t MIN_SETTLE_PROBE = 1u << 20;   /* 1 MB */

    /*
     * Build the whole probe ladder up front and trial-compress every candidate
     * at every rung in parallel (see runInParallel). The ladder doubles the
     * sample from 64 KB until it is exhausted, and never finalizes below
     * MIN_SETTLE_PROBE, so that an adaptive encoder cannot win on a sample where
     * it has not converged yet.
     */
    std::vector<uint32_t> probes;
    {
        uint32_t p = (MIN_QUAL_SAMPLE < sampleLen) ? MIN_QUAL_SAMPLE : sampleLen;
        while (true) {
            probes.push_back(p);
            if (p >= sampleLen) {
                break;
            }
            p = (p > sampleLen / 2) ? sampleLen : (p * 2);
        }
    }

    std::vector<RoundSlots> slots(probes.size());
    for (size_t i = 0; i < probes.size(); i++) {
        slots[i].fcv2Ok.assign(cfgs.size(), false);
        slots[i].fcv2Len.assign(cfgs.size(), 0);
        slots[i].fcv2Us.assign(cfgs.size(), 0);
    }

    std::vector<std::function<void()>> tasks;
    tasks.reserve(probes.size() * (cfgs.size() + (size_t)tryQual + (size_t)tryBwtCm +
                                   (size_t)tryQcm));
    for (size_t i = 0; i < probes.size(); i++) {
        const uint32_t probe = probes[i];
        const size_t count = recordsForBudget(records, probe);
        if (tryQual) {
            tasks.push_back([&records, &freqByByte, &slots, i, count]() {
                slots[i].qualOk = trialQual(records, count, freqByByte, slots[i].qualLen, slots[i].qualUs);
            });
        }
        for (size_t c = 0; c < cfgs.size(); c++) {
            tasks.push_back([&records, &freqByByte, &slots, &cfgs, i, c, count]() {
                uint32_t len = 0;
                uint32_t us = 0;
                if (trialFcv2(records, count, freqByByte, cfgs[c], len, us)) {
                    slots[i].fcv2Ok[c] = true;
                    slots[i].fcv2Len[c] = len;
                    slots[i].fcv2Us[c] = us;
                }
            });
        }
        if (tryBwtCm) {
            const int bwtLevel = bwtLevelFor(probe);
            tasks.push_back([&records, &slots, i, count, bwtLevel]() {
                slots[i].cmOk = trialBwtCm(records, count, bwtLevel, slots[i].cmLen, slots[i].cmUs);
            });
        }
        if (tryQcm) {
            tasks.push_back([&records, &slots, i, count]() {
                slots[i].qcmOk = trialQcm(records, count, slots[i].qcmLen, slots[i].qcmUs);
            });
        }
    }
    runInParallel(tasks);

    QualRoundResult final;
    bool finalSet = false;
    uint32_t settleProbe = sampleLen;

    for (size_t i = 0; i < probes.size(); i++) {
        const uint32_t probe = probes[i];
        QualRoundResult r = assembleRound(slots[i], cfgs);
        final = r;
        finalSet = true;
        settleProbe = probe;
        sel.rounds++;

        if (!r.anyOk) {
            break;
        }
        if (probe >= sampleLen) {
            break;
        }
        if (probe < MIN_SETTLE_PROBE) {
            continue;   /* too early to trust the ranking, just keep doubling */
        }
        /* Only one candidate compresses successfully; more data cannot change
         * the comparison. */
        {
            /*
             * The runner-up is the smallest candidate that is not the leader. This
             * used to be an if-chain per leader; listing the candidates once keeps
             * it correct as candidates are added (coder_qcm joined the set).
             */
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
            /* The leader is far enough ahead that more data will not flip the
             * ranking. */
            if ((double)(runnerUp - r.bestLen) / (double)runnerUp >= SETTLE_MARGIN) {
                break;
            }
        }
    }

    sel.decidedLen = finalSet ? settleProbe : sampleLen;
    sel.trialCount = 0;
    /* Only the coders the config table put in play are reported; a coder that was
     * never trialled is absent rather than listed with a zero length. */
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
        /* fcv2 won; hand the winning tier from the last round to the caller so
         * encoding and prior training stay consistent with it. */
        sel.fcv2Params = toQualParams(final.fcv2Cfg);
    }
    sel.status = FieldStatus::SELECTED;

    /*
     * Name the winning preset, not just the coder: "picked=coder_fcv2" does not
     * say which of the presets was chosen, which is the thing worth seeing when
     * comparing a run against the table.
     */
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