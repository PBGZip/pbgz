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

/*
 * The fcv2 model count is fixed at FCV2_MODEL_COUNT (m0..m4) for every file,
 * short- or long-read: there is no per-read-length model-count rule anymore, so
 * every candidate of a trial encodes with the same full mix. Read length still
 * matters to the trial, but only as the threshold below which fcv2 is worth
 * trialling at all (QualSelector drops the coder above it). QUAL_LONG_READ_LEN is
 * that 500-byte boundary; the measured short- and long-read files sit at ~90 and
 * ~14000 bytes per record, so the exact value only has to land inside that gap.
 */
const uint32_t QUAL_LONG_READ_LEN = 500;

/* Whether a candidate list a caller passed in holds this coder. */
static bool contains(const std::vector<CoderType>& list, CoderType type)
{
    return std::find(list.begin(), list.end(), type) != list.end();
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

/* How far the leader must be ahead for the winner to be settled at the FIRST
 * settle rung (1 MB); the same threshold as the generic-field selectCoder. */
const double SETTLE_MARGIN = 0.03;

/* The margin at every larger rung (>= 2 MB). A lead measured on a small sample is
 * not trustworthy - coder_qual and coder_qcm can trade the lead, and a lead that
 * flips only at the full sample is exactly what the full ladder exists to catch - so
 * the first rung keeps the historical 3%. By 2 MB the lead has held across a full
 * doubling, which is the same evidence the ratio-plateau rule trusts, so a much
 * smaller margin settles it, saving the later rungs (whose samples keep doubling and
 * so dominate the cost) that the verdict would not have changed. */
const double LATER_SETTLE_MARGIN = 0.003;

/*
 * The margin used from the 2 MB rung up: LATER_SETTLE_MARGIN. The 1 MB rung always
 * uses SETTLE_MARGIN, so files that settle there are unaffected.
 */
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
    p.prevShift = cfg.prevShift;
    p.useDedup = cfg.useDedup;
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
 * fcv2 candidate tiers for the given compression level and mean read length.
 *
 * The model count is fixed at FCV2_MODEL_COUNT, so it is not a tier: every
 * candidate encodes with the same full mix, and a trial never mixes two
 * different questions.
 *
 * The presets themselves are the FCV2_PRESETS table in coder_fcv2.h, so they can
 * be printed, measured and compared as data; this function only decides which of
 * them take part. The only pruning left is the block-volume tier
 * (FCV2_PRESETS[i].minVolumeTier); read length is handled one level up, where the
 * whole fcv2 coder is dropped on long reads (see the prune in QualSelector::select):
 *
 *   - ultra (4/ps0) and ultra-80 (80/4/ps0): always kept. They lead the short-
 *     read ranking at every measured block volume.
 *   - fine (24/ps0): needs a per-block QUAL volume large enough for its fine
 *     buckets to stay populated (it leads from a few MB of per-block QUAL up), so
 *     it is tried only at -l 8/9 (100k reads/block).
 *
 * Which of them actually wins is not predicted here; see the table's comment for
 * why a feature rule is not used.
 *
 * Duplicate-read dedup is not a preset: dedup is neutral on typical data and
 * does not depend on the block volume, so leaving it out cannot be recovered by
 * a tier choice at another volume. (The per-read average-quality and transition
 * tiers are gone entirely with their models.)
 */
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

/*
 * Run the trial-compressions on every core that is idle at this point.
 *
 * Codec pre-selection happens before any data block is dispatched, so all
 * worker threads are parked in workStartBarrier and the machine is empty. The
 * quality-value selection dominates the decision's wall time, and it is a ladder
 * of independent measurements: each rung (sample size) x each candidate coder
 * writes only its own slot.
 *
 * What is handed over here is one rung's candidates, not the whole ladder: the
 * ladder is walked lazily (see the walk in selectQualCoder), and the candidates
 * of the rung being trialled are what run side by side.
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

    /* The mean record length of the sample prunes the fcv2 candidate set (see
     * the preset table). The winning preset is handed to the compression side
     * and to prior training through fcv2Params below; the decoder reads it back
     * from the stream header. */
    const uint64_t meanLen = (uint64_t)sampleLen / (uint64_t)records.size();

    /*
     * fcv2 is the short-read quality coder: its context mix (sequencing cycle,
     * base, predecessor quality) pays off on Illumina-style reads. Two cases make
     * trialling it pure cost, so it is dropped from the QUAL candidate set before
     * the trial:
     *
     *   - a quality alphabet wider than its Huffman tree cap (FCV2_TREE_CAP),
     *     where it cannot be built at all; and
     *   - long reads (mean length >= QUAL_LONG_READ_LEN), where the cycle context
     *     stops paying and coder_qual / coder_qcm win instead (measured on the
     *     long-read files here: on ERR11436629 and ERR10386218 every rung put
     *     coder_qcm first and fcv2 last).
     *
     * This trades a small ratio risk on data that straddles the threshold for the
     * fcv2 trial's cost; the candidate set is otherwise untouched.
     */
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

    /*
     * Which quality coders are in play is the config table's decision, not this
     * file's: only what the QUAL row lists is trialled (the caller passes the row's
     * list, so this file stays format-agnostic). An empty row therefore means "no
     * selection" and the column falls through to the field's fallback coder.
     */
    const bool tryQual  = contains(activeCandidates, CoderType::QUAL);
    const bool tryFcv2  = contains(activeCandidates, CoderType::FCV2);
    const bool tryBwtCm = contains(activeCandidates, CoderType::BWT_CM);
    const bool tryQcm   = contains(activeCandidates, CoderType::QCM);

    /* The candidate tiers are only built when fcv2 is in play; the slots below
     * keep the same shape either way, with every tier marked as not-tried. */
    const std::vector<Fcv2Cfg> cfgs =
        tryFcv2 ? candidateFcv2Cfgs(compressLevel) : std::vector<Fcv2Cfg>();

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
     * Enumerate the probe ladder up front (it is cheap - the sample doubles from
     * 64 KB until it is exhausted), then trial one rung at a time, that rung's
     * candidates in parallel (see runInParallel and the walk below). The ladder
     * never finalizes below MIN_SETTLE_PROBE, so that an adaptive encoder cannot
     * win on a sample where it has not converged yet.
     */
    std::vector<uint32_t> probes;
    {
        /*
         * Start the ladder at the oldest rung any settle rule can read, so the rungs
         * below it (whose compressed sizes nothing ever consults) are not computed.
         *
         * A verdict is only ever taken at a probe >= MIN_SETTLE_PROBE (1 MB), and the
         * margin rule reads only that rung. The plateau rule starts at
         * MIN_PLATEAU_PROBE (2 MB) and looks back exactly two rungs, so the oldest
         * rung it reads is one doubling below MIN_PLATEAU_PROBE/2, i.e.
         * MIN_SETTLE_PROBE/2 (512 KB) on the default ladder. The rungs below that -
         * 64 KB, 128 KB, 256 KB - are compressed and then read by nothing. Dropping
         * them keeps the settle points and the two-rung history the plateau rule
         * compares, so the recorded round sizes and the verdict are unchanged; it only
         * removes the work whose result was discarded. (If MIN_PLATEAU_PROBE stops
         * being 2 x MIN_SETTLE_PROBE, this start must be revisited.)
         */
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

    /* Slots for the rung being trialled; fresh per rung, so a trial that failed on this rung cannot
     * inherit the previous rung's success. */

    QualRoundResult final;
    bool finalSet = false;
    uint32_t settleProbe = sampleLen;

    /*
     * Ratio-plateau early exit (see the walk below). The ladder otherwise walks
     * to the full sample whenever the leader's margin never reaches SETTLE_MARGIN;
     * but the adaptive coders being compared (bwt_cm, qcm, fcv2) settle their
     * per-rung compression ratio after only a couple of MB, and once the top two
     * candidates' ratios have both stopped moving across consecutive doubling
     * rungs, more data cannot flip them without un-converging. Stop then, instead
     * of paying for the full-sample re-compressions the verdict does not read.
     */
    const uint32_t MIN_PLATEAU_PROBE = 2u << 20;   /* 2 MB */
    const double plPlateauBand = 0.002;            /* 0.2% absolute ratio movement */
    const bool plateauEnabled = true;
    const double laterSettleMargin = resolveLaterSettleMargin();
    QualRoundResult prevRound;
    QualRoundResult prevPrevRound;
    uint32_t prevProbe = 0;
    bool havePrev = false;
    bool havePrevPrev = false;

    /*
     * Walk the ladder one rung at a time: a rung is trialled only if the rung before it did not
     * settle. Every rung used to be trialled up front, on the argument that they run concurrently
     * and therefore cost no wall time - but this phase is what saturates the cores, so the rungs the
     * verdict never reads are plain wasted work, and they are the expensive ones (the sample doubles
     * from rung to rung, so the cost does too): on con_sorted.bam the decision settles at the 2 MB
     * rung while the 4 MB, 8 MB and full-sample rungs behind it - some 85% of the trialling - were
     * computed and discarded. Walking lazily computes exactly the rungs the decision reads, in the
     * same order and with the same settle rule, so the verdict is unchanged by construction.
     */
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
            const double settleMargin =
                (probe <= (1u << 20)) ? SETTLE_MARGIN : laterSettleMargin;
            if ((double)(runnerUp - r.bestLen) / (double)runnerUp >= settleMargin) {
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