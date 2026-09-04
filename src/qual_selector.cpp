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

#include "coder/coder_io.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_bwt_cm.h"
#include "log/logger.h"

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
 * This candidate was originally not in the quality-value candidate set — generic
 * fields only try bwt_cm through CodecSelector, while the quality-value column
 * chose only between coder_qual and fcv2. Measurement showed this was a real
 * gap: on 4 MB of real quality values, fcv2 25.81%, bwt_cm 26.31%, coder_qual
 * 33.74%; fcv2 still wins, but bwt_cm beats coder_qual by 7.4 percentage
 * points, and coder_qual happens to be the current fallback choice.
 *
 * fcv2 has clearly inapplicable scenarios: it needs each record's length and
 * strand direction, which only the QUAL column of an aligned SAM can provide
 * (see CoderFactory::coderSupports). In those scenarios falling back to bwt_cm
 * instead of coder_qual reduces the cost from 7.4 percentage points to 0.49
 * percentage points.
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
 * shifted more. But the statistical rule is only a prior; measurements (on a
 * synthetic 30 MB of quality values, the fine tier beats the coarse tier by
 * about 0.3 percentage points) show that which tier actually wins depends on
 * the data, so several tiers are put into the candidate set and select()
 * trial-compresses them all, picking the smallest compressed size — parameter
 * choice is driven by measurement.
 *
 * The only data-characteristic decision is the sample size: the fine tier needs
 * enough sample to converge, and with a small sample including it as a candidate
 * only adds trial-compression time and risks being wrongly chosen by a
 * not-yet-converged illusion, so with small samples only the default tier and
 * the coarse tier are kept.
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
 * fcv2 candidate tiers for the given compression level.
 *
 * Selection stays measurement-driven (each candidate is trial-compressed and
 * the smallest wins), but the candidate set is pruned by the block-volume tier
 * so trial time is not spent on tiers that cannot win at that volume. The
 * pruning follows a grid scan on real SAM QUAL (test/qual_tier_scan.cpp,
 * con_sorted.sam, 10k/25k/100k-read blocks):
 *
 *   - default: always kept (the conservative baseline).
 *   - ultra (cycleBucket=4, deltaBucket=2, prevShift=0): best at every block
 *     volume by 0.4%-1.2% over the previous best coarse; the coarse cycle bins
 *     keep every slot populated while prevShift=0 keeps the order-2 context
 *     informative.
 *   - coarse (8/4/ps2): kept as a volume-independent fallback configuration.
 *   - fine (24/12/ps0): only wins once the per-block QUAL volume is large
 *     enough (measured only at ~27 MB, hundreds of thousands of reads) for its
 *     fine buckets to stay populated; tried only at -l 8/9 (100k reads/block).
 *
 * The per-read average-quality tier (qa) and duplicate-read dedup were dropped
 * from the candidates: on the measured data they lose by 1%-1.7% (qa) and are
 * neutral (dedup), and neither is volume-dependent, so pruning them cannot be
 * recovered by a tier choice on other volumes.
 */
std::vector<Fcv2Cfg> candidateFcv2Cfgs(uint8_t compressLevel)
{
    std::vector<Fcv2Cfg> cfgs;
    cfgs.push_back(Fcv2Cfg());   /* default tier */

    Fcv2Cfg ultra;               /* coarsest cycle bins, no predecessor quantization */
    ultra.cycleBucket = 4;
    ultra.deltaBucket = 2;
    ultra.prevShift = 0;
    cfgs.push_back(ultra);

    Fcv2Cfg coarse;              /* coarse tiers, predecessor shifted right two bits */
    coarse.cycleBucket = 8;
    coarse.deltaBucket = 4;
    coarse.prevShift = 2;
    cfgs.push_back(coarse);

    if (levelVolumeTier(compressLevel) >= 3) {
        Fcv2Cfg fine;            /* fine tiers, predecessor not quantized */
        fine.cycleBucket = 24;
        fine.deltaBucket = 12;
        fine.prevShift = 0;
        cfgs.push_back(fine);
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
    uint32_t bestLen = 0;
    CoderType bestCoder = CoderType::QUAL;
    bool     anyOk = false;
    uint32_t fcv2Us = 0;
    uint32_t qualUs = 0;
    uint32_t cmUs = 0;
};

QualRoundResult runRound(const std::vector<QualSampleRecord>& records,
                         const std::vector<uint32_t>& freqByByte,
                         const std::vector<Fcv2Cfg>& cfgs,
                         uint32_t probe)
{
    QualRoundResult r;
    size_t count = recordsForBudget(records, probe);

    if (trialQual(records, count, freqByByte, r.qualLen, r.qualUs)) {
        r.qualOk = true;
        r.anyOk = true;
        r.bestLen = r.qualLen;
        r.bestCoder = CoderType::QUAL;
    }

    /* Pick the smallest among the several fcv2 tiers; this represents fcv2 when
     * compared against coder_qual / bwt_cm. */
    bool fcv2Picked = false;
    for (size_t i = 0; i < cfgs.size(); i++) {
        uint32_t len = 0, us = 0;
        if (trialFcv2(records, count, freqByByte, cfgs[i], len, us)) {
            if (!fcv2Picked || len < r.fcv2Len) {
                r.fcv2Len = len;
                r.fcv2Us = us;
                r.fcv2Cfg = cfgs[i];
                fcv2Picked = true;
            }
        }
    }
    if (fcv2Picked) {
        r.fcv2Ok = true;
        r.anyOk = true;
        if (r.fcv2Len < r.bestLen || !r.qualOk) {
            r.bestLen = r.fcv2Len;
            r.bestCoder = CoderType::FCV2;
        }
    }

    if (trialBwtCm(records, count, bwtLevelFor(probe), r.cmLen, r.cmUs)) {
        r.cmOk = true;
        r.anyOk = true;
        if (r.cmLen < r.bestLen || (!r.qualOk && !r.fcv2Ok)) {
            r.bestLen = r.cmLen;
            r.bestCoder = CoderType::BWT_CM;
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

    const std::vector<Fcv2Cfg> cfgs = candidateFcv2Cfgs(compressLevel);

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
     * that the adaptive encoders have converged. Measured on this file's quality
     * values, bwt_cm leads fcv2 by 2%-5% on small samples (64 KB-512 KB; it
     * compresses unusually well in this range), but the lead narrows as the
     * sample grows and fcv2 overtakes it around 2-4 MB. Finalizing immediately
     * at the 3% lead threshold would wrongly pick bwt_cm in the first round.
     * Hence a minimum finalization sample size is set: below it, only double,
     * never finalize.
     */
    const uint32_t MIN_SETTLE_PROBE = 1u << 20;   /* 1 MB */

    uint32_t probe = (MIN_QUAL_SAMPLE < sampleLen) ? MIN_QUAL_SAMPLE : sampleLen;
    QualRoundResult final;
    bool finalSet = false;

    while (true) {
        QualRoundResult r = runRound(records, freqByByte, cfgs, probe);
        final = r;
        finalSet = true;
        sel.rounds++;

        if (!r.anyOk) {
            break;
        }
        if (probe >= sampleLen) {
            break;
        }
        if (probe < MIN_SETTLE_PROBE) {
            probe = (probe > sampleLen / 2) ? sampleLen : (probe * 2);
            continue;
        }
        /* Only one candidate compresses successfully; adding more data gives
         * nothing to compare against. */
        {
            uint32_t runnerUp = UINT32_MAX;
            if (r.bestCoder == CoderType::QUAL) {
                if (r.fcv2Ok) runnerUp = r.fcv2Len;
                if (r.cmOk && r.cmLen < runnerUp) runnerUp = r.cmLen;
            } else if (r.bestCoder == CoderType::FCV2) {
                runnerUp = r.qualOk ? r.qualLen : UINT32_MAX;
                if (r.cmOk && r.cmLen < runnerUp) runnerUp = r.cmLen;
            } else { /* BWT_CM */
                runnerUp = r.qualOk ? r.qualLen : UINT32_MAX;
                if (r.fcv2Ok && r.fcv2Len < runnerUp) runnerUp = r.fcv2Len;
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
        probe = (probe > sampleLen / 2) ? sampleLen : (probe * 2);
    }

    sel.decidedLen = finalSet ? probe : sampleLen;
    sel.trialCount = 0;
    sel.addTrial(CoderType::QUAL, final.qualOk ? final.qualLen : 0, final.qualUs);
    sel.addTrial(CoderType::FCV2, final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us);
    sel.addTrial(CoderType::BWT_CM, final.cmOk ? final.cmLen : 0, final.cmUs);

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

    LOG_DEBUG("Qual codec trial: coder_qual=%u (%u us), fcv2=%u (%u us), bwt_cm=%u (%u us), picked=%s",
              final.qualOk ? final.qualLen : 0, final.qualUs,
              final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us,
              final.cmOk ? final.cmLen : 0, final.cmUs,
              coderTypeToMagic(sel.selectedCoder));
    return sel;
}
