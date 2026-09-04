/*
 * preprocess_info.h - Data model for the file preprocessing module
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

#pragma once

#include <stdint.h>
#include <atomic>
#include <string>
#include <utility>
#include <vector>

#include "io_block.h"

/*
 * Candidate coders that the pre-selector may trial-compress with.
 * The magic strings match what each coder writes into coder_io::meta,
 * so the decompression side can dispatch on them unchanged.
 */
enum class CoderType : uint8_t {
    BWT_CM = 0,   /* coder_bwt_cm: BWT + context model (current default)   */
    FC,           /* coder_fc:     LZP + BWT + MTF + context range coder   */
    SIMPLE_RC,    /* coder_simple_rc: lossy, removed from the encode/decode path, kept only as an enum placeholder */
    FCV2,         /* coder_fcv2:   context-mixing coder for quality values, only for the SAM QUAL column */
    QUAL,         /* coder_qual:   quality-specific context model coder    */
    AFFIX_MATCH,  /* coder_affix_match: prefix/suffix matching column coder, for regular SAM fields */
    ARITH,        /* coder_arith:  order-0 adaptive arithmetic byte-stream coder */
    COUNT
};

/* Map a CoderType to the magic string stored in block meta. */
static inline const char* coderTypeToMagic(CoderType type)
{
    switch (type) {
    case CoderType::BWT_CM:    return "coder_bwt_cm";
    case CoderType::FC:        return "coder_fc";
    case CoderType::SIMPLE_RC: return "coder_simple_rc";
    case CoderType::FCV2:      return "coder_fcv2";
    case CoderType::QUAL:      return "coder_qual";
    case CoderType::AFFIX_MATCH: return "coder_affix_match";
    case CoderType::ARITH:     return "coder_arith";
    default:                   return "unknown";
    }
}

/* Selection status for a single field. */
enum class FieldStatus : uint8_t {
    SELECTED = 0,  /* a coder was chosen successfully                  */
    FAILED,        /* trial compression failed; use the default coder  */
    SKIPPED        /* field not analyzed (empty / too small / n-a)     */
};

/*
 * Context parameter tier for the fcv2 quality coder. QualSelector picks it based on
 * data characteristics (quality alphabet size, sample size, mean read length) and
 * includes it in trial compression; once selected it travels with PreprocessInfo to
 * the compression side and prior training. The decoder reads the same parameter set
 * back from the stream header and does not rely on this struct. Field semantics are
 * defined in Fcv2Cfg in coder_fcv2.h.
 */
struct QualFcv2Params {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  deltaMax    = 32;
    int  deltaBucket = 8;
    int  prevShift   = 1;
    bool useDelta    = true;
    bool useDedup    = false;
    bool useQa       = false;

    bool operator==(const QualFcv2Params& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               deltaMax == o.deltaMax && deltaBucket == o.deltaBucket &&
               prevShift == o.prevShift && useDelta == o.useDelta &&
               useDedup == o.useDedup && useQa == o.useQa;
    }
    bool operator!=(const QualFcv2Params& o) const { return !(*this == o); }
};

/* Codec selection result for one field. */
struct FieldCodecSelection {
    FieldStatus status;
    CoderType   selectedCoder;
    uint32_t    sampleLen;
    uint32_t    bestCompLen;

    /*
     * Trial-compression results for each candidate coder: which coder, how many bytes
     * it produced, and how many microseconds it took.
     *
     * This used to be two hard-coded field pairs (trialBwtCmLen / trialFcLen) because
     * the generic path happened to have only two candidates, bwt_cm and fc. The QUAL
     * column has a different candidate set, so the same field pair had to be repurposed
     * for other coders and the labels swapped at print time - once there are more than
     * two candidates, that reuse becomes self-contradictory. With an array, each field
     * declares its own candidate set and the print side just reads them off.
     *
     * The trial time exists so the selection policy can trade off compression ratio
     * against speed instead of always picking the minimum: it is common for two coders
     * to differ by under a percentage point in ratio while differing severalfold in
     * throughput. Note however that the sample is only a few hundred KB to 1-2 MB, so
     * the relative error of a single timing is not small (repeat runs of the same binary
     * fluctuate on the order of a few percent). Treat this number as indicative of
     * order-of-magnitude differences, not suitable for comparing two coders that differ
     * by a few percent.
     */
    static const uint32_t TRIAL_MAX = 3;
    CoderType   trialCoder[TRIAL_MAX] = { CoderType::BWT_CM, CoderType::BWT_CM, CoderType::BWT_CM };
    uint32_t    trialLen[TRIAL_MAX] = { 0, 0, 0 };
    uint32_t    trialUs[TRIAL_MAX] = { 0, 0, 0 };
    uint32_t    trialCount = 0;

    /* Record one candidate's trial result; results beyond TRIAL_MAX are dropped outright and do not affect selection itself. */
    void addTrial(CoderType coder, uint32_t len, uint32_t usec)
    {
        if (trialCount >= TRIAL_MAX) {
            return;
        }
        trialCoder[trialCount] = coder;
        trialLen[trialCount] = len;
        trialUs[trialCount] = usec;
        trialCount++;
    }

    /*
     * Number of sample bytes actually used to finalize the decision, and how many
     * rounds were run to get there.
     *
     * Evaluation no longer "compresses the whole sample"; instead it starts from a
     * small sample and doubles it each round, stopping as soon as the leader opens up
     * enough of a gap. Most fields are decided within the first round or two, leaving
     * the remaining budget for the genuinely hard-to-separate fields. decidedLen being
     * smaller than sampleLen means this field was finalized early.
     */
    uint32_t    decidedLen = 0;
    uint32_t    rounds = 0;

    /* Context parameter tier carried when fcv2 is selected; meaningless for other coders. */
    QualFcv2Params fcv2Params;

    FieldCodecSelection()
        : status(FieldStatus::SKIPPED),
          selectedCoder(CoderType::BWT_CM),
          sampleLen(0),
          bestCompLen(0) {}

    /* Compression ratio of the selected coder on the sample (1.0 = no gain). */
    double ratio() const
    {
        return (sampleLen == 0) ? 1.0 : (double)bestCompLen / (double)sampleLen;
    }
};

/* SAM mandatory field indices (0-based, matching the SAM spec columns). */
enum SamField : uint32_t {
    SAM_QNAME = 0,
    SAM_FLAG,
    SAM_RNAME,
    SAM_POS,
    SAM_MAPQ,
    SAM_CIGAR,
    SAM_RNEXT,
    SAM_PNEXT,
    SAM_TLEN,
    SAM_SEQ,
    SAM_QUAL,
    SAM_FIELD_COUNT   /* 11 mandatory fields */
};

/* FASTQ field indices (one record = ID / SEQ / comment / QUAL lines). */
enum FastqField : uint32_t {
    FQ_ID = 0,
    FQ_SEQ,
    FQ_QUAL,
    FQ_COMMENT,
    FQ_FIELD_COUNT    /* 4 fields */
};

/*
 * Result of file preprocessing.
 *
 * Preprocessing runs exactly once, independent of the compression pipeline's
 * block-level concurrency.  The first worker to arrive for a SAM/FASTQ block
 * performs the trial-compression; other workers see "in progress" and simply
 * use their default coder without blocking, so the pipeline never stalls.
 *
 * The publish pattern is lock-free: writers set `state` to DONE with
 * memory_order_release; readers acquire-load `state` and only consult `fields`
 * after seeing DONE, which gives them a consistent snapshot without locks.
 *
 * The struct is intentionally extensible: future preprocessing results (ID
 * split-symbol statistics, quality-score range, read-length histogram, ...)
 * can be added as new members without changing the selection contract.
 */
enum class PreprocessState : uint8_t {
    IDLE = 0,
    RUNNING,
    DONE
};

struct PreprocessInfo {
    BlockType fileType;
    std::atomic<PreprocessState> state;

    /* Sum of the per-field sample lengths, excluding tab separators, newlines, and header lines. */
    uint32_t  sampleBytes;

    /* Raw data bytes actually scanned for analysis, i.e. the size of the slice taken from the head of the data block. */
    uint32_t  scannedBytes;

    std::vector<FieldCodecSelection> fields;

    /*
     * Snapshot of the model trained by fcv2 on QUAL accumulated across blocks
     * (up to 45 MB: first block plus subsequent pre-training blocks), together with
     * the amount of data it was trained on.
     *
     * The prior travels with the preprocessing result rather than being handed to some
     * coder instance: PreprocessInfo is the only hand-off point between the pipeline and
     * the coder decision, so the coder layer need not know about data blocks or file
     * offsets, preserving the layering where coder/ never depends back on src/. An empty
     * snapshot unambiguously means no prior is available this run, letting the caller
     * unconditionally fall back to the existing cold-start behavior.
     */
    std::vector<uint8_t> qualPriorSnapshot;
    uint64_t qualPriorTrainingBytes;

    /*
     * POS delta-varint file-level prior: a packed 256 x uint16 weight table
     * (coder_arith::kPriorBytes bytes) trained by CodecSelector on the POS
     * preprocessing sample. Unlike the QUAL prior it is produced directly
     * during preprocessing (the distribution is stable, no cross-block
     * accumulation is needed), so it travels from analyze to the compression
     * side through this snapshot and is written into the file-level meta for
     * the decoding side. Empty means no prior.
     */
    std::vector<uint8_t> posPriorSnapshot;

    /*
     * analyze decides this file is worth training and publishing a QUAL prior for.
     * The decision is produced during preprocessing (the RUNNING stage) and read only
     * by the reader thread; the actual training is deferred until the reader thread has
     * consumed the pre-training blocks (cross-block accumulation); see
     * CompressEngine::finalizePretrain.
     */
    bool qualPriorRequested;

    PreprocessInfo() : fileType(TYPE_UNKNOW), state(PreprocessState::IDLE), sampleBytes(0), scannedBytes(0),
                       qualPriorTrainingBytes(0), qualPriorRequested(false) {}

    void reset(BlockType type)
    {
        fileType = type;
        state.store(PreprocessState::IDLE, std::memory_order_relaxed);
        sampleBytes = 0;
        scannedBytes = 0;
        fields.clear();
        qualPriorSnapshot.clear();
        qualPriorTrainingBytes = 0;
        qualPriorRequested = false;
        posPriorSnapshot.clear();
    }

    bool isDone() const
    {
        return state.load(std::memory_order_acquire) == PreprocessState::DONE;
    }

    /*
     * Atomically claim the right to run preprocessing.
     * Returns true if this caller won the claim (should run preprocessing),
     * false if another thread is already running or has finished.
     */
    bool tryClaim()
    {
        PreprocessState expected = PreprocessState::IDLE;
        return state.compare_exchange_strong(
            expected, PreprocessState::RUNNING,
            std::memory_order_acq_rel, std::memory_order_acquire);
    }

    void markDone()
    {
        state.store(PreprocessState::DONE, std::memory_order_release);
    }

    const FieldCodecSelection* getField(uint32_t idx) const
    {
        return (idx < fields.size()) ? &fields[idx] : nullptr;
    }

    /*
     * Return the coder to use for a field.  If preprocessing is DONE and
     * selected a coder, return it; otherwise return fallback (the actuator's
     * default).  Safe to call from any worker without blocking.
     */
    CoderType coderFor(uint32_t idx, CoderType fallback) const
    {
        if (!isDone()) {
            return fallback;
        }
        const FieldCodecSelection* sel = getField(idx);
        if (sel != nullptr && sel->status == FieldStatus::SELECTED) {
            return sel->selectedCoder;
        }
        return fallback;
    }

    bool wantQualPrior() const { return qualPriorRequested; }

    void setQualPriorRequested(bool v) { qualPriorRequested = v; }

    const std::vector<uint8_t>& qualPrior() const
    {
        return qualPriorSnapshot;
    }

    uint64_t qualPriorTrainedBytes() const
    {
        return qualPriorTrainingBytes;
    }

    void setQualPrior(std::vector<uint8_t>&& snapshot, uint64_t trainedBytes)
    {
        qualPriorSnapshot = std::move(snapshot);
        qualPriorTrainingBytes = qualPriorSnapshot.empty() ? 0 : trainedBytes;
    }

    const std::vector<uint8_t>& posPrior() const
    {
        return posPriorSnapshot;
    }

    void setPosPrior(std::vector<uint8_t>&& snapshot)
    {
        posPriorSnapshot = std::move(snapshot);
    }
};
