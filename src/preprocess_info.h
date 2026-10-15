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
#include <utility>
#include <vector>

#include "io_block.h"

/*
 * Candidate coders the pre-selector may trial-compress with. Magic strings match what
 * each coder writes into coder_io::meta, so decompression can dispatch on them.
 */
enum class CoderType : uint8_t {
    BWT_CM = 0,   /* coder_bwt_cm: BWT + context model (current default)   */
    FC,           /* coder_fc:     LZP + BWT + MTF + context range coder   */
    SIMPLE_RC,    /* coder_simple_rc: lossy, removed; kept only as an enum placeholder */
    FCV2,         /* coder_fcv2:   context-mixing coder for quality values (SAM QUAL only) */
    QUAL,         /* coder_qual:   quality-specific context model coder    */
    AFFIX_MATCH,  /* coder_affix_match: prefix/suffix matching column coder */
    ARITH,        /* coder_arith:  order-0 adaptive arithmetic byte-stream coder */
    ARITH_K,      /* coder_arith_k: htscodecs order-1 rANS byte coder (see coder_arith_k.h) */
    RANS,         /* coder_rans:   static order-0 entropy coder, fast preset only (never a trial candidate) */
    QCM,          /* coder_qcm:    order-2 context model for the quality column (see coder_qcm.h) */
    QNAME,        /* coder_qname:  the QNAME column's own coder; never a trial candidate */
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
    case CoderType::ARITH_K:   return "coder_arith_k";
    case CoderType::RANS:      return "coder_rans";
    case CoderType::QCM:       return "coder_qcm";
    case CoderType::QNAME:     return "coder_qname";
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
 * Context parameter tier for the fcv2 quality coder, picked by QualSelector's trial
 * compression. Once selected it travels with PreprocessInfo to the compression side and
 * prior training; the decoder reads it back from the stream header. Field semantics are
 * defined in Fcv2Cfg in coder_fcv2.h.
 */
struct QualFcv2Params {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  prevShift   = 1;
    bool useDedup    = false;

    bool operator==(const QualFcv2Params& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               prevShift == o.prevShift && useDedup == o.useDedup;
    }
    bool operator!=(const QualFcv2Params& o) const { return !(*this == o); }
};

/* Codec selection result for one field. */
struct FieldCodecSelection {
    FieldStatus status;
    CoderType   selectedCoder;
    uint32_t    sampleLen;
    uint32_t    bestCompLen;

    /* Trial-compression results per candidate coder: which coder, bytes produced, microseconds
       taken. The candidate set differs per field, so each field declares its list. Timings are
       indicative of order-of-magnitude differences, not small ones. */
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

    /* Sample bytes actually used to finalize the decision, and rounds run to get there.
       Evaluation starts small and doubles until the leader opens enough of a gap. */
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
 * Result of file preprocessing. Preprocessing runs once, independent of block-level
 * concurrency: the first worker to arrive runs the trial compression, others see
 * "in progress" and use their default coder without blocking. Publication is lock-free
 * (writer sets DONE with release, readers acquire-load it before consulting `fields`).
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

    /* Snapshot of the model fcv2 trained on QUAL accumulated across blocks, plus the amount
       of data it was trained on. Empty means no prior is available this run, so the caller
       falls back to cold start. */
    std::vector<uint8_t> qualPriorSnapshot;
    uint64_t qualPriorTrainingBytes;

    /* POS delta-varint file-level prior: a packed 256 x uint16 weight table
       (coder_arith::kPriorBytes bytes), produced directly during preprocessing and written
       into the file-level meta. Empty means no prior. */
    std::vector<uint8_t> posPriorSnapshot;

    /* Which coder codes the SEQ match stream ("m"), set once per file by the codec
       pre-selection (CodecSelector::selectSeqReferenceCoder). A value >= 0 is honored as-is;
       -1 means no verdict, and the coder side uses coder_bwt_cm. */
    std::atomic<int32_t> seqMatchCoder;

    /* Which form the block's N positions travel in: 0 = absolute 4-byte offsets ("npos"),
       1 = varint deltas ("nposd"), 2 = runs ("nposr"), -1 while undecided. All three are
       measured and the smallest wins; measured once on the first block that has positions. */
    std::atomic<int32_t> nposForm;

    /* Which layout writes the QNAME column: 1 = affix segmentation, 0 = coder_qname's
       cross-line deduplication, -1 while undecided. Decided from a sample of the first block
       (see CodecSelector::selectQnameLayout); -1 keeps the affix layout. */
    std::atomic<int32_t> qnameUseAffix;

    /* Whether this file is worth training and publishing a QUAL prior for; training is
       deferred until the pre-training blocks are consumed (see
       CompressEngine::finalizePretrain). */
    bool qualPriorRequested;

    /* Initializers follow the declaration order of the members. */
    PreprocessInfo() : fileType(TYPE_UNKNOW), state(PreprocessState::IDLE), sampleBytes(0), scannedBytes(0),
                       qualPriorTrainingBytes(0), seqMatchCoder(-1), nposForm(-1), qnameUseAffix(-1),
                       qualPriorRequested(false) {}

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
        seqMatchCoder.store(-1, std::memory_order_relaxed);
        nposForm.store(-1, std::memory_order_relaxed);
        qnameUseAffix.store(-1, std::memory_order_relaxed);
    }

    bool isDone() const
    {
        return state.load(std::memory_order_acquire) == PreprocessState::DONE;
    }

    /* Atomically claim the right to run preprocessing. True if this caller won the claim,
       false if another thread is already running or has finished. */
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

    /* Coder to use for a field: the selected one if preprocessing is DONE, else fallback.
       Safe from any worker without blocking. */
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

    /* The two SEQ verdicts, published by the codec pre-selection (see
       CodecSelector::selectSeqReferenceCoder) and read by every block afterwards. A caller may
       also write them before analyze to force one candidate. */
    void setSeqMatchCoder(int32_t coder)
    {
        seqMatchCoder.store(coder, std::memory_order_release);
    }

    void setNposForm(int32_t form)
    {
        nposForm.store(form, std::memory_order_release);
    }

    /* Which QNAME layout the file's names suit; see the field. */
    void setQnameUseAffix(int32_t useAffix)
    {
        qnameUseAffix.store(useAffix, std::memory_order_release);
    }
};
