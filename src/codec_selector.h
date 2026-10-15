/*
 * codec_selector.h - Codec pre-selection for the file preprocessing module
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
#include <string>
#include <vector>

#include "pbgz_types.h"
#include "preprocess_info.h"
#include "qual_selector.h"

class RoughIOBlock;
class Reference;

/* A per-field sample split into lines, fed line by line to line-based coders
   (coder_affix_match); the trailing tab is included. */
struct LineSample {
    const uint8_t* data;
    uint32_t len;
};

/*
 * CodecSelector: codec pre-selection step of file preprocessing. It samples the first
 * data block, splits it into per-field byte streams, trial-compresses each with every
 * candidate coder, and records the best in a PreprocessInfo.
 *
 * Self-contained (parses the block itself) and failure-safe: fields too small or whose
 * trial fails are marked SKIPPED/FAILED and actuators fall back to defaults.
 */
class CodecSelector {
public:
    /* Cap on QUAL bytes used to train the prior; the gain past this point is
       negligible and the extra time delays worker startup. */
    static constexpr uint32_t QUAL_PRIOR_TRAIN_MAX = 45u * 1024u * 1024u; // TEMP

    /*
     * Analyze a sample of block and fill info. Returns 0 on success (including skipped
     * fields), -1 only on a fatal error. inputTotalBytes is the whole input size (0 when
     * unknown); mode selects per-field candidate lists; reference enables the SEQ trial.
     */
    static int32_t analyze(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                           uint8_t compressLevel = 5, uint8_t mode = PBGZ_MODE_ARCHIVE,
                           Reference* reference = nullptr);

    /*
     * Trial-compress one byte stream with every candidate coder and return the result.
     * Candidates are tried in registry trialPriority order, later ties losing. lines
     * supplies per-line samples for line-based coders (empty means they are skipped).
     */
    static FieldCodecSelection selectCoder(const uint8_t* data, uint32_t len,
                                           const std::vector<LineSample>* lines,
                                           const std::vector<CoderType>& candidates);

    /* Convenience overload: rebuilds the historical candidate set (coder_fc; plus
       coder_affix_match when trialAffix; plus coder_bwt_cm when allowBwtCm). */
    static FieldCodecSelection selectCoder(const uint8_t* data, uint32_t len,
                                           bool trialAffix = false,
                                           const std::vector<LineSample>* lines = nullptr,
                                           bool allowBwtCm = false);

    /* Compare coder_bwt_cm and coder_arith on the rebuilt POS delta-varint stream (not
       the raw text). posLines/chrLines supply the per-line samples; returns SKIPPED when
       the stream is too small. */
    static FieldCodecSelection selectPosDeltaCoder(const std::vector<LineSample>& posLines,
                                                   const std::vector<LineSample>& chrLines,
                                                   std::vector<uint64_t>* varintCounts = nullptr);

    /* Cross-block QUAL accumulator for prior training, held by CompressEngine on the
       read thread; collection stops at QUAL_PRIOR_TRAIN_MAX. */
    struct QualPriorAccum {
        std::vector<QualSampleRecord> records;
        std::vector<uint32_t> freqByByte;
        uint64_t collectedBytes = 0;

        bool full() const { return collectedBytes >= QUAL_PRIOR_TRAIN_MAX; }
    };

    /* Append one block's QUAL records to the accumulator (capped at QUAL_PRIOR_TRAIN_MAX). */
    static void accumulateQualPrior(RoughIOBlock* block, QualPriorAccum& acc);

    /* Train fcv2 on the accumulated samples and export a snapshot; trainedBytes reports
       the bytes trained. params must match the compression-side tier. Empty on failure. */
    static std::vector<uint8_t> trainQualPriorModel(const QualPriorAccum& acc,
                                                    const QualFcv2Params& params,
                                                    uint64_t* trainedBytes);

private:
    static int32_t analyzeSam(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                              uint8_t compressLevel, uint8_t mode, Reference* reference);

    /* FASTQ: ID / SEQ / QUAL / comment, trialled with the FASTQ candidate lists under
       the run's mode and level. QUAL uses QualSelector's dedicated record-level path. */
    static int32_t analyzeFastq(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                                uint8_t compressLevel, uint8_t mode);

    /* One record of the SEQ reference trial: its SEQ range in the sample plus views of
       FLAG/RNAME/POS/CIGAR into the block. */
    struct SeqRefRecord {
        uint32_t seqBeg = 0;
        uint32_t seqEnd = 0;
        LineSample flag;
        LineSample rname;
        LineSample pos;
        LineSample cigar;
    };

    /* Build the SEQ match payload, apply the same RLE split as the block pass, and trial
       the candidates; on success fills info.fields[SAM_SEQ] and pins
       seqMatchCoder/nposForm. False when no record can use the reference. */
    static bool selectSeqReferenceCoder(RoughIOBlock* block, Reference* reference, uint8_t level,
                                        PreprocessInfo& info);

    /* Collect SEQ reference records up to sampleBudget bases, plus the N offsets over
       the same sample. */
    static uint32_t extractSeqReferenceSamples(RoughIOBlock* block,
                                               std::vector<SeqRefRecord>& records,
                                               std::vector<uint8_t>& seq,
                                               std::vector<uint32_t>& nPositions,
                                               uint32_t sampleBudget);

    /* Choose between affix segmentation and coder_qname for the QNAME column by running
       both. Fills info.fields[SAM_QNAME] and pins qnameUseAffix; false when there is
       nothing to choose. */
    static bool selectQnameLayout(RoughIOBlock* block, PreprocessInfo& info);

    /* Extract per-field concatenated samples from a block. */
    /* Collect QUAL record samples (not a concatenated column) because the quality
       coders need record boundaries, each record's bases, and the FLAG strand bit. */
    static void extractQualSamples(RoughIOBlock* block,
                                   std::vector<QualSampleRecord>& records,
                                   std::vector<uint32_t>& freqByByte,
                                   uint32_t sampleBudget,
                                   uint64_t qualBudget = UINT64_MAX);

    static uint32_t extractSamFieldSamples(RoughIOBlock* block,
                                       std::vector<std::string>& fieldBufs,
                                       std::vector<std::vector<LineSample>>& fieldLines,
                                       uint32_t sampleBudget,
                                       uint8_t mode);

    /* RNAME+POS line views for the POS delta trial, collected by line count so the
       trial sees block-like volume; only pointers are stored. */
    static uint32_t extractPosDeltaSamples(RoughIOBlock* block,
                                           std::vector<LineSample>& posLines,
                                           std::vector<LineSample>& chrLines,
                                           uint32_t maxLines);
    /* FASTQ sampling in one pass: per-field concatenated bytes, per-line views, and the
       QUAL per-record sample with alphabet counts. */
    static uint32_t extractFastqFieldSamples(RoughIOBlock* block,
                                         std::vector<std::string>& fieldBufs,
                                         std::vector<std::vector<LineSample>>& fieldLines,
                                         std::vector<QualSampleRecord>& qualRecords,
                                         std::vector<uint32_t>& qualFreq,
                                         uint32_t sampleBudget,
                                         uint64_t qualBudget = UINT64_MAX);

    /* Choose the smallest coder_bwt_cm block level whose size fits the sample. */
    static int pickBwtLevel(uint32_t sampleLen);
};
