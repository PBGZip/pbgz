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

#include "preprocess_info.h"
#include "qual_selector.h"

class RoughIOBlock;

/*
 * A per-field sample split into lines. coder_affix_match is a line-based
 * encoder whose prefix/suffix matching gains come from adjacent lines, so it
 * must be fed line by line during trial compression (including the trailing
 * tab); joining the whole column into one byte stream loses line boundaries,
 * and the measured numbers would have nothing in common with real compression.
 */
struct LineSample {
    const uint8_t* data;
    uint32_t len;
};

/*
 * CodecSelector
 *
 * Codec pre-selection step of the file preprocessing module. It scans a sample
 * of the first data block, splits it into per-field byte streams (SAM columns
 * or FASTQ lines), trial-compresses each stream with every candidate coder and
 * records the best one into a PreprocessInfo.
 *
 * The selector is self-contained: it parses the block itself and does not rely
 * on any actuator state, which keeps it independently unit-testable.
 *
 * Failure policy: any field that is too small to give a reliable comparison, or
 * whose trial compression fails, is marked SKIPPED / FAILED. Actuators then fall
 * back to their built-in default coder, so a preprocessing failure never breaks
 * compression.
 */
class CodecSelector {
public:
    /* Prior training uses at most 45 MB of QUAL. The measured gain at this
     * point is about 1.84 percentage points; going to 90 MB adds only about
     * 0.10 percentage points, yet the extra time falls in the serial phase
     * before parallel compression kicks in, directly delaying worker-thread
     * startup, so training must be cut off before the gain curve flattens. */
    static constexpr uint32_t QUAL_PRIOR_TRAIN_MAX = 45u * 1024u * 1024u; // TEMP

    /*
     * Analyze a sample of the given block and fill info.
     * Returns 0 on success (including the "some fields skipped" case),
     * -1 only on a fatal error such as a null block.
     *
     * inputTotalBytes is the byte count of the entire input; 0 means unknown
     * (piped input). It is only used for decisions that pay off when amortized
     * over the whole file — currently whether the QUAL prior should be trained
     * and written. It is a parameter rather than part of PreprocessInfo: the
     * latter is the product of the decision and should not double as an input.
     */
    static int32_t analyze(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                           uint8_t compressLevel = 5);

    /*
     * Trial-compress one byte stream with every candidate coder and return the
     * selection result. Exposed for unit testing.
     *
     * trialAffix controls whether coder_affix_match is included as a candidate.
     * affix is a column-wise prefix/suffix matching encoder that only helps
     * certain SAM regular fields (FLAG/POS/MAPQ/CIGAR/PNEXT/TLEN); for other
     * fields (e.g. SEQ, QNAME) it would only waste time in the comparison and
     * might even be wrongly selected.
     *
     * affix must be fed line by line (see LineSample); lines supplies the
     * per-line sample. When it is empty the trial degrades to one whole-stream
     * encode_line call, whose result is not trustworthy, so affix is excluded
     * from the comparison in that case.
     */
    static FieldCodecSelection selectCoder(const uint8_t* data, uint32_t len,
                                           bool trialAffix = false,
                                           const std::vector<LineSample>* lines = nullptr);

    /*
     * Trial-compress the POS delta-varint stream (the actual byte stream fed to
     * the POS coder, see compressPosFieldDelta) with coder_bwt_cm and
     * coder_arith, and return the better one. The raw column text must not be
     * used for this comparison: the varint delta stream has a very different
     * distribution from the decimal POS text.
     *
     * posLines / chrLines are the per-line POS and RNAME samples (LineSample
     * views including the trailing tab). The delta baseline resets at every
     * RNAME change, mirroring compressPosFieldDelta's chromosome-switch logic.
     * Returns SKIPPED when the rebuilt varint stream is too small to compare
     * reliably; the actuator then falls back to bwt_cm.
     */
    static FieldCodecSelection selectPosDeltaCoder(const std::vector<LineSample>& posLines,
                                                   const std::vector<LineSample>& chrLines,
                                                   std::vector<uint64_t>* varintCounts = nullptr);

    /*
     * Cross-block training sample accumulator for the QUAL prior. Held by
     * CompressEngine on the read thread: QUAL of the first N blocks is
     * accumulated record by record (each block is scanned to its end, but
     * collection stops once QUAL_PRIOR_TRAIN_MAX is reached); after the target
     * number of blocks is read or the accumulator is full, trainQualPriorModel
     * trains once and exports a snapshot.
     *
     * Why across blocks: the first block usually holds only ~9 MB of QUAL, while
     * the prior's gain roughly converges at 45 MB (see QUAL_PRIOR_TRAIN_MAX).
     * Accumulating across multiple blocks lets small blocks approach that
     * training volume too.
     */
    struct QualPriorAccum {
        std::vector<QualSampleRecord> records;
        std::vector<uint32_t> freqByByte;
        uint64_t collectedBytes = 0;

        bool full() const { return collectedBytes >= QUAL_PRIOR_TRAIN_MAX; }
    };

    /* Append one block's QUAL records to the accumulator (cap
     * QUAL_PRIOR_TRAIN_MAX; anything beyond it is dropped). */
    static void accumulateQualPrior(RoughIOBlock* block, QualPriorAccum& acc);

    /* Train fcv2 on the accumulated samples and export a model snapshot;
     * trainedBytes reports the actual number of bytes trained. params is the
     * context parameter tier selected by QualSelector (meaningful only when
     * fcv2 wins) and must match the tier used for real compression, otherwise
     * loading the prior snapshot fails wholesale due to a layout mismatch.
     * Returns an empty vector on failure or when the samples are empty. */
    static std::vector<uint8_t> trainQualPriorModel(const QualPriorAccum& acc,
                                                    const QualFcv2Params& params,
                                                    uint64_t* trainedBytes);

private:
    static int32_t analyzeSam(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info,
                              uint8_t compressLevel);
    static int32_t analyzeFastq(RoughIOBlock* block, PreprocessInfo& info);

    /* Extract per-field concatenated samples from a block. */
    /*
     * Collect samples for the quality-value column only.
     *
     * Generic sampling concatenates the whole column into contiguous bytes,
     * losing record boundaries; but both quality-value candidates need
     * record-level information — fcv2 needs each record's length to recover the
     * sequencing cycle index, and coder_qual needs the corresponding base
     * sequence as context, while the strand direction is taken from bit 0x10 of
     * FLAG. So records are collected here individually, without concatenation.
     */
    static void extractQualSamples(RoughIOBlock* block,
                                   std::vector<QualSampleRecord>& records,
                                   std::vector<uint32_t>& freqByByte,
                                   uint32_t sampleBudget,
                                   uint64_t qualBudget = UINT64_MAX);

    static uint32_t extractSamFieldSamples(RoughIOBlock* block,
                                       std::vector<std::string>& fieldBufs,
                                       std::vector<std::vector<LineSample>>& fieldLines,
                                       uint32_t sampleBudget);

    /*
     * RNAME+POS line views collected by line count instead of the shared byte
     * budget, for the POS delta-varint trial. Only pointers into the block
     * buffer are stored, so sampling a whole block is almost free, and the
     * trial then measures the coders at the volume real compression feeds them
     * (a small sample overstates the cost of a coder that pays a model cold
     * start once per block).
     */
    static uint32_t extractPosDeltaSamples(RoughIOBlock* block,
                                           std::vector<LineSample>& posLines,
                                           std::vector<LineSample>& chrLines,
                                           uint32_t maxLines);
    static uint32_t extractFastqFieldSamples(RoughIOBlock* block,
                                         std::vector<std::string>& fieldBufs,
                                         uint32_t sampleBudget);

    /* Choose the smallest coder_bwt_cm block level whose size fits the sample. */
    static int pickBwtLevel(uint32_t sampleLen);
};
