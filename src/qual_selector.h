/*
 * qual_selector.h - encoder evaluation for the quality-value column
 *
 * The quality-value column cannot go through the generic evaluation path, for
 * three reasons.
 *
 * First, the candidates differ. Generic fields are evaluated with byte-stream
 * compressors such as coder_bwt_cm and coder_fc, while quality values are
 * actually compressed with coder_qual and fcv2. Previously the preprocessing
 * evaluated the quality-value column with bwt_cm and fc, but compressQuality
 * used the other two, so the evaluation result was never actually used.
 *
 * Second, the input shape differs. Generic evaluation concatenates the whole
 * column into one contiguous byte stream, losing record boundaries; but both of
 * these candidates need record-level information — fcv2 needs each record's
 * length to recover the sequencing cycle index, and coder_qual needs the
 * corresponding base sequence as context.
 *
 * Third, a compilation constraint. coder_qual.h brings in clr.h indirectly,
 * whose RangeCoder conflicts with the same-named class in fc/rangecoder.h
 * (brought in by coder_fc.h); the two cannot appear in the same translation
 * unit, so the quality-value evaluation must live in its own file.
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "preprocess_info.h"

/* One record from the sample; both evaluated candidates need this information. */
struct QualSampleRecord {
    std::string qual;   /* quality value */
    std::string seq;    /* corresponding base sequence, used by coder_qual as context */
    bool        rev;    /* strand direction, taken from bit 0x10 of FLAG, used by fcv2 to recover the cycle index */
};

class QualSelector {
public:
    /*
     * Trial-compress coder_qual, fcv2 (over several context parameter tiers)
     * and bwt_cm on the sampled records and return the one with the smaller
     * compressed size. freqByByte counts quality-value occurrences, indexed by
     * raw byte value. All three candidates need it: coder_qual builds its own
     * frequency table from it, and fcv2 builds its Huffman tree from it.
     *
     * Evaluation uses multi-round convergence: start at 64 KB and double the
     * sample each round, finalizing once the leader opens a gap of more than 3%
     * (see SETTLE_MARGIN in codec_selector.cpp). The winning fcv2 parameter tier
     * (if any) is returned via FieldCodecSelection.fcv2Params so the compression
     * side and prior training stay consistent with it; the decoder reads the
     * same parameter set back from the fcv2 bitstream header.
     *
     * If the sample is too small, SKIPPED is returned and the caller keeps its
     * default encoder. If evaluation fails (none of the candidates runs),
     * FAILED is returned and the default is used as well; an evaluation problem
     * must never make compression impossible.
     */
    static FieldCodecSelection select(const std::vector<QualSampleRecord>& records,
                                      const std::vector<uint32_t>& freqByByte);
};
