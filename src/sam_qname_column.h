/*
 * sam_qname_column.h - the QNAME column's layouts: the analysis they need, and the two encoders
 * the file-level choice is between.
 *
 * The column is written one of two ways, and which one is a property of the file's naming scheme
 * rather than of a block: affix segmentation (coder_affix_match over the split segments) wins on a
 * single FASTQ, whose names share a constant prefix and carry locally ordered fragment numbers,
 * while coder_qname's cross-line deduplication wins on a concatenated file, whose alternating
 * prefixes leave nothing for the adjacent-prefix model to hold on to.
 *
 * So it is a file-level verdict like the SEQ column's: the codec pre-selection runs both encoders
 * on a sample of the first block and publishes which won (see PreprocessInfo::qnameUseAffix), and
 * the block pass only encodes the way it says. This module exists so that the trial and the block
 * pass run the *same* encoders - a verdict measured on a different encoder than the one used is
 * worthless - and so that neither of them has to carry the layouts itself.
 *
 * Everything the encoders need is passed in: they read the block's bytes and its parsed tab
 * positions, and write into an output block the caller owns (the trial hands them a scratch one).
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <json/json.h>

class RoughIOBlock;
struct coder_err_sink;

/*
 * A QNAME segment that takes only a handful of distinct values inside a block - a tile number, a
 * lane - is a dictionary that block can carry once. The "dict" layout (see the module source) is
 * built on this cap, and the reader checks it against the count it finds.
 */
inline constexpr uint32_t kMaxDictEntries = 16;

/*
 * Caps for the fixed-alphabet layout ("hexd"): a UUID chunk draws from 16 characters, a base32
 * barcode from 32; past that a uniform character code costs more than the textual layout's model,
 * and the layout is measured against it anyway.
 */
inline constexpr uint32_t kMaxHexAlphabet = 32;
inline constexpr uint32_t kMaxHexLen = 64;

/*
 * What the ID analysis found in a block's QNAME column.
 *
 * The column is split at every separator character a read name is built from (see
 * analyzeQnameFirstLine); each occurrence becomes a boundary between two segments, and the
 * segments are what the layouts code. A segment that some line does not carry at all, or a line
 * whose separators do not line up with the first line's, makes the split unavailable - which is
 * what `posLength == UINT32_MAX` says, and the column is then written whole.
 *
 * The reader needs the same facts (it rebuilds the column segment by segment), so the analysis is
 * a value both sides hold rather than state inside an encoder.
 */
struct IdSplitAnalysis {
    /* The separators, in the order the first line carries them. */
    std::vector<uint8_t> symbols;
    /* Per line, where each of them sits in the QNAME field (see analyzeQnameLine). */
    std::vector<std::vector<int32_t>> positions;
    std::vector<uint32_t> minLen;   /* shortest segment seen, per separator */
    std::vector<uint32_t> maxLen;   /* longest segment seen, per separator */
    /* Segments counted so far; once the split turns out to be unavailable (a line does not carry
       the separators in the same order) every later line leaves it at UINT32_MAX. */
    uint32_t posLength = 0;
};

/*
 * Read one line's QNAME field: the first content line establishes which separators the column
 * uses and where they sit in it, and every line after that only has to confirm it carries them in
 * the same order (see analyzeQnameLine). `length` covers the field up to and including its tab.
 */
int32_t analyzeQnameFirstLine(const uint8_t* buffer, uint32_t length, IdSplitAnalysis& analysis);
int32_t analyzeQnameLine(const uint8_t* buffer, uint32_t length, IdSplitAnalysis& analysis);

/* Where the QNAME column's lines are: the block being compressed, and what its parse recorded. */
struct QnameColumnInput {
    const uint8_t* buffer = nullptr;
    const std::vector<size_t>* npos = nullptr;
    int64_t headEndLine = 0;
    /* The tab positions of each content line (RoughIOBlock's parse, the actuator's contentPos). */
    const std::vector<std::vector<int64_t>>* fieldTabs = nullptr;
};

/*
 * Encode the whole column with the affix-segmentation layout: one sub-stream per separator, each
 * coded with whichever of its layouts came out smallest (see the module source for the catalogue).
 *
 * `trialLines` stops after that many content lines, which is how the trial measures the encoder on
 * a sample; 0 means the whole block. On success `fieldMeta` describes the streams and
 * `fieldSrcLen` is the column's source length. Returns the encoded size, negative on failure.
 */
int32_t encodeQnameSplit(const QnameColumnInput& in, const IdSplitAnalysis& analysis,
                         RoughIOBlock* outBlock, coder_err_sink* sink, uint32_t trialLines,
                         Json::Value& fieldMeta, uint32_t& fieldSrcLen);

/* Same, with coder_qname: the column as one stream, deduplicated across lines. */
int32_t encodeQnameByQname(const QnameColumnInput& in, RoughIOBlock* outBlock, coder_err_sink* sink,
                           uint32_t trialLines, Json::Value& fieldMeta, uint32_t& fieldSrcLen);
