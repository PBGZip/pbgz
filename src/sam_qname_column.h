/*
 * sam_qname_column.h - the QNAME column's layouts and the analysis they need.
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

#include <cstdint>
#include <vector>

#include <json/json.h>

class RoughIOBlock;
struct coder_err_sink;

/* Cap on distinct values the "dict" layout may carry for one segment. */
inline constexpr uint32_t kMaxDictEntries = 16;

/* Caps for the fixed-alphabet layout ("hexd"): alphabet size and segment length. */
inline constexpr uint32_t kMaxHexAlphabet = 32;
inline constexpr uint32_t kMaxHexLen = 64;

/* What the ID analysis found in a block's QNAME column: separators, per-line positions and
   segment lengths. posLength == UINT32_MAX means the split is unavailable. */
struct IdSplitAnalysis {
    /* The separators, in the order the first line carries them. */
    std::vector<uint8_t> symbols;
    /* Per line, where each of them sits in the QNAME field (see analyzeQnameLine). */
    std::vector<std::vector<int32_t>> positions;
    std::vector<uint32_t> minLen;   /* shortest segment seen, per separator */
    std::vector<uint32_t> maxLen;   /* longest segment seen, per separator */
    /* Segments counted so far; UINT32_MAX once the split is unavailable. */
    uint32_t posLength = 0;
    /* Scratch: positions found in the line being analysed, cleared per call. */
    std::vector<int32_t> lineScratch;
};

/* Read one line's QNAME field; the first line establishes the separators, later lines confirm
   the same order. `length` covers the field through its tab. */
int32_t analyzeQnameFirstLine(const uint8_t* buffer, uint32_t length, IdSplitAnalysis& analysis);
int32_t analyzeQnameLine(const uint8_t* buffer, uint32_t length, IdSplitAnalysis& analysis);

/* One line's tab positions in recorded order: entry k is the offset of the (k+1)-th tab, or of
   the line's '\n'. At most CAP entries (the parse stops at the OPTION column). */
struct LineTabs {
    static const uint32_t CAP = 11;
    int64_t pos[CAP];
    uint32_t count = 0;

    uint32_t size() const { return count; }
    bool empty() const { return count == 0; }
    int64_t operator[](size_t k) const { return pos[k]; }
    /* Nothing pushes past CAP: the parse stops at the OPTION column. */
    void push_back(int64_t v) { if (count < CAP) pos[count++] = v; }
};

/* Where the QNAME column's lines are: the block being compressed, and what its parse recorded. */
struct QnameColumnInput {
    const uint8_t* buffer = nullptr;
    const std::vector<size_t>* npos = nullptr;
    int64_t headEndLine = 0;
    /* The tab positions of each content line (RoughIOBlock's parse, the actuator's contentPos). */
    const std::vector<LineTabs>* fieldTabs = nullptr;
};

/* Encode the whole column with the affix-segmentation layout: one sub-stream per separator.
   `trialLines` stops after that many content lines (0 means the whole block). Returns the
   encoded size, negative on failure. */
int32_t encodeQnameSplit(const QnameColumnInput& in, const IdSplitAnalysis& analysis,
                         RoughIOBlock* outBlock, coder_err_sink* sink, uint32_t trialLines,
                         Json::Value& fieldMeta, uint32_t& fieldSrcLen);

/* Same, with coder_qname: the column as one stream, deduplicated across lines. */
int32_t encodeQnameByQname(const QnameColumnInput& in, RoughIOBlock* outBlock, coder_err_sink* sink,
                           uint32_t trialLines, Json::Value& fieldMeta, uint32_t& fieldSrcLen);
