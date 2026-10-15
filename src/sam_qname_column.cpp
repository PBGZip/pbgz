/*
 * sam_qname_column.cpp - the QNAME column's layouts (see sam_qname_column.h).
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

#include "sam_qname_column.h"

#include <cstring>
#include <memory>

#include "coder/coder_affix_match.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_io.h"
#include "coder/coder_qname.h"
#include "coder/id_int_model.h"
#include "io_block.h"
#include "seq_stream_util.h"
#include "log/logger.h"
#include "utils/memory_util.h"

/* Separator characters a read name is built from; each occurring one becomes a sub-stream. */
static const std::string kIdSplitDefault = "/:= _.,-#\r\t\n";

/* Address of the block's output area, with the sink the coders report into. */
static std::shared_ptr<coder_io> makeColumnIo(RoughIOBlock* outBlock, const char* name,
                                              coder_err_sink* sink)
{
    return std::make_shared<coder_io>(outBlock->getCurrent(), (int32_t)outBlock->getRemain(), sink,
                                      name);
}

int32_t analyzeQnameFirstLine(const uint8_t* pBuffer, uint32_t bufLen,
                              IdSplitAnalysis& analysis) {
    if (pBuffer == nullptr || bufLen == 0) {
        return -1;
    }

    std::vector<int32_t> currentLinePos;
    for (uint32_t i = 0; i < bufLen; ++i) {
        char ch = pBuffer[i];
        if (kIdSplitDefault.find(ch) != std::string::npos) {
            analysis.symbols.push_back(ch);
            currentLinePos.push_back(static_cast<int32_t>(i));
        }
    }

    // Initialize max and min length for each separator
    for (size_t i = 0; i < analysis.symbols.size(); ++i) {
        analysis.minLen.push_back(UINT32_MAX);
        analysis.maxLen.push_back(0);
    }

    // Store first line ID analysis information to analysis.positions
    uint32_t lastPos = 0;
    for (size_t idx = 0; idx < currentLinePos.size(); ++idx) {
        uint32_t pos = currentLinePos[idx];
        uint32_t curLen = pos - lastPos + (0 == idx ? 1 : 0);   // First line no needs offset
        if (curLen < analysis.minLen[idx]) {
            analysis.minLen[idx] = curLen;
        }
        if (curLen > analysis.maxLen[idx]) {
            analysis.maxLen[idx] = curLen;
        }
        analysis.posLength++;
        lastPos = pos;
    }

    // Add the current line positions to analysis.positions
    analysis.positions.push_back(currentLinePos);
    return 0;
}

int32_t analyzeQnameLine(const uint8_t* pBuffer, uint32_t bufferLen,
                         IdSplitAnalysis& analysis) {
    std::vector<int32_t>& currentLinePos = analysis.lineScratch;   /* reused: see the note in the header */
    currentLinePos.clear();
    uint32_t lastPos = 0;
    uint32_t lastFindPos = 0;
    for (uint32_t idx = 0; idx < analysis.symbols.size(); ++idx) {
        uint8_t symbol = analysis.symbols[idx];
        uint32_t pos = lastPos;
        for (; pos < bufferLen; ++pos) {
            if (*(pBuffer + pos) == symbol) {
                uint32_t curLen = pos - lastFindPos + (0 == idx ? 1 : 0);
                if (curLen < analysis.minLen[idx]) {
                    analysis.minLen[idx] = curLen;
                }
                if (curLen > analysis.maxLen[idx]) {
                    analysis.maxLen[idx] = curLen;
                }
                currentLinePos.push_back(static_cast<int32_t>(pos));
                analysis.posLength++;
                lastPos = pos + 1;
                lastFindPos = pos;
                break;
            }
            lastPos = pos + 1;
        }

        if (pos >= bufferLen) {
            analysis.posLength = UINT32_MAX; // Mark as unavailable
            break;
        }
    }

    // Add the current line positions to analysis.positions
    analysis.positions.push_back(currentLinePos);
    return 0;
}
/* One split symbol's segment texts (with trailing separator) and bare payloads, from the block. */
struct IdSplitSegments {
    std::vector<std::string> texts;
    std::vector<std::string> payloads;
    bool allNumeric = true;    /* every payload is a plain decimal integer, no leading zero */
    uint32_t srcLength = 0;    /* bytes this segment occupies in the source lines */
};

/* Dictionary layout: the segment's few distinct texts, plus one counter-coded index per line. */
struct IdDictLayout {
    bool ok = false;
    std::vector<uint8_t> buf;
    uint32_t stepWidth = 1;
    uint32_t stepLowBits = 0;
    uint32_t values = 0;       /* indices written */
};

/* Fixed-alphabet layout: payload characters coded uniformly over the segment's alphabet. */
struct IdHexLayout {
    bool ok = false;
    std::string alphabet;
    std::vector<uint8_t> buf;
    uint32_t minLen = 0;
    uint32_t maxLen = 0;
};

/* The textual layout: the segment texts through coder_affix_match. */
static void encodeIdSegmentText(const std::vector<std::string>& texts, coder_io* io)
{
    std::shared_ptr<coder_affix_match> coder = std::make_shared<coder_affix_match>(io);
    for (size_t k = 0; k < texts.size(); ++k) {
        if (texts[k].empty()) {
            continue;
        }
        /* need2hold: the segment texts are temporary, so the coder holds its own context copy. */
        coder->encode_line((uint8_t*)texts[k].data(), (uint32_t)texts[k].size(), true);
    }
    coder->encode_flush();
}

/* Build the dictionary layout, or leave it unset when the segment does not fit the shape. */
static IdDictLayout buildIdDictLayout(const std::vector<std::string>& texts)
{
    IdDictLayout layout;
    std::vector<std::string> entries;
    for (size_t k = 0; k < texts.size(); ++k) {
        if (texts[k].empty()) {
            /* A row that does not carry this segment is not a dictionary entry. */
            return layout;
        }
        bool known = false;
        for (size_t e = 0; e < entries.size(); ++e) {
            if (entries[e] == texts[k]) {
                known = true;
                break;
            }
        }
        if (known) {
            continue;
        }
        if (entries.size() >= kMaxDictEntries) {
            return layout;
        }
        entries.push_back(texts[k]);
    }
    if (entries.size() < 2) {
        return layout;
    }

    uint8_t head[16];
    const uint32_t hn = writeVarint(head, (uint32_t)entries.size());
    layout.buf.assign(head, head + hn);
    for (size_t e = 0; e < entries.size(); ++e) {
        uint8_t lb[8];
        const uint32_t ln = writeVarint(lb, (uint32_t)entries[e].size());
        layout.buf.insert(layout.buf.end(), lb, lb + ln);
        layout.buf.insert(layout.buf.end(), entries[e].begin(), entries[e].end());
    }

    /* Indices go through a counter layout (see id_int::Counter), which codes the runs. */
    std::vector<uint32_t> indexOf;
    indexOf.reserve(texts.size());
    int64_t prevIndex = 0;
    uint64_t maxStep = 0;
    for (size_t k = 0; k < texts.size(); ++k) {
        uint32_t e = 0;
        for (; e + 1 < entries.size(); ++e) {
            if (entries[e] == texts[k]) {
                break;
            }
        }
        indexOf.push_back(e);
        const uint64_t zz = id_int::zigzag32((int64_t)e - prevIndex);
        if (zz > id_int::Counter::kSmallStep && zz > maxStep) {
            maxStep = zz;
        }
        prevIndex = (int64_t)e;
    }
    const uint32_t stepWidth = id_int::widthFor(maxStep);
    const uint32_t stepLowBits = id_int::lowBitsFor(stepWidth);
    std::vector<uint8_t> idx(texts.size() + 256);
    RangeCoder rc;
    rc.output((char*)idx.data(), (char*)idx.data() + idx.size());
    rc.StartEncode();
    id_int::Counter counter;
    counter.reset(stepWidth, stepLowBits);
    prevIndex = 0;
    for (size_t k = 0; k < indexOf.size(); ++k) {
        counter.encode(rc, (int64_t)indexOf[k] - prevIndex);
        prevIndex = (int64_t)indexOf[k];
    }
    rc.FinishEncode();
    if (rc.err != 0) {
        return layout;
    }
    layout.buf.insert(layout.buf.end(), idx.begin(), idx.begin() + rc.size_out());
    layout.ok = true;
    layout.stepWidth = stepWidth;
    layout.stepLowBits = stepLowBits;
    layout.values = (uint32_t)indexOf.size();
    return layout;
}

/* Build the fixed-alphabet layout, or leave it unset when the payloads do not fit the shape. */
static IdHexLayout buildIdHexLayout(const std::vector<std::string>& payloads)
{
    IdHexLayout layout;
    if (payloads.size() < 2) {
        return layout;
    }
    uint32_t minLen = UINT32_MAX;
    uint32_t maxLen = 0;
    std::string alpha;
    for (size_t k = 0; k < payloads.size(); ++k) {
        if (payloads[k].empty()) {
            return layout;
        }
        const uint32_t len = (uint32_t)payloads[k].size();
        if (len < minLen) {
            minLen = len;
        }
        if (len > maxLen) {
            maxLen = len;
        }
        for (size_t c = 0; c < payloads[k].size(); ++c) {
            const char ch = payloads[k][c];
            if (alpha.find(ch) == std::string::npos) {
                if (alpha.size() >= kMaxHexAlphabet) {
                    return layout;
                }
                alpha.push_back(ch);
            }
        }
    }
    if (alpha.size() < 2 || maxLen > kMaxHexLen ||
        alpha.find_first_not_of("0123456789") == std::string::npos) {
        return layout;
    }

    const uint32_t lenRange = (maxLen > minLen) ? (maxLen - minLen + 1) : 1;
    const uint32_t alphaSize = (uint32_t)alpha.size();
    std::vector<uint8_t> hb(payloads.size() * (kMaxHexLen + 8) + 64);
    RangeCoder hrc;
    hrc.output((char*)hb.data(), (char*)hb.data() + hb.size());
    hrc.StartEncode();
    for (size_t k = 0; k < payloads.size(); ++k) {
        if (lenRange > 1) {
            hrc.Encode((uint32_t)payloads[k].size() - minLen, 1, lenRange);
        }
        for (size_t c = 0; c < payloads[k].size(); ++c) {
            hrc.Encode((uint32_t)alpha.find(payloads[k][c]), 1, alphaSize);
        }
    }
    hrc.FinishEncode();
    if (hrc.err != 0 || hrc.size_out() == 0) {
        return layout;
    }
    layout.ok = true;
    layout.alphabet = alpha;
    layout.minLen = minLen;
    layout.maxLen = maxLen;
    layout.buf.assign(hb.begin(), hb.begin() + hrc.size_out());
    return layout;
}

/* What the value-domain layouts (numeric, counter, model, uniform, text) did with an all-digit
   segment; all are encoded into one shared stream. */
struct IdValueLayout {
    bool wrote = false;
    bool numeric = false;
    bool counter = false;
    bool model = false;
    bool uniform = false;
    bool textWroteOut = false;
    uint32_t srcLength = 0;    /* the varint layout's size, which is what the decoder expands to */
    uint32_t dstLength = 0;
    uint32_t width = 0;
    uint32_t lowBits = 0;
    uint32_t values = 0;
    uint32_t uniBase = 0;
    uint32_t uniTotal = 0;
};

/* What the scan of an all-digit segment produced: values, deltas, varints and test facts. */
struct IdValueScan {
    bool ok = false;              /* every payload parsed and every step 32-bit decodable */
    bool correlated = false;      /* enough pure +1 steps for the ordinal layouts to pay */
    uint32_t varintLength = 0;
    uint32_t valueWidth = 0;
    uint32_t valueLowBits = 0;
    uint64_t minValue = 0;
    uint64_t maxValue = 0;
    std::vector<uint64_t> values;
    std::vector<int64_t> deltas;

    /* Whether id_int::Model can carry these values at all (see id_int::widthFits). */
    bool modelFits() const
    {
        return ok && !values.empty() && maxValue > 0 && id_int::widthFits(valueWidth);
    }
};

/* Scan the payloads, building the varints as it goes. */
static IdValueScan scanIdValues(const std::vector<std::string>& payloads, uint8_t* varints,
                                uint32_t varintCapacity)
{
    IdValueScan scan;
    uint32_t pos = 0;
    uint64_t prevVal = 0;
    uint32_t deltaOneCnt = 0;    /* how often the ordinal just increments */
    uint32_t deltaTotal = 0;
    /* UINT64_MAX until the first value arrives (the uniform layout's base). */
    uint64_t minVal = UINT64_MAX;
    scan.ok = true;
    scan.values.reserve(payloads.size());
    scan.deltas.reserve(payloads.size());
    for (size_t k = 0; k < payloads.size() && scan.ok; ++k) {
        if (payloads[k].empty()) {
            continue;               /* keep the row count aligned */
        }
        errno = 0;
        char* endPtr = nullptr;
        unsigned long long v = std::strtoull(payloads[k].c_str(), &endPtr, 10);
        if (errno == ERANGE || endPtr != payloads[k].c_str() + payloads[k].size()) {
            scan.ok = false;
            break;
        }
        const uint64_t cur = (uint64_t)v;
        scan.values.push_back(cur);
        if (cur > scan.maxValue) {
            scan.maxValue = cur;
        }
        if (cur < minVal) {
            minVal = cur;
        }
        const uint64_t delta = (cur >= prevVal) ? (cur - prevVal) : (prevVal - cur);
        const uint64_t zz = (cur >= prevVal) ? (delta << 1) : ((delta << 1) | 1u);
        if (zz > 0xFFFFFFFFull) {   /* keep it 32-bit decodable */
            scan.ok = false;
            break;
        }
        if (pos + 10 > varintCapacity) {
            scan.ok = false;
            break;
        }
        pos += writeVarint(varints + pos, (uint32_t)zz);
        scan.deltas.push_back((int64_t)cur - (int64_t)prevVal);
        if (cur == prevVal + 1) {
            deltaOneCnt++;
        }
        deltaTotal++;
        prevVal = cur;
    }

    /* Varints pay off only when consecutive ordinals are correlated; require many pure +1 steps. */
    scan.varintLength = pos;
    scan.minValue = minVal;
    scan.correlated = (deltaTotal > 0) && (deltaOneCnt * 10 >= deltaTotal);
    scan.valueWidth = id_int::widthFor(scan.maxValue);
    scan.valueLowBits = id_int::lowBitsFor(scan.valueWidth);
    return scan;
}

/* One candidate encoded behind the previous one in the segment's stream; the winner is promoted. */
struct IdValueCandidate {
    bool ok = false;
    uint32_t length = 0;
    uint32_t width = 0;
    uint32_t lowBits = 0;
    uint32_t rows = 0;        /* rows it codes: the counter's deltas, or the values */
    uint32_t base = 0;        /* the uniform layout's range */
    uint32_t total = 0;
    Json::Value coderMeta;    /* the textual candidate's own coder meta, when it has one */
};

/* The zigzag delta varints through coder_bwt_cm. */
static IdValueCandidate buildIdOrdinalVarints(coder_io* io, const uint8_t* varints, uint32_t length)
{
    IdValueCandidate cand;
    std::shared_ptr<coder_bwt_cm> coder = std::make_shared<coder_bwt_cm>(io);
    coder->encode_line(varints, length);
    coder->encode_flush();
    if (io->err == coder_io::IO_OK && io->data_len > 0) {
        cand.ok = true;
        cand.length = (uint32_t)io->data_len;
    }
    return cand;
}

/* The deltas through id_int::Counter, written at `offset`. */
static IdValueCandidate buildIdDeltaCounter(coder_io* io, uint32_t offset,
                                            const std::vector<int64_t>& deltas)
{
    IdValueCandidate cand;
    /* The tree carries only steps past the small ones (see id_int::Counter). */
    uint64_t maxStep = 0;
    for (size_t k = 0; k < deltas.size(); ++k) {
        const uint64_t zz = id_int::zigzag32(deltas[k]);
        if (zz > id_int::Counter::kSmallStep && zz > maxStep) {
            maxStep = zz;
        }
    }
    const uint32_t width = id_int::widthFor(maxStep);
    const uint32_t lowBits = id_int::lowBitsFor(width);
    if (!id_int::widthFits(width)) {
        return cand;
    }

    RangeCoder rc;
    rc.output((char*)(io->data + offset), (char*)(io->data + io->data_capacity));
    rc.StartEncode();
    id_int::Counter counter;
    counter.reset(width, lowBits);
    for (size_t k = 0; k < deltas.size(); ++k) {
        counter.encode(rc, deltas[k]);
    }
    rc.FinishEncode();
    if (rc.err == 0 && rc.size_out() > 0) {
        cand.ok = true;
        cand.length = (uint32_t)rc.size_out();
        cand.width = width;
        cand.lowBits = lowBits;
        cand.rows = (uint32_t)deltas.size();
    }
    return cand;
}

/* The values through id_int::Model: high bits in an adaptive tree, low bits raw. */
static IdValueCandidate buildIdValueModel(coder_io* io, uint32_t offset,
                                          const std::vector<uint64_t>& values, uint32_t width,
                                          uint32_t lowBits)
{
    IdValueCandidate cand;
    id_int::Model model;
    model.reset(width, lowBits);
    RangeCoder rc;
    rc.output((char*)(io->data + offset), (char*)(io->data + io->data_capacity));
    rc.StartEncode();
    for (size_t k = 0; k < values.size(); ++k) {
        model.encode(rc, values[k]);
    }
    rc.FinishEncode();
    if (rc.err == 0 && rc.size_out() > 0) {
        cand.ok = true;
        cand.length = (uint32_t)rc.size_out();
        cand.width = width;
        cand.lowBits = lowBits;
        cand.rows = (uint32_t)values.size();
    }
    return cand;
}

/* The exact range the values occupy, coded with no model at all (see Model::resetUniform). */
static IdValueCandidate buildIdUniformRange(coder_io* io, uint32_t offset,
                                            const std::vector<uint64_t>& values,
                                            uint64_t minValue, uint64_t maxValue)
{
    IdValueCandidate cand;
    const uint64_t total = (maxValue >= minValue) ? (maxValue - minValue + 1) : 0;
    if (total == 0 || total > id_int::kMaxUniformTotal) {
        return cand;
    }
    RangeCoder rc;
    rc.output((char*)(io->data + offset), (char*)(io->data + io->data_capacity));
    rc.StartEncode();
    for (size_t k = 0; k < values.size(); ++k) {
        rc.Encode((uint32_t)(values[k] - minValue), 1, (uint32_t)total);
    }
    rc.FinishEncode();
    if (rc.err == 0 && rc.size_out() > 0) {
        cand.ok = true;
        cand.length = (uint32_t)rc.size_out();
        cand.base = (uint32_t)minValue;
        cand.total = (uint32_t)total;
        cand.rows = (uint32_t)values.size();
    }
    return cand;
}

/* The decimal text through coder_affix_match, encoded last. */
static IdValueCandidate buildIdTextSegment(coder_io* io, uint32_t offset,
                                           const std::vector<std::string>& texts,
                                           coder_err_sink* sink)
{
    IdValueCandidate cand;
    std::shared_ptr<coder_io> textIo =
        std::make_shared<coder_io>(io->data + offset, io->data_capacity - (int32_t)offset,
                                   sink, "QNAME sub-stream");
    encodeIdSegmentText(texts, textIo.get());
    if (textIo->err == coder_io::IO_OK) {
        cand.ok = true;
        cand.length = (uint32_t)textIo->data_len;
        cand.coderMeta = textIo->meta;
    }
    return cand;
}

/* Settle the correlated-ordinal layouts: varints vs. the counter; the smaller wins. */
static void settleIdOrdinalLayouts(coder_io* io, const uint8_t* varints, uint32_t varintLength,
                                   const IdValueScan& scan, IdValueLayout& layout)
{
    const IdValueCandidate numeric = buildIdOrdinalVarints(io, varints, varintLength);
    IdValueCandidate counter;
    if (numeric.ok && !scan.deltas.empty() && (uint32_t)io->data_capacity > numeric.length) {
        counter = buildIdDeltaCounter(io, numeric.length, scan.deltas);
    }

    if (counter.ok && counter.length < numeric.length) {
        memmove(io->data, io->data + numeric.length, counter.length);
        io->data_len = counter.length;
        io->meta["magic"] = "int_field";
        layout.wrote = true;
        layout.counter = true;
        layout.srcLength = varintLength;
        layout.values = counter.rows;
        layout.width = counter.width;
        layout.lowBits = counter.lowBits;
        layout.dstLength = counter.length;
    } else if (numeric.ok) {
        /* The varints stand, and the coder that wrote them left data_len at their length. */
        layout.wrote = true;
        layout.numeric = true;
        layout.srcLength = varintLength;
    }
}

/* Settle the scattered-value layouts: the model, the exact range, and the decimal text. */
static void settleIdValueLayouts(coder_io* io, const std::vector<std::string>& texts,
                                 const IdValueScan& scan, RoughIOBlock* outBlock,
                                 int64_t segBlockStart, coder_err_sink* sink,
                                 IdValueLayout& layout)
{
    const IdValueCandidate model = buildIdValueModel(io, 0, scan.values, scan.valueWidth,
                                                     scan.valueLowBits);
    IdValueCandidate uniform;
    if (model.ok && model.length < (uint32_t)io->data_capacity) {
        uniform = buildIdUniformRange(io, model.length, scan.values, scan.minValue, scan.maxValue);
    }

    if (!model.ok) {
        /* No room for the value layouts: the caller falls through to its textual path. */
        return;
    }
    if (model.length + uniform.length >= (uint32_t)io->data_capacity) {
        /* No room left for the text layout: the value layouts stand. */
        io->data_len = model.length;
        layout.wrote = true;
        layout.model = true;
        layout.width = model.width;
        layout.lowBits = model.lowBits;
        layout.dstLength = model.length;
        return;
    }

    const IdValueCandidate text = buildIdTextSegment(io, model.length + uniform.length, texts, sink);
    /* A tie keeps the earlier layout, for compatibility with older archives. */
    if (uniform.ok && uniform.length < model.length &&
        (!text.ok || uniform.length <= text.length)) {
        memmove(io->data, io->data + model.length, uniform.length);
        io->data_len = uniform.length;
        layout.wrote = true;
        layout.uniform = true;
        layout.uniBase = uniform.base;
        layout.uniTotal = uniform.total;
        layout.dstLength = uniform.length;
    } else if (text.ok && text.length < model.length) {
        memmove(io->data, io->data + model.length + uniform.length, text.length);
        io->data_len = text.length;
        io->meta = text.coderMeta;
        outBlock->setDataLen(segBlockStart + text.length);
        layout.textWroteOut = true;
    } else {
        io->data_len = model.length;
        layout.wrote = true;
        layout.model = true;
        layout.width = model.width;
        layout.lowBits = model.lowBits;
        layout.dstLength = model.length;
    }
}

static IdValueLayout buildIdValueLayouts(const std::vector<std::string>& texts,
                                         const std::vector<std::string>& payloads,
                                         bool allNumeric, bool constMode, coder_io* io,
                                         RoughIOBlock* outBlock, int64_t segBlockStart,
                                         coder_err_sink* sink)
{
    IdValueLayout layout;
    /* An all-digit segment travels as correlated varints, modelled values, or decimal text. */
    if (constMode || !allNumeric || texts.empty()) {
        return layout;
    }
    /* Worst case 10 bytes per varint (64-bit value, 7 payload bits each). */
    const uint32_t cap = (uint32_t)(texts.size() * 10);
    uint8_t* varBuf = MemoryUtil::safeAlloc<uint8_t>(cap);
    if (varBuf == nullptr) {
        return layout;
    }

    const IdValueScan scan = scanIdValues(payloads, varBuf, cap);
    if (scan.ok && scan.correlated && scan.varintLength > 0) {
        settleIdOrdinalLayouts(io, varBuf, scan.varintLength, scan, layout);
    } else if (scan.modelFits()) {
        settleIdValueLayouts(io, texts, scan, outBlock, segBlockStart, sink, layout);
    }
    MemoryUtil::safeFree(varBuf);
    return layout;
}

/* Which layout a segment ended up in, and what it costs. */
struct IdSplitWinner {
    enum Kind {
        Error = -1,   /* the output buffer could not take the segment; the caller gives up */
        Text = 0,     /* the legacy textual layout - no "mode" in the metadata */
        Const,        /* "const" */
        Numeric,      /* "numeric" */
        Counter,      /* "cnt" */
        IntModel,     /* "intd" */
        IntUniform,   /* "intu" */
        Dict,         /* "dict" */
        Hex,          /* "hexd" */
    };
    Kind kind = Text;
    uint32_t dstLength = 0;
};

/* Whether the segment holds the same text in every line of the block. */
static bool isConstantIdSplit(const IdSplitSegments& seg, std::string& constText)
{
    if (seg.texts.empty() || seg.texts[0].empty()) {
        return false;
    }
    constText = seg.texts[0];
    for (size_t k = 1; k < seg.texts.size(); ++k) {
        /* A line missing the segment makes it non-constant. */
        if (seg.texts[k].empty() || seg.texts[k] != constText) {
            return false;
        }
    }
    return true;
}

/* Gather one split symbol's segments from the block, recording texts, payloads and stats. */
static IdSplitSegments collectIdSplitSegments(uint32_t splitIdx, uint32_t trialLines,
                                              const std::vector<size_t>& npos,
                                              const uint8_t* buffer,
                                              const std::vector<std::vector<int32_t>>& splitPos,
                                              const std::vector<uint8_t>& splitSymbols,
                                              int64_t headEndLine)
{
    IdSplitSegments seg;
    const uint32_t lineNum = (uint32_t)npos.size();
    for (uint32_t lineIdx = (uint32_t)headEndLine; lineIdx < lineNum; ++lineIdx) {
        if (trialLines != 0 && lineIdx - (uint32_t)headEndLine >= trialLines) {
            break;
        }
        const uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
        const uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }

        const uint32_t contentId = lineIdx - (uint32_t)headEndLine;
        const uint8_t* segmentStart = nullptr;
        uint32_t segmentLength = 0;
        if (splitIdx == 0) {
            segmentStart = line;
            segmentLength = (uint32_t)splitPos[contentId][0] + 1;
        } else {
            const uint32_t prevPos = (uint32_t)splitPos[contentId][splitIdx - 1];
            const uint32_t currPos = (uint32_t)splitPos[contentId][splitIdx];
            segmentStart = line + prevPos + 1;
            segmentLength = currPos - prevPos;
        }

        /* Keep the full segment for the text path and the bare payload for the numeric one. */
        const bool hasSep = (segmentLength > 0 && segmentStart[segmentLength - 1] == splitSymbols[splitIdx]);
        const uint32_t payloadLen = hasSep ? segmentLength - 1 : segmentLength;
        std::string payload((char*)segmentStart, payloadLen);
        if (!payload.empty()) {
            if (payload.find_first_not_of("0123456789") != std::string::npos) {
                seg.allNumeric = false;
            } else if (payload.size() > 1 && payload[0] == '0') {
                /* A leading zero cannot round-trip through the varint stream; force the textual path. */
                seg.allNumeric = false;
            }
        }
        seg.texts.push_back(std::string((char*)segmentStart, segmentLength));
        seg.payloads.push_back(payload);
        seg.srcLength += segmentLength;
    }
    return seg;
}

/* Settle which layout the segment travels in, leave it in the block, and return the choice. */
static IdSplitWinner settleIdSplitLayout(const IdSplitSegments& seg, bool constMode,
                                         const std::string& constText, const IdValueLayout& value,
                                         const IdDictLayout& dict, const IdHexLayout& hex,
                                         coder_io* io, RoughIOBlock* outBlock,
                                         int64_t segBlockStart)
{
    IdSplitWinner winner;
    if (constMode) {
        /* Stored once, raw: no coder, and nothing per line. */
        if (outBlock->getRemain() < constText.size()) {
            LOG_ERROR("Encode id constant segment overflow: output buffer too small");
            winner.kind = IdSplitWinner::Error;
            return winner;
        }
        memcpy(outBlock->getCurrent(), constText.data(), constText.size());
        outBlock->setDataLen(outBlock->getDataLen() + constText.size());
        winner.dstLength = (uint32_t)constText.size();
    } else if (value.counter) {
        /* The counter layout wrote its payload itself and dropped the varint one it beat. */
        outBlock->setDataLen(segBlockStart + value.dstLength);
        winner.dstLength = value.dstLength;
    } else if (value.model || value.uniform) {
        /* The value layout wrote its payload itself and settled the data length already. */
        outBlock->setDataLen(segBlockStart + value.dstLength);
        winner.dstLength = value.dstLength;
    } else {
        if (!value.textWroteOut) {
            if (!value.numeric) {
                encodeIdSegmentText(seg.texts, io);
            }
            if (io->err != coder_io::IO_OK) {
                LOG_ERROR("Encode id segment overflow: output buffer too small");
                winner.kind = IdSplitWinner::Error;
                return winner;
            }
            /* Advance the block's data length by whatever this segment produced. */
            outBlock->setDataLen(outBlock->getDataLen() + io->data_len);
        }
        winner.dstLength = io->data_len;
    }

    bool dictWon = false;
    if (dict.ok && !dict.buf.empty() && dict.buf.size() <= winner.dstLength) {
        memcpy(outBlock->getBuffer() + segBlockStart, dict.buf.data(), dict.buf.size());
        outBlock->setDataLen(segBlockStart + dict.buf.size());
        winner.dstLength = (uint32_t)dict.buf.size();
        dictWon = true;
    }
    bool hexWon = false;
    if (hex.ok && !hex.buf.empty() && hex.buf.size() <= winner.dstLength) {
        memcpy(outBlock->getBuffer() + segBlockStart, hex.buf.data(), hex.buf.size());
        outBlock->setDataLen(segBlockStart + hex.buf.size());
        winner.dstLength = (uint32_t)hex.buf.size();
        hexWon = true;
    }

    if (hexWon) {
        winner.kind = IdSplitWinner::Hex;
    } else if (dictWon) {
        winner.kind = IdSplitWinner::Dict;
    } else if (value.numeric) {
        winner.kind = IdSplitWinner::Numeric;
    } else if (constMode) {
        winner.kind = IdSplitWinner::Const;
    } else if (value.counter) {
        winner.kind = IdSplitWinner::Counter;
    } else if (value.uniform) {
        winner.kind = IdSplitWinner::IntUniform;
    } else if (value.model) {
        winner.kind = IdSplitWinner::IntModel;
    } else {
        winner.kind = IdSplitWinner::Text;
    }
    return winner;
}

/* Describe the winning layout in the stream's metadata. */
static Json::Value buildIdSplitMeta(const IdSplitWinner& winner, uint32_t splitIdx,
                                    uint32_t srcLength, const IdValueLayout& value,
                                    const IdDictLayout& dict, const IdHexLayout& hex,
                                    const coder_io* io)
{
    Json::Value meta;
    meta["srclen"] = srcLength;
    meta["dstlen"] = winner.dstLength;
    meta["coder"] = io->meta;
    meta["splitidx"] = splitIdx;
    switch (winner.kind) {
    case IdSplitWinner::Hex:
        meta["mode"] = "hexd";
        meta["hexa"] = hex.alphabet;
        meta["hexn"] = (Json::Value::UInt)hex.minLen;
        meta["hexm"] = (Json::Value::UInt)hex.maxLen;
        meta["coder"]["magic"] = "int_field";
        break;
    case IdSplitWinner::Dict:
        meta["mode"] = "dict";
        meta["intw"] = (Json::Value::UInt)dict.stepWidth;
        meta["intk"] = (Json::Value::UInt)dict.stepLowBits;
        meta["numv"] = (Json::Value::UInt)dict.values;
        meta["coder"]["magic"] = "int_field";
        break;
    case IdSplitWinner::Numeric:
        meta["mode"] = "numeric";
        break;
    case IdSplitWinner::Const:
        meta["mode"] = "const";
        meta["coder"]["magic"] = "const";
        break;
    case IdSplitWinner::Counter:
        meta["mode"] = "cnt";
        /* srclen is the size the varint layout would have had. */
        meta["srclen"] = value.srcLength;
        meta["intw"] = (Json::Value::UInt)value.width;
        meta["intk"] = (Json::Value::UInt)value.lowBits;
        meta["numv"] = (Json::Value::UInt)value.values;
        meta["coder"]["magic"] = "int_field";
        break;
    case IdSplitWinner::IntUniform:
        meta["mode"] = "intu";
        meta["intb"] = (Json::Value::UInt)value.uniBase;
        meta["intn"] = (Json::Value::UInt)value.uniTotal;
        meta["coder"]["magic"] = "int_field";
        break;
    case IdSplitWinner::IntModel:
        meta["mode"] = "intd";
        meta["intw"] = (Json::Value::UInt)value.width;
        meta["intk"] = (Json::Value::UInt)value.lowBits;
        meta["coder"]["magic"] = "int_field";
        break;
    case IdSplitWinner::Text:
    default:
        break;   /* the legacy textual layout carries no "mode" */
    }
    return meta;
}


int32_t encodeQnameSplit(const QnameColumnInput& in, const IdSplitAnalysis& analysis,
                         RoughIOBlock* outBlock, coder_err_sink* sink, uint32_t trialLines,
                         Json::Value& fieldMeta, uint32_t& fieldSrcLen)
{
    // One sub-stream per split symbol
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    // Process each split symbol
    for (uint32_t i = 0; i < analysis.symbols.size(); ++i) {
        std::shared_ptr<coder_io> idIo = makeColumnIo(outBlock, "QNAME sub-stream", sink);
        const int64_t segBlockStart = (int64_t)outBlock->getDataLen();

        const IdSplitSegments seg = collectIdSplitSegments(i, trialLines, *in.npos,
                                                           in.buffer, analysis.positions,
                                                           analysis.symbols, in.headEndLine);
        std::string constText;
        const bool constMode = isConstantIdSplit(seg, constText);

        /* The dict and hex layouts are built in their own buffers; the smallest wins. */
        const IdDictLayout dict = constMode ? IdDictLayout{} : buildIdDictLayout(seg.texts);
        const IdHexLayout hex = constMode ? IdHexLayout{} : buildIdHexLayout(seg.payloads);
        /* Value layouts write into the segment's stream and settle the data length themselves. */
        const IdValueLayout value = buildIdValueLayouts(seg.texts, seg.payloads, seg.allNumeric,
                                                        constMode, idIo.get(), outBlock,
                                                        segBlockStart, sink);

        const IdSplitWinner winner = settleIdSplitLayout(seg, constMode, constText, value, dict, hex,
                                                         idIo.get(), outBlock, segBlockStart);
        if (winner.kind == IdSplitWinner::Error) {
            return -1;
        }
        /* The segment's source size, except for the varint layout which reports its own. */
        const uint32_t segSrcLen = value.numeric ? value.srcLength : seg.srcLength;
        streamMeta.append(buildIdSplitMeta(winner, i, segSrcLen, value, dict, hex, idIo.get()));
        totalSrcLength += seg.srcLength;
        totalDstLength += winner.dstLength;
    }

    // Set field metadata
    fieldMeta["totalsrclen"] = totalSrcLength;
    fieldMeta["totaldstlen"] = totalDstLength;
    fieldMeta["splitsym"] = std::string((char*)analysis.symbols.data(), analysis.symbols.size());
    fieldMeta["streams"] = streamMeta;
    fieldMeta["field"] = 0;

    fieldSrcLen = totalSrcLength;

    LOG_DEBUG("SAM ID compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            totalSrcLength, totalDstLength, (double)(totalDstLength * 100)/(double)totalSrcLength);

    return totalDstLength;
}

int32_t encodeQnameByQname(const QnameColumnInput& in, RoughIOBlock* outBlock, coder_err_sink* sink,
                           uint32_t trialLines, Json::Value& fieldMeta, uint32_t& fieldSrcLen)
{
    const std::vector<size_t>& npos = *in.npos;
    const uint32_t lineNum = (uint32_t)npos.size();
    const uint8_t* buffer = in.buffer;

    std::shared_ptr<coder_io> fieldIo = makeColumnIo(outBlock, "QNAME", sink);
    std::shared_ptr<coder_qname> qnameCoder = std::make_shared<coder_qname>(fieldIo.get());

    fieldSrcLen = 0;
    for (uint32_t lineIdx = (uint32_t)in.headEndLine; lineIdx < lineNum; ++lineIdx) {
        if (trialLines != 0 && lineIdx - (uint32_t)in.headEndLine >= trialLines) {
            break;
        }
        uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
        if (buffer[lineStart] == '@') {
            continue;
        }
        uint32_t contentId = lineIdx - (uint32_t)in.headEndLine;
        const uint8_t* idStart = buffer + lineStart;
        /* The QNAME field runs up to the first tab, inclusive. */
        uint32_t idLength = (uint32_t)(*in.fieldTabs)[contentId][0] + 1;
        qnameCoder->encode_line(idStart, idLength);
        fieldSrcLen += idLength;
    }

    qnameCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode QNAME (coder_qname) overflow: output buffer too small");
        return -1;
    }
    outBlock->setDataLen(outBlock->getDataLen() + fieldIo->data_len);

    Json::Value streamMeta;
    Json::Value tmpMeta;
    tmpMeta["srclen"] = fieldSrcLen;
    tmpMeta["dstlen"] = fieldIo->data_len;
    tmpMeta["coder"] = fieldIo->meta;
    tmpMeta["splitidx"] = 0;
    streamMeta.append(tmpMeta);

    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["splitsym"] = "\t";
    fieldMeta["field"] = 0;

    LOG_INFO("SAM QNAME (coder_qname) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}
