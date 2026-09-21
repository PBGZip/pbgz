/*
 * sam_qname_column.cpp - the QNAME column's layouts (see sam_qname_column.h).
 */

#include "sam_qname_column.h"

#include <cstring>
#include <memory>

#include "actg.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_io.h"
#include "coder/coder_qname.h"
#include "coder/id_int_model.h"
#include "io_block.h"
#include "log/logger.h"
#include "sam_field_layout.h"
#include "utils/memory_util.h"

/*
 * The characters the QNAME column can be split on: the separators a read name is built from.
 * Every one of them that occurs in the column becomes a sub-stream of its own (see
 * analyzeQnameFirstLine), which is what lets the layouts code each segment with the alphabet it
 * actually uses.
 */
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
    std::vector<int32_t> currentLinePos;
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
/*
 * One split symbol's segments, gathered from the block: the texts as they appear (trailing
 * split symbol included, which is what the textual layout codes) and the bare payloads
 * (without it, which is what the value-domain layouts code), plus the two facts every
 * layout needs from the scan.
 */
struct IdSplitSegments {
    std::vector<std::string> texts;
    std::vector<std::string> payloads;
    bool allNumeric = true;    /* every payload is a plain decimal integer, no leading zero */
    uint32_t srcLength = 0;    /* bytes this segment occupies in the source lines */
};

/*
 * The dictionary layout: the few texts this segment takes in the block, written once, plus
 * one index per line coded as a counter (see id_int::Counter). Empty when the segment does
 * not fit the shape (a line without it, or too many distinct texts) or when the index
 * stream could not be coded.
 */
struct IdDictLayout {
    bool ok = false;
    std::vector<uint8_t> buf;
    uint32_t stepWidth = 1;
    uint32_t stepLowBits = 0;
    uint32_t values = 0;       /* indices written */
};

/*
 * The fixed-alphabet layout: every character of the payload drawn uniformly from the small
 * alphabet the segment uses, so a random hex chunk costs exactly the four bits a character
 * carries, and the split symbol comes free because the layout implies it.
 */
struct IdHexLayout {
    bool ok = false;
    std::string alphabet;
    std::vector<uint8_t> buf;
    uint32_t minLen = 0;
    uint32_t maxLen = 0;
};

/* The textual layout: the segment texts through coder_affix_match, which is what every
   other layout has to beat. */
static void encodeIdSegmentText(const std::vector<std::string>& texts, coder_io* io)
{
    std::shared_ptr<coder_affix_match> coder = std::make_shared<coder_affix_match>(io);
    for (size_t k = 0; k < texts.size(); ++k) {
        if (texts[k].empty()) {
            continue;
        }
        coder->encode_line((uint8_t*)texts[k].data(), (uint32_t)texts[k].size());
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
    const uint32_t hn = tlenPutVarint(head, (uint32_t)entries.size());
    layout.buf.assign(head, head + hn);
    for (size_t e = 0; e < entries.size(); ++e) {
        uint8_t lb[8];
        const uint32_t ln = tlenPutVarint(lb, (uint32_t)entries[e].size());
        layout.buf.insert(layout.buf.end(), lb, lb + ln);
        layout.buf.insert(layout.buf.end(), entries[e].begin(), entries[e].end());
    }

    /*
     * The indices walk a counter layout (see id_int::Counter): what they do is stay put and
     * step by one when the value changes, and coding that as a run is what keeps a segment
     * that changes twice in a block down to a few dozen bytes. An adaptive model over the
     * index values cannot do it - measured 1844 B for such a sequence against the textual
     * layout's 419 B - because it has to learn a step function one symbol at a time.
     */
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

/*
 * Build the fixed-alphabet layout, or leave it unset when the payloads do not fit the shape.
 *
 * One character is the const layout's business, and past a small alphabet a uniform code
 * costs more per character than the textual layout's model does. All-digit payloads are the
 * value-domain layouts' business: they code a run of digits as one value rather than one
 * character at a time and beat this by construction. Building the candidate for them only
 * buys the scan: measured on a 100 MB Nanopore file whose QNAME is mostly numbers, that was
 * 11 s of a 26 s run.
 */
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

/*
 * What the value-domain layouts did with an all-digit segment. They are built in one helper
 * because they share the stream they write into: each candidate is encoded behind the
 * previous one, so the winner keeps the bytes it already has and none is encoded twice.
 *
 *   numeric  the zigzag deltas as varints, left to the stream's own coder (coder_bwt_cm),
 *            which pays off when consecutive ordinals are correlated;
 *   counter  the deltas themselves through id_int::Counter, which codes the runs rather
 *            than the bytes and beats the varints when the ordinal just increments;
 *   model    the values through id_int::Model: the high bits go through an adaptive tree
 *            over the value domain, the low bits travel raw, which costs about
 *            log2(values in use) per row when the ordinals scatter;
 *   uniform  the exact range the values occupy, coded with no model at all.
 *
 * `wrote` says one of them stands in `io`; `textWroteOut` says the textual comparison that
 * runs last kept its own bytes and has already settled the block's data length.
 */
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

/*
 * What the scan of an all-digit segment produced, and nothing about how a layout encodes it:
 * the values in row order (only the rows that carry a payload), their steps, the zigzag
 * varints the first layout codes, and the two facts the later ones test.
 *
 * The varints are built during the scan because that is where the values are already in hand;
 * the caller sizes the buffer for the 10-bytes-per-row worst case (64-bit value, 7 payload
 * bits per varint byte).
 */
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
    /* UINT64_MAX until the first value arrives: this is the base of the uniform layout,
       and taking it as 0 would make that layout pay for the empty space below the range. */
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
        pos += tlenPutVarint(varints + pos, (uint32_t)zz);
        scan.deltas.push_back((int64_t)cur - (int64_t)prevVal);
        if (cur == prevVal + 1) {
            deltaOneCnt++;
        }
        deltaTotal++;
        prevVal = cur;
    }

    /*
     * Varints only pay off when consecutive ordinals are correlated (successive reads, or
     * mate pairs stored adjacently). When the file is ordered by alignment position and the
     * ordinals scatter (e.g. "ERR031968.<random id>"), the deltas are large and nearly
     * incompressible, and the decimal text wins because affix matching still finds repeated
     * digit prefixes. Require a meaningful share of pure +1 steps before switching away from
     * the text path.
     */
    scan.varintLength = pos;
    scan.minValue = minVal;
    scan.correlated = (deltaTotal > 0) && (deltaOneCnt * 10 >= deltaTotal);
    scan.valueWidth = id_int::widthFor(scan.maxValue);
    scan.valueLowBits = id_int::lowBitsFor(scan.valueWidth);
    return scan;
}

/*
 * One candidate that has just been encoded *behind* the previous one inside the segment's own
 * stream: `length` is what it came to, and the values the metadata needs travel with it, so
 * the settlement below never has to know which builder filled which.
 *
 * Encoding the candidates behind each other - rather than each in a buffer of its own - is
 * what keeps this path free of copies: the winner is promoted to the front of the stream (one
 * memmove of its own length) and the losers are simply overwritten by what follows.
 */
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

/* The zigzag delta varints through coder_bwt_cm: the layout for a correlated ordinal. Its
   stream's coder leaves data_len at the length written, so a losing counter (below) cannot
   disturb it. */
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

/* The deltas themselves through id_int::Counter, which codes the runs rather than the bytes.
   Written at `offset`, so the varints it may beat stay in front of it. */
static IdValueCandidate buildIdDeltaCounter(coder_io* io, uint32_t offset,
                                            const std::vector<int64_t>& deltas)
{
    IdValueCandidate cand;
    /* The tree carries only the steps past the small ones (see id_int::Counter): its width
       must come from those, not from the largest step overall. */
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

/* The values themselves: their high bits go through an adaptive tree over the value domain,
   the low bits travel raw, which costs about log2(values in use) per row. Written at `offset`
   (0 unless something was encoded before it). */
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

/* The decimal text through coder_affix_match: what every other layout has to beat, encoded
   last so all of them can be compared without encoding any of them twice. */
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

/* Settle the two layouts that cover a correlated ordinal: the varints, and the counter that
   codes their deltas as runs. Both are encoded - the counter behind the varints, so the winner
   keeps the bytes it already has - and the smaller wins. */
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

/*
 * Settle the three layouts that cover scattered values: the model, their exact range, and the
 * decimal text. The values are scattered, so neither the deltas nor the decimal digits can
 * reach their entropy; the model codes the values themselves, the uniform layout codes their
 * range with no model at all, and the text is encoded last so all three are compared once.
 */
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
    /* A tie keeps the earlier layout, which is what an archive written before the later one
       existed carries - hence the strict comparisons except where the original said
       "no larger". */
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
    /*
     * An all-digit segment can travel three ways, and the strongest one wins: zigzag delta
     * varints (correlated ordinals), the value itself through id_int::Model (spread over a
     * range), and the decimal text through coder_affix_match (repeated digit prefixes).
     */
    if (constMode || !allNumeric || texts.empty()) {
        return layout;
    }
    /*
     * Worst case 10 bytes per varint (64-bit value, 7 payload bits each). The encoder falls
     * back to the text path if it would not be shorter.
     */
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

/*
 * Which layout a segment ended up in, and what it costs. The order of the kinds is the
 * precedence the metadata has always used, so a switch over them reproduces it.
 */
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

/*
 * A segment holding the same text in every line of the block is stored once and costs
 * nothing per line. That is what the constant parts of a naming scheme are - instrument,
 * run, flowcell, lane - and every layout below still spends a fraction of a bit per line
 * repeating them.
 */
static bool isConstantIdSplit(const IdSplitSegments& seg, std::string& constText)
{
    if (seg.texts.empty() || seg.texts[0].empty()) {
        return false;
    }
    constText = seg.texts[0];
    for (size_t k = 1; k < seg.texts.size(); ++k) {
        /* A segment that some line does not carry at all is not constant: the other lines
           would be given text they never had. */
        if (seg.texts[k].empty() || seg.texts[k] != constText) {
            return false;
        }
    }
    return true;
}

/*
 * Gather one split symbol's segments from the block.
 *
 * Numeric sub-stream mode: when every occurrence of a segment is a plain decimal integer
 * (typical for FASTQ/SEQ read ordinals such as "SRR2769247.<n>"), the segment can be stored
 * as a binary stream of zigzag(delta) LEB128 varints instead of decimal text. Adjacent
 * ordinals are strongly correlated (delta == 1 dominates), and the decimal text wastes that
 * structure across digit boundaries. Non-numeric segments keep the original textual
 * coder_affix_match path unchanged.
 */
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

        /*
         * Keep the full segment (trailing split symbol included) for the text path, and the
         * bare payload for the numeric path / the all-digits test.
         */
        const bool hasSep = (segmentLength > 0 && segmentStart[segmentLength - 1] == splitSymbols[splitIdx]);
        const uint32_t payloadLen = hasSep ? segmentLength - 1 : segmentLength;
        std::string payload((char*)segmentStart, payloadLen);
        if (!payload.empty()) {
            if (payload.find_first_not_of("0123456789") != std::string::npos) {
                seg.allNumeric = false;
            } else if (payload.size() > 1 && payload[0] == '0') {
                /* A numeric segment with a leading zero (e.g. a zero-padded date/month "01"
                   inside a Nanopore QNAME) cannot round-trip through the zigzag-delta varint
                   stream: the decoder re-emits plain decimal text and would silently drop the
                   "0". Force the lossless textual path for such segments. */
                seg.allNumeric = false;
            }
        }
        seg.texts.push_back(std::string((char*)segmentStart, segmentLength));
        seg.payloads.push_back(payload);
        seg.srcLength += segmentLength;
    }
    return seg;
}

/*
 * Settle which layout the segment travels in and leave it in the block, then say which one
 * that was so the metadata can describe it.
 *
 * Everything that writes into the segment's own stream has written by now: the value layouts
 * wrote their payload and moved the data length themselves (they may also have already kept
 * the textual layout - see buildIdValueLayouts), and the dictionary and fixed-alphabet
 * layouts wait in buffers of their own. They are the ones that can still win here: if the
 * smaller of the two is no larger than what stands, its bytes replace it at the segment's own
 * offset, and the room is there by construction.
 */
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
            /*
             * The block's data length advances by whatever this segment produced, whichever
             * layout wrote it: the value layout writes its payload itself, and the textual
             * layout that lost the comparison has already settled the length.
             */
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

/*
 * Describe the layout that won in the stream's metadata. Absent "mode" (or "text") is the
 * legacy textual layout; "numeric" is the zigzag delta varints; "intd" the values themselves
 * through id_int::Model; "intu" their exact range; "cnt" the deltas of a counter; "dict" one
 * index per line into the few texts the block carries for this segment; "const" the single
 * text it carries for all of them; "hexd" one character at a time over a fixed alphabet.
 */
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
        /* srclen is the size the varint layout would have had: the decoder expands the coded
           deltas back into exactly that form, so everything downstream of it is unchanged. */
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
    // Similar to FastqActuator::compressIdInSplit, create multiple streams for each split symbol
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    // Process each split symbol (similar to FastqActuator)
    for (uint32_t i = 0; i < analysis.symbols.size(); ++i) {
        std::shared_ptr<coder_io> idIo = makeColumnIo(outBlock, "QNAME sub-stream", sink);
        const int64_t segBlockStart = (int64_t)outBlock->getDataLen();

        const IdSplitSegments seg = collectIdSplitSegments(i, trialLines, *in.npos,
                                                           in.buffer, analysis.positions,
                                                           analysis.symbols, in.headEndLine);
        std::string constText;
        const bool constMode = isConstantIdSplit(seg, constText);

        /*
         * The dictionary layout (see buildIdDictLayout) and the fixed-alphabet layout (see
         * buildIdHexLayout) are built in buffers of their own, so they can be tried and
         * compared against everything else; the smallest of them wins in settleIdSplitLayout.
         * Both are skipped when the segment is constant, which costs nothing per line.
         */
        const IdDictLayout dict = constMode ? IdDictLayout{} : buildIdDictLayout(seg.texts);
        const IdHexLayout hex = constMode ? IdHexLayout{} : buildIdHexLayout(seg.payloads);
        /*
         * The value-domain layouts are built in one helper (see buildIdValueLayouts); a value
         * layout writes its payload into the segment's own stream and settles the block's data
         * length itself, which is why nothing is copied here.
         */
        const IdValueLayout value = buildIdValueLayouts(seg.texts, seg.payloads, seg.allNumeric,
                                                        constMode, idIo.get(), outBlock,
                                                        segBlockStart, sink);

        const IdSplitWinner winner = settleIdSplitLayout(seg, constMode, constText, value, dict, hex,
                                                         idIo.get(), outBlock, segBlockStart);
        if (winner.kind == IdSplitWinner::Error) {
            return -1;
        }
        /* The segment's source size - except for the varint layout, whose stream is what the
           decoder expands back into, and which the counter layout reports in its place. */
        const uint32_t segSrcLen = value.numeric ? value.srcLength : seg.srcLength;
        streamMeta.append(buildIdSplitMeta(winner, i, segSrcLen, value, dict, hex, idIo.get()));
        totalSrcLength += seg.srcLength;
        totalDstLength += winner.dstLength;
    }

    // Set field metadata with streams (similar to FastqActuator)
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
        /* The QNAME field runs up to the first tab (inclusive, consistent with compressIdFieldInAll). */
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
