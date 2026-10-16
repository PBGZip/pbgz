/*
 * sam_seq_payload.cpp - the reference-coded SEQ payload (see sam_seq_payload.h).
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

#include "sam_seq_payload.h"

#include <atomic>
#include "seq_stream_util.h"   /* appendVarint, used by the payload writers */
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>

#include "actg.h"
#include "log/logger.h"
#include "sam_info.h"

/* The one CIGAR parse: reference span, read span, and operation list. Lower-case letters are
   counted in the spans as they always have been; `ops` is filled only when requested. */
static void parseCigarImpl(const uint8_t* cigar, uint32_t cigarLength, uint32_t* refConsumedOut,
                           uint32_t* seqConsumedOut, std::vector<CigarOp>* ops)
{
    if (ops != nullptr) {
        ops->clear();
    }
    if (refConsumedOut != nullptr) {
        *refConsumedOut = 0;
    }
    if (seqConsumedOut != nullptr) {
        *seqConsumedOut = 0;
    }
    if (cigar == nullptr || cigarLength == 0) {
        return;
    }
    uint32_t refConsumed = 0;
    uint32_t seqConsumed = 0;
    uint32_t currentNumber = 0;
    for (uint32_t i = 0; i < cigarLength; ++i) {
        const char ch = (char)cigar[i];
        if (ch >= '0' && ch <= '9') {
            currentNumber = currentNumber * 10 + (uint32_t)(ch - '0');
            continue;
        }
        if (currentNumber > 0) {
            if (ops != nullptr &&
                (ch == 'M' || ch == 'I' || ch == 'D' || ch == 'N' || ch == 'S' || ch == 'H' ||
                 ch == 'P' || ch == '=' || ch == 'X')) {
                ops->push_back(CigarOp{ch, currentNumber});
            }
            switch (ch) {
                case 'M': case '=': case 'X': case 'm':   /* both spans claim these */
                    refConsumed += currentNumber;
                    seqConsumed += currentNumber;
                    break;
                case 'D': case 'N': case 'd': case 'n':   /* reference only */
                    refConsumed += currentNumber;
                    break;
                case 'I': case 'S': case 'i': case 's': case 'x':   /* read only */
                    seqConsumed += currentNumber;
                    break;
                default:
                    break;
            }
            currentNumber = 0;
        }
    }
    if (refConsumedOut != nullptr) {
        *refConsumedOut = refConsumed;
    }
    if (seqConsumedOut != nullptr) {
        *seqConsumedOut = seqConsumed;
    }
}

uint32_t parseCigarOps(const uint8_t* cigar, uint32_t cigarLength, std::vector<CigarOp>& ops)
{
    uint32_t refConsumed = 0;
    parseCigarImpl(cigar, cigarLength, &refConsumed, nullptr, &ops);
    return refConsumed;
}

uint32_t cigarSeqConsumed(const uint8_t* cigar, uint32_t cigarLength)
{
    uint32_t seqConsumed = 0;
    parseCigarImpl(cigar, cigarLength, nullptr, &seqConsumed, nullptr);
    return seqConsumed;
}

uint32_t cigarRefConsumed(const uint8_t* cigar, uint32_t cigarLength)
{
    uint32_t refConsumed = 0;
    parseCigarImpl(cigar, cigarLength, &refConsumed, nullptr, nullptr);
    return refConsumed;
}

/* splitSeqMatchStream lives in seq_stream_util.h. */

bool seqRecordUsesReference(uint16_t chrId, uint16_t flag, uint64_t startPos, uint32_t refConsumed,
                            Reference* reference, int64_t& refPos)
{
    refPos = 0;
    if (reference == nullptr) {
        /* No reference to address: every record stores its own bases as they stand. */
        return false;
    }
    if (chrId == SEQ_CHR_ID_NONE || chrId == SEQ_CHR_ID_UNKNOWN || (flag & SEQ_FLAG_UNMAPPED)) {
        return false;
    }
    const int64_t chrStartPos = SamInfo::getInstance().getPositionByIndex(chrId);
    if (chrStartPos == -1) {
        return false;
    }
    refPos = chrStartPos + startPos - 1;   /* SAM is 1-based */
    const uint64_t needSquash = (uint64_t)(refConsumed >> 2) + !!(refConsumed & 0x3) + 1;
    if (refPos < 0 ||
        ((uint64_t)(refPos >> 2)) + needSquash > (uint64_t)reference->getSquashLength()) {
        return false;
    }
    return true;
}

bool referenceLinesUpWithHeader(const Reference* reference)
{
    if (reference == nullptr) {
        return false;
    }
    const std::vector<Reference::SequenceSpan>& seqs = reference->getSequences();
    if (seqs.empty()) {
        /* An NI index carries no names; an unverifiable reference is accepted. */
        return true;
    }

    /* Where each name really starts, to compare with the @SQ running-sum offsets. */
    std::unordered_map<std::string, uint64_t> offsetByName;
    offsetByName.reserve(seqs.size() * 2);
    for (const Reference::SequenceSpan& s : seqs) {
        offsetByName.emplace(s.name, s.offset);
    }

    const std::vector<ChromosomeInfo>& chrs = SamInfo::getInstance().getAllChromosomeInfo();
    for (const ChromosomeInfo& c : chrs) {
        const auto it = offsetByName.find(c.name);
        if (it == offsetByName.end() || it->second != c.position) {
            return false;
        }
    }
    return true;
}

void warnReferenceNotUsableForSeq()
{
    static std::atomic<bool> warned{false};
    if (warned.exchange(true, std::memory_order_relaxed)) {
        return;
    }
    /* Nothing else reports this: the reference loads and verifies, the SEQ column just stops
       improving. */
    LOG_WARNING("The reference does not line up with the file's @SQ list: its sequences are not "
                "the ones the header names, in the header's order, so a read's position would "
                "address the wrong bases. SEQ is compressed without it; supply a reference of "
                "the same assembly as the alignment.");
    fprintf(stderr, "warning: reference does not line up with the @SQ list; "
                    "SEQ falls back to compression without a reference\n");
}

/* Walk one record's ops and produce its coded SEQ bases: each base as its 2-bit value, XORed
   with the reference where the op consumes it (M/=/X). `stretch` supplies the reference bases.
   False when the ops do not describe the read exactly. */
template <class StretchFn>
static bool seqCodedBasesWalk(const std::vector<CigarOp>& ops, const uint8_t* seq,
                              uint32_t seqLength, int64_t refPos, StretchFn stretch,
                              uint8_t* ref2bitScratch, uint8_t* out)
{
    const SeqWalkEnd end = seqWalkCigar(ops, seqLength, refPos,
        [&](uint32_t readPos, int64_t pos, uint32_t len) {
            stretch(ref2bitScratch, len, pos);
            for (uint32_t i = 0; i < len; ++i) {
                const uint8_t read2 = (uint8_t)((seq[readPos + i] >> 1) & 0x3);
                out[readPos + i] = (uint8_t)(read2 ^ ref2bitScratch[i]);
            }
            return SeqWalkAction::Continue;
        },
        [&](uint32_t readPos, uint32_t len) {
            for (uint32_t i = 0; i < len; ++i) {
                out[readPos + i] = (uint8_t)((seq[readPos + i] >> 1) & 0x3);
            }
            return SeqWalkAction::Continue;
        },
        [](uint32_t) { return SeqWalkAction::Continue; });
    return end.complete;
}

bool buildSeqReferenceCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                                 uint32_t seqLength, int64_t refPos, Reference* reference,
                                 uint8_t* ref2bitScratch, uint8_t* out)
{
    return seqCodedBasesWalk(ops, seq, seqLength, refPos,
                             [reference](uint8_t* dst, uint32_t len, int64_t pos) {
                                 reference->getStretch2Bits1Char(dst, len, pos);
                             },
                             ref2bitScratch, out);
}

SeqRecordPayload buildSeqRecordPayload(uint16_t chrId, uint16_t flag, uint64_t startPos,
                                       uint32_t refConsumed, const std::vector<CigarOp>& ops,
                                       const uint8_t* seq, uint32_t seqLength,
                                       Reference* reference, uint8_t* coded, uint8_t* ref2bit)
{
    SeqRecordPayload out;
    out.bytes = coded;
    out.length = seqLength;

    int64_t refPos = 0;
    const bool useReference =
        seqRecordUsesReference(chrId, flag, startPos, refConsumed, reference, refPos);
    if (useReference && !ops.empty()) {
        if (buildSeqReferenceCodedBases(ops, seq, seqLength, refPos, reference, ref2bit, coded)) {
            out.usedReference = true;
            out.refPos = refPos;
        } else {
            /* CIGAR and SEQ disagree: store the bases as their own 2-bit codes. */
            actgEncode(seq, coded, seqLength);
        }
    } else {
        if (seqLength > 0) {
            std::memcpy(coded, seq, seqLength);
        }
        out.rawBases = true;
    }
    return out;
}

void seqSquashStretch(const uint8_t* squash, uint8_t* out, uint32_t outLen, uint64_t actgPos)
{
    uint32_t i = 0;
    for (; i < outLen && (actgPos & 0x3u) != 0; ++i, ++actgPos) {
        out[i] = (uint8_t)((squash[actgPos >> 2] >> (6 - 2 * (actgPos & 0x3u))) & 0x3);
    }
    for (; i + 4 <= outLen; i += 4, actgPos += 4) {
        const uint8_t b = squash[actgPos >> 2];
        out[i + 0] = (uint8_t)((b >> 6) & 0x3);
        out[i + 1] = (uint8_t)((b >> 4) & 0x3);
        out[i + 2] = (uint8_t)((b >> 2) & 0x3);
        out[i + 3] = (uint8_t)(b & 0x3);
    }
    for (; i < outLen; ++i, ++actgPos) {
        out[i] = (uint8_t)((squash[actgPos >> 2] >> (6 - 2 * (actgPos & 0x3u))) & 0x3);
    }
}

/* The same walk, against a block's own consensus instead of a file's reference. */
bool buildSeqSquashCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                              uint32_t seqLength, int64_t refPos, const uint8_t* squash,
                              uint8_t* ref2bitScratch, uint8_t* out)
{
    return seqCodedBasesWalk(ops, seq, seqLength, refPos,
                             [squash](uint8_t* dst, uint32_t len, int64_t pos) {
                                 seqSquashStretch(squash, dst, len, (uint64_t)pos);
                             },
                             ref2bitScratch, out);
}

bool seqRecordUsesSquashWindow(uint16_t chrId, uint16_t flag, uint64_t startPos,
                               uint32_t refConsumed, int64_t windowStart, int64_t windowBases,
                               int64_t& refPos)
{
    refPos = 0;
    if (chrId == SEQ_CHR_ID_NONE || chrId == SEQ_CHR_ID_UNKNOWN || (flag & SEQ_FLAG_UNMAPPED)) {
        return false;
    }
    const int64_t chrStartPos = SamInfo::getInstance().getPositionByIndex(chrId);
    if (chrStartPos == -1) {
        return false;
    }
    refPos = chrStartPos + startPos - 1 - windowStart;   /* relative to the window's first base */
    if (refPos < 0) {
        return false;
    }
    const uint64_t needSquash = (uint64_t)(refConsumed >> 2) + !!(refConsumed & 0x3) + 1;
    const uint64_t windowSquashLen = (uint64_t)(windowBases >> 2) + !!(windowBases & 0x3);
    if (((uint64_t)(refPos >> 2)) + needSquash > windowSquashLen) {
        return false;
    }
    return true;
}

SeqRecordPayload buildSeqSquashRecordPayload(uint16_t chrId, uint16_t flag, uint64_t startPos,
                                             uint32_t refConsumed, const std::vector<CigarOp>& ops,
                                             const uint8_t* seq, uint32_t seqLength,
                                             const uint8_t* squash, int64_t windowStart,
                                             int64_t windowBases, uint8_t* coded,
                                             uint8_t* ref2bit)
{
    SeqRecordPayload out;
    out.bytes = coded;
    out.length = seqLength;

    int64_t refPos = 0;
    const bool useReference =
        seqRecordUsesSquashWindow(chrId, flag, startPos, refConsumed, windowStart, windowBases,
                                  refPos);
    if (useReference && !ops.empty()) {
        if (buildSeqSquashCodedBases(ops, seq, seqLength, refPos, squash, ref2bit, coded)) {
            out.usedReference = true;
            out.refPos = refPos + windowStart;
        } else {
            actgEncode(seq, coded, seqLength);
        }
    } else {
        if (seqLength > 0) {
            std::memcpy(coded, seq, seqLength);
        }
        out.rawBases = true;
    }
    return out;
}

void buildRunForm(const SeqExceptionClass& exc, std::vector<uint8_t>& out)
{
    for (uint32_t gap : exc.runGaps) {
        appendVarint(out, gap);
    }
    for (uint32_t len : exc.runLens) {
        appendVarint(out, len);
    }
}
