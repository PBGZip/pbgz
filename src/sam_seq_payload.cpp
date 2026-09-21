/*
 * sam_seq_payload.cpp - the reference-coded SEQ payload (see sam_seq_payload.h).
 */

#include "sam_seq_payload.h"

#include <cstring>

#include "actg.h"
#include "log/logger.h"
#include "sam_field_layout.h"
#include "sam_info.h"

/* The one CIGAR walk; `ops` is filled only when the caller wants the operations. */
static uint32_t parseCigarImpl(const uint8_t* cigar, uint32_t cigarLength, std::vector<CigarOp>* ops)
{
    if (ops != nullptr) {
        ops->clear();
    }
    if (cigar == nullptr || cigarLength == 0) {
        return 0;
    }
    uint32_t refConsumed = 0;
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
                case 'M': case 'D': case 'N': case '=': case 'X':
                case 'm': case 'd': case 'n':
                    refConsumed += currentNumber;
                    break;
                default:
                    break;
            }
            currentNumber = 0;
        }
    }
    return refConsumed;
}

uint32_t parseCigarOps(const uint8_t* cigar, uint32_t cigarLength, std::vector<CigarOp>& ops)
{
    return parseCigarImpl(cigar, cigarLength, &ops);
}

uint32_t cigarRefConsumed(const uint8_t* cigar, uint32_t cigarLength)
{
    return parseCigarImpl(cigar, cigarLength, nullptr);
}

SeqRleSplit splitSeqMatchStream(const uint8_t* match, uint32_t matchLen)
{
    SeqRleSplit split;
    uint32_t nNonZero = 0;
    for (uint32_t i = 0; i < matchLen; i++) {
        if (match[i] != 0) {
            nNonZero++;
        }
    }
    const uint32_t zeroCount = matchLen - nNonZero;
    split.useRle = (matchLen > 0) &&
                   ((uint64_t)zeroCount * 100ull >= (uint64_t)matchLen * 98ull);
    LOG_DEBUG("SEQ match: matchLen=%u zeros=%u (%.4f) nonZero=%u -> RLE split %s",
              matchLen, zeroCount, (double)zeroCount / (double)(matchLen ? matchLen : 1),
              nNonZero, split.useRle ? "on" : "off");
    if (!split.useRle) {
        return split;
    }

    split.run = std::make_unique<uint8_t[]>((nNonZero + 1) * 5 + 16);
    split.val = std::make_unique<uint8_t[]>(nNonZero + 16);
    uint32_t rp = 0;
    uint32_t vp = 0;
    uint32_t run = 0;
    for (uint32_t i = 0; i < matchLen; i++) {
        if (match[i] == 0) {
            run++;
        } else {
            uint32_t v = run;
            while (v >= 0x80) { split.run[rp++] = (uint8_t)(v | 0x80); v >>= 7; }
            split.run[rp++] = (uint8_t)v;
            split.val[vp++] = match[i];
            run = 0;
        }
    }
    uint32_t tail = run;                /* trailing run of zeros */
    while (tail >= 0x80) { split.run[rp++] = (uint8_t)(tail | 0x80); tail >>= 7; }
    split.run[rp++] = (uint8_t)tail;

    split.runLength = rp;               /* length of the run-length segment */
    split.valLength = vp;               /* length of the value segment */
    return split;
}

bool seqRecordUsesReference(uint16_t chrId, uint16_t flag, uint64_t startPos, uint32_t refConsumed,
                            Reference* reference, int64_t& refPos)
{
    refPos = 0;
    if (chrId == 0xFFFF || chrId == 0xFFFE || (flag & 0x04)) {
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

bool buildSeqReferenceCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                                 uint32_t seqLength, int64_t refPos, Reference* reference,
                                 uint8_t* ref2bitScratch, uint8_t* out)
{
    uint32_t readPos = 0;
    int64_t refPosLocal = refPos;
    for (size_t oi = 0; oi < ops.size(); ++oi) {
        const CigarOp& op = ops[oi];
        switch (op.op) {
            case 'M': case '=': case 'X':
                if (readPos + op.len > seqLength) {
                    return false;
                }
                reference->getStretch2Bits1Char(ref2bitScratch, op.len, refPosLocal);
                for (uint32_t i = 0; i < op.len; ++i) {
                    const uint8_t read2 = (seq[readPos + i] >> 1) & 0x3;
                    out[readPos + i] = read2 ^ ref2bitScratch[i];
                }
                readPos += op.len;
                refPosLocal += op.len;
                break;
            case 'I': case 'S':
                if (readPos + op.len > seqLength) {
                    return false;
                }
                for (uint32_t i = 0; i < op.len; ++i) {
                    out[readPos + i] = (seq[readPos + i] >> 1) & 0x3;
                }
                readPos += op.len;
                break;
            case 'D': case 'N':
                refPosLocal += op.len;
                break;
            case 'H': case 'P':
            default:
                break; // consume neither SEQ nor reference
        }
    }
    return readPos == seqLength;
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
            /* CIGAR and SEQ disagree, or the CIGAR runs past the SEQ: the bases go in as their
               2-bit codes anyway. The walk did not run, so the caller must not count this record
               as reference-coded (its bases were not compared against anything). */
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
