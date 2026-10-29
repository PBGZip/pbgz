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
#include <cctype>
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

    /* The reference names (and lengths) are what have to agree with the @SQ list. The order need
       not: both the encoder and the decoder address a chromosome at the running sum of the @SQ
       lengths, so a reordered FASTA keeps decompression correct - only the affected contigs are
       XORed against the wrong bases and compress less well. Requiring the same order made a
       reference whose unplaced contigs are ordered differently fall back to no reference at all,
       which costs far more than it protects. */
    std::unordered_map<std::string, Reference::SequenceSpan> refByName;
    refByName.reserve(seqs.size() * 2);
    for (const Reference::SequenceSpan& s : seqs) {
        refByName.emplace(s.name, s);
    }

    const std::vector<ChromosomeInfo>& chrs = SamInfo::getInstance().getAllChromosomeInfo();
    size_t reordered = 0;
    for (const ChromosomeInfo& c : chrs) {
        const auto it = refByName.find(c.name);
        if (it == refByName.end()) {
            return false;   /* the reference is not the assembly the header names */
        }
        if (it->second.length != 0 && it->second.length != c.length) {
            return false;   /* same name, different sequence length: not this assembly */
        }
        if (it->second.offset != c.position) {
            ++reordered;
        }
    }
    if (reordered > 0) {
        LOG_WARNING("Reference and @SQ list hold the same sequences in a different order "
                    "(%zu of %zu contigs); the reference is used with the @SQ layout, so those "
                    "contigs compress less well (decompression is unaffected).",
                    reordered, chrs.size());
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
                "the ones the header names (a name is missing or its length differs), so a read's "
                "position would address the wrong bases. SEQ is compressed without it; supply a "
                "reference of the same assembly as the alignment.");
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

/* Recompute MD:Z and NM:i exactly as BWA/samtools calmd do (see bam_fillmd1_core): mismatches and
   deletions extend MD, insertions/deletions count toward NM, soft/hard clips and reference skips
   do not, and an N/N or N/base comparison is a mismatch. Both sides run this same function on the
   same record, so the encoder can validate a regenerated value against the stored one and only
   keep an exception when they differ. */
bool deriveMdNm(Reference* reference, int64_t chrStart, uint64_t pos, uint16_t flag,
                const std::vector<CigarOp>& ops, const uint8_t* seq, uint32_t seqLen,
                std::vector<uint8_t>& refScratch, std::string& md, int64_t& nm)
{
    md.clear();
    nm = 0;
    if (reference == nullptr || seq == nullptr || seqLen == 0 ||
        ops.empty() || (flag & SEQ_FLAG_UNMAPPED) || pos == 0) {
        return false;
    }
    int64_t rp = chrStart + (int64_t)pos - 1;
    if (rp < 0) {
        return false;
    }
    const uint64_t squashLen = (uint64_t)reference->getSquashLength();
    uint32_t readPos = 0;
    uint32_t matched = 0;
    char num[16];
    const auto appendNum = [&](uint32_t v) {
        const int n = snprintf(num, sizeof(num), "%u", v);
        md.append(num, (size_t)n);
    };
    /* getStretchActg writes 4-(pos&3) bytes even for a short read, so always fetch from a
       4-aligned start and index past the lead. The scratch grows to the fetched length
       (a deletion may be longer than the read). */
    const auto fetchRef = [&](int64_t at, uint32_t len, uint32_t& lead) -> bool {
        const int64_t aligned = at & ~(int64_t)3;
        lead = (uint32_t)(at - aligned);
        const uint32_t fetchLen = (lead + len + 3u) & ~3u;
        if (aligned < 0 || ((uint64_t)(aligned >> 2)) + (fetchLen >> 2) > squashLen) {
            return false;
        }
        if (refScratch.size() < fetchLen) {
            refScratch.resize(fetchLen);
        }
        reference->getStretchActg(refScratch.data(), fetchLen, (uint64_t)aligned);
        return true;
    };

    for (const CigarOp& op : ops) {
        switch (op.op) {
            case 'M': case '=': case 'X': {
                uint32_t lead = 0;
                if ((uint64_t)readPos + op.len > (uint64_t)seqLen ||
                    !fetchRef(rp, op.len, lead)) {
                    return false;
                }
                /* One base, exactly as calmd compares it. */
                const auto consumeBase = [&](uint32_t k) {
                    const uint8_t rb = (uint8_t)std::toupper(refScratch[lead + k]);
                    const uint8_t cb = (uint8_t)std::toupper(seq[readPos + k]);
                    const bool isMatch = (cb == '=') || (cb == rb && cb != 'N' && rb != 'N');
                    if (isMatch) {
                        ++matched;
                    } else {
                        appendNum(matched);
                        matched = 0;
                        md.push_back((char)rb);
                        ++nm;
                    }
                };
                uint32_t i = 0;
#if defined(__aarch64__) && defined(__ARM_NEON)
                /* Most M blocks match end to end; compare 16 bases at a time and only run the
                   per-base path for a block that has a mismatch.  The canonical form below
                   (lower-case folded to upper, then the calmd rule) is exactly std::toupper
                   followed by that rule. */
                {
                    const uint8x16_t va = vdupq_n_u8('a'), vz = vdupq_n_u8('z');
                    const uint8x16_t vclear = vdupq_n_u8(0xDF);
                    const uint8x16_t vN = vdupq_n_u8('N');
                    const uint8x16_t vEqMark = vdupq_n_u8('=');
                    while (i + 16 <= op.len) {
                        const uint8x16_t cbRaw = vld1q_u8(seq + readPos + i);
                        const uint8x16_t rbRaw = vld1q_u8(refScratch.data() + lead + i);
                        const uint8x16_t cb = vbslq_u8(
                            vandq_u8(vcgeq_u8(cbRaw, va), vcleq_u8(cbRaw, vz)),
                            vandq_u8(cbRaw, vclear), cbRaw);
                        const uint8x16_t rb = vbslq_u8(
                            vandq_u8(vcgeq_u8(rbRaw, va), vcleq_u8(rbRaw, vz)),
                            vandq_u8(rbRaw, vclear), rbRaw);
                        const uint8x16_t isN = vorrq_u8(vceqq_u8(cb, vN), vceqq_u8(rb, vN));
                        const uint8x16_t isMatch = vorrq_u8(
                            vandq_u8(vceqq_u8(cb, rb), vmvnq_u8(isN)),
                            vceqq_u8(cb, vEqMark));
                        if (vminvq_u8(isMatch) == 0xFF) {
                            matched += 16;
                            i += 16;
                            continue;
                        }
                        for (uint32_t k = 0; k < 16; ++k) {
                            consumeBase(i + k);
                        }
                        i += 16;
                    }
                }
#endif
                for (; i < op.len; ++i) {
                    consumeBase(i);
                }
                readPos += op.len;
                rp += op.len;
                break;
            }
            case 'I':
                if ((uint64_t)readPos + op.len > (uint64_t)seqLen) {
                    return false;
                }
                readPos += op.len;
                nm += op.len;
                break;
            case 'S':
                if ((uint64_t)readPos + op.len > (uint64_t)seqLen) {
                    return false;
                }
                readPos += op.len;
                break;
            case 'D': {
                uint32_t lead = 0;
                if (!fetchRef(rp, op.len, lead)) {
                    return false;
                }
                appendNum(matched);
                matched = 0;
                md.push_back('^');
                for (uint32_t i = 0; i < op.len; ++i) {
                    md.push_back((char)std::toupper(refScratch[lead + i]));
                }
                rp += op.len;
                nm += op.len;
                break;
            }
            case 'N':
                rp += op.len;
                break;
            case 'H': case 'P':
                break;
            default:
                return false;
        }
    }
    if (readPos != seqLen) {
        return false;
    }
    appendNum(matched);   /* calmd always appends the trailing match run, even when it is 0 */
    return true;
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
