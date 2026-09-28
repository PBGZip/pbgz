/*
 * sam_seq_payload.cpp - the reference-coded SEQ payload (see sam_seq_payload.h).
 */

#include "sam_seq_payload.h"

#include <atomic>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>

#include "actg.h"
#include "log/logger.h"
#include "sam_field_layout.h"
#include "sam_info.h"

/*
 * The one CIGAR parse: the reference span, the read span, and the operation list all come out of it.
 * `ops` (and the spans) are filled only when the caller wants them.
 *
 * The two spans keep the letter sets they have always had, lower-case forms included, because both
 * feed decisions that the archive records: the read span is what the SEQ payloads are measured
 * against, and the reference span decides whether a record has a window to code against at all. They
 * therefore disagree about exactly one letter - a lower-case 'x' counts towards the read span only -
 * which is a difference to keep rather than to fix quietly, and the comment on each group says so.
 */
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

/*
 * splitSeqMatchStream now lives in seq_stream_util.h, as an inline the FASTQ SEQ column uses too;
 * the 98% threshold it applies and the measurements behind it are documented there.
 */

bool seqRecordUsesReference(uint16_t chrId, uint16_t flag, uint64_t startPos, uint32_t refConsumed,
                            Reference* reference, int64_t& refPos)
{
    refPos = 0;
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
        /*
         * The squash came from an NI index, which carries no names. Its layout has always been
         * taken on trust, and callers relying on that must keep working, so an unverifiable
         * reference is accepted rather than refused.
         */
        return true;
    }

    /*
     * Where each name really starts. The offset a SEQ lookup will use is the running sum of the
     * @SQ lengths in header order (SamInfo::calculateChromosomePositions), so the two agree only
     * when the FASTA presents the same sequences in the same order.
     */
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
    /*
     * Worth saying out loud: the reference is loaded and its checksum verified, so nothing else
     * in the run reports a problem - the SEQ column just quietly stops improving, by 20% of the
     * archive on a single-end BAM measured here (471 MB against 434 MB, once PNEXT's own waste
     * is set aside).
     */
    LOG_WARNING("The reference does not line up with the file's @SQ list: its sequences are not "
                "the ones the header names, in the header's order, so a read's position would "
                "address the wrong bases. SEQ is compressed without it; supply a reference of "
                "the same assembly as the alignment.");
    fprintf(stderr, "warning: reference does not line up with the @SQ list; "
                    "SEQ falls back to compression without a reference\n");
}

/*
 * Walk one record's ops and produce its coded SEQ bases: each read base as its 2-bit value, XORed
 * with the reference's 2-bit base where the op consumes reference (M/=/X) and as it stands where it
 * does not (I/S).
 *
 * Only the reference side varies between the two callers below - the file's genome and a block's own
 * consensus both present their bases as 2-bit, but are fetched differently - so it arrives as a
 * callable. A template parameter rather than std::function: this runs once per read of every block,
 * and an indirect call per op is exactly what the walk should not pay.
 *
 * Answers false when the ops do not describe the read (an op runs past seqLength, or the ops do not
 * consume it exactly), which is what makes the caller store the bases as they stand instead.
 */
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
