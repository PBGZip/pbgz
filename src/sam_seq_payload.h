/*
 * sam_seq_payload.h - the reference-coded SEQ payload: how it is built, split, and how its N
 * positions travel. The encoder-side counterpart of sam_field_layout.h.
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

#include "reference.h"

/* One CIGAR operation: the op character and its length. */
struct CigarOp {
    char op;
    uint32_t len;
};

/* Split a CIGAR string into operations and answer the reference span it consumes (M/D/N/=/X).
   The one CIGAR parser: the span-only entry shares the loop with the operation-list one. */
uint32_t parseCigarOps(const uint8_t* cigar, uint32_t cigarLength, std::vector<CigarOp>& ops);
uint32_t cigarRefConsumed(const uint8_t* cigar, uint32_t cigarLength);

/* The span the CIGAR claims of the read: M/I/S/=/X and lower-case forms. */
uint32_t cigarSeqConsumed(const uint8_t* cigar, uint32_t cigarLength);

/* The run/value split (SeqRleSplit / splitSeqMatchStream) lives in seq_stream_util.h. */

/* Whether this record's bases can be coded against the reference, and where its first base sits
   in it. Requires a mapped chromosome, no unmapped flag, and a span inside the reference. */
/* SAM's placeholder chromosome ids (0xFFFF "*", 0xFFFE unknown) and the FLAG unmapped bit. */
inline constexpr uint16_t SEQ_CHR_ID_NONE = 0xFFFF;
inline constexpr uint16_t SEQ_CHR_ID_UNKNOWN = 0xFFFE;
inline constexpr uint8_t SEQ_FLAG_UNMAPPED = 0x04;

bool seqRecordUsesReference(uint16_t chrId, uint16_t flag, uint64_t startPos, uint32_t refConsumed,
                            Reference* reference, int64_t& refPos);

/* Whether the loaded reference is laid out the way the @SQ list addresses it. Returns true when
   the correspondence holds or cannot be checked (an NI-loaded reference carries no names). */
bool referenceLinesUpWithHeader(const Reference* reference);

/* Tell the user, once per run, that the reference cannot be used for SEQ. */
void warnReferenceNotUsableForSeq();


/* Walk the CIGAR writing the record's per-base 2-bit payload: XOR against the reference on
   M/=/X, the base's own code on I/S, nothing on D/N and clipping. `out` receives `seqLength`
   bytes; false when the walk did not consume the read exactly. */
/* What a callback tells the walk: keep going, or stop here. */
enum class SeqWalkAction { Continue, Stop };

/* Where a walk left off, and whether the ops described the read. */
struct SeqWalkEnd {
    uint32_t readPos = 0;      /* read bases consumed */
    int64_t  refPos = 0;       /* reference bases consumed, from the walk's own starting point */
    bool     complete = false; /* the ops consumed the read exactly, and asked for no stop */
};

/* Walk one record's ops in read order, handing each to onMatch(readPos, refPos, len) for M/=/X,
   onInsert(readPos, len) for I/S, onDelete(len) for D/N. H/P consume neither. An op running past
   seqLength ends the walk as incomplete. */

template <class MatchFn, class InsertFn, class DeleteFn>
static inline SeqWalkEnd seqWalkCigar(const std::vector<CigarOp>& ops, uint32_t seqLength,
                                      int64_t refPos, MatchFn onMatch, InsertFn onInsert,
                                      DeleteFn onDelete)
{
    SeqWalkEnd end;
    end.refPos = refPos;
    for (size_t oi = 0; oi < ops.size(); ++oi) {
        const CigarOp& op = ops[oi];
        switch (op.op) {
            case 'M': case '=': case 'X':
                if (end.readPos + op.len > seqLength) {
                    return end;
                }
                if (onMatch(end.readPos, end.refPos, op.len) == SeqWalkAction::Stop) {
                    return end;
                }
                end.readPos += op.len;
                end.refPos += op.len;
                break;
            case 'I': case 'S':
                if (end.readPos + op.len > seqLength) {
                    return end;
                }
                if (onInsert(end.readPos, op.len) == SeqWalkAction::Stop) {
                    return end;
                }
                end.readPos += op.len;
                break;
            case 'D': case 'N':
                if (onDelete(op.len) == SeqWalkAction::Stop) {
                    return end;
                }
                end.refPos += op.len;
                break;
            case 'H': case 'P':
            default:
                break; // consume neither SEQ nor reference
        }
    }
    end.complete = (end.readPos == seqLength);
    return end;
}

bool buildSeqReferenceCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                                 uint32_t seqLength, int64_t refPos, Reference* reference,
                                 uint8_t* ref2bitScratch, uint8_t* out);

/* The same walk against a caller-owned squash buffer instead of the file's reference, where
   squash[0] stands for base offset `refPos`. */

bool buildSeqSquashCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                              uint32_t seqLength, int64_t refPos, const uint8_t* squash,
                              uint8_t* ref2bitScratch, uint8_t* out);

/* The windowed twin of seqRecordUsesReference: the first base is rebased on [windowStart,
   windowStart + windowBases), or false when the record cannot be coded against it. */

bool seqRecordUsesSquashWindow(uint16_t chrId, uint16_t flag, uint64_t startPos,
                               uint32_t refConsumed, int64_t windowStart, int64_t windowBases,
                               int64_t& refPos);

/*
 * Read `outLen` 2-bit codes out of a squash buffer (four per byte, the first base in the high
 * bits) starting at base offset `actgPos`. The caller has bounds-checked the stretch. Shared by
 * the encoder's walk above and the decoder's, which restores bases the same way.
 */
void seqSquashStretch(const uint8_t* squash, uint8_t* out, uint32_t outLen, uint64_t actgPos);

/*
 * What one record contributes to the match stream, once its bases have been coded.
 */
struct SeqRecordPayload {
    const uint8_t* bytes = nullptr;   /* the bytes to append (the caller's `coded` buffer) */
    uint32_t length = 0;
    /* The payload holds the record's own characters rather than 2-bit codes. The two do not
       overlap (a 2-bit code is 0..3, a SEQ character is at least '='), so the decoder tells them
       apart by asking the same question the encoder asked. */
    bool rawBases = false;
    /* The reference walk ran: the record's bases are XORed against the reference, and the caller
       owes the reference its match statistics (updateMatchedGene). False both for a record that
       cannot consult the reference and for one whose CIGAR and SEQ disagree, since neither ran the
       walk - which is also what the caller's exception loop needs to know, a record that could not
       use the reference carrying its plain N in the payload. */
    bool usedReference = false;
    int64_t refPos = 0;               /* where its first base sits, for that update */
};

/* Code one record's bases into `coded`. A record that can consult the reference gets one 2-bit
   code per base (the XOR against the reference over its ops); one that cannot writes its own
   characters. `coded` must hold `seqLength` bytes. */
SeqRecordPayload buildSeqRecordPayload(uint16_t chrId, uint16_t flag, uint64_t startPos,
                                       uint32_t refConsumed, const std::vector<CigarOp>& ops,
                                       const uint8_t* seq, uint32_t seqLength,
                                       Reference* reference, uint8_t* coded, uint8_t* ref2bit);

/* The same decision against a block-local squash window. `refPos` in the result stays absolute in
   the file's coordinate space; no reference match statistics are updated. */
SeqRecordPayload buildSeqSquashRecordPayload(uint16_t chrId, uint16_t flag, uint64_t startPos,
                                             uint32_t refConsumed, const std::vector<CigarOp>& ops,
                                             const uint8_t* seq, uint32_t seqLength,
                                             const uint8_t* squash, int64_t windowStart,
                                             int64_t windowBases, uint8_t* coded,
                                             uint8_t* ref2bit);

/* Encoder-side accumulator for one exception character: positions as varint deltas, plus - for
   the class the form verdict may ask for - the absolute form and the runs of consecutive
   positions. Stream names: "nposd" (forward varint deltas, default), "npos" (absolute 4-byte),
   "nposr" (runs). */
struct SeqExceptionClass {
    /* Set on the class the form verdict may ask for; positions and runs are kept for it. */
    bool keepForms = false;
    std::vector<uint32_t> abs;
    std::vector<uint8_t> varint;
    /* Runs of consecutive positions: one gap and one length each. */
    std::vector<uint32_t> runGaps;
    std::vector<uint32_t> runLens;
    uint32_t runStart = 0;
    uint32_t runLen = 0;
    uint32_t prevRunEnd = 0;
    uint32_t count = 0;
    uint32_t last = 0;

    void add(uint32_t pos)
    {
        if (keepForms) {
            abs.push_back(pos);
        }
        uint32_t delta = pos - last;   /* strictly increasing within the block */
        if (keepForms) {
            if (delta == 1 && runLen != 0) {
                runLen++;
            } else {
                closeRun();
                runStart = pos;
                runLen = 1;
            }
        }
        last = pos;
        while (delta >= 0x80) {
            varint.push_back((uint8_t)(delta | 0x80));
            delta >>= 7;
        }
        varint.push_back((uint8_t)delta);
        count++;
    }

    /* Ends the run in progress, if any, so runGaps/runLens describe every position added. */
    void closeRun()
    {
        if (runLen != 0) {
            runGaps.push_back(runStart - prevRunEnd);
            runLens.push_back(runLen);
            prevRunEnd = last + 1;
            runLen = 0;
        }
    }
};

/* The run form's payload: every gap, then every length, each as a varint (see above). */
void buildRunForm(const SeqExceptionClass& exc, std::vector<uint8_t>& out);
