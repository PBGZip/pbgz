/*
 * sam_seq_payload.h - the reference-coded SEQ payload: how it is built, how it is split, and the
 * three forms its N positions can travel in.
 *
 * This is the encoder-side counterpart of sam_field_layout.h (which holds the readers). It exists
 * because two stages need the same machinery and neither owns it:
 *
 *   - the codec pre-selection trials the SEQ match stream on a sample of the first block, to pick
 *     the coder the column will be written with and the form its N positions will travel in;
 *   - the block pass walks the CIGAR and writes that payload, block by block.
 *
 * A trial measured on a different payload than the one written is worthless, so the payload's
 * construction is here rather than in either caller: what one record contributes
 * (buildSeqRecordPayload), how the stream splits (splitSeqMatchStream) and what the N positions'
 * three forms look like (SeqExceptionClass, buildRunForm) are single implementations. What stays
 * with the callers is only what genuinely differs - where each gets its record's columns and CIGAR
 * from (the block pass has the field decoders' parsed tables, the trial parses the block's text),
 * and what each does with the result (the block pass updates the reference's match statistics and
 * collects the exception characters; the trial only measures).
 *
 * Nothing here measures anything: choosing between candidates is the codec pre-selection's job,
 * and it is the only stage that runs a trial at all. What the block pass needs is here, what a
 * trial needs to build a comparable payload is here, and nothing else is.
 *
 * Nothing here knows about an actuator, a block or the meta: every function is a pure walk over
 * buffers, the reference and the CIGAR it is given.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "reference.h"

/*
 * One CIGAR operation: the op character and its length. The SEQ reference walk iterates these, and
 * the CIGAR column's own parse produces them.
 */
struct CigarOp {
    char op;
    uint32_t len;
};

/*
 * Split a CIGAR string into operations, and answer the reference span it consumes.
 *
 * Only the operations SAM defines are recorded (MIDNSHP=X); anything else in the field is skipped,
 * which is what the column parse has always done. The span adds up M/D/N/=/X - the operations that
 * consume reference sequence - and, defensively, their lower-case spellings, which the span has
 * always counted even though no valid CIGAR uses them.
 *
 * This is the one CIGAR parser: the CIGAR column's span and operation list, and the SEQ reference
 * walk's copy of them, all come through here. Only the span is wanted on some paths (the span goes
 * into the field decoders' tables whether or not anything walks the operations), so it has an entry
 * of its own rather than forcing a list to be built and thrown away. Both share one loop.
 */
uint32_t parseCigarOps(const uint8_t* cigar, uint32_t cigarLength, std::vector<CigarOp>& ops);
uint32_t cigarRefConsumed(const uint8_t* cigar, uint32_t cigarLength);

/*
 * The two halves of a match stream split into runs and values.
 *
 * The match stream is a sparse 0..3 byte stream in which a block whose reads all take their bases
 * from the reference is ~99% zeros. When the zero share is high enough to pay for it, the stream is
 * replaced by two independent sub-streams: this one (the varint run lengths of the zero runs) and
 * the surviving non-zero values, one byte each. The decision is splitSeqMatchStream's, and its two
 * halves are the "m" and "mval" sub-streams (see writeSeqMatchStreams).
 */
struct SeqRleSplit {
    bool useRle = false;
    std::unique_ptr<uint8_t[]> run;   /* varint run lengths of the zero runs */
    std::unique_ptr<uint8_t[]> val;   /* the surviving non-zero values, one byte each */
    uint32_t runLength = 0;
    uint32_t valLength = 0;
};

/*
 * Whether the match stream is split, and the split itself when it is.
 *
 * The verdict is per block and depends only on the payload: RLE is taken when at least 98% of the
 * bytes are zero. Values are written into the two buffers owned by the result; the caller only
 * reads them.
 */
SeqRleSplit splitSeqMatchStream(const uint8_t* match, uint32_t matchLen);

/*
 * Whether this record's bases can be coded against the reference at all, and where its first base
 * sits in it.
 *
 * Three things have to hold: the record is mapped to a chromosome (0xFFFF and 0xFFFE are the SAM
 * sentinels for "none" and "not in the dictionary"), it is not flagged unmapped, and the span it
 * would consume stays inside the reference - an alignment running past the end would read squashed
 * bytes that are not there. `refConsumed` is the reference span of the record's CIGAR (M/D/N/=/X),
 * which the caller has from the CIGAR parsing.
 *
 * The chromosome's offset in the reference comes from the process-wide table (SamInfo), which is
 * filled from the file's @SQ lines - in the preprocessing stage by the reader thread, and again by
 * whichever actuator pass sees the header.
 */
bool seqRecordUsesReference(uint16_t chrId, uint16_t flag, uint64_t startPos, uint32_t refConsumed,
                            Reference* reference, int64_t& refPos);

/*
 * Walk the CIGAR and write the record's per-base 2-bit payload: XOR against the reference on the
 * operations that consume it (M/=/X), the base's own 2-bit value on the ones that consume SEQ only
 * (I/S), nothing at all on reference-only (D/N) and on clipping. Walking the CIGAR is what keeps an
 * indel from being compared against a reference window the read never aligned to, and the decoder
 * mirrors this same walk.
 *
 * Answers false when the walk did not consume exactly `seqLength` bases - SEQ and CIGAR disagree, or
 * the CIGAR runs past the SEQ - which the caller answers by encoding the bases as they are.
 *
 * `ref2bitScratch` is the caller's scratch for the reference stretch (one op never exceeds the
 * read's length plus the coder's unaligned-write slack); `out` receives `seqLength` bytes.
 */
bool buildSeqReferenceCodedBases(const std::vector<CigarOp>& ops, const uint8_t* seq,
                                 uint32_t seqLength, int64_t refPos, Reference* reference,
                                 uint8_t* ref2bitScratch, uint8_t* out);

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

/*
 * Code one record's bases into `coded`, the payload slot they occupy.
 *
 * This is the decision the trial and the block pass both have to make identically: a trial
 * measured on a different payload than the one the block pass writes is worthless, so the branch
 * itself lives here rather than in each caller. What differs between the two is only what they
 * bring in (the block pass has the field decoders' parsed columns, the trial parses the block's
 * text) and what they do with the result (the block pass updates the reference's statistics and
 * collects the exception characters; the trial only measures).
 *
 * A record that can consult the reference (mapped, with CIGAR operations to walk) gets one 2-bit
 * code per base, the XOR against the reference over those operations; one that cannot writes its
 * own characters instead. `coded` must hold `seqLength` bytes and `ref2bit` the record's length
 * plus the coder's unaligned-write slack.
 */
SeqRecordPayload buildSeqRecordPayload(uint16_t chrId, uint16_t flag, uint64_t startPos,
                                       uint32_t refConsumed, const std::vector<CigarOp>& ops,
                                       const uint8_t* seq, uint32_t seqLength,
                                       Reference* reference, uint8_t* coded, uint8_t* ref2bit);

/*
 * Encoder-side accumulator for one exception character: a strictly increasing list of its positions
 * in the variable-length delta form above, plus - for the one class whose form is a file-level
 * decision - the absolute form and the runs of consecutive positions.
 *
 * The positions are written in one of three forms, and the stream name says which one:
 *   "nposd" - forward deltas, varint (every character; the default)
 *   "npos"  - absolute 4-byte block offsets
 *   "nposr" - runs: one gap and one length per run of consecutive positions
 * Only the 'N' class keeps more than one of them, because which one is written is the file-level
 * verdict PreprocessInfo::nposForm: the accumulator keeps all three so that any of them can be
 * produced, whichever way the pre-selection's trial came out (see
 * SamCodecActuator::writeSeqExceptionStreams). A lone "npos" with no "ch" is the first layout's
 * single list, which held 'N' and 'n' together and is still read; see the decoder in
 * sam_field_layout.h, where the names themselves live too.
 */
struct SeqExceptionClass {
    /* Set on the class the verdict may ask for a form of; positions and runs are kept for it and
       nothing else, since no other class ever needs them. */
    bool keepForms = false;
    std::vector<uint32_t> abs;
    std::vector<uint8_t> varint;
    /* Runs of consecutive positions: one gap (from the end of the previous run) and one length
       each. The N positions of ERR14949932 arrive in runs of 34.3 on average, 98% of them
       exactly 35 long, so describing them as runs makes the source ~1% of the size of the
       one-varint-per-position form before the coder sees either. */
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
