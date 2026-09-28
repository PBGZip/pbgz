/*
 * field_coder_config.h - Per-field encoder configuration table for SAM
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

#include <algorithm>
#include <cstdint>
#include <vector>

#include "pbgz_types.h"
#include "preprocess_info.h"

/*
 * Per-field encoder configuration table for SAM, split by compression profile:
 *
 *   archive (candidates / fallback)
 *       ratio-oriented coders; used by -m archive (the default). These are the
 *       candidates trial-compressed during preprocessing and the actuator's
 *       fallback when nothing is selected / selection failed / selection was
 *       not wired up.
 *
 *   fast (fastCandidates / fastFallback)
 *       speed-oriented coders; used by -m fast. The fast path is meant to stay
 *       quick, so its list holds cheap coders only. With a single entry the
 *       trial is a formality (and the encoder is deterministic); adding a
 *       second fast coder to the list makes the selector compare them without
 *       any actuator change.
 *
 * Changing a field's candidates or default touches only this table. Both the
 * preprocessing trial scope and the actuator defaults read from here.
 *
 * Which coders may actually be selected in a run is the intersection of this
 * table and CoderFactory's descriptor table (the profile and level gate). A
 * coder with no row in the registry, or one whose profile does not match the
 * current -m mode, is filtered out before the trial.
 *
 * Special cases:
 *  - QUAL(10): goes through the QualSelector-specific path, which reads its
 *    candidates from the QUAL row's archive list rather than from the generic
 *    trial (the quality column needs record boundaries and a per-record cycle
 *    index, so it cannot be trialled as a plain byte stream). Dropping a coder
 *    from the QUAL row is how it is taken out of the QUAL decision.
 *  - POS(3) / PNEXT(7) / TLEN(8): always use delta / inferred compression (see
 *    compressPosFieldDelta / compressPNextFieldDelta / compressTLen). POS still
 *    selects its underlying entropy coder, but on the rebuilt varint stream
 *    rather than on raw text, so it keeps its own candidate list in both
 *    profiles; PNEXT/TLEN have no generic candidates at all.
 */

struct FieldCoderConfig {
    /* Archive profile (default -m archive). */
    std::vector<CoderType> candidates;
    CoderType fallback;

    /* Fast profile (-m fast). */
    std::vector<CoderType> fastCandidates;
    CoderType fastFallback;
};

/* Number of SAM fields participating in coder selection: 11 mandatory fields + the 12th OPTION column. */
static const uint32_t SAM_FIELD_COUNT_SELECT = SAM_FIELD_COUNT + 1;

/*
 * The fast list is speed-oriented: coder_rans is a static order-0 entropy coder
 * with per-symbol O(1) cost and no BWT/context model, which is exactly what the
 * fast preset wants. Every field's fast fallback is coder_rans because that is
 * what the fast structured path has always used for generic fields; POS and the
 * inferred fields (PNEXT/TLEN) have no generic fast candidate list because they
 * are not fed to the generic trial (POS is evaluated on its varint delta stream,
 * see selectPosDeltaCoder; PNEXT/TLEN are differenced/inferred). QUAL's fast
 * list is empty because its dedicated selector keeps using the archive row, and
 * its fast fallback here is therefore unused; the FASTQ table's QUAL row below is
 * the one place that reads a fast QUAL default, and it is coder_qual rather than
 * coder_rans.
 */
#define PBGZ_FQ_RANS   std::vector<CoderType>{CoderType::RANS}
#define PBGZ_FQ_NONE   std::vector<CoderType>{}

inline const FieldCoderConfig kSamFieldCoderConfig[SAM_FIELD_COUNT_SELECT] = {
    /* QNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* FLAG  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* RNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* POS   */ {{CoderType::BWT_CM, CoderType::ARITH}, CoderType::BWT_CM, PBGZ_FQ_NONE, CoderType::RANS},
    /* MAPQ  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* CIGAR */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* RNEXT */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* PNEXT */ {{}, CoderType::BWT_CM, PBGZ_FQ_NONE, CoderType::RANS},
    /* TLEN  */ {{}, CoderType::BWT_CM, PBGZ_FQ_NONE, CoderType::RANS},
    /* SEQ   */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::FC, PBGZ_FQ_RANS, CoderType::RANS},
    /* QUAL  */ {{CoderType::QUAL, CoderType::FCV2, CoderType::BWT_CM, CoderType::QCM}, CoderType::QUAL, PBGZ_FQ_NONE, CoderType::RANS},
    /* OPTION */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
};

/*
 * The FASTQ table: one record is four lines - ID / SEQ / QUAL / comment (see FastqField).
 *
 * It is a table of its own because the two formats' columns are not the same material even
 * where they are analogues: a FASTQ ID is a whole '@...' line rather than one tab-delimited
 * field, its SEQ has no CIGAR beside it to be walked against a reference, and its comment line
 * is usually not written at all. So their candidates and defaults are declared per format
 * rather than shared.
 *
 *   ID      the '@...' line, written line by line and segmented at its separators (see
 *           FastqCodecActuator::compressIdInSplit) - the layout the segmentation coder gains
 *           on, which is why affix is a candidate. coder_fc is not: it wants the column as one
 *           block, so a verdict for it could not be honored here.
 *   SEQ     written whole when there is no reference, so coder_fc is a candidate - the actuator
 *           takes its whole-block path then - subject to its size window (FC_MIN_LEN/
 *           FC_MAX_LEN), outside which the column falls back to bwt_cm.
 *   QUAL    goes through the dedicated QualSelector path, which reads its candidates from this
 *           row - dropping a coder here is how it is taken out of the QUAL decision. The fast
 *           profile leaves the list empty, so no trial runs, and keeps its default at
 *           coder_qual: a quality column is largely a function of the read's bases, which only
 *           a record-level coder can see, so coding QUAL as a plain byte stream (coder_rans)
 *           can come out larger than gzip on QUAL-heavy short reads. coder_qual is the
 *           column's original coder and stays its default in fast mode.
 *   COMMENT only written when the line is neither '+' alone nor a copy of the ID line, so most
 *           files never encode this column and its verdict is simply unused. Written line by
 *           line, like the ID, so its candidates match the ID's.
 */
inline const FieldCoderConfig kFastqFieldCoderConfig[FQ_FIELD_COUNT] = {
    /* ID      */ {{CoderType::BWT_CM, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* SEQ     */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* QUAL    */ {{CoderType::QUAL, CoderType::FCV2, CoderType::BWT_CM, CoderType::QCM}, CoderType::QUAL, PBGZ_FQ_NONE, CoderType::QUAL},
    /* COMMENT */ {{CoderType::BWT_CM, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
};

#undef PBGZ_FQ_RANS
#undef PBGZ_FQ_NONE

inline const FieldCoderConfig* samFieldCoderConfig(uint32_t fieldIdx)
{
    if (fieldIdx >= SAM_FIELD_COUNT_SELECT) {
        return nullptr;
    }
    return &kSamFieldCoderConfig[fieldIdx];
}

/* Empty list returned for out-of-range fields, so callers can always iterate. */
inline const std::vector<CoderType>& emptyCoderList()
{
    static const std::vector<CoderType> kEmpty;
    return kEmpty;
}

/* Candidate list for a field under the given compression mode. */
inline const std::vector<CoderType>& samFieldCandidates(uint32_t fieldIdx, uint8_t mode)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return emptyCoderList();
    }
    return (mode == PBGZ_MODE_FAST) ? cfg->fastCandidates : cfg->candidates;
}

/* Whether a field lists the given coder as a candidate in the given mode. */
inline bool samFieldCandidate(uint32_t fieldIdx, CoderType type, uint8_t mode)
{
    const std::vector<CoderType>& list = samFieldCandidates(fieldIdx, mode);
    return std::find(list.begin(), list.end(), type) != list.end();
}

/* Archive-profile membership (historical accessor; the QUAL path uses it). */
inline bool samFieldCandidate(uint32_t fieldIdx, CoderType type)
{
    return samFieldCandidate(fieldIdx, type, PBGZ_MODE_ARCHIVE);
}

/*
 * Whether the QUAL column may use this coder: the QUAL selector trials exactly
 * the coders the QUAL row above lists, and the row is the only place that decides
 * it.
 *
 * bwt_cm is in the set because it is the best of the three on part of the data.
 * On aligned -l 7 BAM, short-read QUAL codes to coder_qual 34.11% @7MB/s, fcv2
 * 26.21% @4MB/s, bwt_cm 26.56% @6MB/s; long-read QUAL to 50.16% @13MB/s, 48.78%
 * @5MB/s, 48.57% @7MB/s. So fcv2 wins on short reads, and bwt_cm wins on long ones
 * where it is also the faster of the two. It is likewise the coder the QUAL column
 * falls back to where fcv2 is not applicable (an unaligned or single-record block,
 * see CoderFactory::coderSupports).
 *
 * The row is consulted, so there is no separate switch for this: taking a coder
 * out of the QUAL decision means removing it from the QUAL row above.
 */
inline bool qualCoderCandidate(CoderType type)
{
    return samFieldCandidate(SAM_QUAL, type);
}

/*
 * The FASTQ accessors, mirroring the SAM ones: the candidate list for a field under the current
 * -m mode, membership in it, and the mode's default (see qualCoderCandidate above for what the
 * QUAL row decides).
 */
inline const FieldCoderConfig* fastqFieldCoderConfig(uint32_t fieldIdx)
{
    if (fieldIdx >= FQ_FIELD_COUNT) {
        return nullptr;
    }
    return &kFastqFieldCoderConfig[fieldIdx];
}

inline const std::vector<CoderType>& fastqFieldCandidates(uint32_t fieldIdx, uint8_t mode)
{
    const FieldCoderConfig* cfg = fastqFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return emptyCoderList();
    }
    return (mode == PBGZ_MODE_FAST) ? cfg->fastCandidates : cfg->candidates;
}

inline bool fastqFieldCandidate(uint32_t fieldIdx, CoderType type, uint8_t mode)
{
    const std::vector<CoderType>& list = fastqFieldCandidates(fieldIdx, mode);
    return std::find(list.begin(), list.end(), type) != list.end();
}

inline CoderType fastqFieldDefaultCoder(uint32_t fieldIdx, uint8_t mode, CoderType fallback)
{
    const FieldCoderConfig* cfg = fastqFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return fallback;
    }
    return (mode == PBGZ_MODE_FAST) ? cfg->fastFallback : cfg->fallback;
}

/* Archive-profile default coder; returns the caller-supplied fallback when nothing is registered. */
inline CoderType samFieldDefaultCoder(uint32_t fieldIdx, CoderType fallback)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    return (cfg != nullptr) ? cfg->fallback : fallback;
}

/* Default coder for the given compression mode. */
inline CoderType samFieldDefaultCoder(uint32_t fieldIdx, uint8_t mode, CoderType fallback)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return fallback;
    }
    return (mode == PBGZ_MODE_FAST) ? cfg->fastFallback : cfg->fallback;
}
