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
 * archive (candidates/fallback) for -m archive, fast (fastCandidates/fastFallback)
 * for -m fast. Changing a field's candidates or default touches only this table;
 * actual selectability is the intersection with CoderFactory's descriptor table.
 *
 * Special cases: QUAL goes through QualSelector (which reads its candidates from the
 * QUAL row); POS/PNEXT/TLEN always use delta/inferred compression.
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

/* The fast lists are speed-oriented (coder_rans). Every field's fast fallback is coder_rans;
   POS and the inferred fields (PNEXT/TLEN) have no generic fast candidate list, and QUAL's
   fast list is empty (its dedicated selector uses the archive row). */
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
    /* SEQ   */ {{CoderType::FC}, CoderType::FC, PBGZ_FQ_RANS, CoderType::RANS},
    /* QUAL  */ {{CoderType::QUAL, CoderType::FCV2, CoderType::QCM}, CoderType::QUAL, PBGZ_FQ_NONE, CoderType::RANS},
    /* OPTION */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
};

/*
 * The FASTQ table: one record is four lines - ID / SEQ / QUAL / comment (see FastqField).
 * A separate table because the two formats' columns are not the same material even where
 * analogous (a FASTQ ID is a whole '@...' line, its SEQ has no CIGAR, its comment is
 * usually absent).
 *   ID      '@...' line, written line by line and segmented at separators; affix is a
 *           candidate, coder_fc is not (it wants the column as one block).
 *   SEQ     written whole when there is no reference, so coder_fc is a candidate, subject
 *           to its size window.
 *   QUAL    goes through the dedicated QualSelector path; fast leaves the list empty and
 *           keeps coder_qual as default.
 *   COMMENT written only when it is neither '+' alone nor a copy of the ID line; its
 *           candidates match the ID's.
 */
inline const FieldCoderConfig kFastqFieldCoderConfig[FQ_FIELD_COUNT] = {
    /* ID      */ {{CoderType::BWT_CM, CoderType::AFFIX_MATCH}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* SEQ     */ {{CoderType::FC}, CoderType::BWT_CM, PBGZ_FQ_RANS, CoderType::RANS},
    /* QUAL    */ {{CoderType::QUAL, CoderType::FCV2, CoderType::QCM}, CoderType::QUAL, PBGZ_FQ_NONE, CoderType::QUAL},
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

/* Extra candidates a field takes on at the top compression level (-l9, the extreme-ratio path):
   coder_bwt_cm for SEQ and QUAL, which do not list it in the default archive row. */
inline std::vector<CoderType> samFieldExtremeCandidates(uint32_t fieldIdx, uint8_t mode, uint8_t level)
{
    std::vector<CoderType> extra;
    if (level < 9 || mode != PBGZ_MODE_ARCHIVE || (fieldIdx != SAM_SEQ && fieldIdx != SAM_QUAL)) {
        return extra;
    }
    const std::vector<CoderType>& base = samFieldCandidates(fieldIdx, mode);
    if (std::find(base.begin(), base.end(), CoderType::BWT_CM) == base.end()) {
        extra.push_back(CoderType::BWT_CM);
    }
    return extra;
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

/* Whether the QUAL column may use this coder: the QUAL selector trials exactly the coders
   the QUAL row lists, so removing a coder from that row takes it out of the QUAL decision. */
inline bool qualCoderCandidate(CoderType type)
{
    return samFieldCandidate(SAM_QUAL, type);
}

/* FASTQ accessors mirroring the SAM ones: candidate list, membership, and mode default. */
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

/* The FASTQ twin of samFieldExtremeCandidates: SEQ and QUAL take on coder_bwt_cm at -l9. */
inline std::vector<CoderType> fastqFieldExtremeCandidates(uint32_t fieldIdx, uint8_t mode, uint8_t level)
{
    std::vector<CoderType> extra;
    if (level < 9 || mode != PBGZ_MODE_ARCHIVE || (fieldIdx != FQ_SEQ && fieldIdx != FQ_QUAL)) {
        return extra;
    }
    const std::vector<CoderType>& base = fastqFieldCandidates(fieldIdx, mode);
    if (std::find(base.begin(), base.end(), CoderType::BWT_CM) == base.end()) {
        extra.push_back(CoderType::BWT_CM);
    }
    return extra;
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
