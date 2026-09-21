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
 * list is empty because its dedicated selector keeps using the archive row.
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
