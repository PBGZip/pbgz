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

#include "preprocess_info.h"

/*
 * Per-field encoder configuration table for SAM: candidate coders (trial-compressed
 * and compared per field during preprocessing) and the default coder (the actuator's
 * fallback when nothing is selected / selection failed / selection was not wired up).
 * Changing a field's candidates or default touches only this table; both the
 * preprocessing trial scope and the actuator defaults are read from here, so there is
 * no need to edit two places.
 *
 * Special cases:
 *  - QUAL(10): goes through the QualSelector-specific path (coder_qual / fcv2 /
 *    bwt_cm), does not participate in generic trial compression; only a default is
 *    registered here.
 *  - POS(3) / PNEXT(7) / TLEN(8): always use delta / inferred compression (see
 *    compressPosFieldDelta / compressPNextFieldDelta / compressTLen); no generic
 *    candidates, the default is the underlying bwt_cm.
 */

struct FieldCoderConfig {
    std::vector<CoderType> candidates; /* candidate coders, trial-compressed in order */
    CoderType fallback;                /* default coder */
};

/* Number of SAM fields participating in coder selection: 11 mandatory fields + the 12th OPTION column. */
static const uint32_t SAM_FIELD_COUNT_SELECT = SAM_FIELD_COUNT + 1;

inline const FieldCoderConfig kSamFieldCoderConfig[SAM_FIELD_COUNT_SELECT] = {
    /* QNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* FLAG  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* RNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* POS   */ {{}, CoderType::BWT_CM},
    /* MAPQ  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* CIGAR */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* RNEXT */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* PNEXT */ {{}, CoderType::BWT_CM},
    /* TLEN  */ {{}, CoderType::BWT_CM},
    /* SEQ   */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::FC},
    /* QUAL  */ {{}, CoderType::QUAL},
    /* OPTION */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
};

inline const FieldCoderConfig* samFieldCoderConfig(uint32_t fieldIdx)
{
    if (fieldIdx >= SAM_FIELD_COUNT_SELECT) {
        return nullptr;
    }
    return &kSamFieldCoderConfig[fieldIdx];
}

/* Whether a field lists the given coder as a candidate. */
inline bool samFieldCandidate(uint32_t fieldIdx, CoderType type)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return false;
    }
    return std::find(cfg->candidates.begin(), cfg->candidates.end(), type) != cfg->candidates.end();
}

/* Default coder; returns the caller-supplied fallback when nothing is registered. */
inline CoderType samFieldDefaultCoder(uint32_t fieldIdx, CoderType fallback)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    return (cfg != nullptr) ? cfg->fallback : fallback;
}
