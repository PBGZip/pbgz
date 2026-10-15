/*
 * sam_field_rules.h - rules a SAM field's value follows, shared by the actuators that read the
 * same field from different representations (SamCodecActuator's text vs BamCodecActuator's BAM
 * columns), so the two cannot drift apart.
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

/* The template length a pair implies, in one of two aligner conventions differing by 1: bwa
   subtracts one (minusOne = true), minimap2 does not. The caller chooses the convention per
   block and decides whether the pair has a template length at all. */
static inline int32_t samTemplateLen(int64_t pos, int64_t nextPos, uint32_t selfSpan,
                                     uint32_t mateSpan, bool minusOne)
{
    const int64_t conv = minusOne ? 1 : 0;
    if (pos < nextPos) {
        return (int32_t)(nextPos + (int64_t)mateSpan - pos - conv);
    }
    if (pos > nextPos) {
        return (int32_t)-(pos + (int64_t)selfSpan - nextPos - conv);
    }
    /* The two ends coincide: the mate's side of the template decides, as in the positive branch. */
    return (int32_t)(nextPos + (int64_t)mateSpan - pos - conv);
}
