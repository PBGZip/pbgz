/*
 * sam_field_rules.h - the rules a SAM field's value follows, shared by the actuators that read the
 * same field from a different representation.
 *
 * SamCodecActuator parses SAM text and BamCodecActuator consumes the columns of a BAM record, but
 * the values they reconstruct are the same values, so the rules behind them have to be the same
 * too. The template length of a pair is one of them: TLEN is encoded against it, and both actuators
 * used to spell the rule out separately, which is exactly how two copies drift apart.
 */

#pragma once

#include <cstdint>

/*
 * The template length a pair implies, in one of the two conventions aligners write. They differ by
 * exactly 1:
 *   bwa       right-end position + right-read length - left-end position - 1   (minusOne = true)
 *   minimap2  right-end position + right-read length - left-end position       (minusOne = false)
 * Neither is "standard", so the convention is chosen per block by counting which one matches more
 * records; the wrong one marks 99%+ of records as exceptions and degenerates the exception stream
 * into storing full values.
 *
 * The caller decides whether the pair has a template length at all - unpaired, either end unmapped,
 * or a mate on another reference gives 0 - and what the two spans are (0 is allowed, for a mate
 * whose span is not in this block).
 */
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
