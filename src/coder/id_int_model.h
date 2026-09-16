/*
 * id_int_model.h - value-domain coding for an all-digit ID segment.
 *
 * The split ID path carries a segment whose text is all digits in one of two layouts today: zigzag
 * delta varints (they win when successive values are correlated - ordinals, mate pairs - because
 * then the deltas are small) and the decimal text through coder_affix_match (which wins when digit
 * prefixes repeat). Neither codes a value that is spread over a range at its own entropy: varint
 * deltas are larger than the values themselves for scattered input, and text pays a character model
 * per decimal digit, so a value the aligner picked essentially at random costs more than its
 * information content.
 *
 * This codes the value itself instead: the high bits walk an adaptive binary tree over the value
 * domain, the low bits are coded as one uniform symbol. For a value uniform over [0, max] the cost
 * lands on log2(max + 1) - the entropy - and for any other distribution the tree follows the mass
 * wherever it sits, so it needs no assumption about how the segment is ordered.
 *
 * The layout is fixed by two parameters the stream meta carries (see compressIdFieldSplit), so the
 * payload needs no header of its own:
 *   "intw" - bits per value (1..32), from the segment's maximum
 *   "intk" - how many low bits are coded raw; the remaining width-k bits walk the tree
 *
 * Encoder and decoder share this type, so the two walks cannot drift apart.
 */

#pragma once

#include <cstdint>
#include <vector>

#include "clr.h"
#include "simple_model.h"

namespace id_int {

/*
 * Bits modelled as a tree; anything wider is split into a modelled high part and raw low bits.
 * 1 << 10 nodes keeps a segment's model at a few tens of KB while still letting the tree follow any
 * distribution over the value range.
 */
static const uint32_t kMaxTreeBits = 6;

/* Bits needed to hold maxValue (at least 1). */
inline uint32_t widthFor(uint64_t maxValue)
{
    uint32_t w = 0;
    while (maxValue != 0) {
        w++;
        maxValue >>= 1;
    }
    return (w == 0) ? 1 : w;
}

/*
 * How many low bits are coded raw. The cap is not a tuning knob: the raw part travels as one symbol
 * whose frequency total is 2^lowBits, and the range coder needs its range (at least 2^24) to cover
 * that total, so past 16 bits the arithmetic would divide its range down to nothing. That in turn
 * bounds the layout to kMaxTreeBits + kMaxRawBits bits per value; anything wider is not carried by
 * this layout at all (see idIntLayoutFits).
 */
static const uint32_t kMaxRawBits = 16;

inline uint32_t lowBitsFor(uint32_t width)
{
    const uint32_t raw = (width > kMaxTreeBits) ? (width - kMaxTreeBits) : 0;
    return (raw > kMaxRawBits) ? kMaxRawBits : raw;
}

/* Whether a value of that width can travel through Model (and Counter) at all. */
inline bool widthFits(uint32_t width)
{
    return width >= 1 && width <= (kMaxTreeBits + kMaxRawBits);
}

/*
 * The uniform layout (see Model::resetUniform) spends the whole value on one symbol, so its
 * frequency total is the segment's range. The range coder needs a total it can divide its range by
 * (see kMaxRawBits), and keeping it inside 16 bits leaves the same margin the raw low bits do.
 */
static const uint64_t kMaxUniformTotal = (uint64_t)1 << 16;

class Model {
public:
    /* size_t so that a segment with fewer lines than the tree has nodes is not resized by accident. */
    void reset(uint32_t width, uint32_t lowBits)
    {
        uniform_ = false;
        width_ = width;
        lowBits_ = lowBits;
        const uint32_t treeBits = (width_ > lowBits_) ? (width_ - lowBits_) : 0;
        tree_.assign(((size_t)1 << treeBits) + 1, SIMPLE_MODEL<2>());
    }

    /*
     * Exact range: every value in [base, base + total) is equally likely and nothing is modelled.
     *
     * This is the layout for a segment that sits on its range without structure - an aligner's
     * choice of coordinate, say. The tree above can follow any distribution, but following a flat
     * one costs it (K-1)/2 * log2(N)/N bits per line in redundancy alone, which for 64 leaves over
     * a block is more than a sixth of the value's own entropy; and the modelled high bits cannot
     * help when the low bits are the ones that carry the range. Paying log2(total) per value with
     * no model at all is then simply smaller, and the encoder picks it by size like any other layout.
     */
    void resetUniform(uint64_t base, uint64_t total)
    {
        uniform_ = true;
        base_ = base;
        total_ = (total == 0) ? 1 : total;
        width_ = 0;
        lowBits_ = 0;
        tree_.clear();
    }

    bool isUniform() const { return uniform_; }

    /* Writes one value; the caller owns the RangeCoder's start/finish. */
    void encode(RangeCoder& rc, uint64_t value)
    {
        if (uniform_) {
            rc.Encode((uint32_t)(value - base_), 1, (uint32_t)total_);
            return;
        }
        const uint32_t treeBits = (width_ > lowBits_) ? (width_ - lowBits_) : 0;
        uint32_t node = 1;
        for (int32_t b = (int32_t)treeBits - 1; b >= 0; --b) {
            const uint32_t bit = (uint32_t)((value >> (uint32_t)(b + (int32_t)lowBits_)) & 1u);
            tree_[node].encodeSymbol(&rc, (uint16_t)bit);
            node = node * 2 + bit;
        }
        if (lowBits_ > 0) {
            rc.Encode((uint32_t)(value & (((uint64_t)1 << lowBits_) - 1)), 1, (uint32_t)1 << lowBits_);
        }
    }

    /* Reads one value back. */
    uint64_t decode(RangeCoder& rc)
    {
        if (uniform_) {
            const uint32_t total = (uint32_t)total_;
            const uint32_t low = rc.GetFreq(total);
            rc.Decode(low, 1);
            return base_ + low;
        }
        const uint32_t treeBits = (width_ > lowBits_) ? (width_ - lowBits_) : 0;
        uint64_t value = 0;
        uint32_t node = 1;
        for (uint32_t i = 0; i < treeBits; ++i) {
            const uint32_t bit = tree_[node].decodeSymbol(&rc);
            node = node * 2 + bit;
            value = (value << 1) | bit;
        }
        value <<= lowBits_;
        if (lowBits_ > 0) {
            const uint32_t total = (uint32_t)1 << lowBits_;
            const uint32_t low = rc.GetFreq(total);
            rc.Decode(low, 1);
            value |= low;
        }
        return value;
    }

    uint32_t width() const { return width_; }
    uint32_t lowBits() const { return lowBits_; }

private:
    bool uniform_ = false;
    uint64_t base_ = 0;
    uint64_t total_ = 0;
    uint32_t width_ = 0;
    uint32_t lowBits_ = 0;
    std::vector<SIMPLE_MODEL<2>> tree_;
};

/* Zigzag over 32 bits: the sign of a small difference in the low bit, so a step of -1 stays as
   short as a step of +1. */
inline uint32_t zigzag32(int64_t value)
{
    const int32_t v = (int32_t)value;
    return (uint32_t)((v << 1) ^ (v >> 31));
}

inline int64_t unzigzag32(uint32_t value)
{
    return (int64_t)(int32_t)((value >> 1) ^ (uint32_t)(-(int32_t)(value & 1)));
}

/*
 * The delta sequence of a counter that mostly stays put and then steps by a little - a read
 * ordinal, a flowcell row number. What has to be coded is "how long does it stay put" and "how
 * far does it step", not the values themselves.
 *
 * The varint layout spends one byte per line on the zero steps and leaves the run structure to the
 * byte model behind it, which cannot see it: measured on a coordinate-sorted BAM it costs 0.984
 * bit/line, where coding the runs costs 0.934 (runs 0.833 + steps 0.101, both of them the measured
 * entropy of the sequence). That difference is the whole point of this type.
 *
 * One bit per line says whether the run ends here; the model behind it is picked by how long the
 * run has lasted, which is what turns a fixed per-line flag (0.891 bit, worse) into the run-length
 * distribution. When the run ends, the step travels as zigzag: the small ones - which is nearly all
 * of them - through a symbol model of their own, and anything past that through a value-domain tree
 * (see Model) that also absorbs the rare backwards jump. The split matters: a single stray step of
 * -19363 would otherwise set the tree's width to 16 bits and make every ordinary +1 pay ten raw
 * bits of it as well.
 */
class Counter {
public:
    static const uint32_t kRunCtx = 16;
    /* Steps 1..kSmallStep are coded straight; kSmallStep + 1 is the escape to the tree. */
    static const uint32_t kSmallStep = 7;

    void reset(uint32_t escapeWidth, uint32_t escapeLowBits)
    {
        runCtx_.assign(kRunCtx, SIMPLE_MODEL<2>());
        stepSmall_.reset();
        stepBig_.reset(escapeWidth, escapeLowBits);
        run_ = 0;
    }

    void encode(RangeCoder& rc, int64_t delta)
    {
        const uint32_t ctx = (run_ < kRunCtx) ? run_ : (kRunCtx - 1);
        if (delta == 0) {
            runCtx_[ctx].encodeSymbol(&rc, 0);
            run_++;
            return;
        }
        runCtx_[ctx].encodeSymbol(&rc, 1);
        const uint32_t zz = zigzag32(delta);
        if (zz >= 1 && zz <= kSmallStep) {
            stepSmall_.encodeSymbol(&rc, (uint16_t)(zz - 1));
        } else {
            stepSmall_.encodeSymbol(&rc, (uint16_t)kSmallStep);
            stepBig_.encode(rc, zz);
        }
        run_ = 0;
    }

    int64_t decode(RangeCoder& rc)
    {
        const uint32_t ctx = (run_ < kRunCtx) ? run_ : (kRunCtx - 1);
        if (runCtx_[ctx].decodeSymbol(&rc) == 0) {
            run_++;
            return 0;
        }
        run_ = 0;
        const uint16_t sym = stepSmall_.decodeSymbol(&rc);
        const uint32_t zz = (sym < kSmallStep) ? (uint32_t)(sym + 1) : (uint32_t)stepBig_.decode(rc);
        return unzigzag32(zz);
    }

private:
    std::vector<SIMPLE_MODEL<2>> runCtx_;
    SIMPLE_MODEL<kSmallStep + 1> stepSmall_;
    Model stepBig_;
    uint32_t run_ = 0;
};

}  // namespace id_int
