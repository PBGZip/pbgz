/*
 * id_int_model.h - value-domain coding for an all-digit ID segment.
 *
 * Codes the value itself rather than varint deltas or decimal text: the high
 * bits walk an adaptive binary tree over the value domain, the low bits are
 * coded as one uniform symbol.
 *
 * Layout is fixed by the stream meta: "intw" (bits per value, 1..32) and "intk"
 * (number of low bits coded raw).
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

#include "clr.h"
#include "simple_model.h"

namespace id_int {

/* Bits modelled as a tree; wider values split into modelled high and raw low bits. */
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

/* Max raw low bits: the raw part is one symbol whose total the range coder must cover. */
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

/* Max total for the uniform layout: the range coder must be able to divide by it. */
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

    /* Uniform layout: every value in [base, base + total) is equally likely and nothing is modelled. */
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
 * Codes a mostly-constant counter's delta sequence: one bit per line for
 * "run continues", then the step via zigzag (small steps through a symbol model,
 * larger ones through the value-domain tree).
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
