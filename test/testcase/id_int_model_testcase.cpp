/*
 * id_int_model_testcase.cpp - value-domain coding of an all-digit ID segment (see id_int_model.h).
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * The layout exists to code a value spread over a range at its own entropy instead of paying for
 * its decimal digits (see the file header of id_int_model.h). Both halves are pinned here: the
 * round trip for the shapes a QNAME segment actually takes - scattered coordinates, a constant, a
 * tiny range, values that need the raw low bits - and the cost, which is the reason to have it at
 * all: a uniform range must land near log2(range) per value, well below what the digits cost.
 */

#include <gtest/gtest.h>

#include <cstdint>
#include <vector>

#include "coder/id_int_model.h"

namespace {

/* Values that look like a flowcell coordinate: spread over [0, 30000], nothing sorted. */
std::vector<uint64_t> scatteredValues(size_t count, uint64_t range = 30001)
{
    std::vector<uint64_t> values;
    values.reserve(count);
    uint64_t seed = 0x9E3779B97F4A7C15ull;
    for (size_t i = 0; i < count; ++i) {
        seed = seed * 6364136223846793005ull + 1442695040888963407ull;
        values.push_back((seed >> 33) % range);
    }
    return values;
}

std::vector<uint8_t> encodeAll(const std::vector<uint64_t>& values, uint32_t width, uint32_t lowBits)
{
    std::vector<uint8_t> out(values.size() * 8 + 64);
    RangeCoder rc;
    rc.output((char*)out.data(), (char*)out.data() + out.size());
    rc.StartEncode();
    id_int::Model model;
    model.reset(width, lowBits);
    for (size_t i = 0; i < values.size(); ++i) {
        model.encode(rc, values[i]);
    }
    rc.FinishEncode();
    EXPECT_EQ(rc.err, 0);
    const int n = rc.size_out();
    out.resize((n > 0) ? (size_t)n : 0);
    return out;
}

std::vector<uint64_t> decodeAll(const std::vector<uint8_t>& bytes, uint32_t width, uint32_t lowBits,
                                size_t count)
{
    std::vector<uint8_t> buf = bytes;
    RangeCoder rc;
    rc.input((char*)buf.data(), (char*)buf.data() + buf.size());
    rc.StartDecode();
    EXPECT_EQ(rc.err, 0);
    id_int::Model model;
    model.reset(width, lowBits);
    std::vector<uint64_t> out;
    out.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        out.push_back(model.decode(rc));
    }
    return out;
}

/* The layout the encoder picks for a segment: width from the maximum, raw low bits beyond the tree. */
void layoutFor(uint64_t maxValue, uint32_t& width, uint32_t& lowBits)
{
    width = id_int::widthFor(maxValue);
    lowBits = id_int::lowBitsFor(width);
}

}  // namespace

TEST(IdIntModelTest, RoundTripsScatteredValues) {
    const std::vector<uint64_t> values = scatteredValues(5000);
    uint32_t width = 0, lowBits = 0;
    layoutFor(30000, width, lowBits);
    /* 30000 needs 15 bits; how many of them stay raw is what kMaxTreeBits decides. */
    EXPECT_EQ(width, 15u);
    EXPECT_EQ(lowBits, id_int::lowBitsFor(width));
    EXPECT_LT(lowBits, width);
    EXPECT_GT(width - lowBits, 0u);

    const std::vector<uint8_t> bytes = encodeAll(values, width, lowBits);
    EXPECT_EQ(decodeAll(bytes, width, lowBits, values.size()), values);
}

TEST(IdIntModelTest, RoundTripsConstantTinyAndWideValues) {
    /* A constant segment (a tile that never changes in the block). */
    const std::vector<uint64_t> constant(500, 11102);
    uint32_t w = 0, k = 0;
    layoutFor(11102, w, k);
    EXPECT_EQ(decodeAll(encodeAll(constant, w, k), w, k, constant.size()), constant);

    /* A tiny range: every value fits in the tree, so nothing is left raw. */
    std::vector<uint64_t> tiny;
    for (int i = 0; i < 500; ++i) {
        tiny.push_back((uint64_t)(i % 7));
    }
    layoutFor(6, w, k);
    EXPECT_EQ(w, 3u);
    EXPECT_EQ(k, 0u);
    EXPECT_EQ(decodeAll(encodeAll(tiny, w, k), w, k, tiny.size()), tiny);

    /* Wide values: the low bits travel raw, up to the width the layout can carry (see
       kMaxRawBits) - beyond that the encoder refuses the layout rather than losing precision. */
    std::vector<uint64_t> wide = {0, 1, 0x3FFFFFull, 0x200000ull, 0x123456ull};
    layoutFor(0x3FFFFFull, w, k);
    EXPECT_EQ(w, 22u);                       /* the widest a value can be and still fit the layout */
    EXPECT_EQ(k, id_int::lowBitsFor(w));
    EXPECT_TRUE(id_int::widthFits(w));
    EXPECT_EQ(decodeAll(encodeAll(wide, w, k), w, k, wide.size()), wide);

    /* And the bound itself: a 32-bit value cannot be carried, so the encoder must not try. */
    layoutFor(0xFFFFFFFFull, w, k);
    EXPECT_EQ(w, 32u);
    EXPECT_FALSE(id_int::widthFits(w));
}

/*
 * The reason the layout exists: a value uniform over [0, max] must cost about log2(max + 1) bits.
 * Over 30001 values that is 14.87 bits each, while the decimal text - what the segment cost before -
 * needs a character model per digit and lands around 15.8. The margin asserted here is wide enough
 * to survive the models' start-up cost over 200k values.
 */
TEST(IdIntModelTest, CostsAboutTheEntropyOfAUniformRange) {
    const std::vector<uint64_t> values = scatteredValues(200000);
    uint32_t width = 0, lowBits = 0;
    layoutFor(30000, width, lowBits);
    const std::vector<uint8_t> bytes = encodeAll(values, width, lowBits);

    const double bitsPerValue = 8.0 * (double)bytes.size() / (double)values.size();
    const double entropy = 14.87;   /* log2(30001) */
    EXPECT_LT(bitsPerValue, 15.6);
    EXPECT_GT(bitsPerValue, entropy - 0.2);
}

/*
 * The counter layout that the run-length mode rests on: a sequence that stays put and steps, with
 * the occasional backwards jump (the case that would be lost if the steps were coded at one width).
 */
TEST(IdIntModelTest, CounterRoundTripsRunsAndSteps) {
    std::vector<int64_t> deltas;
    uint64_t seed = 0xC0FFEEull;
    for (int i = 0; i < 5000; ++i) {
        seed = seed * 6364136223846793005ull + 1442695040888963407ull;
        const uint64_t r = seed >> 33;
        if (r % 10 < 7) {
            deltas.push_back(0);                 /* the run continues */
        } else {
            deltas.push_back((int64_t)(r % 3) + 1);   /* a small step */
        }
    }
    deltas.push_back(-19363);                     /* the stray backwards jump */
    deltas.push_back(1);

    uint64_t maxEscape = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
        const uint64_t zz = id_int::zigzag32(deltas[i]);
        if (zz > id_int::Counter::kSmallStep && zz > maxEscape) {
            maxEscape = zz;
        }
    }
    const uint32_t width = id_int::widthFor(maxEscape);
    const uint32_t lowBits = id_int::lowBitsFor(width);

    std::vector<uint8_t> bytes(deltas.size() * 6 + 4096);
    RangeCoder rc;
    rc.output((char*)bytes.data(), (char*)bytes.data() + bytes.size());
    rc.StartEncode();
    id_int::Counter enc;
    enc.reset(width, lowBits);
    for (size_t i = 0; i < deltas.size(); ++i) {
        enc.encode(rc, deltas[i]);
    }
    rc.FinishEncode();
    EXPECT_EQ(rc.err, 0);
    bytes.resize((size_t)rc.size_out());

    std::vector<uint8_t> in = bytes;
    RangeCoder drc;
    drc.input((char*)in.data(), (char*)in.data() + in.size());
    drc.StartDecode();
    id_int::Counter dec;
    dec.reset(width, lowBits);
    for (size_t i = 0; i < deltas.size(); ++i) {
        EXPECT_EQ(dec.decode(drc), deltas[i]) << "step " << i;
    }
}
