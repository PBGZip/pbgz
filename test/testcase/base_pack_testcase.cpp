/*
 * base_pack_testcase.cpp - SIMD vs scalar packing equivalence
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

#include <gtest/gtest.h>

#include <cstring>
#include <random>
#include <vector>

#include "coder/base_pack.h"

namespace {

std::vector<uint8_t> randomBases(uint32_t len, uint32_t seed)
{
    static const char kBases[] = "ACGT";
    std::mt19937 rng(seed);
    std::vector<uint8_t> in(len);
    for (uint32_t i = 0; i < len; ++i) {
        in[i] = (uint8_t)kBases[rng() & 3u];
    }
    return in;
}

}  // namespace

/* The SIMD path (SSE on x86, NEON on ARM) must produce exactly the scalar
   bytes for every length, including the tails. */
TEST(BasePackTest, MatchesScalar) {
    for (uint32_t len = 0; len <= 600; ++len) {
        const std::vector<uint8_t> in = randomBases(len, 20261028u + len);
        const size_t cap = (len + 3) / 4 + 4;
        std::vector<uint8_t> packed(cap, 0xAA);
        std::vector<uint8_t> scalarOut(cap, 0xAA);

        const bool simdOk = pbgz_base_pack::packBases(in.data(), len, packed.data());
        const bool scalarOk = pbgz_base_pack::scalar(in.data(), len, scalarOut.data());

        ASSERT_EQ(simdOk, scalarOk) << "len=" << len;
        ASSERT_TRUE(simdOk) << "len=" << len;
        ASSERT_EQ(0, memcmp(packed.data(), scalarOut.data(), (len + 3) / 4))
            << "len=" << len;
    }
}

TEST(BasePackTest, RejectsNonAcgt) {
    const uint8_t bad[] = {'A', 'C', 'G', 'T', 'N'};
    std::vector<uint8_t> out(4, 0);
    EXPECT_FALSE(pbgz_base_pack::packBases(bad, sizeof(bad), out.data()));
    EXPECT_FALSE(pbgz_base_pack::scalar(bad, sizeof(bad), out.data()));

    const uint8_t lower[] = {'A', 'C', 'G', 't'};
    EXPECT_FALSE(pbgz_base_pack::packBases(lower, sizeof(lower), out.data()));
    EXPECT_FALSE(pbgz_base_pack::scalar(lower, sizeof(lower), out.data()));
}
