/*
 * qname_coder_testcase.cpp - Test cases for the coder_qname coder
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
#include <cstdlib>
#include <string>
#include <vector>
#include <memory>

#include "coder/coder.h"
#include "coder/coder_io.h"
#include "coder/coder_qname.h"
#include "utils/memory_util.h"

class QNameCoderTest : public ::testing::Test
{
public:
    static const int BUFFER_SIZE = 1 << 24;

    void SetUp() override
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        compressed_data = (uint8_t*) malloc(BUFFER_SIZE);
        decompressed_data = (uint8_t*) malloc(BUFFER_SIZE);
        ASSERT_NE(compressed_data, nullptr);
        ASSERT_NE(decompressed_data, nullptr);
    }

    void TearDown() override
    {
        free(compressed_data);
        free(decompressed_data);
    }

    uint8_t* compressed_data = nullptr;
    uint8_t* decompressed_data = nullptr;

    /* Encode all QNAMEs (including the trailing '\t') and return the encoder's output length. */
    int32_t encodeAll(std::vector<std::string>& lines)
    {
        coder_io io(compressed_data, BUFFER_SIZE);
        coder_qname enc(&io);
        for (auto& s : lines) {
            enc.encode_line((const uint8_t*)s.data(), (uint32_t)s.size());
        }
        enc.encode_flush();
        return io.data_len;
    }

    bool roundtrip(std::vector<std::string>& lines)
    {
        int32_t packed = encodeAll(lines);
        coder_io io(compressed_data, packed);
        coder_qname dec(&io);
        for (size_t i = 0; i < lines.size(); i++) {
            int32_t len = dec.decode_line(decompressed_data, BUFFER_SIZE, UINT8_MAX, false);
            if (len < 0 || (size_t)len != lines[i].size() ||
                memcmp(decompressed_data, lines[i].data(), len) != 0) {
                return false;
            }
        }
        return true;
    }
};

/* Realistic QNAMEs: a fixed 9-char prefix (5 values) + '.' + a 7~8 digit number, repeated in pairs but never adjacent. */
static std::vector<std::string> makeRealisticNames(int pairs, unsigned seed)
{
    const char* prefixes[5] = {"ERR013136", "ERR015528", "ERR156633", "ERR162873", "ERR162876"};
    std::vector<std::string> names;
    names.reserve(pairs * 2);
    srand(seed);
    for (int i = 0; i < pairs; i++) {
        int p = rand() % 5;
        int id = 1000000 + rand() % 20000000;
        std::string n = std::string(prefixes[p]) + "." + std::to_string(id) + "\t";
        names.push_back(n);
        /* R1/R2 of the same fragment share a name but are separated by about 1~60 lines after coordinate sorting. */
        int gap = 1 + rand() % 60;
        for (int g = 0; g < gap && (int)names.size() < pairs * 2; g++) {
            int p2 = rand() % 5;
            int id2 = 1000000 + rand() % 20000000;
            names.push_back(std::string(prefixes[p2]) + "." + std::to_string(id2) + "\t");
        }
        if ((int)names.size() < pairs * 2) {
            names.push_back(n);
        }
    }
    return names;
}

TEST_F(QNameCoderTest, BasicRoundTrip)
{
    std::vector<std::string> lines = {"ERR015528.21801860\t", "ERR162873.22213189\t",
        "ERR013136.12036163\t", "ERR013136.10734940\t", "ERR013136.4208740\t"};
    EXPECT_TRUE(roundtrip(lines));
}

/* Round-trip consistency across many records (including cross-line duplicates) and a compression ratio below 50% of the raw size. */
TEST_F(QNameCoderTest, RealisticRoundTripAndRatio)
{
    std::vector<std::string> lines = makeRealisticNames(4000, 7);
    ASSERT_GT(lines.size(), 1000u);
    EXPECT_TRUE(roundtrip(lines));

    size_t raw = 0;
    for (auto& s : lines) raw += s.size();
    int32_t packed = encodeAll(lines);
    EXPECT_LT((size_t)packed, raw / 2) << "Fixed prefix + paired-redundant data must compress below 50%";
}

/* Widely separated duplicates beyond the small window: the content must still be lossless (a copy miss degrades to NEW). */
TEST_F(QNameCoderTest, FarApartDuplicateStillLossless)
{
    /* The line distance must exceed RING_SIZE (4096); generate enough lines first, then insert a duplicate pair. */
    std::vector<std::string> lines = makeRealisticNames(6000, 3);
    ASSERT_GT(lines.size(), 5000u);
    /* Create a duplicate pair whose distance far exceeds RING_SIZE. */
    std::string first = lines[0];
    lines.insert(lines.begin() + 5000, first);
    EXPECT_TRUE(roundtrip(lines));
}

/* Many unique names (no copy possible): round-trip is consistent and no infinite loop occurs. */
TEST_F(QNameCoderTest, ManyUniqueNames)
{
    std::vector<std::string> lines;
    lines.reserve(200000);
    srand(99);
    for (int i = 0; i < 200000; i++) {
        int p = rand() % 5;
        int id = 1 + rand() % 40000000;
        lines.push_back(std::string("ERR") + std::to_string(100000 + p) + "." + std::to_string(id) + "\t");
    }
    EXPECT_TRUE(roundtrip(lines));
}
