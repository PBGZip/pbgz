/*
 * qname_dedup_coder_testcase.cpp - the QNAME dedup layout's flag/distance stream
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
#include <vector>

#include "coder/coder_io.h"
#include "coder/coder_qname_dedup.h"
#include "utils/memory_util.h"

class QNameDedupCoderTest : public ::testing::Test
{
public:
    static const int BUFFER_SIZE = 1 << 20;

    void SetUp() override
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        packed_ = (uint8_t*)malloc(BUFFER_SIZE);
        roundtrip_ = (uint8_t*)malloc(BUFFER_SIZE);
        ASSERT_NE(packed_, nullptr);
        ASSERT_NE(roundtrip_, nullptr);
    }

    void TearDown() override
    {
        free(packed_);
        free(roundtrip_);
    }

    /* flags + three distance planes, as encodeQnameDedup builds it */
    static std::vector<uint8_t> makePayload(uint32_t lines, uint32_t reps)
    {
        std::vector<uint8_t> p((size_t)lines + (size_t)reps * 3, 0);
        uint32_t written = 0;
        for (uint32_t i = 0; i < lines; ++i) {
            if (written < reps && (i % 3) == 2) {
                p[i] = 1;
                const uint32_t dist = 1 + (i * 7919u) % (i + 1);
                p[lines + written] = (uint8_t)dist;
                p[lines + reps + written] = (uint8_t)(dist >> 8);
                p[lines + 2 * reps + written] = (uint8_t)(dist >> 16);
                ++written;
            }
        }
        return p;
    }

    int32_t encode(const std::vector<uint8_t>& payload, uint32_t lines, uint32_t reps)
    {
        coder_io io(packed_, BUFFER_SIZE);
        io.meta["lines"] = (Json::Value::UInt)lines;
        io.meta["reps"] = (Json::Value::UInt)reps;
        coder_qname_dedup enc(&io);
        enc.encode_line(payload.data(), (uint32_t)payload.size());
        enc.encode_flush();
        return io.data_len;
    }

    uint8_t* packed_ = nullptr;
    uint8_t* roundtrip_ = nullptr;
};

TEST_F(QNameDedupCoderTest, RoundtripsFlagsAndDistancePlanes)
{
    const uint32_t lines = 1000;
    const uint32_t reps = 300;
    const std::vector<uint8_t> payload = makePayload(lines, reps);
    const int32_t packed = encode(payload, lines, reps);
    ASSERT_GT(packed, 0);

    coder_io io(packed_, packed);
    io.meta["lines"] = (Json::Value::UInt)lines;
    io.meta["reps"] = (Json::Value::UInt)reps;
    coder_qname_dedup dec(&io);
    ASSERT_EQ(dec.decode_line(roundtrip_, (uint32_t)payload.size(), UINT8_MAX, false),
              (int32_t)payload.size());
    EXPECT_EQ(memcmp(roundtrip_, payload.data(), payload.size()), 0);

    /* A skewed payload must be smaller than its raw form. */
    EXPECT_LT((uint32_t)packed, (uint32_t)payload.size());
}

TEST_F(QNameDedupCoderTest, HandlesANoRepeatBlock)
{
    const uint32_t lines = 4096;
    const uint32_t reps = 0;
    const std::vector<uint8_t> payload = makePayload(lines, reps);
    const int32_t packed = encode(payload, lines, reps);
    ASSERT_GT(packed, 0);

    coder_io io(packed_, packed);
    io.meta["lines"] = (Json::Value::UInt)lines;
    io.meta["reps"] = (Json::Value::UInt)reps;
    coder_qname_dedup dec(&io);
    ASSERT_EQ(dec.decode_line(roundtrip_, (uint32_t)payload.size(), UINT8_MAX, false),
              (int32_t)payload.size());
    EXPECT_EQ(memcmp(roundtrip_, payload.data(), payload.size()), 0);
}

TEST_F(QNameDedupCoderTest, RejectsAMismatchedShape)
{
    const uint32_t lines = 64;
    const uint32_t reps = 8;
    const std::vector<uint8_t> payload = makePayload(lines, reps);
    const int32_t packed = encode(payload, lines, reps);
    ASSERT_GT(packed, 0);

    coder_io io(packed_, packed);
    io.meta["lines"] = (Json::Value::UInt)lines;
    io.meta["reps"] = (Json::Value::UInt)reps;
    coder_qname_dedup dec(&io);
    /* one byte short of the advertised shape */
    EXPECT_LT(dec.decode_line(roundtrip_, (uint32_t)payload.size() - 1, UINT8_MAX, false), 0);
}
