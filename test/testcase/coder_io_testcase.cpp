/*
 * coder_io_testcase.cpp - Cpp file for pbgz project
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

#include <cstdint>
#include <vector>

#include "coder/coder_io.h"

/*
 * Behavior of the out-of-bounds error sink.
 *
 * These cases guard one property: "any stream that overflows during a single block process must
 * be queryable at the exit point." Previously this relied on every call site checking its own err;
 * the result was that SAM checked all 12 streams while FASTQ and the index checked none, so data
 * corrupted by an overflow was written out as if it had succeeded.
 */
class CoderIoErrSinkTest : public ::testing::Test {};

/* Without a sink attached, behavior is unchanged: only the local err is set and nothing else is affected. */
TEST_F(CoderIoErrSinkTest, StandaloneIoKeepsErrorLocal)
{
    std::vector<uint8_t> buf(2);
    coder_io io(buf.data(), (int32_t)buf.size());

    io.putc('a');
    io.putc('b');
    EXPECT_EQ(io.err, coder_io::IO_OK);

    io.putc('c');
    EXPECT_EQ(io.err, coder_io::IO_BUF_FULL);
}

/* A write overflow is reported to the sink immediately, together with which stream overflowed. */
TEST_F(CoderIoErrSinkTest, WriteOverflowReachesSink)
{
    coder_err_sink sink;
    std::vector<uint8_t> buf(1);
    coder_io io(buf.data(), (int32_t)buf.size(), &sink, "QUAL");

    io.putc('a');
    EXPECT_TRUE(sink.ok());

    io.putc('b');
    EXPECT_FALSE(sink.ok());
    EXPECT_EQ(sink.err, coder_io::IO_BUF_FULL);
    EXPECT_STREQ(sink.what, "QUAL");
}

/* Read exhaustion is reported as well. '\0' itself is valid data, so a flag alone distinguishes a real 0 from nothing left to read. */
TEST_F(CoderIoErrSinkTest, ReadExhaustionReachesSink)
{
    coder_err_sink sink;
    std::vector<uint8_t> buf{0x00};
    coder_io io(buf.data(), (int32_t)buf.size(), &sink, "SEQ");

    EXPECT_EQ(io.getc(), 0x00);
    EXPECT_TRUE(sink.ok()) << "A real read 0 must not be treated as exhaustion";

    EXPECT_EQ(io.getc(), 0x00);
    EXPECT_FALSE(sink.ok());
    EXPECT_EQ(sink.err, coder_io::IO_READ_EMPTY);
    EXPECT_STREQ(sink.what, "SEQ");
}

/*
 * A block opens over a dozen streams; the sink keeps only the first error: subsequent errors are
 * usually chain reactions of it, and the first one has the most diagnostic value.
 */
TEST_F(CoderIoErrSinkTest, SinkKeepsFirstErrorAcrossStreams)
{
    coder_err_sink sink;
    std::vector<uint8_t> bufA(1);
    std::vector<uint8_t> bufB(1);

    coder_io first(bufA.data(), (int32_t)bufA.size(), &sink, "QNAME");
    coder_io second(bufB.data(), (int32_t)bufB.size(), &sink, "QUAL");

    first.putc('a');
    first.putc('b');
    ASSERT_EQ(sink.err, coder_io::IO_BUF_FULL);
    ASSERT_STREQ(sink.what, "QNAME");

    second.putc('a');
    second.getc();
    second.getc();
    EXPECT_EQ(sink.err, coder_io::IO_BUF_FULL) << "The first error must be preserved";
    EXPECT_STREQ(sink.what, "QNAME");
}

/* err is sticky: repeated overflows on the same stream never overwrite the recorded error. */
TEST_F(CoderIoErrSinkTest, ErrorIsSticky)
{
    coder_err_sink sink;
    std::vector<uint8_t> buf(1);
    coder_io io(buf.data(), (int32_t)buf.size(), &sink, "CIGAR");

    io.getc();
    io.getc();
    ASSERT_EQ(io.err, coder_io::IO_READ_EMPTY);

    io.putc('x');
    EXPECT_EQ(io.err, coder_io::IO_READ_EMPTY);
    EXPECT_EQ(sink.err, coder_io::IO_READ_EMPTY);
}
