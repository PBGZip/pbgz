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
 * 越界错误汇聚点的行为。
 *
 * 这些用例守的是一条性质："一次块处理里任何一路流越界，都必须能在出口被问到"。
 * 原来它靠每个调用点自觉去查各自的 err，结果 SAM 查了 12 路、FASTQ 和索引一路
 * 都没查，越界写坏的数据被当成成功写了出去。
 */
class CoderIoErrSinkTest : public ::testing::Test {};

/* 不接汇聚点时行为与原来一致：只置本地 err，不影响别人 */
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

/* 写越界立刻上报汇聚点，并带上是哪一路流 */
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

/* 读端耗尽同样上报。'\0' 本身是合法数据，只能靠标志区分"真读到 0"和"没得读" */
TEST_F(CoderIoErrSinkTest, ReadExhaustionReachesSink)
{
    coder_err_sink sink;
    std::vector<uint8_t> buf{0x00};
    coder_io io(buf.data(), (int32_t)buf.size(), &sink, "SEQ");

    EXPECT_EQ(io.getc(), 0x00);
    EXPECT_TRUE(sink.ok()) << "读到真实的 0 不该被当成耗尽";

    EXPECT_EQ(io.getc(), 0x00);
    EXPECT_FALSE(sink.ok());
    EXPECT_EQ(sink.err, coder_io::IO_READ_EMPTY);
    EXPECT_STREQ(sink.what, "SEQ");
}

/*
 * 一个块要开十几路流，汇聚点只留第一个错误：后续错误多半是它的连锁反应，
 * 第一个最有定位价值。
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
    EXPECT_EQ(sink.err, coder_io::IO_BUF_FULL) << "首个错误应保持不变";
    EXPECT_STREQ(sink.what, "QNAME");
}

/* err 是粘性的：同一路流反复越界不会把错误刷成别的 */
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
