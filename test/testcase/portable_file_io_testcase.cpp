/*
 * portable_file_io_testcase.cpp - Tests for the streaming file I/O backend
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

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "portable_file_io.h"

namespace {

const char* kPortableIoFile = "ut_portable_io.bin";

std::string readWholeFile(BufferedFileReader& reader) {
    std::string out;
    char buf[64];
    size_t n = 0;
    while ((n = reader.readIO(buf, sizeof(buf))) > 0) {
        out.append(buf, n);
    }
    return out;
}

}  // namespace

TEST(PortableFileIoTest, ReadWriteRoundTrip) {
    std::remove(kPortableIoFile);
    const std::string payload = "portable streaming io payload";

    BufferedFileWriter writer(kPortableIoFile);
    ASSERT_EQ(writer.openIO(), 0);
    EXPECT_EQ(writer.writeIO(payload.data(), payload.size()), payload.size());
    writer.closeIO();
    EXPECT_EQ(writer.getWriteError(), 0);

    BufferedFileReader reader(kPortableIoFile);
    ASSERT_EQ(reader.openIO(), 0);
    EXPECT_EQ(reader.getFileSize(), payload.size());
    EXPECT_EQ(readWholeFile(reader), payload);
    reader.closeIO();

    std::remove(kPortableIoFile);
}

TEST(PortableFileIoTest, SeekAndReadLine) {
    std::remove(kPortableIoFile);
    const std::string text = "alpha\r\nbeta\n\ngamma";

    BufferedFileWriter writer(kPortableIoFile);
    ASSERT_EQ(writer.openIO(), 0);
    ASSERT_EQ(writer.writeIO(text.data(), text.size()), text.size());
    writer.closeIO();

    BufferedFileReader reader(kPortableIoFile);
    ASSERT_EQ(reader.openIO(), 0);

    std::string line;
    EXPECT_EQ(reader.readLine(line), 5u);
    EXPECT_EQ(line, "alpha");
    EXPECT_EQ(reader.readLine(line), 4u);
    EXPECT_EQ(line, "beta");
    EXPECT_EQ(reader.readLine(line), 0u);
    EXPECT_EQ(line, "");
    EXPECT_EQ(reader.readLine(line), 5u);
    EXPECT_EQ(line, "gamma");
    EXPECT_TRUE(reader.isEOF());

    EXPECT_EQ(reader.seekIO(0), 0);
    EXPECT_EQ(reader.getCurrentPos(), 0u);
    EXPECT_FALSE(reader.isEOF());
    EXPECT_EQ(reader.readIO(nullptr, 0), 0u);

    reader.closeIO();
    std::remove(kPortableIoFile);
}

TEST(PortableFileIoTest, WriteAtThenAppend) {
    std::remove(kPortableIoFile);
    const std::string tail = " 0123456789";

    BufferedFileWriter writer(kPortableIoFile);
    ASSERT_EQ(writer.openIO(), 0);
    ASSERT_EQ(writer.writeIO(tail.data(), tail.size()), tail.size());

    const uint64_t patch = 0x1122334455667788ull;
    EXPECT_GT(writer.writeIOAt(1, &patch, sizeof(patch)), 0);
    ASSERT_EQ(writer.seekToEnd(), 0);
    ASSERT_EQ(writer.writeIO("END", 3), 3u);
    writer.closeIO();
    EXPECT_EQ(writer.getWriteError(), 0);

    BufferedFileReader reader(kPortableIoFile);
    ASSERT_EQ(reader.openIO(), 0);
    std::string content = readWholeFile(reader);
    reader.closeIO();

    ASSERT_GE(content.size(), sizeof(patch) + 1 + 3);
    EXPECT_EQ(content[0], ' ');
    EXPECT_EQ(memcmp(content.data() + 1, &patch, sizeof(patch)), 0);
    EXPECT_EQ(content.substr(content.size() - 3), "END");
    EXPECT_EQ(reader.getFileSize(), content.size());

    std::remove(kPortableIoFile);
}

TEST(PortableFileIoTest, LargeSequentialIo) {
    std::remove(kPortableIoFile);
    const size_t total = 8u << 20;
    std::vector<uint8_t> payload(total);
    for (size_t i = 0; i < total; ++i) {
        payload[i] = (uint8_t)(i * 131u + 17u);
    }

    BufferedFileWriter writer(kPortableIoFile);
    ASSERT_EQ(writer.openIO(), 0);
    size_t written = 0;
    while (written < total) {
        const size_t chunk = (total - written > (1u << 20)) ? (1u << 20) : (total - written);
        ASSERT_EQ(writer.writeIO(payload.data() + written, chunk), chunk);
        written += chunk;
    }
    writer.closeIO();

    BufferedFileReader reader(kPortableIoFile);
    ASSERT_EQ(reader.openIO(), 0);
    EXPECT_EQ(reader.getFileSize(), total);
    std::vector<uint8_t> got(total);
    size_t read = 0;
    while (read < total) {
        const size_t chunk = (total - read > (1u << 20)) ? (1u << 20) : (total - read);
        ASSERT_EQ(reader.readIO(got.data() + read, chunk), chunk);
        read += chunk;
    }
    reader.closeIO();
    EXPECT_EQ(memcmp(payload.data(), got.data(), total), 0);

    std::remove(kPortableIoFile);
}
