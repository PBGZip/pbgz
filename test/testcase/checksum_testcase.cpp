/*
 * checksum_testcase.cpp - tests for checksum algorithms and metadata markers
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

#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "farmhash/src/farmhash.h"
#include "utils/checksum.h"

TEST(ChecksumTest, Md5Compatibility) {
    EXPECT_EQ(checksumHex(ChecksumType::MD5, "abc"),
              "900150983cd24fb0d6963f7d28e17f72");
}

TEST(ChecksumTest, FarmHashStreamMatchesFingerprint64) {
    std::vector<uint8_t> data(4096);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 37u + 11u);
    }

    for (size_t len = 0; len <= 300; ++len) {
        const std::string direct = checksumHex(
            ChecksumType::FarmHash64, data.data(), len);
        for (size_t chunk : {size_t(1), size_t(2), size_t(7), size_t(63),
                             size_t(64), size_t(65), size_t(127)}) {
            ChecksumStream stream(ChecksumType::FarmHash64);
            for (size_t off = 0; off < len; off += chunk) {
                const size_t n = (off + chunk <= len) ? chunk : (len - off);
                stream.update(data.data() + off, n);
            }
            EXPECT_EQ(stream.finish(), direct) << "len=" << len << " chunk=" << chunk;
        }
    }

    std::string direct = checksumHex(ChecksumType::FarmHash64, data.data(), data.size());
    ChecksumStream stream(ChecksumType::FarmHash64);
    for (size_t off = 0; off < data.size(); off += 13) {
        const size_t n = (off + 13 <= data.size()) ? 13 : (data.size() - off);
        stream.update(data.data() + off, n);
    }
    EXPECT_EQ(stream.finish(), direct);
}

TEST(ChecksumTest, FarmHashStreamLargeInput) {
    std::vector<uint8_t> data(8u << 20);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 29u + (i >> 7));
    }
    const std::string direct = checksumHex(
        ChecksumType::FarmHash64, data.data(), data.size());

    for (size_t chunk : {size_t(1), size_t(63), size_t(64), size_t(65),
                         size_t(4096), size_t(1u << 20), data.size()}) {
        ChecksumStream stream(ChecksumType::FarmHash64);
        for (size_t off = 0; off < data.size(); off += chunk) {
            const size_t n = (off + chunk <= data.size()) ? chunk : (data.size() - off);
            stream.update(data.data() + off, n);
        }
        EXPECT_EQ(stream.finish(), direct) << "chunk=" << chunk;
    }
}

TEST(ChecksumTest, FarmHashMatchesFarmHashApi) {
    const std::string data = "pbgz farmhash fingerprint64";
    char expected[17];
    snprintf(expected, sizeof(expected), "%016llx",
             (unsigned long long)util::Fingerprint64(data.data(), data.size()));
    EXPECT_EQ(checksumHex(ChecksumType::FarmHash64, data), expected);
}

TEST(ChecksumTest, MetadataMarkerCompatibility) {
    Json::Value meta;

    setChecksumMeta(meta, ChecksumType::FarmHash64, "0011223344556677");
    ChecksumType type = ChecksumType::MD5;
    EXPECT_EQ(getChecksumMeta(meta, type), "0011223344556677");
    EXPECT_EQ(type, ChecksumType::FarmHash64);
    EXPECT_EQ(meta["hashmethod"].asString(), "farmhash64");

    /* A file produced by an older writer: only md5, no marker. */
    Json::Value old;
    old["md5"] = "900150983cd24fb0d6963f7d28e17f72";
    type = ChecksumType::FarmHash64;
    EXPECT_EQ(getChecksumMeta(old, type), "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(type, ChecksumType::MD5);

    /* Explicit md5 marker. */
    Json::Value marked;
    marked["hashmethod"] = "md5";
    marked["md5"] = "900150983cd24fb0d6963f7d28e17f72";
    type = ChecksumType::FarmHash64;
    EXPECT_EQ(getChecksumMeta(marked, type), "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(type, ChecksumType::MD5);
}
