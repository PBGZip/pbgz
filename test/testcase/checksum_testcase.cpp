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

TEST(ChecksumTest, FarmHash128StreamMatchesApi) {
    std::vector<uint8_t> data(4096);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 37u + 11u);
    }

    /* 0..600 covers the seed/144-byte/short-tail boundaries of CityHash128. */
    for (size_t len = 0; len <= 600; ++len) {
        const std::string direct = checksumHex(
            ChecksumType::FarmHash128, data.data(), len);
        for (size_t chunk : {size_t(1), size_t(2), size_t(7), size_t(63),
                             size_t(64), size_t(65), size_t(127), size_t(128),
                             size_t(129), size_t(255)}) {
            ChecksumStream stream(ChecksumType::FarmHash128, len);
            for (size_t off = 0; off < len; off += chunk) {
                const size_t n = (off + chunk <= len) ? chunk : (len - off);
                stream.update(data.data() + off, n);
            }
            EXPECT_EQ(stream.finish(), direct) << "len=" << len << " chunk=" << chunk;
        }
    }

    const std::string direct =
        checksumHex(ChecksumType::FarmHash128, data.data(), data.size());
    ChecksumStream stream(ChecksumType::FarmHash128, data.size());
    for (size_t off = 0; off < data.size(); off += 13) {
        const size_t n = (off + 13 <= data.size()) ? 13 : (data.size() - off);
        stream.update(data.data() + off, n);
    }
    EXPECT_EQ(stream.finish(), direct);
}

TEST(ChecksumTest, FarmHash128StreamLargeInput) {
    std::vector<uint8_t> data(8u << 20);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 29u + (i >> 7));
    }
    const std::string direct = checksumHex(
        ChecksumType::FarmHash128, data.data(), data.size());

    for (size_t chunk : {size_t(63), size_t(64), size_t(65), size_t(4096),
                         size_t(1u << 20), data.size()}) {
        ChecksumStream stream(ChecksumType::FarmHash128, data.size());
        for (size_t off = 0; off < data.size(); off += chunk) {
            const size_t n = (off + chunk <= data.size()) ? chunk : (data.size() - off);
            stream.update(data.data() + off, n);
        }
        EXPECT_EQ(stream.finish(), direct) << "chunk=" << chunk;
    }
}

TEST(ChecksumTest, FarmHash128ShortFeedIsSafe) {
    std::vector<uint8_t> data(256);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 17u + 3u);
    }
    for (size_t fed : {size_t(0), size_t(1), size_t(15), size_t(16), size_t(100),
                       size_t(119), size_t(120), size_t(127)}) {
        for (size_t declared : {fed, fed + 100, size_t(200)}) {
            ChecksumStream stream(ChecksumType::FarmHash128, declared);
            stream.update(data.data(), fed);
            EXPECT_EQ(stream.finish(),
                      checksumHex(ChecksumType::FarmHash128, data.data(), fed))
                << "fed=" << fed << " declared=" << declared;
        }
    }
}

TEST(ChecksumTest, FarmHash128MatchesFarmHashApi) {
    const std::string data = "pbgz farmhash fingerprint128";
    auto value = util::Fingerprint128(data.data(), data.size());
    char expected[33];
    snprintf(expected, sizeof(expected), "%016llx%016llx",
             (unsigned long long)util::Uint128High64(value),
             (unsigned long long)util::Uint128Low64(value));
    const std::string actual = checksumHex(ChecksumType::FarmHash128, data);
    EXPECT_EQ(actual, expected);
    EXPECT_EQ(actual.size(), 32u);
}

TEST(ChecksumTest, MetadataMarkerCompatibility) {
    Json::Value meta;

    setChecksumMeta(meta, ChecksumType::FarmHash128, "00112233445566778899aabbccddeeff");
    ChecksumType type = ChecksumType::MD5;
    EXPECT_EQ(getChecksumMeta(meta, type), "00112233445566778899aabbccddeeff");
    EXPECT_EQ(type, ChecksumType::FarmHash128);
    EXPECT_EQ(meta[kChecksumMethodKey].asString(), "f");
    EXPECT_EQ(meta["hash"].asString(), "00112233445566778899aabbccddeeff");

    /* A file produced by an older writer: only md5, no marker. */
    Json::Value old;
    old["md5"] = "900150983cd24fb0d6963f7d28e17f72";
    type = ChecksumType::FarmHash128;
    EXPECT_EQ(getChecksumMeta(old, type), "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(type, ChecksumType::MD5);

    /* Legacy marker as written by the previous release. */
    Json::Value marked;
    marked["hashmethod"] = "md5";
    marked["md5"] = "900150983cd24fb0d6963f7d28e17f72";
    type = ChecksumType::FarmHash128;
    EXPECT_EQ(getChecksumMeta(marked, type), "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(type, ChecksumType::MD5);

    /* MD5 writes the short tag too. */
    Json::Value md5Meta;
    setChecksumMeta(md5Meta, ChecksumType::MD5, "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(md5Meta[kChecksumMethodKey].asString(), "m");
    type = ChecksumType::FarmHash128;
    EXPECT_EQ(getChecksumMeta(md5Meta, type), "900150983cd24fb0d6963f7d28e17f72");
    EXPECT_EQ(type, ChecksumType::MD5);
}

TEST(ChecksumTest, ShortTagNames) {
    EXPECT_EQ(checksumTypeFromName("f"), ChecksumType::FarmHash128);
    EXPECT_EQ(checksumTypeFromName("m"), ChecksumType::MD5);
    EXPECT_EQ(checksumTypeFromName("md5"), ChecksumType::MD5);
    EXPECT_EQ(checksumTypeFromName(""), ChecksumType::MD5);
}
