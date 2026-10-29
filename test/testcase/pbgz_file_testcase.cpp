/*
 * pbgz_file_testcase.cpp - Test case for PBGZ file operations
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
#include "pbgz_file.h"
#include <coder_json.h>
#include "farmhash/src/farmhash.h"

#include <vector>

// Test cases for PbgzFileHeader
TEST(PbgzFileHeader, PbgzFileHeaderInit) {
    PbgzFileHeader header;
    EXPECT_EQ(header.getBlockType(), FILE_HEADER);
    EXPECT_EQ(header.getVersionStr(), "2.2.0");
}

// Test cases for PbgzFileMeta
TEST(PbgzFileMeta, PbgzFileMetaInit) {
    PbgzFileMeta meta;
    EXPECT_EQ(meta.getBlockType(), FILE_META);
    EXPECT_TRUE(meta.getMetaData().empty());
}   

TEST(PbgzFileMeta, PbgzFileMetaSetKey) {
    PbgzFileMeta meta;
    std::string key = "testKey";
    std::string value = "testValue";
    meta.setMetaData(key, value);
    EXPECT_EQ(meta.getMetaData(key).asString(), value);
}

TEST(PbgzDataBlock, PbgzDataBlockSetMetaData) {
    PbgzDataBlock dataBlock;
    std::string key = "testKey";
    Json::Value value;
    value["test"] = "value";
    dataBlock.setMetaData(key, value);
    EXPECT_EQ(dataBlock.getMetaData(key), value);
}

TEST(PbgzDataBlock, PbgzDataBlockSetData) {
    PbgzDataBlock dataBlock;
    const char* testData = "This is test data.";
    uint32_t dataLength = strlen(testData);
    EXPECT_EQ(dataBlock.setBlockData((uint8_t*)testData, dataLength), 0);
    EXPECT_EQ(dataBlock.getDataLength(), dataLength);
    EXPECT_NE(dataBlock.getDataPtr(), nullptr);
    EXPECT_EQ(memcmp(dataBlock.getDataPtr(), testData, dataLength), 0);
}


/*
 * Container checksums switched from util::Hash64() to the portable
 * Fingerprint64 so files written on x86 and ARM can be verified on either
 * platform.  Files written by older versions still carry util::Hash64() and
 * must keep verifying on the platform that wrote them.
 */
TEST(PbgzDataBlock, PortableChecksumAndLegacyFile) {
    std::vector<uint8_t> data(256);
    for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (uint8_t)(i * 29u + 7u);
    }

    PbgzDataBlock block;
    ASSERT_EQ(block.setBlockData(data.data(), (uint32_t)data.size()), 0);
    block.setMetaData("field", "value");
    block.calcChecksum();
    EXPECT_EQ(block.verifyCheckSum(true), 0);

    Json::StreamWriterBuilder writer;
    const std::string meta = Json::writeString(writer, block.getMetaData());

    /* Legacy files stored util::Hash64(), which is CPU-specific and is no
       longer accepted.  Such files skip the container check (their content
       still carries MD5), so they can cross architectures. */
    block.setDataCheckSum(util::Hash64((const char*)data.data(), data.size()));
    block.setMetaCheckSum(util::Hash64(meta));
    EXPECT_EQ(block.verifyCheckSum(false), 0);

    /* A corrupted portable checksum must still fail. */
    block.setDataCheckSum(util::Fingerprint64((const char*)data.data(), data.size()) ^ 1ULL);
    block.setMetaCheckSum(util::Fingerprint64(meta));
    EXPECT_NE(block.verifyCheckSum(true), 0);
}
