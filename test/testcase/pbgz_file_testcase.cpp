#include <gtest/gtest.h>
#include "pbgz_file.h"
#include "coder.h"
#include "utils/memory_util.h"
#include <coder_json.h>

// Test cases for PbgzFileHeader
TEST(PbgzFileHeader, PbgzFileHeaderInit) {
    PbgzFileHeader header;
    EXPECT_EQ(header.getBlockType(), FILE_HEADER);
    EXPECT_EQ(header.getVersionStr(), "2.0.0");
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

