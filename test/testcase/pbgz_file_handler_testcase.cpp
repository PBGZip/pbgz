#include <gtest/gtest.h>
#include "pbgz_file_handler.h"

TEST(PbgzFileWriteTest, WriteNewFile) {
    PbgzFileWriter writer("pbgz_write_test.pbgz");
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    writer.close();
}


TEST(PbgzFileWriteTest, WriteFileMeta) {
    PbgzFileWriter writer("pbgz_write_file_meta_test.pbgz");
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    PbgzFileMeta fileMeta;
    std::string key = "testKey";
    Json::Value value;
    value["test"] = "value";
    fileMeta.setMetaData(key, value);
    writer.setFileMeta(fileMeta);
    writer.writeFileMeta();
    writer.close();
}

TEST(PbgzFileWriteTest, WriteBlockData) {
    PbgzFileWriter writer("pbgz_write_data_block_test.pbgz");
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    PbgzFileMeta fileMeta;
    std::string key = "fileTestKey";
    Json::Value value;
    value["test"] = "value";
    fileMeta.setMetaData(key, value);
    writer.setFileMeta(fileMeta);
    writer.writeFileMeta();

    PbgzDataBlock blockData;
    key = "metaTestKey";
    blockData.setMetaData(key, value);
    const char* testData = "This is test data.";
    uint32_t dataLength = strlen(testData);
    EXPECT_EQ(blockData.setBlockData((uint8_t*)testData, dataLength), 0);
    writer.writeBlockData(blockData);
    writer.close();
}