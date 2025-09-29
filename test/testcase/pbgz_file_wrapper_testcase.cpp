#include <gtest/gtest.h>
#include "pbgz_file_wrapper.h"

TEST(PbgzFileWriteTest, WriteNewFile) {
    IOWriter* iowriter = new FileWriter("pbgz_write_test.pbgz");
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    writer.close();
    delete iowriter;
}


TEST(PbgzFileWriteTest, WriteFileMeta) {
    IOWriter* iowriter = new FileWriter("pbgz_write_file_meta_test.pbgz");
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    PbgzFileMeta fileMeta;
    std::string key = "testKey";
    Json::Value value;
    value["test"] = "value";
    fileMeta.setMetaData(key, value);
    writer.setFileMeta(fileMeta);
    writer.writeFileMeta();
    writer.close();
    delete iowriter;
}

TEST(PbgzFileWriteTest, WriteBlockData) {
     IOWriter* iowriter = new FileWriter("pbgz_write_data_block_test.pbgz");
    PbgzFileWriter writer(iowriter);
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
    delete iowriter;
    
}

TEST(PbgzFileReadTest, ReadFileMeta) {
    IOReader* ioreader = new FileReader("pbgz_write_file_meta_test.pbgz");
    PbgzFileReader reader(ioreader);
    EXPECT_EQ(reader.open(), 0) << "Failed to open file for reading";
    PbgzFileHeader fileHeader = reader.getFileHeader();
    EXPECT_EQ(fileHeader.getBlockType(), FILE_HEADER) << "Invalid file magic";
    PbgzFileMeta fileMeta = reader.getFileMeta();
    EXPECT_EQ(fileMeta.getBlockType(), FILE_META) << "Failed to read file metadata";
    Json::Value value = fileMeta.getMetaData("testKey");
    ASSERT_FALSE(value.isNull()) << "Metadata 'testKey' not found";
    ASSERT_EQ(value["test"].asString(), "value") << "Metadata 'testKey' has incorrect value";
    reader.close();
    delete ioreader;
}

TEST(PbgzFileReadTest, ReadBlockData) {
    IOReader* ioreader = new FileReader("pbgz_write_data_block_test.pbgz");
    PbgzFileReader reader(ioreader);
    EXPECT_EQ(reader.open(), 0) << "Failed to open file for reading";
    PbgzFileHeader fileHeader = reader.getFileHeader();
    EXPECT_EQ(fileHeader.getBlockType(), FILE_HEADER) << "Invalid file magic";
    PbgzFileMeta fileMeta = reader.getFileMeta();
    EXPECT_EQ(fileMeta.getBlockType(), FILE_META) << "Failed to read file metadata";
    Json::Value value = fileMeta.getMetaData("fileTestKey");
    ASSERT_FALSE(value.isNull()) << "Metadata 'fileTestKey' not found";
    ASSERT_EQ(value["test"].asString(), "value") << "Metadata 'fileTestKey' has incorrect value";

    PbgzDataBlock dataBlock;
    EXPECT_EQ(reader.readDataBlock(dataBlock), 0) << "Failed to read data block";
    EXPECT_EQ(dataBlock.getBlockType(), FILE_DATA) << "Invalid block type";
    value = dataBlock.getMetaData("metaTestKey");
    ASSERT_FALSE(value.isNull()) << "Metadata 'metaTestKey' not found";
    EXPECT_EQ(value["test"].asString(), "value") << "Metadata 'metaTestKey' has incorrect value";
    const uint8_t* dataPtr = dataBlock.getDataPtr();
    ASSERT_NE(dataPtr, nullptr) << "Data pointer is null";
    std::string dataStr((const char*)dataPtr, dataBlock.getDataLength());
    EXPECT_EQ(dataStr, "This is test data.") << "Data content is incorrect";
    reader.close();
    delete ioreader;
}