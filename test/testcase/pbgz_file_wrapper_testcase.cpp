#include <gtest/gtest.h>
#include "pbgz_file_wrapper.h"

TEST(PbgzFileWriteTest, WriteNewFile) {
    std::string fileName = "pbgz_write_test.pbgz";
    (void)remove(fileName.c_str());
    IOWriter* iowriter = new FileWriter(fileName);
    iowriter->openIO();
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    writer.close();
    iowriter->closeIO();
    delete iowriter;
}


TEST(PbgzFileWriteTest, WriteFileMeta) {
    std::string fileName = "pbgz_write_file_meta_test.pbgz";
    (void)remove(fileName.c_str());
    IOWriter* iowriter = new FileWriter(fileName);
    iowriter->openIO();
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    PbgzFileMeta fileMeta;
    std::string key = "testKey";
    Json::Value value;
    value["test"] = "value";
    value["test2"] = "value";
    value["test3"] = "value";
    fileMeta.setMetaData(key, value);
    writer.setFileMeta(fileMeta);
    writer.writeFileMeta();
    writer.close();
    iowriter->closeIO();
    delete iowriter;
}

TEST(PbgzFileWriteTest, WriteBlockData) {
    std::string fileName = "pbgz_write_data_block_test.pbgz";
    (void)remove(fileName.c_str());
    IOWriter* iowriter = new FileWriter(fileName);
    iowriter->openIO();
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    PbgzFileMeta fileMeta;
    std::string key = "fileTestKey";
    Json::Value value;
    value["test"] = "value";
    value["test2"] = "value";
    value["test3"] = "value";
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
    iowriter->closeIO();
    delete iowriter;
    
}

TEST(PbgzFileReadTest, ReadFileMeta) {
    std::string fileName = "pbgz_write_file_meta_test.pbgz";
    IOReader* ioreader = new FileReader(fileName);
    ioreader->openIO();
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
    ioreader->closeIO();
    delete ioreader;
}

TEST(PbgzFileReadTest, ReadBlockData) {
    std::string fileName = "pbgz_write_data_block_test.pbgz";
    IOReader* ioreader = new FileReader(fileName);
    ioreader->openIO();
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
    ioreader->closeIO();
    delete ioreader;
}