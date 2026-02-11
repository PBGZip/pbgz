#include <gtest/gtest.h>
#include <coder.h>

#include "pbgz_file_wrapper.h"
#include "utils/memory_util.h"
#include "utils/path_util.h"

class PbgzFileWriteTest : public ::testing::Test {
public:
    void SetUp() override { 
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    }

    void TearDown() override { 

    }
};

class PbgzFileReadTest : public ::testing::Test {
public:
    void SetUp() override { 
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    }

    void TearDown() override { 

    }
};

TEST_F(PbgzFileWriteTest, WriteNewFile) {
    std::string fileName = "pbgz_write_test.pbgz";
    (void)remove(fileName.c_str());
    IOWriter* iowriter = new FileWriter(fileName);
    iowriter->openIO();
    PbgzFileWriter writer(iowriter);
    ASSERT_EQ(writer.open(), 0) << "Failed to open file for writing";
    writer.close();
    iowriter->closeIO();
    delete iowriter;
    
    // Clean up test generated files
    PathUtil::removeFile(fileName);
}

TEST_F(PbgzFileWriteTest, WriteFileMeta) {
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
    writer.setBaseFileMeta(fileMeta);
    writer.writeBaseFileMeta();
    writer.close();
    iowriter->closeIO();
    delete iowriter;
    
    // Clean up test generated files
    PathUtil::removeFile(fileName);
}

TEST_F(PbgzFileWriteTest, WriteBlockData) {
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
    writer.setBaseFileMeta(fileMeta);
    writer.writeBaseFileMeta();
    
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
    
    // Clean up test generated files
    PathUtil::removeFile(fileName);
}

TEST_F(PbgzFileReadTest, ReadFileMeta) {
    std::string fileName = "pbgz_write_file_meta_test.pbgz";
    
    // First create test file
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
    writer.setBaseFileMeta(fileMeta);
    writer.writeBaseFileMeta();
    writer.close();
    iowriter->closeIO();
    delete iowriter;
    
    // Now read file
    IOReader* ioreader = new FileReader(fileName);
    ioreader->openIO();
    PbgzFileReader reader(ioreader);
    EXPECT_EQ(reader.open(), 0) << "Failed to open file for reading";
    PbgzFileHeader fileHeader = reader.getFileHeader();
    EXPECT_EQ(fileHeader.getBlockType(), FILE_HEADER) << "Invalid file magic";
    PbgzFileMeta fileMetaRead = reader.getBaseFileMeta();
    EXPECT_EQ(fileMetaRead.getBlockType(), FILE_META) << "Failed to read file metadata";
    Json::Value valueRead = fileMetaRead.getMetaData("testKey");
    ASSERT_FALSE(valueRead.isNull()) << "Metadata 'testKey' not found";
    ASSERT_EQ(valueRead["test"].asString(), "value") << "Metadata 'testKey' has incorrect value";
    reader.close();
    ioreader->closeIO();
    delete ioreader;
    
    // Clean up test generated files
    PathUtil::removeFile(fileName);
}

TEST_F(PbgzFileReadTest, ReadBlockData) {
    std::string fileName = "pbgz_write_data_block_test.pbgz";
    
    // First create test file
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
    writer.setBaseFileMeta(fileMeta);
    writer.writeBaseFileMeta();
    
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
    
    // Now read file
    IOReader* ioreader = new FileReader(fileName);
    ioreader->openIO();
    PbgzFileReader reader(ioreader);
    EXPECT_EQ(reader.open(), 0) << "Failed to open file for reading";
    PbgzFileHeader fileHeader = reader.getFileHeader();
    EXPECT_EQ(fileHeader.getBlockType(), FILE_HEADER) << "Invalid file magic";
    PbgzFileMeta fileMetaRead = reader.getBaseFileMeta();
    EXPECT_EQ(fileMetaRead.getBlockType(), FILE_META) << "Failed to read file metadata";
    Json::Value valueRead = fileMetaRead.getMetaData("fileTestKey");
    ASSERT_FALSE(valueRead.isNull()) << "Metadata 'fileTestKey' not found";
    ASSERT_EQ(valueRead["test"].asString(), "value") << "Metadata 'fileTestKey' has incorrect value";

    PbgzDataBlock dataBlock;
    uint8_t readerBuff[2048];
    dataBlock.setDataPtr(readerBuff);
    EXPECT_EQ(reader.readDataBlock(dataBlock), 0) << "Failed to read data block";
    EXPECT_EQ(dataBlock.getBlockType(), FILE_DATA) << "Invalid block type";
    valueRead = dataBlock.getMetaData("metaTestKey");
    ASSERT_FALSE(valueRead.isNull()) << "Metadata 'metaTestKey' not found";
    EXPECT_EQ(valueRead["test"].asString(), "value") << "Metadata 'metaTestKey' has incorrect value";
    const uint8_t* dataPtr = dataBlock.getDataPtr();
    ASSERT_NE(dataPtr, nullptr) << "Data pointer is null";
    std::string dataStr((const char*)dataPtr, dataBlock.getDataLength());
    EXPECT_EQ(dataStr, "This is test data.") << "Data content is incorrect";
    reader.close();
    ioreader->closeIO();
    delete ioreader;
    
    // Clean up test generated files
    PathUtil::removeFile(fileName);
}