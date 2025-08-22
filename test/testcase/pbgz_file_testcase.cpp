#include <gtest/gtest.h>
#include "pbgz_file.h"

// Test cases for PbgzFileHeader
TEST(PbgzFileHeader, PbgzFileHeaderInit) {
    PbgzFileHeader header;
    EXPECT_EQ(header.getBlockType(), FILE_HEADER);
    EXPECT_EQ(header.getVersionStr(), "2.0.0");
}

TEST(PbgzFileHeader, PbgzFileHeaderSerialize) {
    PbgzFileHeader header;
    uint8_t buffer[7];
    uint32_t dataLength = 0;    
    header.serialize(buffer, sizeof(buffer), dataLength);
    EXPECT_EQ(dataLength, 7);
    EXPECT_EQ(memcmp(buffer, PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH), 0);
    char version[3] = {2, 0, 0};
    EXPECT_EQ(memcmp(buffer + PBGZ_FILE_MAGIC_LENGTH, version, 3), 0);   
}

TEST(PbgzFileHeader, PbgzFileHeaderUnserialize) {
    PbgzFileHeader header;
    uint8_t buffer[7];
    uint32_t dataLength = 0;    
    header.serialize(buffer, sizeof(buffer), dataLength);
    EXPECT_EQ(dataLength, 7);
    
    PbgzFileHeader newHeader;
    EXPECT_EQ(newHeader.unserialize(buffer, dataLength), 0);
    EXPECT_EQ(newHeader.getBlockType(), FILE_HEADER);
    EXPECT_EQ(newHeader.getVersionStr(), "2.0.0");
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

TEST(PbgzFileMeta, PbgzFileMetaSerialize) {
    PbgzFileMeta meta;
    std::string key = "testKey";
    std::string value = "testValue";
    meta.setMetaData(key, value);
    uint8_t buffer[1024];
    uint32_t dataLength = 0;
    EXPECT_EQ(meta.serialize(buffer, sizeof(buffer), dataLength), 0);
    EXPECT_GT(dataLength, 0);
    EXPECT_EQ(memcmp(buffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH), 0);

    uint32_t metaLength;
    memcpy(&metaLength, buffer + PBGZ_FILE_META_MAGIC_LENGTH, PBGZ_FILE_META_SIZE_LENGTH);
    EXPECT_GT(metaLength, 0);   
    std::string jsonString(reinterpret_cast<char*>(buffer + PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH), metaLength);
    Json::CharReaderBuilder readerBuilder;
    Json::Value jsonData;   
    std::istringstream jsonStream(jsonString);
    std::string errs;
    EXPECT_TRUE(Json::parseFromStream(readerBuilder, jsonStream, &jsonData, &errs));
    EXPECT_EQ(jsonData[key].asString(), value); 
    uint64_t checksum = *(uint64_t*)(buffer + PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH + metaLength);
    EXPECT_EQ(checksum, 0);    
}

TEST(PbgzFileMeta, PbgzFileMetaUnserialize) {
    PbgzFileMeta meta;
    std::string key = "testKey";    
    std::string value = "testValue";
    meta.setMetaData(key, value);
    uint8_t buffer[1024];
    uint32_t dataLength = 0;
    EXPECT_EQ(meta.serialize(buffer, sizeof(buffer), dataLength), 0);
    EXPECT_GT(dataLength, 0);
    PbgzFileMeta newMeta;
    EXPECT_EQ(newMeta.unserialize(buffer, dataLength), 0);
    EXPECT_EQ(newMeta.getBlockType(), FILE_META);
    EXPECT_EQ(newMeta.getMetaData(key).asString(), value);
    EXPECT_EQ(newMeta.getMetaChecksum(), 0); // Assuming checksum is not set in this test
}

TEST(PbgzDataBlock, PbgzDataBlockInit) {
    PbgzDataBlock dataBlock;
    EXPECT_EQ(dataBlock.getBlockType(), FILE_DATA);
    EXPECT_TRUE(dataBlock.getMetaData().empty());
    EXPECT_EQ(dataBlock.getDataLength(), 0);
    EXPECT_EQ(dataBlock.getDataPtr(), nullptr);
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

TEST(PbgzDataBlock, PbgzDataBlockSerialize) {
    PbgzDataBlock dataBlock;
    std::string key = "testKey";        
    Json::Value value;
    value["test"] = "value";
    dataBlock.setMetaData(key, value);
    const char* testData = "This is test data.";
    uint32_t dataLength = strlen(testData);
    EXPECT_EQ(dataBlock.setBlockData((uint8_t*)testData, dataLength), 0);   

    uint8_t buffer[1024];
    uint32_t dataLengthOut = 0;
    EXPECT_EQ(dataBlock.serialize(buffer, sizeof(buffer), dataLengthOut), 0);
    EXPECT_GT(dataLengthOut, 0);
    EXPECT_EQ(memcmp(buffer, &PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH), 0);

    uint32_t metaLength;
    memcpy(&metaLength, buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH, PBGZ_DATA_BLOCK_META_SIZE_LENGTH);
    EXPECT_GT(metaLength, 0);
    std::string jsonString(reinterpret_cast<char*>(buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH), metaLength);
    Json::CharReaderBuilder readerBuilder;
    Json::Value jsonData;
    std::istringstream jsonStream(jsonString);
    std::string errs;
    EXPECT_TRUE(Json::parseFromStream(readerBuilder, jsonStream, &jsonData, &errs));
    EXPECT_EQ(jsonData[key]["test"].asString(), value["test"].asString());
    uint64_t metaChecksum = (*(uint64_t*)(buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + metaLength));
    EXPECT_EQ(metaChecksum, 0); // Assuming checksum is not set in this test

    uint32_t dataLengthRead;
    memcpy(&dataLengthRead, buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + metaLength + PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH, 
          sizeof(dataLengthRead));
    EXPECT_EQ(dataLengthRead, dataLength);

    EXPECT_EQ(memcmp(buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + metaLength + PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH + PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH,
                   dataBlock.getDataPtr(), dataLength), 0);

    uint64_t dataBlockChecksum = *(uint64_t*)(buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + metaLength + PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH + PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH + dataLength);
    EXPECT_EQ(dataBlockChecksum, 0); // Assuming checksum is not set in this test
    uint64_t originDataChecksum = *(uint64_t*)(buffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + metaLength + PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH + PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH + dataLength + PBGZ_DATA_BLOCK_CHECKSUM_LENGTH);
    EXPECT_EQ(originDataChecksum, 0); // Assuming origin checksum is not set in this test
}

TEST(PbgzDataBlock, PbgzDataBlockUnserialize) {
    PbgzDataBlock dataBlock;
    std::string key = "testKey";
    Json::Value value;
    value["test"] = "value";
    dataBlock.setMetaData(key, value);
    const char* testData = "This is test data.";
    uint32_t dataLength = strlen(testData);
    EXPECT_EQ(dataBlock.setBlockData((uint8_t*)testData, dataLength), 0);
    uint8_t buffer[1024];
    uint32_t dataLengthOut = 0;
    EXPECT_EQ(dataBlock.serialize(buffer, sizeof(buffer), dataLengthOut), 0);
    EXPECT_GT(dataLengthOut, 0);
    PbgzDataBlock newDataBlock;
    EXPECT_EQ(newDataBlock.unserialize(buffer, dataLengthOut), 0);
    EXPECT_EQ(newDataBlock.getBlockType(), FILE_DATA);
    EXPECT_EQ(newDataBlock.getMetaData(key)["test"].asString(), value["test"].asString());
    EXPECT_EQ(newDataBlock.getDataLength(), dataLength);
    EXPECT_NE(newDataBlock.getDataPtr(), nullptr);
    EXPECT_EQ(memcmp(newDataBlock.getDataPtr(), testData, dataLength), 0);
    EXPECT_EQ(newDataBlock.getOriginDataChecksum(), 0); // Assuming origin checksum is not set in this test
}