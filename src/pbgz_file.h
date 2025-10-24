/*
 * pbgz_file.h - Header file for pbgz project
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

#pragma once

#include <string>
#include <json/json.h>

#include "pbgz_types.h"

/// @brief Magic value for block in PBGZ file format, used to identify block type
typedef enum {
    INVALID = 0,  // unknown
    FILE_HEADER =  1, // File header information
    FILE_META = 2,    // File meta information
    FILE_DATA = 3,    // File data
}PbgzBlockType;

const std::string PBGZ_FILE_MAGIC = "PBGZ"; // Magic value for PBGZ file format
const uint32_t PBGZ_FILE_META_MAGIC = 0x000000EB; // Magic value for PBGZ file meta information
const uint32_t PBGZ_DATA_BLOCK_MAGIC = 0x000000DB; // Magic value for PBGZ file meta information


/// @brief File type in PBGZ file format
typedef enum {
    BINARY_FILE = 0, // Binary file
    FASTQ_FILE = 1,  // FASTQ format file
} FileFormat;

const uint32_t PBGZ_FILE_MAGIC_LENGTH = 4; // Length of PBGZ file magic value
const uint32_t PBGZ_FILE_VERSION_LENGTH = 3; // Length of file version number

const uint32_t PBGZ_FILE_META_MAGIC_LENGTH = 4; // Length of PBGZ file meta magic value
const uint32_t PBGZ_FILE_META_SIZE_LENGTH = 4; // Length of PBGZ file meta information length
const uint32_t PBGZ_FILE_META_CHECKSUM_LENGTH = 8; // Length of PBGZ data block magic value

const uint32_t PBGZ_DATA_BLOCK_MAGIC_LENGTH = 4; // Length of PBGZ data block magic value
const uint32_t PBGZ_DATA_BLOCK_META_SIZE_LENGTH = 4; // Length of PBGZ data block meta information length
const uint32_t PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH = 4;  // Length of PBGZ data block data length
const uint32_t PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH = 8; // Length of PBGZ data block meta information checksum
const uint32_t PBGZ_DATA_BLOCK_CHECKSUM_LENGTH = 8; // Length of PBGZ data block checksum


class Serializable {
public:
    /**
     * @brief Serialize the object to a buffer.
     * @param buffer Destination buffer.
     * @param bufferLength Length of the buffer.
     * @param dataLength Actual length of serialized data.
     * @return Status code.
     */
    virtual int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) = 0;

    /**
     * @brief Deserialize the object from a buffer.
     * @param buffer Source buffer.
     * @param bufferLength Length of the buffer.
     * @return Status code.
     */
    virtual int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) = 0;
};

/// @brief PBGZ file header
class PbgzFileHeader {
public:
    // /**
    //  * @brief Serialize the file header to a buffer.
    //  */
    // int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override ;

    // /**
    //  * @brief Deserialize the file header from a buffer.
    //  */
    // int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    /**
     * @brief Get the version string of the file header.
     */
    std::string getVersionStr();

    /**
     * @brief Get the block type of the file header.
     */
    PbgzBlockType getBlockType() const {
        return blockType;
    }

    void setBlockType(PbgzBlockType type) {
        blockType = type;
    }

    uint8_t* getVerion() {
        return version;
    }

    void setVersion(uint8_t* buffer, uint32_t size) {
        if (size != PBGZ_FILE_VERSION_LENGTH) {
            return;
        }
        memcpy(version, buffer, size);
    }

    PbgzFileHeader() : blockType(FILE_HEADER) {
        version[0] = PBGZ_VERSION_MAJOR;
        version[1] = PBGZ_VERSION_MINOR;
        version[2] = PBGZ_VERSION_PATCH;
    }

private:
    PbgzBlockType blockType;  // Block type
    uint8_t version[3];  // Version number
};

/// @brief PBGZ file meta information
class PbgzFileMeta {
public:
    // /**
    //  * @brief Serialize the file meta information to a buffer.
    //  */
    // int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override ;

    // /**
    //  * @brief Unserialize the file meta information from a buffer.
    //  */
    // int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    /**
     * @brief Check if the meta data checksum is valid.
     * @return True if checksum is valid, false otherwise.
     */
    uint64_t getMetaChecksum() const{
        return metaChecksum;
    }

    void setMetaChecksum(uint64_t checksum) {
        metaChecksum = checksum;
    }

    Json::Value& getMetaData() {
        return metaData;
    }

    Json::Value getMetaData(const std::string& key) const {
        return metaData[key];
    }

    void setMetaData(const std::string &key, const Json::Value& value) {
        metaData[key] = value;
    }

    PbgzBlockType getBlockType() const {
        return blockType;   
    }

    void setBlockType(PbgzBlockType type) {
        blockType = type;
    }
 
    PbgzFileMeta():blockType(FILE_META) {
        metaData.clear();
        metaChecksum = 0;
    }
    
    virtual ~PbgzFileMeta() {
        metaData.clear();
    }

private:
    PbgzBlockType blockType;   // Block type
    Json::Value metaData; // Meta information in JSON format
    uint64_t metaChecksum; // Checksum for meta data
};


/// @brief  Pbgz file block data
class PbgzDataBlock {
public:
    PbgzBlockType getBlockType() {
        return blockType;
    }

    void setBlockType(PbgzBlockType type) {
        blockType = type;
    }

    Json::Value& getMetaData() {
        return dataMetaInfo;
    }

    Json::Value& getMetaData(const std::string &key) {
        return dataMetaInfo[key];
    }

    void setMetaData(const Json::Value& value) {
        dataMetaInfo = value;
    }

    void setMetaData(const std::string& key, const Json::Value& value) {
        dataMetaInfo[key] = value;
    }

    uint32_t getDataLength() {
        return blockDataLength;
    }

    void setDataLength(uint32_t length) {
        blockDataLength = length;
    }

    uint8_t* getDataPtr() {
        return pBlockData;
    }

    void setDataPtr(uint8_t* dataPtr) {
        pBlockData = dataPtr;
    }

    int32_t setBlockData(uint8_t* data, uint32_t length);

    // /**
    //  * @brief Serialize the data block to a buffer.
    //  */
    // int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override;

    // /**
    //  * @brief Unserialize the data block from a buffer.
    //  */
    // int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    PbgzDataBlock():blockType(FILE_DATA), blockDataLength(0), pBlockData(nullptr){
        dataMetaInfo.clear();
        metaChecksum = 0;
        dataBlockChecksum = 0; 
    }

    virtual ~PbgzDataBlock() {
        dataMetaInfo.clear();
        // 外部地址，不能在此释放
        pBlockData = nullptr;
        
    }

    void calcChecksum();

    void setMetaCheckSum(uint64_t checkSum) {
        metaChecksum = checkSum;
    }

    void setDataCheckSum(uint64_t checkSum) {
        dataBlockChecksum = checkSum;
    }

    uint64_t getMetaCheckSum() {
        return metaChecksum;
    }

    uint64_t getDataCheckSum() {
        return dataBlockChecksum;
    }

    int32_t verifyCheckSum();

private:
    PbgzBlockType blockType;       // block type
    Json::Value dataMetaInfo;      // block meta infomation
    uint64_t metaChecksum;         // Checksum for meta data
    uint32_t blockDataLength;      // block Data length
    uint8_t *pBlockData;           // block data 
    uint64_t dataBlockChecksum;    // block data checksum
};
