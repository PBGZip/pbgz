/*
 * pbgz_file.h - Header file for pbgz_v2 project
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

#ifndef PBGZ_FILE_H
#define PBGZ_FILE_H

#include <string>
#include <json/json.h>

/// @brief Magic value for block in PBGZ file format, used to identify block type
typedef enum {
    INVALID = 0,  // unknown
    FILE_HEADER =  1, // File header information
    FILE_META = 2,    // File meta information
    FILE_DATA = 3,    // File data
}PbgzBlockType;

const std::string PBGZ_FILE_MAGIC = "PBGZ"; // Magic value for PBGZ file format
const uint32_t PBGZ_FILE_META_MAGIC = 0x0000EB; // Magic value for PBGZ file meta information
const uint32_t PBGZ_DATA_BLOCK_MAGIC = 0x0000DB; // Magic value for PBGZ file meta information


/// @brief File type in PBGZ file format
typedef enum {
    BINARY = 0, // Binary file
    FASTQ = 1,  // FASTQ format file
} FileFormat;

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
class PbgzFileHeader : public Serializable {
public:
    /**
     * @brief Serialize the file header to a buffer.
     */
    int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override ;

    /**
     * @brief Deserialize the file header from a buffer.
     */
    int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    /**
     * @brief Get the version string of the file header.
     */
    std::string getVersionStr();

    PbgzFileHeader() :blockType(INVALID) {
        memset(version, 0, sizeof(version));
    }

private:
    PbgzBlockType blockType;  // Block type
    uint8_t version[3];  // Version number
};

/// @brief PBGZ file meta information
class PbgzFileMeta : public Serializable{
public:
    /**
     * @brief Serialize the file meta information to a buffer.
     */
    int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override ;

    /**
     * @brief Unserialize the file meta information from a buffer.
     */
    int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    Json::Value getMetaData() const {
        return metaData;
    }

    Json::Value getMetaData(std::string& key) {
        return metaData[key];
    }

    void setMetaData(std::string &key, Json::Value &value) {
        metaData[key] = value;
    }

    void setMetaData(std::string &key, std::string&value) {
        metaData[key] = value;
    }

    PbgzFileMeta():blockType(INVALID), metaLength(0) {
        metaData.clear();
        memset(metaChecksum, 0, sizeof(metaChecksum));
    }

    virtual ~PbgzFileMeta() {
        metaData.clear();
    }

private:
    PbgzBlockType blockType;   // Block type
    uint32_t metaLength; // Length of meta data
    Json::Value metaData; // Meta information in JSON format
    uint8_t metaChecksum[8]; // Checksum for meta data
};


/// @brief  Pbgz file block data
class PbgzDataBlock : public Serializable {
public:
    PbgzBlockType getBlockType() const {
        return blockType;
    }

    void setBlockType(PbgzBlockType type) {
        blockType = type;
    }

    uint32_t getMetaLength() const {
        return metaLength;
    }

    void setMetaLength(uint32_t length) {
        metaLength = length;
    }

    Json::Value getMetaData() const {
        return metaData;
    }

    Json::Value getMetaData(const std::string &key) const {
        return metaData[key];
    }

    void setMetaData(const Json::Value& value) {
        metaData = value;
    }

    void setMetaData(const std::string& key, const Json::Value& value) {
        metaData[key] = value;
    }

    uint32_t getDataLength() const {
        return dataLength;
    }

    void setDataLength(uint32_t length) {
        dataLength = length;
    }

    const uint8_t* getDataPtr() const {
        return pData;
    }

    int32_t setBlockData(uint8_t* data, uint32_t length);

    /**
     * @brief Serialize the data block to a buffer.
     */
    int32_t serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) override;

    /**
     * @brief Unserialize the data block from a buffer.
     */
    int32_t unserialize(uint8_t* buffer, uint32_t bufferLength) override;

    PbgzDataBlock():blockType(INVALID), metaLength(0), dataLength(0), pData(nullptr){
        metaData.clear();
        memset(dataBlockChecksum, 0 ,sizeof(dataBlockChecksum));
        memset(originDataChecksum, 0, sizeof(originDataChecksum));
    }

protected:
    virtual ~PbgzDataBlock() {
        metaData.clear();
        if (!pData) {
            free(pData);
            pData = nullptr;
        }
    }

private:
    PbgzBlockType blockType;   // block type
    uint32_t metaLength;           // block meta infomation length
    Json::Value metaData;      // block meta infomation
    uint32_t dataLength;           // block Data length
    uint8_t *pData;                 // block data 
    uint8_t dataBlockChecksum[8];   // block data checksum
    uint8_t originDataChecksum[8];   // origin data checksum
};

#endif