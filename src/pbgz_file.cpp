/*
 * pbgz_file.cpp - cpp file for pbgz project
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

#include <iostream>

#include "pbgz_file.h"
#include "log/logger.h"


std::string PbgzFileHeader::getVersionStr() {
    // Convert version array to string format "x.x.x"
    // Assuming version is in the format [major, minor, patch]
    // where each element is a single character representing a digit.
    // Adjust the conversion if the version format is different.
    if (version[0] == '\0') {
        return "Unknown Version";
    }
    // Ensure that version elements are valid characters    
    char buffer[12];
    snprintf(buffer, sizeof(buffer), "%d.%d.%d", version[0], version[1], version[2]);
    return std::string(buffer);
}

int32_t PbgzFileHeader::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    if (blockType != FILE_HEADER) {
        LOG_STDOUT(LOG_ERROR, "Invalid block type for serialization: %d", blockType);
        return -1;
    }
    if (bufferLength < PBGZ_FILE_MAGIC_LENGTH + sizeof(version)) {
        return -1; // Buffer too small
    }

    // Serialize block type, 4 byte without '\0'
    memcpy(buffer, PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH);
    // Serialize version
    memcpy(buffer + PBGZ_FILE_MAGIC_LENGTH, version, sizeof(version));

    dataLength = PBGZ_FILE_MAGIC_LENGTH + sizeof(version);
    return 0; 
}

int32_t PbgzFileHeader::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < PBGZ_FILE_MAGIC_LENGTH + sizeof(version)) {
        return -1; // Buffer too small
    }

    // Unserialize block type
    if (strncmp(reinterpret_cast<const char*>(buffer), PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH) != 0) {
        LOG_STDOUT(LOG_ERROR, "Invalid magic value for PBGZ file header.");
        return -1; // Invalid magic
    }
    blockType = FILE_HEADER;

    // Unserialize version
    memcpy(version, buffer + PBGZ_FILE_MAGIC_LENGTH, sizeof(version));

    return 0; 
}

int32_t PbgzFileMeta::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    if (blockType != FILE_META) {
        LOG_STDOUT(LOG_ERROR, "Invalid block type for serialization: %d", blockType);
        return -1;
    }   

    Json::StreamWriterBuilder writer;
    std::string jsonString =  Json::writeString(writer, metaData);
    metaLength = jsonString.size();
    if (bufferLength < PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH + metaLength + PBGZ_FILE_META_CHECKSUM_LENGTH) {
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Serialize block type
    memcpy(buffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH);
    dataOffset += PBGZ_FILE_META_MAGIC_LENGTH;

    // Serialize meta data length
    memcpy(buffer + dataOffset, &metaLength, PBGZ_FILE_META_SIZE_LENGTH);
    dataOffset += PBGZ_FILE_META_SIZE_LENGTH;

    // Serialize meta data
    // 暂时没有加密，待实现
    /// TODO: Implement meta data encryption if needed
    memcpy(buffer + dataOffset, jsonString.c_str(), metaLength);
    dataOffset += metaLength;

    // Serialize checksum
    memcpy(buffer + dataOffset, &metaChecksum, PBGZ_FILE_META_CHECKSUM_LENGTH);

    dataLength = dataOffset + PBGZ_FILE_META_CHECKSUM_LENGTH;
    return 0; 
}

int32_t PbgzFileMeta::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH + sizeof(metaChecksum)) {
        return -1; // Buffer too small
    }

    uint dataOffset = 0;
    // Unserialize block type
    if (memcmp(buffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH) != 0) {
        LOG_STDOUT(LOG_ERROR, "Invalid magic value for PBGZ file meta.");
        return -1; // Invalid magic
    }
    blockType = FILE_META;
    dataOffset += PBGZ_FILE_META_MAGIC_LENGTH;

    // Unserialize meta data length
    if (bufferLength < dataOffset + PBGZ_FILE_META_SIZE_LENGTH) {
        return -   1 ; // Buffer too small for metaLength   
    }
    memcpy(&metaLength, buffer + dataOffset, sizeof(metaLength));
    dataOffset += PBGZ_FILE_META_SIZE_LENGTH;   

    // Unserialize meta data
    std::string jsonString(reinterpret_cast<char*>(buffer + dataOffset), metaLength);
    Json::CharReaderBuilder readerBuilder;
    std::istringstream jsonStream(jsonString);
    std::string errs;
    if (!Json::parseFromStream(readerBuilder, jsonStream, &metaData, &errs)) {
        LOG_STDOUT(LOG_ERROR, "Failed to parse JSON: %s", errs);
        return -1; // JSON parse error
    }
    dataOffset += metaLength;

    // Unserialize checksum
    metaChecksum = *(uint64_t*)(buffer + dataOffset);

    return 0; 
}

int32_t PbgzDataBlock::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    if (blockType != FILE_DATA) {
        LOG_STDOUT(LOG_ERROR, "Invalid block type for serialization: %d", blockType);
        return -1;
    }

    Json::StreamWriterBuilder writer;
    std::string jsonString = Json::writeString(writer, dataMetaInfo);
    dataMetaLength = jsonString.size();
    if (bufferLength < PBGZ_DATA_BLOCK_MAGIC_LENGTH +  PBGZ_DATA_BLOCK_META_SIZE_LENGTH + dataMetaLength + 
        PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH + blockDataLength + PBGZ_DATA_BLOCK_CHECKSUM_LENGTH + PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Buffer too small for serialization");
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Serialize block type
    memcpy(buffer, &PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_MAGIC_LENGTH;

    // Serialize meta length
    memcpy(buffer + dataOffset, &dataMetaLength, PBGZ_DATA_BLOCK_META_SIZE_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_META_SIZE_LENGTH;

    // Serialize meta data
    memcpy(buffer + dataOffset, jsonString.c_str(), dataMetaLength);
    dataOffset += dataMetaLength;

    // Serialize meta checksum
    memcpy(buffer + dataOffset, &metaChecksum, PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH;

    // Serialize data length
    memcpy(buffer + dataOffset, &blockDataLength, PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH;

    // Serialize block data
    if (pBlockData && blockDataLength > 0) {
        memcpy(buffer + dataOffset, pBlockData, blockDataLength);
        dataOffset += blockDataLength;
    }

    // Serialize checksums
    memcpy(buffer + dataOffset, &dataBlockChecksum, PBGZ_DATA_BLOCK_CHECKSUM_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_CHECKSUM_LENGTH;
    
    memcpy(buffer + dataOffset, &originDataChecksum, PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH );
    dataLength = dataOffset + PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH;
    
    return 0; 
}

int32_t PbgzDataBlock::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < PBGZ_DATA_BLOCK_MAGIC_LENGTH + PBGZ_DATA_BLOCK_META_SIZE_LENGTH + PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH
         + PBGZ_DATA_BLOCK_CHECKSUM_LENGTH + PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH) {
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Unserialize block type
    if (memcmp(buffer, &PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH) != 0) {
        LOG_STDOUT(LOG_ERROR, "Invalid magic value for PBGZ data block.");
        return -1; // Invalid magic
    }

    blockType = FILE_DATA;
    dataOffset += PBGZ_DATA_BLOCK_MAGIC_LENGTH;

    // Unserialize meta length
    memcpy(&dataMetaLength, buffer + dataOffset, PBGZ_DATA_BLOCK_META_SIZE_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_META_SIZE_LENGTH;

    // Unserialize meta data
    if (bufferLength < dataOffset + dataMetaLength) {
        return -1; // Buffer too small for meta data
    }
    std::string jsonString(reinterpret_cast<char*>(buffer + dataOffset), dataMetaLength);
    Json::CharReaderBuilder readerBuilder;
    std::istringstream jsonStream(jsonString);
    std::string errs;
    if (!Json::parseFromStream(readerBuilder, jsonStream, &dataMetaInfo, &errs)) {
        LOG_STDOUT(LOG_FATAL, "Failed to parse JSON: %s", errs.c_str());
        return -1; // JSON parse error
    }
    dataOffset += dataMetaLength;

    // Unserialize meta checksum
    if (bufferLength < dataOffset + PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH) {
        return -1; // Buffer too small for meta checksum
    }
    metaChecksum = *(uint64_t*)(buffer + dataOffset);
    dataOffset += PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH; 

    // Unserialize data length
    memcpy(&blockDataLength, buffer + dataOffset, PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH);
    dataOffset += PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH;

    if (bufferLength < dataOffset + blockDataLength ) {
        return -1; // Buffer too small for data
    }

    // Unserialize block data
    if (blockDataLength > 0) {
        pBlockData = static_cast<uint8_t*>(malloc(blockDataLength));
        if (!pBlockData) {
            LOG_STDOUT(LOG_FATAL, "Failed to allocate memory for block data.");
            return -1; // Memory allocation failed
        }
        memcpy(pBlockData, buffer + dataOffset, blockDataLength);
        dataOffset += blockDataLength;
    }

    // Unserialize checksums
    dataBlockChecksum = *(uint64_t*)(buffer + dataOffset);
    dataOffset += PBGZ_DATA_BLOCK_CHECKSUM_LENGTH;
    
    originDataChecksum = *(uint64_t*)(buffer + dataOffset);
    return 0; 
}

int32_t PbgzDataBlock::setBlockData(uint8_t* data, uint32_t length) {
    if (pBlockData) {
        free(pBlockData);
        pBlockData = nullptr;
    }
    
    pBlockData = static_cast<uint8_t*>(malloc(length));
    if (!pBlockData) {
        std::cerr << "Failed to allocate memory for block data." << std::endl;
        return -1;
    }

    blockDataLength = length;
    memcpy(pBlockData, data, length);
    return 0;
}



