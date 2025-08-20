/*
 * pbgz_file.cpp - cpp file for pbgz_v2 project
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
    if (version[0] == '\0' || version[1] == '\0' || version[2] == '\0') {
        return "Unknown Version";
    }
    // Ensure that version elements are valid characters    
    char buffer[12];
    snprintf(buffer, sizeof(buffer), "%d.%d.%d", version[0], version[1], version[2]);
    return std::string(buffer);
}

int32_t PbgzFileHeader::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    if (blockType != FILE_HEADER) {
        
        return -1;
    }
    if (bufferLength < PBGZ_FILE_MAGIC.length() + sizeof(version)) {
        return -1; // Buffer too small
    }

    // Serialize block type, 4 byte without '\0'
    memcpy(buffer, &PBGZ_FILE_MAGIC, PBGZ_FILE_MAGIC.length());
    // Serialize version
    memcpy(buffer + PBGZ_FILE_MAGIC.length(), version, sizeof(version));

    dataLength = PBGZ_FILE_MAGIC.length() + sizeof(version);
    return 0; 
}

int32_t PbgzFileHeader::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < PBGZ_FILE_MAGIC.length() + sizeof(version)) {
        return -1; // Buffer too small
    }

    // Unserialize block type
    memcpy(&blockType, buffer, PBGZ_FILE_MAGIC.length());
    // Unserialize version
    memcpy(version, buffer + PBGZ_FILE_MAGIC.length(), sizeof(version));

    return 0; 
}

int32_t PbgzFileMeta::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    Json::FastWriter writer;
    std::string jsonString = writer.write(metaData);
    metaLength = jsonString.size();
    if (bufferLength < sizeof(PBGZ_FILE_META_MAGIC) + metaLength + sizeof(metaChecksum)) {
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Serialize block type
    memcpy(buffer, &blockType, sizeof(PBGZ_FILE_META_MAGIC));
    dataOffset += sizeof(PBGZ_FILE_META_MAGIC);

    // Serialize meta data length
    memcpy(buffer + dataOffset, &metaLength, sizeof(metaLength));
    dataOffset += sizeof(metaLength);

    // Serialize meta data
    memcpy(buffer + dataOffset, jsonString.c_str(), metaLength);
    dataOffset += metaLength;

    // Serialize checksum
    memcpy(buffer + dataOffset, metaChecksum, sizeof(metaChecksum));

    dataLength = sizeof(PBGZ_FILE_META_MAGIC) + jsonString.size() + sizeof(metaChecksum);
    return 0; 
}

int32_t PbgzFileMeta::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < sizeof(PBGZ_FILE_META_MAGIC) + sizeof(metaChecksum)) {
        return -1; // Buffer too small
    }

    uint dataOffset = 0;
    // Unserialize block type
    memcpy(&blockType, buffer, sizeof(PBGZ_FILE_META_MAGIC));
    dataOffset += sizeof(PBGZ_FILE_META_MAGIC);

    // Unserialize meta data length
    if (bufferLength < dataOffset + sizeof(metaLength)) {
        return -   1 ; // Buffer too small for metaLength   
    }
    memcpy(&metaLength, buffer + dataOffset, sizeof(metaLength));
    dataOffset += sizeof(metaLength);   

    // Unserialize meta data
    std::string jsonString(reinterpret_cast<char*>(buffer + dataOffset), metaLength);
    Json::CharReaderBuilder readerBuilder;
    std::istringstream jsonStream(jsonString);
    std::string errs;
    if (!Json::parseFromStream(readerBuilder, jsonStream, &metaData, &errs)) {
        std::cerr << "Failed to parse JSON: " << errs << std::endl;
        return -1; // JSON parse error
    }
    dataOffset += metaLength;

    // Unserialize checksum
    memcpy(metaChecksum, buffer + dataOffset, sizeof(metaChecksum));

    return 0; 
}


int32_t PbgzDataBlock::serialize(uint8_t* buffer, uint32_t bufferLength, uint32_t& dataLength) {
    Json::FastWriter writer;
    std::string jsonString = writer.write(metaData);
    metaLength = jsonString.size();
    if (bufferLength < sizeof(PBGZ_DATA_BLOCK_MAGIC) + sizeof(metaLength) + metaLength + sizeof(dataLength) + dataLength + sizeof(dataBlockChecksum) + sizeof(originDataChecksum)) {
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Serialize block type
    memcpy(buffer, &blockType, sizeof(PBGZ_DATA_BLOCK_MAGIC));
    dataOffset += sizeof(PBGZ_DATA_BLOCK_MAGIC);

    // Serialize meta length
    memcpy(buffer + dataOffset, &metaLength, sizeof(metaLength));
    dataOffset += sizeof(metaLength);

    // Serialize meta data
    memcpy(buffer + dataOffset, jsonString.c_str(), metaLength);
    dataOffset += metaLength;

    // Serialize data length
    memcpy(buffer + dataOffset, &dataLength, sizeof(dataLength));
    dataOffset += sizeof(dataLength);

    // Serialize block data
    if (pData && dataLength > 0) {
        memcpy(buffer + dataOffset, pData, dataLength);
        dataOffset += dataLength;
    }

    // Serialize checksums
    memcpy(buffer + dataOffset, dataBlockChecksum, sizeof(dataBlockChecksum));
    dataOffset += sizeof(dataBlockChecksum);
    
    memcpy(buffer + dataOffset, originDataChecksum, sizeof(originDataChecksum));
    dataLength = sizeof(PBGZ_DATA_BLOCK_MAGIC) + sizeof(metaLength) + metaLength + sizeof(dataLength) + dataLength + sizeof(dataBlockChecksum) + sizeof(originDataChecksum);
    
    return 0; 
}


int32_t PbgzDataBlock::unserialize(uint8_t* buffer, uint32_t bufferLength) {
    if (bufferLength < sizeof(PBGZ_DATA_BLOCK_MAGIC) + sizeof(metaLength) + sizeof(dataLength) + sizeof(dataBlockChecksum) + sizeof(originDataChecksum)) {
        return -1; // Buffer too small
    }

    uint32_t dataOffset = 0;
    // Unserialize block type
    memcpy(&blockType, buffer, sizeof(PBGZ_DATA_BLOCK_MAGIC));
    dataOffset += sizeof(PBGZ_DATA_BLOCK_MAGIC);

    // Unserialize meta length
    memcpy(&metaLength, buffer + dataOffset, sizeof(metaLength));
    dataOffset += sizeof(metaLength);

    // Unserialize meta data
    if (bufferLength < dataOffset + metaLength) {
        return -1; // Buffer too small for meta data
    }
    std::string jsonString(reinterpret_cast<char*>(buffer + dataOffset), metaLength);
    Json::CharReaderBuilder readerBuilder;
    std::istringstream jsonStream(jsonString);
    std::string errs;
    if (!Json::parseFromStream(readerBuilder, jsonStream, &metaData, &errs)) {
        LOG_STDOUT(LOG_FATAL, "Failed to parse JSON: %s", errs.c_str());
        return -1; // JSON parse error
    }
    dataOffset += metaLength;

    // Unserialize data length
    memcpy(&dataLength, buffer + dataOffset, sizeof(dataLength));
    dataOffset += sizeof(dataLength);

    if (bufferLength < dataOffset + dataLength ) {
        return -1; // Buffer too small for data
    }

    // Unserialize block data
    if (dataLength > 0) {
        pData = static_cast<uint8_t*>(malloc(dataLength));
        if (!pData) {
            return -1; // Memory allocation failed
        }
        memcpy(pData, buffer + dataOffset, dataLength);
        dataOffset += dataLength;
    }

    // Unserialize checksums
    memcpy(dataBlockChecksum, buffer + dataOffset, sizeof(dataBlockChecksum));
    dataOffset += sizeof(dataBlockChecksum);
    
    memcpy(originDataChecksum, buffer + dataOffset, sizeof(originDataChecksum));

    return 0; 
}

int32_t PbgzDataBlock::setBlockData(uint8_t* data, uint32_t length) {
    if (pData) {
        free(pData);
        pData = nullptr;
    }
    
    pData = static_cast<uint8_t*>(malloc(length));
    if (!pData) {
        std::cerr << "Failed to allocate memory for block data." << std::endl;
        return -1;
    }

    memcpy(pData, data, length);
    return 0;
}



