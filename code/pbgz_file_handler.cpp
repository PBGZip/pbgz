/*
 * pbgz_file_handler.cpp - CPP file for pbgz_v2 project
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

#include "pbgz_file_handler.h"
#include "log/logger.h"

// 2MB read buffer size
const uint32_t PBGZ_FILE_READ_BUFFER_LENGTH = 2 * 1024 * 1024; 

uint8_t* PbgzFileReader::getFileReadBuffer(){
    thread_local static uint8_t* pReadBuffer = nullptr;
    if (!pReadBuffer) {
        pReadBuffer = new uint8_t[PBGZ_FILE_READ_BUFFER_LENGTH];
        if (!pReadBuffer) {
            LOG_STDOUT(LOG_ERROR, "Failed to allocate memory for read buffer.");
            return nullptr; // Memory allocation error
        }
    }
    return pReadBuffer;
}

int32_t PbgzFileReader::open() {
    // Implement the read logic for PBGZ file format
    // This is a placeholder implementation
    LOG_STDOUT(LOG_INFO, "Reading PBGZ file...");
    pFile = fopen(fileName.c_str(), "rb");
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "Failed to open file: %s", fileName.c_str());
        return -1; // File open error
    }

    if (0 != initFileHeadAndMeta()) {
        fclose(pFile);
        return -1; // Initialization error
    }

    return 0;
}

int32_t PbgzFileReader::initFileHeadAndMeta(bool isCheckMagic) {
    // Read the file header
    uint8_t* pReadBuffer = getFileReadBuffer();
    if (!pReadBuffer) {
        LOG_STDOUT(LOG_ERROR, "Failed to get read buffer.");   
        return -1; 
    }

    // Read the magic value, 4 byte
    size_t readLen = 0;
    // If isCheckMagic is true, we will not read the magic value again
    if (isCheckMagic) {
        memcpy(pReadBuffer, PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH);
    } else {
        readLen = fread(pReadBuffer, PBGZ_FILE_MAGIC_LENGTH, 1, pFile);
        if (readLen !=  PBGZ_FILE_MAGIC_LENGTH) {
            LOG_STDOUT(LOG_ERROR, "Failed to read file header from: %s", fileName.c_str());
            return -1; // File read error
        }

        if (memcmp(pReadBuffer, PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH) != 0) {
            // 安照16进制打印，自行解析
            LOG_STDOUT(LOG_ERROR, "%s is not a valid pbgz file, magic no is %X", fileName.c_str(), 
                (uint32_t)(*(uint32_t*)pReadBuffer));
            return -1; // Invalid magic
        }
    }

    // Read the version number, 3 byte
    readLen = fread(pReadBuffer + PBGZ_FILE_MAGIC_LENGTH, 1, PBGZ_FILE_VERSION_LENGTH, pFile);
    if (readLen != PBGZ_FILE_VERSION_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read file version from: %s", fileName.c_str());
        return -1; // File read error
    }   

    LOG_STDOUT(LOG_INFO, "PBGZ file opened successfully: %s", fileName.c_str());

    PbgzFileHeader header;
    header.unserialize(reinterpret_cast<uint8_t*>(pReadBuffer), PBGZ_FILE_MAGIC_LENGTH + PBGZ_FILE_VERSION_LENGTH);
    fileHeaderMap[currentFileIndex++] = header;

    // Read file meta information
    // read file meta magic
    readLen = fread(pReadBuffer, PBGZ_FILE_META_MAGIC_LENGTH, 1, pFile);
    if (readLen != PBGZ_FILE_META_MAGIC_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "%s is not a valid pbgz file.", fileName.c_str());
        return -1;
    }

    if (0 != memcmp(pReadBuffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH)) {
        LOG_STDOUT(LOG_ERROR, "%s is not a valid pbgz file, file meta magic no is %X", fileName.c_str(), 
            (uint32_t)(*(uint32_t*)pReadBuffer));
        return -1;
    }

    // read file meta length
    readLen = fread(pReadBuffer + PBGZ_FILE_META_MAGIC_LENGTH, 1, PBGZ_FILE_META_SIZE_LENGTH, pFile);
    if (readLen != PBGZ_FILE_META_SIZE_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read file meta length from: %s", fileName.c_str());
        return -1; // File read error       
    }
    const uint32_t metaLength = (uint32_t)(*(uint32_t*)(pReadBuffer + PBGZ_FILE_META_MAGIC_LENGTH));

    // read the rest of the meta data and checksum
    readLen = fread(pReadBuffer + PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH, 1, metaLength + PBGZ_FILE_META_CHECKSUM_LENGTH, pFile);
    if (readLen != metaLength + PBGZ_FILE_META_CHECKSUM_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read file meta data from: %s", fileName.c_str());
        return -1; // File read error   
    }

    PbgzFileMeta fileMeta;
    // file meta 的总长度：maigic头（4）+ metaLength（4）+ metaData + checksum（8）
    uint32_t fileMetaDataLength = PBGZ_FILE_META_MAGIC_LENGTH + PBGZ_FILE_META_SIZE_LENGTH + metaLength + PBGZ_FILE_META_CHECKSUM_LENGTH;
    int ret = fileMeta.unserialize(pReadBuffer, fileMetaDataLength);
    if (ret != 0) {
        LOG_STDOUT(LOG_ERROR, "Failed to unserialize file meta data from: %s", fileName.c_str());
        return -1; // Unserialization error
    }

    // 保持和header相同的序列号
    fileMetaMap[currentFileIndex] = fileMeta;
    return 0;
}       

void PbgzFileReader::close() {
    if (pFile) {
        fclose(pFile);
        pFile = nullptr;
    }

    return;
}

PbgzFileReader::~PbgzFileReader() {
    if (pFile) {
        fclose(pFile);
    }
}

const PbgzFileHeader& PbgzFileReader::getFileHeader() {
    if (currentFileIndex < 0) {
        LOG_STDOUT(LOG_ERROR, "No file has been read yet.");
        throw std::runtime_error("No file has been read yet.");
    }

    if (fileHeaderMap.find(currentFileIndex) == fileHeaderMap.end()) {
        LOG_STDOUT(LOG_ERROR, "File header not found for current file index: %d", currentFileIndex);
        throw std::runtime_error("File header not found for current file index.");
    }
    
    return fileHeaderMap[currentFileIndex];
}

const PbgzFileMeta& PbgzFileReader::getFileMeta() {
    // if currentFileIndex is -1, it means no file has been read yet
    if (currentFileIndex < 0) {
        LOG_STDOUT(LOG_ERROR, "No file has been read yet.");
        throw std::runtime_error("No file has been read yet.");
    }

    // Return the file meta information for the current file index
    if (fileMetaMap.find(currentFileIndex) == fileMetaMap.end()) {
        LOG_STDOUT(LOG_ERROR, "File meta information not found for current file index: %", currentFileIndex);   
        throw std::runtime_error("File meta information not found for current file index.");
    }

    return fileMetaMap[currentFileIndex];
}

int32_t PbgzFileReader::readDataBlock(PbgzDataBlock& dataBlock) {
    LOG_STDOUT(LOG_INFO, "Reading PBGZ data block...");
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "File is not open.");
        return -1; // File not open
    }   

    uint8_t* pReadBuffer = getFileReadBuffer();
    if (!pReadBuffer) {
        LOG_STDOUT(LOG_ERROR, "Failed to get read buffer.");        
        return -1; // Memory allocation error
    }

    // Read the data block magic value, 4byte
    size_t readLen = fread(pReadBuffer, PBGZ_DATA_BLOCK_MAGIC_LENGTH, 1, pFile);
    if (readLen != PBGZ_DATA_BLOCK_MAGIC_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block header.");
        return -1; // File read error
    }

    if (memcmp(pReadBuffer, &PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH) != 0) {
        LOG_STDOUT(LOG_ERROR, "Invalid magic value for PBGZ data block. %X", (uint32_t)(*(uint32_t*)pReadBuffer));
        // It maybe a new file, try to parse the file header and meta information
         if (memcmp(pReadBuffer, &PBGZ_FILE_MAGIC, PBGZ_FILE_MAGIC_LENGTH) == 0) {
            LOG_STDOUT(LOG_INFO, "Detected a new PBGZ file format, reinitializing file header and meta.");
            // 由于已经读取了文件的Maigc值,不需要在读取了
            if (0 != initFileHeadAndMeta(true)) {
                LOG_STDOUT(LOG_ERROR, "Failed to initialize file header and meta.");
                return -1; // Initialization error
            }

            // continue to read the data block from new file
            return readDataBlock(dataBlock); 
        }
        return -1; // Invalid magic
    }

    uint32_t offset = PBGZ_DATA_BLOCK_MAGIC_LENGTH;
    // Read the meta length, 4byte
    readLen = fread(pReadBuffer + offset, 1, PBGZ_DATA_BLOCK_META_SIZE_LENGTH, pFile);
    if (readLen != PBGZ_DATA_BLOCK_META_SIZE_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block meta size.");
        return -1; // File read error           
    }

    uint32_t metaLength = (uint32_t)(*(uint32_t*)(pReadBuffer + PBGZ_DATA_BLOCK_MAGIC_LENGTH));
    if (metaLength == 0) {
        LOG_STDOUT(LOG_ERROR, "Meta length is zero, invalid data block.");
        return -1; // Invalid meta length       
    }
    offset += PBGZ_DATA_BLOCK_META_SIZE_LENGTH;

    // Read the meta data
    readLen = fread(pReadBuffer + offset, 1, metaLength, pFile);
    if (readLen != metaLength) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block meta data."); 
        return -1; // File read error   
    }
    offset += metaLength;

    // Read the meta checksum, 8byte
    readLen = fread(pReadBuffer + offset , 1, PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH, pFile);
    if (readLen != PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block meta checksum.");
        return -1; // File read error
    }
    offset += PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH;

    // Read the data length, 4byte
    readLen = fread(pReadBuffer + offset, 1, PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH, pFile);
    if (readLen != PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block data length.");
        return -1; // File read error
    }   
    
    uint32_t dataLength = (uint32_t)(*(uint32_t*)(pReadBuffer + offset));
    if (dataLength == 0) {
        LOG_STDOUT(LOG_ERROR, "Data length is zero, invalid data block.");  
        return -1; // Invalid data length
    }

    // Read the block data
    offset += PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH;
    readLen = fread(pReadBuffer + offset, 1, dataLength, pFile);
    if (readLen != dataLength) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block data.");
        return -1; // File read error
    }

    // Read the data block checksum, 8byte
    offset += dataLength;
    readLen = fread(pReadBuffer + offset, 1, PBGZ_DATA_BLOCK_CHECKSUM_LENGTH, pFile);
    if (readLen != PBGZ_DATA_BLOCK_CHECKSUM_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block checksum.");
        return -1; // File read error
    }

    // Read the original data checksum, 8byte
    offset += PBGZ_DATA_BLOCK_CHECKSUM_LENGTH;
    readLen = fread(pReadBuffer + offset, 1, PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH, pFile);
    if (readLen != PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH) {
        LOG_STDOUT(LOG_ERROR, "Failed to read PBGZ data block original checksum.");     
        return -1; // File read error
    }

    // Unserialize the data block
    int ret = dataBlock.unserialize(pReadBuffer, offset + PBGZ_DATA_BLOCK_ORIGIN_CHECKSUM_LENGTH);
    if (ret != 0) {
        LOG_STDOUT(LOG_ERROR, "Failed to unserialize PBGZ data block.");
        return -1; // Unserialization error 
    }
        
    return 0;
}   


uint8_t* PbgzFileWriter::getFileWriteBuffer() {
    thread_local static uint8_t* pWriteBuffer = nullptr;
    if (!pWriteBuffer) {
        pWriteBuffer = new uint8_t[PBGZ_FILE_READ_BUFFER_LENGTH];
        if (!pWriteBuffer) {
            LOG_STDOUT(LOG_ERROR, "Failed to allocate memory for write buffer.");
            return nullptr; // Memory allocation error
        }
    }
    return pWriteBuffer;
}

int32_t PbgzFileWriter::open(){
    LOG_STDOUT(LOG_INFO, "Reading PBGZ file...");
    if (fileName.length() == 0) {
        LOG_STDOUT(LOG_ERROR, "File name is empty.");
        return -1; // Invalid file name
    }
    
    pFile = fopen(fileName.c_str(), "w+");
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "Failed to open file: %s", fileName.c_str());
        return -1; // File open error
    }

    if (0 != initFileHead()) {
        fclose(pFile);
        pFile = nullptr;
        return -1; // Initialization error
    }

    return 0;
}

int32_t PbgzFileWriter::close(){
    if (pFile) {
        fclose(pFile);
        pFile = nullptr;
    }

    return 0;
}

int32_t PbgzFileWriter::initFileHead() {
     if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "File is not open.");
        return -1; // File not open
    }

    uint8_t* pWriteBuffer = getFileWriteBuffer();
    if (!pWriteBuffer) {
        LOG_STDOUT(LOG_ERROR, "Failed to get write buffer.");
        return -1; // Memory allocation error   
    }

    uint32_t dataLength = 0;
    // Serialize the file header to the write buffer
    fileHeader.serialize(pWriteBuffer, PBGZ_FILE_READ_BUFFER_LENGTH, dataLength);
    if (dataLength == 0) {
        LOG_STDOUT(LOG_ERROR, "Failed to serialize file header.");
        return -1; // Serialization error
    }

    fwrite(pWriteBuffer, 1, dataLength, pFile);
    fflush(pFile);
    LOG_STDOUT(LOG_INFO, "File header initialized successfully.");
    return 0;
}

int32_t PbgzFileWriter::writeFileMeta() {
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "File is not open.");
        return -1; // File not open
    }

    uint8_t* pWriteBuffer = getFileWriteBuffer();
    if (!pWriteBuffer) {
        LOG_STDOUT(LOG_ERROR, "Failed to get write buffer.");
        return -1; // Memory allocation error
    }

    uint32_t dataLength = 0;
    // Serialize the file meta to the write buffer
    int32_t ret = fileMeta.serialize(pWriteBuffer, PBGZ_FILE_READ_BUFFER_LENGTH, dataLength);
    if (ret != 0) {
        LOG_STDOUT(LOG_ERROR, "Failed to serialize file meta.");
        return -1; // Serialization error
    }

    if (dataLength == 0) {
        LOG_STDOUT(LOG_ERROR, "No data to serialize file meta.");
        return 0; // No data to write
    }

    fwrite(pWriteBuffer, 1, dataLength, pFile);
    fflush(pFile);
    LOG_STDOUT(LOG_INFO, "File meta written successfully.");
    
    return 0;
}

int32_t PbgzFileWriter::writeBlockData(PbgzDataBlock& dataBlock) {
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "File is not open.");
        return -1; // File not open
    }

    uint8_t* pWriteBuffer = getFileWriteBuffer();
    if (!pWriteBuffer) {
        LOG_STDOUT(LOG_ERROR, "Failed to get write buffer.");
        return -1; // Memory allocation error
    }

    uint32_t dataLength = 0;
    int32_t ret = dataBlock.serialize(pWriteBuffer, PBGZ_FILE_READ_BUFFER_LENGTH, dataLength);
    if (ret != 0) {
        LOG_STDOUT(LOG_ERROR, "Failed to serialize block data.");
        return -1; // Serialization error
    }

    if (dataLength == 0) {
        LOG_STDOUT(LOG_ERROR, "No data to serialize block data.");
        return 0; // No data to write
    }

    fwrite(pWriteBuffer, 1, dataLength, pFile);
    fflush(pFile);
    LOG_STDOUT(LOG_INFO, "Block data written successfully.");
    
    return 0;
}

