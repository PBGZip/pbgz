/*
 * pbgz_file_handler.cpp - CPP file for pbgz project
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

#include "pbgz_file_wrapper.h"
#include "log/logger.h"
#include "coder_json.h"
#include "utils/memory_util.h"

// Buffer as temporary storage for parsing meta and other information, no need to be too large
const uint32_t PBGZ_FILE_READ_BUFFER_LENGTH = 16 * 1024 * 1024;

uint8_t* PbgzFileReader::getFileReadBuffer(){
	// 线程局部智能指针：自动管理内存，线程退出时析构释放
    thread_local static std::unique_ptr<uint8_t[]> pReadBuffer;
    
    if (!pReadBuffer) {
        // 用 unique_ptr 分配内存（无需手动 delete）
        pReadBuffer = std::make_unique<uint8_t[]>(PBGZ_FILE_READ_BUFFER_LENGTH);
        if (!pReadBuffer) { // make_unique 失败时返回空指针（C++17+）
            LOG_ERROR("Failed to allocate memory for write buffer.");
        }
    }

    // 返回智能指针指向的原始缓冲区地址（供业务使用）
    return pReadBuffer.get();
}

int32_t PbgzFileReader::open() {
     // Implement the read logic for PBGZ file format
    if (0 != initFileHeadAndMeta()) {
        LOG_ERROR("Create PbgzFileReader, load head and meta failed");
        return -1;
    }
    return 0;
}

int32_t PbgzFileReader::initFileHeadAndMeta(bool isCheckMagic) {
    if (ioReader == nullptr) {
        LOG_ERROR("IO reader is NULL.");
        return -1;
    }
    // Read the file header
    uint8_t* pReadBuffer = getFileReadBuffer();
    if (pReadBuffer == nullptr) {
        LOG_ERROR("Failed to get read buffer.");   
        return -1; 
    }

    PbgzFileHeader fileHeader;
    // Read the magic value, 4 byte
    size_t readLen = 0;
    // If isCheckMagic is true, we will not read the magic value again
    if (!isCheckMagic) {
        readLen = ioReader->readIO(pReadBuffer, PBGZ_FILE_MAGIC_LENGTH);
        if (readLen !=  PBGZ_FILE_MAGIC_LENGTH) {
            LOG_ERROR("Failed to read file header from io, readLen=%d", readLen);
            return -1; // File read error
        }
        if (memcmp(pReadBuffer, PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH) != 0) {
            // Print in hexadecimal format, parse manually
            LOG_ERROR("IO is not a valid pbgz format, magic no is %X", (uint32_t)(*(uint32_t*)pReadBuffer));
            return -1; // Invalid magic
        }
    }
    fileHeader.setBlockType(FILE_HEADER);

    // Read the version number, 3 byte
    readLen = ioReader->readIO(pReadBuffer, PBGZ_FILE_VERSION_LENGTH);
    if (readLen != PBGZ_FILE_VERSION_LENGTH) {
        LOG_ERROR("Failed to read file version from IO");
        return -1; // File read error
    }   
    // Unserialize version
    fileHeader.setVersion(pReadBuffer + PBGZ_FILE_MAGIC_LENGTH, PBGZ_FILE_VERSION_LENGTH);
    
    // 读取文件头扩展部分，
    readLen = ioReader->readIO(pReadBuffer, sizeof(uint32_t));
    if (readLen > 0) {    
        uint32_t fileHeadSize = *(uint32_t*)pReadBuffer;
        if (fileHeadSize > 0) {
            readLen = ioReader->readIO(pReadBuffer, fileHeadSize);
            if (readLen >= sizeof(uint64_t)) {
                fileHeader.setDynamicMetaOffset(*(uint64_t*)pReadBuffer);
            }
        }
    }
    
    fileHeaderMap[++currentFileIndex] = fileHeader; 

    // Read file meta information
    PbgzFileMeta baseFileMeta;
    baseFileMeta.setMetaType(BASE_FILE_META);
    if (0 != readFileMeta(baseFileMeta)) {
        LOG_ERROR("Read base file meta failed.");
        return -1;
    }
    // Keep the same serial number as header
    baseFileMetaMap[currentFileIndex] = baseFileMeta;

    if (fileHeader.getDynamicMetaOffset() == 0) {
        return 0;
    }
    uint64_t dynamicOffset = fileHeader.getDynamicMetaOffset();
    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    if (fileReader == nullptr) {
        return 0;   // 如果从管道读取，转换会失败，属于正常的
    }

    size_t readOffset = fileReader->getCurrentPos(); // 备份当前读取的位置
    fileReader->seekIO(dynamicOffset);

    PbgzFileMeta dyncFileMeta;
    dyncFileMeta.setMetaType(DYNAMIC_FILE_META);
    if (0 != readFileMeta(dyncFileMeta)) {
        LOG_ERROR("Read dynamic file meta failed. offset = %llu", dynamicOffset);
        return -1;
    }
    dynamicFileMetaMap[currentFileIndex] = dyncFileMeta;
    fileReader->seekIO(readOffset);

    return 0;
}       

int32_t PbgzFileReader::readFileMeta(PbgzFileMeta& fileMeta, bool isCheckMagic) {
     if (ioReader == nullptr) {
        LOG_ERROR("IO reader is NULL.");
        return -1;
    }
    // Read the file header
    uint8_t* pReadBuffer = getFileReadBuffer();
    if (pReadBuffer == nullptr) {
        LOG_ERROR("Failed to get read buffer.");   
        return -1; 
    }

    size_t readLen = 0;

    // read file meta magic
    if (isCheckMagic) {
        readLen = ioReader->readIO(pReadBuffer, PBGZ_FILE_META_MAGIC_LENGTH);
        if (readLen == 0) {
            LOG_INFO("No file meta information found in IO");
            return 0; // No file meta information, not an error
        }
        if (readLen != PBGZ_FILE_META_MAGIC_LENGTH) {
            LOG_ERROR("IO is not a valid pbgz file.");
            return -1;
        }
        if (0 != memcmp(pReadBuffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH)) {
            LOG_ERROR("IO is not a valid pbgz format, file meta magic no is %X, expect %X", 
                (*(uint32_t*)pReadBuffer), PBGZ_FILE_META_MAGIC);
            return -1;
        }
    }
    fileMeta.setBlockType(FILE_META);

    // read file meta length
    readLen = ioReader->readIO(pReadBuffer, PBGZ_FILE_META_SIZE_LENGTH);
    if (readLen != PBGZ_FILE_META_SIZE_LENGTH) {
        LOG_ERROR("Failed to read file meta length from IO.");
        return -1; // File read error       
    }

    const uint32_t metaLength = (uint32_t)(*(uint32_t*)(pReadBuffer));
    // read the the meta data 
    readLen = ioReader->readIO(pReadBuffer, metaLength);
    if (readLen != metaLength) {
        LOG_ERROR("Failed to read file meta data from IO.");
        return -1; // File read error   
    }
    coder_json fileMetaCoder;
    fileMetaCoder.decoder(pReadBuffer, metaLength, fileMeta.getMetaData());

    readLen = ioReader->readIO(pReadBuffer, PBGZ_FILE_META_CHECKSUM_LENGTH);
    if (readLen != PBGZ_FILE_META_CHECKSUM_LENGTH) {
        LOG_ERROR("Failed to read file meta data from IO.");
        return -1; // File read error   
    }
    fileMeta.setMetaChecksum(*(uint64_t*)pReadBuffer);
    return 0;
}

void PbgzFileReader::close() {
    if (ioReader) {
        ioReader->closeIO();
    }

    return;
}

PbgzFileHeader& PbgzFileReader::getFileHeader() {
    if (currentFileIndex < 0) {
        LOG_ERROR("No file has been read yet.");
        throw std::runtime_error("No file has been read yet.");
    }

    if (fileHeaderMap.find(currentFileIndex) == fileHeaderMap.end()) {
        LOG_ERROR("File header not found for current file index: %d", currentFileIndex);
        throw std::runtime_error("File header not found for current file index.");
    }
    
    return fileHeaderMap[currentFileIndex];
}

PbgzFileMeta& PbgzFileReader::getBaseFileMeta() {
    // if currentFileIndex is -1, it means no file has been read yet
    if (currentFileIndex < 0) {
        LOG_ERROR("No file has been read yet.");
        throw std::runtime_error("No file has been read yet.");
    }

    return baseFileMetaMap[currentFileIndex];
}

PbgzFileMeta& PbgzFileReader::getDynamicFileMeta() {
     // if currentFileIndex is -1, it means no file has been read yet
    if (currentFileIndex < 0) {
        LOG_ERROR("No file has been read yet.");
        throw std::runtime_error("No file has been read yet.");
    }

    return dynamicFileMetaMap[currentFileIndex];
}

int32_t PbgzFileReader::readDataBlock(PbgzDataBlock& dataBlock) {
    if (!ioReader) {
        LOG_ERROR("File is not open.");
        return -1; // File not open
    }   

    uint8_t* pReadBuffer = getFileReadBuffer();
    if (!pReadBuffer) {
        LOG_ERROR("Failed to get read buffer.");        
        return -1; // Memory allocation error
    }

    // Read the data block magic value, 4byte
    size_t readLen = ioReader->readIO(pReadBuffer, PBGZ_DATA_BLOCK_MAGIC_LENGTH);
    if (readLen == 0) { // Indicates read completion
        return 0;
    }

    if (readLen != PBGZ_DATA_BLOCK_MAGIC_LENGTH) {
        LOG_ERROR("Failed to read PBGZ data block header.");
        return -1; // File read error
    }

    if (memcmp(pReadBuffer, &PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH) != 0) {
        LOG_DEBUG("Not a magic value for PBGZ data block. %X", (uint32_t)(*(uint32_t*)pReadBuffer));
        // It maybe a new file, try to parse the file header and meta information
        if (memcmp(pReadBuffer, &PBGZ_FILE_MAGIC, PBGZ_FILE_MAGIC_LENGTH) == 0) {
            LOG_INFO("Detected a new PBGZ file format, reinitializing file header and meta.");
            // Since the file Magic value has already been read, no need to read again
            if (0 != initFileHeadAndMeta(true)) {
                LOG_ERROR("Failed to initialize file header and meta.");
                return -1; // Initialization error
            }

            // continue to read the data block from new file
            return readDataBlock(dataBlock); 
        } else if (memcmp(pReadBuffer, &PBGZ_FILE_META_MAGIC, PBGZ_FILE_MAGIC_LENGTH) == 0) {
            PbgzFileMeta tmpMete;
            if (0 != readFileMeta(tmpMete, false)) {
                LOG_INFO("Failed to read file dynamic file meta.");
            }
            return readDataBlock(dataBlock);
        }
        return -1; // Unexpect block, maybe meta, should ignore
    }

    dataBlock.setBlockType(FILE_DATA);

    // Read the meta length, 4byte
    readLen = ioReader->readIO(pReadBuffer, PBGZ_DATA_BLOCK_META_SIZE_LENGTH);
    if (readLen != PBGZ_DATA_BLOCK_META_SIZE_LENGTH) {
        LOG_ERROR("Failed to read PBGZ data block meta size.");
        return -1; // File read error           
    }
    uint32_t metaLength = (uint32_t)(*(uint32_t*)pReadBuffer);
    if (metaLength == 0) {
        LOG_ERROR("Meta length is zero, invalid data block.");
        return -1; // Invalid meta length       
    }
    // Read the meta data
    readLen = ioReader->readIO(pReadBuffer, metaLength);
    if (readLen != metaLength) {
        LOG_ERROR("Failed to read PBGZ data block meta data."); 
        return -1; // File read error   
    }
    coder_json blockMetaCoder;
    blockMetaCoder.decoder(pReadBuffer, metaLength, dataBlock.getMetaData());
   
    // Read the meta checksum, 8byte
    readLen = ioReader->readIO(pReadBuffer , PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH);
    if (readLen != PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH) {
        LOG_ERROR("Failed to read PBGZ data block meta checksum.");
        return -1; // File read error
    }
    dataBlock.setMetaCheckSum(*(uint64_t*)pReadBuffer);

    // Read the data length, 4byte
    readLen = ioReader->readIO(pReadBuffer, PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH);
    if (readLen != PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH) {
        LOG_ERROR("Failed to read PBGZ data block data length.");
        return -1; // File read error
    }   
    uint32_t dataLength = (uint32_t)(*(uint32_t*)(pReadBuffer));
    if (dataLength == 0) {
        LOG_ERROR("Data length is zero, invalid data block.");  
        return -1; // Invalid data length
    }
    dataBlock.setDataLength(dataLength);

    // Read the block data
    // PbgzDataBlock data storage is not memory allocated by itself, the address is passed in from external, reducing copying
    readLen = ioReader->readIO(dataBlock.getDataPtr(), dataLength);
    if (readLen != dataLength) {
        LOG_ERROR("Failed to read PBGZ data block data. expect %d, actual %d", dataLength, readLen);
        return -1; // File read error
    }

    // Read the data block checksum, 8byte
    readLen = ioReader->readIO(pReadBuffer, PBGZ_DATA_BLOCK_CHECKSUM_LENGTH);
    if (readLen != PBGZ_DATA_BLOCK_CHECKSUM_LENGTH) {
        LOG_ERROR("Failed to read PBGZ data block checksum.");
        return -1; // File read error
    }
    dataBlock.setDataCheckSum(*(uint64_t*)pReadBuffer);
    return 0;
}   


uint8_t* PbgzFileWriter::getFileWriteBuffer() {
    // 线程局部智能指针：自动管理内存，线程退出时析构释放
    thread_local static std::unique_ptr<uint8_t[]> pWriteBuffer;
    
    if (!pWriteBuffer) {
        // 用 unique_ptr 分配内存（无需手动 delete）
        pWriteBuffer = std::make_unique<uint8_t[]>(PBGZ_FILE_READ_BUFFER_LENGTH);
        if (!pWriteBuffer) { // make_unique 失败时返回空指针（C++17+）
            LOG_ERROR("Failed to allocate memory for write buffer.");
        }
    }

    // 返回智能指针指向的原始缓冲区地址（供业务使用）
    return pWriteBuffer.get();
}

int32_t PbgzFileWriter::open(){
    if (ioWriter == nullptr) {
        return -1;
    }
    if (0 != initFileHead()) {
        return -1; // Initialization error
    }
    return 0;
}

int32_t PbgzFileWriter::close(){
    return 0;
}

int32_t PbgzFileWriter::initFileHead() {
    if (!ioWriter) {
        LOG_ERROR("IO is not open.");
        return -1; // File not open
    }

    uint32_t writeLen = 0;
    // Serialize block type, 4 byte without '\0'
    writeLen += ioWriter->writeIO(PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC_LENGTH);
    // Serialize version
    writeLen += ioWriter->writeIO(fileHeader.getVerion(), PBGZ_FILE_VERSION_LENGTH);
    uint32_t fileHeadExtLength = sizeof(PbgzFileHeaderExt);
    writeLen += ioWriter->writeIO(&fileHeadExtLength, sizeof(uint32_t));
    writeLen += ioWriter->writeIO(&fileHeader.getFileHeaderExt(), fileHeadExtLength);
    if (writeLen > 0) {
        return 0;
    }
    return -1;
}

int32_t PbgzFileWriter::writeFileMeta(PbgzFileMeta& fileMeta){
    if (!ioWriter) {
        LOG_ERROR("IO is not open.");
        return -1; // File not open
    }

    coder_json fileMetaCoder;
    std::string fileMetaEncodeStr;
    fileMetaCoder.encoder(fileMeta.getMetaData(), fileMetaEncodeStr);
    uint32_t metaLength = fileMetaEncodeStr.length();

    uint32_t writeLen = 0;
    // write block type
    writeLen += ioWriter->writeIO(&PBGZ_FILE_META_MAGIC, PBGZ_FILE_META_MAGIC_LENGTH);

    // write meta data length
    writeLen += ioWriter->writeIO( &metaLength, PBGZ_FILE_META_SIZE_LENGTH);

    // write meta data
    writeLen += ioWriter->writeIO(fileMetaEncodeStr.c_str(), metaLength);

    // write checksum
    uint64_t checksum = fileMeta.getMetaChecksum();
    writeLen += ioWriter->writeIO(&checksum, PBGZ_FILE_META_CHECKSUM_LENGTH);

    return writeLen;
}

int32_t PbgzFileWriter::writeBlockData(PbgzDataBlock& dataBlock) {
    if (!ioWriter) {
        LOG_ERROR("IO is not open.");
        return -1; // File not open
    }

    uint8_t* pWriteBuffer = getFileWriteBuffer();
    if (!pWriteBuffer) {
        LOG_ERROR("Failed to get write buffer.");
        return -1; // Memory allocation error
    }

    uint32_t writeLen = 0;
    // Write block type
    writeLen += ioWriter->writeIO(&PBGZ_DATA_BLOCK_MAGIC, PBGZ_DATA_BLOCK_MAGIC_LENGTH);

    // Write meta length
    coder_json blockMetaCoder;
    std::string blockMetaOut;
    blockMetaCoder.encoder(dataBlock.getMetaData(), blockMetaOut);
    uint32_t dataMetaLength = blockMetaOut.length();
    writeLen += ioWriter->writeIO(&dataMetaLength, PBGZ_DATA_BLOCK_META_SIZE_LENGTH);
    
    // Write meta data
    writeLen += ioWriter->writeIO(blockMetaOut.c_str(), dataMetaLength);

    // Write meta checksum
    uint64_t checksum = dataBlock.getMetaCheckSum();
    writeLen += ioWriter->writeIO(&checksum, PBGZ_DATA_BLOCK_META_CHECKSUM_LENGTH);

    // Serialize data length
    uint32_t blockDataLength = dataBlock.getDataLength();
    writeLen += ioWriter->writeIO(&blockDataLength, PBGZ_DATA_BLOCK_DATA_SIZE_LENGTH);

    // Write block data
    void* pBlockData = dataBlock.getDataPtr();
    if (pBlockData && blockDataLength > 0) {
        writeLen += ioWriter->writeIO(pBlockData, blockDataLength);
    }

    // Write checksums
    uint64_t dataChecksum = dataBlock.getDataCheckSum();
    writeLen += ioWriter->writeIO(&dataChecksum, PBGZ_DATA_BLOCK_CHECKSUM_LENGTH);
    
    return writeLen;
}

void PbgzFileWriter::updateMetaOffset(uint64_t dynamicMetaOffset) {
    FileWriter* pFileWrite = dynamic_cast<FileWriter*>(ioWriter);
    if (pFileWrite == nullptr) {
        return;
    }

    LOG_INFO("Dynamic file offset = %llu", dynamicMetaOffset);

    pFileWrite->writeIOAt(11, &dynamicMetaOffset, sizeof(dynamicMetaOffset));
    return;
}
