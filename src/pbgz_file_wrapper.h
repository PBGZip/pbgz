/*
 * pbgz_file_handler.h - Header file for pbgz project
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

#include <stdint.h>

#include "pbgz_file.h"
#include "io_wrapper.h"
#include "io_block.h"

class RoughIOBlock;


/// @brief Reader for pbgz format files
class PbgzFileReader {
public:
    int32_t open();

    PbgzFileHeader& getFileHeader();

    PbgzFileMeta& getBaseFileMeta();

    PbgzFileMeta& getDynamicFileMeta();

    /* Prefetched format-probing bytes, consumed in order before regular reads. */
    void setPreReadData(const uint8_t* data, size_t len);

    /* Start byte offset of the currently read pbgz package in the input stream. */
    uint64_t getCurrentFileStart() const;

    /* Absolute start of the block most recently read by readDataBlock. */
    uint64_t getCurrentBlockStart() const { return currentBlockStart; }

    /* Index of the pbgz package currently being read. */
    int32_t getCurrentFileIndex() const { return currentFileIndex; }

    /* Read one data block; if dst is non-null, resize it and read the data into it. */
    int32_t readDataBlock(PbgzDataBlock& dataBlock, RoughIOBlock* dst = nullptr);

    PbgzFileReader(IOReader* pReader) : ioReader(pReader) {
        currentFileIndex = -1;  // -1 indicates no file has been read yet
        currentBlockStart = 0;
        preReadSize = 0;
        preReadPos = 0;
    }

    void close();

    virtual ~PbgzFileReader() { }

private:
    uint64_t currentBlockStart;

    int32_t initFileHeadAndMeta(bool isCheckMagic = false);

    static uint8_t* getFileReadBuffer();

    int32_t readFileMeta(PbgzFileMeta& fileMeta, bool isCheckMagic = true, bool usePreRead = true);

    /* Consume the prefetch buffer first, then read directly from ioReader. */
    size_t readFromSource(void* pBuffer, size_t n);

private:
    std::map<int32_t, PbgzFileHeader> fileHeaderMap;  // File headers; a file may hold several concatenated packages
    std::map<int32_t, PbgzFileMeta> baseFileMetaMap;      // Base file metadata
    std::map<int32_t, PbgzFileMeta> dynamicFileMetaMap;   // Dynamic file metadata (reader side only)
    std::map<int32_t, uint64_t> fileStartMap;             // Start position of each concatenated package, for offset translation
    int32_t currentFileIndex; // Current file sequence number
    IOReader* ioReader;

    /* Prefetched format-probing bytes, consumed in order first when reading back */
    uint8_t preReadBuf[BLOCK_TYPE_DETECT_SIZE];
    size_t preReadSize;   // prefetch bytes not yet consumed
    size_t preReadPos;    // Position within the prefetch buffer
};

/// @brief Writer for pbgz format files
/// @note: 
class PbgzFileWriter {
public:
    int32_t open();

    int32_t close();

    int32_t writeFileMeta(PbgzFileMeta& fileMeta);

    int32_t writeBaseFileMeta() {
        return writeFileMeta(baseFileMeta);
    }

    int32_t writeDynamicFileMeta() {
        return writeFileMeta(dynamicFileMeta);
    }

    int32_t writeBlockData(PbgzDataBlock& dataBlock);

    int32_t setBaseFileMeta(const PbgzFileMeta& metaInfo) {
        baseFileMeta = metaInfo;
        return 0;
    }

     int32_t setDynamicFileMeta(const PbgzFileMeta& metaInfo) {
        dynamicFileMeta = metaInfo;
        return 0;
    }

    PbgzFileWriter(IOWriter* pWriter ) : ioWriter(pWriter) {}

    virtual ~PbgzFileWriter() { 
        ioWriter = nullptr;
    }

    PbgzFileHeader& getFileHeader() {
        return fileHeader;
    }

    PbgzFileMeta& getBaseFileMeta() {
        return baseFileMeta;
    }

    PbgzFileMeta& getDynamicFileMeta() {
        return dynamicFileMeta;
    }

    void updateMetaOffset(uint64_t dynamicMetaOffset);

private:
    static uint8_t* getFileWriteBuffer();

    int32_t initFileHead();

private:
    PbgzFileHeader fileHeader; // File header
    PbgzFileMeta baseFileMeta; // Base file meta
    PbgzFileMeta dynamicFileMeta; // Dynamic file meta 
    IOWriter* ioWriter;
};
