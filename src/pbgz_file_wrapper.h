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
#include <string>

#include "pbgz_file.h"
#include "io_wrapper.h"


/// @brief Reader for pbgz format files
class PbgzFileReader {
public:
    int32_t open();

    const PbgzFileHeader& getFileHeader();

    const PbgzFileMeta& getFileMeta();

    int32_t readDataBlock(PbgzDataBlock& dataBlock);

    PbgzFileReader(IOReader* pReader) : ioReader(pReader) {
        currentFileIndex = -1;  // Initialize file index to -1, indicating no file has been read yet
        // Implement the read logic for PBGZ file format
        if (0 != initFileHeadAndMeta()) {
            throw std::runtime_error("Create PbgzFileReader, load head and meta failed");
        }
    }

    void close();

    virtual ~PbgzFileReader() { }

private:
    int32_t initFileHeadAndMeta(bool isCheckMagic = false);

    static uint8_t* getFileReadBuffer();

private:
    std::map<int32_t, PbgzFileHeader> fileHeaderMap;  // File headers, a file may consist of multiple compressed packages concatenated together
    std::map<int32_t, PbgzFileMeta> fileMetaMap;      // File metadata
    int32_t currentFileIndex; // Current file sequence number, indicating which file is currently being read
    IOReader* ioReader;
};

/// @brief Writer for pbgz format files
/// @note: 
class PbgzFileWriter {
public:
    int32_t open();

    int32_t close();

    int32_t writeFileMeta();

    int32_t writeBlockData(PbgzDataBlock& dataBlock);

    int32_t setFileMeta(const PbgzFileMeta& metaInfo) {
        fileMeta = metaInfo;
        return 0;
    }

    PbgzFileWriter(IOWriter* pWriter ) : ioWriter(pWriter) {}

    virtual ~PbgzFileWriter() { 
        ioWriter = nullptr;
    }

    PbgzFileHeader& getFileHeader() {
        return fileHeader;
    }

    PbgzFileMeta& getFileMeta() {
        return fileMeta;
    }

private:
    static uint8_t* getFileWriteBuffer();

    int32_t initFileHead();

private:
    PbgzFileHeader fileHeader; // File header
    PbgzFileMeta fileMeta; // File metadata
    IOWriter* ioWriter;
};
