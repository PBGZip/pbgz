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
#include "io_block.h"

class RoughIOBlock;


/// @brief Reader for pbgz format files
class PbgzFileReader {
public:
    int32_t open();

    PbgzFileHeader& getFileHeader();

    PbgzFileMeta& getBaseFileMeta();

    PbgzFileMeta& getDynamicFileMeta();

    /*
     * Prefetched format-probing data (bytes BlockFactory read from ioReader before
     * creating the reader).
     *
     * Pipe input is not seekable, so the bytes consumed by probing must be handed back
     * to the reader through this interface or the file header would be lost. These
     * bytes are consumed in order ahead of regular reads; seek scenarios (dynamic meta
     * reread) read ioReader directly, bypassing the prefetch buffer, and the prefetch
     * buffer remains in effect once the seek position is restored.
     */
    void setPreReadData(const uint8_t* data, size_t len);

    /*
     * Start byte position, in the whole input stream, of the pbgz package currently
     * being read.
     *
     * It exists because of cat concatenation: once several independently compressed
     * packages are joined end-to-end into one file, from the second package onward no
     * offset recorded inside a package equals its real position in the concatenated
     * file. During compression output always starts writing at 0, so an offset stored
     * in a package is essentially "the distance from its own header"; when reading a
     * single package the start happens to be 0, relative and absolute values coincide,
     * and the problem is masked - only concatenation exposes it.
     *
     * The read side simply takes this start and reconstructs the real position as
     * "package start + in-package offset"; the written format needs no byte changed,
     * and compatibility with existing files is unaffected.
     */
    uint64_t getCurrentFileStart() const;

    /*
     * Absolute start of the block most recently read by readDataBlock (the byte where
     * the magic number sits).
     *
     * An auxiliary block's identity is its absolute address: the package start only
     * suffices while "exactly one of each kind per package" holds, and breaks once
     * priors are sharded. The cache keys on this address, so re-encountering the same
     * block hits directly without re-reading.
     */
    uint64_t getCurrentBlockStart() const { return currentBlockStart; }

    /*
     * Which pbgz package is currently being read (increments in order after cat
     * concatenation). Unlike getCurrentFileStart, this sequence number works under
     * pipe input too - it depends only on the number of package headers parsed, not on
     * any file position.
     */
    int32_t getCurrentFileIndex() const { return currentFileIndex; }

    /*
     * Read one data block. When dst is non-null, the block data is first resized into
     * dst to dataLength (the input block buffer is preallocated at the parameter-level
     * block_size, but a file's block may be larger, e.g. a -l 9 512MB block still
     * reaches ~72MB compressed), then the data is read into dst's buffer.
     */
    int32_t readDataBlock(PbgzDataBlock& dataBlock, RoughIOBlock* dst = nullptr);

    PbgzFileReader(IOReader* pReader) : ioReader(pReader) {
        currentFileIndex = -1;  // Initialize file index to -1, indicating no file has been read yet
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

    /* Consume the prefetch buffer first, then read ioReader; dynamic meta rereads (already seeked) go through direct ioReader reads */
    size_t readFromSource(void* pBuffer, size_t n);

private:
    std::map<int32_t, PbgzFileHeader> fileHeaderMap;  // File headers, a file may consist of multiple compressed packages concatenated together
    std::map<int32_t, PbgzFileMeta> baseFileMetaMap;      // Base file metadata
    std::map<int32_t, PbgzFileMeta> dynamicFileMetaMap;   // Dynamic file meta data, only exits in file reader
    std::map<int32_t, uint64_t> fileStartMap;             // start position of each concatenated package in the input stream, used to translate in-package relative offsets back into real positions
    int32_t currentFileIndex; // Current file sequence number, indicating which file is currently being read
    IOReader* ioReader;

    /* Prefetched format-probing data: bytes already consumed from ioReader by the Creator; consumed in order first when reading back */
    uint8_t preReadBuf[BLOCK_TYPE_DETECT_SIZE];
    size_t preReadSize;   // prefetch bytes not yet consumed
    size_t preReadPos;    // consumption position within the prefetch buffer
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
