/*
 * io_block.h - Header file for pbgz project
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
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <memory>
#include <string>

#include "utils/memory_util.h"
#include "bam_columns.h"

/* One BAM reference-sequence entry (name + length) from the BAM header. */
struct BamRef {
    std::string name;
    int32_t len;
};

const size_t BLOCK_SIZE = 268435456;

/* Bytes read ahead to detect the file format. */
const size_t BLOCK_TYPE_DETECT_SIZE = 1 << 20;

/* Initial input block allocation size. */
const size_t FIXED_INPUT_BLOCK_SIZE = 1 << 20;

typedef enum
{
    TYPE_UNKNOW = 0,                // Unknown block type (0)
    GZIP = (1 << 0),                // GZIP compressed block (1)
    BINARY = (1 << 1),              // Binary block (2)
    FASTQ_GEN2 = (1 << 2),          // FASTQ generation 2 block (4)
    FASTQ_GEN3 = (1 << 3),          // FASTQ generation 3 block (8)
    BINARY_GZIP = (BINARY | GZIP),  // Binary block with GZIP compression (3)
    FASTQ_GEN2_GZIP = (FASTQ_GEN2 | GZIP), // FASTQ gen2 block with GZIP compression (5)
    FASTQ_GEN3_GZIP = (FASTQ_GEN3 | GZIP), // FASTQ gen3 block with GZIP compression (9)
    BAM = (1 << 4),                 // BAM format block (16)
    REFERENCE = (1 << 5),           // 
    REFERENCE_INDEX = 1 << 6,       //
    SAM = 1 << 7, 
    SAM_GZIP = (SAM | GZIP),        // SAM block with GZIP compression (129)
    QUAL_PRIOR = 1 << 8,            // Pre-trained QUAL model snapshot (256)
    PBGZFILE = UINT32_MAX           // Special value for PBGZ file type (4294967295)
} BlockType;


class RoughIOBlock /* Rough block */
{
public:
    RoughIOBlock(size_t len = BLOCK_SIZE) {
        blockSize = len;
        bufferSize = blockSize * 2; // Extra space to handle cases where compressed data is larger than original
        buffer =  MemoryUtil::safeAlloc<uint8_t>(bufferSize);
        reset();
    }

    ~RoughIOBlock() {
        MemoryUtil::safeFree(buffer);
    }

    void reset() {
        blockId = -1;
        maxLineLen = 0;
        blockType = TYPE_UNKNOW;
        npos.clear();
        dataLen = 0;
        metaLen = 0;
        syncAux = false;
        packageStart = 0;
        packageIndex = -1;
        md5Hex.clear();
        md5Ready = false;
        bamCols.reset();
        bamPrebuilt = false;
        rawBam = false;
        bamRefs.reset();
    }

    /* Grow the buffer on demand; returns 0 on success, -1 on realloc failure. */
    int32_t ensureCapacity(size_t neededSize) {
        if (bufferSize >= neededSize) {
            return 0;
        }
        size_t newSize = neededSize;
        if (newSize < bufferSize * 2) {
            newSize = bufferSize * 2;
        }
        if (newSize > ((size_t)2 << 30) || newSize > bufferSize * 16) {
            LOG_WARNING("ensureCapacity: target %zu bytes (have %zu), abnormally large; "
                        "block_size may be wrong or coder keeps overflowing", newSize, bufferSize);
        }
        uint8_t* newBuffer = static_cast<uint8_t*>(realloc(buffer, newSize));
        if (newBuffer == nullptr) {
            LOG_ERROR("ensureCapacity: realloc failed, need %zu, have %zu", neededSize, bufferSize);
            return -1;
        }
        buffer = newBuffer;
        bufferSize = newSize;
        return 0;
    }

    /* Which pbgz package this block belongs to; also the identity of an auxiliary block. */
    int64_t getPackageIndex() const {
        return packageIndex;
    }

    void setPackageIndex(int64_t index) {
        packageIndex = index;
    }

    /* Start position of this block's pbgz package in the file. */
    int64_t getPackageStart() const {
        return packageStart;
    }

    void setPackageStart(int64_t start) {
        packageStart = start;
    }

    /* True for synchronously emitted auxiliary blocks; cleared by reset(). */
    bool isSyncAux() const {
        return syncAux;
    }

    void setSyncAux(bool flag) {
        syncAux = flag;
    }

    uint8_t* getBuffer() const {
        return buffer;
    }

    size_t getBufferSize() {
        return bufferSize;
    }

    void setBlockId(int64_t id) {
        blockId = id;
    }

    int64_t getBlockId() const {
        return blockId;
    }

    void setDataLen(int64_t len) {
        dataLen = len;
    }

    int64_t getDataLen() const {
        return dataLen;
    }

    std::vector<size_t>& getNpos() {
        return npos;
    }

    void setBlockType(BlockType type) {
        blockType = type;
    }

    BlockType getBlockType() const {
        return blockType;
    }

    void setMaxLineLen(size_t len) {
        maxLineLen = len;
    }   

    size_t getMaxLineLen() const {
        return maxLineLen;
    }

    uint8_t* getCurrent() {
        return buffer + dataLen + metaLen;
    } 

    size_t getRemain() {
        return bufferSize - dataLen - metaLen;
    }

    uint8_t* getMetaBuffer() {
        return buffer + dataLen;
    }

    size_t getMetaLen(){   
        return metaLen;
    }

    void setMetaLen(size_t len){
        metaLen = len;
    }

    size_t  getTotalDataLen() {
        return dataLen + metaLen;
    }

    size_t getBlockSize() const {
        return blockSize;
    }

    /* MD5 of the raw block text, computed by the reader thread. */
    bool hasMd5() const {
        return md5Ready;
    }

    std::string getMd5() const {
        return md5Hex;
    }

    void setMd5(const std::string& value) {
        md5Hex = value;
        md5Ready = true;
    }

    /* Structured BAM columns carried alongside the raw record bytes. */
    std::shared_ptr<BamColumns> getBamColumns() const {
        return bamCols;
    }

    void setBamColumns(const std::shared_ptr<BamColumns>& c) {
        bamCols = c;
    }

    /* True when the payload already is the BAM byte stream of its records. */
    bool isBamPrebuilt() const {
        return bamPrebuilt;
    }

    void setBamPrebuilt(bool value) {
        bamPrebuilt = value;
    }

    /* Archive-mode BAM block: payload is raw records, worker renders to SAM text. */
    bool isRawBam() const {
        return rawBam;
    }

    void setRawBam(bool value) {
        rawBam = value;
    }

    /* Shared read-only BAM reference-sequence table for worker-side rendering. */
    std::shared_ptr<const std::vector<BamRef>> getBamRefs() const {
        return bamRefs;
    }

    void setBamRefs(const std::shared_ptr<const std::vector<BamRef>>& refs) {
        bamRefs = refs;
    }

private:
    uint8_t *buffer;
    size_t bufferSize;             /* Size of buffer */
    BlockType blockType;             /* Current block type */
    size_t blockSize;              /* Block size */
    std::vector<size_t> npos;      /* Positions of newline characters in buffer, starting from 0 */
    int64_t blockId;
    bool syncAux = false;
    int64_t packageStart = 0;
    int64_t packageIndex = -1;
    int64_t dataLen;
    size_t maxLineLen;
    size_t metaLen;
    std::string md5Hex;            /* Cached MD5 of the raw block, filled on the read side */
    bool md5Ready = false;
    std::shared_ptr<BamColumns> bamCols;
    bool bamPrebuilt = false;      /* payload already is a BAM byte stream (see setBamPrebuilt) */
    bool rawBam = false;           /* archive-mode block carries raw BAM records, worker renders to SAM text */
    std::shared_ptr<const std::vector<BamRef>> bamRefs; /* BAM ref table for worker-side rendering */
};
