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
#include <algorithm>
#include <string>
#include <string.h>

#include "utils/memory_util.h"

const uint32_t BLOCK_SIZE = 268435456;

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
    QUAL_PRIOR = 1 << 8,            // Pre-trained QUAL model snapshot (256)
    PBGZFILE = UINT32_MAX           // Special value for PBGZ file type (4294967295)
} BlockType;


class RoughIOBlock /* Rough block */
{
public:
    RoughIOBlock(uint32_t len = BLOCK_SIZE) {
        blockSize = len;
        bufferSize = blockSize * 2; // Extra space to handle cases where compressed data is larger than original
        buffer =  MemoryUtil::safeAlloc<uint8_t>(bufferSize); //static_cast<uint8_t*>(calloc(bufferSize, sizeof(int8_t)));
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
    }

    uint8_t* getBuffer() const {
        return buffer;
    }

    uint32_t getBufferSize() {
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

    std::vector<uint32_t>& getNpos() {
        return npos;
    }

    void setBlockType(BlockType type) {
        blockType = type;
    }

    BlockType getBlockType() const {
        return blockType;
    }

    void setMaxLineLen(uint32_t len) {
        maxLineLen = len;
    }   

    uint32_t getMaxLineLen() const {
        return maxLineLen;
    }

    uint8_t* getCurrent() {
        return buffer + dataLen + metaLen;
    } 

    uint32_t getRemain() {
        return bufferSize - dataLen - metaLen;
    }

    uint8_t* getMetaBuffer() {
        return buffer + dataLen;
    }

    uint32_t getMetaLen(){   
        return metaLen;
    }

    void setMetaLen(uint32_t len){
        metaLen = len;
    }

    uint32_t  getTotalDataLen() {
        return dataLen + metaLen;
    }

    uint32_t getBlockSize() {
        return blockSize;
    }
    
private:
    uint8_t *buffer;
    uint32_t bufferSize;             /* Size of buffer */
    BlockType blockType;             /* Current block type */
    uint32_t blockSize;              /* Block size */
    std::vector<uint32_t> npos;      /* Positions of newline characters in buffer, starting from 0 */
    int64_t blockId;
    int64_t dataLen;
    uint32_t maxLineLen;
    uint32_t metaLen;
};
