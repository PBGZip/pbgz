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
    PBGZFILE = UINT32_MAX           // Special value for PBGZ file type (4294967295)
} BlockType;


class RoughIOBlock /* 粗糙的块*/
{
public:
    RoughIOBlock(uint32_t len = BLOCK_SIZE) {
        blockSize = len;
        bufferSize = blockSize * 2; // 多部分空间, 用于解决压缩后变大的问题
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
    uint32_t bufferSize;             /* buffer的size */
    BlockType blockType;             /* 当前块对应块类型 */
    uint32_t blockSize;              /* 块大小 */
    std::vector<uint32_t> npos;      /* buffer中换行符的位置，从0开始 */
    int64_t blockId;
    int64_t dataLen;
    uint32_t maxLineLen;
    uint32_t metaLen;
};


