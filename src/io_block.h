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

const size_t BLOCK_SIZE = 268435456;

/* Number of bytes the Creator reads ahead to determine the file format. This data is later merged into the first block / first package header, so nothing is wasted. */
const size_t BLOCK_TYPE_DETECT_SIZE = 1 << 20;

/* Initial allocation size for an input block: no longer allocates one large buffer for -l up front; starts at a fixed 1MB and grows on demand with realloc while reading. */
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
    }

    /*
     * Grows the buffer on demand. Returns 0 on success, -1 on failure (realloc
     * failed).
     *
     * Primary defense: the actuator's block entry pre-allocates block_size (the
     * upper bound determined at compression time) * 2, plugging all field
     * overflows. coder_io's putc check and the decode error return chain serve
     * as a fallback.
     *
     * Abnormal warning: when the target size exceeds the usual upper bound (2GB)
     * or growing to 16x the current size is still not enough—the former usually
     * means block_size is anomalous, the latter usually means a coder keeps
     * overflowing (the fallback chain also calls this)—both indicate a deeper
     * problem. It logs a warning but does not abort: the data may simply be
     * large, and the decision is left to the caller based on the error code.
     *
     * Called after reset() and before any data is written: at this point
     * dataLen == 0, the buffer holds no valid data, and there are no dangling
     * pointers into the buffer, so realloc is safe. After growing, bufferSize is
     * updated, and later getRemain/getCurrent/getBufferSize all see the new
     * value. When a block is returned to the freeOutputPool and reused, it
     * carries the larger buffer, and subsequent blocks keep using it without a
     * per-block realloc.
     */
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

    /*
     * Which pbgz package this block belongs to. It is also used as the identity
     * of an auxiliary block—under pipe input an absolute file offset degenerates
     * to 0 and is meaningless, so it cannot serve as identity.
     */
    int64_t getPackageIndex() const {
        return packageIndex;
    }

    void setPackageIndex(int64_t index) {
        packageIndex = index;
    }

    /*
     * Start position of this block's pbgz package in the file.
     *
     * The auxiliary-block addresses recorded in the block meta are **relative to
     * the package**—at compression time it is impossible to know where the file
     * will later be cat'ed to. The decompression side must add the start position
     * of the package this block belongs to in order to reconstruct absolute
     * addresses; this is the same semantics as getCurrentFileStart() + offset
     * used for the reference genome. In the single-package case this value is 0,
     * so relative degenerates to absolute.
     */
    int64_t getPackageStart() const {
        return packageStart;
    }

    void setPackageStart(int64_t start) {
        packageStart = start;
    }

    /*
     * Position-addressing flag: this block is a "synchronously emitted auxiliary
     * block"; the writer thread writes it immediately on sight, it does not take
     * part in blockId reordering, and it does not consume the data block id
     * sequence. Blocks come from a shared pool and are reused repeatedly, so this
     * flag must be cleared in reset(); otherwise a stale flag from the previous
     * round would cause an ordinary data block to be written out early as an
     * auxiliary block.
     */
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
};
