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

/* Creator 预读用于判定文件格式的字节数。该数据随后并入首个块/首个包头部，不浪费。 */
const size_t BLOCK_TYPE_DETECT_SIZE = 1 << 20;

/* 输入块初始分配大小：不再按 -l 一次性分配大缓冲，固定 1MB 起步，读取时按需 realloc 扩容。 */
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
     * 按需扩容。返回 0 成功，-1 失败（realloc 失败）。
     *
     * 主防线：actuator 块入口按 block_size（压缩时确定的上界）×2 预分配，堵所有
     * 字段越界。coder_io 的 putc 检查与 decode 错误返回链作兜底。
     *
     * 异常告警：目标 size 超过常规上界（2GB）或相对当前已扩到 16 倍还不够——
     * 前者多半是 block_size 异常，后者多半是 coder 反复越界（兜底链也调它），
     * 都说明有更深层问题。打警告但不中断：数据可能就是大，留给上层据错误码决策。
     *
     * 在 reset() 之后、任何数据写入之前调用：此刻 dataLen == 0，缓冲里没有有效
     * 数据，不存在指向缓冲内部的悬空指针，realloc 安全。扩容后 bufferSize 更新，
     * 后续 getRemain/getCurrent/getBufferSize 都看到新值。块归还 freeOutputPool
     * 后复用带着更大的缓冲，后续块沿用，不会每块 realloc。
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
     * 本块属于第几个 pbgz 包。辅助块的身份用的就是它——
     * 绝对文件偏移在管道输入下退化为 0，根本不成立，不能当身份。
     */
    int64_t getPackageIndex() const {
        return packageIndex;
    }

    void setPackageIndex(int64_t index) {
        packageIndex = index;
    }

    /*
     * 本块所属 pbgz 包在文件中的起点。
     *
     * 块 meta 里记录的辅助块地址是**包内相对**的——压缩时无从预知自己将来被 cat 到
     * 哪个位置。解压侧必须加上本块所属包的起点才还原成绝对地址，与参考基因组走的
     * getCurrentFileStart() + offset 是同一套语义。单包场景该值为 0，退化成相对即绝对。
     */
    int64_t getPackageStart() const {
        return packageStart;
    }

    void setPackageStart(int64_t start) {
        packageStart = start;
    }

    /*
     * 位置寻址标记：本块是"同步发射的辅助块"，写线程见到即写，
     * 不参与 blockId 顺序重排，也不消耗数据块的 id 序列。
     * 块取自公共池并被反复复用，故必须在 reset() 里清掉，
     * 否则上一轮的标记会让一个普通数据块被当成辅助块提前写出。
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

    size_t getBlockSize() {
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
