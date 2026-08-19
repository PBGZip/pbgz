/*
 * block_wrapper.cpp - Source file for pbgz project
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

#include "block_wrapper.h"
#include "log/logger.h"
#include "pbgz_types.h"
#include <cstring>
#include <set>
#include "utils/md5_util.h"
#include "pbgz_manager.h"
#include "config_manager.h"
#include <zlib.h>

#if defined(__x86_64__)
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

BlockType BlockReader::analyzeBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    return BINARY;
}

int64_t BlockReader::readBlock(RoughIOBlock* blockPtr, BlockType fileType) {
    if (ioReader == nullptr || blockPtr == nullptr) {
        LOG_ERROR("IO reader or block pointer is null.");
        return -1;
    }

    uint8_t* buffer = blockPtr->getBuffer();
    size_t target = readTargetBytes(blockPtr);
    size_t totalLen = 0;

    blockPtr->setBlockId(blockId++);

    prependBufferedData(blockPtr, totalLen);

    /* 输入块按需扩容：初始只分配 1MB，读取目标（-l）更大时 realloc */
    buffer = blockPtr->getBuffer();
    if (blockPtr->getBufferSize() < target) {
        if (0 != blockPtr->ensureCapacity(target)) {
            LOG_ERROR("Ensure input block capacity failed.");
            return -1;
        }
        buffer = blockPtr->getBuffer();
    }

    while (totalLen < target) {
        size_t readLen = ioReader->readIO(buffer + totalLen, target - totalLen);
        if (readLen == 0) { 
            break; // EOF
        }
        totalLen += readLen;
    }

    if (totalLen == 0) {
        return 0; // EOF
    }

    blockPtr->setDataLen(static_cast<int64_t>(totalLen));

    // 子类完成格式预分析（换行符位置记录、Fastq/Sam actuator 预分析等）
    BlockType type = analyzeBlock(blockPtr, fileType);
    blockPtr->setBlockType(type);

    LOG_DEBUG("Read One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    
    return blockPtr->getDataLen();
}

size_t BlockReader::prependBufferedData(RoughIOBlock* blockPtr, size_t& totalLen, bool capAtBlockSize) {
    uint8_t* buffer = blockPtr->getBuffer();
    /* 封顶按读取目标（-l）而非块自身 blockSize（输入块只分配 1M，读取目标可能更大/更小） */
    const size_t capSize = capAtBlockSize ? readTargetBytes(blockPtr) : (size_t)-1;
    size_t placed = 0;

    /*
     * 上一块对齐遗留的 cache 先并入：它在文件里的位置先于剩余的预读数据（detect）。
     * 若先并入 detect，当 detect 超过读取目标（如 -l 1 的 512KB < 1MB 预读）时，
     * 上一块拆到 cache 的跨块记录会被丢在块边界后面，造成数据错乱。
     */
    if (cacheLen > 0 && (!capAtBlockSize || totalLen < capSize)) {
        size_t toCopy = cacheLen;
        if (capAtBlockSize && toCopy > capSize - totalLen) {
            toCopy = capSize - totalLen;
        }
        if (totalLen + toCopy > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(totalLen + toCopy)) {
                return placed;
            }
            buffer = blockPtr->getBuffer();
        }
        memcpy(buffer + totalLen, cache, toCopy);
        totalLen += toCopy;
        placed += toCopy;
        cacheLen -= toCopy;
        if (cacheLen > 0) {
            memmove(cache, cache + toCopy, cacheLen);
        }
    }

    // Creator 预读的格式探测数据并入（排在 cache 之后）
    if (detectLen > 0) {
        size_t toCopy = detectLen;
        if (capAtBlockSize && toCopy > capSize - totalLen) {
            toCopy = capSize - totalLen;
        }
        if (toCopy > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(toCopy)) {
                return placed;
            }
            buffer = blockPtr->getBuffer();
        }
        memcpy(buffer + totalLen, detectBuf, toCopy);
        totalLen += toCopy;
        placed += toCopy;
        if (toCopy < detectLen) {
            memmove(detectBuf, detectBuf + toCopy, detectLen - toCopy);
            detectLen -= toCopy;
        } else {
            detectLen = 0;
        }
    }
    return placed;
}

bool BlockReader::alignToRecordBoundary(RoughIOBlock* blockPtr, bool isFastq) {
    const int64_t totalLen = blockPtr->getDataLen();
    std::vector<size_t>& npos = blockPtr->getNpos();
    const int32_t lineNum = static_cast<int32_t>(npos.size());
    uint8_t* buffer = blockPtr->getBuffer();
    const size_t readTarget = readTargetBytes(blockPtr);

    int64_t remainLen = 0;
    if (lineNum > 0 && totalLen >= (int64_t)readTarget) {
        if (isFastq) {
            const int32_t completeLines = (lineNum >> 2) << 2;
            if (completeLines == 0) {
                // 整块连一条完整 FASTQ 记录都不够，无法对齐
                return false;
            }
            remainLen = totalLen - npos[completeLines - 1] - 1;
        } else {
            remainLen = totalLen - npos[lineNum - 1] - 1;
        }
    }

    if (remainLen > 0) {
        // If there's still data in cache, do memmove first
        if (cacheLen > 0) {
            memmove(cache + remainLen, cache, cacheLen);
        }
        // Put remaining data into cache
        memcpy(cache, buffer + totalLen - remainLen, remainLen);
        cacheLen = cacheLen + remainLen;
        const int64_t newTotalLen = totalLen - remainLen;
        blockPtr->setDataLen(newTotalLen);
        if (isFastq) {
            for (int t = 0; t < lineNum - ((lineNum >> 2) << 2); ++t) {
                npos.pop_back(); // Remove last newline position
            }
        }
        return true;
    }

    /*
     * 末尾块（totalLen < 读取目标，EOF）：结构化编码要求块恰好落在记录
     * 边界上——SAM 以 \n 收尾，FASTQ 还需行数为 4 的倍数（完整记录）。
     * 不满足即数据不合规（缺结尾换行、记录被截断），绝不能篡改数据去凑
     * 格式（旧实现会补一个 \n，还原后多出一字节），而是把本块降级为
     * BINARY 通用压缩：压缩率让位于忠实还原。
     */
    bool clean = false;
    if (lineNum > 0 && npos[lineNum - 1] == (size_t)(totalLen - 1)) {
        clean = !isFastq || (lineNum % 4 == 0);
    }
    if (!clean) {
        LOG_INFO("Tail block not aligned to record boundary, fallback to binary codec (blockId=%ld).",
                 (long)blockPtr->getBlockId());
        return false;
    }
    return true;
}

BlockType FastqBlockReader::analyzeBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    uint8_t* buffer = blockPtr->getBuffer();
    const int64_t dataLen = blockPtr->getDataLen();

    if (dataLen == 0 || buffer[0] != '@') {
        return BINARY;
    }

    // Build lookup table for O(1) base character validation
    static bool isValidBase[256];
    static bool lookupTableInitialized = false;
    if (!lookupTableInitialized) {
        for (int i = 0; i < 256; i++) {
            isValidBase[i] = false;
        }
        isValidBase[(uint8_t)'A'] = true;
        isValidBase[(uint8_t)'C'] = true;
        isValidBase[(uint8_t)'G'] = true;
        isValidBase[(uint8_t)'N'] = true;
        isValidBase[(uint8_t)'T'] = true;
        isValidBase[(uint8_t)'a'] = true;
        isValidBase[(uint8_t)'c'] = true;
        isValidBase[(uint8_t)'g'] = true;
        isValidBase[(uint8_t)'n'] = true;
        isValidBase[(uint8_t)'t'] = true;
        lookupTableInitialized = true;
    }

    // Prefetch optimization for large FASTQ blocks
    if (dataLen >= 4096) {
        for (int64_t i = 0; i < dataLen; i += 64) {
#if defined(__x86_64__)
            _mm_prefetch((char*)(buffer + i + 256), _MM_HINT_T0);
#elif defined(__aarch64__)
            __builtin_prefetch(buffer + i + 256, 0, 3);
#endif
        }
    }

    int64_t lineNum = 0;        // Line count
    int32_t baseLen = 0;        // Base length
    int32_t maxBaseLen = 0;     // Maximum base length
    int64_t lastEndlinePos = 0; // Previous newline position
    int64_t endlinePos = 0;     // Current newline position
    int64_t linePos = 0;        // Position within current line
    const char* bufPtr = reinterpret_cast<const char*>(buffer);

    // Use memchr to find newlines efficiently, then do per-character validation
    const char* searchStart = bufPtr;
    while (true) {
        const char* newlinePtr = static_cast<const char*>(memchr(searchStart, '\n',
                                                      dataLen - (searchStart - bufPtr)));
        if (newlinePtr == nullptr) {
            break;
        }

        // Validate characters between last newline and current newline
        int64_t newlinePos = newlinePtr - bufPtr;
        int lineMod = lineNum % 4;

        for (int64_t pos = (searchStart - bufPtr); pos < newlinePos; ++pos) {
            if (lineMod == 0) {
                if (linePos == 0 && buffer[pos] != '@') {
                    return BINARY;
                }
            } else if (lineMod == 1) {
                if (!isValidBase[static_cast<uint8_t>(buffer[pos])]) {
                    return BINARY;
                }
            } else if (lineMod == 2) {
                if (linePos == 0 && buffer[pos] != '+') {
                    return BINARY;
                }
            }
            ++linePos;
        }

        // Process newline
        lastEndlinePos = endlinePos;
        endlinePos = newlinePos;

        if (lineMod == 1) {
            baseLen = endlinePos - lastEndlinePos - 1;
            if (baseLen > maxBaseLen) {
                maxBaseLen = baseLen;
            }
        } else if (lineMod == 3) {
            if (endlinePos - lastEndlinePos - 1 != baseLen) {
                return BINARY;
            }
        }

        ++lineNum;
        blockPtr->getNpos().push_back(static_cast<uint32_t>(newlinePos));
        linePos = 0;
        searchStart = newlinePtr + 1;
    }

    if (maxBaseLen == 0) {
        return BINARY;
    }

    blockPtr->setMaxLineLen(maxBaseLen);
    BlockType baseType = (maxBaseLen > GENE2_MAX_BASE) ? FASTQ_GEN3 : FASTQ_GEN2;

    // 块尾记录对齐（把不完整记录移入缓存）
    if (!alignToRecordBoundary(blockPtr, true)) {
        return BINARY;
    }

    /*
     * 预分析在 coder 侧 actuator 中执行并持有其状态（压缩需要），这里只做换行符位置
     * 记录与结构校验，不重复运行 FastqCodecActuator::preAnalysis。
     * GZ 输入时块类型仍按非 GZ 的 FASTQ 类型设置，由 io 层透明解压负责还原。
     */
    return baseType;
}

/*
 * SAM/SAM-GZ 块读取：
 *   1. 头部（@ 行）独立成一个块返回；
 *   2. 数据区按 readsPerBlock 条 read 分块（由压缩级别决定 10000/25000/100000）。
 * 块大小由 read 行数决定，不再受 byte blockSize 约束，缓冲不足时按需扩容。
 */
int64_t SamBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType fileType) {
    if (ioReader == nullptr || blockPtr == nullptr) {
        LOG_ERROR("IO reader or block pointer is null.");
        return -1;
    }

    blockPtr->setBlockId(blockId++);
    std::vector<size_t>& npos = blockPtr->getNpos();
    npos.clear();

    size_t totalLen = 0;
    /* SAM 块按 read 行数分块，不受 byte blockSize 约束，一次性并入全部预读数据与缓存 */
    prependBufferedData(blockPtr, totalLen, false);

    uint8_t* buffer = blockPtr->getBuffer();
    const size_t chunkSize = 1 << 20;

    /* 缓冲为空（无预读数据，如测试直接构造）时先读一段，才能判断是否以 @ 头部开头 */
    if (totalLen == 0) {
        if (totalLen + chunkSize > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(totalLen + chunkSize)) {
                return -1;
            }
            buffer = blockPtr->getBuffer();
        }
        const size_t readLen = ioReader->readIO(buffer, chunkSize);
        if (readLen == 0) {
            return 0;   // EOF
        }
        totalLen = readLen;
    }

    const bool startsWithHeader = (totalLen > 0 && buffer[0] == '@');
    /* 只有"纯头部块"不含数据行；其余块（数据块、降级块、未拆分头部的块）都携带数据 */
    lastBlockHasData = !(startsWithHeader && splitHeader);

    /*
     * splitHeader 且块以 @ 开头（头部块）：读到第一条数据行即停，本块只返回头部行。
     * 否则：读到 readsPerBlock 条数据行（或 EOF）为止，头部行并入首个数据块。
     */
    const size_t readTarget = (splitHeader && startsWithHeader) ? 1 : readsPerBlock;

    size_t dataLineCount = 0;
    size_t lineStart = 0;
    size_t scanPos = 0;

    while (dataLineCount < readTarget) {
        // 先扫描当前缓冲里已有的字节
        size_t i = scanPos;
        while (i < totalLen) {
            if (buffer[i] == '\n') {
                if (buffer[lineStart] != '@') {
                    ++dataLineCount;
                }
                npos.push_back(i);
                lineStart = i + 1;
                if (dataLineCount >= readTarget) {
                    break;
                }
            }
            ++i;
        }
        scanPos = i;
        if (dataLineCount >= readTarget) {
            break;
        }

        // 需要更多数据：按需扩容后继续读
        if (totalLen + chunkSize > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(totalLen + chunkSize)) {
                break;
            }
            buffer = blockPtr->getBuffer();
        }
        const size_t readLen = ioReader->readIO(buffer + totalLen, chunkSize);
        if (readLen == 0) {
            break; // EOF
        }
        totalLen += readLen;
    }

    if (totalLen == 0) {
        return 0; // EOF
    }

    // ---- 块尾是否落在完整的行边界上 ----
    const bool stoppedAtTarget = (dataLineCount >= readTarget);
    const bool endsWithNewline = (buffer[totalLen - 1] == '\n');
    if (!stoppedAtTarget && !endsWithNewline) {
        /* 末尾记录被截断（无结尾换行）：整块按二进制压缩，保证数据不丢 */
        LOG_INFO("SAM tail block not aligned to record boundary, fallback to binary codec (blockId=%ld).",
                 (long)blockPtr->getBlockId());
        blockPtr->setDataLen((int64_t)totalLen);
        blockPtr->setBlockType(BINARY);
        return (int64_t)totalLen;
    }

    // ---- 确定本块保留到第几行 ----
    size_t keepLines = 0;
    if (startsWithHeader && splitHeader) {
        /* 头部块：只保留所有 @ 头部行，数据行留给后续块 */
        size_t lineStartPos = 0;
        for (size_t idx = 0; idx < npos.size(); ++idx) {
            if (buffer[lineStartPos] != '@') {
                break;
            }
            keepLines = idx + 1;
            lineStartPos = npos[idx] + 1;
        }
    } else {
        /* 数据块：保留头部行（若有）+ readsPerBlock 条数据行（不足则全保留） */
        size_t lineStartPos = 0;
        size_t dataLines = 0;
        for (size_t idx = 0; idx < npos.size(); ++idx) {
            keepLines = idx + 1;
            if (buffer[lineStartPos] != '@') {
                ++dataLines;
                if (dataLines >= readsPerBlock) {
                    break;
                }
            }
            lineStartPos = npos[idx] + 1;
        }
    }

    if (keepLines == 0 || keepLines > npos.size()) {
        /* 防御：没有完整行，整块按二进制 */
        LOG_INFO("SAM block has no complete line, fallback to binary codec (blockId=%ld).",
                 (long)blockPtr->getBlockId());
        blockPtr->setDataLen((int64_t)totalLen);
        blockPtr->setBlockType(BINARY);
        return (int64_t)totalLen;
    }

    // ---- 块尾对齐：保留行之后的内容移入缓存，供下一块使用 ----
    const size_t keepLen = npos[keepLines - 1] + 1;
    if (keepLen < totalLen) {
        const size_t tailLen = totalLen - keepLen;
        if (cacheLen > 0) {
            memmove(cache + tailLen, cache, cacheLen);
        }
        memcpy(cache, buffer + keepLen, tailLen);
        cacheLen += tailLen;
        totalLen = keepLen;
    }
    if (npos.size() > keepLines) {
        npos.resize(keepLines);
    }
    blockPtr->setDataLen((int64_t)totalLen);

    /*
     * 预分析在 coder 侧 actuator 中执行并持有其状态（压缩需要），这里只做换行符位置
     * 记录与结构校验，不重复运行 SamCodecActuator::preAnalysis。
     * GZ 输入时块类型仍按非 GZ 的 SAM 类型设置，由 io 层透明解压负责还原。
     */
    blockPtr->setBlockType(SAM);
    return (int64_t)totalLen;
}

BlockType SamBlockReader::analyzeBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    uint8_t* buffer = blockPtr->getBuffer();
    const int64_t dataLen = blockPtr->getDataLen();
    const char* bufPtr = reinterpret_cast<const char*>(buffer);

    // 记录换行符位置（memchr 自带 SIMD 优化）
    const char* searchStart = bufPtr;
    while (true) {
        const char* newlinePtr = static_cast<const char*>(memchr(searchStart, '\n',
                                                      dataLen - (searchStart - bufPtr)));
        if (newlinePtr == nullptr) {
            break;
        }
        blockPtr->getNpos().push_back(static_cast<uint32_t>(newlinePtr - bufPtr));
        searchStart = newlinePtr + 1;
    }

    if (blockPtr->getNpos().empty()) {
        return BINARY;
    }

    // 块尾对齐到完整行
    if (!alignToRecordBoundary(blockPtr, false)) {
        return BINARY;
    }

    /*
     * 结构校验兜底路径（SamBlockReader::readBlock 已覆写，通常不会走到这里）：
     * 不做全量 actuator 预分析，SAM 深校验与 BINARY 回落由 coder 侧 actuatorPreProc 负责。
     */
    return SAM;
}

namespace {

    /* ---- BAM 二进制字段解析（小端） ---- */
    inline int32_t bamI32(const uint8_t*& p) {
        int32_t v;
        memcpy(&v, p, 4);
        p += 4;
        return v;
    }
    inline uint16_t bamU16(const uint8_t*& p) {
        uint16_t v;
        memcpy(&v, p, 2);
        p += 2;
        return v;
    }

    const char* const BAM_CIGAR_OPS = "MIDNSHP=XB";
    /* BAM 4-bit 碱基编码：索引 0..15 -> "=ACMGRSVTWYHKDBN" */
    const char BAM_BASE_MAP[16] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                   'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

    void buildBamCigar(const uint8_t* data, uint16_t nCigarOp, std::string& out) {
        for (uint16_t i = 0; i < nCigarOp; ++i) {
            uint32_t op;
            memcpy(&op, data + i * 4, 4);
            const uint32_t len = op >> 4;      /* 高 28 位是长度 */
            const unsigned code = op & 0xF;    /* 低 4 位是操作码 */
            out += std::to_string(len);
            out += (code < 11) ? BAM_CIGAR_OPS[code] : 'M';
        }
    }

    void buildBamSeq(const uint8_t* data, int32_t lSeq, std::string& out) {
        for (int32_t i = 0; i < lSeq; ++i) {
            const uint8_t nibble = (data[i >> 1] >> (4 - 4 * (i & 1))) & 0xF;
            out += BAM_BASE_MAP[nibble];
        }
    }

    void buildBamQual(const uint8_t* data, int32_t lSeq, std::string& out) {
        bool missing = true;
        for (int32_t i = 0; i < lSeq; ++i) {
            if (data[i] != 0xFF) {
                missing = false;
                break;
            }
        }
        if (missing) {
            return;  /* 调用方统一输出 '*' */
        }
        /*
         * BAM 里 QUAL 存的是 phred 分数（htslib/samtools 写 BAM 时对 SAM 的 phred+33
         * 减 33 再存），转回 SAM 需加回 33。
         */
        for (int32_t i = 0; i < lSeq; ++i) {
            out += (char)(data[i] + 33);
        }
    }

    void appendBamAux(const uint8_t*& p, const uint8_t* end, std::string& line) {
        while (p + 3 <= end) {
            char tag[3] = {(char)p[0], (char)p[1], 0};
            const char type = (char)p[2];
            p += 3;
            line += "\t";
            line += tag;
            line += ":";
            /*
             * SAM 的辅助类型只有 A/i/f/Z/H/B。BAM 里整数子类型 c/C/s/S/i/I 一律按
             * 'i' 输出（与 samtools view 一致），数值原样保留。
             */
            line += 'i';
            line += ":";
            switch (type) {
            case 'A':
                /* 字符型，SAM 用 'A' */
                line[line.size() - 2] = 'A';
                line += (char)*p;
                p += 1;
                break;
            case 'c':
                line += std::to_string((int)(int8_t)*p);
                p += 1;
                break;
            case 'C':
                line += std::to_string((int)*p);
                p += 1;
                break;
            case 's': {
                int16_t v;
                memcpy(&v, p, 2);
                line += std::to_string(v);
                p += 2;
                break;
            }
            case 'S': {
                uint16_t v;
                memcpy(&v, p, 2);
                line += std::to_string(v);
                p += 2;
                break;
            }
            case 'i': {
                int32_t v;
                memcpy(&v, p, 4);
                line += std::to_string(v);
                p += 4;
                break;
            }
            case 'I': {
                uint32_t v;
                memcpy(&v, p, 4);
                line += std::to_string(v);
                p += 4;
                break;
            }
            case 'f': {
                float v;
                memcpy(&v, p, 4);
                char buf[64];
                snprintf(buf, sizeof(buf), "%g", (double)v);
                line[line.size() - 2] = 'f';
                line += buf;
                p += 4;
                break;
            }
            case 'Z':
            case 'H': {
                const size_t n = strnlen((const char*)p, end - p);
                line[line.size() - 2] = type;
                line.append((const char*)p, n);
                p += n + 1;
                break;
            }
            case 'B': {
                const char sub = (char)*p++;
                int32_t count;
                memcpy(&count, p, 4);
                p += 4;
                line[line.size() - 2] = 'B';
                line += sub;
                for (int32_t i = 0; i < count && p + 1 <= end; ++i) {
                    line += ',';
                    switch (sub) {
                    case 'c': line += std::to_string((int)(int8_t)*p); p += 1; break;
                    case 'C': line += std::to_string((int)*p); p += 1; break;
                    case 's': { int16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
                    case 'S': { uint16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
                    case 'i': { int32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
                    case 'I': { uint32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
                    case 'f': { float v; memcpy(&v, p, 4); char bf[64]; snprintf(bf, sizeof(bf), "%g", (double)v); line += bf; p += 4; } break;
                    default: return; break;
                    }
                }
                break;
            }
            default:
                return;   /* 未知类型，跳过剩余 */
            }
        }
    }

}  // namespace

size_t BamBlockReader::readRawFromSource(void* dst, size_t n) {
    size_t got = 0;
    if (detectLen > 0) {
        const size_t toCopy = (n < detectLen) ? n : detectLen;
        memcpy(dst, detectBuf, toCopy);
        detectLen -= toCopy;
        if (detectLen > 0) {
            memmove(detectBuf, detectBuf + toCopy, detectLen);
        }
        got = toCopy;
    }
    if (got < n) {
        got += ioReader->readIO((uint8_t*)dst + got, n - got);
    }
    return got;
}

size_t BamBlockReader::readBamBytes(void* dst, size_t n) {
    return readRawFromSource(dst, n);
}

size_t BamGzBlockReader::readBamBytes(void* dst, size_t n) {
    if (!inflateReady) {
        memset(&inflateState, 0, sizeof(inflateState));
        inflateInit2(&inflateState, 32 + MAX_WBITS);
        inflateReady = true;
        gzInLen = 0;
        gzInEof = false;
    }

    size_t got = 0;
    while (got < n) {
        inflateState.next_out = (Bytef*)((uint8_t*)dst + got);
        inflateState.avail_out = (uInt)(n - got);
        const int rc = inflate(&inflateState, Z_NO_FLUSH);
        got = n - inflateState.avail_out;
        if (rc == Z_STREAM_END) {
            /* BGZF 是多个 gzip 成员拼接，重置后继续 inflate 剩余输入 */
            if (inflateState.avail_in > 0) {
                inflateReset(&inflateState);
                continue;
            }
            break;
        }
        if (rc != Z_OK && rc != Z_BUF_ERROR) {
            break;   /* 数据损坏 */
        }
        if (inflateState.avail_in == 0) {
            if (gzInEof) {
                break;
            }
            const size_t inRead = readRawFromSource(gzInBuf, sizeof(gzInBuf));
            if (inRead == 0) {
                gzInEof = true;
                break;
            }
            inflateState.next_in = gzInBuf;
            inflateState.avail_in = (uInt)inRead;
        }
    }
    return got;
}

int32_t BamBlockReader::parseBamHeader() {
    uint8_t magic[4];
    if (readBamBytes(magic, 4) != 4 || memcmp(magic, "BAM\x1", 4) != 0) {
        LOG_ERROR("BAM magic not found.");
        return -1;
    }

    int32_t lText = 0;
    if (readBamBytes(&lText, 4) != 4 || lText < 0 || lText > (int32_t)(1 << 24)) {
        LOG_ERROR("Invalid BAM header text length %d", lText);
        return -1;
    }
    std::string text;
    text.resize(lText);
    if (lText > 0 && readBamBytes(&text[0], lText) != (size_t)lText) {
        LOG_ERROR("Failed to read BAM header text.");
        return -1;
    }

    int32_t nRef = 0;
    if (readBamBytes(&nRef, 4) != 4 || nRef < 0 || nRef > 100000) {
        LOG_ERROR("Invalid BAM ref count %d", nRef);
        return -1;
    }
    refs.clear();
    for (int32_t i = 0; i < nRef; ++i) {
        int32_t lName = 0;
        if (readBamBytes(&lName, 4) != 4 || lName <= 0 || lName > 4096) {
            LOG_ERROR("Invalid BAM ref name length %d", lName);
            return -1;
        }
        std::string name;
        name.resize(lName);
        if (readBamBytes(&name[0], lName) != (size_t)lName) {
            return -1;
        }
        while (!name.empty() && name.back() == '\0') {
            name.pop_back();
        }
        int32_t refLen = 0;
        if (readBamBytes(&refLen, 4) != 4) {
            return -1;
        }
        refs.push_back({name, refLen});
    }

    /* 以 BAM 头文本为基础，缺 @SQ 的参考序列补上 @SQ 行 */
    std::set<std::string> sqRefs;
    headerText.clear();
    {
        size_t start = 0;
        while (start < text.size()) {
            size_t end = text.find('\n', start);
            if (end == std::string::npos) {
                end = text.size();
            }
            std::string line = text.substr(start, end - start);
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) {
                line.pop_back();
            }
            if (!line.empty()) {
                if (line.rfind("@SQ", 0) == 0) {
                    size_t sn = line.find("SN:");
                    if (sn != std::string::npos) {
                        size_t p = sn + 3;
                        size_t e = line.find('\t', p);
                        if (e == std::string::npos) {
                            e = line.size();
                        }
                        sqRefs.insert(line.substr(p, e - p));
                    }
                }
                headerText += line;
                headerText += '\n';
            }
            if (end == text.size()) {
                break;
            }
            start = end + 1;
        }
    }
    for (const BamRef& ref : refs) {
        if (sqRefs.find(ref.name) == sqRefs.end()) {
            headerText += "@SQ\tSN:";
            headerText += ref.name;
            headerText += "\tLN:";
            headerText += std::to_string(ref.len);
            headerText += '\n';
        }
    }
    return 0;
}

int32_t BamBlockReader::parseBamRecord(const uint8_t* data, int32_t size, std::string& line) {
    const uint8_t* p = data;
    const uint8_t* end = data + size;

    const int32_t refID = bamI32(p);
    const int32_t pos = bamI32(p);
    const uint8_t lReadName = *p++;
    const uint8_t mapq = *p++;
    (void)bamU16(p);              // bin
    const uint16_t nCigarOp = bamU16(p);
    const uint16_t flag = bamU16(p);
    const int32_t lSeq = bamI32(p);
    const int32_t nextRefID = bamI32(p);
    const int32_t nextPos = bamI32(p);
    const int32_t tlen = bamI32(p);

    if (lReadName == 0 || p > end) {
        return -1;
    }
    std::string qname((const char*)p, lReadName - 1);
    p += lReadName;

    const uint8_t* cigarData = p;
    p += (size_t)nCigarOp * 4;
    const uint8_t* seqData = p;
    p += (size_t)(lSeq + 1) / 2;
    const uint8_t* qualData = p;
    p += (size_t)lSeq;

    if (p > end) {
        return -1;
    }

    const std::string rname = (refID >= 0 && refID < (int32_t)refs.size()) ? refs[refID].name : "*";
    const std::string rnext = (nextRefID >= 0 && nextRefID < (int32_t)refs.size())
        ? ((nextRefID == refID) ? "=" : refs[nextRefID].name) : "*";

    std::string cigar;
    buildBamCigar(cigarData, nCigarOp, cigar);
    if (cigar.empty()) {
        cigar = "*";
    }

    std::string seq;
    buildBamSeq(seqData, lSeq, seq);
    if (seq.empty()) {
        seq = "*";
    }

    std::string qual;
    buildBamQual(qualData, lSeq, qual);
    if (qual.empty()) {
        qual = "*";
    }

    line = qname;
    line += "\t"; line += std::to_string(flag);
    line += "\t"; line += rname;
    line += "\t"; line += std::to_string((pos >= 0) ? (pos + 1) : 0);
    line += "\t"; line += std::to_string(mapq);
    line += "\t"; line += cigar;
    line += "\t"; line += rnext;
    line += "\t"; line += std::to_string((nextPos >= 0) ? (nextPos + 1) : 0);
    line += "\t"; line += std::to_string(tlen);
    line += "\t"; line += seq;
    line += "\t"; line += qual;

    appendBamAux(p, end, line);
    return 0;
}

/*
 * BAM 读取：头部作为独立块返回（转成 SAM 头），数据区按 readsPerBlock 条 read
 * 解压成 SAM 行组成一个 SAM 块，交由 SAM 压缩器压缩。块类型统一标记为 BAM，
 * 打印 FileType 显示 BAM，内容与处理路径与 SAM 一致。
 */
int64_t BamBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    if (ioReader == nullptr || blockPtr == nullptr) {
        LOG_ERROR("IO reader or block pointer is null.");
        return -1;
    }
    blockPtr->setBlockId(blockId++);
    std::vector<size_t>& npos = blockPtr->getNpos();
    npos.clear();

    uint8_t* out = blockPtr->getBuffer();
    size_t outLen = 0;

    if (!headerParsed) {
        if (0 != parseBamHeader()) {
            LOG_ERROR("Parse BAM header failed.");
            blockPtr->setBlockType(BINARY);
            return -1;
        }
        headerParsed = true;
    }

    // 写入头部（若尚未写）
    if (!headerWritten) {
        if (headerText.empty()) {
            headerText = "@HD\tVN:1.6\n";
        }
        if (headerText.size() > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(headerText.size())) {
                return -1;
            }
            out = blockPtr->getBuffer();
        }
        memcpy(out, headerText.data(), headerText.size());
        outLen = headerText.size();
        for (size_t i = 0; i < headerText.size(); ++i) {
            if (headerText[i] == '\n') {
                npos.push_back(i);
            }
        }
        headerWritten = true;
        lastBlockHasData = !splitHeader;
        if (splitHeader) {
            blockPtr->setDataLen((int64_t)outLen);
            blockPtr->setBlockType(BAM);
            return (int64_t)outLen;
        }
    }

    // 数据区：逐条 read 解压成 SAM 行，凑满 readsPerBlock 条
    uint32_t reads = 0;
    while (reads < readsPerBlock) {
        int32_t blockSize = 0;
        if (readBamBytes(&blockSize, 4) != 4) {
            break;   // EOF
        }
        if (blockSize <= 0 || blockSize > (1 << 30)) {
            break;
        }
        std::vector<uint8_t> rec((size_t)blockSize);
        if (readBamBytes(rec.data(), rec.size()) != rec.size()) {
            break;
        }
        std::string line;
        if (0 != parseBamRecord(rec.data(), blockSize, line)) {
            break;
        }
        if (outLen + line.size() + 1 > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(outLen + line.size() + 1)) {
                break;
            }
            out = blockPtr->getBuffer();
        }
        memcpy(out + outLen, line.data(), line.size());
        outLen += line.size();
        out[outLen] = '\n';
        npos.push_back(outLen);
        outLen += 1;
        ++reads;
        lastBlockHasData = true;
    }

    if (outLen == 0) {
        return 0;   // EOF（头部已在之前块输出）
    }

    blockPtr->setDataLen((int64_t)outLen);
    blockPtr->setBlockType(BAM);
    return (int64_t)outLen;
}

namespace BlockUtil {
    bool isFastqBlock(BlockType type) {
        return (type == FASTQ_GEN2 || type == FASTQ_GEN3 || type == FASTQ_GEN2_GZIP || type == FASTQ_GEN3_GZIP);
    }

    bool isSAMBlock(BlockType type) {
        /* BAM 块内容是 SAM 文本，当作 SAM 处理，因此 SAM 支持的能力 BAM 一并支持 */
        return (type == SAM || type == SAM_GZIP || type == BAM);
    }

    bool isBAMBlock(BlockType type) {
        return (type == BAM);
    }

    bool isAuxiliaryBlock(BlockType type) {
        return (type == REFERENCE || type == REFERENCE_INDEX || type == QUAL_PRIOR);
    }

    std::string getBlockTypeName(BlockType type) {
        switch(type) {
            case TYPE_UNKNOW:
                return "TYPE_UNKNOW";
            case GZIP:
                return "GZIP";
            case BINARY:
                return "BINARY";
            case FASTQ_GEN2:
                return "FASTQ_GEN2";
            case FASTQ_GEN3:
                return "FASTQ_GEN3";
            case BINARY_GZIP:
                return "BINARY_GZIP";
            case FASTQ_GEN2_GZIP:
                return "FASTQ_GEN2";
            case FASTQ_GEN3_GZIP:
                return "FASTQ_GEN3";
            case BAM:
                return "BAM";
            case SAM:
                return "SAM";
            case SAM_GZIP:
                return "SAM";
            case PBGZFILE:
                return "PBGZFILE";
            default:
                return "UNKNOWN_TYPE";
        }
    }
}

namespace {

    bool startsWithPbgz(const uint8_t* buf, size_t len) {
        return len >= 4 && memcmp(buf, "PBGZ", 4) == 0;
    }

    bool startsWithBamMagic(const uint8_t* buf, size_t len) {
        return len >= 4 && buf[0] == 'B' && buf[1] == 'A' && buf[2] == 'M' && buf[3] == 1;
    }

    /* BGZF 魔数：1f 8b 08 04（BAM 本身就是 BGZF 流） */
    bool isBgzfStream(const uint8_t* buf, size_t len) {
        return len >= 4 && buf[0] == 0x1f && buf[1] == 0x8b && buf[2] == 0x08 && (buf[3] & 0x04) == 0x04;
    }

    /* 解压首个 gzip/BGZF 块；成功时 out 为解压后的字节，返回 true */
    bool inflateFirstBlock(const uint8_t* buf, size_t len, std::vector<uint8_t>& out) {
        if (len < 20) {
            return false;
        }
        z_stream zs;
        memset(&zs, 0, sizeof(zs));
        if (inflateInit2(&zs, 32 + MAX_WBITS) != Z_OK) {
            return false;
        }
        out.resize(1 << 20);
        zs.next_in = const_cast<Bytef*>(buf);
        zs.avail_in = (uInt)len;
        zs.next_out = out.data();
        zs.avail_out = (uInt)out.size();
        int rc = inflate(&zs, Z_SYNC_FLUSH);
        inflateEnd(&zs);
        out.resize(out.size() - zs.avail_out);
        return (rc == Z_OK || rc == Z_STREAM_END) && !out.empty();
    }

    /* 检查解压后是否以 BAM 魔数开头 */
    bool isBamByInflate(const uint8_t* buf, size_t len) {
        std::vector<uint8_t> out;
        return inflateFirstBlock(buf, len, out) && out.size() >= 4 &&
               out[0] == 'B' && out[1] == 'A' && out[2] == 'M' && out[3] == 1;
    }

    bool isBamSample(const uint8_t* buf, size_t len) {
        if (startsWithBamMagic(buf, len)) {
            return true;
        }
        if (isBgzfStream(buf, len)) {
            return isBamByInflate(buf, len);
        }
        return false;
    }

    /* 校验整个探测长度内的每一行，防止"首个记录像 FASTQ/SAM、其后全不是"的文件误入对应 reader。
     * 若样本内出现任何一行不合规，即整体判为 BINARY，交给 analyzeBlock 兜底的场景不复存在。 */
    bool isSamSample(const uint8_t* buf, size_t len) {
        if (len < 4 || buf[0] != '@') {
            return false;
        }
        static const char* const samHeads[] = {"@HD", "@SQ", "@RG", "@PG", "@CO"};

        size_t lineCount = 0;
        size_t lineStart = 0;
        for (size_t i = 0; i < len; ++i) {
            if (buf[i] != '\n') {
                continue;
            }
            const size_t lineLen = i - lineStart;
            if (lineLen == 0) {
                return false;   // 空行不是合法的 SAM 行
            }
            if (buf[lineStart] == '@') {
                // 头部行：必须是已知的 SAM 头部类型
                bool known = false;
                for (const char* head : samHeads) {
                    if (lineLen >= 3 && memcmp(buf + lineStart, head, 3) == 0) {
                        known = true;
                        break;
                    }
                }
                if (!known) {
                    return false;
                }
            } else {
                // 数据行：SAM 至少 11 个必选字段（QNAME..QUAL），即至少 10 个制表符
                int tabs = 0;
                for (size_t k = lineStart; k < i; ++k) {
                    if (buf[k] == '\t') {
                        ++tabs;
                    }
                }
                if (tabs < 10) {
                    return false;
                }
            }
            ++lineCount;
            lineStart = i + 1;
        }
        return lineCount >= 1;
    }

    /* 校验整个探测长度内的每一行 FASTQ 结构（@ID/碱基/+/质量行），任一行不合规即判 BINARY。
     * 要求首条完整记录存在（line >= 4），防止首个记录像 FASTQ、其后全不是的文件误入 Fastq reader。 */
    bool isFastqSample(const uint8_t* buf, size_t len) {
        if (len < 2 || buf[0] != '@') {
            return false;
        }
        static bool isValidBase[256] = { false };
        static bool tableInit = false;
        if (!tableInit) {
            isValidBase[(uint8_t)'A'] = true;
            isValidBase[(uint8_t)'C'] = true;
            isValidBase[(uint8_t)'G'] = true;
            isValidBase[(uint8_t)'N'] = true;
            isValidBase[(uint8_t)'T'] = true;
            isValidBase[(uint8_t)'a'] = true;
            isValidBase[(uint8_t)'c'] = true;
            isValidBase[(uint8_t)'g'] = true;
            isValidBase[(uint8_t)'n'] = true;
            isValidBase[(uint8_t)'t'] = true;
            tableInit = true;
        }

        size_t line = 0;
        size_t start = 0;
        size_t baseLen = 0;
        for (size_t i = 1; i < len; ++i) {
            if (buf[i] != '\n') {
                continue;
            }
            const size_t lineLen = i - start;
            switch (line & 0x3) {
            case 0: {
                if (lineLen < 1 || buf[start] != '@') {
                    return false;
                }
                break;
            }
            case 1: {
                if (lineLen == 0) {
                    return false;
                }
                for (size_t k = start; k < i; ++k) {
                    if (!isValidBase[(uint8_t)buf[k]]) {
                        return false;
                    }
                }
                baseLen = lineLen;
                break;
            }
            case 2: {
                if (lineLen < 1 || buf[start] != '+') {
                    return false;
                }
                break;
            }
            case 3: {
                if (lineLen != baseLen) {
                    return false;
                }
                break;
            }
            default:
                return false;
            }
            ++line;
            start = i + 1;
        }
        return line >= 4;
    }

}  // namespace

/*
 * 判定文件是否为 BAM：原始 BAM 魔数，或 BGZF/gzip 流解压后以 BAM 魔数开头。
 * 供引擎在创建 ioReader 时区分"透明 gz 解压"与"交给 FileReader + BamGzBlockReader"。
 */
bool BlockUtil::isBamFile(const std::string& fileName) {
    FILE* fp = fopen(fileName.c_str(), "rb");
    if (fp == nullptr) {
        return false;
    }
    std::vector<uint8_t> buf(BLOCK_TYPE_DETECT_SIZE);
    const size_t len = fread(buf.data(), 1, buf.size(), fp);
    fclose(fp);
    return isBamSample(buf.data(), len);
}

/* 按文件内容探测输入格式：gz/BGZF 自动解压后再判。仅用于普通文件，管道/STDIN 勿用。 */
BlockType BlockUtil::detectInputFileType(const std::string& fileName) {
    FILE* fp = fopen(fileName.c_str(), "rb");
    if (fp == nullptr) {
        return TYPE_UNKNOW;
    }
    std::vector<uint8_t> buf(BLOCK_TYPE_DETECT_SIZE);
    const size_t len = fread(buf.data(), 1, buf.size(), fp);
    fclose(fp);
    if (len == 0) {
        return TYPE_UNKNOW;
    }
    if (startsWithPbgz(buf.data(), len)) {
        return PBGZFILE;
    }
    if (startsWithBamMagic(buf.data(), len)) {
        return BAM;
    }
    /* 通用 gzip 魔数（不要求 BGZF 的 FEXTRA 标志，普通 .gz 也要能解） */
    if (len >= 2 && buf[0] == 0x1f && buf[1] == 0x8b) {
        std::vector<uint8_t> plain;
        if (!inflateFirstBlock(buf.data(), len, plain)) {
            return BINARY;
        }
        /* 外层 gzip 解开后可能仍是 BGZF（如 gzip 过的 BAM），再解一层 */
        if (plain.size() >= 2 && plain[0] == 0x1f && plain[1] == 0x8b) {
            std::vector<uint8_t> plain2;
            if (inflateFirstBlock(plain.data(), plain.size(), plain2)) {
                if (plain2.size() >= 4 && plain2[0] == 'B' && plain2[1] == 'A' &&
                    plain2[2] == 'M' && plain2[3] == 1) {
                    return BAM;
                }
                if (isSamSample(plain2.data(), plain2.size())) {
                    return SAM;
                }
                if (isFastqSample(plain2.data(), plain2.size())) {
                    return FASTQ_GEN2;
                }
                return BINARY;
            }
        }
        if (plain.size() >= 4 && plain[0] == 'B' && plain[1] == 'A' &&
            plain[2] == 'M' && plain[3] == 1) {
            return BAM;
        }
        if (isSamSample(plain.data(), plain.size())) {
            return SAM;
        }
        if (isFastqSample(plain.data(), plain.size())) {
            return FASTQ_GEN2;
        }
        return BINARY;
    }
    if (isSamSample(buf.data(), len)) {
        return SAM;
    }
    if (isFastqSample(buf.data(), len)) {
        return FASTQ_GEN2;
    }
    return BINARY;
}

/*
 * 读取一小部分数据解析文件格式，按格式创建对应的 BlockReader 子类。
 * 预读数据交给子类并入首个块；PBGZ 需要回到文件起点由 PbgzBlockReader 自行解析。
 */
BlockReader* BlockFactory::createBlockReader(IOReader* ioReader, uint8_t compressLevel, bool splitSamHeader) {
    if (ioReader == nullptr) {
        LOG_ERROR("Create block reader failed: io reader is null.");
        return nullptr;
    }

    /* SAM 数据块按 read 行数分块：1-5 -> 10000，6-7 -> 25000，8-9 -> 100000 */
    uint32_t samReadsPerBlock = 10000;
    if (compressLevel >= 8) {
        samReadsPerBlock = 100000;
    } else if (compressLevel >= 6) {
        samReadsPerBlock = 25000;
    }

    uint8_t* detectBuf = MemoryUtil::safeAlloc<uint8_t>(BLOCK_TYPE_DETECT_SIZE);
    if (detectBuf == nullptr) {
        LOG_ERROR("Create block reader failed: alloc detect buffer failed.");
        return nullptr;
    }
    const size_t detectLen = ioReader->readIO(detectBuf, BLOCK_TYPE_DETECT_SIZE);
    if (detectLen == (size_t)-1) {
        LOG_ERROR("Create block reader failed: read detect data failed.");
        MemoryUtil::safeFree(detectBuf);
        return nullptr;
    }

    BlockReader* reader = nullptr;
    if (startsWithPbgz(detectBuf, detectLen)) {
        /*
         * PBGZ：预读字节交还给 PbgzBlockReader（管道输入不可 seek 回退），
         * 头部解析由 PbgzBlockReader::init() 从预读缓冲消费。
         */
        reader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader, detectBuf, detectLen);
    } else if (isBamSample(detectBuf, detectLen)) {
        if (startsWithBamMagic(detectBuf, detectLen)) {
            /* 已是原始 BAM 字节流（io 层已解 BGZF，如 .bam） */
            reader = MemoryUtil::safeNewClass<BamBlockReader>(ioReader, detectBuf, detectLen, samReadsPerBlock, splitSamHeader);
        } else {
            /* 内层仍是 BGZF（如 .bam.gz），需先 inflate 成原始 BAM */
            reader = MemoryUtil::safeNewClass<BamGzBlockReader>(ioReader, detectBuf, detectLen, samReadsPerBlock, splitSamHeader);
        }
    } else if (isSamSample(detectBuf, detectLen)) {
        reader = MemoryUtil::safeNewClass<SamBlockReader>(ioReader, detectBuf, detectLen, samReadsPerBlock, splitSamHeader);
    } else if (isFastqSample(detectBuf, detectLen)) {
        reader = MemoryUtil::safeNewClass<FastqBlockReader>(ioReader, detectBuf, detectLen);
    } else {
        reader = MemoryUtil::safeNewClass<BinaryBlockReader>(ioReader, detectBuf, detectLen);
    }

    MemoryUtil::safeFree(detectBuf);

    if (reader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return nullptr;
    }
    /* 读取目标按 -l 决定（输入块只分配 1MB，读时按此目标按需扩容） */
    reader->setReadBlockBytes(ConfigManager::getInstance().getBlockSizeByCompressLevel(compressLevel));
    if (0 != reader->init()) {
        LOG_ERROR("Create block reader init failed.");
        MemoryUtil::safeDeleteClass(reader);
        return nullptr;
    }
    return reader;
}

BlockWriter* BlockFactory::createBlockWriter(IOWriter* ioWriter) {
    if (ioWriter == nullptr) {
        LOG_ERROR("Create block writer failed: io writer is null.");
        return nullptr;
    }
    BlockWriter* blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
    if (blockWriter == nullptr || 0 != blockWriter->init()) {
        LOG_ERROR("Create block writer failed.");
        MemoryUtil::safeDeleteClass(blockWriter);
        return nullptr;
    }
    return blockWriter;
}

PbgzBlockWriter* BlockFactory::createPbgzBlockWriter(IOWriter* ioWriter) {
    if (ioWriter == nullptr) {
        LOG_ERROR("Create pbgz block writer failed: io writer is null.");
        return nullptr;
    }
    PbgzBlockWriter* blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
    if (blockWriter == nullptr || 0 != blockWriter->init()) {
        LOG_ERROR("Create pbgz block writer failed.");
        MemoryUtil::safeDeleteClass(blockWriter);
        return nullptr;
    }
    return blockWriter;
}

int64_t PbgzBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType __attribute__ ((unused)) fileType) {
    if (pbgzFileReader == nullptr || blockPtr == nullptr) {
        return -1;
    }

    PbgzDataBlock pbgzDataBlock;
    pbgzDataBlock.setDataPtr(blockPtr->getBuffer());
    if (0 != pbgzFileReader->readDataBlock(pbgzDataBlock, blockPtr)) {
        LOG_ERROR("Read Pbgz data block failed.");
        return -1;
    }

    if (0 != pbgzDataBlock.verifyCheckSum()) {
        LOG_ERROR("Verify pbgz block checksum failed, block id = %ld", pbgzDataBlock.getMetaData("blockid").asInt64());
        return -1;
    }
    
    blockPtr->setDataLen(pbgzDataBlock.getMetaData("datalen").asInt64());
    blockPtr->setMetaLen(pbgzDataBlock.getMetaData("metalen").asInt64());
    blockPtr->setBlockId(pbgzDataBlock.getMetaData("blockid").asInt64());
    /*
     * "本块属于哪个包"是块的固有属性，必须由唯一的生产者写入。
     * 放在各调用点去设，只要漏掉一处（区域查询、头部预读等），
     * 该路径上的块就会拿 0 当包起点，把辅助块的相对地址错译成别包的地址。
     */
    blockPtr->setPackageStart((int64_t)pbgzFileReader->getCurrentFileStart());
    blockPtr->setPackageIndex(pbgzFileReader->getCurrentFileIndex());
    std::string blockType = pbgzDataBlock.getMetaData("blocktype").asString();
    if (blockType == "fastq_gen2" || blockType == "fastq_gen2_gzip") {
        /* 旧文件里 *_gzip 块类型是"原始输入为 GZ"的标记；统一按非 GZ 类型处理 */
        blockPtr->setBlockType(FASTQ_GEN2);
    } else if (blockType == "fastq_gen3" || blockType == "fastq_gen3_gzip") {
        blockPtr->setBlockType(FASTQ_GEN3);
    } else if (blockType == "binary" || blockType == "binary_gzip") {
        blockPtr->setBlockType(BINARY);
    } else if (blockType == "refe_gene") {
        blockPtr->setBlockType(REFERENCE);
    } else if (blockType == "refe_gene_index") {
        blockPtr->setBlockType(REFERENCE_INDEX);
    } else if (blockType == "sam" || blockType == "sam_gzip") {
        blockPtr->setBlockType(SAM);
    } else if (blockType == "bam") {
        blockPtr->setBlockType(BAM);
    } else if (blockType == "qual_prior") {
        blockPtr->setBlockType(QUAL_PRIOR);
    }
    // Copy entire block information
    memcpy(blockPtr->getBuffer(), pbgzDataBlock.getDataPtr(), pbgzDataBlock.getDataLength());

    LOG_DEBUG("Read One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    
    return pbgzDataBlock.getDataLength();
}   

int32_t PbgzBlockReader::init() {
    if (ioReader == nullptr) {
        return -1;
    }
    if (pbgzFileReader == nullptr) {
        LOG_ERROR("Create PbgzFileReader failed");
        return -1;
    }
    return pbgzFileReader->open();
}

int32_t BlockWriter::writeBlock(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr || ioWriter == nullptr) {
        return -1;
    }
    if (blockPtr->getDataLen() != 0 || blockPtr->getMetaLen() != 0) {
        ioWriter->writeIO(blockPtr->getBuffer(), blockPtr->getDataLen());
    }
    LOG_DEBUG("Write One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    return 0;
}

int32_t PbgzBlockWriter::init() {
   if (ioWriter == nullptr || pbgzFileWriter == nullptr) {
        return -1;
    }
    pbgzFileWriter->open();
    return 0;
}


int32_t PbgzBlockWriter::writeBaseFileMeta() {
    if (ioWriter == nullptr || pbgzFileWriter == nullptr) {
        return -1;
    }
    pbgzFileWriter->getBaseFileMeta().setMetaData("writer", "pbgz_writer_v" + PbgzManager::getInstance().getVersion());
    pbgzFileWriter->getBaseFileMeta().setMetaData("hashmethod", "md5");
    pbgzFileWriter->writeBaseFileMeta();
    return 0;
}

int32_t PbgzBlockWriter::writeDynamicFileMeta() {
    if (ioWriter == nullptr || pbgzFileWriter == nullptr) {
        return -1;
    }
    pbgzFileWriter->writeDynamicFileMeta();
    return 0;
}

void PbgzBlockWriter::updateHeadExt(){ 
    if (ioWriter == nullptr || pbgzFileWriter == nullptr) {
        return;
    }

    FileWriter* pFileWrite = dynamic_cast<FileWriter*>(ioWriter);
    if (pFileWrite == nullptr) {
        return;
    }

    pbgzFileWriter->updateMetaOffset(pFileWrite->getCurrentPos());
}

int32_t PbgzBlockWriter::writeBlock(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr) {
        return -1;
    }

    if (pbgzFileWriter == nullptr) {
        return -1;
    }

    LOG_DEBUG("Write One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    
    if (blockPtr->getDataLen() == 0 && blockPtr->getMetaLen() == 0) {
        return 0;
    }

    PbgzDataBlock dataBlock;    
    dataBlock.setBlockData(blockPtr->getBuffer(), blockPtr->getTotalDataLen());
    dataBlock.setMetaData("datalen", (Json::Value::UInt64)blockPtr->getDataLen());
    dataBlock.setMetaData("metalen", (Json::Value::UInt64)blockPtr->getMetaLen());
    dataBlock.setMetaData("blockid", blockPtr->getBlockId());
    if (blockPtr->getBlockType() == FASTQ_GEN2) {
        dataBlock.setMetaData("blocktype", "fastq_gen2");
    } else if (blockPtr->getBlockType() == FASTQ_GEN3) {
        dataBlock.setMetaData("blocktype", "fastq_gen3");
    } else if (blockPtr->getBlockType() == FASTQ_GEN2_GZIP) {
        dataBlock.setMetaData("blocktype", "fastq_gen2_gzip");
    } else if (blockPtr->getBlockType() == FASTQ_GEN3_GZIP) {
        dataBlock.setMetaData("blocktype", "fastq_gen3_gzip");
    } else if (blockPtr->getBlockType() == BINARY) {
        dataBlock.setMetaData("blocktype", "binary");
    } else if (blockPtr->getBlockType() == REFERENCE) {
        dataBlock.setMetaData("blocktype", "refe_gene");
    } else if (blockPtr->getBlockType() == REFERENCE_INDEX) {
        dataBlock.setMetaData("blocktype", "refe_gene_index");
    } else if (blockPtr->getBlockType() == SAM || blockPtr->getBlockType() == SAM_GZIP) {
        dataBlock.setMetaData("blocktype", "sam");
    } else if (blockPtr->getBlockType() == BAM) {
        dataBlock.setMetaData("blocktype", "bam");
    } else if (blockPtr->getBlockType() == QUAL_PRIOR) {
        dataBlock.setMetaData("blocktype", "qual_prior");
    }

    /// Calculate checksum of Meta and data
    dataBlock.calcChecksum();
    return pbgzFileWriter->writeBlockData(dataBlock);
}

/* ==================== BamWriter（解压 -b：SAM -> BAM） ==================== */

namespace {

    /* ---- 小端字节流写出 ---- */
    inline void putU8(std::vector<uint8_t>& out, uint8_t v) {
        out.push_back(v);
    }
    inline void putU16(std::vector<uint8_t>& out, uint16_t v) {
        out.push_back((uint8_t)(v & 0xFF));
        out.push_back((uint8_t)((v >> 8) & 0xFF));
    }
    inline void putI32(std::vector<uint8_t>& out, int32_t v) {
        uint32_t u = (uint32_t)v;
        out.push_back((uint8_t)(u & 0xFF));
        out.push_back((uint8_t)((u >> 8) & 0xFF));
        out.push_back((uint8_t)((u >> 16) & 0xFF));
        out.push_back((uint8_t)((u >> 24) & 0xFF));
    }
    inline void putU32(std::vector<uint8_t>& out, uint32_t v) {
        putI32(out, (int32_t)v);
    }

    /* ---- BAM CIGAR 操作码（与 BAM_CIGAR_OPS "MIDNSHP=XB" 的索引一致） ---- */
    int cigarOpCode(char c) {
        switch (c) {
        case 'M': return 0;
        case 'I': return 1;
        case 'D': return 2;
        case 'N': return 3;
        case 'S': return 4;
        case 'H': return 5;
        case 'P': return 6;
        case '=': return 7;
        case 'X': return 8;
        default:  return -1;
        }
    }

    /* 参考序列上消耗长度的 CIGAR 操作：M/D/N/=/X */
    bool cigarConsumesRef(int code) {
        return code == 0 || code == 2 || code == 3 || code == 7 || code == 8;
    }

    /* 解析 CIGAR 字符串（"*" 或空时返回 true 且无操作）；refSpan 为参考序列消耗长度 */
    bool parseCigar(const char* s, size_t n, std::vector<uint32_t>& ops, int64_t& refSpan) {
        ops.clear();
        refSpan = 0;
        size_t i = 0;
        while (i < n) {
            size_t j = i;
            while (j < n && s[j] >= '0' && s[j] <= '9') {
                ++j;
            }
            if (j == i) {
                return false;   // 缺少长度
            }
            uint32_t len = 0;
            for (size_t k = i; k < j; ++k) {
                len = len * 10 + (uint32_t)(s[k] - '0');
            }
            if (j >= n) {
                return false;   // 缺少操作符
            }
            const int code = cigarOpCode(s[j]);
            if (code < 0) {
                return false;
            }
            ops.push_back((len << 4) | (uint32_t)code);
            if (cigarConsumesRef(code)) {
                refSpan += len;
            }
            i = j + 1;
        }
        return true;
    }

    /* BAM 4-bit 碱基编码：字符 -> 索引（与读取侧 BAM_BASE_MAP 互逆） */
    int baseNibble(char c) {
        switch (c) {
        case '=': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'M': return 3;
        case 'G': return 4;
        case 'R': return 5;
        case 'S': return 6;
        case 'V': return 7;
        case 'T': return 8;
        case 'W': return 9;
        case 'Y': return 10;
        case 'H': return 11;
        case 'K': return 12;
        case 'D': return 13;
        case 'B': return 14;
        case 'N': return 15;
        default:  return 15;
        }
    }

    /* 经典 reg2bin：beg 0-based 起点，end 0-based 终点（不含） */
    uint16_t samReg2Bin(int64_t beg, int64_t end) {
        const int64_t e = end - 1;
        if ((beg >> 14) == (e >> 14)) {
            return (uint16_t)(((1 << 15) - 1) / 7 + (beg >> 14));
        }
        if ((beg >> 17) == (e >> 17)) {
            return (uint16_t)(((1 << 12) - 1) / 7 + (beg >> 17));
        }
        if ((beg >> 20) == (e >> 20)) {
            return (uint16_t)(((1 << 9) - 1) / 7 + (beg >> 20));
        }
        if ((beg >> 23) == (e >> 23)) {
            return (uint16_t)(((1 << 6) - 1) / 7 + (beg >> 23));
        }
        if ((beg >> 26) == (e >> 26)) {
            return (uint16_t)(((1 << 3) - 1) / 7 + (beg >> 26));
        }
        return 0;
    }

    /* 从 @SQ 行里取 TAG 字段值（"SN:" / "LN:"） */
    std::string sqFieldValue(const std::string& line, const char* tag) {
        const size_t pos = line.find(tag);
        if (pos == std::string::npos) {
            return "";
        }
        size_t valStart = pos + strlen(tag);
        size_t end = line.find('\t', valStart);
        if (end == std::string::npos) {
            end = line.size();
        }
        return line.substr(valStart, end - valStart);
    }

    /*
     * 把一条 SAM 可选字段 "TAG:TYPE:VALUE" 转成 BAM aux 字节（含 TAG 与类型），
     * 失败返回 -1（调用方跳过该字段）。命名区别于读取侧的 appendBamAux（BAM->SAM）。
     */
    int appendSamAuxToBam(const char* opt, size_t n, std::vector<uint8_t>& out) {
        size_t p = 0;
        while (p < n && opt[p] != ':') {
            ++p;
        }
        if (p != 2) {
            return -1;   // TAG 必须是 2 字符
        }
        if (p + 3 > n || opt[p + 2] != ':') {
            return -1;   // 期待 "TAG:TYPE:VALUE"
        }
        const char type = opt[p + 1];
        const char* value = opt + p + 3;
        const size_t vlen = n - (p + 3);

        putU8(out, (uint8_t)opt[0]);
        putU8(out, (uint8_t)opt[1]);
        switch (type) {
        case 'A': {
            if (vlen != 1) {
                return -1;
            }
            putU8(out, (uint8_t)'A');
            putU8(out, (uint8_t)value[0]);
            break;
        }
        case 'i': {
            long long v = 0;
            try {
                v = std::stoll(std::string(value, vlen));
            } catch (...) {
                return -1;
            }
            if (v >= INT8_MIN && v <= INT8_MAX) {
                putU8(out, (uint8_t)'c');
                putU8(out, (uint8_t)(int8_t)v);
            } else if (v >= 0 && v <= UINT8_MAX) {
                putU8(out, (uint8_t)'C');
                putU8(out, (uint8_t)v);
            } else if (v >= INT16_MIN && v <= INT16_MAX) {
                putU8(out, (uint8_t)'s');
                putU16(out, (uint16_t)(int16_t)v);
            } else if (v >= 0 && v <= UINT16_MAX) {
                putU8(out, (uint8_t)'S');
                putU16(out, (uint16_t)v);
            } else if (v >= INT32_MIN && v <= INT32_MAX) {
                putU8(out, (uint8_t)'i');
                putI32(out, (int32_t)v);
            } else {
                putU8(out, (uint8_t)'I');
                putU32(out, (uint32_t)v);
            }
            break;
        }
        case 'f': {
            double d = strtod(std::string(value, vlen).c_str(), nullptr);
            float f = (float)d;
            putU8(out, (uint8_t)'f');
            uint8_t raw[4];
            memcpy(raw, &f, 4);
            out.insert(out.end(), raw, raw + 4);
            break;
        }
        case 'Z':
        case 'H': {
            putU8(out, (uint8_t)type);
            out.insert(out.end(), value, value + vlen);
            putU8(out, 0);
            break;
        }
        case 'B': {
            /* VALUE 形如 "SUBTYPE,V1,V2,..." */
            if (vlen < 3 || value[1] != ',') {
                return -1;
            }
            const char sub = value[0];
            std::vector<double> vals;
            {
                size_t i = 2;
                while (i <= vlen) {
                    size_t j = i;
                    while (j < vlen && value[j] != ',') {
                        ++j;
                    }
                    try {
                        vals.push_back(std::stod(std::string(value + i, j - i)));
                    } catch (...) {
                        return -1;
                    }
                    i = j + 1;
                }
            }
            putU8(out, (uint8_t)'B');
            putU8(out, (uint8_t)sub);
            const size_t countPos = out.size();
            putI32(out, 0);   // 计数占位，随后回填
            for (double dv : vals) {
                switch (sub) {
                case 'c': putU8(out, (uint8_t)(int8_t)dv); break;
                case 'C': putU8(out, (uint8_t)dv); break;
                case 's': putU16(out, (uint16_t)(int16_t)dv); break;
                case 'S': putU16(out, (uint16_t)dv); break;
                case 'i': putI32(out, (int32_t)dv); break;
                case 'I': putU32(out, (uint32_t)dv); break;
                case 'f': {
                    float f = (float)dv;
                    uint8_t raw[4];
                    memcpy(raw, &f, 4);
                    out.insert(out.end(), raw, raw + 4);
                    break;
                }
                default:
                    return -1;
                }
            }
            const int32_t cnt = (int32_t)vals.size();
            out[countPos] = (uint8_t)(cnt & 0xFF);
            out[countPos + 1] = (uint8_t)((cnt >> 8) & 0xFF);
            out[countPos + 2] = (uint8_t)((cnt >> 16) & 0xFF);
            out[countPos + 3] = (uint8_t)((cnt >> 24) & 0xFF);
            break;
        }
        default:
            return -1;
        }
        return 0;
    }

}  // namespace

BamWriter::BamWriter(IOWriter* pIoWriter) : BlockWriter(pIoWriter) {
    headerWritten = false;
    passThrough = false;
    finished = false;
    bgzfLen = 0;
    bgzfReady = false;
    memset(&bgzfZs, 0, sizeof(bgzfZs));
}

BamWriter::~BamWriter() {
    if (bgzfReady) {
        (void)finish();   // 兜底：残留块与 EOF 标记写出
        deflateEnd(&bgzfZs);
    }
}

int32_t BamWriter::writeRaw(const void* data, size_t len) {
    if (ioWriter == nullptr || len == 0) {
        return 0;
    }
    const size_t n = ioWriter->writeIO(data, len);
    if (n != len) {
        ioWriter->latchWriteError(-1);
        return -1;
    }
    return 0;
}

int32_t BamWriter::bgzfFlushBlock() {
    if (bgzfLen == 0) {
        return 0;
    }
    if (!bgzfReady) {
        memset(&bgzfZs, 0, sizeof(bgzfZs));
        /* 原始 deflate 流：gzip 头/尾由 BGZF 块自己拼 */
        if (deflateInit2(&bgzfZs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            LOG_ERROR("BamWriter: deflateInit2 failed.");
            return -1;
        }
        bgzfReady = true;
    }

    const size_t uncomprLen = bgzfLen;
    std::vector<uint8_t> comp(deflateBound(&bgzfZs, (uLong)uncomprLen) + 64);
    bgzfZs.next_in = bgzfBuf;
    bgzfZs.avail_in = (uInt)uncomprLen;
    bgzfZs.next_out = comp.data() + 18;   // 18 字节 gzip 头留白
    bgzfZs.avail_out = (uInt)(comp.size() - 18);
    const int rc = deflate(&bgzfZs, Z_FINISH);
    if (rc != Z_STREAM_END) {
        LOG_ERROR("BamWriter: deflate block failed, rc=%d.", rc);
        return -1;
    }
    const size_t compLen = (comp.size() - 18) - bgzfZs.avail_out;
    if (deflateReset(&bgzfZs) != Z_OK) {
        return -1;
    }

    const uint32_t crc = crc32(0L, bgzfBuf, (uInt)uncomprLen);
    const uint32_t total = (uint32_t)(18 + compLen + 8);
    const uint8_t header[18] = {
        0x1f, 0x8b, 0x08, 0x04,          /* gzip magic + deflate + FEXTRA */
        0, 0, 0, 0,                       /* mtime */
        0, 0xff,                          /* XFL / OS */
        0x06, 0x00,                       /* XLEN = 6 */
        0x42, 0x43,                       /* "BC" */
        0x02, 0x00,                       /* SLEN = 2 */
        (uint8_t)((total - 1) & 0xFF), (uint8_t)(((total - 1) >> 8) & 0xFF)   /* BSIZE */
    };
    uint8_t trailer[8];
    memcpy(trailer, &crc, 4);
    const uint32_t isize = (uint32_t)uncomprLen;
    memcpy(trailer + 4, &isize, 4);

    if (0 != writeRaw(header, sizeof(header)) ||
        0 != writeRaw(comp.data() + 18, compLen) ||
        0 != writeRaw(trailer, sizeof(trailer))) {
        return -1;
    }
    bgzfLen = 0;
    return 0;
}

int32_t BamWriter::bgzfWrite(const void* data, size_t len) {
    const uint8_t* p = (const uint8_t*)data;
    while (len > 0) {
        size_t n = sizeof(bgzfBuf) - bgzfLen;
        if (n > len) {
            n = len;
        }
        memcpy(bgzfBuf + bgzfLen, p, n);
        bgzfLen += n;
        p += n;
        len -= n;
        if (bgzfLen == sizeof(bgzfBuf)) {
            if (0 != bgzfFlushBlock()) {
                return -1;
            }
        }
    }
    return 0;
}

int32_t BamWriter::finish() {
    if (finished) {
        return 0;
    }
    finished = true;
    if (passThrough) {
        return 0;   // 透传模式没有 BGZF 状态
    }
    if (0 != bgzfFlushBlock()) {
        return -1;
    }
    /* BGZF EOF 标记（空块），标准 BAM 的收尾块 */
    static const uint8_t kEofBlock[28] = {
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff,
        0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
        0x03, 0x00, 0, 0, 0, 0, 0, 0, 0, 0
    };
    return writeRaw(kEofBlock, sizeof(kEofBlock));
}

int32_t BamWriter::writeBamHeader(const uint8_t* data, size_t len, size_t& dataStart) {
    refs.clear();
    refIndex.clear();
    std::string headerText;

    size_t pos = 0;
    while (pos < len && data[pos] == '@') {
        size_t nl = pos;
        while (nl < len && data[nl] != '\n') {
            ++nl;
        }
        size_t lineLen = nl - pos;
        if (lineLen > 0 && data[pos + lineLen - 1] == '\r') {
            --lineLen;
        }
        if (lineLen >= 3 && memcmp(data + pos, "@SQ", 3) == 0) {
            const std::string line((const char*)data + pos, lineLen);
            const std::string sn = sqFieldValue(line, "SN:");
            const std::string ln = sqFieldValue(line, "LN:");
            if (!sn.empty()) {
                int32_t refLen = 0;
                if (!ln.empty()) {
                    try {
                        refLen = std::stoi(ln);
                    } catch (...) {
                        refLen = 0;
                    }
                }
                refs.push_back({sn, refLen});
            }
        }
        headerText.append((const char*)data + pos, lineLen);
        headerText += '\n';
        pos = nl + 1;
    }
    if (headerText.empty()) {
        headerText = "@HD\tVN:1.6\n";
    }
    for (size_t i = 0; i < refs.size(); ++i) {
        refIndex[refs[i].name] = (int32_t)i;
    }

    std::vector<uint8_t> hdr;
    hdr.reserve(64 + headerText.size() + refs.size() * 32);
    hdr.push_back('B');
    hdr.push_back('A');
    hdr.push_back('M');
    hdr.push_back(1);
    putI32(hdr, (int32_t)headerText.size());
    hdr.insert(hdr.end(), headerText.begin(), headerText.end());
    putI32(hdr, (int32_t)refs.size());
    for (const BamRef& r : refs) {
        putI32(hdr, (int32_t)(r.name.size() + 1));
        hdr.insert(hdr.end(), r.name.begin(), r.name.end());
        hdr.push_back(0);
        putI32(hdr, r.len);
    }

    dataStart = pos;
    return bgzfWrite(hdr.data(), hdr.size());
}

int32_t BamWriter::writeBamRecord(const uint8_t* line, size_t len) {
    /* 去掉行尾换行/回车 */
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        --len;
    }
    if (len == 0) {
        return 0;
    }

    /* 按 \t 切字段 */
    std::vector<const char*> fields;
    std::vector<size_t> flens;
    {
        size_t start = 0;
        for (size_t i = 0; i <= len; ++i) {
            if (i == len || line[i] == '\t') {
                fields.push_back((const char*)line + start);
                flens.push_back(i - start);
                start = i + 1;
            }
        }
    }
    if (fields.size() < 11) {
        LOG_WARNING("BamWriter: SAM record with too few fields, skipped.");
        return 0;
    }

    const char* qname = fields[0];
    const size_t qlen = flens[0];
    const char* rname = fields[2];
    const size_t rlen = flens[2];
    const char* rnext = fields[6];
    const size_t rnextLen = flens[6];
    const char* cigarStr = fields[5];
    const size_t clen = flens[5];
    const char* seqStr = fields[9];
    const size_t slen = flens[9];
    const char* qualStr = fields[10];
    const size_t qlen2 = flens[10];

    if (qlen == 0) {
        return 0;   // QNAME 不能为空，非法记录跳过
    }

    const uint16_t flag = (uint16_t)strtoul(fields[1], nullptr, 10);
    const uint8_t mapq = (uint8_t)strtoul(fields[4], nullptr, 10);
    const int64_t pos1 = strtoll(fields[3], nullptr, 10);      // 1-based
    const int64_t pnext1 = strtoll(fields[7], nullptr, 10);    // 1-based
    const int32_t btlen = (int32_t)strtol(fields[8], nullptr, 10);

    /* refID：RNAME '*' / 空 -> -1，否则查参考序列表 */
    int32_t refId = -1;
    if (rlen > 0 && !(rlen == 1 && rname[0] == '*')) {
        const std::string nm(rname, rlen);
        std::map<std::string, int32_t>::const_iterator it = refIndex.find(nm);
        if (it != refIndex.end()) {
            refId = it->second;
        } else {
            LOG_WARNING("BamWriter: RNAME '%s' not in @SQ, treat as unmapped.", nm.c_str());
        }
    }
    int32_t bpos = (refId >= 0 && pos1 > 0) ? (int32_t)(pos1 - 1) : -1;

    int32_t nextRefId = -1;
    if (rnextLen == 1 && rnext[0] == '=') {
        nextRefId = refId;   // '=' 表示与 RNAME 同参考序列
    } else if (rnextLen > 0 && !(rnextLen == 1 && rnext[0] == '*')) {
        const std::string nm(rnext, rnextLen);
        std::map<std::string, int32_t>::const_iterator it = refIndex.find(nm);
        if (it != refIndex.end()) {
            nextRefId = it->second;
        }
    }
    int32_t bnext = (nextRefId >= 0 && pnext1 > 0) ? (int32_t)(pnext1 - 1) : -1;

    /* CIGAR */
    std::vector<uint32_t> cigarOps;
    int64_t refSpan = 0;
    if (!(clen == 1 && cigarStr[0] == '*')) {
        if (!parseCigar(cigarStr, clen, cigarOps, refSpan)) {
            LOG_WARNING("BamWriter: invalid CIGAR, treat as no CIGAR.");
            cigarOps.clear();
            refSpan = 0;
        }
    }
    const uint16_t nCigar = (uint16_t)cigarOps.size();

    /* bin */
    uint16_t bin = 0;
    if (refId >= 0 && bpos >= 0) {
        const int64_t end = (int64_t)bpos + (refSpan > 0 ? refSpan : 1);
        bin = samReg2Bin(bpos, end);
    }

    /* SEQ */
    int32_t lSeq = 0;
    std::vector<uint8_t> packedSeq;
    if (slen > 0 && !(slen == 1 && seqStr[0] == '*')) {
        lSeq = (int32_t)slen;
        packedSeq.resize((size_t)((slen + 1) / 2), 0);
        for (size_t i = 0; i < slen; ++i) {
            const int nib = baseNibble(seqStr[i]);
            if (i & 1) {
                packedSeq[i >> 1] |= (uint8_t)(nib & 0xF);
            } else {
                packedSeq[i >> 1] = (uint8_t)((nib & 0xF) << 4);
            }
        }
    }

    /* QUAL：缺失（*）时全 0xFF；否则逐字节 ascii-33 */
    std::vector<uint8_t> quals;
    const bool qualMissing = (qlen2 == 0) || (qlen2 == 1 && qualStr[0] == '*');
    if (!qualMissing) {
        size_t qn = (qlen2 < (size_t)lSeq) ? qlen2 : (size_t)lSeq;
        quals.resize(qn);
        for (size_t i = 0; i < qn; ++i) {
            int phred = (int)(uint8_t)qualStr[i] - 33;
            if (phred < 0) {
                phred = 0;
            }
            if (phred > 93) {
                phred = 93;
            }
            quals[i] = (uint8_t)phred;
        }
    }

    /* read_name（BAM 规定 l_read_name <= 255，含结尾 \0） */
    size_t qn = (qlen < 254) ? qlen : 254;
    const uint8_t lReadName = (uint8_t)(qn + 1);

    /* 组装记录 */
    std::vector<uint8_t> rec;
    rec.reserve(32 + qn + cigarOps.size() * 4 + packedSeq.size() + quals.size() + 16);
    putI32(rec, 0);   // block_size 占位
    putI32(rec, refId);
    putI32(rec, bpos);
    putU8(rec, lReadName);
    putU8(rec, mapq);
    putU16(rec, bin);
    putU16(rec, nCigar);
    putU16(rec, flag);
    putI32(rec, lSeq);
    putI32(rec, nextRefId);
    putI32(rec, bnext);
    putI32(rec, btlen);
    rec.insert(rec.end(), qname, qname + qn);
    rec.push_back(0);
    for (size_t i = 0; i < cigarOps.size(); ++i) {
        putU32(rec, cigarOps[i]);
    }
    rec.insert(rec.end(), packedSeq.begin(), packedSeq.end());
    /* QUAL：有值时逐字节追加；缺失（*）或短于 SEQ 时补 0xFF，保证记录长度与 BAM 布局一致 */
    if (!quals.empty()) {
        rec.insert(rec.end(), quals.begin(), quals.end());
    }
    for (int32_t i = (int32_t)quals.size(); i < lSeq; ++i) {
        rec.push_back(0xFF);
    }
    for (size_t i = 11; i < fields.size(); ++i) {
        if (0 != appendSamAuxToBam(fields[i], flens[i], rec)) {
            LOG_WARNING("BamWriter: skip invalid SAM option field.");
        }
    }

    const int32_t blockSize = (int32_t)(rec.size() - 4);
    memcpy(rec.data(), &blockSize, 4);

    return bgzfWrite(rec.data(), rec.size());
}

int32_t BamWriter::writeDataLines(const uint8_t* buffer, size_t start, size_t end) {
    size_t pos = start;
    while (pos < end) {
        size_t nl = pos;
        while (nl < end && buffer[nl] != '\n') {
            ++nl;
        }
        if (nl > pos && buffer[pos] != '@') {
            if (0 != writeBamRecord(buffer + pos, nl - pos)) {
                return -1;
            }
        }
        pos = nl + 1;
    }
    return 0;
}

int32_t BamWriter::writeBlock(RoughIOBlock* blockPtr) {
    if (ioWriter == nullptr || blockPtr == nullptr) {
        return -1;
    }
    const uint8_t* buffer = blockPtr->getBuffer();
    const int64_t dataLen = blockPtr->getDataLen();
    if (dataLen <= 0) {
        return 0;
    }

    if (!passThrough && !headerWritten) {
        const BlockType type = blockPtr->getBlockType();
        if (BlockUtil::isSAMBlock(type) || BlockUtil::isBAMBlock(type)) {
            size_t dataStart = (size_t)dataLen;
            if (0 != writeBamHeader(buffer, (size_t)dataLen, dataStart)) {
                return -1;
            }
            headerWritten = true;
            /* 头部行之后可能还跟着数据行（头部未独立成块的块） */
            if (dataStart < (size_t)dataLen) {
                return writeDataLines(buffer, dataStart, (size_t)dataLen);
            }
            return 0;
        }
        LOG_WARNING("BamWriter: block type %d is not SAM/BAM, pass through raw output.", type);
        passThrough = true;
        /* 落到下面的透传分支，把当前块也原样写出 */
    }

    if (passThrough) {
        return writeRaw(buffer, (size_t)dataLen);
    }

    return writeDataLines(buffer, 0, (size_t)dataLen);
}



