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
#include "profile_stats.h"
#include "log/logger.h"
#include "pbgz_types.h"
#include <cstring>
#include <set>
#include <atomic>
#include <thread>
#include "pbgz_manager.h"
#include "config_manager.h"
#include <zlib.h>

#include "bam_record.h"

#if defined(__x86_64__)
#include <immintrin.h>
#include "utils/md5_util.h"
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

BlockType BlockReader::analyzeBlock(RoughIOBlock* /*blockPtr*/, BlockType /*fileType*/) {
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

    /* Grow the input block on demand: it initially allocates only 1MB, and is realloc'd when the read target (-l) is larger */
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

    // Subclass performs format pre-analysis (newline position recording, Fastq/Sam actuator pre-analysis, etc.)
    BlockType type = analyzeBlock(blockPtr, fileType);
    blockPtr->setBlockType(type);

    LOG_DEBUG("Read One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    
    return blockPtr->getDataLen();
}

size_t BlockReader::prependBufferedData(RoughIOBlock* blockPtr, size_t& totalLen, bool capAtBlockSize) {
    uint8_t* buffer = blockPtr->getBuffer();
    /* Cap by the read target (-l) rather than the block's own blockSize (the input block only allocates 1M; the read target may be larger or smaller) */
    const size_t capSize = capAtBlockSize ? readTargetBytes(blockPtr) : (size_t)-1;
    size_t placed = 0;

    /* Merge the previous block's leftover cache first: it precedes the remaining detect data in the file. */
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

    // Merge in the format-detection data prefetched by the Creator (appended after cache)
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
                // The whole block cannot hold even one complete FASTQ record; cannot align
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
     * Tail block (EOF): structured coding needs a record-aligned end (\n-terminated, FASTQ line count a multiple of 4).
     * Malformed data is never tampered with; the block degrades to BINARY instead.
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

    // Align the tail of the block to a record boundary (move the incomplete record into cache)
    if (!alignToRecordBoundary(blockPtr, true)) {
        return BINARY;
    }

    /* Only record newline positions and validate structure; pre-analysis runs in the coder-side actuator. */
    return baseType;
}

/*
 * CRAM-style block slicing for SAM/BAM: a block ends at the read-count limit (readsPerBlock) or the base cap
 * (blockMaxBases = readsPerBlock * 500), whichever comes first.
 */

/* Sequenced bases in one SAM record line (SEQ column, field 10); 0 for headers, '*' and malformed rows. */
static size_t samRowSeqLength(const uint8_t* row, size_t len) {
    int tabs = 0;
    size_t seqStart = 0;
    size_t i = 0;
    for (; i < len; ++i) {
        if (row[i] == '\t') {
            ++tabs;
            if (tabs == 9) {
                seqStart = i + 1;
            } else if (tabs == 10) {
                break;   /* i now points at the tab that terminates SEQ */
            }
        }
    }
    if (tabs < 9) {
        return 0;
    }
    const size_t seqEnd = (tabs >= 10) ? i : len;   /* no 10th tab -> SEQ runs to the line end */
    if (seqEnd <= seqStart) {
        return 0;
    }
    if (seqEnd - seqStart == 1 && row[seqStart] == '*') {
        return 0;
    }
    return seqEnd - seqStart;
}

/*
 * SAM/SAM-GZ block reading: the @ header is its own block; the data region is split into blocks of readsPerBlock reads.
 */
int64_t SamBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    if (ioReader == nullptr || blockPtr == nullptr) {
        LOG_ERROR("IO reader or block pointer is null.");
        return -1;
    }

    blockPtr->setBlockId(blockId++);
    std::vector<size_t>& npos = blockPtr->getNpos();
    npos.clear();

    size_t totalLen = 0;
    /* SAM blocks are split by read line count, unconstrained by byte blockSize; merge all prefetch data and cache at once */
    prependBufferedData(blockPtr, totalLen, false);

    uint8_t* buffer = blockPtr->getBuffer();
    const size_t chunkSize = 1 << 20;

    /* Empty buffer: read a chunk first so we can tell whether it starts with @ header lines. */
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
    /* Only a "pure header block" carries no data lines; all other blocks (data blocks, degraded blocks, blocks with an unsplit header) carry data */
    lastBlockHasData = !(startsWithHeader && splitHeader);

    /* Header block (splitHeader and @ start): stop at the first data line. Otherwise read up to readsPerBlock data lines. */
    const size_t readTarget = (splitHeader && startsWithHeader) ? 1 : readsPerBlock;

    size_t dataLineCount = 0;
    size_t blockBases = 0;    /* sequenced bases already kept in this block */
    size_t lineStart = 0;
    size_t scanPos = 0;

    /* Slice on reads or bases reached; a line overshooting the cap starts the next block, and reads are never split. */
    bool reachedBaseCap = false;
    /* Only EOF marks a truncated tail; a buffer-cap stop is aligned to the last line below. */
    bool eofReached = false;

    while (dataLineCount < readTarget) {
        // First scan the bytes already present in the current buffer
        size_t i = scanPos;
        while (i < totalLen) {
            if (buffer[i] == '\n') {
                const bool isDataLine = (buffer[lineStart] != '@');
                if (isDataLine) {
                    const size_t rowBases = samRowSeqLength(buffer + lineStart, i - lineStart);
                    if (dataLineCount > 0 && blockMaxBases > 0 && blockBases + rowBases > blockMaxBases) {
                        reachedBaseCap = true; /* this line and everything after it belongs to the next block */
                        break;
                    }
                    blockBases += rowBases;
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
        if (reachedBaseCap || dataLineCount >= readTarget) {
            break;
        }

        // Need more data: grow the buffer on demand and continue reading
        if (totalLen + chunkSize > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(totalLen + chunkSize)) {
                break;
            }
            buffer = blockPtr->getBuffer();
        }
        const size_t readLen = ioReader->readIO(buffer + totalLen, chunkSize);
        if (readLen == 0) {
            eofReached = true;
            break; // EOF
        }
        totalLen += readLen;
    }

    if (totalLen == 0) {
        return 0; // EOF
    }

    // ---- Whether the end of the block falls on a complete line boundary ----
    const bool stoppedAtTarget = (dataLineCount >= readTarget);
    const bool endsWithNewline = (buffer[totalLen - 1] == '\n');
    if (!stoppedAtTarget && !endsWithNewline && eofReached) {
        /* Truncated trailing record (no ending newline): compress the whole block as binary. */
        LOG_INFO("SAM tail block not aligned to record boundary, fallback to binary codec (blockId=%ld).",
                 (long)blockPtr->getBlockId());
        blockPtr->setDataLen((int64_t)totalLen);
        blockPtr->setBlockType(BINARY);
        return (int64_t)totalLen;
    }

    // ---- Determine up to which line this block keeps ----
    size_t keepLines = 0;
    if (startsWithHeader && splitHeader) {
        /* Header block: keep only all @ header lines; data lines are left for subsequent blocks */
        size_t lineStartPos = 0;
        for (size_t idx = 0; idx < npos.size(); ++idx) {
            if (buffer[lineStartPos] != '@') {
                break;
            }
            keepLines = idx + 1;
            lineStartPos = npos[idx] + 1;
        }
    } else {
        /* Data block: keep header lines (if any) + readsPerBlock data lines (keep all if fewer remain) */
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
        /* No complete line found; treat the whole block as binary. */
        LOG_INFO("SAM block has no complete line, fallback to binary codec (blockId=%ld).",
                 (long)blockPtr->getBlockId());
        blockPtr->setDataLen((int64_t)totalLen);
        blockPtr->setBlockType(BINARY);
        return (int64_t)totalLen;
    }

    // ---- Align the block tail: move the content after the kept lines into cache for the next block ----
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

    /* Only record newline positions and validate structure; pre-analysis runs in the coder-side actuator. */
    blockPtr->setBlockType(SAM);
    return (int64_t)totalLen;
}

BlockType SamBlockReader::analyzeBlock(RoughIOBlock* blockPtr, BlockType /*fileType*/) {
    uint8_t* buffer = blockPtr->getBuffer();
    const int64_t dataLen = blockPtr->getDataLen();
    const char* bufPtr = reinterpret_cast<const char*>(buffer);

    // Record newline positions (memchr has built-in SIMD optimization)
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

    // Align the block tail to a complete line
    if (!alignToRecordBoundary(blockPtr, false)) {
        return BINARY;
    }

    /* Fallback for structural validation; deep SAM validation and BINARY fallback happen in the coder-side actuator. */
    return SAM;
}

namespace {

    /* ---- BAM binary field parsing (little-endian) ---- */
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
    /* BAM 4-bit base encoding: index 0..15 -> "=ACMGRSVTWYHKDBN" */
    const char BAM_BASE_MAP[16] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                   'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

    void buildBamCigar(const uint8_t* data, uint16_t nCigarOp, std::string& out) {
        for (uint16_t i = 0; i < nCigarOp; ++i) {
            uint32_t op;
            memcpy(&op, data + i * 4, 4);
            const uint32_t len = op >> 4;      /* Upper 28 bits are the length */
            const unsigned code = op & 0xF;    /* Lower 4 bits are the opcode */
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
            return;  /* The caller uniformly outputs '*' */
        }
        /* BAM QUAL stores raw phred scores; adding 33 recovers SAM's phred+33. */
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
            /* SAM auxiliary types are only A/i/f/Z/H/B; BAM integer subtypes print as 'i'. */
            line += 'i';
            line += ":";
            switch (type) {
            case 'A':
                /* Character type; SAM uses 'A' */
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
                return;   /* Unknown type; skip the rest */
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
    /* Single non-blocking read; callers needing a full read must loop (see readBamBytes). */
    if (got < n) {
        got += ioReader->readIO((uint8_t*)dst + got, n - got);
    }
    return got;
}

size_t BamBlockReader::readBamBytes(void* dst, size_t n) {
    size_t got = 0;
    uint8_t* out = (uint8_t*)dst;
    while (got < n) {
        const size_t nRead = readRawFromSource(out + got, n - got);
        if (nRead == 0) {
            break; // EOF
        }
        got += nRead;
    }
    return got;
}

size_t BamGzBlockReader::readBamBytes(void* dst, size_t n) {
    PBGZ_PROF_SCOPE(pbgzprof::READ_INFLATE);
    if (!parStarted) {
        parStarted = true;
        /* Probe the stream head once; the bytes are replayed by the parallel or serial decoder. */
        probeBuf.resize(1 << 16);
        size_t got = 0;
        while (got < probeBuf.size()) {
            const size_t r = readRawFromSource(probeBuf.data() + got, probeBuf.size() - got);
            if (r == 0) {
                break;
            }
            got += r;
        }
        probeBuf.resize(got);

        if (got >= 12) {
            auto rawFn = [](void* c, void* d, size_t len) -> size_t {
                return static_cast<BamGzBlockReader*>(c)->readRawFromSource(d, len);
            };
            if (bgzfPar.start(this, rawFn, probeBuf.data(), probeBuf.size(), parThreads)) {
                parMode = true;
            }
        }
    }

    if (parMode) {
        return bgzfPar.read(dst, n);
    }

    /* Serial fallback for plain gzip; the compressed input first comes from the probe bytes above. */
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
            /* gzip concatenates multiple members; reset and continue inflating. */
            if (inflateState.avail_in > 0) {
                inflateReset(&inflateState);
                continue;
            }
            if (gzInEof) {
                break;
            }
            const size_t inRead = readRawBam(gzInBuf, sizeof(gzInBuf));
            if (inRead == 0) {
                gzInEof = true;
                break;
            }
            inflateReset(&inflateState);
            inflateState.next_in = gzInBuf;
            inflateState.avail_in = (uInt)inRead;
            continue;
        }
        if (rc != Z_OK && rc != Z_BUF_ERROR) {
            break;   /* Corrupted data */
        }
        if (inflateState.avail_in == 0) {
            if (gzInEof) {
                break;
            }
            const size_t inRead = readRawBam(gzInBuf, sizeof(gzInBuf));
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

/* Raw source for the serial fallback: drain the probe bytes first, then the underlying source. */
size_t BamGzBlockReader::readRawBam(void* dst, size_t n) {
    uint8_t* out = (uint8_t*)dst;
    size_t got = 0;
    while (got < n && !probeBuf.empty()) {
        const size_t c = (probeBuf.size() < n - got) ? probeBuf.size() : (n - got);
        memcpy(out + got, probeBuf.data(), c);
        probeBuf.erase(probeBuf.begin(), probeBuf.begin() + (ptrdiff_t)c);
        got += c;
    }
    if (got < n) {
        got += readRawFromSource(out + got, n - got);
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

    /* Build on the BAM header text, appending @SQ lines for reference sequences that lack them */
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

int32_t bamRecordToSamLine(const uint8_t* data, int32_t size,
                           const std::vector<BamRef>& refs, std::string& line) {
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
 * Worker-side render of a raw-record BAM block into byte-identical SAM text, recomputing npos and the MD5.
 */
int32_t bamRenderRawBlockToText(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr) {
        return -1;
    }
    const std::shared_ptr<const std::vector<BamRef>> refs = blockPtr->getBamRefs();
    if (refs == nullptr) {
        LOG_ERROR("Raw BAM block(%ld) has no reference table for rendering.",
                  (long)blockPtr->getBlockId());
        return -1;
    }

    const uint8_t* src = blockPtr->getBuffer();
    const size_t srcLen = (size_t)blockPtr->getDataLen();

    /* Render into a grow-on-demand buffer, then swap it in as the block content. */
    std::string text;
    text.reserve(srcLen * 3 / 2 + 64);

    std::vector<size_t> newNpos;
    size_t pos = 0;
    while (pos + 4 <= srcLen) {
        int32_t size = 0;
        memcpy(&size, src + pos, 4);
        if (size <= 0 || (size_t)size > srcLen - pos - 4 || size > (1 << 30)) {
            LOG_ERROR("Raw BAM block(%ld): bad record length %d at offset %zu "
                      "(srcLen=%zu, first8=%02x %02x %02x %02x %02x %02x %02x %02x).",
                      (long)blockPtr->getBlockId(), size, pos, srcLen,
                      (unsigned)src[0], (unsigned)src[1], (unsigned)src[2], (unsigned)src[3],
                      (unsigned)src[4], (unsigned)src[5], (unsigned)src[6], (unsigned)src[7]);
            return -1;
        }
        std::string line;
        if (0 != bamRecordToSamLine(src + pos + 4, size, *refs, line)) {
            LOG_ERROR("Raw BAM block(%ld): failed to render record at offset %zu.",
                      (long)blockPtr->getBlockId(), pos);
            return -1;
        }
        const size_t lineStart = text.size();
        text.append(line);
        text.push_back('\n');
        newNpos.push_back(lineStart + line.size());
        pos += 4 + (size_t)size;
    }
    if (pos != srcLen) {
        LOG_ERROR("Raw BAM block(%ld): trailing %zu bytes after the last record.",
                  (long)blockPtr->getBlockId(), srcLen - pos);
        return -1;
    }

    /* Replace the block payload with the rendered text. */
    if (text.size() > (size_t)blockPtr->getBufferSize()) {
        if (0 != blockPtr->ensureCapacity(text.size())) {
            return -1;
        }
    }
    memcpy(blockPtr->getBuffer(), text.data(), text.size());
    blockPtr->setDataLen((int64_t)text.size());
    blockPtr->getNpos().swap(newNpos);

    /* MD5 of the rendered text: the reader skipped hashing raw blocks, so it is computed here. */
    std::string md5;
    calcMd5sum(md5, blockPtr->getBuffer(), (uint32_t)blockPtr->getDataLen());
    blockPtr->setMd5(md5);
    blockPtr->setRawBam(false);

    /* Count the rendered (text) size for the file-level ratio. */
    PbgzManager::getInstance().updateReadDataLen(blockPtr);
    return 0;
}

int32_t BamBlockReader::parseBamRecord(const uint8_t* data, int32_t size, std::string& line) {
    return bamRecordToSamLine(data, size, refs, line);
}



/*
 * BAM reading: the header is its own block (as SAM text); reads are decompressed into SAM lines and grouped into readsPerBlock blocks. The block type stays BAM.
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

    // Write the header (if not yet written)
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

    if (structMode) {
        colsCur = std::make_shared<BamColumns>();
        colsCur->recordStreamOffset = (uint32_t)outLen;
    }

    // Data region: decompress each read into a SAM line until the read or base cap is reached.
    uint32_t reads = 0;
    size_t blockBases = 0;   /* sequenced bases already kept in this block */
    bool blockIsRaw = false; /* this block streamed raw records (worker renders) */
    /* Whether this block is past the text-rendered prefix, i.e. streams raw records. */
    const bool rawThisBlock = rawStreamMode && (int64_t)blockId > rawStreamTextLimit;

    /* A read that would push the block over the base cap is parked for the next block; reads are never split. */
    if (!pendingSamLine.empty() && !rawThisBlock) {
        if (outLen + pendingSamLine.size() + 1 > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(outLen + pendingSamLine.size() + 1)) {
                return -1;
            }
            out = blockPtr->getBuffer();
        }
        memcpy(out + outLen, pendingSamLine.data(), pendingSamLine.size());
        outLen += pendingSamLine.size();
        out[outLen] = '\n';
        npos.push_back(outLen);
        outLen += 1;
        blockBases += samRowSeqLength((const uint8_t*)pendingSamLine.data(), pendingSamLine.size());
        pendingSamLine.clear();
        ++reads;
        lastBlockHasData = true;
    }

    if (structMode && !pendingBamRecord.empty()) {
        /* A record parked to keep the previous block inside its base cap; streams in like any other. */
        std::vector<uint8_t> pend;
        pend.swap(pendingBamRecord);
        if (outLen + 4 + pend.size() > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(outLen + 4 + pend.size())) {
                return -1;
            }
            out = blockPtr->getBuffer();
        }
        int32_t pendSize = (int32_t)pend.size();
        memcpy(out + outLen, &pendSize, 4);
        outLen += 4;
        memcpy(out + outLen, pend.data(), pend.size());
        outLen += pend.size();
        ++reads;
        lastBlockHasData = true;
    } else if (rawThisBlock && !pendingBamRecord.empty()) {
        /* Same parked record in raw-stream mode: append it raw and count its bases like the text flush. */
        blockIsRaw = true;
        std::vector<uint8_t> pend;
        pend.swap(pendingBamRecord);
        if (outLen + 4 + pend.size() > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(outLen + 4 + pend.size())) {
                return -1;
            }
            out = blockPtr->getBuffer();
        }
        int32_t pendSize = (int32_t)pend.size();
        memcpy(out + outLen, &pendSize, 4);
        outLen += 4;
        memcpy(out + outLen, pend.data(), pend.size());
        outLen += pend.size();
        int32_t pendSeqLen = 0;
        if (pend.size() >= 20) {
            const uint8_t* pCore = pend.data() + 16;
            pendSeqLen = bamI32(pCore);
            if (pendSeqLen < 0) {
                pendSeqLen = 0;
            }
        }
        blockBases += (unsigned long)(pendSeqLen > 0 ? (unsigned long)pendSeqLen : 0);
        ++reads;
        lastBlockHasData = true;
    }

    while (reads < readsPerBlock) {
        int32_t blockSize = 0;
        if (readBamBytes(&blockSize, 4) != 4) {
            break;   // EOF
        }
        if (blockSize <= 0 || blockSize > (1 << 30)) {
            break;
        }
        recBuf.resize((size_t)blockSize);
        if (readBamBytes(recBuf.data(), recBuf.size()) != recBuf.size()) {
            break;
        }
        /* l_seq (SEQ length) is the fixed field of the 32-byte BAM core at offset 16. */
        int32_t seqLen = 0;
        if (recBuf.size() >= 20) {
            const uint8_t* pCore = recBuf.data() + 16;
            seqLen = bamI32(pCore);
            if (seqLen < 0) {
                seqLen = 0;
            }
        }
        if (structMode) {
            if (reads > 0 && blockMaxBases > 0 &&
                (uint64_t)blockBases + (uint64_t)(seqLen > 0 ? (uint32_t)seqLen : 0) > blockMaxBases) {
                pendingBamRecord.assign(recBuf.begin(), recBuf.end());
                break;
            }
            /* Raw record bytes only ([u32 size][record]); column parsing is deferred to the worker. */
            if (outLen + 4 + (size_t)blockSize > blockPtr->getBufferSize()) {
                if (0 != blockPtr->ensureCapacity(outLen + 4 + (size_t)blockSize)) {
                    break;
                }
                out = blockPtr->getBuffer();
            }
            memcpy(out + outLen, &blockSize, 4);
            outLen += 4;
            memcpy(out + outLen, recBuf.data(), (size_t)blockSize);
            outLen += (size_t)blockSize;
            blockBases += (uint64_t)(seqLen > 0 ? (uint32_t)seqLen : 0);
            ++reads;
            lastBlockHasData = true;
            continue;
        }

        if (rawThisBlock) {
            /* Raw-stream mode: after the leading text blocks, stream raw records; the worker renders them. */
            if (reads > 0 && blockMaxBases > 0 &&
                (uint64_t)blockBases + (uint64_t)(seqLen > 0 ? (uint32_t)seqLen : 0) > blockMaxBases) {
                pendingBamRecord.assign(recBuf.begin(), recBuf.end());
                break;
            }
            blockIsRaw = true;
            if (outLen + 4 + (size_t)blockSize > blockPtr->getBufferSize()) {
                if (0 != blockPtr->ensureCapacity(outLen + 4 + (size_t)blockSize)) {
                    break;
                }
                out = blockPtr->getBuffer();
            }
            memcpy(out + outLen, &blockSize, 4);
            outLen += 4;
            memcpy(out + outLen, recBuf.data(), (size_t)blockSize);
            outLen += (size_t)blockSize;
            blockBases += (uint64_t)(seqLen > 0 ? (uint32_t)seqLen : 0);
            ++reads;
            lastBlockHasData = true;
            continue;
        }

        lineBuf.clear();
        if (0 != parseBamRecord(recBuf.data(), blockSize, lineBuf)) {
            break;
        }
        /* Slice on reads or bases reached; park this read rather than overshoot the cap. */
        if (reads > 0 && blockMaxBases > 0 && (uint64_t)blockBases + (uint64_t)seqLen > blockMaxBases) {
            if (rawStreamMode && ((int64_t)blockId + 1 > rawStreamTextLimit)) {
                /* Next block is raw-stream: park the BAM record itself so raw framing stays intact. */
                pendingBamRecord.assign(recBuf.begin(), recBuf.end());
            } else {
                pendingSamLine.assign(lineBuf);   // copy keeps lineBuf capacity reusable
            }
            break;
        }
        if (outLen + lineBuf.size() + 1 > blockPtr->getBufferSize()) {
            if (0 != blockPtr->ensureCapacity(outLen + lineBuf.size() + 1)) {
                break;
            }
            out = blockPtr->getBuffer();
        }
        memcpy(out + outLen, lineBuf.data(), lineBuf.size());
        outLen += lineBuf.size();
        out[outLen] = '\n';
        npos.push_back(outLen);
        outLen += 1;
        blockBases += (uint64_t)(seqLen > 0 ? (uint32_t)seqLen : 0);
        ++reads;
        lastBlockHasData = true;
    }

    if (outLen == 0) {
        return 0;   // EOF (the header was already emitted in a previous block)
    }

    blockPtr->setDataLen((int64_t)outLen);
    blockPtr->setBlockType(BAM);
    if (blockIsRaw) {
        /* Raw records: npos stays empty and the worker hashes the text after rendering. */
        blockPtr->setRawBam(true);
        if (!rawRefs) {
            rawRefs = std::make_shared<const std::vector<BamRef>>(refs);
        }
        blockPtr->setBamRefs(rawRefs);
    }
    if (structMode && colsCur != nullptr) {
        colsCur->nRecords = reads;   /* payload parsed later, by the worker */
        colsCur->columnsBuilt = false;
        blockPtr->setBamColumns(colsCur);
    }
    return (int64_t)outLen;
}

namespace BlockUtil {
    bool isFastqBlock(BlockType type) {
        return (type == FASTQ_GEN2 || type == FASTQ_GEN3 || type == FASTQ_GEN2_GZIP || type == FASTQ_GEN3_GZIP);
    }

    bool isSAMBlock(BlockType type) {
        /* BAM blocks contain SAM text and are treated as SAM, so every capability supported for SAM is also supported for BAM */
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

    /* BGZF magic: 1f 8b 08 04 (BAM itself is a BGZF stream) */
    bool isBgzfStream(const uint8_t* buf, size_t len) {
        return len >= 4 && buf[0] == 0x1f && buf[1] == 0x8b && buf[2] == 0x08 && (buf[3] & 0x04) == 0x04;
    }

    /* Inflate the first gzip/BGZF block; on success out holds the inflated bytes and true is returned */
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

    /* Check whether the inflated data starts with the BAM magic number */
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

    /* Validate every line in the sample; any malformed line makes the whole sample BINARY. */
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
                return false;   // An empty line is not a valid SAM line
            }
            if (buf[lineStart] == '@') {
                // Header line: must be a known SAM header type
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
                // Data line: SAM requires at least 11 mandatory fields (QNAME..QUAL), i.e. at least 10 tabs
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

    /* Validate the FASTQ structure of every line; any malformed line makes the sample BINARY, and one full record is required. */
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
 * Whether a file is BAM: raw BAM magic, or a BGZF/gzip stream inflating to BAM.
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

/* Detect the input format from file content, inflating gz/BGZF first. Regular files only; not for pipes/STDIN. */
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
    /* Generic gzip magic (does not require BGZF's FEXTRA flag; plain .gz must also be inflatable) */
    if (len >= 2 && buf[0] == 0x1f && buf[1] == 0x8b) {
        std::vector<uint8_t> plain;
        if (!inflateFirstBlock(buf.data(), len, plain)) {
            return BINARY;
        }
        /* After inflating the outer gzip the result may still be BGZF (e.g. a gzipped BAM); inflate one more layer */
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
 * SAM block granularity per level (1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000 reads), used by the readers and the codec pre-selection trial.
 */
uint32_t BlockFactory::samReadsPerBlockOfLevel(uint8_t compressLevel)
{
    if (compressLevel >= 8) {
        return 100000;
    }
    if (compressLevel >= 6) {
        return 25000;
    }
    return 10000;
}

BlockReader* BlockFactory::createBlockReader(IOReader* ioReader, uint8_t compressLevel, bool splitSamHeader,
                                              bool bamStructMode, bool pbgzAsOpaque, uint32_t bgzfThreads) {
    if (ioReader == nullptr) {
        LOG_ERROR("Create block reader failed: io reader is null.");
        return nullptr;
    }

    /* SAM data blocks are split by read line count: 1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000 */
    uint32_t samReadsPerBlock = samReadsPerBlockOfLevel(compressLevel);

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
        if (pbgzAsOpaque) {
            /*
             * Compressing an archive: read it as opaque binary, since its blocks hold already-coded payloads.
             */
            reader = MemoryUtil::safeNewClass<BinaryBlockReader>(ioReader, detectBuf, detectLen);
        } else {
            /*
             * PBGZ: prefetched bytes are returned to PbgzBlockReader (piped input cannot be seeked back).
             */
            reader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader, detectBuf, detectLen);
        }
    } else if (isBamSample(detectBuf, detectLen)) {
        if (startsWithBamMagic(detectBuf, detectLen)) {
            /* Already a raw BAM byte stream (the io layer has inflated BGZF, e.g. .bam) */
            BamBlockReader* bamReader = MemoryUtil::safeNewClass<BamBlockReader>(ioReader, detectBuf, detectLen, samReadsPerBlock, splitSamHeader);
            if (bamReader != nullptr && bamStructMode) {
                bamReader->setStructMode(true);
            }
            reader = bamReader;
        } else {
            /* The inner layer is still BGZF (e.g. a plain .bam, or .bam.gz); inflate it to raw BAM first */
            BamGzBlockReader* bamReader = MemoryUtil::safeNewClass<BamGzBlockReader>(ioReader, detectBuf, detectLen, samReadsPerBlock, splitSamHeader);
            if (bamReader != nullptr && bamStructMode) {
                bamReader->setStructMode(true);
            }
            if (bamReader != nullptr && bgzfThreads > 0) {
                /* Parallel inflate workers, derived from -t by the caller (see bamBgzfThreadNum). */
                bamReader->setParThreads((int)bgzfThreads);
            }
            reader = bamReader;
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
    /* The read target is set by -l (the input block only allocates 1MB and is grown to this target on demand) */
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
    if (0 != pbgzFileReader->readDataBlock(pbgzDataBlock, blockPtr)) {        LOG_ERROR("Read Pbgz data block failed.");
        return -1;
    }

    if (0 != pbgzDataBlock.verifyCheckSum()) {
        LOG_ERROR("Verify pbgz block checksum failed, block id = %ld", pbgzDataBlock.getMetaData("blockid").asInt64());
        return -1;
    }
    
    blockPtr->setDataLen(pbgzDataBlock.getMetaData("datalen").asInt64());
    blockPtr->setMetaLen(pbgzDataBlock.getMetaData("metalen").asInt64());
    blockPtr->setBlockId(pbgzDataBlock.getMetaData("blockid").asInt64());
    /* Package identity is intrinsic to the block and written by its producer, not at each call site. */
    blockPtr->setPackageStart((int64_t)pbgzFileReader->getCurrentFileStart());
    blockPtr->setPackageIndex(pbgzFileReader->getCurrentFileIndex());
    std::string blockType = pbgzDataBlock.getMetaData("blocktype").asString();
    if (blockType == "fastq_gen2" || blockType == "fastq_gen2_gzip") {
        /* In old files the *_gzip block type marks "the raw input was GZ"; uniformly treat it as the non-GZ type */
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

/* ==================== BamWriter (decompress -b: SAM -> BAM) ==================== */

namespace {

    /* Extract the value of a TAG field from an @SQ line ("SN:" / "LN:") */
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
     * Parse the leading SAM header lines and append the BAM header byte stream to `out`.
     * dataStart is the offset of the first data line; refIndex is filled with the @SQ order.
     */
    int32_t buildBamHeaderPayload(const uint8_t* data, size_t len, size_t& dataStart,
                                  std::map<std::string, int32_t>& refIndex,
                                  std::vector<uint8_t>& out)
    {
        refIndex.clear();
        std::vector<BamRef> refs;
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

        out.clear();
        out.reserve(64 + headerText.size() + refs.size() * 32);
        out.push_back('B');
        out.push_back('A');
        out.push_back('M');
        out.push_back(1);
        bamrec::putI32(out, (int32_t)headerText.size());
        out.insert(out.end(), headerText.begin(), headerText.end());
        bamrec::putI32(out, (int32_t)refs.size());
        for (const BamRef& r : refs) {
            bamrec::putI32(out, (int32_t)(r.name.size() + 1));
            out.insert(out.end(), r.name.begin(), r.name.end());
            out.push_back(0);
            bamrec::putI32(out, r.len);
        }

        dataStart = pos;
        return 0;
    }

}  // namespace

namespace {

/* Uncompressed payload of one BGZF block; the format caps a block at 64KB. */
const size_t kBgzfBlockBytes = 65536;

/*
 * Deflate one chunk into a complete BGZF block (gzip header with "BC" extra field, deflate payload, crc32/isize trailer).
 * Every call owns its z_stream, so a batch of blocks can be deflated concurrently.
 */
int32_t bgzfCompressBlock(const uint8_t* src, size_t srcLen, std::vector<uint8_t>& out)
{
    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    /* Raw deflate stream: the gzip header/trailer is assembled by the BGZF block itself */
    if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        LOG_ERROR("BamWriter: deflateInit2 failed.");
        return -1;
    }

    out.resize(deflateBound(&zs, (uLong)srcLen) + 26);   /* + 18 header + 8 trailer */
    zs.next_in = (Bytef*)src;
    zs.avail_in = (uInt)srcLen;
    zs.next_out = out.data() + 18;   /* leave room for the 18-byte gzip header */
    zs.avail_out = (uInt)(out.size() - 18);
    const int rc = deflate(&zs, Z_FINISH);
    if (rc != Z_STREAM_END) {
        LOG_ERROR("BamWriter: deflate block failed, rc=%d.", rc);
        deflateEnd(&zs);
        return -1;
    }
    const size_t compLen = (out.size() - 18) - zs.avail_out;
    deflateEnd(&zs);

    const uint32_t crc = crc32(0L, src, (uInt)srcLen);
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
    /* Shrink to [header][deflate payload][trailer]. */
    out.resize(18 + compLen + 8);
    memcpy(out.data(), header, sizeof(header));
    const uint32_t isize = (uint32_t)srcLen;
    memcpy(out.data() + 18 + compLen, &crc, 4);
    memcpy(out.data() + 18 + compLen + 4, &isize, 4);
    return 0;
}

/*
 * Run fn(i) for every i in [0, n) on the calling thread plus up to maxThreads-1 helpers; fn must only touch storage it owns.
 * Order is preserved (each task owns its slot), so the output is identical to a serial run.
 */
template <typename Fn>
void parallelFor(size_t n, uint32_t maxThreads, const Fn& fn)
{
    if (n == 0) {
        return;
    }
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) {
        hw = 2;
    }
    if (maxThreads == 0 || maxThreads > hw) {
        maxThreads = hw;
    }
    if (maxThreads == 0) {
        maxThreads = 1;
    }
    size_t helpers = (size_t)maxThreads - 1;
    if (helpers > n - 1) {
        helpers = n - 1;
    }
    if (helpers == 0) {
        for (size_t i = 0; i < n; ++i) {
            fn(i);
        }
        return;
    }
    std::atomic<size_t> next(0);
    auto work = [&]() {
        size_t i;
        while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n) {
            fn(i);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(helpers);
    for (size_t k = 0; k < helpers; ++k) {
        pool.emplace_back(work);
    }
    work();
    for (size_t k = 0; k < pool.size(); ++k) {
        pool[k].join();
    }
}

/*
 * Compress `total` bytes into consecutive BGZF blocks in parallel, preserving order; the byte stream is identical to serial.
 */
int32_t bgzfCompressRange(const uint8_t* src, size_t total, uint32_t maxThreads,
                          std::vector<std::vector<uint8_t>>& out)
{
    const size_t nBlocks = (total + kBgzfBlockBytes - 1) / kBgzfBlockBytes;
    if (nBlocks == 0) {
        return 0;
    }
    out.resize(nBlocks);

    std::atomic<bool> failed(false);
    parallelFor(nBlocks, maxThreads, [&](size_t i) {
        const size_t off = i * kBgzfBlockBytes;
        size_t len = kBgzfBlockBytes;
        if (off + len > total) {
            len = total - off;
        }
        if (0 != bgzfCompressBlock(src + off, len, out[i])) {
            failed.store(true, std::memory_order_relaxed);
        }
    });
    return failed.load() ? -1 : 0;
}

/*
 * Convert a span of SAM text lines into the raw BAM record stream, preserving line order.
 * Chunks are aligned to line boundaries and converted in parallel; `out` is appended to.
 */
int32_t convertSamDataToBam(const uint8_t* buffer, size_t start, size_t end,
                            const std::map<std::string, int32_t>& refIndex,
                            uint32_t maxThreads, std::vector<uint8_t>& out)
{
    if (start >= end) {
        return 0;
    }
    const size_t len = end - start;

    uint32_t budget = maxThreads;
    if (budget == 0) {
        budget = std::thread::hardware_concurrency();
        if (budget == 0) {
            budget = 2;
        }
    }
    size_t nChunks = (len + (512 * 1024) - 1) / (512 * 1024);
    if (nChunks > (size_t)budget) {
        nChunks = (size_t)budget;
    }
    if (nChunks < 1) {
        nChunks = 1;
    }

    /* Line-aligned chunk boundaries. */
    std::vector<size_t> bounds;
    bounds.reserve(nChunks + 1);
    bounds.push_back(start);
    for (size_t c = 1; c < nChunks; ++c) {
        size_t p = start + (len * c) / nChunks;
        while (p < end && buffer[p] != '\n') {
            ++p;
        }
        if (p < end) {
            ++p;   /* the next chunk starts just after the newline */
        }
        if (p > bounds.back()) {
            bounds.push_back(p);
        }
    }
    bounds.push_back(end);

    const size_t nParts = bounds.size() - 1;
    std::vector<std::vector<uint8_t>> parts(nParts);
    std::atomic<bool> failed(false);
    parallelFor(nParts, maxThreads, [&](size_t i) {
        bamrec::BamRecordScratch scratch;
        const uint8_t* p = buffer + bounds[i];
        const uint8_t* e = buffer + bounds[i + 1];
        std::vector<uint8_t>& part = parts[i];
        while (p < e) {
            const uint8_t* nl = (const uint8_t*)memchr(p, '\n', (size_t)(e - p));
            const size_t lineLen = (nl != nullptr) ? (size_t)(nl - p) : (size_t)(e - p);
            if (lineLen > 0 && *p != '@') {
                const int32_t built = bamrec::buildBamRecordFromSamLine(p, lineLen, refIndex, scratch);
                if (built < 0) {
                    failed.store(true, std::memory_order_relaxed);
                    return;
                }
                if (built > 0) {
                    part.insert(part.end(), scratch.rec.begin(), scratch.rec.end());
                }
            }
            p += lineLen + 1;
        }
    });
    if (failed.load()) {
        return -1;
    }

    size_t total = 0;
    for (size_t i = 0; i < parts.size(); ++i) {
        total += parts[i].size();
    }
    if (total == 0) {
        return 0;
    }
    out.reserve(out.size() + total);
    for (size_t i = 0; i < parts.size(); ++i) {
        out.insert(out.end(), parts[i].begin(), parts[i].end());
    }
    return 0;
}

}  // namespace

int32_t SamToBamConverter::convert(RoughIOBlock* blockPtr)
{
    if (blockPtr == nullptr) {
        return -1;
    }
    if (blockPtr->getDataLen() <= 0 || blockPtr->isBamPrebuilt()) {
        return 0;
    }

    const BlockType type = blockPtr->getBlockType();
    const bool isAlignment =
        BlockUtil::isSAMBlock(type) || BlockUtil::isBAMBlock(type);
    if (!isAlignment) {
        /* Non-alignment input (-b on FASTQ/binary): let BamWriter pass it through. */
        if (!headerWritten) {
            passThrough = true;
        }
        return 0;
    }
    if (passThrough) {
        return 0;
    }

    const uint8_t* data = blockPtr->getBuffer();
    const size_t dataLen = (size_t)blockPtr->getDataLen();
    size_t dataStart = 0;
    std::vector<uint8_t> bam;

    if (!headerWritten) {
        if (0 != buildBamHeaderPayload(data, dataLen, dataStart, refIndex, bam)) {
            return -1;
        }
        headerWritten = true;
    }
    if (dataStart < dataLen &&
        0 != convertSamDataToBam(data, dataStart, dataLen, refIndex, maxThreads, bam)) {
        return -1;
    }

    if (bam.empty()) {
        blockPtr->setDataLen(0);
        return 0;
    }
    if (0 != blockPtr->ensureCapacity(bam.size())) {
        LOG_ERROR("SamToBamConverter: output buffer grow failed");
        return -1;
    }
    memcpy(blockPtr->getBuffer(), bam.data(), bam.size());
    blockPtr->setDataLen((int64_t)bam.size());
    blockPtr->setBamPrebuilt(true);
    return 0;
}

BamWriter::BamWriter(IOWriter* pIoWriter, uint32_t maxThreads)
    : BlockWriter(pIoWriter), bgzfEncoderReady(false), maxThreads(maxThreads) {
    headerWritten = false;
    passThrough = false;
    finished = false;
    /* Mirror the read-side BgzfParallelDecoder: deflate runs on its own worker pool. */
    bgzfEncoderReady = bgzfEncoder.start((int)maxThreads, bgzfCompressBlock,
                                         this, &BamWriter::bgzfSink);
    if (!bgzfEncoderReady) {
        LOG_ERROR("BamWriter: failed to start BGZF deflate pool, using serial fallback.");
    }
}

BamWriter::~BamWriter() {
    if (!finished && !passThrough) {
        (void)finish();   // Fallback: flush the trailing block and the EOF marker
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

int32_t BamWriter::bgzfSink(void* ctx, const uint8_t* data, size_t len) {
    BamWriter* writer = static_cast<BamWriter*>(ctx);
    return (writer != nullptr) ? writer->writeRaw(data, len) : -1;
}

int32_t BamWriter::bgzfWrite(const void* data, size_t len) {
    if (len == 0) {
        return 0;
    }
    /* Normal path: hand the bytes to the independent deflate pool; it owns the 64KB
       membering and writes completed BGZF blocks back through bgzfSink. */
    if (bgzfEncoderReady) {
        return bgzfEncoder.write(data, len);
    }
    /* Serial fallback when the pool could not be started. */
    /* Accumulate, then frame every complete 64KB block into one parallel batch; the remainder stays buffered. */
    const size_t base = bgzfPending.size();
    bgzfPending.resize(base + len);
    memcpy(bgzfPending.data() + base, data, len);

    const size_t consumed = (bgzfPending.size() / kBgzfBlockBytes) * kBgzfBlockBytes;
    if (consumed == 0) {
        return 0;
    }
    std::vector<std::vector<uint8_t>> blocks;
    if (0 != bgzfCompressRange(bgzfPending.data(), consumed, maxThreads, blocks)) {
        return -1;
    }
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (0 != writeRaw(blocks[i].data(), blocks[i].size())) {
            return -1;
        }
    }
    bgzfPending.erase(bgzfPending.begin(), bgzfPending.begin() + (std::ptrdiff_t)consumed);
    return 0;
}

int32_t BamWriter::bgzfFlushRemaining() {
    if (bgzfEncoderReady) {
        return bgzfEncoder.finish();
    }
    if (bgzfPending.empty()) {
        return 0;
    }
    /* The final block may be shorter than 64KB, allowed only at end of stream. */
    std::vector<std::vector<uint8_t>> blocks;
    if (0 != bgzfCompressRange(bgzfPending.data(), bgzfPending.size(), maxThreads, blocks)) {
        return -1;
    }
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (0 != writeRaw(blocks[i].data(), blocks[i].size())) {
            return -1;
        }
    }
    bgzfPending.clear();
    return 0;
}

int32_t BamWriter::finish() {
    if (finished) {
        return 0;
    }
    finished = true;
    if (passThrough) {
        return 0;   // Pass-through mode has no BGZF state
    }
    if (0 != bgzfFlushRemaining()) {
        return -1;
    }
    /* BGZF EOF marker (empty block), the terminating block of a standard BAM */
    static const uint8_t kEofBlock[28] = {
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff,
        0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
        0x03, 0x00, 0, 0, 0, 0, 0, 0, 0, 0
    };
    return writeRaw(kEofBlock, sizeof(kEofBlock));
}

int32_t BamWriter::writeBamHeader(const uint8_t* data, size_t len, size_t& dataStart) {
    std::vector<uint8_t> hdr;
    if (0 != buildBamHeaderPayload(data, len, dataStart, refIndex, hdr)) {
        return -1;
    }
    return bgzfWrite(hdr.data(), hdr.size());
}

int32_t BamWriter::writeDataLines(const uint8_t* buffer, size_t start, size_t end) {
    std::vector<uint8_t> joined;
    if (0 != convertSamDataToBam(buffer, start, end, refIndex, maxThreads, joined)) {
        return -1;
    }
    if (joined.empty()) {
        return 0;
    }
    return bgzfWrite(joined.data(), joined.size());
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

    /*
     * Prebuilt block: the conversion stage (or a structured-BAM worker) already produced
     * BAM bytes, so only frame them into BGZF here. The first prebuilt block carries the
     * BAM header built by SamToBamConverter.
     */
    if (!passThrough && blockPtr->isBamPrebuilt()) {
        headerWritten = true;
        return bgzfWrite(buffer, (size_t)dataLen);
    }

    if (!passThrough && !headerWritten) {
        const BlockType type = blockPtr->getBlockType();
        if (BlockUtil::isSAMBlock(type) || BlockUtil::isBAMBlock(type)) {
            size_t dataStart = (size_t)dataLen;
            if (0 != writeBamHeader(buffer, (size_t)dataLen, dataStart)) {
                return -1;
            }
            headerWritten = true;
            /* Data lines may follow the header lines (blocks whose header was not split into its own block) */
            if (dataStart < (size_t)dataLen) {
                return writeDataLines(buffer, dataStart, (size_t)dataLen);
            }
            return 0;
        }
        LOG_WARNING("BamWriter: block type %d is not SAM/BAM, pass through raw output.", type);
        passThrough = true;
        /* Fall through to the pass-through branch below to write the current block as-is */
    }

    if (passThrough) {
        return writeRaw(buffer, (size_t)dataLen);
    }

    return writeDataLines(buffer, 0, (size_t)dataLen);
}



