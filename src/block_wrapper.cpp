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
#include "utils/md5_util.h"
#include "pbgz_manager.h"

#if defined(__x86_64__)
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

BlockType BlockReader::constructBlock(RoughIOBlock* blockPtr) {
    /// Get file type by analyzing block content
    uint8_t* buffer = const_cast<uint8_t*>(blockPtr->getBuffer());
    // if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 2 && buffer[0] == 0x1f && buffer[1] == 0x8b) {
    //     return GZIP;
    // }

    // if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 4  && buffer[0] == 'B' && buffer[1] == 'A' 
    //     && buffer[2] == 'M' && buffer[3] == '\1') {
    //     return BAM;
    // }

    // Check if it's a SAM file by detecting headers starting with @HD, @SQ, @RG, @PG, @CO
    bool isNeedCheckSam = false;
    // First block checks file header, subsequent blocks are processed based on first block's result
    if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 4 && buffer[0] == '@') {
        std::string head((char*)buffer, 3);
        if (head == "@HD" || head == "@SQ" || head == "@RG" || head == "@PG" || head == "@CO") {
            isNeedCheckSam = true;
        }
    }
    if (isNeedCheckSam || blockPtr->getBlockType() == SAM) {
        // Use memchr for optimized newline finding - highly optimized SIMD implementation
        const int64_t dataLen = blockPtr->getDataLen();
        const char* bufPtr = reinterpret_cast<const char*>(buffer);
        const char* searchStart = bufPtr;
        
        // memchr is highly optimized with SIMD assembly in standard libraries
        // Iterative search using memchr to find all newline positions
        while (true) {
            const char* newlinePtr = static_cast<const char*>(memchr(searchStart, '\n', 
                                                      dataLen - (searchStart - bufPtr)));
            if (newlinePtr == nullptr) {
                break;
            }
            blockPtr->getNpos().push_back(static_cast<uint32_t>(newlinePtr - bufPtr));
            searchStart = newlinePtr + 1;
        }
        return SAM;
    }

    if (buffer[0] == '@') {   // FASTQ format
        // Determine if it's GEN2 or GEN3
        int64_t lineNum = 0;   // Line count
        int32_t baseLen = 0;   // Base length
        int32_t maxBaseLen = 0; // Maximum base length
        int64_t lastEndlinePos = 0;  // Previous newline position
        int64_t endlinePos = 0;  // Current newline position
        
        // Build lookup table for O(1) base character validation
        // Much faster than std::string::find() which does linear search
        static bool isValidBase[256];
        static bool lookupTableInitialized = false;
        
        // Initialize lookup table only once (thread-safe by C++11 static initialization)
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
        
        int64_t linePos = 0;   // Position within current line
        
        const int64_t dataLen = blockPtr->getDataLen();
        const char* bufPtr = reinterpret_cast<const char*>(buffer);
        
        // For FASTQ, use memchr for fast newline finding but still validate per character
        // because format correctness is critical for FASTQ
        if (dataLen >= 4096) {
            // Prefetch optimization for large FASTQ blocks
            for (int64_t i = 0; i < dataLen; i += 64) {
#if defined(__x86_64__)
                _mm_prefetch((char*)(buffer + i + 256), _MM_HINT_T0);
#elif defined(__aarch64__)
                __builtin_prefetch(buffer + i + 256, 0, 3);
#endif
            }
        }
        
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
            
            // Validate the characters in the current line segment
            for (int64_t pos = (searchStart - bufPtr); pos < newlinePos; ++pos) {
                if (lineMod == 0) {
                    if (linePos == 0 && buffer[pos] != '@') {
                        return BINARY;
                    }
                } else if (lineMod == 1) {
                    // O(1) lookup table replace O(n) string search
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
        if (maxBaseLen > GENE2_MAX_BASE) {
            return FASTQ_GEN3;
        } else {
            return FASTQ_GEN2;
        }
    } else {
        return BINARY;
    }
    return BINARY;
}


int64_t BlockReader::readBlock(RoughIOBlock* blockPtr, BlockType fileType) {
    if (ioReader == nullptr || blockPtr == nullptr) {
        LOG_ERROR("IO reader or block pointer is null.");
        return -1;
    }

    uint8_t* buffer = blockPtr->getBuffer();
    size_t blockSize = blockPtr->getBlockSize();
    size_t totalLen = 0;

    blockPtr->setBlockId(blockId++);

    // First copy data from cache to buffer
    if (cacheLen > 0) {
        size_t toCopy = (cacheLen < blockSize) ? cacheLen : blockSize;
        memcpy(buffer, cache, toCopy);
        totalLen += toCopy;
        cacheLen -= toCopy;
        if (cacheLen > 0) {
            memmove(cache, cache + toCopy, cacheLen);
        }
    }

    while (totalLen < blockSize) {
        size_t readLen = ioReader->readIO(buffer + totalLen, blockSize - totalLen);
        if (readLen == 0) { 
            break; // EOF
        }
        totalLen += readLen;
    }

    if (totalLen == 0) {
        return 0; // EOF
    }

    blockPtr->setDataLen(static_cast<int64_t>(totalLen));

    // First block or already determined as FASTQ format based on first block
    if (fileType == TYPE_UNKNOW || BlockUtil::isFastqBlock(fileType)) {
        BlockType type = constructBlock(blockPtr);
        if (blockPtr->getBlockId() == 0) {
            blockPtr->setBlockType(type);
        } else {
            if (type == BINARY && (fileType == BAM || fileType == GZIP || fileType == SAM)) {
                blockPtr->setBlockType(fileType);
            } else {
                blockPtr->setBlockType(type);
            }
        }
    } else {
        blockPtr->setBlockType(fileType);
        if (BlockUtil::isSAMBlock(fileType)) {
            (void)constructBlock(blockPtr);
        }
    }

    if (BlockUtil::isFastqBlock(blockPtr->getBlockType()) || BlockUtil::isSAMBlock(blockPtr->getBlockType())) {
        // If FASTQ format, ensure block integrity
        int32_t lineNum = static_cast<int32_t>(blockPtr->getNpos().size());
        int64_t remainLen = 0;
        if (totalLen >= blockSize) {
            if (BlockUtil::isFastqBlock(blockPtr->getBlockType())) { 
                remainLen = totalLen - blockPtr->getNpos()[((lineNum >> 2) << 2) - 1] - 1;      
            } else if (BlockUtil::isSAMBlock(blockPtr->getBlockType())) {
                remainLen = totalLen - blockPtr->getNpos()[lineNum - 1] - 1;
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
            totalLen -= remainLen;
        } else {
            // if not end with \n, add one
            uint32_t lastNPos  = blockPtr->getNpos()[lineNum - 1];
            if (lastNPos < totalLen - 1) {
                blockPtr->getBuffer()[totalLen] = '\n';
                blockPtr->getNpos().push_back(totalLen);
                totalLen += 1;
                lineNum += 1;
            }
        }

        blockPtr->setDataLen(static_cast<int64_t>(totalLen));
        if (BlockUtil::isFastqBlock(blockPtr->getBlockType())) {
            for (int t = 0; t < lineNum - ((lineNum >> 2) << 2); ++t)  {
                blockPtr->getNpos().pop_back(); // Remove last newline position
            }
        }
    }

    LOG_DEBUG("Read One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d,lineNum=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen(),blockPtr->getNpos().size());
    
    return totalLen;
}

namespace BlockUtil {
    bool isFastqBlock(BlockType type) {
        return (type == FASTQ_GEN2 || type == FASTQ_GEN3 || type == FASTQ_GEN2_GZIP || type == FASTQ_GEN3_GZIP);
    }

    bool isSAMBlock(BlockType type) {
        return (type == SAM);
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
                return "FASTQ_GEN2_GZIP";
            case FASTQ_GEN3_GZIP:
                return "FASTQ_GEN3_GZIP";
            case BAM:
                return "BAM";
            case SAM:
                return "SAM";
            case PBGZFILE:
                return "PBGZFILE";
            default:
                return "UNKNOWN_TYPE";
        }
    }
}

int64_t PbgzBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType __attribute__ ((unused)) fileType) {
    if (pbgzFileReader == nullptr || blockPtr == nullptr) {
        return -1;
    }

    PbgzDataBlock pbgzDataBlock;
    pbgzDataBlock.setDataPtr(blockPtr->getBuffer());
    if (0 != pbgzFileReader->readDataBlock(pbgzDataBlock)) {
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
    if (blockType == "fastq_gen2") {
        blockPtr->setBlockType(FASTQ_GEN2);
    } else if (blockType == "fastq_gen3") {
        blockPtr->setBlockType(FASTQ_GEN3);
    } else if (blockType == "binary") {
        blockPtr->setBlockType(BINARY);
    } else if (blockType == "refe_gene") {
        blockPtr->setBlockType(REFERENCE);
    } else if (blockType == "refe_gene_index") {
        blockPtr->setBlockType(REFERENCE_INDEX);
    } else if (blockType == "sam") {
        blockPtr->setBlockType(SAM);
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
    dataBlock.setMetaData("datalen", blockPtr->getDataLen());
    dataBlock.setMetaData("metalen", blockPtr->getMetaLen());
    dataBlock.setMetaData("blockid", blockPtr->getBlockId());
    if (blockPtr->getBlockType() == FASTQ_GEN2) {
        dataBlock.setMetaData("blocktype", "fastq_gen2");
    } else if (blockPtr->getBlockType() == FASTQ_GEN3) {
        dataBlock.setMetaData("blocktype", "fastq_gen3");
    } else if (blockPtr->getBlockType() == BINARY) {
        dataBlock.setMetaData("blocktype", "binary");
    } else if (blockPtr->getBlockType() == REFERENCE) {
        dataBlock.setMetaData("blocktype", "refe_gene");
    } else if (blockPtr->getBlockType() == REFERENCE_INDEX) {
        dataBlock.setMetaData("blocktype", "refe_gene_index");
    } else if (blockPtr->getBlockType() == SAM) {
        dataBlock.setMetaData("blocktype", "sam");
    } else if (blockPtr->getBlockType() == QUAL_PRIOR) {
        dataBlock.setMetaData("blocktype", "qual_prior");
    }

    /// Calculate checksum of Meta and data
    dataBlock.calcChecksum();
    return pbgzFileWriter->writeBlockData(dataBlock);
}



