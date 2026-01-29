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

BlockType BlockReader::constructBlock(RoughIOBlock* blockPtr) {
    /// Get file type by analyzing block content
    uint8_t* buffer = const_cast<uint8_t*>(blockPtr->getBuffer());
    if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 2 && buffer[0] == 0x1f && buffer[1] == 0x8b) {
        return GZIP;
    }

    if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 4  && buffer[0] == 'B' && buffer[1] == 'A' 
        && buffer[2] == 'M' && buffer[3] == '\1') {
        return BAM;
    }

    if (buffer[0] == '@') {   // FASTQ format
        // Determine if it's GEN2 or GEN3
        int64_t lineNum = 0;   // Line count
        int32_t baseLen = 0;   // Base length
        int32_t maxBaseLen = 0; // Maximum base length
        int64_t lastEndlinePos = 0;  // Previous newline position
        int64_t endlinePos = 0;  // Current newline position
        const std::string validBase = "ACGTNactgn";
        int64_t linePos = 0;   // Position within current line
        for (int64_t pos = 0; pos < blockPtr->getDataLen(); ++pos) {
            int lineMod = lineNum % 4;  
            if (buffer[pos] == '\n') {
                lastEndlinePos = endlinePos;
                endlinePos = pos;
                // End of base line
                if (lineMod == 1) {
                    baseLen = endlinePos - lastEndlinePos - 1;
                    // Update maximum base count
                    if (baseLen > maxBaseLen) {
                        maxBaseLen = baseLen;
                    }
                } else if (lineMod == 3) {  
                    // Quality line and base line must have same length
                    if (endlinePos - lastEndlinePos - 1 != baseLen) {
                        return BINARY;
                    }
                }

                ++lineNum;
                blockPtr->getNpos().push_back(static_cast<uint32_t>(pos));
                linePos = 0;
                continue;
            }
            
            if (lineMod == 0) {
                if (linePos == 0 && buffer[pos] != '@') { // Line ID
                    return BINARY;
                }
            } else if (lineMod == 1) {
                if (validBase.find(buffer[pos]) == std::string::npos) {
                    return BINARY;
                }
            } else if (lineMod == 2) {
                if (linePos == 0 && buffer[pos] != '+') {
                    return BINARY;
                }
            }
            ++linePos;
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
            if (type == BINARY && (fileType == BAM || fileType == GZIP)) {
                blockPtr->setBlockType(fileType);
            } else {
                blockPtr->setBlockType(type);
            }
        }
    } else {
        blockPtr->setBlockType(fileType);
    }

    if (BlockUtil::isFastqBlock(blockPtr->getBlockType())) {
        // If FASTQ format, ensure block integrity
        int32_t lineNum = static_cast<int32_t>(blockPtr->getNpos().size());
        int64_t remainLen = totalLen - blockPtr->getNpos()[((lineNum >> 2) << 2) - 1] - 1;
        if (remainLen > 0) {
            // If there's still data in cache, do memmove first
            if (cacheLen > 0) {
                memmove(cache + remainLen, cache, cacheLen);
            }
            // Put remaining data into cache
            memcpy(cache, buffer + totalLen - remainLen, remainLen);
            cacheLen = cacheLen + remainLen;
            totalLen -= remainLen;
        }

        blockPtr->setDataLen(static_cast<int64_t>(totalLen));
        for (int t = 0; t < lineNum - ((lineNum >> 2) << 2); ++t)  {
            blockPtr->getNpos().pop_back(); // Remove last newline position
        }
    }
    
    return totalLen;
}

namespace BlockUtil {
    bool isFastqBlock(BlockType type) {
        return (type == FASTQ_GEN2 || type == FASTQ_GEN3 || type == FASTQ_GEN2_GZIP || type == FASTQ_GEN3_GZIP);
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
        LOG_ERROR("Read Pbgz data block faild.");
        return -1;
    }

    if (0 != pbgzDataBlock.verifyCheckSum()) {
        LOG_ERROR("Verify pbgz block checksum failed");
        return -1;
    }
    
    blockPtr->setDataLen(pbgzDataBlock.getMetaData("datalen").asInt64());
    blockPtr->setMetaLen(pbgzDataBlock.getMetaData("metalen").asInt64());
    blockPtr->setBlockId(pbgzDataBlock.getMetaData("blockid").asInt64());
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
    }
    // Copy entire block information
    memcpy(blockPtr->getBuffer(), pbgzDataBlock.getDataPtr(), pbgzDataBlock.getDataLength());

    LOG_DEBUG("Read One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen());
    
    return pbgzDataBlock.getDataLength();
}   

int32_t PbgzBlockReader::init() {
    if (ioReader == nullptr) {
        return -1;
    }
    if (pbgzFileReader == nullptr) {
        LOG_ERROR("Creat PbgzFileReader failed");
        return -1;
    }
    return pbgzFileReader->open();
}

int32_t BlockWriter::writeBlock(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr || ioWriter == nullptr) {
        return -1;
    }
    ioWriter->writeIO(blockPtr->getBuffer(), blockPtr->getDataLen());
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

    LOG_DEBUG("Write One Block, blockId=%d,blockType=%d,dataLen=%d,metalen=%d.", 
        blockPtr->getBlockId(), blockPtr->getBlockType(), blockPtr->getDataLen(),blockPtr->getMetaLen());

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
    } 

    /// Calculate checksum of Meta and data
    dataBlock.calcChecksum();
    return pbgzFileWriter->writeBlockData(dataBlock);
}



