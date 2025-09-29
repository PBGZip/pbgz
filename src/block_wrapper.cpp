#include "block_wrapper.h"
#include "log/logger.h"
#include "pbgz_types.h"
#include <cstring>
#include "utils/md5util.h"

BlockType BlockReader::constructBlock(RoughIOBlock* blockPtr) {
    /// 通过分析block内容获取文件类型
    uint8_t* buffer = const_cast<uint8_t*>(blockPtr->getBuffer());
    if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 2 && buffer[0] == 0x1f && buffer[1] == 0x8b) {
        return GZIP;
    }

    if (blockPtr->getBlockId() == 0 && blockPtr->getDataLen() >= 4  && buffer[0] == 'B' && buffer[1] == 'A' 
        && buffer[2] == 'M' && buffer[3] == '\1') {
        return BAM;
    }

    if (buffer[0] == '@') {   // FASTQ格式
        // 判断是GEN2还是GEN3
        int64_t lineNum = 0;   // 行数
        int32_t baseLen = 0;   // 碱基长度
        int32_t maxBaseLen = 0; // 最大碱基长度
        int64_t lastEndlinePos = 0;  // 上一个换行符位置
        int64_t endlinePos = 0;  // 当前换行符位置
        const std::string validBase = "ACGTNactgn";
        int64_t linePos = 0;   // 当前行内的位置
        for (int64_t pos = 0; pos < blockPtr->getDataLen(); ++pos) {
            int lineMod = lineNum % 4;  
            if (buffer[pos] == '\n') {
                lastEndlinePos = endlinePos;
                endlinePos = pos;
                // 碱基行的结束符
                if (lineMod == 1) {
                    baseLen = endlinePos - lastEndlinePos - 1;
                    // 更新最大碱基数
                    if (baseLen > maxBaseLen) {
                        maxBaseLen = baseLen;
                    }
                } else if (lineMod == 3) {  
                    // 质量值行和碱基行长度必须一致
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
                if (linePos == 0 && buffer[pos] != '@') { // 行ID
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
    size_t bufferSize = blockPtr->getBufferSize();
    size_t totalLen = 0;

    blockPtr->setBlockId(blockId++);

    // 先将缓存中的数据拷贝到buffer中
    if (cacheLen > 0) {
        size_t toCopy = (cacheLen < bufferSize) ? cacheLen : bufferSize;
        memcpy(buffer, cache, toCopy);
        totalLen += toCopy;
        cacheLen -= toCopy;
        if (cacheLen > 0) {
            memmove(cache, cache + toCopy, cacheLen);
        }
    }

    while (totalLen < bufferSize) {
        size_t readLen = ioReader->readIO(buffer + totalLen, bufferSize - totalLen);
        if (readLen == 0) { 
            break; // EOF
        }
        totalLen += readLen;
    }

    if (totalLen == 0) {
        return 0; // EOF
    }

    blockPtr->setDataLen(static_cast<int64_t>(totalLen));

    // 首块或者根据首块已经确定为FASTQ格式
    if (fileType == TYPE_UNKNOW || BlockUtil::isFastqBlock(fileType)) {
        BlockType type = constructBlock(blockPtr);
        if (blockPtr->getBlockId() == 0) {
            blockPtr->setBlockType(type);
        } else {
            if (type == BINARY && (fileType == BAM || fileType == GZIP)) {
                blockPtr->setBlockType(fileType);
            }
            else {
                blockPtr->setBlockType(type);
            }
        }
    } else {
        blockPtr->setBlockType(fileType);
    }

    if (BlockUtil::isFastqBlock(blockPtr->getBlockType())) {
        // 如果是FASTQ格式，确保块的完整性
        int32_t lineNum = static_cast<int32_t>(blockPtr->getNpos().size());
        int64_t remainLen = totalLen - blockPtr->getNpos()[((lineNum >> 2) << 2) - 1] - 1;
        if (remainLen > 0) {
            // 如果缓存中还有数据，先memmove
            if (cacheLen > 0) {
                memmove(cache + remainLen, cache, cacheLen);
            }
            // 将剩余的数据放到缓存中
            memcpy(cache, buffer + totalLen - remainLen, remainLen);
            cacheLen = cacheLen + remainLen;
            totalLen -= remainLen;
        }

        blockPtr->setDataLen(static_cast<int64_t>(totalLen));
        for (int t = 0; t < lineNum - ((lineNum >> 2) << 2); ++t)  {
            blockPtr->getNpos().pop_back(); // 去掉最后一个换行符位置
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

int64_t PbgzBlockReader::readBlock(RoughIOBlock* blockPtr, BlockType fileType) {
    if (pbgzFileReader == nullptr || blockPtr == nullptr) {
        return -1;
    }

    PbgzDataBlock pbgzDataBlock;
    if (0 != pbgzFileReader->readDataBlock(pbgzDataBlock)) {
        LOG_ERROR("Read Pbgz data block faild.");
        return -1;
    }

    // 
    if (0 != pbgzDataBlock.verifyCheckSum()) {
        LOG_ERROR("Verify pbgz block checksum failed");
        return -1;
    }
    
    blockPtr->setBlockId(blockId++);
    blockPtr->setDataLen(pbgzDataBlock.getMetaData("datalen").asInt64());
    blockPtr->setMetaLen(pbgzDataBlock.getMetaData("metalen").asInt64());
    std::string blockType = pbgzDataBlock.getMetaData("blocktype").asString();
    if (blockType == "fastq_gen2") {
        blockPtr->setBlockType(FASTQ_GEN2);
    } else if (blockType == "fastq_gen3") {
        blockPtr->setBlockType(FASTQ_GEN3);
    } else if (blockType == "binary") {
        blockPtr->setBlockType(BINARY);
    } 
    // 前面放的是data, 后面放的是meta
    memcpy(blockPtr->getBuffer(), pbgzDataBlock.getDataPtr(), pbgzDataBlock.getDataLength());
    
    return pbgzDataBlock.getDataLength();
}   

int32_t PbgzBlockReader::init() {
    if (pbgzFileReader != nullptr) {
        return 0;
    }
    if (ioReader == nullptr) {
        return -1;
    }

    pbgzFileReader = new PbgzFileReader(ioReader);
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
    if (pbgzFileWriter != nullptr) {
        return 0;
    }

    if (ioWriter == nullptr) {
        return -1;
    }

    pbgzFileWriter = new PbgzFileWriter(ioWriter);
    if (pbgzFileWriter == nullptr) {
        return -1;
    }

    pbgzFileWriter->open();

    PbgzFileMeta pgzeFileMeta;
    pgzeFileMeta.setMetaData("writer", "pbgz_writer_v2.0.0");
    pgzeFileMeta.setMetaData("hashmethod", "md5");
    pgzeFileMeta.setMetaData("blocksize", BLOCK_SIZE);
    pbgzFileWriter->setFileMeta(pgzeFileMeta);
    pbgzFileWriter->writeFileMeta();

    return 0;
}

int32_t PbgzBlockWriter::writeBlock(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr) {
        return -1;
    }

    if (pbgzFileWriter == nullptr) {
        return -1;
    }

    PbgzDataBlock dataBlock;    
    dataBlock.setBlockData(blockPtr->getBuffer(), blockPtr->getTotalDataLen());
    dataBlock.setMetaData("datalen", blockPtr->getDataLen());
    dataBlock.setMetaData("metalen", blockPtr->getMetaLen());
    if (blockPtr->getBlockType() == FASTQ_GEN2) {
        dataBlock.setMetaData("blocktype", "fastq_gen2");
    } else if (blockPtr->getBlockType() == FASTQ_GEN3) {
        dataBlock.setMetaData("blocktype", "fastq_gen3");
    } else if (blockPtr->getBlockType() == BINARY) {
        dataBlock.setMetaData("blocktype", "binary");
    } 
    /// 计算Meta和数据的校验和
    dataBlock.calcChecksum();
    pbgzFileWriter->writeBlockData(dataBlock);
    return 0;
}



