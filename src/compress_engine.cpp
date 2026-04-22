/*
 * compress_engine.cpp - Source file for pbgz project
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


#include "compress_engine.h"
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"
#include "codec_actuator_adapter.h"

int32_t CompressEngine::init() {
    if (0 != CodecEngine::init()) {
        return -1;
    }

    indexBlockQueue->setCapility(2);
    freeIndexBlockQueue->setCapility(2);
    for (uint32_t i = 0; i < freeIndexBlockQueue->getCapility(); ++i) {
        RoughIOBlock* ptr = MemoryUtil::safeNewClass<RoughIOBlock>();
        if (ptr == nullptr) {
            LOG_ERROR("Failed to create RoughIOBlock for freeIndexBlockQueue");
            return -1;
        }
        freeIndexBlockQueue->push(ptr);
    }
    
    return 0;
}

CompressEngine::~CompressEngine() {
    while(!freeIndexBlockQueue->empty()) {
        RoughIOBlock* ptr = freeIndexBlockQueue->get();
        MemoryUtil::safeDeleteClass(ptr);
    }

    while(!indexBlockQueue->empty()) {
        RoughIOBlock* ptr = indexBlockQueue->get();
        MemoryUtil::safeDeleteClass(ptr);
    }
}

BlockReader* CompressEngine::createBlockReader() {
    BlockReader* blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return nullptr;
    }
    if (0 != blockReader->init()) {
        LOG_ERROR("BlockReader init failed");
        MemoryUtil::safeDeleteClass(blockReader);
        return nullptr;
    }

    return blockReader;
}

BlockWriter* CompressEngine::createBlockWriter() {
   PbgzBlockWriter* blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("pbgzWriter is NULL.");
        return nullptr;
    }
    blockWriter->init();
    blockWriter->setBaseFileMeta(baseFileMeta);
    blockWriter->writeBaseFileMeta();
    return blockWriter;
}

void CompressEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void CompressEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    MemoryUtil::safeDeleteClass(blockWriter);
}

Actuator* CompressEngine::actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    FastqCompressActuator* fastqActuator = dynamic_cast<FastqCompressActuator*>(actuator);
    if (fastqActuator != nullptr) {
        FastqCodecActuator* fastqCodecActuator = fastqActuator->getCodecActuator();
        if (0 != fastqCodecActuator->preAnalysis()) {
            LOG_INFO("Fastq preAnalysis failed");
            MemoryUtil::safeDeleteClass(fastqActuator);
            return MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr);
        }
    }

    SamCompressActuator* samActuator = dynamic_cast<SamCompressActuator*>(actuator);
    if (samActuator != nullptr) {
        SamCodecActuator* samCodecActuator = samActuator->getCodecActuator();
        if (0 != samCodecActuator->preAnalysis()) {
            LOG_INFO("sam preAnalysis failed");
            MemoryUtil::safeDeleteClass(samActuator);
            return MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr);
        }
    }
    
    return actuator;
}



int32_t CompressEngine::startEnginePreProc() {
    if (!initReference()) {
        LOG_ERROR("Engine init reference failed.");
        return -1;
    }
    return 0;
}

bool CompressEngine::initReference() {
    if (!parameter.referenceGenic.empty()){
        pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    }
    if (pRefGene) {
        if (!pRefGene->makeIndex()) {
            MemoryUtil::safeDeleteClass(pRefGene);
        }
    }
    return true;
}

uint32_t CompressEngine::calcPackRefeBlockSize() {
    int64_t each = (16 << 20);  
    int64_t total = pRefGene->getSquashLength();

    uint32_t result = total / each;
    int64_t remain = total % each;
    if (remain > 0) {
        result++;
    }
    return result;
}


/*  Save reference genome */
int64_t CompressEngine::packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef) {
    int64_t block = 0, offset = 0;
    std::vector<std::thread> tpools;
    
    uint32_t pcnt = parameter.threadNum;
    int64_t each = (16 << 20);
    int64_t total = pRefGene->getSquashLength();
    int64_t remain = total;
    maxBlockLen =0;
    totalEncLen = 0;
    
    BlockingQueue<RefeInfo> inputPool(pcnt);
    Reference* refe = pRefGene;
    int64_t current;
    std::mutex m;
    BlockingQueueType* freeOutput = freeOutputPool.get();
    BlockingQueueType* outputPool = outputDataPool.get();
    uint32_t count = blockCount;

    // Start threads
    for (uint32_t idx = 0; idx < pcnt; ++idx) {
        tpools.push_back(std::thread([&inputPool, &m, &freeOutput, &outputPool, &refe, &each, &maxBlockLen, &totalEncLen, &isSanitizeRef, &count]() {
            uint8_t* output = MemoryUtil::safeAlloc<uint8_t>(each);
            while(true) {
                RefeInfo refe2do = inputPool.get();
                int64_t plen = refe2do.second.second;
                if (plen == 0) {
                    break;
                }
                const uint8_t *p = refe->getSquash() + refe2do.second.first;
                if (isSanitizeRef) {
                    refe->sanitizeRefSquash(refe2do.second.first, plen);
                }

                coder_io refeIo(output, each);
                coder_ppmd cppmd(&refeIo);
                cppmd.encode(p, plen);
                cppmd.encode_flush();

                Json::Value meta;
                meta["block"] = (Json::Value::Int)refe2do.first.first;
                meta["offset"] = (Json::Value::Int64)refe2do.second.first;
                meta["srclen"] = (Json::Value::Int)plen;
                meta["dstlen"] = (Json::Value::Int)refeIo.data_len;
                meta["coder"] = refeIo.meta;
                coder_json cmeta;
                std::string metaString;
                cmeta.encoder(meta, metaString); /*  Compress block meta */
            
                /* Write compressed stream of current block's meta */
                RoughIOBlock* outBlock = freeOutput->get();
                outBlock->reset();
                 int64_t currBlockLen = metaString.length() + refeIo.data_len;
                if (currBlockLen >= outBlock->getRemain()) {
                    LOG_ERROR("reference block pack failed");
                    break;
                }
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);
                outBlock->setBlockId(count + refe2do.first.first);
                outputPool->push(outBlock);

                m.lock();
                maxBlockLen = (currBlockLen > maxBlockLen) ? currBlockLen : maxBlockLen;
                totalEncLen += currBlockLen;
                m.unlock();
            }
            MemoryUtil::safeFree(output);
        }));
    }

    while (remain > 0) {
        current = std::min(each, remain);
        inputPool.push(std::make_pair(std::make_pair(block, nullptr), std::make_pair(offset, current)));
        block++;
        offset += current;
        remain -= current;
    }
    // Write end marker
    for (uint32_t i = 0; i < pcnt; ++i) {
        std::pair<std::pair<int64_t, uint8_t *>, std::pair<int64_t, int64_t>> refe2do;
        refe2do.second.second = 0;
        inputPool.push(refe2do);
    }
    for (auto &t : tpools) {
        if (t.joinable()) {
            t.join();
        }
    }

    return block;
}

int32_t CompressEngine::startEnginePostProc() {
    bool isPackRefeInTail = pRefGene && !parameter.isUnpackRef && parameter.outputFile != STDOUT;
    if (isPackRefeInTail) {
        int64_t maxRefLen = 0;
        int64_t totolEncLen = 0;
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath());
        refeMeta["offset"]= 0;
        dynamicFileMeta.setMetaData("refe", refeMeta);
        FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
        if (fileWriter != nullptr) {
            int64_t refBlockCount = packReference(maxRefLen, totolEncLen);
            dynamicFileMeta.getMetaData("refe")["max_block_len"] = maxRefLen;
            dynamicFileMeta.getMetaData("refe")["blocks"] = refBlockCount;
        }
    }
    
    return CodecEngine::engineStartAfterProc();
}

void CompressEngine::setDataBlockPosition(uint32_t blockId) {
    if (!parameter.isMakeIndex) {
        return;
    }
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }
    pbgzIndex.setBlockPosition(blockId, fileWriter->getCurrentPos());
    return;
 }

 void CompressEngine::startWriteIndexTask() {
    if (!parameter.isMakeIndex) {
        return;
    }

    return;
 }

 int32_t CompressEngine::startWorkPreProc() {
    // Refresh file meta
    if (pRefGene) {
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath()); /* Contains md5 information, used for decompression verification */
        refeMeta["max_block_len"] = 16 << 20;
        bool isPackRefeInHeader = !parameter.isUnpackRef && parameter.outputFile == STDOUT;
        if (isPackRefeInHeader) {
            refeMeta["blocks"] = calcPackRefeBlockSize();
            baseFileMeta.setMetaData("refe",refeMeta);
            int64_t maxRefLen = 0;
            int64_t totolEncLen = 0;
            (void)packReference(maxRefLen, totolEncLen, false);
        } else {
            refeMeta["blocks"] = 0;
            baseFileMeta.setMetaData("refe",refeMeta);
        }
    }

    return 0;
 }

 Actuator* CompressEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    if (parameter.isMakeIndex) {
        fprintf(stderr, "Fastq file will not make index.");
        parameter.isMakeIndex = false;
    }

    Actuator* pActuator = nullptr;
    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<FastqCompressActuator>(inBlockPtr, outBlockPtr, pRefGene);
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamCompressActuator>(inBlockPtr, outBlockPtr, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr);
    } 
    
    if (pActuator == nullptr) {
        LOG_ERROR("Not support block type: %d, blockId=%d", inBlockPtr->getBlockType(), inBlockPtr->getBlockId());
        freeInputPool->push(inBlockPtr);
        outBlockPtr->reset();
        outBlockPtr->setBlockId(inBlockPtr->getBlockId());
        // When an error occurs, push a block with length 0 but correct ID, the write thread ignores blocks with length 0 to prevent thread waiting
        outputDataPool->push(outBlockPtr);
        return nullptr;
    }

    if (pActuator->initial() != 0) {
        LOG_ERROR("actuator initial failed");
        MemoryUtil::safeDeleteClass(pActuator);
        return nullptr;
    }  
    
    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType()) || BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
    }

    return pActuator;
}
