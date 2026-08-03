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
#include "codec_selector.h"

int64_t CompressEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    int64_t ret = CodecEngine::readOneBlock(blockReader, fileType);
    
    // Initialize statistics when file type is determined (first block)
    if (ret > 0 && !statsInitialized && fileType != TYPE_UNKNOW) {
        initStatsBasedOnFileType(fileType);
    }
    
    return ret;
}

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

void CompressEngine::initStatsBasedOnFileType(BlockType fileType) {
    if (statsInitialized) {
        return;
    }
    
    if (BlockUtil::isSAMBlock(fileType)) {
        stats = std::make_unique<SamStat>();
        if (stats->init() != 0) {
            LOG_ERROR("Failed to initialize SamStat");
        } else {
            LOG_INFO("Initialized SamStat for SAM file compression");
        }
    } else if (BlockUtil::isFastqBlock(fileType)) {
        stats = std::make_unique<FastqStat>();
        if (stats->init() != 0) {
            LOG_ERROR("Failed to initialize FastqStat");
        } else {
            LOG_INFO("Initialized FastqStat for FASTQ file compression");
        }
    } else {
        LOG_DEBUG("File type %d does not require stats initialization", fileType);
    }
    
    statsInitialized = true;
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
            BinaryCompressActuator* binaryActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
            if (binaryActuator == nullptr) {
                return nullptr;
            }
            binaryActuator->initial();
            return binaryActuator;
        }
    }

    SamCompressActuator* samActuator = dynamic_cast<SamCompressActuator*>(actuator);
    if (samActuator != nullptr) {
        SamCodecActuator* samCodecActuator = samActuator->getCodecActuator();
        if (0 != samCodecActuator->preAnalysis()) {
            LOG_INFO("sam preAnalysis failed");
            MemoryUtil::safeDeleteClass(samActuator);
            BinaryCompressActuator* binaryActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
            if (binaryActuator == nullptr) {
                return nullptr;
            }
            binaryActuator->initial();
            return binaryActuator;
        }
    }

    runFilePreprocessOnce(inBlockPtr);

    return actuator;
}

void CompressEngine::runFilePreprocessOnce(RoughIOBlock* inBlockPtr) {
    /*
     * Lock-free, non-blocking claim: the first worker to arrive (state==IDLE)
     * atomically transitions to RUNNING and performs the trial-compression.
     * All other workers see RUNNING/DONE and proceed immediately with their
     * default coder -- the compression pipeline never stalls on preprocessing.
     *
     * Preprocessing's natural concurrency is "test different coders on one
     * sample" (independent, parallelizable), which is logically distinct from
     * the pipeline's "different blocks on different workers" concurrency.
     */
    if (!preprocessInfo.tryClaim()) {
        return;
    }

    if (0 != CodecSelector::analyze(inBlockPtr, preprocessInfo)) {
        LOG_INFO("File preprocessing (codec selection) failed, using default coders");
    } else if (parameter.verbose) {
        static const char* samFieldNames[] = {
            "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
            "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
        };
        static const char* fqFieldNames[] = {
            "ID", "SEQ", "QUAL", "COMMENT"
        };

        bool isSam = (preprocessInfo.fileType == SAM);
        uint32_t fieldCount = isSam ? SAM_FIELD_COUNT : FQ_FIELD_COUNT;
        const char* const* names = isSam ? samFieldNames : fqFieldNames;

        fprintf(stderr, "\n[preprocess] codec selection"
                " (scanned %u raw bytes, field samples total %u):\n",
                preprocessInfo.scannedBytes, preprocessInfo.sampleBytes);
        for (uint32_t i = 0; i < fieldCount && i < preprocessInfo.fields.size(); ++i) {
            const auto& sel = preprocessInfo.fields[i];
            const char* status = (sel.status == FieldStatus::SELECTED) ? "selected" :
                                 (sel.status == FieldStatus::SKIPPED) ? "skipped" : "failed";
            if (sel.status == FieldStatus::SELECTED) {
                /* 吞吐按 样本字节 / 试压耗时 折算成 MB/s；耗时为 0 时不显示速度 */
                double bwtMBps = sel.trialBwtCmUs ? (double)sel.decidedLen / sel.trialBwtCmUs : 0.0;
                double fcMBps  = sel.trialFcUs    ? (double)sel.decidedLen / sel.trialFcUs    : 0.0;
                /*
                 * 质量值列走的是专用评估路径，两个试压字段承载的是 coder_qual 和
                 * fcv2 的结果，而不是通用路径的 bwt_cm 和 fc。标签必须跟着换，
                 * 否则读的人会把 coder_qual 的数字当成 coder_bwt_cm 的。
                 */
                const bool isQualRow = isSam && (i == (uint32_t)SAM_QUAL);
                const char* candA = isQualRow ? "coder_qual" : "bwt_cm";
                const char* candB = isQualRow ? "fcv2"       : "fc";
                fprintf(stderr, "  %-6s -> %-14s  %u -> %u (%.2f%%)"
                        "  [%s %.2f%% @%.0fMB/s | %s %.2f%% @%.0fMB/s] %u轮\n",
                        names[i], coderTypeToMagic(sel.selectedCoder),
                        sel.decidedLen, sel.bestCompLen,
                        sel.decidedLen ? 100.0 * sel.bestCompLen / sel.decidedLen : 0.0,
                        candA, sel.decidedLen ? 100.0*sel.trialBwtCmLen/sel.decidedLen : 0, bwtMBps,
                        candB, sel.decidedLen ? 100.0*sel.trialFcLen/sel.decidedLen : 0, fcMBps,
                        sel.rounds);
            } else {
                fprintf(stderr, "  %-6s -> %-14s  (sample %u bytes)\n",
                        names[i], status, sel.sampleLen);
            }
        }
        fprintf(stderr, "\n");
    }
    preprocessInfo.markDone();;
}

void CompressEngine::writeFilePostProc(BlockWriter* blockWriter) {
    CodecEngine::writeFilePostProc(blockWriter);
    // make index
    if (parameter.isMakeIndex) {
        std::string indexFileName = parameter.outputFile + ".pbgzi";
        if (PathUtil::fileExists(indexFileName)) {
            PathUtil::removeFile(indexFileName);
        }
        SamIndex::getInstance().updateFileOffsetsFromBlockPosition();
        SamIndex::getInstance().dumpToFile(indexFileName);
    }
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
            LOG_ERROR("makeIndex failed, reference will not be used");
            MemoryUtil::safeDeleteClass(pRefGene);
            pRefGene = nullptr;
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

                RoughIOBlock* outBlock = freeOutput->get();
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
        refeMeta["ni_name"] = "";
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
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }
    BlockPosition::getInstance().setBlockPosition(blockId, fileWriter->getCurrentPos());
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
        refeMeta["ni_name"] = "";
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
    Actuator* pActuator = nullptr;
    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType())) {
        if (parameter.isMakeIndex) {
            fprintf(stderr, "Fastq file will not make index.");
            parameter.isMakeIndex = false;
        }
        pActuator = MemoryUtil::safeNewClass<FastqCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
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
