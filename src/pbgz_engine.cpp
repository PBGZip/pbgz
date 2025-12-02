/*
 * pbgz_engine.cpp - Cpp file for pbgz project
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

#include "pbgz_engine.h"
#include "log/logger.h"
#include "coder/coder.h"
#include "utils/memory_util.h"
#include "pbgz_manager.h"
#include "block_wrapper.h"
#include "fastq_actuator.h"
#include "binary_actuator.h"
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"
#include "config_manager.h"
#ifdef __SSE4_2__ 
#include "hardware.h"
#endif
#include <set>

int PbgzEngine::init() {
    // Register allocation functions required by coder
    coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
    coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
    coder_ns::register_exit_proc(pbgzExitProc);
    coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    coder_ns::resister_logger_proc(coderLog);

    // Create queues
    freeInputPool.setCapility(parameter.threadNum);
    inputDataPool.setCapility(parameter.threadNum);
    freeOutputPool.setCapility(parameter.threadNum << 1);
    outputDataPool.setCapility(parameter.threadNum << 1);

    uint32_t blockBufferSize = ConfigManager::getInstance().getBlockSizeByCompressLevel(parameter.compressLevel);
    // First push empty blocks to free queue
    for (uint32_t i = 0; i < freeInputPool.getCapility(); ++i) {
        RoughIOBlock* inPtr = MemoryUtil::safeNewClass<RoughIOBlock>(blockBufferSize);
        if (inPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }
        freeInputPool.push(inPtr);
    }

    for (uint32_t j = 0; j < freeOutputPool.getCapility(); ++j) {
        RoughIOBlock* outPtr = MemoryUtil::safeNewClass<RoughIOBlock>(blockBufferSize);
        if (outPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }
        freeOutputPool.push(outPtr);
    }

    if (parameter.inputFile == "/dev/stdin") {
        ioReader = MemoryUtil::safeNewClass<PipeReader>();
        LOG_INFO("Create PipeReader.");
    } else {
        if (PathUtil::isGzFile(parameter.inputFile)) {
            bool isSupportSimd = false;
#ifdef __SSE4_2__
            isSupportSimd = Hardware().isSupportSimd();
#endif
            if (isSupportSimd) {
                ioReader = MemoryUtil::safeNewClass<FastGzFileReader>(parameter.inputFile);
                LOG_INFO("Create FastGzFileReader.");
            } else {
                ioReader = MemoryUtil::safeNewClass<GzFileReader>(parameter.inputFile, parameter.threadNum);
                LOG_INFO("Create GzFileReader.");
            }
        } else {
            ioReader = MemoryUtil::safeNewClass<FileReader>(parameter.inputFile);
            LOG_INFO("Create FileReader.");
            
        }
    }
    if (ioReader == nullptr) {
        LOG_ERROR("Create IO reader failed.");
        return -1;
    }
    ioReader->openIO();

    if (parameter.outputFile == "/dev/stdout") {
        if(parameter.isDecToGZ) {
            ioWriter = MemoryUtil::safeNewClass<GzPipeWriter>(parameter.threadNum);
            LOG_INFO("Create GzPipeWriter.");
        } else {
            ioWriter = MemoryUtil::safeNewClass<PipeWriter>();
            LOG_INFO("Create PipeWriter.");
        }
    } else {
        if(parameter.isDecToGZ) {
            ioWriter = MemoryUtil::safeNewClass<GzFileWriter>(parameter.outputFile, parameter.threadNum);
            LOG_INFO("Create GzFileWriter.");
        } else {
            ioWriter = MemoryUtil::safeNewClass<FileWriter>(parameter.outputFile);
            LOG_INFO("Create FileWriter.");
        }
    }
    if (ioWriter == nullptr) {
        LOG_ERROR("Create IO reader failed.");
        return -1;
    }
    ioWriter->openIO();

    return 0;
}

PbgzEngine::~PbgzEngine() {
    for (auto& th : coderThreads) {
        if (th.joinable()) {
            th.join();
        }
    }
    
    if (writeThread.joinable()) {
        writeThread.join();
    }

    // Release resources
    if (ioReader) {
        ioReader->closeIO();
        MemoryUtil::safeDeleteClass(ioReader);
    }

    if (ioWriter) {
        ioWriter->closeIO();
        MemoryUtil::safeDeleteClass(ioWriter);
    }

    while(!freeInputPool.empty()) {
        RoughIOBlock* inPtr = freeInputPool.get();
        MemoryUtil::safeDeleteClass(inPtr);
    }

    while(!inputDataPool.empty()) {
        RoughIOBlock* inPtr = inputDataPool.get();
        MemoryUtil::safeDeleteClass(inPtr);
    }

    while(!freeOutputPool.empty()) {
        RoughIOBlock* outPtr = freeOutputPool.get();
        MemoryUtil::safeDeleteClass(outPtr);
    }

    while(!outputDataPool.empty()) {
        RoughIOBlock* outPtr = outputDataPool.get();
        MemoryUtil::safeDeleteClass(outPtr);
    }

    MemoryUtil::safeDeleteClass(pRefGene);
}

int32_t PbgzEngine::start() {
    PbgzManager::getInstance().printHeadInfo(parameter);

    Timer costTimer(true);
    if (!parameter.isDecompress) {
        if (!initReferenceForCompress()) {
            LOG_ERROR("Init refrence for compress failed.");
            return -1;
        }
    }
    int32_t ret = startWriteTask();
    if (ret != 0) {
        LOG_ERROR("Start write task failed.");
        return -1;
    }

    ret = startCoderTask();
    if (ret != 0) {
        LOG_ERROR("Start coder task failed.");  
        return -1;
    }

    ret = startReadTask();
    if (ret != 0) {
        LOG_ERROR("Start read task failed.");
        // 编码的终止符号
        for (uint32_t i = 0; i < parameter.threadNum; ++i) {
            inputDataPool.push(nullptr);
        }
        // 输出终止符号
        outputDataPool.push(nullptr);
        return -1;
    }

    for (auto& th : coderThreads) {
        if (th.joinable()) {
            th.join();
        }
    }

    bool isPackRefeInTail = !parameter.isDecompress && pRefGene && !parameter.isUnpackRef && parameter.outputFile != "/dev/stdout";
    if (isPackRefeInTail) {
        std::vector<RoughIOBlock*> blockVec;
        int64_t maxRefLen = 0;
        int64_t totolEncLen = 0;
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath());
        FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
        if (fileWriter != nullptr) {
            int64_t refBlockCount = packReference(blockVec, maxRefLen, totolEncLen);
            refeMeta["max_block_len"] = maxRefLen;
            refeMeta["blocks"] = refBlockCount;
            // Throw to write file thread to write
            for (RoughIOBlock* item : blockVec) {
                item->setBlockId(blockCount++);
                outputDataPool.pushForce(item);
            }
            dynamicFileMeta.setMetaData("refe", refeMeta); 
        }
        blockVec.clear();
    }

    if (!parameter.isDecompress && pRefGene && parameter.outputFile != "/dev/stdout") {
        packBlockRefePosIdx();
    }

    // Write end marker
    outputDataPool.push(nullptr);
    writeThread.join();
    PbgzManager::getInstance().printTailInfo(costTimer, parameter);
    return 0;
}

int64_t PbgzEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    RoughIOBlock* blockPtr = freeInputPool.get();
    if (blockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }
    blockPtr->reset();

    int64_t ret = blockReader->readBlock(blockPtr, fileType);
    if (blockPtr->getBlockType() == REFERENCE || blockPtr->getBlockType() == REFERENCE_INDEX) {
        freeInputPool.push(blockPtr);
        return -2;
    }

    if (blockPtr->getBlockId() == 0){ 
        fileType = blockPtr->getBlockType();
        PbgzManager::getInstance().printFileType(fileType);
    }

    PbgzManager::getInstance().updateReadDataLen(blockPtr);
    inputDataPool.push(blockPtr);
    blockCount++;
    return ret;
}

int32_t PbgzEngine::startReadTask() {
    pthread_setname_np(pthread_self(), "readtask");
    BlockReader* blockReader = nullptr;
    if (parameter.isDecompress) {   // Decompression mode, read content from pbgz file
        blockReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
        PbgzBlockReader* pbgzReader = dynamic_cast<PbgzBlockReader*>(blockReader);
        if (pbgzReader == nullptr) {
            MemoryUtil::safeDeleteClass(blockReader);
            return -1;
        }
        pbgzReader->init();
        baseFileMeta = pbgzReader->getBaseFileMeta();
        dynamicFileMeta = pbgzReader->getDynamicFileMeta();
        if (!initRefGeneForDecomress(pbgzReader)) {
            LOG_INFO("Init reference for decompress failed");
            MemoryUtil::safeDeleteClass(blockReader);
            return -1;
        }

        if (!initRefeIndexForDecompress(pbgzReader)) {
            LOG_INFO("Init reference index for decompress failed");
            MemoryUtil::safeDeleteClass(blockReader);
            return -1;
        }
    } else {  // Compression mode, read content from non-pbgz file
        blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
        if (blockReader == nullptr) {
            LOG_ERROR("Create block reader failed.");
            return -1;
        }
        if (0 != blockReader->init()) {
            LOG_ERROR("BlockReader init failed");
            MemoryUtil::safeDeleteClass(blockReader);
            return -1;
        }
    }

    if (parameter.isDecompress && !parameter.refeGenePos.empty() && pRefGene) {
        readBlockByPostition(blockReader);
    } else {
        BlockType fileType = TYPE_UNKNOW;
        int64_t ret = 0;
        do {
            ret = readOneBlock(blockReader, fileType);
        } while(ret > 0 || ret == -2);
    }

    // Reached end of file or error, insert empty block as end marker to data queue
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        inputDataPool.push(nullptr);
    }

    delete blockReader;
    blockReader = nullptr;
    return 0;
}


void PbgzEngine::readBlockByPostition(BlockReader* blockReader) {
    if (blockReader == nullptr) {
        return;
    }
    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    if (fileReader == nullptr) {
        return;
    }
    if (0 != pRefGene->calcRefPosByDesc(parameter.refeGenePos, refeBeginPos, refeEndPos)) {
        return;
    }

    LOG_INFO("RefeBeginPos=%ld, RefeEndPos=%ld", refeBeginPos, refeEndPos);

    std::set<uint32_t> blockIds; 
    auto lowerIt = blockRefeIndex.lower_bound(refeBeginPos);
    auto upperIt = blockRefeIndex.upper_bound(refeEndPos);
    for (auto it = lowerIt; it != upperIt; ++it) {
        if (it->first >= refeBeginPos && it->first <= refeEndPos) {
            blockIds.insert(it->second);
        }
    }
    
    BlockType fileType = TYPE_UNKNOW;
    size_t offsetBak = fileReader->getCurrentPos();
    for (auto blockId : blockIds) {
        int64_t offset = dynamicFileMeta.getMetaData("block_pos")[blockId].asInt64();
        fileReader->seekIO(offset);
        int64_t ret = readOneBlock(blockReader, fileType);
        if (ret == 2) {
            continue;
        }
    }
    fileReader->seekIO(offsetBak);
    return;
}

int32_t PbgzEngine::startCoderTask() {
    auto coderTask = [this](int32_t id) {
        pthread_setname_np(pthread_self(), std::string("codertask_").append(std::to_string(id)).c_str());
        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool.get();
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
                break;
            }

            RoughIOBlock* outBlockPtr = freeOutputPool.get();
            if (outBlockPtr == nullptr) {
                LOG_ERROR("Get free output block failed.");
                freeInputPool.push(inBlockPtr);
                break;
            }
            outBlockPtr->reset();
            outBlockPtr->setBlockId(inBlockPtr->getBlockId());
            outBlockPtr->setBlockType(inBlockPtr->getBlockType());

            Actuator* pActuator = nullptr;
            if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType())) {
                pActuator = MemoryUtil::safeNewClass<FastqActuator>(inBlockPtr, outBlockPtr, pRefGene);
                FastqActuator* fastqActuator = dynamic_cast<FastqActuator*>(pActuator);
                if (fastqActuator != nullptr) {
                    if (parameter.isDecompress && !parameter.refeGenePos.empty() && pRefGene) {
                        fastqActuator->setRefeBeginEndPos(refeBeginPos, refeEndPos);
                    }
                    if (!parameter.isDecompress && 0 != fastqActuator->preAnalysis()) { // Only need to analyze block content in compression scenario
                        pActuator = MemoryUtil::safeNewClass<BinaryActuator>(inBlockPtr, outBlockPtr);
                    }
                }
            } else if (inBlockPtr->getBlockType() == BINARY) {
                pActuator = MemoryUtil::safeNewClass<BinaryActuator>(inBlockPtr, outBlockPtr);
            } else {
                freeInputPool.push(inBlockPtr);
                freeOutputPool.push(outBlockPtr);
                // LOG_ERROR("Not support block type: %d", inBlockPtr->getBlockType());
                continue;
            }

            if (pActuator == nullptr) {
                freeInputPool.push(inBlockPtr);
                freeOutputPool.push(outBlockPtr);
                LOG_ERROR("Create actuator failed.");
                break;
            }

            int32_t ret = 0;
            if (parameter.isDecompress) {
                ret = pActuator->decompress();
            } else {
                ret = pActuator->compress();
                FastqActuator* fastqActuator = dynamic_cast<FastqActuator*>(pActuator);
                if (fastqActuator != nullptr) {
                    blockRefePos[outBlockPtr->getBlockId()] = fastqActuator->getRefeMappedPositons();
                    LOG_INFO("BlockId:%d, MappedLength:%d", outBlockPtr->getBlockId(), fastqActuator->getRefeMappedPositons().size());
                }
            }

            if (ret != 0) {
                LOG_ERROR("Coder task failed.");
                freeInputPool.push(inBlockPtr);
                freeOutputPool.push(outBlockPtr);
                delete pActuator;
                pActuator = nullptr;
                fprintf(stderr, "Warning: block(%ld) process failed.\n", inBlockPtr->getBlockId());
                continue;
            }
            delete pActuator;
            pActuator = nullptr;
            freeInputPool.push(inBlockPtr);
            outputDataPool.push(outBlockPtr);
        }
    };
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        coderThreads.emplace_back(std::thread(coderTask, i));
    }
    return 0;
}   

int32_t PbgzEngine::startWriteTask() {
    auto writerTask = [this]() -> int32_t {
        pthread_setname_np(pthread_self(), "writetask");
        BlockWriter* blockWriter = nullptr;
        if (parameter.isDecompress) {  // Decompression mode, file write in non-pbgz format
            blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
            if (blockWriter == nullptr) {
                LOG_ERROR("Failed to create block writer.");
                return -1;
            }
            blockWriter->init();
        } else {   // Compression mode, file write in pbgz format
            blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
            PbgzBlockWriter* pbgzWriter =  dynamic_cast<PbgzBlockWriter*>(blockWriter);
            if (pbgzWriter == nullptr) {
                LOG_ERROR("pbgzWriter is NULL.");
                return -1;
            }
            pbgzWriter->init();
            pbgzWriter->setBaseFileMeta(baseFileMeta);
            pbgzWriter->writeBaseFileMeta();
        }
        
        int64_t blockId2Write = 0; 
        while (true) {
            RoughIOBlock* outBlockPtr = outputDataPool.get();
            if (outBlockPtr == nullptr) {     // Got end marker
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outblockPtr = outputSortedCache.front();
                    setDataBlockPosition(outblockPtr->getBlockId());
                    blockWriter->writeBlock(outblockPtr);
                    PbgzManager::getInstance().updateWriteDataLen(outblockPtr);
                    freeOutputPool.push(outblockPtr);
                    outputSortedCache.pop_front();
                }

                // 待所有数据写入完成之后更新扩展头，并写入动态文件meta
                PbgzBlockWriter* pbgzWriter =  dynamic_cast<PbgzBlockWriter*>(blockWriter);
                if (pbgzWriter != nullptr) {
                    pbgzWriter->updateHeadExt();
                    pbgzWriter->setDynamicFileMeta(dynamicFileMeta);
                    pbgzWriter->writeDynamicFileMeta();
                    resetReferenceOffset();
                }
                break;
            } else {
                if (outBlockPtr->getBlockType() == REFERENCE) {
                    /// Writing reference blocks is one-time, release after writing
                    FileWriter* fileWriter =  dynamic_cast<FileWriter*>(ioWriter);
                    if (fileWriter != nullptr) {
                        updateReferenceOffset(fileWriter->getCurrentPos());
                    }
                    blockWriter->writeBlock(outBlockPtr);
                    PbgzManager::getInstance().updateWriteDataLen(outBlockPtr);
                    MemoryUtil::safeDeleteClass(outBlockPtr);
                    continue;
                }

                if (outBlockPtr->getBlockType() == REFERENCE_INDEX) {
                    FileWriter* fileWriter =  dynamic_cast<FileWriter*>(ioWriter);
                    if (fileWriter != nullptr) {
                        LOG_INFO("Reference index offset is %ld", fileWriter->getCurrentPos());
                        dynamicFileMeta.getMetaData("refe_index")["offset"] = fileWriter->getCurrentPos();
                    }
                    blockWriter->writeBlock(outBlockPtr);
                    PbgzManager::getInstance().updateWriteDataLen(outBlockPtr);
                    MemoryUtil::safeDeleteClass(outBlockPtr);
                    continue;
                }

                outputSortedCache.push_back(outBlockPtr);
                outputSortedCache.sort([](const RoughIOBlock* p1, RoughIOBlock* p2) {
                    return p1->getBlockId() <= p2->getBlockId();
                });
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outblockPtr = outputSortedCache.front();
                    if (outblockPtr->getBlockId() == blockId2Write) {
                        setDataBlockPosition(outblockPtr->getBlockId());
                        blockWriter->writeBlock(outblockPtr);
                        blockId2Write++;
                        PbgzManager::getInstance().updateWriteDataLen(outblockPtr);
                        freeOutputPool.push(outblockPtr);
                        outputSortedCache.pop_front();
                    } else {
                        break;
                    }
                } 
            }
        }
        MemoryUtil::safeDeleteClass(blockWriter);
        return 0;
    };

    writeThread = std::thread(writerTask);
    return 0;
}

bool PbgzEngine::initRefGeneForDecomress(PbgzBlockReader* blockReader) {
    if (blockReader == nullptr) {
        return false;
    }
    if (!baseFileMeta.getMetaData("refe").isObject()) {
        LOG_INFO("No reference");
        return true;
    }
    Json::Value& metaRefe = baseFileMeta.getMetaData("refe");
    if (dynamicFileMeta.getMetaData().isMember("refe")) {
        LOG_INFO("Read meta data from dynamic file meta");
        metaRefe = dynamicFileMeta.getMetaData("refe");
    }

    std::string fastaName = metaRefe["fasta_name"].asString();
    int64_t fastaLength = metaRefe["fasta_len"].asInt64();
    std::string niName = metaRefe["ni_name"].asString();

    if (metaRefe["blocks"].asInt64() > 0) {
        /*  There is pack reference in compression package, directly unpack */
        if(unpackReference(blockReader, metaRefe)) {
            return true;
        }
    }
    
    std::string fastaNameInput = parameter.referenceGenic;
    if (fastaNameInput.empty()) { /* No reference file specified */
        fprintf(stderr, "need to specify the following FASTA file:\n\n");
        fprintf(stderr, "\t%-12s : %s\n", "File Name", metaRefe["fasta_name"].asString().c_str());
        fprintf(stderr, "\t%-12s : %ld\n", "File Length", metaRefe["fasta_len"].asInt64());
        fprintf(stderr, "\t%-12s : %s\n", "File MD5", metaRefe["fasta_md5"].asString().c_str());
        LOG_ERROR("reference file needs to be specified to complete decompression");
        return false;
    }

    fastaNameInput = PathUtil::getAbsPath(fastaNameInput);
    int64_t fastqFileLenInput = PathUtil::getFileSize(fastaNameInput);
    /* check whether fasta is matched */
    if (PathUtil::getFileName(fastaNameInput) != fastaName) {
        fastaNameInput = PathUtil::getAbsPath(fastaNameInput);
        LOG_ERROR("initialize reference failed: used fasta %s, should be %s", fastaNameInput.c_str(), fastaName.c_str());
        return false;
    }
    if (fastaLength != fastqFileLenInput) {
        LOG_ERROR("initialize reference failed: used fasta file len %ld, should be %ld", fastqFileLenInput, fastaLength);
        return false;
    }

    Reference refeCheck(parameter.referenceGenic, parameter.threadNum);
    std::string niNameInput;
    refeCheck.getNiFileFromReference(niNameInput);
    niNameInput = PathUtil::getAbsPath(niNameInput);
    if (PathUtil::getFileName(niNameInput) != niName)  {
        LOG_ERROR("initialize reference failed: used ni file %s, should be %s", niNameInput.c_str(), niName.c_str());
        return false;
    }
    /* matched, do make index */
    pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    if (pRefGene == nullptr) {
        return false;
    }
    if (!pRefGene->initSquashByNiFile()) {
        LOG_ERROR("initialize reference failed");
        return false;
    }
    
    return true;
}

bool PbgzEngine::initRefeIndexForDecompress(PbgzBlockReader* blockReader) {
    if (blockReader == nullptr) {
        return false;
    }

    if (!dynamicFileMeta.getMetaData().isMember("refe_index")) {
        return true;
    }

    FileReader* pFileReader = dynamic_cast<FileReader*>(ioReader);
    if (pFileReader == nullptr) {
        return false;
    }

    size_t readOffset = pFileReader->getCurrentPos(); // 备份读取位置
    
    Json::Value refeIndexMeta = dynamicFileMeta.getMetaData("refe_index");
    int64_t offset = refeIndexMeta["offset"].asInt64();
    LOG_INFO("Reference index offset = %ld", offset);
    pFileReader->seekIO(offset);
    RoughIOBlock* refeIndexBlock = MemoryUtil::safeNewClass<RoughIOBlock>();
    blockReader->readBlock(refeIndexBlock);
    coder_json cmeta;
    Json::Value refIndexBlockMeta;
    cmeta.decoder(refeIndexBlock->getMetaBuffer(), refeIndexBlock->getMetaLen(), refIndexBlockMeta);
    uint32_t srcLen = refIndexBlockMeta["srclen"].asUInt();
    uint32_t dstLen = refIndexBlockMeta["dstlen"].asUInt(); 
    if (dstLen != refeIndexBlock->getDataLen()) {
        LOG_ERROR("RefeIndex invalid, dstlen = %d, expecet = %d.", refeIndexBlock->getDataLen(), dstLen);
        pFileReader->seekIO(readOffset);
        return false;
    }
    Json::Value& refePosDesc = refIndexBlockMeta["refe_desc"];
    Json::Value::iterator iter = refePosDesc.begin();
    for(; iter != refePosDesc.end(); ++iter) {
        std::string posPair = iter->asString();
        size_t pos = 0;
        if ((pos = posPair.find(",")) != std::string::npos)  {
            int64_t offset = std::stoll(posPair.substr(0, pos));
            int64_t length = std::stoll(posPair.substr(pos + 1));
            pRefGene->getRefeDescPos()[iter.key().asString()] = std::make_pair(offset, length);
        }
    }

    if (refIndexBlockMeta["coder"]["magic"].asString() != "coder_ppmd") {
        LOG_ERROR("coder name invalid.", refIndexBlockMeta["coder"]["magic"].asCString());
        pFileReader->seekIO(readOffset);
        return false;
    }
    uint8_t* outBuffer = MemoryUtil::safeAlloc<uint8_t>(srcLen);
    coder_io refIndexIO(refeIndexBlock->getBuffer(), refeIndexBlock->getDataLen());
    refIndexIO.meta = refIndexBlockMeta["coder"];
    coder_ppmd coder(&refIndexIO);
    coder.decode(outBuffer, srcLen);
    uint32_t readLength = 0;
    while(readLength < srcLen) {
        uint32_t blockId = *(uint32_t*)(outBuffer + readLength);
        readLength += sizeof(uint32_t);
        uint32_t refeIndexSize = *(uint32_t*)(outBuffer + readLength);
        readLength += sizeof(uint32_t);
        for (uint32_t idx = 0; idx < refeIndexSize; ++idx) {
            int64_t refeIndexPos = *(int64_t*)(outBuffer + readLength);
            blockRefeIndex[refeIndexPos] = blockId;
            readLength += sizeof(int64_t);
        }
    }
    pFileReader->seekIO(readOffset);

    // Json::Value blockPos = dynamicFileMeta.getMetaData("block_pos");
    // Json::Value::iterator jsonIter = blockPos.begin();
    // for (; jsonIter!= blockPos.end(); ++jsonIter) {
    //     LOG_INFO("BlockId:%d, FileOffset=%d", jsonIter.key().asUInt(), jsonIter->asInt64());
    // }
    // pRefGene->showRefeDesPosInfo();
    // auto mapItem = blockRefeIndex.begin();
    // for (; mapItem != blockRefeIndex.end(); ++mapItem) {
    //     LOG_INFO("MappedPosition:%ld, blockId:%d", mapItem->first, mapItem->second);
    // }
    
    return true;
}

bool PbgzEngine::initReferenceForCompress() {
    if (!parameter.referenceGenic.empty()){
        pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    }
    if (pRefGene) {
        if (!pRefGene->makeIndex()) {
            MemoryUtil::safeDeleteClass(pRefGene);
        }
    }
    // Refresh file meta
    if (pRefGene) {
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath()); /* Contains md5 information, used for decompression verification */
        std::vector<RoughIOBlock*> blockVec;
        bool isPackRefeInHeader = !parameter.isUnpackRef && parameter.outputFile == "/dev/stdout";
        if (isPackRefeInHeader) {
            int64_t maxRefLen = 0;
            int64_t totolEncLen = 0;
            int64_t refBlockCount = packReference(blockVec, maxRefLen, totolEncLen, false);
            refeMeta["max_block_len"] = maxRefLen;
            refeMeta["blocks"] = refBlockCount;
            // Throw to write file thread to write
            for (RoughIOBlock* item : blockVec) {
                item->setBlockId(blockCount++);
                outputDataPool.pushForce(item);
            }
            blockVec.clear();
        }
        baseFileMeta.setMetaData("refe",refeMeta);
    }
    return true;
}

bool PbgzEngine::unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta) {
    int64_t refeSquashLen = refeMeta["squash_len"].asInt64();
    std::string fastaName = refeMeta["fasta_name"].asString();
    int64_t maxLen = refeMeta["max_block_len"].asInt64();
    int64_t blocks = refeMeta["blocks"].asInt64();
    std::string md5Packed = refeMeta["md5"].asString();
    int64_t offset = 0;
    if (refeMeta.isMember("offset")) {
        offset = refeMeta["offset"].asInt64();
        LOG_INFO("Reference offset = %d", offset);
    }

    FileReader* fileReader = nullptr;
    size_t readPos = 0;
    if (offset > 0) {
        fileReader = dynamic_cast<FileReader*>(ioReader);
        if (fileReader == nullptr) {
            return false;
        }
        readPos = fileReader->getCurrentPos();
        fileReader->seekIO(offset);
    }
    int64_t pcnt = parameter.threadNum;
    RoughIOBlock* block[pcnt];
    for (int64_t n = 0; n < pcnt; n++) {
        block[n] = MemoryUtil::safeNewClass<RoughIOBlock>(maxLen);
    }
    
    BlockingQueue<RoughIOBlock*> inputPool(pcnt);
    BlockingQueue<RoughIOBlock*> inputBlock(1);
    pRefGene =  MemoryUtil::safeNewClass<Reference>(fastaName, parameter.threadNum);
    uint8_t* refeSquash = pRefGene->initSquashByStream(refeSquashLen);

    std::vector<std::thread> tpools;
    for (int64_t n = 0; n < blocks; n++) {
        if (n < pcnt) {
            blockReader->readBlock(block[n]);
            inputPool.push(block[n]);
            tpools.push_back(std::thread([&inputPool, &inputBlock, &refeSquash, &n]() {
                for (;;) {
                    RoughIOBlock* currBlock = inputPool.get();
                    if (currBlock == nullptr) {
                        return;
                    }
                    /* First parse out current block's meta information, get corresponding line's stream information and encoder information */
                    coder_json cmeta;
                    Json::Value refBlockMeta;
                    cmeta.decoder(currBlock->getMetaBuffer(), currBlock->getMetaLen(), refBlockMeta);
                    if (refBlockMeta["coder"]["magic"].asString() == "coder_ppmd") {
                        int64_t offset = refBlockMeta["offset"].asInt64();
                        int32_t srcLen = refBlockMeta["srclen"].asInt();
                        int32_t dstLen = refBlockMeta["dstlen"].asInt();

                        /* Point to the position of data to be decompressed in current block of reference */
                        if (currBlock->getDataLen() != dstLen) {
                            LOG_ERROR("reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                refBlockMeta["block"].asInt(), dstLen, currBlock->getDataLen());
                            return;
                        }
                        coder_io refeIo(currBlock->getBuffer(), currBlock->getDataLen());
                        refeIo.meta = refBlockMeta["coder"];
                        coder_ppmd cppmd(&refeIo);
                        int32_t len = cppmd.decode(refeSquash + offset, srcLen);
                        if (srcLen != len) {
                            LOG_ERROR("reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                    refBlockMeta["block"].asInt(), srcLen, len);
                            return;
                        } 
                    } else {
                        LOG_ERROR("undefined coder %s", refBlockMeta["coder"]["magic"].asString().c_str());
                        return;
                    }
                    inputBlock.push(currBlock);
                }
            }));
        }
        else {
            RoughIOBlock* current = inputBlock.get(); 
            blockReader->readBlock(current);
            inputPool.push(current);
        }
    }

    if (fileReader != nullptr && readPos != 0) {
        fileReader->seekIO(readPos);
    }

    for (int i= 0; i < pcnt; ++i) {
        (void)inputBlock.get();
        inputPool.push(nullptr);
    }
    for (auto &t : tpools) {
        if (t.joinable()) {
            t.join();
        }
    }
    for (int64_t n = 0; n < pcnt; n++) {
        MemoryUtil::safeDeleteClass(block[n]);
    }
    return true;
}

/*  Save reference genome */
int64_t PbgzEngine::packReference(std::vector<RoughIOBlock*>& blockVec, int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef) {
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

    // Start threads
    for (uint32_t idx = 0; idx < pcnt; ++idx) {
        tpools.push_back(std::thread([&inputPool, &m, &blockVec, &refe, &each, &maxBlockLen, &totalEncLen, &isSanitizeRef]() {
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
            
                int64_t currBlockLen = metaString.length() + refeIo.data_len;
                if (currBlockLen >= plen) {
                    LOG_ERROR("reference block pack failed");
                    break;
                }
                /* Write compressed stream of current block's meta */
                RoughIOBlock* outBlock = MemoryUtil::safeNewClass<RoughIOBlock>();
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);

                m.lock();
                blockVec.push_back(outBlock);
                maxBlockLen = (currBlockLen > maxBlockLen) ? currBlockLen : maxBlockLen;
                totalEncLen += currBlockLen;
                m.unlock();
            }
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

    std::sort(blockVec.begin(), blockVec.end(), [](RoughIOBlock* a, RoughIOBlock* b) {
        return a->getBlockId() < b->getBlockId(); // 降序条件
    });

    return block;
}

void PbgzEngine::updateReferenceOffset(int64_t offset) {
    if (refeOffsetFLag) {
        return;
    }
    LOG_INFO("Reference offset is %ld", offset);
    refeOffsetFLag = true;
    dynamicFileMeta.getMetaData("refe")["offset"] = offset;
    return;
}

void PbgzEngine::resetReferenceOffset() {
    refeOffsetFLag = false;
}

 void PbgzEngine::setDataBlockPosition(uint32_t blockId) {
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }
    if (!dynamicFileMeta.getMetaData().isMember("block_pos")) {
        Json::Value meta;
        meta[blockId] = fileWriter->getCurrentPos();
        dynamicFileMeta.setMetaData("block_pos", meta);
    } else {
        dynamicFileMeta.getMetaData("block_pos")[blockId] = fileWriter->getCurrentPos();
    }
    return;
 }

int32_t PbgzEngine::packBlockRefePosIdx() {
    uint32_t packLength = 0;
    packLength += blockRefePos.size() * (sizeof(uint32_t) + sizeof(uint32_t));
    auto iter = blockRefePos.begin();
    for(; iter != blockRefePos.end(); ++iter) {
        packLength += iter->second.size() * sizeof(int64_t);
    }

    uint8_t* buffer = MemoryUtil::safeAlloc<uint8_t>(packLength);
    if (buffer == nullptr) {
        return -1;
    }

    RoughIOBlock* outBlock = MemoryUtil::safeNewClass<RoughIOBlock>();
    if (outBlock == nullptr) {
        MemoryUtil::safeFree(buffer);
        return -1;
    }

    iter = blockRefePos.begin();
    uint32_t offset = 0;
    for(; iter != blockRefePos.end(); ++iter) {
        uint32_t blockId = iter->first;
        memcpy(buffer + offset, &blockId, sizeof(uint32_t));
        offset += 4;
        uint32_t length = iter->second.size();
        memcpy(buffer + offset, &length, sizeof(uint32_t));
        offset += 4;
        memcpy(buffer + offset, iter->second.data(), length * sizeof(int64_t));
        offset += length * sizeof(int64_t);
    }

    coder_io idxIo(outBlock->getBuffer(), outBlock->getBufferSize());
    coder_ppmd cppmd(&idxIo);
    cppmd.encode(buffer, packLength);
    cppmd.encode_flush();
    outBlock->setDataLen(idxIo.data_len);

    Json::Value meta;
    meta["srclen"] = (Json::Value::Int)packLength;
    meta["dstlen"] = (Json::Value::Int)idxIo.data_len;
    meta["coder"] = idxIo.meta;
    
    Json::Value refePosDesMeta;
    RefDescInfo& refDesc = pRefGene->getRefeDescPos();
    RefDescInfoIter refDescIter = refDesc.begin();
    for (; refDescIter != refDesc.end(); ++refDescIter ) {
        refePosDesMeta[refDescIter->first] = std::to_string(refDescIter->second.first) + "," + std::to_string(refDescIter->second.second);
    }
    meta["refe_desc"] = refePosDesMeta;

    coder_json cmeta;
    std::string metaString;
    cmeta.encoder(meta, metaString);
    memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
    outBlock->setMetaLen(metaString.length());
    outBlock->setBlockId(blockCount++);
    outBlock->setBlockType(REFERENCE_INDEX);
    
    Json::Value refeIndexMeta;
    refeIndexMeta["blocks"] = 1;
    dynamicFileMeta.setMetaData("refe_index", refeIndexMeta); 

    outputDataPool.push(outBlock);
    MemoryUtil::safeFree(buffer);
    return 0;
}