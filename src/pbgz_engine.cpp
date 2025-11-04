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
#include "hardware.h"

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
    } else {
        if (PathUtil::isGzFile(parameter.inputFile)) {
            bool isSupportSimd = Hardware().isSupportSimd();
            if (isSupportSimd) {
                ioReader = MemoryUtil::safeNewClass<FastGzFileReader>(parameter.inputFile);
            } else {
                ioReader = MemoryUtil::safeNewClass<GzFileReader>(parameter.inputFile, parameter.threadNum);
            }
        } else {
            ioReader = MemoryUtil::safeNewClass<FileReader>(parameter.inputFile);
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
        } else {
            ioWriter = MemoryUtil::safeNewClass<PipeWriter>();
        }
    } else {
        if(parameter.isDecToGZ) {
            ioWriter = MemoryUtil::safeNewClass<GzFileWriter>(parameter.outputFile, parameter.threadNum);
        } else {
            ioWriter = MemoryUtil::safeNewClass<FileWriter>(parameter.outputFile);
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
    // Release resources
    if (ioReader) {
        ioReader->closeIO();
        delete ioReader;
        ioReader = nullptr;
    }

    if (ioWriter) {
        ioWriter->closeIO();
        delete ioWriter;
        ioWriter = nullptr;
    }

    while(!freeInputPool.empty()) {
        RoughIOBlock* inPtr = freeInputPool.get();
        if (inPtr) {
            delete inPtr;
            inPtr = nullptr;
        }
    }

    while(!inputDataPool.empty()) {
        RoughIOBlock* inPtr = inputDataPool.get();
        if (inPtr) {
            delete inPtr;
            inPtr = nullptr;
        }
    }

    while(!freeOutputPool.empty()) {
        RoughIOBlock* outPtr = freeOutputPool.get();
        if (outPtr) {
            delete outPtr;
            outPtr = nullptr;
        }
    }

    while(!outputDataPool.empty()) {
        RoughIOBlock* outPtr = outputDataPool.get();
        if (outPtr) {
            delete outPtr;
            outPtr = nullptr;
        }
    }

    if (pRefGene) {
        delete pRefGene;
        pRefGene = nullptr;
    }
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
        return -1;
    }

    for (auto& th : coderThreads) {
        if (th.joinable()) {
            th.join();
        }
    }

    // Write end marker
    outputDataPool.push(nullptr);
    writeThread.join();
    PbgzManager::getInstance().printTailInfo(costTimer, parameter);
    return 0;
}

int32_t PbgzEngine::startReadTask() {
    pthread_setname_np(pthread_self(), "readtask");
    BlockReader* blockReader = nullptr;
    if (parameter.isDecompress) {   // Decompression mode, read content from pbgz file
        blockReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
        PbgzBlockReader* pbgzReader = dynamic_cast<PbgzBlockReader*>(blockReader);
        if (!initRefGeneForDecomress(pbgzReader)) {
            LOG_INFO("Init reference for decompress failed");
        }
        fileMeta = pbgzReader->getFileMeta();
    } else {  // Compression mode, read content from non-pbgz file
        blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
    }

    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return -1;
    }
    if (0 != blockReader->init()) {
        LOG_ERROR("BlockReader init failed");
        delete blockReader;
        blockReader = nullptr;
        return -1;
    }

    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        RoughIOBlock* blockPtr = freeInputPool.get();
        if (blockPtr == nullptr) {
            LOG_ERROR("Get free block failed.");
            return -1;
        }
        blockPtr->reset();

        ret = blockReader->readBlock(blockPtr, fileType);
        if (ret <= 0) {
            if (ret < 0) {
                LOG_ERROR("Read block failed.");
            }
            // Reached end of file or error, insert empty block as end marker to data queue
            for (uint32_t i = 0; i < parameter.threadNum; ++i) {
                inputDataPool.push(nullptr);
            }
            break;
        } 

        if (blockPtr->getBlockId() == 0){ 
            fileType = blockPtr->getBlockType();
            PbgzManager::getInstance().printFileType(fileType);
        }

        PbgzManager::getInstance().updateReadDataLen(blockPtr);
        inputDataPool.push(blockPtr);
    } while(ret > 0);

    delete blockReader;
    blockReader = nullptr;
    return 0;
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
                    if (!parameter.isDecompress && 0 != fastqActuator->preAnalysis()) { // Only need to analyze block content in compression scenario
                        pActuator = MemoryUtil::safeNewClass<BinaryActuator>(inBlockPtr, outBlockPtr);
                    }
                }
            } else if (inBlockPtr->getBlockType() == BINARY) {
                pActuator = MemoryUtil::safeNewClass<BinaryActuator>(inBlockPtr, outBlockPtr);
            } else {
                freeInputPool.push(inBlockPtr);
                LOG_ERROR("Not support block type: %d", inBlockPtr->getBlockType());
                break;
            }

            if (pActuator == nullptr) {
                freeInputPool.push(inBlockPtr);
                LOG_ERROR("Create actuator failed.");
                break;
            }

            int32_t ret = 0;
            if (parameter.isDecompress) {
                ret = pActuator->decompress();
            } else {
                ret = pActuator->compress();
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
        } else {   // Compression mode, file write in pbgz format
            blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
            PbgzBlockWriter* pbgzWriter =  dynamic_cast<PbgzBlockWriter*>(blockWriter);
            if (pbgzWriter == nullptr) {
                LOG_ERROR("pbgzWriter is NULL.");
                return -1;
            }
            pbgzWriter->setFileMeta(fileMeta);
        }
        if (blockWriter == nullptr) {
            LOG_ERROR("Failed to create block writer.");
            return -1;
        }
        
        blockWriter->init();

        int64_t blockId2Write = 0; 
        while (true) {
            RoughIOBlock* outBlockPtr = outputDataPool.get();
            if (outBlockPtr == nullptr) {     // Got end marker
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outblockPtr = outputSortedCache.front();
                    blockWriter->writeBlock(outblockPtr);
                    PbgzManager::getInstance().updateWriteDataLen(outblockPtr);
                    freeOutputPool.push(outblockPtr);
                    outputSortedCache.pop_front();
                }
                break;
            } else {
                if (outBlockPtr->getBlockType() == REFERENCE) {
                    /// Writing reference blocks is one-time, release after writing
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
    const PbgzFileMeta& fileMeta = blockReader->getFileMeta();
    if (!fileMeta.getMetaData("refe").isObject()) {
        LOG_INFO("No reference");
        return true;
    }
    
    Json::Value metaRefe = fileMeta.getMetaData("refe");
    std::string fastaName = metaRefe["fasta_name"].asString();
    int64_t fastaLength = metaRefe["fasta_len"].asInt64();
    std::string niName = metaRefe["ni_name"].asString();

    if (metaRefe["blocks"].asInt64() > 0) {
        /*  There is pack reference in compression package, directly unpack */
        unpackReference(blockReader);
    } else {
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
        if (fastaNameInput != fastaName) {
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
        if (niNameInput != niName)  {
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
    }
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
        refeMeta["fasta_name"] = PathUtil::getAbsPath(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getAbsPath(pRefGene->getNiFilePath()); /* Contains md5 information, used for decompression verification */
        std::vector<RoughIOBlock*> blockVec;
        if (!parameter.isUnpackRef) {
            int64_t maxRefLen = 0;
            int64_t totolEncLen = 0;
            int64_t refBlockCount = packReference(blockVec, maxRefLen, totolEncLen);
            refeMeta["max_block_len"] = maxRefLen;
            refeMeta["blocks"] = refBlockCount;
        }
        fileMeta.setMetaData("refe",refeMeta);
        // Throw to write file thread to write
        if (!parameter.isUnpackRef) {
            for (auto item : blockVec) {
                outputDataPool.pushForce(item);
            }
            blockVec.clear();
        }
    }
    return true;
}

void PbgzEngine::unpackReference(PbgzBlockReader* blockReader) {
    Json::Value metaRefe = blockReader->getFileMeta().getMetaData("refe");
    int64_t refeSquashLen = metaRefe["squash_len"].asInt64();
    std::string fastaName = metaRefe["fasta_name"].asString();
    int64_t maxLen = metaRefe["max_block_len"].asInt64();
    int64_t blocks = metaRefe["blocks"].asInt64();
    std::string md5Packed = metaRefe["md5"].asString();
    
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
    return;
}

/*  Save reference genome */
int64_t PbgzEngine::packReference(std::vector<RoughIOBlock*>& blockVec, int64_t &maxBlockLen, int64_t &totalEncLen) {
    int64_t block = 0, offset = 0;
    std::vector<std::thread> tpools;
    
    uint32_t pcnt = parameter.threadNum;
    int64_t each = (16 << 20);
    int64_t total = pRefGene->getSquashLength();
    int64_t remain = total;
    maxBlockLen = totalEncLen = 0;
    
    uint8_t *output = MemoryUtil::safeAlloc<uint8_t>(pcnt * each);
    
    BlockingQueue<RefeInfo> inputPool(pcnt * 2);
    Reference* refe = pRefGene;
    int64_t current;
    std::mutex m;

    // Start threads
    for (uint32_t idx = 0; idx < pcnt; ++idx) {
        tpools.push_back(std::thread([&inputPool, &output, &m, &blockVec, &refe, &each, &maxBlockLen, &totalEncLen]() {
            while(true) {
                RefeInfo refe2do = inputPool.get();
                int64_t plen = refe2do.second.second;
                if (plen == 0) {
                    break;
                }
                const uint8_t *p = refe->getSquash() + refe2do.second.first;
                refe->sanitizeRefSquash(refe2do.second.first, refe2do.second.second);

                coder_io refeIo(refe2do.first.second, each);
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
            
                m.lock();
                int64_t currBlockLen = metaString.length() + refeIo.data_len;
                if (currBlockLen >= plen) {
                    LOG_ERROR("reference block pack failed");
                    m.unlock();
                    break;
                }
                /* Write compressed stream of current block's meta */
                RoughIOBlock* outBlock = MemoryUtil::safeNewClass<RoughIOBlock>();
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);
                blockVec.push_back(outBlock);
                maxBlockLen = (currBlockLen > maxBlockLen) ? currBlockLen : maxBlockLen;
                totalEncLen += currBlockLen;
                m.unlock();
            }
        }));
    }

    while (remain > 0) {
        current = std::min(each, remain);
        inputPool.push(std::make_pair(std::make_pair(block, output + ((block % pcnt) * each)), std::make_pair(offset, current)));
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

    MemoryUtil::safeFree(output);
    return block;
}
