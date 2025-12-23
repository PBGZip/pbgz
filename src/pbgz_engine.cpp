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

int32_t PbgzEngine::init() {
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

    int32_t ret = engineStartPreProc();
    if (ret != 0) {
        LOG_ERROR("engine pre proc failed.");
        return -1;
    }

    ret = startWriteTask();
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

    ret = engineStartAfterProc();
    if (ret != 0) {
        LOG_ERROR("engine after proc failed.");
        return -1;
    }

    // Write end marker
    outputDataPool.push(nullptr);
    writeThread.join();
    PbgzManager::getInstance().printTailInfo(costTimer, isPrintRatio());
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

    if (blockPtr->getBlockId() == 0 && blockPtr->getTotalDataLen() > 0){ 
        fileType = blockPtr->getBlockType();
        PbgzManager::getInstance().printFileType(fileType);
    }

    PbgzManager::getInstance().updateReadDataLen(blockPtr);
    inputDataPool.push(blockPtr);
    if (ret > 0) {
        blockCount++;
    }
    return ret;
}

void PbgzEngine::readBlocks(BlockReader* blockReader) {
    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        ret = readOneBlock(blockReader, fileType);
    } while(ret > 0 || ret == -2);

    return;
}

int32_t PbgzEngine::startReadTask() {
    pthread_setname_np(pthread_self(), "readtask");
    BlockReader* blockReader = createBlockReader();
    if (blockReader == nullptr) {
        return -1;
    }
    readBlocks(blockReader);
    // Reached end of file or error, insert empty block as end marker to data queue
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        inputDataPool.push(nullptr);
    }

    MemoryUtil::safeDeleteClass(blockReader);
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
                if (parameter.isMakeIndex) {
                    fprintf(stderr, "Fastq file will not make index.");
                    parameter.isMakeIndex = false;
                }
                pActuator = MemoryUtil::safeNewClass<FastqActuator>(inBlockPtr, outBlockPtr, pRefGene);
                pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
                if (pActuator == nullptr) {
                    break;
                }
            } else if (inBlockPtr->getBlockType() == BINARY) {
                if (parameter.isMakeIndex) {
                    fprintf(stderr, "Binary file will not make index.");
                    parameter.isMakeIndex = false;
                }
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

            int32_t ret = actuatorProc(pActuator, inBlockPtr, outBlockPtr);
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
        BlockWriter* blockWriter = createBlockWriter();
        if (blockWriter == nullptr) {
            PbgzManager::getInstance().exitProc(-1, "Inner error");
            return -1;
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


void PbgzEngine::updateReferenceOffset(int64_t offset) {
    if (refeOffsetFLag) {
        return;
    }
    LOG_DEBUG("Reference offset is %ld", offset);
    refeOffsetFLag = true;
    dynamicFileMeta.getMetaData("refe")["offset"] = offset;
    return;
}

void PbgzEngine::resetReferenceOffset() {
    refeOffsetFLag = false;
}
