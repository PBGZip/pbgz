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

int PbgzEngine::init() {
    // 注册coder需要的注册函数
    coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
    coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
    coder_ns::register_exit_proc(pbgzExitProc);
    coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    coder_ns::resister_logger_proc(coderLog);

    // 创建队列
    freeInputPool.setCapility(parameter.threadNum);
    inputDataPool.setCapility(parameter.threadNum);
    freeOutputPool.setCapility(parameter.threadNum << 1);
    outputDataPool.setCapility(parameter.threadNum << 1);

    // 首先往空闲队列压入空的block
    for (uint32_t i = 0; i < freeInputPool.getCapility(); ++i) {
        RoughIOBlock* inPtr = new RoughIOBlock();
        if (inPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }

        freeInputPool.push(inPtr);
    }

    for (uint32_t j = 0; j < freeOutputPool.getCapility(); ++j) {
        RoughIOBlock* outPtr = new RoughIOBlock();
        if (outPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }

        freeOutputPool.push(outPtr);
    }

    if (parameter.inputFile == "/dev/stdin") {
        ioReader = MemoryUtil::safeNewClass<PipeReader>();
    } else {
        ioReader = MemoryUtil::safeNewClass<FileReader>(parameter.inputFile);
    }

    if (ioReader == nullptr) {
        LOG_ERROR("Create IO reader failed.");
        return -1;
    }

    ioReader->openIO();

    if (parameter.outputFile == "/dev/stdout") {
        ioWriter = new PipeWriter();
    } else {
        ioWriter = new FileWriter(parameter.outputFile);
    }

    if (ioWriter == nullptr) {
        LOG_ERROR("Create IO reader failed.");
        return -1;
    }

    ioWriter->openIO();

    return 0;
}

PbgzEngine::~PbgzEngine() {
    // 释放资源
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
}

int32_t PbgzEngine::start() {
    PbgzManager::getInstance().printHeadInfo(parameter);

    Timer costTimer(true);
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

    // 写入结束标记
    outputDataPool.push(nullptr);

    writeThread.join();

    PbgzManager::getInstance().printTailInfo(costTimer, parameter);

    return 0;
}

int32_t PbgzEngine::startReadTask() {
    BlockReader* blockReader = nullptr;
    if (parameter.isDecompress) {   // 解压模式，从pbgz文件读取内容
        blockReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
    } else {  // 压缩模式，从非pbgz文件读取内
        blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
    }

    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return -1;
    }
    if (0 != blockReader->init()) {
        LOG_ERROR("BlockReader init failed");
        return -1;
    }

    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        RoughIOBlock* blockPtr = freeInputPool.get();
        if (blockPtr == nullptr) {
            LOG_ERROR("Get free block failed。");
            return -1;
        }
        blockPtr->reset();

        ret = blockReader->readBlock(blockPtr, fileType);
        if (ret <= 0) {
            if (ret < 0) {
                LOG_ERROR("Read block failed.");
            }
            // 读到文件结尾或者出错, 往数据队列插入一个空的block作为结束标志
            for (int i = 0; i < parameter.threadNum; ++i) {
                usleep(20);
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
    auto coderTask = [this]() {
        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool.get();
            if (inBlockPtr == nullptr) {  // 读到空指针，表示拿到了结束标志
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
                pActuator = MemoryUtil::safeNewClass<FastqActuator>(inBlockPtr, outBlockPtr);
                FastqActuator* fastqActuator = dynamic_cast<FastqActuator*>(pActuator);
                if (fastqActuator != nullptr) {
                    if (!parameter.isDecompress && 0 != fastqActuator->preAnalysis()) { // 压缩场景才需要对块的内容进行分析
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
                delete pActuator;
                pActuator = nullptr;
                break;
            }
            delete pActuator;
            pActuator = nullptr;
            freeInputPool.push(inBlockPtr);
            outputDataPool.push(outBlockPtr);
        }
    };

    for (int8_t i = 0; i < parameter.threadNum; ++i) {
        coderThreads.emplace_back(std::thread(coderTask));
    }

    return 0;
}   

int32_t PbgzEngine::startWriteTask() {
    auto writerTask = [this]() {
        BlockWriter* blockWriter = nullptr;
        if (parameter.isDecompress) {  // 解压模式，文件写入为非pbgz格式
            blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
        } else {   // 压缩模式，写文件格式为pbgz格式
            blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
        }

        if (blockWriter == nullptr) {
            LOG_ERROR("Failed to create block writer.");
            return -1;
        }
        blockWriter->init();

        int64_t blockId2Write = 0; 
        while (true) {
            RoughIOBlock* outBlockPtr = outputDataPool.get();
            if (outBlockPtr == nullptr) {     // 拿到了结束标记
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outblockPtr = outputSortedCache.front();
                    blockWriter->writeBlock(outblockPtr);
                    PbgzManager::getInstance().updateWriteDataLen(outblockPtr);
                    outputSortedCache.pop_front();
                }
                break;
            } else {
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

        delete blockWriter;
        blockWriter = nullptr;
    };

    writeThread = std::thread(writerTask);
    return 0;
}





