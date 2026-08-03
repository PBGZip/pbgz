/*
 * pbgz_engine.cpp - Source file for pbgz project
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
#include "utils/timer.h"
#include "log/logger.h"
#include "utils/path_util.h"
#include "hardware.h"
#include "coder.h"
#include "pbgz_manager.h"
#include "config_manager.h"
#include "block_wrapper.h"
#include "actuator.h"


PbgzEngine::PbgzEngine(const PbgzParameter&  para) {
    parameter = para;
    ioReader = nullptr;
    ioWriter = nullptr;
    syncFlag = false;
    blockId2Write = 0;
    blockCount = 0;

    freeInputPool = std::make_unique<BlockingQueueType>();
    inputDataPool = std::make_unique<BlockingQueueType>();
    freeOutputPool = std::make_unique<BlockingQueueType>();
    outputDataPool = std::make_unique<BlockingQueueType>();
}

PbgzEngine::~PbgzEngine() {
    for (auto& th : workThreads) {
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

    while(!freeInputPool->empty()) {
        RoughIOBlock* inPtr = freeInputPool->get();
        MemoryUtil::safeDeleteClass(inPtr);
    }

    while(!inputDataPool->empty()) {
        RoughIOBlock* inPtr = inputDataPool->get();
        MemoryUtil::safeDeleteClass(inPtr);
    }

    while(!freeOutputPool->empty()) {
        RoughIOBlock* outPtr = freeOutputPool->get();
        MemoryUtil::safeDeleteClass(outPtr);
    }

    while(!outputDataPool->empty()) {
        RoughIOBlock* outPtr = outputDataPool->get();
        MemoryUtil::safeDeleteClass(outPtr);
    }
}

uint32_t PbgzEngine::getBlockSize() {
    return ConfigManager::getInstance().getBlockSizeByCompressLevel(parameter.compressLevel);
}

int32_t PbgzEngine::init() {
    // Register allocation functions required by coder
    coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
    coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
    coder_ns::register_exit_proc(pbgzExitProc);
    coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    coder_ns::resister_logger_proc(coderLog);
    coder_ns::initFcCoder();

    // Create queues
    freeInputPool->setCapility(parameter.threadNum);
    inputDataPool->setCapility(parameter.threadNum);
    freeOutputPool->setCapility(parameter.threadNum << 1);
    outputDataPool->setCapility(parameter.threadNum << 1);

    uint32_t blockBufferSize = getBlockSize();
    // First push empty blocks to free queue
    for (uint32_t i = 0; i < freeInputPool->getCapility(); ++i) {
        RoughIOBlock* inPtr = MemoryUtil::safeNewClass<RoughIOBlock>(blockBufferSize);
        if (inPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }
        freeInputPool->push(inPtr);
    }

    for (uint32_t j = 0; j < freeOutputPool->getCapility(); ++j) {
        RoughIOBlock* outPtr = MemoryUtil::safeNewClass<RoughIOBlock>(blockBufferSize);
        if (outPtr == nullptr) {
            LOG_ERROR("PbgzEngine init failed.");
            return -1;
        }
        freeOutputPool->push(outPtr);
    }

    if (parameter.inputFile == STDIN) {
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

    if (parameter.outputFile == STDOUT) {
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

int32_t PbgzEngine::start() {
    printHeadInfo();
    Timer costTimer(true);

    int32_t ret = startEnginePreProc();
    if (ret != 0) {
        LOG_ERROR("call startEnginePreProc failed, ret = %d", ret);
        return ret;
    }

    ret = prepareFileMeta();
    if (ret != 0) {
        LOG_ERROR("call prepareFileMeta failed, ret = %d", ret);
        return ret;
    }

    ret = startWriteTask();
    if (ret != 0) {
        LOG_ERROR("call startWriteTask failed, ret = %d", ret);
        return ret;
    }

    ret = startWorkPreProc(); 
    if (ret != 0) {
        LOG_ERROR("call startWorkPreProc failed, ret = %d", ret);
        return ret;
    }

    ret =startWorkTask();
    if (ret != 0) {
        LOG_ERROR("call startWorkTask failed, ret = %d", ret);
        return ret;
    }

    ret =startReadPreProc();
    if (ret != 0) {
        LOG_ERROR("call startReadPreProc failed, ret = %d", ret);
        return ret;
    }

    ret =startReadTask();
    if (ret != 0) {
        LOG_ERROR("call startReadTask failed, ret = %d", ret);

        // Termination symbol for encoding
        for (uint32_t i = 0; i < parameter.threadNum; ++i) {
            inputDataPool->push(nullptr);
        }
        // Output termination symbol
        outputDataPool->push(nullptr);
        return ret;
    }

    for (auto& th : workThreads) {
        if (th.joinable()) {
            th.join();
        }
    }

    ret = startEnginePostProc();
    if (ret != 0) {
        LOG_ERROR("call startEnginePostProc failed.");
        return -1;
    }

    // Write end marker
    outputDataPool->push(nullptr);
    writeThread.join();
    
    printTailInfo(costTimer);

    return 0;
}

int64_t PbgzEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    RoughIOBlock* blockPtr = freeInputPool->get();
    if (blockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }
    blockPtr->reset();

    int64_t ret = blockReader->readBlock(blockPtr, fileType);
    if (ret <= 0) {
        freeInputPool->push(blockPtr);
        return ret;
    }
    
    if (blockPtr->getBlockId() == 0 && blockPtr->getTotalDataLen() > 0){ 
        fileType = blockPtr->getBlockType();
    }

    updateInputStatics(blockPtr);
    inputDataPool->push(blockPtr);
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
    pthread_setname_np("readtask");
    BlockReader* blockReader = createBlockReader();
    if (blockReader == nullptr) {
        return -1;
    }
    readBlocks(blockReader);
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        inputDataPool->push(nullptr);
    }
    releaseBlockReader(blockReader);
    return 0;
}

int32_t PbgzEngine::startWriteTask() {
    auto writerTask = [this]() -> int32_t {
        pthread_setname_np("writetask");
        BlockWriter* blockWriter = createBlockWriter();
        if (blockWriter == nullptr) {
            PbgzManager::getInstance().exitProc(-1, "Inner error");
            return -1;
        }

        while (true) {
            RoughIOBlock* outBlockPtr = outputDataPool->get();
            if (outBlockPtr == nullptr) {     // Got end marker
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outBlockTmp = outputSortedCache.front();
                    if (outBlockTmp != nullptr) {
                        writeOneBlock(blockWriter, outBlockTmp);
                    }
                }
                writeFilePostProc(blockWriter);
                break;
            } else {
                writeBlockPreProc(blockWriter, outBlockPtr);
                outputSortedCache.push_back(outBlockPtr);
                outputSortedCache.sort([](const RoughIOBlock* p1, RoughIOBlock* p2) {
                    return p1->getBlockId() <= p2->getBlockId();
                });
                while (!outputSortedCache.empty()) {
                    RoughIOBlock* outblockPtr = outputSortedCache.front();
                    if (outblockPtr->getBlockId() == blockId2Write) {
                        writeOneBlock(blockWriter, outblockPtr);
                    } else {
                        break;
                    }
                } 
            }
        }
        releaseBlockWriter(blockWriter);
        return 0;
    };

    writeThread = std::thread(writerTask);
    return 0;
}

void PbgzEngine::writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr) {
    blockWriter->writeBlock(outblockPtr);
    blockId2Write++;
    updateOutputStatics(outblockPtr);
    freeOutputPool->push(outblockPtr);
    outputSortedCache.pop_front();
    return;
}

int32_t PbgzEngine::startWorkTask() {
    auto coderTask = [this](int32_t id) {
        pthread_setname_np(std::string("codertask_").append(std::to_string(id)).c_str());
        if (id > 0 ) {
            std::unique_lock<std::mutex> lock(mutex);
            conditionVar.wait(lock);
        } else if (parameter.verbose) {
            /*
             * 只有 id == 0 的线程会不经等待地跳过上面的 wait，直接进入
             * while 循环处理第一个数据块；此时其余 id > 0 的工作线程都还
             * 阻塞在 conditionVar 上。也就是说本分支在整次运行里天然
             * 只会被执行一次，无须再额外加去重标志。
             * 这里输出一行提示，让 -v 的使用者看得见"第一块串行、编码器
             * 决策尚未拍板"这个短暂但关键的阶段——它是首块试压结果决定
             * 全文件后续 coder 选择、进而保证压缩产物可复现的窗口期。
             */
            fprintf(stderr, "[流水线] 串行压缩开始，编码器决策未完成，其余线程等待中\n");
        }

        LOG_INFO("Coder task (%d) begin to running!", id);

        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool->get();
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
                if (isNeedNotify(true)) {
                    conditionVar.notify_all();
                }
                break;
            }
            RoughIOBlock* outBlockPtr = freeOutputPool->get();
            if (outBlockPtr == nullptr) {
                LOG_ERROR("Get free output block failed.");
                freeInputPool->push(inBlockPtr);
                break;
            }
            outBlockPtr->reset();
            outBlockPtr->setBlockId(inBlockPtr->getBlockId());

            Actuator* pActuator = createActuator(inBlockPtr, outBlockPtr);
            if (pActuator == nullptr) {
                /*
                 * 单个块建不出执行器，绝不能退出循环：工作线程一旦提前结束，
                 * 剩余线程要独自消化整条输入队列，全部退出则读线程会永远阻塞在
                 * 有界队列上。这里与下面 actuatorProc 失败走完全相同的收尾——
                 * 归还输入块、推一个同 id 的零长块让写线程继续推进，然后继续取下一块。
                 */
                LOG_ERROR("Create actuator failed for block(%ld)", inBlockPtr->getBlockId());
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(inBlockPtr->getBlockId());
                outputDataPool->push(outBlockPtr);
                continue;
            }

            int32_t ret = actuatorProc(pActuator, inBlockPtr, outBlockPtr);
            if (ret != 0) {
                LOG_ERROR("Coder task failed for block(%d)", inBlockPtr->getBlockId());
                fprintf(stderr, "Warning: block(%ld) process failed.\n", inBlockPtr->getBlockId());
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(inBlockPtr->getBlockId());
                // When an error occurs, push a block with length 0 but correct ID, the write thread ignores blocks with length 0 to prevent thread waiting
                outputDataPool->push(outBlockPtr);
                MemoryUtil::safeDeleteClass(pActuator);
                continue;
            }
            
            outBlockPtr->setBlockType(inBlockPtr->getBlockType());

            if (isNeedNotify(pActuator->getNotifyFlag())) {
                /*
                 * 走到这里意味着：首块已经完整跑完（含 CodecSelector 的
                 * 试压/挑选），编码器决策终于"拍板"，可以叫醒其余工作线程
                 * 进入并行阶段。isNeedNotify 内部由 syncFlag 自锁——首次
                 * 传入 true 才返回 true，之后永久返回 false——所以这段
                 * 分支在整次运行里天然只会走一次，无须外层再加一层去重。
                 * 打印必须是单次 fprintf 保证原子性，因为此刻其他线程即将
                 * 被唤醒并可能立刻向 stderr 写日志。
                 */
                if (parameter.verbose) {
                    fprintf(stderr, "[流水线] 拍板完成，并行压缩启动（%u 线程）\n", parameter.threadNum);
                }
                conditionVar.notify_all();
            }
    
            MemoryUtil::safeDeleteClass(pActuator);
            
            freeInputPool->push(inBlockPtr);
            outputDataPool->push(outBlockPtr);
        }
    };
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        workThreads.emplace_back(std::thread(coderTask, i));
    }
    return 0;
}

bool PbgzEngine::isNeedNotify(bool flag) {
    if (flag && !syncFlag) {
        syncFlag = flag;
        return true;
    }

    return false;
}

int32_t PbgzEngine::actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*) {
    return actuator->process();
} 
