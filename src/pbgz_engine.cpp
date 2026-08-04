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
#include <cstring>
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

    /*
     * 落盘必须在这里完成并检查。closeIO 里的 msync 是 mmap 写出真正回写磁盘的时机，
     * 原来只在析构里做，比 start() 返回还晚，失败无处可报。重复调用无害：
     * closeIO 会把 mappedAddress 和 fd 置空。
     */
    if (ioWriter != nullptr) {
        ioWriter->closeIO();
        int32_t writeErr = ioWriter->getWriteError();
        if (writeErr != 0) {
            LOG_ERROR("Flush output failed: %s", strerror(writeErr));
            fprintf(stderr, "Error: flushing output failed (%s), output is incomplete.\n", strerror(writeErr));
            return -1;
        }
    }

    /*
     * 块级错误的最终处理：任何一块压/解失败，整体结果都不可信。
     * 零长块只是保流水线不死锁的过渡，绝不能让进程以成功退出、留下残缺文件。
     */
    if (taskFailed.load()) {
        LOG_ERROR("Some blocks failed during processing; result is incomplete.");
        fprintf(stderr, "Error: some blocks failed, output is incomplete.\n");
        return -1;
    }

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

    readBlockPreProc(blockReader);

    int64_t ret = blockReader->readBlock(blockPtr, fileType);
    if (ret <= 0) {
        freeInputPool->push(blockPtr);
        return ret;
    }

    /*
     * 辅助块按定义就不属于数据流，对所有引擎都一样，所以规则只能写在这里。
     * 没人认领不是错误：那是旧版本读新格式时需要的前向兼容行为。
     */
    if (BlockUtil::isAuxiliaryBlock(blockPtr->getBlockType())) {
        PbgzBlockReader* pbgzReader = dynamic_cast<PbgzBlockReader*>(blockReader);
        (void)offerAuxBlock(blockPtr, pbgzReader ? pbgzReader->getCurrentFileIndex() : 0);
        freeInputPool->push(blockPtr);
        return -2;
    }

    const BlockIntake intake = intakeBlock(blockReader, blockPtr);
    if (intake != BlockIntake::DISPATCH) {
        freeInputPool->push(blockPtr);
        return (intake == BlockIntake::SKIP) ? -2 : -1;
    }

    if (blockPtr->getBlockId() == 0 && blockPtr->getTotalDataLen() > 0){
        fileType = blockPtr->getBlockType();
        /* 必须在 push 之前：入队即意味着某个工作线程可能立刻开始处理这个块。 */
        fileDecisionProc(blockPtr);
    }

    updateInputStatics(blockPtr);
    inputDataPool->push(blockPtr);
    if (ret > 0) {
        blockCount++;
    }

    return ret;
}

/*
 * 返回最后一次 readOneBlock 的结果：> 0 或 -2 不会出现（循环会继续），
 * 0 是干净 EOF，-1 是读错误。调用方必须据此区分"读完"与"读坏"——
 * 截断/损坏的输入绝不允许被当成正常 EOF 静默收尾。
 */
int64_t PbgzEngine::readBlocks(BlockReader* blockReader) {
    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        ret = readOneBlock(blockReader, fileType);
    } while(ret > 0 || ret == -2);

    return ret;
}

int32_t PbgzEngine::startReadTask() {
    pthread_setname_np("readtask");
    BlockReader* blockReader = createBlockReader();
    if (blockReader == nullptr) {
        return -1;
    }
    int64_t readRet = readBlocks(blockReader);
    for (uint32_t i = 0; i < parameter.threadNum; ++i) {
        inputDataPool->push(nullptr);
    }
    releaseBlockReader(blockReader);
    if (readRet < 0 && readRet != -2) {
        LOG_ERROR("Read task ended with an IO error: %ld", readRet);
        return -1;
    }
    return 0;
}

bool PbgzEngine::offerAuxBlock(RoughIOBlock* blockPtr, int64_t packageIndex) {
    for (size_t i = 0; i < auxConsumers.size(); ++i) {
        if (auxConsumers[i]->claim(blockPtr, packageIndex)) {
            return true;
        }
    }
    LOG_DEBUG("No consumer claimed auxiliary block, type=%d", blockPtr->getBlockType());
    return false;
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
            } else if (outBlockPtr->isSyncAux()) {
                /*
                 * 辅助块位置寻址，必须抢在 writeBlockPreProc / 顺序重排 / blockId2Write
                 * 之前处理：它不属于数据块 id 序列，一旦落进下面那套排序逻辑，
                 * 就会顶掉一个数据块的位置并让写指针永久错位。
                 */
                completeSyncAuxBlock(blockWriter, outBlockPtr);
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
        if (blockWriter->getWriteError() != 0) {
            taskFailed.store(true);
        }
        releaseBlockWriter(blockWriter);
        return 0;
    };

    writeThread = std::thread(writerTask);
    return 0;
}

void PbgzEngine::writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr) {
    if (blockWriter->writeBlock(outblockPtr) < 0) {
        LOG_ERROR("Write block %u failed.", outblockPtr->getBlockId());
        taskFailed.store(true);
    }
    blockId2Write++;
    updateOutputStatics(outblockPtr);
    freeOutputPool->push(outblockPtr);
    outputSortedCache.pop_front();
    return;
}

void PbgzEngine::completeSyncAuxBlock(BlockWriter* blockWriter, RoughIOBlock* block) {
    /*
     * 取位置必须紧贴 writeBlock 之前：writeBlock 写出的第一个字节就是容器魔数，
     * 所以此刻的写指针正是解压侧 readBlock 需要的那个偏移。
     */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    const int64_t offset = (fileWriter != nullptr) ? (int64_t)fileWriter->getCurrentPos() : -1;

    blockWriter->writeBlock(block);
    updateOutputStatics(block);
    /* 归还必须早于唤醒：发射者被唤醒后可能立刻再申请空闲块。 */
    freeOutputPool->push(block);

    {
        std::lock_guard<std::mutex> lock(auxEmitMutex);
        auxEmitOffset = offset;
        auxEmitDone = true;
    }
    auxEmitCond.notify_all();
}

int64_t PbgzEngine::emitSyncAuxBlock(RoughIOBlock* block) {
    if (block == nullptr) {
        return -1;
    }
    {
        std::lock_guard<std::mutex> lock(auxEmitMutex);
        auxEmitDone = false;
        auxEmitOffset = -1;
    }

    block->setSyncAux(true);
    outputDataPool->push(block);

    std::unique_lock<std::mutex> lock(auxEmitMutex);
    auxEmitCond.wait(lock, [this]() { return auxEmitDone; });
    return auxEmitOffset;
}

int32_t PbgzEngine::startWorkTask() {
    auto coderTask = [this](int32_t id) {
        pthread_setname_np(std::string("codertask_").append(std::to_string(id)).c_str());

        LOG_INFO("Coder task (%d) begin to running!", id);

        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool->get();
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
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
                taskFailed.store(true);
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(inBlockPtr->getBlockId());
                outputDataPool->push(outBlockPtr);
                continue;
            }

            /*
             * coder 层的 check_exit/coder_exit 现在抛 coder_exception（以前是在库里
             * 直接 _Exit 杀进程，不留日志）。这里是它的最终处理地：接住并转成
             * 与其它失败完全一致的返回值，汇入 taskFailed。catch(...) 也一并拦住：
             * 工作线程里漏出异常会直接 std::terminate，比丢块更难查。
             */
            int32_t ret = 0;
            try {
                ret = actuatorProc(pActuator, inBlockPtr, outBlockPtr);
            } catch (const coder_exception& e) {
                LOG_ERROR("Coder error on block(%ld): [%d] %s", inBlockPtr->getBlockId(), e.getCode(), e.what());
                ret = e.getCode() != 0 ? e.getCode() : -1;
            } catch (const std::exception& e) {
                LOG_ERROR("Unexpected exception on block(%ld): %s", inBlockPtr->getBlockId(), e.what());
                ret = -1;
            } catch (...) {
                LOG_ERROR("Unknown exception on block(%ld)", inBlockPtr->getBlockId());
                ret = -1;
            }
            if (ret != 0) {
                LOG_ERROR("Coder task failed for block(%d)", inBlockPtr->getBlockId());
                fprintf(stderr, "Warning: block(%ld) process failed.\n", inBlockPtr->getBlockId());
                taskFailed.store(true);
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(inBlockPtr->getBlockId());
                // When an error occurs, push a block with length 0 but correct ID, the write thread ignores blocks with length 0 to prevent thread waiting
                outputDataPool->push(outBlockPtr);
                MemoryUtil::safeDeleteClass(pActuator);
                continue;
            }

            outBlockPtr->setBlockType(inBlockPtr->getBlockType());

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

int32_t PbgzEngine::actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*) {
    return actuator->process();
}
