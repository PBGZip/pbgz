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
    /*
     * 输入侧容量取 threadNum + 1：SAM 头部独立成块后不参与 QUAL 先验预训练计数，
     * 若容量恰为 threadNum，读线程在读满队列后阻塞于取空闲块，永远凑不满预训练目标
     * 块数，finalizePretrain 不触发、coder 线程等不到 priorSettled -> 死锁。
     * 多出一个槽位让读线程能比 coder 多读一块，预训练在队列满前即可收尾。
     */
    freeInputPool->setCapility(parameter.threadNum + 1);
    inputDataPool->setCapility(parameter.threadNum + 1);
    freeOutputPool->setCapility(parameter.threadNum << 1);
    outputDataPool->setCapility(parameter.threadNum << 1);

    uint32_t blockBufferSize = getBlockSize();
    // First push empty blocks to free queue
    /*
     * 输入块初始只分配固定 1MB；读取时按 -l 决定的读取目标（getBlockSize()）读，
     * 容量不足由 reader 的 ensureCapacity 按需 realloc 扩容，避免按 -l 一次性分配
     * 大内存（如 -l 9 的 512MB×N）。输出块与文件头 block_size 仍按 -l。
     */
    for (uint32_t i = 0; i < freeInputPool->getCapility(); ++i) {
        RoughIOBlock* inPtr = MemoryUtil::safeNewClass<RoughIOBlock>(FIXED_INPUT_BLOCK_SIZE);
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
            /*
             * 标准 BAM 是 BGZF（gzip）流，逐字节看和 .gz 文件无法区分。但它的 gzip 是
             * 格式内层，不该走"透明 gz 解压"：透明解压用的是 Gz/FastGz reader，拿不到
             * 输入文件总长（fileDecisionProc 里 dynamic_cast<FileReader*> 失败），QUAL
             * 先验的"值不值得写"判据会退化，BAM 输入就会白白多写一个先验辅助块。
             * 因此 BAM 一律交给 FileReader + BamGzBlockReader（后者内部自行 inflate），
             * 块类型仍是 BAM，压缩行为与 SAM 输入对齐。
             */
            if (BlockUtil::isBamFile(parameter.inputFile)) {
                ioReader = MemoryUtil::safeNewClass<FileReader>(parameter.inputFile);
                LOG_INFO("Create FileReader (BAM).");
            } else {
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

    /* 注册输出文件：异常退出（如缺参考基因）时由 PbgzManager 清理删除。 */
    if (parameter.outputFile != STDOUT && parameter.outputFile != "-") {
        PbgzManager::getInstance().addOutputFile(parameter.outputFile);
    }

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
    }

    /*
     * 文件级决策必须发生在 push 之前（入队即意味着某个工作线程可能立刻开始处理）。
     * SAM 头部块（只含 @ 行、无数据行）不触发决策，推迟到首个数据块，否则编码器选型
     * 采不到字段样本只能退到默认编码器。必须在 push 之前调用，其余引擎行为不变。
     */
    if (!fileDecisionInvoked && blockReader->blockHasData(blockPtr)) {
        fileDecisionInvoked = true;
        fileDecisionProc(blockPtr);
    }

    /* 跨块预训练等按块累积工作（压缩引擎用），也在入队前执行。 */
    pretrainBlockProc(blockPtr);

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

    /* EOF/错误后收尾未完成的预处理（如跨块先验训练），必须在推终止符之前。 */
    readLoopPostProc();
    return ret;
}

int32_t PbgzEngine::startReadTask() {
    pthread_setname_np(pthread_self(), "readtask");
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
        pthread_setname_np(pthread_self(), "writetask");
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
        pthread_setname_np(pthread_self(), std::string("codertask_").append(std::to_string(id)).c_str());

        LOG_INFO("Coder task (%d) begin to running!", id);

        /* 压缩引擎的先验发布门闩：发布前 coder 线程等待通知，发布后才开始拉块压缩。 */
        workStartBarrier();

        /*
         * 首块串行化（perf 分支机制，见 firstCoderNotify）：id==0 的线程直接开始，
         * 其余线程等第 0 块处理完（preAnalysis 填充 SamInfo 等共享状态）再放行，
         * 避免并发 preAnalysis 在染色体表未就绪时把块误降级为 BINARY。
         *
         * 等待必须带谓词：第 0 块如果建不出执行器（如格式不支持的块类型），
         * 处理会立刻失败、放行通知可能早于其它线程开始等待而错过——不带谓词的裸等待
         * 会永久阻塞，让整个流水线死锁。谓词在持锁下复查 coderStartSync，先置位先得。
         */
        if (id > 0) {
            std::unique_lock<std::mutex> lock(coderStartMutex);
            coderStartCond.wait(lock, [this] { return coderStartSync; });
        }

        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool->get();
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
                /* 第 0 块未能完成（读错/全错）时也要放行，否则其余线程永久阻塞。 */
                if (firstCoderNotify(true)) {
                    coderStartCond.notify_all();
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

            /*
             * 越界错误的唯一查处点。coder_io 越界只置标志、不中断，执行器不一定把它
             * 转成返回值——原来 SAM 查了 12 处、FASTQ 和索引一处没查，越界写坏的数据
             * 就这样被当成"成功"写了出去。现在所有 coder_io 的错误都汇到执行器身上，
             * 这里问一次就覆盖全部流、全部执行器，新增的流也不会漏。
             *
             * 放在 ret 判断之前而不是合并进去：ret 已经是 0 时才需要补这一问，
             * ret 非 0 时失败原因更具体，不该被覆盖。
             */
            if (ret == 0 && !pActuator->ioError().ok()) {
                LOG_ERROR("Stream '%s' out of bounds on block(%ld): %s", pActuator->ioError().what,
                          inBlockPtr->getBlockId(),
                          pActuator->ioError().err == coder_io::IO_BUF_FULL ? "output buffer too small"
                                                                            : "input stream exhausted");
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

            /* 首个完成的块放行其余 coder 线程（通常是第 0 块）。 */
            if (firstCoderNotify(pActuator->getNotifyFlag())) {
                coderStartCond.notify_all();
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

bool PbgzEngine::firstCoderNotify(bool flag) {
    std::lock_guard<std::mutex> lock(coderStartMutex);
    if (flag && !coderStartSync) {
        coderStartSync = true;
        return true;
    }
    return false;
}

int32_t PbgzEngine::actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*) {
    return actuator->process();
}
