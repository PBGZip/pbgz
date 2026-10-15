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
#include "profile_stats.h"


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
    /* Input capacity is threadNum + 1, so the reader can finish pre-training before the queue fills. */
    freeInputPool->setCapility(parameter.threadNum + 1);
    inputDataPool->setCapility(parameter.threadNum + 1);
    freeOutputPool->setCapility(parameter.threadNum << 1);
    outputDataPool->setCapility(parameter.threadNum << 1);

    uint32_t blockBufferSize = getBlockSize();
    // First push empty blocks to free queue
    /* Input blocks start at 1MB and grow on demand; output blocks follow -l. */
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
            /* BAM is a BGZF stream, but must bypass transparent gz decompression so the block type stays BAM. */
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

    /* Register the output file for cleanup on abnormal exit. */
    if (parameter.outputFile != STDOUT && parameter.outputFile != "-") {
        PbgzManager::getInstance().addOutputFile(parameter.outputFile);
    }

    return 0;
}

int32_t PbgzEngine::start() {
    printHeadInfo();
    Timer costTimer(true);
    const std::chrono::steady_clock::time_point profWall0 = std::chrono::steady_clock::now();

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

    /* Flush and check the output here; closeIO is idempotent. */
    if (ioWriter != nullptr) {
        ioWriter->closeIO();
        int32_t writeErr = ioWriter->getWriteError();
        if (writeErr != 0) {
            LOG_ERROR("Flush output failed: %s", strerror(writeErr));
            fprintf(stderr, "Error: flushing output failed (%s), output is incomplete.\n", strerror(writeErr));
            return -1;
        }
    }

    /* Fail the run if any block failed, so a truncated output is never reported as success. */
    if (taskFailed.load()) {
        LOG_ERROR("Some blocks failed during processing; result is incomplete.");
        fprintf(stderr, "Error: some blocks failed, output is incomplete.\n");
        return -1;
    }

    printTailInfo(costTimer);
    
    {
        const double wallSec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - profWall0).count();
        pbgzprof::dump(stderr, wallSec);
    }

    return 0;
}

int64_t PbgzEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    RoughIOBlock* blockPtr = nullptr;
    {
        PBGZ_PROF_SCOPE(pbgzprof::READ_QUEUE_WAIT);
        blockPtr = freeInputPool->get();
    }
    if (blockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }
    blockPtr->reset();

    readBlockPreProc(blockReader);

    int64_t ret = 0;
    {
        PBGZ_PROF_SCOPE(pbgzprof::READ_BLOCK);
        ret = blockReader->readBlock(blockPtr, fileType);
    }
    if (ret <= 0) {
        freeInputPool->push(blockPtr);
        return ret;
    }

    /* An auxiliary block is not part of the data stream and is handled here for every engine. */
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

    /* The file-level decision must happen before the block is enqueued; SAM header blocks are skipped. */
    if (!fileDecisionInvoked && blockReader->blockHasData(blockPtr)) {
        fileDecisionInvoked = true;
        PBGZ_PROF_SCOPE(pbgzprof::READ_DECISION);
        fileDecisionProc(blockPtr);
    }

    /* Per-block accumulation (e.g. cross-block pre-training) runs before enqueueing. */
    {
        PBGZ_PROF_SCOPE(pbgzprof::READ_PRETRAIN);
        pretrainBlockProc(blockPtr);
    }

    /* Input-only work (block checksum) is done here, not on a worker thread. */
    preDispatchBlock(blockPtr);

    updateInputStatics(blockPtr);
    dispatchInputBlock(blockPtr);
    if (ret > 0) {
        blockCount++;
    }

    return ret;
}

/* Returns the last readOneBlock result: > 0 or -2 continue, 0 is EOF, -1 is an error. */
int64_t PbgzEngine::readBlocks(BlockReader* blockReader) {
    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        ret = readOneBlock(blockReader, fileType);
    } while(ret > 0 || ret == -2);

    /* Finish incomplete preprocessing before pushing the termination markers. */
    readLoopPostProc();
    return ret;
}

int32_t PbgzEngine::startReadTask() {
    /* Deliberately unnamed: this runs on the main thread, whose comm names the process. */
    PBGZ_PROF_SCOPE(pbgzprof::READ_TOTAL);
    BlockReader* blockReader = createBlockReader();
    if (blockReader == nullptr) {
        return -1;
    }
    engineReaderReady(blockReader);
    int64_t readRet = readBlocks(blockReader);
    /* Drain reader-side staging before the end markers. */
    renderStageDrain();
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
        PBGZ_PROF_SCOPE(pbgzprof::WRITE_TOTAL);
        BlockWriter* blockWriter = createBlockWriter();
        if (blockWriter == nullptr) {
            PbgzManager::getInstance().exitProc(-1, "Inner error");
            return -1;
        }

        while (true) {
            RoughIOBlock* outBlockPtr = nullptr;
            {
                PBGZ_PROF_SCOPE(pbgzprof::WRITE_WAIT);
                /* Take blocks as they arrive; reordering happens below. */
                outBlockPtr = outputDataPool->get();
            }
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
                /* Auxiliary blocks are position-addressed and bypass the reordering and blockId2Write. */
                completeSyncAuxBlock(blockWriter, outBlockPtr);
            } else {
                writeBlockPreProc(blockWriter, outBlockPtr);
                outputSortedCache.push_back(outBlockPtr);
                { PBGZ_PROF_SCOPE(pbgzprof::WRITE_SORT);
                /* Strict weak ordering: false for equal ids. */
                outputSortedCache.sort([](const RoughIOBlock* p1, RoughIOBlock* p2) {
                    return p1->getBlockId() < p2->getBlockId();
                }); }
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
    PBGZ_PROF_SCOPE(pbgzprof::WRITE_BLOCK);
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
    /* Capture the write position before writeBlock: it is the container-header offset. */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    const int64_t offset = (fileWriter != nullptr) ? (int64_t)fileWriter->getCurrentPos() : -1;

    blockWriter->writeBlock(block);
    updateOutputStatics(block);
    /* Return to the free pool before waking the emitter. */
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
        PBGZ_PROF_SCOPE(pbgzprof::CODER_TOTAL);

        LOG_INFO("Coder task (%d) begin to running!", id);

        /* Coder threads wait for the prior to be released before pulling blocks. */
        {
            PBGZ_PROF_SCOPE(pbgzprof::CODER_BARRIER);
            workStartBarrier();
            if (id > 0) {
                std::unique_lock<std::mutex> lock(coderStartMutex);
                coderStartCond.wait(lock, [this] { return coderStartSync; });
            }
        }

        /* First-block serialization: id==0 starts; others wait until block 0 is processed. */

        while (true) {
            RoughIOBlock* inBlockPtr = nullptr;
            {
                PBGZ_PROF_SCOPE(pbgzprof::CODER_WAIT_INPUT);
                inBlockPtr = inputDataPool->get();
            }
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
                /* Release the others even if block 0 failed, to avoid a deadlock. */
                if (firstCoderNotify(true)) {
                    coderStartCond.notify_all();
                }
                break;
            }
            RoughIOBlock* outBlockPtr = nullptr;
            {
                PBGZ_PROF_SCOPE(pbgzprof::CODER_WAIT_OUTPUT);
                outBlockPtr = freeOutputPool->get();
            }
            if (outBlockPtr == nullptr) {
                LOG_ERROR("Get free output block failed.");
                freeInputPool->push(inBlockPtr);
                break;
            }
            outBlockPtr->reset();
            outBlockPtr->setBlockId(inBlockPtr->getBlockId());

            Actuator* pActuator = nullptr;
            {
                PBGZ_PROF_SCOPE(pbgzprof::CODER_CREATE);
                pActuator = createActuator(inBlockPtr, outBlockPtr);
            }
            if (pActuator == nullptr) {
                /* A failed actuator must not exit the loop; do the same cleanup as actuatorProc failure. */
                LOG_ERROR("Create actuator failed for block(%ld)", inBlockPtr->getBlockId());
                taskFailed.store(true);
                /* Capture the id before returning the block; push() lets the reader rewrite it. */
                const int64_t failedBlockId = inBlockPtr->getBlockId();
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(failedBlockId);
                outputDataPool->push(outBlockPtr);
                continue;
            }

            /* Convert coder exceptions into return values for uniform handling. */
            int32_t ret = 0;
            try {
                PBGZ_PROF_SCOPE(pbgzprof::CODER_COMPRESS);
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

            /* Check coder_io out-of-bounds once, covering every stream and actuator. */
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
                /* Capture the id before push(); the reader rewrites it after. */
                const int64_t failedBlockId = inBlockPtr->getBlockId();
                freeInputPool->push(inBlockPtr);
                outBlockPtr->reset();
                outBlockPtr->setBlockId(failedBlockId);
                // Push a zero-length block with the correct id to keep the writer advancing
                outputDataPool->push(outBlockPtr);
                MemoryUtil::safeDeleteClass(pActuator);
                continue;
            }

            outBlockPtr->setBlockType(inBlockPtr->getBlockType());

            /* Release the other coder threads on first completion. */
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
