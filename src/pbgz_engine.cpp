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
     * Input-side capacity is threadNum + 1: once the SAM header is isolated as its
     * own block it no longer contributes to the QUAL prior pre-training count. If the
     * capacity were exactly threadNum, the reader thread would block on the free-block
     * queue once it filled the queue and could never accumulate the target number of
     * pre-training blocks; finalizePretrain would never fire and the coder threads
     * would wait forever for priorSettled -> deadlock. The extra slot lets the reader
     * get one block ahead of the coders, so pre-training finishes before the queue fills.
     */
    freeInputPool->setCapility(parameter.threadNum + 1);
    inputDataPool->setCapility(parameter.threadNum + 1);
    freeOutputPool->setCapility(parameter.threadNum << 1);
    outputDataPool->setCapility(parameter.threadNum << 1);

    uint32_t blockBufferSize = getBlockSize();
    // First push empty blocks to free queue
    /*
     * Input blocks are initially allocated with a fixed 1MB; reads target the size
     * chosen by -l (getBlockSize()), and when capacity is insufficient the reader's
     * ensureCapacity reallocs on demand. This avoids allocating large memory up front
     * based on -l (e.g. 512MB x N at -l 9). Output blocks and the header block_size
     * still follow -l.
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
             * Standard BAM is a BGZF (gzip) stream that is indistinguishable from a .gz
             * file byte by byte. But its gzip is the inner layer of the format and must
             * not go through "transparent gz decompression": transparent decompression
             * uses the Gz/FastGz readers, which cannot obtain the total length of the
             * input file (the dynamic_cast<FileReader*> in fileDecisionProc fails), so
             * the "is it worth writing" criterion of the QUAL prior degrades and a BAM
             * input ends up writing a needless prior auxiliary block. Therefore BAM is
             * always handed to FileReader + BamGzBlockReader (the latter inflates
             * internally), the block type remains BAM, and compression behavior is
             * aligned with SAM input.
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

    /* Register the output file so PbgzManager can clean it up on abnormal exit (e.g. missing reference genome). */
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
     * Writing to disk must be completed and checked here. The msync inside closeIO is
     * the point where mmap output is actually flushed back to disk; it used to run only
     * in the destructor, later than start() returning, leaving no way to report failure.
     * Calling it again is harmless: closeIO nulls out mappedAddress and fd.
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
     * Final handling of block-level errors: if any block fails to compress/decompress,
     * the overall result cannot be trusted. Zero-length blocks are only a stopgap to
     * keep the pipeline from deadlocking; the process must never exit successfully and
     * leave a truncated file behind.
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
     * By definition an auxiliary block is not part of the data stream and behaves the
     * same for every engine, so this rule can only live here. An auxiliary block that
     * no one claims is not an error: it is the forward-compatible behavior an old
     * version needs when reading a newer format.
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
     * The file-level decision must happen before push (enqueueing means some worker
     * thread may start processing immediately). A SAM header block (only @ lines, no
     * data lines) does not trigger the decision; it is deferred to the first data block,
     * otherwise codec selection would collect no field samples and fall back to the
     * default codec. It must be called before push; behavior of other engines is unchanged.
     */
    if (!fileDecisionInvoked && blockReader->blockHasData(blockPtr)) {
        fileDecisionInvoked = true;
        fileDecisionProc(blockPtr);
    }

    /* Per-block accumulation work such as cross-block pre-training (used by the compression engine) is also executed before enqueueing. */
    pretrainBlockProc(blockPtr);

    updateInputStatics(blockPtr);
    inputDataPool->push(blockPtr);
    if (ret > 0) {
        blockCount++;
    }

    return ret;
}

/*
 * Returns the result of the last readOneBlock: > 0 or -2 never reach here (the loop
 * continues in those cases), 0 is a clean EOF, and -1 is a read error. Callers must
 * use this to distinguish "read to the end" from "read corrupted data" - a truncated
 * or damaged input must never be silently treated as a normal EOF.
 */
int64_t PbgzEngine::readBlocks(BlockReader* blockReader) {
    BlockType fileType = TYPE_UNKNOW;
    int64_t ret = 0;
    do {
        ret = readOneBlock(blockReader, fileType);
    } while(ret > 0 || ret == -2);

    /* Finish any incomplete preprocessing after EOF/error (e.g. cross-block prior training); this must happen before pushing the termination markers. */
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
                 * Auxiliary blocks are position-addressed and must be handled ahead of
                 * writeBlockPreProc / the reordering / blockId2Write: they do not belong
                 * to the data-block id sequence, and once they slip into the sorting
                 * logic below they would displace a data block and permanently misalign
                 * the write pointer.
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
     * The position must be captured immediately before writeBlock: the first byte
     * written by writeBlock is the container magic number, so the write pointer at this
     * instant is exactly the offset that the decompression side's readBlock needs.
     */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    const int64_t offset = (fileWriter != nullptr) ? (int64_t)fileWriter->getCurrentPos() : -1;

    blockWriter->writeBlock(block);
    updateOutputStatics(block);
    /* Return to the free pool must precede the wake-up: the emitter, once woken, may immediately request another free block. */
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

        /* Prior-release latch of the compression engine: coder threads wait for notification before the prior is released, and only then start pulling blocks to compress. */
        workStartBarrier();

        /*
         * First-block serialization (a perf-branch mechanism, see firstCoderNotify):
         * the thread with id==0 starts directly; the others wait until block 0 has been
         * processed (preAnalysis populates shared state such as SamInfo) before being
         * released, so that concurrent preAnalysis cannot downgrade a block to BINARY
         * while the chromosome table is not yet ready.
         *
         * The wait must use a predicate: if no actuator can be built for block 0 (e.g. an
         * unsupported block type), processing fails immediately and the release
         * notification may be emitted before the other threads start waiting, thus being
         * missed - a bare wait without a predicate would then block forever and deadlock
         * the whole pipeline. The predicate re-checks coderStartSync while holding the
         * lock; whoever sets it first wins.
         */
        if (id > 0) {
            std::unique_lock<std::mutex> lock(coderStartMutex);
            coderStartCond.wait(lock, [this] { return coderStartSync; });
        }

        while (true) {
            RoughIOBlock* inBlockPtr = inputDataPool->get();
            if (inBlockPtr == nullptr) {  // Got null pointer, indicating end marker
                /* Even when block 0 was not completed (read error / all-failed), release the others, otherwise they would block forever. */
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
                 * Failing to build an actuator for a single block must never exit the
                 * loop: if a worker thread ends early, the remaining threads have to
                 * drain the entire input queue by themselves, and if all of them exit,
                 * the reader thread blocks forever on the bounded queue. This path takes
                 * exactly the same cleanup as an actuatorProc failure below - return the
                 * input block, push a zero-length block with the same id so the writer
                 * can keep advancing, then continue with the next block.
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
             * The coder layer's check_exit/coder_exit now throw coder_exception (they
             * used to kill the process directly with _Exit inside the library, leaving
             * no log). This is their final handling point: catch and convert them into
             * return values identical to any other failure, funneling into taskFailed.
             * catch(...) is also included: an exception escaping a worker thread would
             * invoke std::terminate directly, which is harder to debug than a dropped block.
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
             * The single place that checks out-of-bounds errors. A coder_io overflow
             * only sets a flag without aborting, and the actuator does not necessarily
             * turn it into a return value - previously SAM was checked in 12 places
             * while FASTQ and the index were never checked, so data corrupted by an
             * out-of-bounds write was silently written out as "success". Now all
             * coder_io errors are funneled onto the actuator, and asking once here
             * covers every stream and every actuator, so newly added streams cannot
             * be missed either.
             *
             * This sits before the ret check rather than merged into it: this extra
             * question is only needed when ret is already 0; when ret is non-zero the
             * failure cause is more specific and must not be overwritten.
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

            /* The first block that completes releases the other coder threads (usually block 0). */
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
