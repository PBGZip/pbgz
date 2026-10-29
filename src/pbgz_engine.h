/*
 * pbgz_engine.h - Head file for pbgz project
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

#pragma once

#include <thread>
#include <list>
#include <atomic>

#include "pbgz_types.h"
#include "blocking_queue.h"
#include "io_block.h"
#include "io_wrapper.h"
#include "block_wrapper.h"
#include "actuator.h"
#include "utils/timer.h"
#include "reference.h"
#include "aux_block_consumer.h"
#include "sam_option_plan.h"

struct PreprocessInfo;

using BlockingQueueType = BlockingQueue<RoughIOBlock*>;

class PbgzEngine {
public:
    PbgzEngine(const PbgzParameter&  para);

    virtual ~PbgzEngine();

    virtual int32_t init();

    virtual int32_t start();

protected:
    virtual void printHeadInfo() { };

    virtual void printTailInfo(Timer&) { }

    virtual uint32_t getBlockSize();

    virtual int32_t startEnginePreProc() { return 0; }

    /* Fill the base file metadata; must run before the writer thread starts. */
    virtual int32_t prepareFileMeta() { return 0; }

    virtual int32_t startWriteTask();

    virtual int32_t startWorkPreProc() { return 0; }

    virtual int32_t startWorkTask();

    virtual int32_t startReadPreProc() { return 0; }

    virtual int32_t startReadTask();

    virtual int32_t startEnginePostProc() { return 0; }

    virtual BlockReader* createBlockReader() = 0;

    virtual BlockWriter* createBlockWriter() = 0;

    virtual void releaseBlockReader(BlockReader* &blockReader) = 0;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) = 0;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) = 0;

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*);

    virtual int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType);

    virtual int64_t readBlocks(BlockReader* blockReader);

    /* Callback after a block is read and before it is enqueued; no-op by default. */
    virtual void pretrainBlockProc(RoughIOBlock* /*blockPtr*/) {}

    /* Called once after the read loop ends; no-op by default. */
    virtual void readLoopPostProc() {}

    /* Synchronization point before workers start pulling blocks; no-op by default. */
    virtual void workStartBarrier() {}

    /* Disposition of a block after it is read and before it enters inputDataPool. */
    enum class BlockIntake {
        DISPATCH,   /* Normal data block; enqueue it for the worker threads */
        SKIP,       /* Data stream that does not belong to this engine; already consumed in place, the read loop continues */
        ABORT       /* Input unusable; terminate the read loop */
    };

    virtual BlockIntake intakeBlock(BlockReader* /*blockReader*/, RoughIOBlock* /*blockPtr*/) {
        return BlockIntake::DISPATCH;
    }

    /* Called just before a block is read, for pre-read state capture. */
    virtual void readBlockPreProc(BlockReader* /*blockReader*/) { }

    /* File-level decision, called once before the first data block is dispatched. */
    virtual void fileDecisionProc(RoughIOBlock* /*firstBlock*/) { }

    /* Called after a block is classified but before it is pushed to the workers. */
    virtual void preDispatchBlock(RoughIOBlock* /*blockPtr*/) { }

    /* Called right after the block reader is created, before any data is read. */
    virtual void engineReaderReady(BlockReader* /*blockReader*/) { }

    /* Reader thread's final hand-off of a block to the workers. */
    virtual void dispatchInputBlock(RoughIOBlock* blockPtr) {
        inputDataPool->push(blockPtr);
    }

    /* Coder thread's hand-off of a finished block; override to insert a post-process stage. */
    virtual void dispatchOutputBlock(RoughIOBlock* blockPtr) {
        outputDataPool->push(blockPtr);
    }

    /* Start an optional post-processing stage between the coder workers and the writer. */
    virtual int32_t startOutputStage() { return 0; }

    /* Stop and join the optional post-processing stage; must be safe to call twice. */
    virtual void drainOutputStage() { }

    /* Called once after the read loop ends, before the coder end markers are pushed. */
    virtual void renderStageDrain() { }

    virtual void updateInputStatics(RoughIOBlock*) { }

    virtual void updateOutputStatics(RoughIOBlock*) { }

    virtual void writeBlockPreProc(BlockWriter*, RoughIOBlock*) { }

    virtual void writeBlockPostProc() { }

    virtual void writeFilePreProc() { }

    virtual void writeFilePostProc(BlockWriter*) { }

    virtual void writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr);

    /* Writer-thread side: persist a synchronous auxiliary block and return its offset. */
    void completeSyncAuxBlock(BlockWriter* blockWriter, RoughIOBlock* block);

    /* Ask each consumer whether it claims an auxiliary block; result only used for logging. */
    bool offerAuxBlock(RoughIOBlock* blockPtr, int64_t packageIndex);

    /* Synchronously emit an auxiliary block and return its container-header offset. */
    int64_t emitSyncAuxBlock(RoughIOBlock* block);

    void registerAuxConsumer(AuxBlockConsumer* consumer) {
        if (consumer != nullptr) {
            auxConsumers.push_back(consumer);
        }
    }

public:
    std::unique_ptr<BlockingQueueType> freeInputPool;   // Free blocks for the read task
    std::unique_ptr<BlockingQueueType> inputDataPool;   // Completed input blocks for compression/decompression tasks
    std::unique_ptr<BlockingQueueType> freeOutputPool;  // Free blocks for compression/decompression tasks
    std::unique_ptr<BlockingQueueType> outputDataPool;  // Processed blocks for the output task
    PbgzParameter parameter;
    /* SAM QNAME coder decision, shared by every block: -1 undecided, 0 = coder_qname, 1 = affix split. */
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;
    const PbgzParameter& getParameter() const { return parameter; }
    virtual Reference* getReference() { return nullptr; }

    /* Block size upper bound from the file header; 0 if not yet read. */
    uint32_t fileBlockSize = 0;
    uint32_t getFileBlockSize() const { return fileBlockSize; }

    /* File preprocessing result (codec pre-selection); null unless compression. */
    virtual const PreprocessInfo* getPreprocessInfo() { return nullptr; }

    /* Mutable view of getPreprocessInfo()'s object; null wherever that is null. */
    virtual PreprocessInfo* getPreprocessInfoMut() { return nullptr; }

    /* Fetch the QUAL prior model snapshot; empty if no prior is available. */
    virtual AuxPayloadPtr getQualPrior(int64_t /*packageIndex*/) { return AuxPayloadPtr(); }

    /* Absolute offset of the prior block's container header; -1 means none. */
    virtual int64_t getQualPriorAddress() const { return -1; }

    /* File-level POS delta prior, parsed from file meta; empty if none. */
    virtual AuxPayloadPtr getPosPrior() { return AuxPayloadPtr(); }

    /* File-level SAM OPTION layout plan (see buildSamOptionPlan). */
    SamOptionPlan& getSamOptionPlan() { return samOptionPlan; }


    std::vector<AuxBlockConsumer*> auxConsumers;

    std::vector<std::thread> workThreads;
    std::thread writeThread;
    int64_t blockId2Write;

    uint32_t blockCount;

    /* Set on block-level failure; checked during engine cleanup. */
    std::atomic<bool> taskFailed{false};

    /* Whether the file-level decision has been executed. */
    bool fileDecisionInvoked = false;

    /* Handoff slot for synchronous auxiliary block emission, guarded by auxEmitMutex. */
    mutable std::mutex auxEmitMutex;
    mutable std::condition_variable auxEmitCond;
    int64_t auxEmitOffset = -1;
    bool auxEmitDone = false;

    /* First-block serialization: coder id==0 finishes block 0 before releasing others. */
    mutable std::mutex coderStartMutex;
    std::condition_variable coderStartCond;
    bool coderStartSync = false;

    /* Returns true the first time, false afterwards; used to release once. */
    bool firstCoderNotify(bool flag);

    /* Written by the reader before the coder workers need it, read-only afterwards. */
    SamOptionPlan samOptionPlan;
};
