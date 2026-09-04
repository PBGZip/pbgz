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

    /*
     * Prepare the base file metadata; must be called before startWriteTask.
     *
     * Why this hook exists separately: as soon as the writer thread starts, it runs
     * createBlockWriter, which immediately calls writeBaseFileMeta to persist
     * baseFileMeta. If the main thread is still filling baseFileMeta at that point, the
     * two sides race on the same JSON - whoever finishes first decides whether that
     * portion ends up in the metadata written to the file header.
     *
     * The practical consequence of this race is that compression becomes
     * non-reproducible: compressing the same input twice may leave the file-header
     * metadata differing by a whole JSON member, shifting everything that follows and
     * producing entirely different output bytes. Both results still decompress
     * correctly (the decompression side tolerates missing members), which is why the
     * problem went unnoticed for a long time, but it breaks the basic property that
     * "the same input always yields the same output" and makes any byte-level
     * regression comparison impossible.
     *
     * Convention therefore: anything that must be written into the base file metadata
     * has to be filled in within this hook; work that requires the writer thread to be
     * ready (e.g. writing the reference block downstream) stays in startWorkPreProc.
     */
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

    /*
     * Callback after a data block is read and before it is enqueued, for engine-side
     * per-block accumulation and similar preprocessing.
     *
     * The base class is a no-op by default. The compression engine uses it to
     * accumulate QUAL prior training samples across blocks (first-block decision plus
     * appending on subsequent blocks); once the target block count or cap is reached it
     * trains and releases the prior.
     */
    virtual void pretrainBlockProc(RoughIOBlock* /*blockPtr*/) {}

    /*
     * Called once after the read loop ends (EOF or read error), for the engine to wrap
     * up any incomplete preprocessing. The base class is a no-op by default.
     */
    virtual void readLoopPostProc() {}

    /*     * Synchronization point before worker threads start pulling blocks to compress.
     * The base class is a no-op by default.
     *
     * The compression engine uses it to enforce "coder threads must not start before the
     * prior is released": the reader thread releases the cross-block trained prior
     * (persisting the auxiliary block and setting the flag) and then notifies, which is
     * when the worker threads are released. Other engines have no such requirement.
     */
    virtual void workStartBarrier() {}

    /*     * Disposition of a block after it is read and before it goes into inputDataPool,
     * asked uniformly by readOneBlock.
     *
     * This is a return value rather than letting each engine override the whole
     * readOneBlock because the skeleton of the read loop (acquire free block / read /
     * first-block decision / enqueue / count) is identical for all engines; the only
     * step that genuinely differs is "is this block mine". Previously each of the three
     * engines duplicated the skeleton, so any change to the skeleton missed two copies.
     */
    enum class BlockIntake {
        DISPATCH,   /* Normal data block; enqueue it for the worker threads */
        SKIP,       /* Data stream that does not belong to this engine; already consumed in place, the read loop continues */
        ABORT       /* Input unusable; terminate the read loop */
    };

    virtual BlockIntake intakeBlock(BlockReader* /*blockReader*/, RoughIOBlock* /*blockPtr*/) {
        return BlockIntake::DISPATCH;
    }

    /* A moment before a block is read, for capturing quantities that are only valid before the read, such as "the block's starting position in the source file". */
    virtual void readBlockPreProc(BlockReader* /*blockReader*/) { }

    /*     * File-level decision: the reader thread calls it once before handing the first
     * data block to inputDataPool.
     *
     * Placement is the whole point of this hook. Decisions such as codec selection and
     * prior training are made exactly once per file and must precede the processing of
     * any block. Placing it "before enqueue" guarantees this ordering by the data-flow
     * position itself: worker threads can only take blocks from the queue, and without
     * a block there is nothing to compress. The queue's enqueue/dequeue provide
     * release/acquire semantics, so visibility of the decision result is covered too;
     * no waiting between worker threads or extra synchronization flags are needed.
     *
     * A property gained along the way: the decision always happens on the reader thread
     * and is always based on block 0, so the sample no longer drifts with scheduling;
     * and when a synchronous auxiliary block (such as the QUAL prior) is emitted the
     * writer thread has not yet received any data block, turning "the auxiliary block
     * physically precedes all data blocks" from a timing coincidence into a positional fact.
     */
    virtual void fileDecisionProc(RoughIOBlock* /*firstBlock*/) { }

    virtual void updateInputStatics(RoughIOBlock*) { }

    virtual void updateOutputStatics(RoughIOBlock*) { }

    virtual void writeBlockPreProc(BlockWriter*, RoughIOBlock*) { }

    virtual void writeBlockPostProc() { }

    virtual void writeFilePreProc() { }

    virtual void writeFilePostProc(BlockWriter*) { }

    virtual void writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr);

    /* Writer-thread side: persist a synchronous auxiliary block and hand its container-header offset back to the emitter. */
    void completeSyncAuxBlock(BlockWriter* blockWriter, RoughIOBlock* block);

    /*
     * When an auxiliary block passes by, ask each consumer in turn whether it claims it.
     * The return value is only used for logging: an auxiliary block claimed by no one is
     * silently skipped, which is the forward-compatible behavior an old version needs
     * when reading a newer format and is not an error.
     */
    bool offerAuxBlock(RoughIOBlock* blockPtr, int64_t packageIndex);

    /*
     * Synchronously emit an auxiliary block: push it to the writer thread, block until
     * it is persisted, and return the absolute file offset of the block's **container
     * header** (i.e. the position from which the block-read manager can parse the block
     * JSON, not the start of the block's payload).
     *
     * Auxiliary blocks are position-addressed: they do not occupy a data-block blockId
     * and do not participate in the writer thread's reordering; the writer writes them
     * on sight. The write must be performed by the writer thread itself because only it
     * holds the BlockWriter and exclusively owns the file write pointer - a caller on
     * any other thread reading getCurrentPos() would get a racy value.
     *
     * May only be called from within fileDecisionProc: at that point the reader thread
     * has not dispatched any data block, the writer thread is necessarily idle, and at
     * most one emission is in flight globally, so a single slot below is enough to carry
     * the result with no per-block mapping. After return the block has already been
     * returned to freeOutputPool by the writer thread; the caller must not touch it.
     */
    int64_t emitSyncAuxBlock(RoughIOBlock* block);

    void registerAuxConsumer(AuxBlockConsumer* consumer) {
        if (consumer != nullptr) {
            auxConsumers.push_back(consumer);
        }
    }

public:
    std::unique_ptr<BlockingQueueType> freeInputPool;   // Free queue, file reading tasks get blocks from here to read data
    std::unique_ptr<BlockingQueueType> inputDataPool;   // Data reading tasks write completed data to this queue, compression/decompression tasks get data from this queue
    std::unique_ptr<BlockingQueueType> freeOutputPool;  // Compression/decompression tasks get blocks from this queue to write processed data, output tasks write free blocks to this queue after processing
    std::unique_ptr<BlockingQueueType> outputDataPool;  // After compression/decompression is completed, write to this queue, output tasks get data from this queue
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;
    const PbgzParameter& getParameter() const { return parameter; }
    virtual Reference* getReference() { return nullptr; }

    /*
     * Upper bound of the block size written to the file header when compressing. On the
     * decompression side createBlockReader reads it back from baseFileMeta, and the
     * actuator uses it to pre-allocate the output buffer with ensureCapacity(block_size*2)
     * - the primary defense against out-of-bounds access in every field. 0 means it has
     * not been read yet (older files do not write it, or createBlockReader has not run);
     * the actuator then falls back to the default getBlockSize(). The coder_io putc
     * check and the decode error-return chain serve as a backstop.
     */
    uint32_t fileBlockSize = 0;
    uint32_t getFileBlockSize() const { return fileBlockSize; }

    /* File preprocessing result (codec pre-selection). Only the compression
       engine populates this; other engines return nullptr. */
    virtual const PreprocessInfo* getPreprocessInfo() { return nullptr; }

    /*
     * Fetch the QUAL prior model snapshot. The two sides obtain it from different
     * sources but through the same interface: the compression side returns the one
     * trained during first-block preprocessing; the decompression side retrieves the
     * claimed one by the auxiliary block's absolute address. Returning empty means no
     * prior is available this run, and callers must fall back to the fixed
     * initial-value model.
     *
     * The decompression side does not seek back to re-read: the prior block now
     * physically precedes all data blocks, so a sequential stream is guaranteed to pass
     * it before reaching any data block; the region-query path prefetches it at engine
     * initialization according to the file metadata. Both paths guarantee the cache is
     * in place by the time a worker thread uses it.
     */
    virtual AuxPayloadPtr getQualPrior(int64_t /*packageIndex*/) { return AuxPayloadPtr(); }

    /* Absolute offset of the prior block's container header; the compression side writes it into the block meta and the decompression side uses it to look up the block. -1 means no prior. */
    virtual int64_t getQualPriorAddress() const { return -1; }

    /* File-level POS delta prior (parsed from the file meta by the decoding
     * engine; empty when the file has none). Unlike the QUAL prior it is a
     * single file-wide table, so no package index is involved. */
    virtual AuxPayloadPtr getPosPrior() { return AuxPayloadPtr(); }


    std::vector<AuxBlockConsumer*> auxConsumers;

    std::vector<std::thread> workThreads;
    std::thread writeThread;
    int64_t blockId2Write;

    uint32_t blockCount;

    /*
     * Final funnel for block-level failures: set when actuatorProc returns non-zero on
     * a worker thread; engine cleanup (after startEnginePostProc) checks it uniformly
     * and returns an error, so the process exits non-zero - a block failure is no
     * longer dismissed with a warning and a zero-length block. atomic because worker
     * threads set it concurrently.
     */
    std::atomic<bool> taskFailed{false};

    /* Whether the file-level decision (codec selection, etc.) has been executed; a SAM header block carries no data, so the decision is deferred to the first data block. */
    bool fileDecisionInvoked = false;

    /* Handoff slot for emitting synchronous auxiliary blocks, shared by the emitter and the writer thread and guarded by auxEmitMutex. */
    mutable std::mutex auxEmitMutex;
    mutable std::condition_variable auxEmitCond;
    int64_t auxEmitOffset = -1;
    bool auxEmitDone = false;

    /*
     * First-block serialization: the coder thread with id==0 finishes processing block 0
     * before releasing the other threads, so concurrent preAnalysis cannot misclassify a
     * block as binary while SamInfo is not yet populated. This mechanism originally
     * existed in the perf branch and was accidentally removed when the prior latch
     * (workStartBarrier) was integrated; restored here. The release happens only once.
     */
    mutable std::mutex coderStartMutex;
    std::condition_variable coderStartCond;
    bool coderStartSync = false;

    /* Returns true the first time it is called and records it; always returns false afterwards, working with notify_all to release exactly once. */
    bool firstCoderNotify(bool flag);
};
