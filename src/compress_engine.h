/*
 * compress_engine.h - Header file for pbgz project
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

#include "codec_engine.h"
#include "actuator.h"
#include <mutex>
#include <condition_variable>
#include "pbgz_stat.h"
#include "block_wrapper.h"
#include "preprocess_info.h"
#include "codec_selector.h"
#include <atomic>

class CompressEngine : public CodecEngine {
public:
    CompressEngine(PbgzParameter& para) : CodecEngine(para) {
        indexBlockQueue = std::make_unique<BlockingQueueType>();
        freeIndexBlockQueue = std::make_unique<BlockingQueueType>();
    }

    virtual ~CompressEngine();

    virtual int32_t init();

    PbgzStat* getStats() { return stats.get(); }

    virtual const PreprocessInfo* getPreprocessInfo() override { return &preprocessInfo; }

    void initStatsBasedOnFileType(BlockType fileType);

    virtual void fileDecisionProc(RoughIOBlock* firstBlock) override;

protected:
    /*
     * The reader thread accumulates QUAL prior training samples block by block; once
     * the target block count (= concurrency) or the 45 MB cap is reached, the prior is
     * trained and published in one pass.
     */
    virtual void pretrainBlockProc(RoughIOBlock* blockPtr) override;

    /* Finalize any incomplete pretraining when the read loop ends (EOF/error). */
    virtual void readLoopPostProc() override;

    /* coder thread startup latch: starts pulling and compressing blocks only after being notified that the prior is published (or decided absent). */
    virtual void workStartBarrier() override;

    /* Train on the accumulated QUAL, publish the prior auxiliary block, set DONE, and notify the waiting coder threads. */
    void finalizePretrain();
    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual int32_t startEnginePreProc() override;

    virtual int32_t startEnginePostProc() override;

    virtual  Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr);


    /*
     * Emit the trained QUAL prior as a QUAL_PRIOR auxiliary block and record the
     * absolute offset of its container header.
     *
     * The call site is finalizePretrain (after cross-block pretraining completes), not
     * file end: because data blocks reference the prior, the prior must be fully
     * compressed before any data block. Since coder threads are held at the latch
     * until the prior is published, it is safe to write the prior to disk before
     * releasing them — even block 0 can use the prior, with no special case needed for
     * it.
     *
     * The actual write to disk is done by the writer thread via emitSyncAuxBlock; this
     * function only builds the block and records the metadata.
     */
    void emitQualPrior();

    virtual void writeFilePostProc(BlockWriter* blockWriter) override;

    virtual void printTailInfo(Timer& costTimer) override {
        PbgzManager::getInstance().printTailInfo(costTimer, true);
        if (parameter.verbose) {
            double secs = costTimer.elapsedSeconds();
            int64_t readLen = PbgzManager::getInstance().getTotalReadLen();
            int64_t writeLen = PbgzManager::getInstance().getTotalWriteLen();
            if (secs > 0) {
                fprintf(stderr, "Compress speed: %.2f MB/s in (%.2f MB/s out), %.3f s\n",
                        (readLen / (1024.0 * 1024.0)) / secs,
                        (writeLen / (1024.0 * 1024.0)) / secs,
                        secs);
            }
        }
        if (parameter.showStat && stats) {
            stats->printStats();
        }
    }

    virtual void setDataBlockPosition(uint32_t blockId) override;

    virtual int32_t prepareFileMeta() override;

    virtual int32_t startWorkPreProc() override;

    virtual Reference* getReference() override { return pRefGene; }

private:
    bool initReference();

    int64_t packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true);

    uint32_t calcPackRefeBlockSize();

    /* Whether the reference genome is packed with the file header: true when output goes to a pipe and an external reference is not requested. */
    bool isPackRefeInHeader() const {
        return !parameter.isUnpackRef && parameter.outputFile == STDOUT;
    }

private:
    std::map<uint32_t, std::vector<int64_t>> blockRefePos;
    std::map<int64_t, uint32_t> blockRefeIndex;
    std::unique_ptr<BlockingQueueType> indexBlockQueue;
    std::unique_ptr<BlockingQueueType> freeIndexBlockQueue;
    std::unique_ptr<PbgzStat> stats;
    bool statsInitialized = false;

    PreprocessInfo preprocessInfo;

    /*
     * Absolute file offset of the QUAL prior block's container header; -1 means this
     * compression has no usable prior. It is written by the reader thread before the
     * first block is dispatched and then read by all worker threads. The queue's
     * enqueue/dequeue already establish happens-before; the atomic is used only to
     * make this publication relationship self-documenting at the type level.
     */
    std::atomic<int64_t> qualPriorOffset{-1};

    /* The trained prior snapshot; copied once when building the block, then shared read-only so all worker threads can use it with zero copies. */
    AuxPayloadPtr qualPriorBlob;

    /*
     * Cross-block pretraining state (accessed only by the reader thread): the
     * accumulator, the count of blocks fed, the target block count, and whether
     * pretraining is active. The target block count is the concurrency; the 45 MB cap
     * is enforced by the accumulator itself.
     */
    CodecSelector::QualPriorAccum qualPriorAccum;
    uint32_t pretrainBlockCount = 0;
    uint32_t pretrainTarget = 0;
    bool pretraining = false;

    /*
     * Publication latch for coder threads: set and signaled by the reader thread after
     * it trains and publishes the prior; coder threads wait in workStartBarrier. For
     * files without a prior it is set immediately in fileDecisionProc, so no non-prior
     * path is slowed down; for empty files it is set as a fallback in
     * readLoopPostProc.
     */
    std::mutex priorMutex;
    std::condition_variable priorCv;
    bool priorSettled = false;

    void signalPriorSettled()
    {
        std::lock_guard<std::mutex> lk(priorMutex);
        priorSettled = true;
        priorCv.notify_all();
    }

public:
    int64_t getQualPriorAddress() const override {
        return qualPriorOffset.load(std::memory_order_acquire);
    }

    AuxPayloadPtr getQualPrior(int64_t /*packageIndex*/) override {
        /* Written before the first block is dispatched, then read-only; the acquire load of the address atomic also establishes its visibility. */
        return (getQualPriorAddress() < 0) ? AuxPayloadPtr() : qualPriorBlob;
    }
};
