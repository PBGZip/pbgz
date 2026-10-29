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
#include <vector>
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

    virtual int32_t init() override;

    PbgzStat* getStats() { return stats.get(); }

    virtual const PreprocessInfo* getPreprocessInfo() override { return &preprocessInfo; }

    virtual PreprocessInfo* getPreprocessInfoMut() override { return &preprocessInfo; }

    void initStatsBasedOnFileType(BlockType fileType);

    virtual void fileDecisionProc(RoughIOBlock* firstBlock) override;

    /* Computes the block MD5 on the reader thread (see PbgzEngine::readOneBlock). */
    virtual void preDispatchBlock(RoughIOBlock* blockPtr) override;

    /* Build the file-level OPTION layout plan on its own thread (see buildSamOptionPlan). */
    void startSamOptionPlan(RoughIOBlock* blockPtr);

    /* Archive-mode BAM: launch the SAM render (conversion) thread pool. */
    virtual void engineReaderReady(BlockReader* blockReader) override;

    /* Archive-mode BAM: divert raw-record blocks into the render stage. */
    virtual void dispatchInputBlock(RoughIOBlock* blockPtr) override;

    /* Archive-mode BAM: flush + join the render stage before the coder end markers. */
    virtual void renderStageDrain() override;

protected:

    /* Render worker: BAM records -> SAM text + block MD5 for one raw block (see bamRenderRawBlockToText). */
    void renderWorker();

    /* True while raw BAM blocks are routed to the render stage (set on the reader thread). */
    bool renderStageEnabled = false;
    uint32_t renderThreadNum = 0;
    std::unique_ptr<BlockingQueueType> renderInputPool;
    std::vector<std::thread> renderThreads;

    /* File-level OPTION layout plan builder; joined before the engine goes away. */
    std::thread optionPlanThread;
    bool optionPlanStarted = false;

    /* Reader thread accumulates QUAL prior samples block by block until the target block
       count or the 45 MB cap is reached. */
    virtual void pretrainBlockProc(RoughIOBlock* blockPtr) override;

    /* Finalize any incomplete pretraining when the read loop ends (EOF/error). */
    virtual void readLoopPostProc() override;

    /* Coder thread startup latch: waits until the prior is published or decided absent. */
    virtual void workStartBarrier() override;

    /* Train on the accumulated QUAL, publish the prior, and notify waiting coder threads. */
    void finalizePretrain();
    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual int32_t startEnginePreProc() override;

    virtual int32_t startEnginePostProc() override;

    virtual  Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;


    /* Build the QUAL_PRIOR auxiliary block and record its header offset. The writer
       thread performs the actual write via emitSyncAuxBlock. */
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

    int64_t packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true,
                          bool syncAux = false);

    uint32_t calcPackRefeBlockSize();

    /* Whether the reference genome is packed with the file header (stdout output, unpacked reference not requested). */
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

    /* Absolute offset of the QUAL prior block header; -1 means no prior. Written by the
       reader thread before the first block is dispatched, then read by all workers. */
    std::atomic<int64_t> qualPriorOffset{-1};

    /* Trained prior snapshot, shared read-only by all worker threads. */
    AuxPayloadPtr qualPriorBlob;

    /* Cross-block pretraining state (reader thread only); the 45 MB cap lives in the
       accumulator. */
    CodecSelector::QualPriorAccum qualPriorAccum;
    uint32_t pretrainBlockCount = 0;
    uint32_t pretrainTarget = 0;
    bool pretraining = false;

    /* Publication latch: signaled once the prior is published or decided absent; coder
       threads wait in workStartBarrier. */
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
        /* Written before the first block is dispatched, then read-only. */
        return (getQualPriorAddress() < 0) ? AuxPayloadPtr() : qualPriorBlob;
    }
};
