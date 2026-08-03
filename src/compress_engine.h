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
#include "pbgz_stat.h"
#include "block_wrapper.h"
#include "preprocess_info.h"
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
    
    virtual int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType) override;
    
protected:
    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual int32_t startEnginePreProc() override;

    virtual int32_t startEnginePostProc() override;

    virtual  Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr);


    /*
     * 把首块训练出的 QUAL 先验发射成一个 QUAL_PRIOR 辅助块，并记下其容器头绝对偏移。
     *
     * 调用点在串行首块窗口内（runFilePreprocessOnce 之后、第 0 块编码之前），
     * 而不是文件收尾：数据块的分片头要写下先验的绝对地址，地址就必须先于一切数据块存在。
     * 放在这里还有个直接好处——连第 0 块自己都能用上先验，不必为它开特例。
     *
     * 实际落盘由写线程经 emitSyncAuxBlock 完成，本函数只负责构块与登记元信息。
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

    void runFilePreprocessOnce(RoughIOBlock* inBlockPtr);

    int64_t packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true);

    uint32_t calcPackRefeBlockSize();

    /* 参考基因组是否随文件头一起打包：输出到管道且未要求外挂参考时为真。 */
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
     * QUAL 先验块容器头的绝对文件偏移，-1 表示本次压缩没有可用先验。
     * 由 0 号工作线程在串行窗口内写入，随后被所有并行工作线程读取，故用原子量发布：
     * conditionVar 的 notify_all 只保证唤醒，不足以为非原子共享量建立可见性。
     */
    std::atomic<int64_t> qualPriorOffset{-1};

    /* 训练出的先验快照，构块时做一次拷贝后只读共享，供全部工作线程零拷贝取用。 */
    AuxPayloadPtr qualPriorBlob;

public:
    int64_t getQualPriorAddress() const override {
        return qualPriorOffset.load(std::memory_order_acquire);
    }

    AuxPayloadPtr getQualPrior(int32_t /*packageIndex*/) override {
        /* 在串行首块窗口内写入，随后只读；地址原子量的 acquire 同时为它建立可见性。 */
        return (getQualPriorAddress() < 0) ? AuxPayloadPtr() : qualPriorBlob;
    }
};
