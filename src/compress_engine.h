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

    virtual void writeFilePostProc(BlockWriter* blockWriter) override ;

    /*
     * 把首块训练出的 QUAL 先验写成一个 QUAL_PRIOR 块，并记下它的文件级偏移。
     *
     * 走的是 refe.offset 那个已经验证过的形状：解压侧靠偏移 seek，不靠块号定位。
     * 之所以必须由写线程在数据块全部落盘之后调用，是因为先验只有读完第 0 块才存在，
     * 没法像参考基因组那样在开工前就打进文件头；而写线程此刻已是单线程，
     * getCurrentPos() 取到的就是这个块的真实起点，不存在竞争，也无需跨线程回填。
     */
    void packQualPrior(BlockWriter* blockWriter);

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
};
