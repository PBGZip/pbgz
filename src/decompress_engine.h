/*
 * decompress_engine.h - Header file for pbgz project
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

#include <map>
#include <mutex>
#include <vector>

#include "codec_engine.h"
#include "actuator.h"


/*
 * QUAL 先验的认领者。读线程路过 QUAL_PRIOR 块时把它解出来存下，工作线程随后取用。
 *
 * blob 由读线程写、工作线程读，两者之间没有其它同步点，所以这里自带一把锁；
 * 读只发生在解码器构造时，不在热路径上。
 */
class QualPriorConsumer : public AuxBlockConsumer {
public:
    bool claim(RoughIOBlock* blockPtr, int32_t packageIndex) override;

    /*
     * 按包序号取该包的先验，没有则返回空。
     * 返回 shared_ptr 是为了让取用者把这份数据一直握到自己解完，
     * 不受认领者后续动作影响。
     */
    AuxPayloadPtr forPackage(int32_t packageIndex) const {
        std::lock_guard<std::mutex> lock(mutex);
        std::map<int32_t, AuxPayloadPtr>::const_iterator it = byPackage.find(packageIndex);
        return (it == byPackage.end()) ? AuxPayloadPtr() : it->second;
    }

private:
    mutable std::mutex mutex;
    std::map<int32_t, AuxPayloadPtr> byPackage;
};

class DecompressEngine : public CodecEngine {
public:
    DecompressEngine(PbgzParameter& para) : CodecEngine(para) {
        readHeadBlockFlag = true;
        registerAuxConsumer(&qualPriorConsumer);
     }

    virtual ~DecompressEngine() { }

    AuxPayloadPtr getQualPrior(int32_t packageIndex) override {
        return qualPriorConsumer.forPackage(packageIndex);
    }

    void setReadHeadFlag(bool flag) {
        readHeadBlockFlag = flag;
    }

protected:
    BlockReader* createBlockReader() override;

    BlockWriter* createBlockWriter() override;

    void releaseBlockReader(BlockReader* &blockReader) override;

    void releaseBlockWriter(BlockWriter* &BlockWriter) override;

    virtual void readBlocks(BlockReader* blockReader);

    virtual void printTailInfo(Timer& costTimer) override {
        PbgzManager::getInstance().printTailInfo(costTimer, false);
        if (parameter.verbose) {
            double secs = costTimer.elapsedSeconds();
            int64_t readLen = PbgzManager::getInstance().getTotalReadLen();
            int64_t writeLen = PbgzManager::getInstance().getTotalWriteLen();
            if (secs > 0) {
                fprintf(stderr, "Decompress speed: %.2f MB/s in (%.2f MB/s out), %.3f s\n",
                        (readLen / (1024.0 * 1024.0)) / secs,
                        (writeLen / (1024.0 * 1024.0)) / secs,
                        secs);
            }
        }
    }

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

private:
    bool initRefGene(PbgzBlockReader* blockReader);

    bool initRefeIndex();

    void readBlockByPostition(BlockReader* blockReader);

    bool unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta);

    bool unpackQualPrior(PbgzBlockReader* blockReader, Json::Value& priorMeta);

    void printFastqFileNotMatchInfo(const Json::Value& metaRefe); 

    virtual Reference* getReference() override { return pRefGene; }

private:
    bool readHeadBlockFlag;
    QualPriorConsumer qualPriorConsumer;
    std::string refPosChrName;
    uint32_t refPosBegin;
    uint32_t refPosEnd;
};
