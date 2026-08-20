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
 * Claimant of the QUAL prior. When the reader thread passes a QUAL_PRIOR block, it
 * decompresses and stores it; worker threads use it afterwards.
 *
 * The blob is written by the reader thread and read by worker threads with no other
 * synchronization point between them, so this carries its own lock; reads happen only
 * when decoders are constructed, not on the hot path.
 */
class QualPriorConsumer : public AuxBlockConsumer {
public:
    bool claim(RoughIOBlock* blockPtr, int64_t packageIndex) override;

    /*
     * Return the prior for the given package index, or an empty payload if there is
     * none.
     * A shared_ptr is returned so the consumer can hold the data until it finishes
     * decompressing, unaffected by the claimant's subsequent actions.
     */
    AuxPayloadPtr forPackage(int64_t packageIndex) const {
        std::lock_guard<std::mutex> lock(mutex);
        std::map<int64_t, AuxPayloadPtr>::const_iterator it = byPackage.find(packageIndex);
        return (it == byPackage.end()) ? AuxPayloadPtr() : it->second;
    }

private:
    mutable std::mutex mutex;
    std::map<int64_t, AuxPayloadPtr> byPackage;
};

class DecompressEngine : public CodecEngine {
public:
    DecompressEngine(PbgzParameter& para) : CodecEngine(para) {
        readHeadBlockFlag = true;
        registerAuxConsumer(&qualPriorConsumer);
     }

    virtual ~DecompressEngine() { }

    AuxPayloadPtr getQualPrior(int64_t packageIndex) override {
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

    int64_t readBlocks(BlockReader* blockReader) override;

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

    int64_t readBlockByPostition(BlockReader* blockReader);

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
