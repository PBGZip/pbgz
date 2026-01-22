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

#include "pbgz_engine.h"
#include "actuator.h"
#include <mutex>

class CompressEngine : public PbgzEngine {
public:
    CompressEngine(PbgzParameter& para) : PbgzEngine(para) {
        
    }

    virtual ~CompressEngine();

    virtual int32_t init();

protected:
    virtual BlockReader* createBlockReader() override;
    
    virtual BlockWriter* createBlockWriter() override;

    virtual int32_t engineStartPreProc() override;

    virtual int32_t engineStartAfterProc() override;

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr);

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock* outBlockPtr);

    virtual bool isPrintRatio() {
        return true;
    }

    virtual void setDataBlockPosition(uint32_t blockId);

    void startWriteIndexTask();

    virtual int32_t beforeCoderProc();

private:
    bool initReference();

    int64_t packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true);

    uint32_t calcPackRefeBlockSize();

private:
    std::map<uint32_t, std::vector<int64_t>> blockRefePos;
    std::map<int64_t, uint32_t> blockRefeIndex;
    BlockingQueue<RoughIOBlock*> indexBlockQueue;
    BlockingQueue<RoughIOBlock*> freeIndexBlockQueue;  
};
