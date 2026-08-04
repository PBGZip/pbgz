/*
 * index_engine.h - Header file for pbgz project
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

#include "pbgz_types.h"
#include "pbgz_engine.h"

class IndexEngine : public PbgzEngine {
public:
    IndexEngine(const PbgzParameter& para) : PbgzEngine(para) {
    }

    virtual ~IndexEngine();

protected:
    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual void readBlockPreProc(BlockReader* blockReader) override;

    virtual BlockIntake intakeBlock(BlockReader* blockReader, RoughIOBlock* blockPtr) override;

    virtual int32_t startEnginePostProc() override;

private:
    std::map<uint32_t, int64_t> blockFileOffset;

    /* 当前这一块在源文件里的起始位置，只能在读之前取得；只由读线程访问。 */
    int64_t pendingBlockOffset = 0;
};
