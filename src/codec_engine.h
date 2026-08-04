/*
 * codec_engine.h - Header file for pbgz project
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

#include <stdint.h>
#include <list>
#include <vector>
#include <thread>
#include <memory>

#include "io_block.h"
#include "blocking_queue.h"
#include "pbgz_types.h"
#include "io_wrapper.h"
#include "reference.h"
#include "block_wrapper.h"
#include "actuator.h"
#include "pbgz_index.h"
#include "pbgz_engine.h"
#include "utils/path_util.h"
#include "pbgz_manager.h"

class CodecEngine : public PbgzEngine {
public:
    CodecEngine(const PbgzParameter& para) : PbgzEngine(para) {
        pRefGene = nullptr;
        refeOffsetFLag = false;
        blockCount = 0;
    }

    virtual ~CodecEngine();

protected:

    void updateReferenceOffset(int64_t offset);

    void resetReferenceOffset();

    virtual void setDataBlockPosition(uint32_t) {
        return;
    }

    virtual void printHeadInfo() override {
        PbgzManager::getInstance().printHeadInfo(parameter);
    }

    virtual void fileDecisionProc(RoughIOBlock* firstBlock) override {
        PbgzManager::getInstance().printFileType(firstBlock->getBlockType());
    }

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*) {
        return actuator;
    }

    virtual int32_t engineStartAfterProc() {
        if (parameter.isRemoveOriginFile) {
            PathUtil::removeFile(parameter.inputFile);
        }

        return 0;
    }

    virtual void updateInputStatics(RoughIOBlock* inBlockPtr) { 
        PbgzManager::getInstance().updateReadDataLen(inBlockPtr);
    }

    virtual void updateOutputStatics(RoughIOBlock* outBlockPtr) {
        PbgzManager::getInstance().updateWriteDataLen(outBlockPtr);
    }

    virtual void writeFilePostProc(BlockWriter* blockWriter) override;

    virtual void writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr) override;

    virtual void writeBlockPreProc(BlockWriter*, RoughIOBlock* outblockPtr) override;

protected:
    Reference* pRefGene;
    PbgzFileMeta baseFileMeta;
    PbgzFileMeta dynamicFileMeta;
    // 写线程补 refe.offset 与尾部线程补参考统计会并发改同一 JSON，必须共用此锁。
    std::mutex dynamicFileMetaMutex;
    bool refeOffsetFLag;

    virtual Reference* getReference() { return pRefGene; }
};
