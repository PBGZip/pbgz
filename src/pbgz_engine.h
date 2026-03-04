/*
 * pbgz_engine.h - Header file for pbgz project
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
#include "utils/path_util.h"

class PbgzEngine {
public:
    PbgzEngine(const PbgzParameter& para) {
        parameter = para;
        ioReader = nullptr;
        ioWriter = nullptr;
        pRefGene = nullptr;
        refeOffsetFLag = false;
        blockCount = 0;
    }

    virtual int32_t init();

    virtual int32_t start();

    virtual ~PbgzEngine();

protected:
    int32_t startReadTask();

    int32_t startCoderTask();

    int32_t startWriteTask();

    void updateReferenceOffset(int64_t offset);

    void resetReferenceOffset();

    void setDataBlockPosition(uint32_t) {
        return;
    }

    virtual int32_t beforeCoderProc() {
        return 0;
    }

    virtual BlockReader* createBlockReader() = 0;

    virtual BlockWriter* createBlockWriter() = 0;

    int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType);

    virtual void readBlocks(BlockReader* blockReader);

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*) {
        return actuator;
    }

    virtual int32_t actuatorProc(Actuator*, RoughIOBlock*, RoughIOBlock*) {
        return 0;
    }

    virtual int32_t engineStartPreProc() {
        return 0;
    }

    virtual int32_t engineStartAfterProc() {
        if (parameter.isRemoveOriginFile) {
            PathUtil::removeFile(parameter.inputFile);
        }

        return 0;
    }

    virtual bool isPrintRatio() {
        return false;
    }

protected:
    BlockingQueue<RoughIOBlock*> freeInputPool;   // Free queue, file reading tasks get blocks from here to read data
    BlockingQueue<RoughIOBlock*> inputDataPool;   // Data reading tasks write completed data to this queue, compression/decompression tasks get data from this queue
    BlockingQueue<RoughIOBlock*> freeOutputPool;  // Compression/decompression tasks get blocks from this queue to write processed data, output tasks write free blocks to this queue after processing
    BlockingQueue<RoughIOBlock*> outputDataPool;  // After compression/decompression is completed, write to this queue, output tasks get data from this queue
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;

    std::vector<std::thread> coderThreads;
    std::thread writeThread;
    Reference* pRefGene;

    PbgzFileMeta baseFileMeta;
    PbgzFileMeta dynamicFileMeta;
    bool refeOffsetFLag;
    uint32_t blockCount;

    PbgzIndex pbgzIndex;
};
