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

#include "io_block.h"
#include "blocking_queue.h"
#include "pbgz_types.h"
#include "io_wrapper.h"
#include "reference.h"
#include "block_wrapper.h"

class PbgzEngine {
public:
    PbgzEngine(const PbgzParameter& para) {
        parameter = para;
        ioReader = nullptr;
        ioWriter = nullptr;
        pRefGene = nullptr;
        refeOffsetFLag = false;
        blockCount = 0;
        refeBeginPos = 0;
        refeEndPos = 0;
    }

    int32_t init();

    int32_t start();

    ~PbgzEngine();

private:
    int32_t startReadTask();

    int32_t startCoderTask();

    int32_t startWriteTask();

    bool initRefGeneForDecomress(PbgzBlockReader* blockReader);

    bool initRefeIndexForDecompress(PbgzBlockReader* blockReader);

    bool unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta);

    int64_t packReference(std::vector<RoughIOBlock*>& blockVec, int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true);

    int32_t packBlockRefePosIdx();

    bool initReferenceForCompress();

    void updateReferenceOffset(int64_t offset);

    void resetReferenceOffset();

    void setDataBlockPosition(uint32_t blockId);

    void readBlockByPostition(BlockReader* blockReader);

    int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType);

private:
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
    std::map<uint32_t, std::vector<int64_t>> blockRefePos;
    std::map<int64_t, uint32_t> blockRefeIndex;
    int64_t refeBeginPos;
    int64_t refeEndPos;
};
