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
    }

    int32_t init();

    int32_t start();

    ~PbgzEngine();

private:
    int32_t startReadTask();

    int32_t startCoderTask();

    int32_t startWriteTask();

    bool initRefGeneForDecomress(PbgzBlockReader* blockReader);

    void unpackReference(PbgzBlockReader* blockReader);

    int64_t packReference(std::vector<RoughIOBlock*>& blockVec, int64_t &maxBlockLen, int64_t &totalEncLen);

    bool initReferenceForCompress();

private:
    BlockingQueue<RoughIOBlock*> freeInputPool;   // 空闲的队列，文件读取任务从这里获取block块去读取数据
    BlockingQueue<RoughIOBlock*> inputDataPool;   // 数据读取任务获取完数据写入到此队列，压缩/解压任务从此队列获取数据进行数据
    BlockingQueue<RoughIOBlock*> freeOutputPool;  // 压缩/解压任务从此队列中获取block写入处理后的数据，输出任务处理后空闲块写入此队列
    BlockingQueue<RoughIOBlock*> outputDataPool;  // 压缩/解压完成之后写入此队列，输出任务从此队列中获取数据
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;

    std::vector<std::thread> coderThreads;
    std::thread writeThread;
    Reference* pRefGene;

    PbgzFileMeta fileMeta;
};
