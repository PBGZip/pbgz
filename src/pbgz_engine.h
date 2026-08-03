/*
 * pbgz_engine.h - Head file for pbgz project
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

#include <thread>
#include <list>

#include "pbgz_types.h"
#include "blocking_queue.h"
#include "io_block.h"
#include "io_wrapper.h"
#include "block_wrapper.h"
#include "actuator.h"
#include "utils/timer.h"
#include "reference.h"

struct PreprocessInfo;

using BlockingQueueType = BlockingQueue<RoughIOBlock*>;

class PbgzEngine {
public:
    PbgzEngine(const PbgzParameter&  para);

    virtual ~PbgzEngine();

    virtual int32_t init();

    virtual int32_t start();

protected:
    virtual void printHeadInfo() { };

    virtual void printTailInfo(Timer&) { }

    virtual uint32_t getBlockSize();

    virtual int32_t startEnginePreProc() { return 0; }

    /*
     * 准备基础文件元信息，必须在 startWriteTask 之前调用。
     *
     * 为什么要单独开这个钩子：写线程一启动就会执行 createBlockWriter，里面立刻
     * 调用 writeBaseFileMeta 把 baseFileMeta 落盘。如果此时主线程还在往
     * baseFileMeta 里塞内容，两边就在同一份 JSON 上赛跑——谁先跑完决定了写进
     * 文件头的元信息里有没有那部分内容。
     *
     * 这个竞争的实际后果是压缩结果不可复现：同一份输入压两次，文件头元信息可能
     * 相差一整个 JSON 成员，后续所有内容随之整体偏移，两次输出的字节完全不同。
     * 两个结果都能正确解压（解压侧对成员缺失做了兼容），所以问题长期没有暴露，
     * 但它让"同一输入必得同一输出"这条基本性质不成立，也让任何字节级回归对比
     * 都无法进行。
     *
     * 因此约定：凡是要写进基础文件元信息的内容，都必须在本钩子里填好；
     * 需要写线程就绪之后才能做的事（例如往下游写参考块）留在 startWorkPreProc。
     */
    virtual int32_t prepareFileMeta() { return 0; }

    virtual int32_t startWriteTask();

    virtual int32_t startWorkPreProc() { return 0; }

    virtual int32_t startWorkTask();

    virtual int32_t startReadPreProc() { return 0; }

    virtual int32_t startReadTask();

    virtual int32_t startEnginePostProc() { return 0; }

    virtual BlockReader* createBlockReader() = 0;

    virtual BlockWriter* createBlockWriter() = 0;

    virtual void releaseBlockReader(BlockReader* &blockReader) = 0;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) = 0;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) = 0;

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*); 
    
    virtual int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType);

    virtual void readBlocks(BlockReader* blockReader);

    virtual void updateInputStatics(RoughIOBlock*) { }

    virtual void updateOutputStatics(RoughIOBlock*) { }

    virtual void writeBlockPreProc(BlockWriter*, RoughIOBlock*) { }

    virtual void writeBlockPostProc() { }

    virtual void writeFilePreProc() { }

    virtual void writeFilePostProc(BlockWriter*) { }

    virtual void writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr);

private:
    bool isNeedNotify(bool flag);

public:
    std::unique_ptr<BlockingQueueType> freeInputPool;   // Free queue, file reading tasks get blocks from here to read data
    std::unique_ptr<BlockingQueueType> inputDataPool;   // Data reading tasks write completed data to this queue, compression/decompression tasks get data from this queue
    std::unique_ptr<BlockingQueueType> freeOutputPool;  // Compression/decompression tasks get blocks from this queue to write processed data, output tasks write free blocks to this queue after processing
    std::unique_ptr<BlockingQueueType> outputDataPool;  // After compression/decompression is completed, write to this queue, output tasks get data from this queue
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;
    const PbgzParameter& getParameter() const { return parameter; }
    virtual Reference* getReference() { return nullptr; }

    /* File preprocessing result (codec pre-selection). Only the compression
       engine populates this; other engines return nullptr. */
    virtual const PreprocessInfo* getPreprocessInfo() { return nullptr; }

    std::vector<std::thread> workThreads;
    std::thread writeThread;
    int64_t blockId2Write; 

    uint32_t blockCount;

    mutable std::mutex mutex;
    mutable std::condition_variable conditionVar;
    bool syncFlag;
};
