

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


class PbgzEngine {

public:
    PbgzEngine(const PbgzParameter&  para);

    virtual ~PbgzEngine();

    virtual int32_t init();

    virtual int32_t start();

protected:
    virtual void printHeadInfo() { };

    virtual void printTailInfo(Timer&) { }

    virtual int32_t startEnginePreProc() { return 0; }

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

    virtual int32_t actuatorProc(Actuator*, RoughIOBlock*, RoughIOBlock*) = 0; 
    
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

protected:
    BlockingQueue<RoughIOBlock*> freeInputPool;   // Free queue, file reading tasks get blocks from here to read data
    BlockingQueue<RoughIOBlock*> inputDataPool;   // Data reading tasks write completed data to this queue, compression/decompression tasks get data from this queue
    BlockingQueue<RoughIOBlock*> freeOutputPool;  // Compression/decompression tasks get blocks from this queue to write processed data, output tasks write free blocks to this queue after processing
    BlockingQueue<RoughIOBlock*> outputDataPool;  // After compression/decompression is completed, write to this queue, output tasks get data from this queue
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;

    std::vector<std::thread> codecThreads;
    std::thread writeThread;
    int64_t blockId2Write; 

    mutable std::mutex mutex;
    mutable std::condition_variable conditionVar;
    bool syncFlag;
};