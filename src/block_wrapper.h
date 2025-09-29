#pragma once

#include <cstring>

#include "io_block.h"
#include "io_wrapper.h"
#include "pbgz_file_wrapper.h"

namespace BlockUtil {
    bool isFastqBlock(BlockType type);

    std::string getBlockTypeName(BlockType type);
}

class BlockReader {
public:
    BlockReader(IOReader* reader): ioReader(reader) {
        cacheLen = 0;
        memset(cache, 0, sizeof(cache));
        blockId = 0;
    }

    virtual ~BlockReader() {
        ioReader = nullptr;
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW);

    virtual int32_t init() { return 0; }

private:
    BlockType constructBlock(RoughIOBlock* blockPtr);
    
protected:
    int64_t blockId;
    IOReader* ioReader;
    uint8_t cache[BLOCK_SIZE];  // 用于存储上次读取后剩余的数据
    uint64_t cacheLen;        // cache中数据的长度
};


class PbgzBlockReader : public BlockReader {
public:
    PbgzBlockReader(IOReader* reader): BlockReader(reader) {
        pbgzFileReader = nullptr;
    }

    ~PbgzBlockReader() {
        if (pbgzFileReader != nullptr) {
            delete pbgzFileReader;
            pbgzFileReader = nullptr;
        }
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType) override;

    virtual int32_t init() override;

private:
    PbgzFileReader* pbgzFileReader;
};


class BlockWriter {
public:
    BlockWriter(IOWriter* pIoWriter): ioWriter(pIoWriter) { }

    virtual int32_t writeBlock(RoughIOBlock* blockPtr);

    virtual int32_t init() { return 0; } 

    virtual ~BlockWriter() {
        ioWriter = nullptr;
    }

protected:
    IOWriter* ioWriter;
};


class PbgzBlockWriter : public BlockWriter {
public:
    PbgzBlockWriter(IOWriter* pIoWriter) : BlockWriter(pIoWriter) {
        pbgzFileWriter = nullptr;
    }
    virtual int32_t writeBlock(RoughIOBlock* blockPtr);

    virtual int32_t init();

    virtual ~PbgzBlockWriter() {
        if (pbgzFileWriter != nullptr) {
            delete pbgzFileWriter ;
            pbgzFileWriter = nullptr;
        }
    }
    
private:
    PbgzFileWriter* pbgzFileWriter;
};