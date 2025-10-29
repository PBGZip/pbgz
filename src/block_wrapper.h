#pragma once

#include <cstring>

#include "io_block.h"
#include "io_wrapper.h"
#include "pbgz_file_wrapper.h"
#include "utils/memory_util.h"

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
    uint8_t cache[BLOCK_SIZE];  // Used to store remaining data from previous read
    uint64_t cacheLen;        // Length of data in cache
};


class PbgzBlockReader : public BlockReader {
public:
    PbgzBlockReader(IOReader* reader): BlockReader(reader) {
        pbgzFileReader = MemoryUtil::safeNewClass<PbgzFileReader>(ioReader);;
    }

    ~PbgzBlockReader() {
       MemoryUtil::safeDeleteClass(pbgzFileReader);
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual int32_t init() override;

    const PbgzFileMeta& getFileMeta() {
        return pbgzFileReader->getFileMeta();
    }

    const PbgzFileHeader& getFileHeader() {
         return pbgzFileReader->getFileHeader();
    }

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
        pbgzFileWriter = MemoryUtil::safeNewClass<PbgzFileWriter>(ioWriter);;
    }
    virtual int32_t writeBlock(RoughIOBlock* blockPtr);

    void setFileMeta(PbgzFileMeta& fileMeta) {
        pbgzFileWriter->setFileMeta(fileMeta);
    }

    virtual int32_t init();

    virtual ~PbgzBlockWriter() {
        MemoryUtil::safeDeleteClass(pbgzFileWriter);
    }
    
private:
    PbgzFileWriter* pbgzFileWriter;
};