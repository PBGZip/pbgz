/*
 * block_wrapper.h - Header file for pbgz project
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

#include <cstring>

#include "io_block.h"
#include "io_wrapper.h"
#include "pbgz_file_wrapper.h"
#include "utils/memory_util.h"

namespace BlockUtil {
    bool isFastqBlock(BlockType type);
    bool isSAMBlock(BlockType type);

    /*
     * 辅助块：不携带用户数据，不进"读→并行压缩→按序写"这条流水线，
     * 由文件级偏移按需 seek 加载（参考基因组、其索引、QUAL 先验都属此类）。
     * 流水线入口只问这一个谓词，新增辅助块类型不必再去每个入口补分支。
     */
    bool isAuxiliaryBlock(BlockType type);

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

    PbgzFileMeta& getBaseFileMeta() {
        return pbgzFileReader->getBaseFileMeta();
    }

    PbgzFileMeta& getDynamicFileMeta() {
        return pbgzFileReader->getDynamicFileMeta();
    }

    PbgzFileHeader& getFileHeader() {
         return pbgzFileReader->getFileHeader();
    }

    uint64_t getCurrentFileStart() const {
        return pbgzFileReader->getCurrentFileStart();
    }

    uint64_t getCurrentBlockStart() const {
        return pbgzFileReader->getCurrentBlockStart();
    }

    int32_t getCurrentFileIndex() const {
        return pbgzFileReader->getCurrentFileIndex();
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

    void setBaseFileMeta(PbgzFileMeta& fileMeta) {
        pbgzFileWriter->setBaseFileMeta(fileMeta);
    }

    void setDynamicFileMeta(PbgzFileMeta& fileMeta) {
        pbgzFileWriter->setDynamicFileMeta(fileMeta);
    }

    virtual int32_t init();

    virtual ~PbgzBlockWriter() {
        MemoryUtil::safeDeleteClass(pbgzFileWriter);
    }

    int32_t writeBaseFileMeta();

    int32_t writeDynamicFileMeta();

    void updateHeadExt();
    
private:
    PbgzFileWriter* pbgzFileWriter;
};
