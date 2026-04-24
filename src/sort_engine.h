/*
 * sort_engin.h - Head file for pbgz project
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

#include "pbgz_engine.h"
#include "utils/timer.h"
#include "pbgz_manager.h"


class SamCombineOutputWriter {
public:
    virtual int32_t writerSam(int64_t samSortPos, std::string& samFileLine) = 0;

    virtual void flush() { return; }

    virtual void close() { return; }
};

class SamCombineOutputFileWriter : public SamCombineOutputWriter {
public:
    SamCombineOutputFileWriter(std::string& outputFile) : fileWriter(outputFile) {
        fileWriter.openIO();
    }

    virtual ~SamCombineOutputFileWriter() {
        fileWriter.closeIO();
    }

    virtual int32_t writerSam(int64_t samSortPos, std::string& samFileLine) override;

    virtual void close() { fileWriter.closeIO(); }

private:
    FileWriter fileWriter;
};

class SamCombineOutputBlockWriter : public SamCombineOutputWriter {
public:
    SamCombineOutputBlockWriter(BlockingQueueType* freeOutputPool, BlockingQueueType* outputDataPool) 
      : freePool(freeOutputPool), outputPool(outputDataPool) {
        writeBlockId = 0;
    }

    int32_t initial(uint32_t beginBlockId);

    virtual ~SamCombineOutputBlockWriter() {
    }

    virtual void flush() override;

    virtual int32_t writerSam(int64_t, std::string& samFileLine) override;

private:
    RoughIOBlock* outBlock;
    uint32_t writeBlockId;
    BlockingQueueType* freePool;
    BlockingQueueType* outputPool;
};

class SortEngine : public PbgzEngine {
public:
   
    SortEngine(const PbgzParameter& para) : PbgzEngine(para) {
    }

    virtual ~SortEngine();

protected:
    virtual void printHeadInfo() override;

    virtual void printTailInfo(Timer& costTimer) override;

    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual int32_t startEnginePostProc() override;

    virtual void updateInputStatics(RoughIOBlock* inBlockPtr) { 
        PbgzManager::getInstance().updateReadDataLen(inBlockPtr);
    }

    virtual void updateOutputStatics(RoughIOBlock* outBlockPtr) {
        PbgzManager::getInstance().updateWriteDataLen(outBlockPtr);
    }

    virtual uint32_t getBlockSize() override;

    int32_t combineAllSamFile(uint16_t inputLevel, uint32_t fileNumber, std::vector<std::string>& outputFiles);

    int32_t combineSamFile(std::vector<std::string> fileList, SamCombineOutputWriter* outputWriter);

};