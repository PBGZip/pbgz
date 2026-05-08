/*
 * decompress_engine.h - Header file for pbgz project
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

#include "codec_engine.h"
#include "actuator.h"


class DecompressEngine : public CodecEngine {
public:
    DecompressEngine(PbgzParameter& para) : CodecEngine(para) {
        readHeadBlockFlag = true;
     }

    virtual ~DecompressEngine() { }

protected:
    BlockReader* createBlockReader() override;

    BlockWriter* createBlockWriter() override;

    void releaseBlockReader(BlockReader* &blockReader) override;

    void releaseBlockWriter(BlockWriter* &BlockWriter) override;

    virtual void readBlocks(BlockReader* blockReader);

    virtual void printTailInfo(Timer& costTimer) override { 
        PbgzManager::getInstance().printTailInfo(costTimer, false);
    }

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

private:
    bool initRefGene(PbgzBlockReader* blockReader);

    bool initRefeIndex();

    void readBlockByPostition(BlockReader* blockReader);

    bool unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta);

    void printFastqFileNotMatchInfo(const Json::Value& metaRefe); 

    void setReadHeafFlag(bool flag) {
        readHeadBlockFlag = flag;
    }

private:
    bool readHeadBlockFlag;

    virtual Reference* getReference() override { return pRefGene; }
};
