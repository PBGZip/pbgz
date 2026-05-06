/*
 * index_engine.cpp - Source file for pbgz project
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

#include <cstdio>
#include <tuple>

#include "index_engine.h"
#include "index_actuator.h"
#include "pbgz_index.h"
#include "utils/memory_util.h"
#include "log/logger.h"

IndexEngine::~IndexEngine() {
}

BlockReader* IndexEngine::createBlockReader() {
    PbgzBlockReader* pbgzReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
    if (pbgzReader == nullptr) {
        LOG_ERROR("Failed to create block reader.");
        return nullptr;
    }
    if (0 != pbgzReader->init()) {
        LOG_ERROR("Block reader init fail.");
        MemoryUtil::safeDeleteClass(pbgzReader);
        return nullptr;
    }
    return pbgzReader;
}

BlockWriter* IndexEngine::createBlockWriter() {
    BlockWriter* blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("Failed to create block writer.");
        return nullptr;
    }
    if (0 != blockWriter->init()) {
        LOG_ERROR("Block writer init fail.");
        MemoryUtil::safeDeleteClass(blockWriter);
        return nullptr;
    }
    return blockWriter;
}

void IndexEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void IndexEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    MemoryUtil::safeDeleteClass(blockWriter);
}

Actuator* IndexEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    IndexActuator* indexActuator = MemoryUtil::safeNewClass<IndexActuator>(inBlockPtr, outBlockPtr, parameter);
    if (indexActuator == nullptr) {
        LOG_ERROR("Failed to create IndexActuator.");
        return nullptr;
    }
    if (indexActuator->initial() != 0) {
        LOG_ERROR("IndexActuator initial failed.");
        MemoryUtil::safeDeleteClass(indexActuator);
        return nullptr;
    }
    return indexActuator;
}


int64_t IndexEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    RoughIOBlock* blockPtr = freeInputPool->get();
    if (blockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }
    blockPtr->reset();

    int64_t ret = blockReader->readBlock(blockPtr, fileType);
    if (ret <= 0) {
        freeInputPool->push(blockPtr);
        return ret;
    }

    if (blockPtr->getBlockType() == REFERENCE) {
        freeInputPool->push(blockPtr);
        return -2;
    }
    
    if (blockPtr->getBlockType() != SAM) {
        fprintf(stderr, "The index command is only valid for SAM  pbgz files.");
        freeInputPool->push(blockPtr);
        return -1;
    }

    if (blockPtr->getBlockId() == 0 && blockPtr->getTotalDataLen() > 0){ 
        fileType = blockPtr->getBlockType();
    }

    updateInputStatics(blockPtr);
    inputDataPool->push(blockPtr);
    if (ret > 0) {
        blockCount++;
    }

    return ret;
}

int32_t IndexEngine::startEnginePostProc() {
    const auto& samIndexList = SamIndex::getInstance().getSamIndexList();

    // Collect all SamIndexItems and sort by blockId
    std::vector<std::tuple<uint16_t, SamIndexItem>> sortedByBlockId;
    for (const auto& chrIndexPair : samIndexList) {
        uint16_t chrIndex = chrIndexPair.first;
        const auto& items = chrIndexPair.second;
        for (const auto& item : items) {
            sortedByBlockId.emplace_back(chrIndex, item);
        }
    }

    std::sort(sortedByBlockId.begin(), sortedByBlockId.end(),
        [](const auto& a, const auto& b) {
            const auto& itemA = std::get<1>(a);
            const auto& itemB = std::get<1>(b);
            uint16_t chrA = std::get<0>(a);
            uint16_t chrB = std::get<0>(b);

            if (itemA.blockId != itemB.blockId) {
                return itemA.blockId < itemB.blockId;
            }
            if (chrA != chrB) {
                return chrA < chrB;
            }
            return itemA.referenceMapPos < itemB.referenceMapPos;
        });

    // Get output block and format output
    RoughIOBlock* outputPtr = freeOutputPool->get();
    if (outputPtr == nullptr) {
        LOG_ERROR("Failed to get free output block for index output");
        return -1;
    }
    outputPtr->reset();
    outputPtr->setBlockId(blockCount);

    uint8_t* current = outputPtr->getCurrent();
    uint32_t written = 0;
    uint32_t remain = outputPtr->getRemain();

    for (const auto& entry : sortedByBlockId) {
        uint16_t chrIndex = std::get<0>(entry);
        const SamIndexItem& item = std::get<1>(entry);

        int len = snprintf((char*)current, remain, "%hu\t%ld\t%u\t%u\n",
                          chrIndex, item.referenceMapPos, item.readNumber,
                          item.blockId);

        if (len <= 0 || len >= (int)remain) {
            LOG_ERROR("Index output buffer overflow");
            freeOutputPool->push(outputPtr);
            return -1;
        }

        current += len;
        written += len;
        remain -= len;
    }

    outputPtr->setDataLen(written);
    outputDataPool->push(outputPtr);

    SamIndex::getInstance().clear();

    return 0;
}
