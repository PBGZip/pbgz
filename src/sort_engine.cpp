/*
 * sort_engine.cpp - Source file for pbgz project
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

#include "sort_engine.h"
#include "pbgz_manager.h"
#include "sam_sort_actuator.h"

SortEngine::~SortEngine() {

}

void SortEngine::printHeadInfo() {
    PbgzManager::getInstance().printHeadInfo(parameter);
}

BlockReader* SortEngine::createBlockReader() {
    BlockReader* blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return nullptr;
    }

    if (blockReader->init() != 0) {
        LOG_ERROR("BlockReader init failed.");
        MemoryUtil::safeDeleteClass(blockReader);
        return nullptr;
    }

    return blockReader;
}

BlockWriter* SortEngine::createBlockWriter() {
    BlockWriter* blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("Create block writer failed.");
        return nullptr;
    }

    if (blockWriter->init() != 0) {
        LOG_ERROR("BlockWriter init failed.");
        MemoryUtil::safeDeleteClass(blockWriter);
        return nullptr;
    }
    return blockWriter;
}

void SortEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void SortEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    MemoryUtil::safeDeleteClass(blockWriter);
}

Actuator* SortEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    SAMSortActuator* sortActutor = MemoryUtil::safeNewClass<SAMSortActuator>(inBlockPtr, outBlockPtr);
    return sortActutor;
}

