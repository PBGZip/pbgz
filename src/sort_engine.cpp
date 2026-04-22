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

#include <algorithm>
#include <set>

#include "sort_engine.h"
#include "pbgz_manager.h"
#include "sam_info.h"
#include "sam_sort_actuator.h"
#include "utils/path_util.h"
#include "sam_sort.h"

SortEngine::~SortEngine() {

}

uint32_t SortEngine::getBlockSize() {
    return 256 << 20; /// 512M 
}

void SortEngine::printHeadInfo() {
    PbgzManager::getInstance().printHeadInfo(parameter);
}

void SortEngine::printTailInfo(Timer& costTimer) {
    fprintf(stderr, "\nSort finish, cost %um%us.\n", (costTimer.elapsed() / 1000) / 60, (costTimer.elapsed() / 1000) % 60);
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
    if (sortActutor->initial() != 0) {
        LOG_ERROR("SortActutor initial failed.");
        MemoryUtil::safeDeleteClass(sortActutor);
        return nullptr;
    }
    return sortActutor;
}


int32_t SortEngine::startEnginePostProc() {
    /// 将生成的排序好的文件进行归并
    /// 先将头部写入文件
    uint32_t writeBlockId = blockCount;
    std::string headFileName = getSortedHeadFileName();
    RoughIOBlock* outBlock = freeOutputPool->get();
    outBlock->reset();
    outBlock->setBlockId(writeBlockId++);
    if (PathUtil::fileExists(headFileName)) {
        FileReader headReader(headFileName);
        headReader.openIO();
        std::string readLine;
        
        while (size_t readed = headReader.readLine(readLine)) {
            if (readed == 0) {
                break;
            }
            if (outBlock->getDataLen() + readed + 1 > outBlock->getBlockSize()) {
                outputDataPool->push(outBlock);
                outBlock = freeOutputPool->get();
                outBlock->reset();
                outBlock->setBlockId(writeBlockId++);
            }
            memcpy(outBlock->getCurrent(), readLine.c_str(), readed);
            outBlock->setDataLen(outBlock->getDataLen() + readed);
            *outBlock->getCurrent() = '\n';
            outBlock->setDataLen(outBlock->getDataLen() + 1);
        }

        if (outBlock->getDataLen() > 0) {
            outputDataPool->push(outBlock);
            outBlock = freeOutputPool->get();
            outBlock->reset();
            outBlock->setBlockId(writeBlockId++);
        }
        headReader.closeIO();
        PathUtil::removeFile(headFileName);
    }

    LOG_DEBUG("Block count: %d", blockCount);

    std::map<uint32_t, FileReader*> sortedFileReader;
    for (uint32_t fileIdx = 0; fileIdx < blockCount; ++fileIdx) {
        std::string fileName = getSortedSamFileName(fileIdx);
        sortedFileReader[fileIdx] = MemoryUtil::safeNewClass<FileReader>(fileName);
        sortedFileReader[fileIdx]->openIO();
        LOG_DEBUG("Open file : %s", fileName.c_str());
    }

    SamFileSort<SortedSamItem> samSorter;
    std::map<uint32_t, SortedSamItem> fisrtUnMapSamLines;
    std::map<uint32_t, std::pair<int64_t, SortedSamItem>> samLineCache;
    uint32_t index = 0;
    for (uint32_t fileIdx = 0; fileIdx < blockCount; ++fileIdx) {
        SortedSamItem oneItem;
        oneItem.blockId = fileIdx;
        std::string readLine;
        int64_t lastSortPos = -1;
        while (true) {
            size_t readedLine = sortedFileReader[fileIdx]->readLine(readLine);
            if (readedLine > 0) {
                int64_t sortPos = oneItem.parseData(readLine);
                if (sortPos >= 0) {
                    // 第一个读取
                    if (lastSortPos == -1) {
                        samSorter.push(sortPos, oneItem, (fileIdx << 16) + index++);
                        lastSortPos = sortPos;
                    } else {
                        // 位置相同，继续加入排序
                        if (lastSortPos == sortPos) {
                            samSorter.push(sortPos, oneItem, (fileIdx << 16) + index++);
                        } else {
                            // 位置不同，进入缓存，结束这个块的读取
                            samLineCache[fileIdx] = std::make_pair(sortPos, oneItem);
                            break;
                        }
                    }
                } else {
                    fisrtUnMapSamLines[fileIdx] = oneItem;
                    break;
                }
            } else {
                break;
            }
        }
    }

    while(!samSorter.empty()) {
        std::set<uint32_t> fileIdxSet;
        int64_t lastSortPos = -1;
        while (!samSorter.empty()) {
            if (lastSortPos != -1 && lastSortPos != samSorter.topSortPos()) {
                break;
            }
            lastSortPos = samSorter.topSortPos();
            SortedSamItem& fisrtItem = samSorter.top();
            if (outBlock->getDataLen() + fisrtItem.samLine.length() > outBlock->getBlockSize()) {
                outputDataPool->push(outBlock);
                outBlock = freeOutputPool->get();
                outBlock->reset();
                outBlock->setBlockId(writeBlockId++);
            }

            memcpy(outBlock->getCurrent(), fisrtItem.samLine.c_str(), fisrtItem.samLine.length());
            outBlock->setDataLen(outBlock->getDataLen() + fisrtItem.samLine.length());
            *outBlock->getCurrent() = '\n';
            outBlock->setDataLen(outBlock->getDataLen() + 1);
            fileIdxSet.insert(fisrtItem.blockId);
            samSorter.pop();
        }

        index = 0;
        for (auto nextFileIdx : fileIdxSet) {
            int64_t lastReadSortPos = -1;
            while (true) {
                SortedSamItem nextItem;
                int64_t nextRefPos = -1;
                if (samLineCache.find(nextFileIdx) != samLineCache.end()) {
                    nextRefPos = samLineCache[nextFileIdx].first;
                    nextItem = samLineCache[nextFileIdx].second;
                    samLineCache.erase(nextFileIdx);
                } else {
                    if (fisrtUnMapSamLines.find(nextFileIdx) != fisrtUnMapSamLines.end()) {
                        break;
                    }
                    std::string nextLine;
                    size_t readedLine = sortedFileReader[nextFileIdx]->readLine(nextLine);
                    if (readedLine == 0) {
                        break;
                    }
                    nextItem.blockId = nextFileIdx;
                    nextRefPos = nextItem.parseData(nextLine);
                }
                if (nextRefPos >= 0) {
                    if (lastReadSortPos == -1) {
                        samSorter.push(nextRefPos, nextItem, (nextFileIdx << 16) + index++);
                        lastReadSortPos = nextRefPos;
                    } else {
                        if (lastReadSortPos == nextRefPos) {
                            samSorter.push(nextRefPos, nextItem, (nextFileIdx << 16) + index++);
                        } else {
                            samLineCache[nextFileIdx] = std::make_pair(nextRefPos, nextItem);
                            break;
                        }
                    }
                } else {
                    fisrtUnMapSamLines[nextFileIdx] = nextItem;
                    break;
                }
            }
        }
        fileIdxSet.clear();
    }

    if (outBlock->getDataLen() > 0) {
        outputDataPool->push(outBlock);
        outBlock = freeOutputPool->get();
        outBlock->reset();
        outBlock->setBlockId(writeBlockId++);
    }

    /// 按照blockid顺序，将未匹配的数据写入sam
    for (uint32_t fileIdx = 0; fileIdx < blockCount; ++fileIdx) {
        if (fisrtUnMapSamLines.find(fileIdx) != fisrtUnMapSamLines.end()) {
            SortedSamItem& unmapItem = fisrtUnMapSamLines[fileIdx];
            memcpy(outBlock->getCurrent(), unmapItem.samLine.c_str(), unmapItem.samLine.length());
            outBlock->setDataLen(outBlock->getDataLen() + unmapItem.samLine.length());
            *outBlock->getCurrent() = '\n';
            outBlock->setDataLen(outBlock->getDataLen() + 1);
        }

        std::string line;
        while (uint32_t readedLength = sortedFileReader[fileIdx]->readLine(line)) {
            if (readedLength == 0) {
                break;
            }
            SortedSamItem item;
            item.parseData(line);
            if (outBlock->getDataLen() + item.samLine.length() > outBlock->getBlockSize()) {
                outputDataPool->push(outBlock);
                outBlock = freeOutputPool->get();
                outBlock->reset();
                outBlock->setBlockId(writeBlockId++);
            }

            memcpy(outBlock->getCurrent(), item.samLine.c_str() , item.samLine.length());
            outBlock->setDataLen(outBlock->getDataLen() + item.samLine.length());
            *outBlock->getCurrent() = '\n';
            outBlock->setDataLen(outBlock->getDataLen() + 1);
        }
        sortedFileReader[fileIdx]->closeIO();
    }

    if (outBlock->getDataLen() > 0) {
        outputDataPool->push(outBlock);
    }

    for (auto rederItem : sortedFileReader) {
        MemoryUtil::safeDeleteClass(rederItem.second);
    }

    for (uint32_t fileIdx = 0; fileIdx < blockCount; ++fileIdx) {
        PathUtil::removeFile(getSortedSamFileName(fileIdx));
    }

    return 0;
}

