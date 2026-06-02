/*
 * pbgz_index.cpp - Source file for pbgz project
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


#include "pbgz_index.h"
#include "io_wrapper.h"
#include "log/logger.h"


void BlockPosition::setBlockPosition(uint32_t blockId, int64_t position) {
    blockPositions[blockId] = position;
}

BlockPosition& BlockPosition::getInstance() {
    static BlockPosition instance;
    return instance;
}


SamIndex& SamIndex::getInstance() {
    static SamIndex instance;
    return instance;
}

void SamIndex::addSamIndex(uint16_t chrIndex, int64_t refPos, uint32_t readNumber, uint32_t blockId, int64_t offset) {
    SamIndexItem item;
    item.referenceMapPos = refPos;
    item.readNumber = readNumber;
    item.blockId = blockId;
    item.fileOffset = offset;
    samIndexList[chrIndex].push_back(item);
}

int32_t SamIndex::getSamBlockByRef(uint16_t chrIndex, int64_t beginRefPos, int64_t endRefPos,
                                   std::set<std::pair<uint32_t, int64_t>>& outBlockList) {
    outBlockList.clear();
    if (samIndexList.find(chrIndex) == samIndexList.end()) {
        return 0;
    }

    auto& items = samIndexList[chrIndex];
    if (beginRefPos == 0 && endRefPos == 0) {
        for (auto it = items.begin(); it != items.end(); ++it) {
            outBlockList.insert(std::make_pair(it->blockId, it->fileOffset));
        }
    } else {
        SamIndexItem lowerBound, upperBound;
        lowerBound.referenceMapPos = beginRefPos;
        upperBound.referenceMapPos = endRefPos;

        auto itLow = std::lower_bound(items.begin(), items.end(), lowerBound);
        auto itHigh = std::upper_bound(items.begin(), items.end(), upperBound);

        if (itLow != items.begin() || itLow->referenceMapPos != beginRefPos) {
            itLow = std::prev(itLow);
        }
        for (auto it = itLow; it != itHigh; ++it) {
            LOG_DEBUG("BlockId = %d, offset = %d", it->blockId, it->fileOffset);
            outBlockList.insert(std::make_pair(it->blockId, it->fileOffset));
        }
    }

    return outBlockList.size();
}

void SamIndex::dumpToFile(std::string& fileName) {
    FileWriter fileWriter(fileName);
    fileWriter.openIO();

    char buffer[256];
    for (const auto& samIndex : samIndexList) {
        uint16_t chrIndex = samIndex.first;
        for (const auto& samItem : samIndex.second) {
            int len = snprintf(buffer, sizeof(buffer), "%hu\t%ld\t%u\t%u\t%ld\n",
                              chrIndex, samItem.referenceMapPos, samItem.readNumber,
                              samItem.blockId, samItem.fileOffset);
            fileWriter.writeIO(buffer, len);
        }
    }
    fileWriter.closeIO();
}

void SamIndex::loadFromFile(std::string& fileName) {
    LOG_DEBUG("Load index file : %s", fileName.c_str());
    FileReader fileReader(fileName);
    fileReader.openIO();

    samIndexList.clear();

    while(!fileReader.isEOF()) {
        std::string line;
        uint32_t readLen = fileReader.readLine(line);
        if (readLen == 0 || line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }

        if (tokens.size() == 5) {
            uint16_t chrIndex = static_cast<uint16_t>(std::stoul(tokens[0]));
            int64_t referenceMapPos = std::stoll(tokens[1]);
            uint32_t readNumber = static_cast<uint32_t>(std::stoul(tokens[2]));
            uint32_t blockId = static_cast<uint32_t>(std::stoul(tokens[3]));
            uint32_t fileOffset = static_cast<uint32_t>(std::stoul(tokens[4]));
            addSamIndex(chrIndex, referenceMapPos, readNumber, blockId, fileOffset);
        }
    }
    fileReader.closeIO();
}

void SamIndex::updateFileOffsetsFromBlockPosition() {
    BlockPosition& blockPosition = BlockPosition::getInstance();
    for (auto& samIndexPair : samIndexList) {
        for (auto& samIndexItem : samIndexPair.second) {
            int64_t position = blockPosition.getBlockPosition(samIndexItem.blockId);
            if (position != -1) {
                samIndexItem.fileOffset = position;
            }
        }
    }
    return;
}

