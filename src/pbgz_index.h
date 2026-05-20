/*
 * pbgz_index.h - Header file for pbgz project
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

#include <stdint.h>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <set>


class BlockPosition {
public:
    void setBlockPosition(uint32_t blockId, int64_t position);

    static BlockPosition& getInstance();

    const std::map<uint32_t, int64_t>& getBlockPosition() {
        return blockPositions;
    }

    int64_t getBlockPosition(uint32_t blockId) {
        if (blockPositions.find(blockId) == blockPositions.end()) {
            return -1;
        }
        return blockPositions[blockId];
    }

private:
    std::map<uint32_t, int64_t> blockPositions;
};

struct SamIndexItem{
    int64_t referenceMapPos;
    uint32_t readNumber;
    uint32_t blockId;
    int64_t fileOffset;

    SamIndexItem() {
        referenceMapPos = -1;
        readNumber = 0;
        blockId = 0;
        fileOffset = 0;
    }

    bool operator<(const SamIndexItem& other) const {
        return referenceMapPos < other.referenceMapPos;
    }
};

class SamIndex {
public:
    SamIndex() {
    }

    static SamIndex& getInstance();

    void addSamIndex(uint16_t chrIndex, int64_t refPos, uint32_t readNumber, uint32_t blockId, int64_t offset = 0); 

    int32_t getSamBlockByRef(uint16_t chrIndex, int64_t beginRefPos, int64_t endRefPos, 
        std::set<std::pair<uint32_t, int64_t>>& outBlockList);

    void dumpToFile(std::string& fileName);

    void loadFromFile(std::string& fileName);

    const std::map<uint16_t, std::vector<SamIndexItem>>& getSamIndexList() const {
        return samIndexList;
    }

    void clear() {
        samIndexList.clear();
    }

    void updateFileOffsetsFromBlockPosition();

private:
    std::map<uint16_t, std::vector<SamIndexItem>> samIndexList;
};