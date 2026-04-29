/*
 * sam_index.cpp - Implementation of sam_index
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
 * FITNESS FOR A PARTICULAR PURPOSE OR NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#if defined(__has_feature)
  #if __has_feature(address_sanitizer)
    #include <sanitizer/lsan_interface.h>
  #endif
#endif

#include "sam_index.h"
#include "io_wrapper.h"


SamIndex& SamIndex::getInstance() {
    static SamIndex instance;
    return instance;
}

void SamIndex::addSamIndex(uint16_t chrIndex, int64_t refPos, uint32_t readNumber, uint32_t blockId) {
    SamIndexItem item;
    item.referenceMapPos = refPos;
    item.readNumber = readNumber;
    item.blockId = blockId;
    samIndexList[chrIndex].push_back(item);
}

int32_t SamIndex::getSamBlockByRef(uint16_t chrIndex, int64_t beginRefPos, int64_t endRefPos,
                                   std::vector<uint32_t>& outBlockList) {
    outBlockList.clear();
    if (samIndexList.find(chrIndex) == samIndexList.end()) {
        return 0;
    }

    auto& items = samIndexList[chrIndex];
    SamIndexItem lowerBound, upperBound;
    lowerBound.referenceMapPos = beginRefPos;
    upperBound.referenceMapPos = endRefPos;

    auto itLow = std::lower_bound(items.begin(), items.end(), lowerBound);
    auto itHigh = std::upper_bound(items.begin(), items.end(), upperBound);

    for (auto it = itLow; it != itHigh; ++it) {
        outBlockList.push_back(it->blockId);
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
            int len = snprintf(buffer, sizeof(buffer), "%hu\t%ld\t%u\t%u\n",
                              chrIndex, samItem.referenceMapPos, samItem.readNumber,
                              samItem.blockId);
            fileWriter.writeIO(buffer, len);
        }
    }
    fileWriter.closeIO();
}

void SamIndex::loadFromFile(std::string& fileName) {
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

        if (tokens.size() == 4) {
            uint16_t chrIndex = static_cast<uint16_t>(std::stoul(tokens[0]));
            int64_t referenceMapPos = std::stoll(tokens[1]);
            uint32_t readNumber = static_cast<uint32_t>(std::stoul(tokens[2]));
            uint32_t blockId = static_cast<uint32_t>(std::stoul(tokens[3]));
            addSamIndex(chrIndex, referenceMapPos, readNumber, blockId);
        }
    }
    fileReader.closeIO();
}
