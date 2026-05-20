/*
 * sam_sort.h - Header file for pbgz project
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

#include <vector>

#include "io_wrapper.h"
#include "sam_info.h"


struct SortedSamItem {
    uint32_t fileIndex;
    std::string samLine;

    SortedSamItem() {
        fileIndex = 0;
    }

    int64_t parseData(std::string& readLine) {
        size_t colonPos = readLine.find(':');
        if (colonPos != std::string::npos) {
            // Parse the string after ":"
            samLine = readLine.substr(colonPos + 1);
            // Parse the number before ":"
            return std::stoll(readLine.substr(0, colonPos));
        } 

        // No ":" found, treat entire line as samLine and set referencePos to -1
        samLine = readLine;
        return  -1;
        
    }
};


template <typename Data>
class SamFileSort {
protected:
    struct CompareItem {
        // Reference genome position
        int64_t dataPos;

        // Index of data in dataItem
        uint32_t dataIndex;

        // 
        uint32_t dataExtPos;


        bool operator>(const CompareItem& other) const {
            if (dataPos != other.dataPos) {
                return dataPos > other.dataPos;
            }
            
            return dataExtPos > other.dataExtPos;
        }

        bool operator<(const CompareItem& other) const {
            if (dataPos == other.dataPos) {
                return dataExtPos < other.dataExtPos;
            }
            return dataPos < other.dataPos;
        }
    };

public:
    Data& top() {
        return dataItem[minHeap.top().dataIndex];
    }

    int64_t topSortPos() {
        return minHeap.top().dataPos;
    }

    bool empty() const {
        return minHeap.empty();
    }

    size_t size() const {
        return minHeap.size();
    }

    void push(int64_t compareValue,  Data& data, uint32_t extCompareValue = 0) {
        uint32_t index;
        if (freeIndices.empty()) {
            index = dataItem.size();
            usedFlags.resize(index + 1, true);
            dataItem.push_back(data);
        } else {
            index = freeIndices.back();
            freeIndices.pop_back();
            usedFlags[index] = true;
            dataItem[index] = data;
        }
        minHeap.push({compareValue, index, extCompareValue});
    }

    void pop() {
        if (minHeap.empty()) return;
        uint32_t index = minHeap.top().dataIndex;
        minHeap.pop();
        usedFlags[index] = false;
        freeIndices.push_back(index);
    }

private:
    std::priority_queue<CompareItem, std::vector<CompareItem>, std::greater<CompareItem>> minHeap;
    std::vector<Data> dataItem;
    std::vector<bool> usedFlags;
    std::vector<uint32_t> freeIndices;
};