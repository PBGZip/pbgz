/*
 * sam_sort_actuator.cpp - SAM file sorting actuator implementation file
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

#include <vector>
#include <cstdio>

#include "sam_sort_actuator.h"
#include "sam_info.h"
#include "io_wrapper.h"

SAMSortActuator::~SAMSortActuator() {
    
}

int32_t SAMSortActuator::initial() {
    if (inBlockPtr == nullptr) {
        return -1;
    }

    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    for (uint32_t idx = 0; idx < lineNum; ++idx) {
        uint32_t begin = (idx == 0 ? 0 : npos[idx - 1]  + 1);
        uint32_t end = (idx < npos.size() ? npos[idx] : strlen((char*)buffer + begin));
        if (begin >= end) {
            return -1;
        }
        std::string line((char*)buffer + begin, end - begin);
        if (buffer[begin] == '@') {
            if (line.length() < 3) {
                return -1;
            }
            headLineNum++;

            if (line.substr(0, 3) == "@SQ") {
                if (SamUtil::parseChromosomeInfo(line) != 0) {
                    LOG_ERROR("Parse chromosome info failed");
                    return -1;
                }
            }
        } else {
            uint32_t lastTabPosInlie = 0;
            uint32_t tabCount = 0;
            uint16_t flag = 0xFFFF;
            
            int64_t chrStart = -1;
            for (uint32_t i = 0; i < line.length(); ++i) {
                if (line.at(i) == '\t' || line.at(i) == '\n') {
                    if (tabCount == 1) {
                        flag = std::stoul(line.substr(lastTabPosInlie + 1, i - lastTabPosInlie - 1));
                    } else if (tabCount == 2) {   /// Chromonome Name 
                        std::string chrName = line.substr(lastTabPosInlie + 1, i - lastTabPosInlie - 1);
                        if (chrName != "*") {
                            uint16_t chrIndex = SamInfo::getInstance().getChrNameIndex(chrName);
                            if (chrIndex != 0xFFFF) {
                                chrStart =  SamInfo::getInstance().getPositionByIndex(chrIndex);
                            }
                        }
                    } else if (tabCount == 3) {  /// 
                        SamSortItem samItem;
                        if ((flag & 0x04) == 0 && chrStart != -1) {
                            int64_t mapPos = std::stoull(line.substr(lastTabPosInlie + 1, i - lastTabPosInlie - 1));
                            samItem.referencePos = ((chrStart + mapPos) << 16 ) + (flag & ~0x800);
                            samItem.lineId = idx;
                            mappedSamItem.push_back(samItem);
                        } else {
                            samItem.lineId = idx;
                            unmappedSamItem.push_back(samItem);
                        }
                    }
                    tabCount++;
                    lastTabPosInlie = i;
                } 
            }
        }

        if (inBlockPtr->getNpos().size() > (size_t)headLineNum) {
           notifyFlag = true;
        }
    }

    return 0;
}

int32_t SAMSortActuator::process() {
    std::stable_sort(mappedSamItem.begin(), mappedSamItem.end());

    std::vector<uint32_t>& npos =  inBlockPtr->getNpos();
    // 先将头部单独写入文件
    if (headLineNum > 0) {
        std::string headName = getSortedHeadFileName();
        if (inBlockPtr->getBlockId() == 0) {
            std::remove(headName.c_str());
        }
        FileWriter headFileWrite(headName);
        headFileWrite.openIO();
        headFileWrite.seekToEnd();
        for (uint32_t lineId = 0; lineId < headLineNum; ++lineId) {
            uint32_t begin = (lineId == 0) ? 0 : npos[lineId - 1] + 1;
            uint32_t lineLength = npos[lineId] - begin + 1;
            headFileWrite.writeIO(inBlockPtr->getBuffer() + begin,  lineLength);
        }
        headFileWrite.closeIO();
    }

    /// 已经排序的写入临时文件
    std::string sortSamName = getSortedSamFileName(inBlockPtr->getBlockId());
    std::remove(sortSamName.c_str());
    FileWriter fileWrite(sortSamName);
    fileWrite.openIO();
    fileWrite.seekToEnd();
    for (auto item = mappedSamItem.begin(); item < mappedSamItem.end(); ++item) {
        uint32_t lineId = item->lineId;
        uint32_t begin = (lineId == 0) ? 0 : npos[lineId - 1] + 1;
        uint32_t lineLength = npos[lineId] - begin + 1;
        std::string writeBuff = std::to_string(item->referencePos) + ":";
        fileWrite.writeIO(writeBuff.c_str(), writeBuff.length());
        fileWrite.writeIO(inBlockPtr->getBuffer() + begin,  lineLength);
    }

    /// 未匹配的部分写入文件
    for (auto item = unmappedSamItem.begin(); item < unmappedSamItem.end(); ++item) {
        uint32_t lineId = item->lineId;
        uint32_t begin = (lineId == 0) ? 0 : npos[lineId - 1] + 1;
        uint32_t lineLength = npos[lineId] - begin + 1;
        std::string writeBuff = "-1:";
        fileWrite.writeIO(writeBuff.c_str(), writeBuff.length());
        fileWrite.writeIO(inBlockPtr->getBuffer() + begin,  lineLength);
    }

    fileWrite.closeIO();
    return 0;
}

std::string getSortedHeadFileName() {
    return "sorted_head.sam";
}

std::string getSortedSamFileName(uint32_t blockId) {
    std::string fileName = "sorted_sam_";
    fileName += std::to_string(blockId) + ".sam";
    return fileName;
}