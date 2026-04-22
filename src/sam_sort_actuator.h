/*
 * sam_sort_actuator.h - SAM file sorting actuator header file
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
#include <string>
#include "actuator.h"

class SamSortItem {
public:
    int64_t referencePos;
    uint32_t lineId;

    SamSortItem() {
        referencePos = -1;
        lineId = 0;
    }

    bool operator<(const SamSortItem& item) const {
        return referencePos < item.referencePos;
    }
};

std::string getSortedHeadFileName();

std::string getSortedSamFileName(uint32_t blockId);


class SAMSortActuator : public Actuator {
public:
    SAMSortActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr) : Actuator(inPtr, outPtr) {
        headLineNum = 0;
        notifyFlag = false;
    }

    virtual ~SAMSortActuator();

    virtual int32_t initial() override;
    
    virtual int32_t process() override;

    virtual bool getNotifyFlag() {
        return notifyFlag;
    } 

private:
    uint32_t headLineNum;
    std::vector<SamSortItem> mappedSamItem;
    std::vector<SamSortItem> unmappedSamItem;
    bool notifyFlag;
};