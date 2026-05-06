/*
 * codec_actuator.h - Head file for pbgz project
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
#include <json/json.h>

#include "io_block.h"
#include "pbgz_types.h"

class CodecActuator {

public:
    virtual int32_t decompress() = 0;
    virtual int32_t compress() = 0;

    CodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzParameter& para): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzPara(para) {};
    
    virtual ~CodecActuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    virtual bool getNotifyFlag() {
        return true;
    } 

protected: 
    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;  
    Json::Value meta;
    PbgzParameter& pbgzPara;
};


