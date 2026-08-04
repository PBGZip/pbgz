/*
 * actuator.h - Header file for pbgz project
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

#include <json/json.h>

#include "io_block.h"
#include "pbgz_types.h"
#include "coder/coder_io.h"

#include <memory>

class PbgzEngine;

class Actuator {
public:
    virtual int32_t process() = 0;

    virtual int32_t initial() { return 0; }

    virtual int32_t cleanup() { return 0; }

    Actuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzEngine(engine) {};

    virtual ~Actuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    /*
     * 本次块处理期间所有 coder_io 的越界错误都汇到这里。引擎在 actuatorProc 返回后
     * 问一次即可——那里本来就是所有执行器共用的收口，与 taskFailed、零长块补位这些
     * 既有的失败处置在同一处，不必每种执行器各自再造一遍。
     *
     * 声明为虚函数是因为 SAM/FASTQ/二进制这三种走的是 CodecActuator 那套平行基类，
     * 由适配器桥接过来；真正开 coder_io 的是被适配者，所以适配器要把这一问转下去。
     */
    virtual const coder_err_sink& ioError() const { return ioErrSink; }

protected:
    /*
     * 构造 coder_io 的唯一入口。走这里出来的视图自带汇聚点，越界必然被引擎看到；
     * 直接 make_shared<coder_io> 则不带，那是留给试压和单测的——试压时装不下只是
     * "这个编码器不合适"的选择依据，不是故障，不该让整块失败。
     */
    std::shared_ptr<coder_io> makeCoderIo(const uint8_t* buff, int32_t len, const char* name) {
        return std::make_shared<coder_io>(buff, len, &ioErrSink, name);
    }

    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
