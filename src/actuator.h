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
#include "coder/coder_io.h"

#include <memory>

class PbgzEngine;

class Actuator {
public:
    virtual int32_t process() = 0;

    virtual int32_t initial() { return 0; }

    virtual int32_t cleanup() { return 0; }

    /*
     * First-block serialization flag: coder thread id==0 finishes block 0 (including
     * preAnalysis) before other threads are released. Defaults to true; SAM overrides
     * it via notifyFlag. See PbgzEngine::firstCoderNotify.
     */
    virtual bool getNotifyFlag() { return true; }

    Actuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzEngine(engine) {};

    virtual ~Actuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    /*
     * Sink for all coder_io out-of-range errors during this block's processing. The
     * engine queries it once after actuatorProc returns. Virtual so adapters forward
     * the query to the adaptee that actually opens coder_io.
     */
    virtual const coder_err_sink& ioError() const { return ioErrSink; }

protected:
    /* Sole entry point for constructing a coder_io; views carry the error sink.
       Direct make_shared<coder_io> is reserved for trial compression and unit tests. */
    std::shared_ptr<coder_io> makeCoderIo(const uint8_t* buff, int32_t len, const char* name) {
        return std::make_shared<coder_io>(buff, len, &ioErrSink, name);
    }

    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
