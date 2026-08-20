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

    /*
     * First-block serialization flag: the engine lets the coder thread with id==0
     * finish processing block 0 first (including preAnalysis, which fills
     * process-wide shared state such as SamInfo) before releasing the other threads.
     * Defaults to true, meaning the block can be released as soon as it is done; SAM
     * overrides it to return based on notifyFlag (set halfway through compression, so
     * release is guaranteed to happen after preAnalysis). Release happens only once;
     * see PbgzEngine::firstCoderNotify.
     */
    virtual bool getNotifyFlag() { return true; }

    Actuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzEngine(engine) {};

    virtual ~Actuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    /*
     * All coder_io out-of-range errors during this block's processing funnel here. The
     * engine asks once after actuatorProc returns - that is already the shared choke
     * point for every actuator, sitting alongside existing failure handling like
     * taskFailed and zero-length block padding, so each actuator type does not need to
     * reinvent it.
     *
     * It is virtual because SAM/FASTQ/binary all go through the parallel CodecActuator
     * base classes and are bridged in by adapters; the entity that actually opens
     * coder_io is the adaptee, so the adapter must forward this query down.
     */
    virtual const coder_err_sink& ioError() const { return ioErrSink; }

protected:
    /*
     * The only entry point for constructing a coder_io. Views created through here
     * carry the error sink, so out-of-range access is guaranteed to be seen by the
     * engine; calling make_shared<coder_io> directly does not, and that is reserved
     * for trial compression and unit tests - in trial compression "does not fit" is
     * just a selection signal that this coder is unsuitable, not a failure, and must
     * not abort the whole block.
     */
    std::shared_ptr<coder_io> makeCoderIo(const uint8_t* buff, int32_t len, const char* name) {
        return std::make_shared<coder_io>(buff, len, &ioErrSink, name);
    }

    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
