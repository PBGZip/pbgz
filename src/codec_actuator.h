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
#include "pbgz_engine.h"
#include "coder_factory.h"
#include "field_coder_config.h"

class CodecActuator {

public:
    virtual int32_t decompress() = 0;
    virtual int32_t compress() = 0;

    CodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzEngine(engine) {};

    virtual ~CodecActuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    /* Aggregation point for this block's out-of-bounds errors; see Actuator::ioError(). */
    const coder_err_sink& ioError() const { return ioErrSink; }

protected:
    /* Compression level from engine parameters; 0 when there is no engine, so each
       encoder keeps its own default. */
    uint8_t engineCompressLevel() const {
        return (pbgzEngine != nullptr) ? pbgzEngine->getParameter().compressLevel : 0;
    }

    /* Active -m mode; archive when there is no engine (e.g. unit tests). */
    uint8_t compressMode() const {
        return (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                       : (uint8_t)PBGZ_MODE_ARCHIVE;
    }

    /* Create a SAM field encoder, type taken from preprocessing's trial result and
       falling back to the field's registered default when there is no usable result.
       Level comes from engineCompressLevel. */
    std::shared_ptr<coder> makeFieldEncoder(uint32_t fieldIdx, CoderType fallback,
                                            coder_io* io, bool lineMode)
    {
        const uint8_t mode = (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                                    : (uint8_t)PBGZ_MODE_ARCHIVE;

        /* -m mode picks the field default from the config table; the caller's
           fallback applies only to fields the table does not cover. */
        const CoderType modeFallback = samFieldDefaultCoder(fieldIdx, mode, fallback);
        return buildFieldEncoder(modeFallback, fieldIdx, io, lineMode, mode);
    }

    /* FASTQ counterpart of makeFieldEncoder; uses kFastqFieldCoderConfig and
       FastqField indices (which overlap SAM's, so the tables must stay separate). */
    std::shared_ptr<coder> makeFastqFieldEncoder(uint32_t fieldIdx, CoderType fallback,
                                                 coder_io* io, bool lineMode)
    {
        const uint8_t mode = (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                                    : (uint8_t)PBGZ_MODE_ARCHIVE;
        const CoderType modeFallback = fastqFieldDefaultCoder(fieldIdx, mode, fallback);
        return buildFieldEncoder(modeFallback, fieldIdx, io, lineMode, mode);
    }

    /* FASTQ field's decided coder without building it; callers need it to choose how
       to feed the column (SEQ goes whole-block to coder_fc, line by line otherwise). */
    CoderType pickedFastqCoder(uint32_t fieldIdx, CoderType fallback)
    {
        const uint8_t mode = (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                                    : (uint8_t)PBGZ_MODE_ARCHIVE;
        return pickedCoder(fieldIdx, fastqFieldDefaultCoder(fieldIdx, mode, fallback), mode);
    }

    /* Build an encoder of an already-decided type, applying the level. In lineMode, a
       type that cannot be fed line by line is replaced by lineModeFallback (level
       reapplied, since the two types accept different ranges). */
    std::shared_ptr<coder> makeEncoderOfType(CoderType type, coder_io* io, bool lineMode,
                                             CoderType lineModeFallback)
    {
        const uint8_t level = engineCompressLevel();
        CoderFactory::applyLevel(io, type, level);
        std::shared_ptr<coder> enc = CoderFactory::makeEncoder(type, io);
        if (lineMode && !enc->supportsLineMode()) {
            CoderFactory::applyLevel(io, lineModeFallback, level);
            enc = CoderFactory::makeEncoder(lineModeFallback, io);
        }
        return enc;
    }

    /* Shared body of the field-encoder entries: pick the coder (mode default, or
       preprocessing's verdict when usable) and build it. A verdict from another -m
       mode is ignored (see eligibleProfile); a picked coder that cannot be fed line
       by line is replaced by modeFallback. */
    std::shared_ptr<coder> buildFieldEncoder(CoderType modeFallback, uint32_t fieldIdx,
                                             coder_io* io, bool lineMode, uint8_t mode)
    {
        return makeEncoderOfType(pickedCoder(fieldIdx, modeFallback, mode), io, lineMode, modeFallback);
    }

    /* Coder preprocessing decided for a field, or modeFallback when it decided nothing
       usable (no info, not finished, failed, skipped sample, or wrong -m mode). */
    CoderType pickedCoder(uint32_t fieldIdx, CoderType modeFallback, uint8_t mode)
    {
        CoderType picked = modeFallback;
        const PreprocessInfo* preInfo = (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
        if (preInfo != nullptr) {
            picked = preInfo->coderFor(fieldIdx, modeFallback);
        }
        if (picked != modeFallback && !CoderFactory::eligibleProfile(picked, mode)) {
            picked = modeFallback;
        }
        return picked;
    }

    /* Sole entry point for constructing coder_io; views carry the aggregation sink.
       Direct make_shared<coder_io> is reserved for trial compression and unit tests. */
    std::shared_ptr<coder_io> makeCoderIo(const uint8_t* buff, int32_t len, const char* name) {
        return std::make_shared<coder_io>(buff, len, &ioErrSink, name);
    }

    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;
    Json::Value meta;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
