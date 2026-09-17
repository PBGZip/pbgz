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

    /* See Actuator::ioError(); the aggregation point for out-of-bounds errors
     * of this block's processing */
    const coder_err_sink& ioError() const { return ioErrSink; }

protected:
    /*
     * The compression level from the engine parameters; returns 0 when there is
     * no engine (decompression/test engine), in which case the caller skips the
     * level setting and each encoder keeps its own default level.
     */
    uint8_t engineCompressLevel() const {
        return (pbgzEngine != nullptr) ? pbgzEngine->getParameter().compressLevel : 0;
    }

    /*
     * Create an encoder for a field, with the type taken from the
     * trial-compression result of the preprocessing phase.
     *
     * fallback is the field's fixed encoder (the one its row in
     * kSamFieldCoderConfig registers as the default). All three of the following
     * cases fall back to it, so there is always something to compress with:
     *   1. The engine provides no preprocessing info (both decompression and
     *      test engines return a null pointer)
     *   2. Preprocessing has not finished yet, or it failed
     *   3. The field's sample is too small and was skipped, so there is no
     *      trustworthy trial-compression result
     *
     * In other words, this path can only make the encoder better; a
     * preprocessing problem can never leave nothing to compress with. The level
     * comes from the engine parameters (engineCompressLevel), is converted
     * according to the finally used encoder type, and is written into coder_io
     * for the encoder to read at construction / first encoding.
     */
    std::shared_ptr<coder> makeFieldEncoder(uint32_t fieldIdx, CoderType fallback,
                                            coder_io* io, bool lineMode)
    {
        const uint8_t mode = (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                                    : (uint8_t)PBGZ_MODE_ARCHIVE;
        const uint8_t level = engineCompressLevel();

        /*
         * The -m mode picks the field's default coder from the config table:
         * archive uses the ratio-oriented fallback, fast the speed-oriented one.
         * The caller's fallback only applies to fields the table does not cover.
         *
         * For fast, the fast fallback is coder_rans for every field, which is
         * what the fast path has always used for generic fields - deterministic
         * and quick (per-symbol O(1), no BWT/context model). The block format is
         * self-describing (the magic travels in meta), so the decoding side
         * needs no preset information.
         */
        const CoderType modeFallback = samFieldDefaultCoder(fieldIdx, mode, fallback);

        CoderType picked = modeFallback;
        const PreprocessInfo* preInfo = (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
        if (preInfo != nullptr) {
            picked = preInfo->coderFor(fieldIdx, modeFallback);
        }

        /*
         * Guard against a selection result from another profile - e.g. a stale
         * archive verdict in a fast run. A coder this mode may not select is
         * replaced by the mode's own default. The default itself is always
         * allowed: it is the field's registered coder for this mode.
         */
        if (picked != modeFallback && !CoderFactory::eligibleProfile(picked, mode)) {
            picked = modeFallback;
        }

        CoderFactory::applyLevel(io, picked, level);
        std::shared_ptr<coder> enc = CoderFactory::makeEncoder(picked, io);

        /*
         * Preprocessing trial compression always measures whole-block, so the
         * picked encoder may not be usable line by line. When this field must
         * be fed line by line but the picked encoder does not support it, fall
         * back to the field's fixed encoder, whose line-mode usage is verified.
         * Better a slightly worse compression ratio than data that cannot be
         * decompressed. The level must also be reapplied for the fallback type
         * (the two encoder types accept different level ranges).
         */
        if (lineMode && !enc->supportsLineMode()) {
            CoderFactory::applyLevel(io, modeFallback, level);
            enc = CoderFactory::makeEncoder(modeFallback, io);
        }
        return enc;
    }

    /*
     * The only entry point for constructing coder_io. Views created here carry
     * the aggregation sink, so any out-of-bounds access is guaranteed to be
     * seen by the engine; a direct make_shared<coder_io> does not, and that is
     * reserved for trial compression and unit tests — in a trial, running out of
     * capacity is only evidence that "this encoder is unsuitable", not a
     * failure, and should not fail the whole block.
     */
    std::shared_ptr<coder_io> makeCoderIo(const uint8_t* buff, int32_t len, const char* name) {
        return std::make_shared<coder_io>(buff, len, &ioErrSink, name);
    }

    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;
    Json::Value meta;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
