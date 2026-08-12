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

class CodecActuator {

public:
    virtual int32_t decompress() = 0;
    virtual int32_t compress() = 0;

    CodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr): inBlockPtr(inPtr), outBlockPtr(outPtr), pbgzEngine(engine) {};

    virtual ~CodecActuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

    /* 见 Actuator::ioError()，本次块处理的越界错误汇聚点 */
    const coder_err_sink& ioError() const { return ioErrSink; }

#if defined(TEST_MODE) || defined(GTEST_ENABLED)
    // Test-only constructor: creates a dummy engine wrapper
    static PbgzEngine* createTestEngine(const PbgzParameter& para);
#endif

protected:
    /*
     * 引擎参数里的压缩级别；无引擎（解压/测试引擎）时返回 0，调用方据此跳过
     * level 设置，编码器保持各自的默认 level。
     */
    uint8_t engineCompressLevel() const {
        return (pbgzEngine != nullptr) ? pbgzEngine->getParameter().compressLevel : 0;
    }

    /*
     * 为某个字段创建编码器，类型取自预处理阶段的试压结果。
     *
     * fallback 传本字段原先写死的那个编码器。以下三种情况都会原样退回 fallback，
     * 使行为和接通预处理选择之前完全一致：
     *   1. 引擎不提供预处理信息（解压引擎、测试引擎都返回空指针）
     *   2. 预处理还没跑完，或者跑失败了
     *   3. 该字段样本太小被跳过，没有可信的试压结果
     *
     * 也就是说这条路径只会让编码器变得更好，不会因为预处理出问题而压不出东西。
     * level 取自引擎参数（engineCompressLevel），按最终使用的编码器类型换算后写入
     * coder_io，编码器构造/首次编码时读取。
     */
    std::shared_ptr<coder> makeFieldEncoder(uint32_t fieldIdx, CoderType fallback,
                                            coder_io* io, bool lineMode)
    {
        CoderType picked = fallback;
        const PreprocessInfo* preInfo = (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
        if (preInfo != nullptr) {
            picked = preInfo->coderFor(fieldIdx, fallback);
        }

        const uint8_t level = engineCompressLevel();
        CoderFactory::applyLevel(io, picked, level);
        std::shared_ptr<coder> enc = CoderFactory::makeEncoder(picked, io);

        /*
         * 预处理试压一律按整块方式测量，选出来的编码器未必能按逐行方式使用。
         * 本字段要逐行喂数据、而选中的编码器不支持时，退回该字段原本写死的那个，
         * 它的用法是经过验证的。宁可压缩率差一点，也不能压出解不开的数据。
         * 退回的类型也要重设 level（两种编码器能接受的级别范围不同）。
         */
        if (lineMode && !enc->supportsLineMode()) {
            CoderFactory::applyLevel(io, fallback, level);
            enc = CoderFactory::makeEncoder(fallback, io);
        }
        return enc;
    }

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
    Json::Value meta;
    PbgzEngine* pbgzEngine;
    coder_err_sink ioErrSink;
};
