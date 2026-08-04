/*
 * index_actuator.h - Index actuator header file
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
 * FITNESS FOR A PARTICULAR PURPOSE OR NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <json/json.h>
#include <memory>
#include <vector>

#include "actuator.h"
#include "coder/coder_io.h"
#include "coder/coder_bwt_cm.h"
#include "pbgz_index.h"

class IndexActuator : public Actuator {
public:
    IndexActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr);

    virtual ~IndexActuator();

    /*
     * 能建索引的块类型，唯一声明处。
     *
     * 索引要的是"每条 read 落在参考基因组的哪个位置"：SAM 的 RNAME/POS 是原文自带的，
     * FASTQ 则要从 compressBaseWithRef 产出的 mpos 子流里取，取法完全不同，所以支持
     * 哪些类型取决于本执行器实现了哪些解码分支。引擎侧只查询，不重复判断。
     */
    static bool supports(BlockType type) {
        return type == SAM;
    }

    virtual int32_t initial() override;

    virtual int32_t process() override;

private:
    int32_t parseHeader(Json::Value& meta);

    int32_t initDecoders(const Json::Value& streams, int32_t& readOffset);

    int32_t skipUnneededFields(const Json::Value& streams, int32_t& readOffset);

    int32_t decodeAndBuildIndex(uint32_t lineNum);

    int64_t headEndLine;

    // Decoder objects
    std::shared_ptr<coder_io> flagIo, chrIo, posIo;
    coder_bwt_cm* flagDecoder;
    coder_bwt_cm* chrDecoder;
    coder_bwt_cm* posDecoder;

    // For sorting validation
    struct SortKey {
        uint16_t chrIndex;
        uint32_t mapPos;

        bool operator<(const SortKey& other) const {
            if (chrIndex != other.chrIndex) {
                return chrIndex < other.chrIndex;
            }
            return mapPos < other.mapPos;
        }
    };

    std::vector<SortKey> sortKeys;
};
