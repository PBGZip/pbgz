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
     * Block types for which an index can be built; the single place this is declared.
     *
     * An index needs "which position in the reference genome each read lands on": SAM
     * carries RNAME/POS natively, while FASTQ must read them from the mpos substream
     * produced by compressBaseWithRef - the extraction differs completely, so which
     * types are supported depends on which decode branches this actuator implements.
     * The engine side only queries; it does not duplicate the decision.
     */
    static bool supports(BlockType type) {
        /* BAM block content is SAM text, so it can be indexed too */
        return type == SAM || type == SAM_GZIP || type == BAM;
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
