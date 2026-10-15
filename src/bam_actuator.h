/*
 * bam_actuator.h - column-wise BAM compression without SAM text.
 *
 * BamCodecActuator consumes the BamColumns produced by BamBlockReader (see bam_columns.h), encoding each column with the coder preprocessing selected for the matching SAM field.
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

#include "bam_columns.h"
#include "codec_actuator.h"
#include "reference.h"

class BamCodecActuator : public CodecActuator {
public:
    BamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr,
                     Reference* pRef = nullptr);
    virtual ~BamCodecActuator() override {}

    int32_t compress() override;
    int32_t decompress() override;

private:
    /* Encode one serialised column with the coder picked for fieldIdx and append its meta to streamMeta. */
    int32_t encodeColumn(uint32_t fieldIdx, CoderType fallback, const uint8_t* data, size_t size,
                         Json::Value& streamMeta, const char* name, const char* mode);

    /* QNAME must be fed line by line, or the prefix/suffix-matching coder degrades to a context model. */
    int32_t encodeQnameColumn(const uint8_t* qnameStream, size_t qnameStreamSize,
                              uint32_t nRecords, Json::Value& streamMeta);

    /* Whole-column coding with a plain bwt_cm block, for columns whose preselected coder has a fragile whole-block decoder. */
    int32_t encodeRawBwtColumn(uint32_t fieldIdx, const uint8_t* data, size_t size,
                               Json::Value& streamMeta, const char* name);

    /* LEB128 helpers for the numeric columns. */
    static void putUvarint(std::vector<uint8_t>& out, uint64_t v);
    static void putSvarint(std::vector<uint8_t>& out, int64_t v);

    /* Columns that need per-record handling or a transform. */
    int32_t encodePosColumn(Json::Value& streamMeta);
    int32_t encodePnextColumn(Json::Value& streamMeta);
    int32_t encodeTlenColumn(Json::Value& streamMeta);
    int32_t encodeCigarColumn(Json::Value& streamMeta);
    int32_t encodeQualColumn(Json::Value& streamMeta);

    /* SEQ against a reference: per-base diff stream, or plain text SEQ when no reference is usable. */
    int32_t encodeSeqColumnRef(Json::Value& streamMeta);

    uint32_t pickedFor(uint32_t fieldIdx, CoderType fallback) const;

    /* Whether the loaded FASTA is laid out the way the @SQ list addresses it; cached per actuator. */
    bool seqReferenceUsable();

    BamColumns* cols;
    Reference* pRefeGene;

    /* seqReferenceUsable's cached verdict; false until it has been taken. */
    bool refeUsableChecked = false;
    bool refeUsableForSeq = false;
};
