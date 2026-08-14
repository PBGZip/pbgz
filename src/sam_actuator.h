/*
 * sam_actuator.h - Header file for sam_actuator
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

#include <map>
#include <vector>
#include <unordered_map>
#include <functional>

#include "codec_actuator.h"
#include "coder.h"
#include "coder_io.h"
#include "coder_qual.h"
#include "coder_fcv2.h"
#include "reference.h"
#include "sam_info.h"
#include "coder_bwt_cm.h"
#include "coder_qname.h"
#include "pbgz_stat.h"

// Forward declaration
class CompressEngine;

/* TLEN 伙伴索引用 (pos, pnext) 作键的哈希函数。 */
struct PairInt64Hash {
    std::size_t operator()(const std::pair<int64_t, int64_t>& key) const {
        return std::hash<int64_t>()(key.first) ^ (std::hash<int64_t>()(key.second) << 1);
    }
};

class SamCodecActuator : public CodecActuator {
public:
    SamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr, Reference* pRefeGene = nullptr);
    virtual ~SamCodecActuator() override;

    int32_t preAnalysis();

    int32_t compress() override;
    int32_t decompress() override;

    int32_t decompressHeader(RoughIOBlock* outputBlock);

    // Field-by-field decompression
    int32_t decompressSamByFields(RoughIOBlock* outputBlock);

    int32_t initDecoder(RoughIOBlock* outputBlock);

    /*
     * 先解码全块的 POS/CIGAR/PNEXT，把完整伙伴索引与参考跨度建好，主循环再逐行
     * 输出时 TLEN 才能用 computeTLEN 还原，异常值因此可以压到接近零。
     */
    int32_t preDecodeForTLEN();

    int32_t decompressRegularField(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    /* OPTION 字段（第 12 列起的所有 tag）按 CRAM 式 tag 列化压缩/解压。当前未启用（OPTION 走 affix），保留代码。 */
    int32_t compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    int32_t decompressOptionField(uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock,
                                  const Json::Value& fieldMeta);

    int32_t decompressPNextFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    int32_t decompressPosFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    int32_t decompressIdField(uint32_t fieldIdx, Json::Value& fieldMeta, RoughIOBlock* outputBlock);

    int32_t decompressChrName(uint32_t fieldIdx, uint32_t lineNo, RoughIOBlock* outputBlock);

    int32_t decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& nposOffset, uint32_t& totalBaseLen, RoughIOBlock* outputBlock);

    int32_t decompressQuality(uint8_t* basePtr, uint32_t actualBaseLen, RoughIOBlock* outputBlock);

    int32_t decompressTLen(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock, const Json::Value& fieldMeta);

    int32_t computeTLEN(uint32_t lineIdx, bool minusOne);

    template<typename T>
    int32_t decompressNumber(uint32_t fieldIdx, uint32_t lineNo, RoughIOBlock* outputBlock) {
        uint32_t outLen = sizeof(T);
        int32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), outLen, UINT8_MAX, false);
        if (fieldLen < 0 || (uint32_t)fieldLen != outLen) {
            LOG_ERROR("Decode failed, filed = %u, lineNo = %u", fieldIdx, lineNo);
            return -1;
        }

        T val = *(T*)outputBlock->getCurrent();
        std::string strVal = std::to_string(val);
        if (fieldIdx == 1) {
            mappedFlag[lineNo] = val;
        } else if (fieldIdx == 3) {
            mappedPos[lineNo] = val;
        } else if (fieldIdx == 7) {
            nextMappedPos[lineNo] = val;
        }
        memcpy(outputBlock->getCurrent(), strVal.c_str(), strVal.length());
        outputBlock->setDataLen(outputBlock->getDataLen() + strVal.length());
        *outputBlock->getCurrent() = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return strVal.length() + 1;
    }

    int32_t decompressCigar(uint32_t fieldIdx, uint8_t splitFlag, uint32_t lineIdx, RoughIOBlock* outputBlock);

    int64_t getHeadLineNumber() {
        return headEndLine;
    }

    int64_t getSamLineNumber() {
        return samLine;
    }

    void initMetaInfo();

private:
    int32_t preAnalysisIdLine(uint8_t* buffer, uint32_t length);

    int32_t preAnalysisIdFirstLine(uint8_t* buffer, uint32_t length);

    // compress SAM Head
    int32_t compressSamHeader();

    // Field-by-field compression
    int32_t compressSamByFields();

    // ID field split compression
    int32_t compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                                 uint32_t trialLines = 0);
    /* QNAME 专用：跨行去重 + 按位置建模（见 coder_qname.h）。trialLines>0 时只压前 N 行（选型试压用）。 */
    int32_t compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                                 uint32_t trialLines = 0);

    // ID field whole compression
    int32_t compressIdFieldInAll(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // Regular field compression
    int32_t compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* 预处理选出的编码器类型；引擎不提供或尚未决出时返回 fallback。 */
    CoderType pickedCoderFor(uint32_t fieldIdx, CoderType fallback) const;

    uint32_t parseCigar(uint8_t* cigarString, uint32_t cigarLength);

    /* 只统计消耗参考序列的 CIGAR 操作（M/D/N/=/X），用于 TLEN 推算。 */
    uint32_t parseCigarRefConsumed(uint8_t* cigarString, uint32_t cigarLength);

    /* PNEXT 按 (PNEXT - POS) 差值文本压缩；连续行的差值远小于原始值，bwt_cm 压得更小。 */
    template<typename CoderType>
    int32_t compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* POS 按与上一行 POS 的差值文本压缩；实测优于定宽二进制和文本 affix。 */
    template<typename CoderType>
    int32_t compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* TLEN 不存原值，解压时按 POS/PNEXT/CIGAR 推算，只对推算不上的行存异常。 */
    template<typename CoderType>
    int32_t compressTLen(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    template<typename T>
    int32_t compressNumber(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
        std::vector<size_t>& npos = inBlockPtr->getNpos();
        uint32_t lineNum = npos.size();
        uint8_t* buffer = inBlockPtr->getBuffer();

        // Create encoder for regular field compression
        std::shared_ptr<coder_io> numberIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        std::shared_ptr<coder_bwt_cm> numberCoder = std::make_shared<coder_bwt_cm>(numberIo.get());
        CoderFactory::applyLevel(numberIo.get(), CoderType::BWT_CM, engineCompressLevel());

        fieldSrcLen = 0;
        uint32_t srcLen = 0;
        // Process each line and extract the current field
        for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
            uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
            uint32_t lineEnd = npos[lineIdx] - lineStart;

            uint8_t* line = buffer + lineStart;
            // Skip header lines (starting with @)
            if (*line == '@') {
                continue;
            }

            // Middle fields: between tabs
            uint32_t contentIdx = lineIdx - headEndLine;
            uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
            uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
            uint8_t* fieldStart = line + prevTabPos + 1;
            uint32_t fieldLength = currTabPos - prevTabPos - 1;

            std::string str = std::string((char*)fieldStart, fieldLength);
            T value = (T)std::stoll(str);
            /*
             * fieldSrcLen 记录的是原始文本占用（含行尾 tab，即 currTabPos - prevTabPos）；
             * 转成数字后压缩流里只有定宽二进制，其总空间记在 meta["srclen"]（srcLen）。
             */
            fieldSrcLen += fieldLength + 1;
            // Encode the field data
            numberCoder->encode_line(reinterpret_cast<const uint8_t*>(&value), sizeof(T));
            srcLen += sizeof(T);

            if (fieldIdx == 3) {
                mappedPos[lineIdx] = value;
            } else if (fieldIdx == 1) {
                mappedFlag[lineIdx] = value;
            } else if (fieldIdx == 7) {
                nextMappedPos[lineIdx] = value;
            }
        }

        // Flush the encoder for this field
        numberCoder->encode_flush();

        // Update output block data length
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + numberIo->data_len);

        // Set field metadata
        fieldMeta["srclen"] = srcLen;
        fieldMeta["dstlen"] = numberIo->data_len;
        fieldMeta["coder"] = numberIo->meta;
        fieldMeta["field"] = fieldIdx;
        /*
         * 数值字段有两种压缩形态：这里按定宽二进制编码；当预处理选出 affix 时，
         * compressSamByFields 会改走文本形态（compressRegularField，mode="string"）。
         * mode 写进 meta，解压侧据此选择解码路径。
         */
        fieldMeta["mode"] = "number";

        LOG_INFO("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldIdx, fieldSrcLen, numberIo->data_len, (double)(numberIo->data_len * 100)/(double)fieldSrcLen);

        return numberIo->data_len;
    }

    int32_t compressChrName(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressBaseWithoutRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // Helper methods
    void setReference(Reference* ref) { pRefeGene = ref; }

    int32_t compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t buildSamIndex();

 private:
    int64_t headEndLine;
    uint32_t samLine;
    std::vector<std::vector<int64_t>> contentPos;

    Reference* pRefeGene;

    // ID analysis related members (similar to FastqActuator)
    uint32_t idPosLength;
    std::vector<uint8_t> idSplitSymbols;
    std::vector<std::vector<int32_t>> idSplitPos; // Position of each separator in ID field for each line
    std::vector<uint32_t> idSplitMinLen; // Minimum length of each separator
    std::vector<uint32_t> idSplitMaxLen; // Maximum length of each separator
    const std::string idSplitDefault = "/:= _.,-#\r\t\n";
    uint16_t maxFieldSize = 0;
    std::vector<std::pair<int64_t, uint16_t>> lineFiledCount;
    std::map<uint32_t, int64_t> mappedPos;
    std::map<uint32_t, uint16_t> mappedChr;
    std::map<uint32_t, int64_t> nextMappedPos;
    std::map<uint32_t, uint16_t> nextMappedChr;
    std::map<uint32_t, uint16_t> mappedFlag;
    /* 每条记录消耗的参考序列长度（CIGAR 的 M/D/N/=/X 之和），TLEN 推算用。 */
    std::map<uint32_t, uint32_t> cigarReadLen;
    /* TLEN 伙伴索引，键为 (pos, pnext)，值为行号。 */
    std::unordered_map<std::pair<int64_t, int64_t>, uint32_t, PairInt64Hash> tlenMateIndex;
    /* initDecoder 记下的各字段压缩流位置，供 preDecodeForTLEN 重建解码器。 */
    std::map<uint32_t, uint32_t> fieldIoStart;
    std::map<uint32_t, uint32_t> fieldIoDstLen;
    std::map<uint32_t, int32_t> fieldIoLevel;
    /* preDecodeForTLEN 预解码的字段内容（POS/CIGAR/PNEXT），主循环直接拷贝。 */
    std::map<uint32_t, std::vector<std::string>> tlenPreDecodedFields;
    /* 解压侧 TLEN 推算异常的缓存行。 */
    std::map<uint32_t, int32_t> tlenCache;
    /* POS 差值解码时的上一行 POS（主循环兜底路径用，preDecodeForTLEN 用局部变量）。 */
    int64_t posDeltaPrev = 0;
    std::vector<std::pair<uint32_t, uint32_t>> unmapedReadLength;

    uint32_t baseNCount;
    uint32_t* baseNPosBuffer;
    uint32_t* baseLengthBuffer;
    uint32_t minBaseLength = UINT32_MAX;
    uint32_t maxBaseLength = 0;

    // SAM header compression related members
    uint32_t headerSrcLen; // Original length of SAM file header
    uint32_t headerDstLen; // Compressed length of SAM file header

    // Quality compression related members
    std::vector<std::pair<uint16_t, uint16_t>> qualFreqTable;

    uint32_t readOffset;

    std::vector<std::shared_ptr<coder_io>> ioVector;
    std::vector<std::shared_ptr<coder>> idDecoders;
    /* QNAME 用 coder_qname 压缩/解压时为 true（见 decompressIdField）。 */
    bool idUsesQnameCoder = false;
    std::map<uint32_t, std::shared_ptr<coder>> fieldDecoders;
    std::shared_ptr<coder_qual> qualCoder;
    /* QUAL 用 fcv2 压缩时的解码器；用其他编码器压缩时保持为空。 */
    std::shared_ptr<coder_fcv2> qualFcv2Decoder;
    /* QUAL 用 bwt_cm 压缩时的解码器；用其他编码器压缩时保持为空。 */
    std::shared_ptr<coder_bwt_cm> qualCmDecoder;

    uint8_t* baseSquashBuffer;
    uint8_t* baseDiffSquashBuffer;
    uint8_t* refeStrecchBuffer;

    /* OPTION tag-split 列化：整块解一次，逐行从缓存取。当前未启用，保留代码。 */
    int32_t decodeOptionColumn(const Json::Value& fieldMeta);
    std::vector<std::string> optionRecLines;
    bool optionCacheEmpty = true;

    const uint8_t atcg4[4] = {'A', 'C', 'T', 'G'};

    uint16_t refPosChrIndex;
    uint32_t refPosBegin;
    uint32_t refPosEnd;
};
