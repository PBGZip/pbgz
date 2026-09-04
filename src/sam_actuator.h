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
#include "coder_bwt_cm.h"

// Forward declaration
class CompressEngine;

/* Hash function for the TLEN mate index keyed by (pos, pnext). */
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
     * Pre-decode the whole block's POS/CIGAR/PNEXT so that the full mate index
     * and reference spans are ready; only then can the main loop reconstruct
     * TLEN via computeTLEN while emitting lines, allowing exceptional values to
     * compress to near zero.
     */
    int32_t preDecodeForTLEN();

    int32_t decompressRegularField(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    /* The OPTION field (all tags from column 12 on) is compressed/decompressed by CRAM-style tag columnization. Currently disabled (OPTION goes through affix); code kept for future use. */
    int32_t compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    int32_t decompressOptionField(uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock,
                                  const Json::Value& fieldMeta);

    int32_t decompressPNextFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    int32_t decompressPosFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    int32_t decompressIdField(uint32_t fieldIdx, Json::Value& fieldMeta, RoughIOBlock* outputBlock);

    int32_t decompressChrName(uint32_t fieldIdx, uint32_t lineNo, RoughIOBlock* outputBlock);

    int32_t decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& nposOffset, uint32_t& totalBaseLen, RoughIOBlock* outputBlock);

    /* Read one record's worth of bytes from the SEQ match stream.
       When matchBlockDecode is set (coder_fc, block-only) the bytes are sliced from
       the pre-decoded matchBlockBuffer; otherwise they are decoded line by line.
       Returns the number of bytes produced, or -1 on failure. */
    int32_t readMatchLine(uint32_t fieldIdx, uint8_t* dst, uint32_t len, uint32_t lineNo);

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
    /* QNAME-specific: cross-line deduplication + position-based modeling (see coder_qname.h). When trialLines > 0, only the first N lines are compressed (used for trial-based selection). */
    int32_t compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                                 uint32_t trialLines = 0);

    // ID field whole compression
    int32_t compressIdFieldInAll(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // Regular field compression
    int32_t compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* The coder type chosen by preprocessing; returns fallback when the engine provides none or the decision is not yet made. */
    CoderType pickedCoderFor(uint32_t fieldIdx, CoderType fallback) const;

    /* A single CIGAR operation (op char + length). */
    struct CigarOp { char op; uint32_t len; };

    uint32_t parseCigar(uint8_t* cigarString, uint32_t cigarLength);

    /* Parses a CIGAR string into an ordered list of (op, len) operations. */
    void parseCigarOps(uint8_t* cigarString, uint32_t cigarLength, std::vector<CigarOp>& ops);

    /* Only counts CIGAR operations that consume reference sequence (M/D/N/=/X); used for TLEN reconstruction. */
    uint32_t parseCigarRefConsumed(uint8_t* cigarString, uint32_t cigarLength);

    /* PNEXT is compressed as (PNEXT - POS) delta text; the deltas of consecutive lines are far smaller than the raw values, so bwt_cm compresses them better. */
    template<typename CoderType>
    int32_t compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* POS is compressed as unsigned varint (LEB128) deltas against the previous line's POS; the baseline resets at each chromosome switch and the decoder reads byte-by-byte. */
    template<typename CoderType>
    int32_t compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* TLEN is not stored verbatim; it is reconstructed from POS/PNEXT/CIGAR on decompression, storing exceptions only for lines that cannot be reconstructed. */
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
             * fieldSrcLen records the raw text size (including the trailing tab,
             * i.e. currTabPos - prevTabPos); once converted to numbers, the stream
             * holds only fixed-width binary, whose total size is recorded in
             * meta["srclen"] (srcLen).
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
         * Numeric fields have two compression forms: fixed-width binary here;
         * when preprocessing selects affix, compressSamByFields switches to the
         * textual form (compressRegularField, mode="string"). mode is written
         * into meta, and the decompression side selects the decode path
         * accordingly.
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
    /* Reference sequence length consumed per record (sum of CIGAR M/D/N/=/X); used for TLEN reconstruction. */
    std::map<uint32_t, uint32_t> cigarReadLen;
    /* TLEN mate index; key is (pos, pnext), value is the line number. */
    std::unordered_map<std::pair<int64_t, int64_t>, uint32_t, PairInt64Hash> tlenMateIndex;
    /* Per-field compressed stream positions recorded by initDecoder, used by preDecodeForTLEN to rebuild decoders. */
    std::map<uint32_t, uint32_t> fieldIoStart;
    std::map<uint32_t, uint32_t> fieldIoDstLen;
    std::map<uint32_t, int32_t> fieldIoLevel;
    /* Field contents pre-decoded by preDecodeForTLEN (POS/CIGAR/PNEXT); the main loop copies them directly. */
    std::map<uint32_t, std::vector<std::string>> tlenPreDecodedFields;
    /* Cached lines whose TLEN reconstruction on the decompression side failed. */
    std::map<uint32_t, int32_t> tlenCache;
    /* PNEXT exception lines (lineNo -> reconstructed PNEXT value), populated by preDecodeForTLEN for the pnext_qname_rebuild mode. */
    std::map<uint32_t, int64_t> pnextCache;
    /* QNAME per data line, populated by preDecodeForTLEN to rebuild PNEXT from mate pairing. */
    std::vector<std::string> decodedQnames;
    /* Previous line's POS for delta decoding (used by the main-loop fallback path; preDecodeForTLEN uses a local variable). */
    int64_t posDeltaPrev = 0;
    std::vector<std::pair<uint32_t, uint32_t>> unmapedReadLength;

    /* Per-line CIGAR operation list (op char, length), parsed once by
       compressCigar / decompressCigar / preDecodeForTLEN and reused by the
       SEQ reference rebuild (M segments from reference, I/S stored separately)
       and by TLEN reconstruction. Indexed by lineNo - headEndLine. */
    std::vector<std::vector<CigarOp>> cigarOpList;

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
    /* Per-sub-stream offset/length of the QNAME compressed stream, recorded by
       initDecoder so preDecodeForTLEN can rebuild independent decoders and
       pre-decode QNAME without consuming idDecoders used by the main loop. */
    std::vector<uint32_t> idStreamOffsets;
    std::vector<uint32_t> idStreamDstLens;
    /* Coder type name per QNAME sub-stream, for preDecodeForTLEN to rebuild decoders. */
    std::vector<std::string> idStreamCoders;
    /* True when QNAME is compressed/decompressed with coder_qname (see decompressIdField). */
    bool idUsesQnameCoder = false;
    /* State for numeric QNAME sub-streams (compressIdFieldSplit "mode":"numeric").
       Such a stream carries no split-symbol terminator, so it cannot be decoded
       one segment per line the way the textual path does; instead the whole varint
       stream is decoded once per block and one value is served per line. Indexes
       run parallel to idDecoders; entries for textual sub-streams stay null/zero. */
    std::vector<uint8_t*> idNumericBufs;
    std::vector<uint32_t> idNumericLens;
    std::vector<uint32_t> idNumericPos;
    std::vector<uint64_t> idNumericAcc;
    /* Releases the buffers above; called before (re)populating them per block. */
    void clearIdNumericState();
    std::map<uint32_t, std::shared_ptr<coder>> fieldDecoders;
    std::shared_ptr<coder_qual> qualCoder;
    /* Decoder for QUAL when compressed with fcv2; kept null when another coder was used. */
    std::shared_ptr<coder_fcv2> qualFcv2Decoder;
    /* Decoder for QUAL when compressed with bwt_cm; kept null when another coder was used. */
    std::shared_ptr<coder_bwt_cm> qualCmDecoder;

    uint8_t* baseSquashBuffer;
    uint8_t* baseDiffSquashBuffer;
    uint8_t* refeStrecchBuffer;
    /* SEQ match stream decoded as a whole block (coder_fc only supports block
       decompression, unlike coder_bwt_cm which is decoded line by line).
       When matchBlockDecode is true, the whole match stream is decoded up front
       into matchBlockBuffer and each record is then served by slicing it with its
       own length (see decompressBase). */
    bool matchBlockDecode = false;
    uint8_t* matchBlockBuffer = nullptr;
    uint32_t matchBlockLength = 0;
    uint32_t matchBlockOffset = 0;
    /* Byte length of the extra "mval" sub-stream that follows the "m" sub-stream when
       RLE is in use; used to advance readOffset past both sub-streams. */
    uint32_t matchExtraOffset = 0;
    /* Reads with missing QUAL (*): expanded into seqLen '*' on compression so that
       each record's length in the quality stream matches the decompression side
       (which fetches each record by SEQ/CIGAR length). */
    std::vector<uint8_t> qualMissingBuf;

    /* OPTION tag-split columnization: decodes the whole block once and serves lines from cache. Currently disabled; code kept for future use. */
    int32_t decodeOptionColumn(const Json::Value& fieldMeta);
    std::vector<std::string> optionRecLines;
    bool optionCacheEmpty = true;

    const uint8_t atcg4[4] = {'A', 'C', 'T', 'G'};

    uint16_t refPosChrIndex;
    uint32_t refPosBegin;
    uint32_t refPosEnd;
};
