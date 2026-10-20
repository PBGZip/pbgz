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

#include <deque>
#include <map>
#include <vector>
#include <unordered_map>
#include <functional>

#include "codec_actuator.h"
#include "coder.h"
#include "coder_io.h"
#include "reference.h"
#include "coder_bwt_cm.h"
#include "coder_factory.h"
/* The reference-coded SEQ payload and CIGAR walk, shared with the codec pre-selection (see sam_seq_payload.h). */
#include "sam_seq_payload.h"
#include "seq_stream_util.h"   /* SeqRleSplit, used by the SEQ payload helpers below */
/* The QNAME column's layouts and analysis, shared with the codec pre-selection. */
#include "sam_qname_column.h"

/* The QNAME column's pre-decode state (see preDecodeForTLEN). */
struct IdPredecodeState {
    std::vector<std::shared_ptr<coder_io>> streams;
    std::vector<std::shared_ptr<coder>> decoders;
    std::vector<uint8_t*> numericBufs;
    std::vector<uint32_t> numericLens;
    std::vector<uint32_t> numericPos;
    std::vector<uint64_t> numericAcc;
};

/* One TLEN-reconstruction column (POS/CIGAR/PNEXT) and its decoder (see TlenColumnPredecoder). */
struct TlenColumnPredecoder {
    std::shared_ptr<coder_io> stream;
    std::shared_ptr<coder> decoder;
    std::string coderName;
    std::string mode;
    /* Column source length: bounds the widest line (sizes the decode buffer). */
    uint32_t srcLen = 0;
};


#include "coder/id_int_model.h"

// Forward declaration
class CompressEngine;

/* Hash function for the TLEN mate index keyed by (pos, pnext). */
struct PairInt64Hash {
    std::size_t operator()(const std::pair<int64_t, int64_t>& key) const {
        return std::hash<int64_t>()(key.first) ^ (std::hash<int64_t>()(key.second) << 1);
    }
};

/* What decoding one line's fields leaves for the fields after it (see decodeSamFields). */
struct SamLineDecodeState {
    uint8_t* basePtr = nullptr;
    uint32_t actualBaseLen = 0;
};

/* What reconstructing one line's TLEN needs (see gatherTLEN). */
struct TlenInput {
    int64_t  pos = 0;
    int64_t  pnext = 0;
    uint32_t refSpan = 0;
    uint32_t mateRefSpan = 0;
};


/* The OPTION column's parse result (see collectOptionTags). */
struct OptionParseState {
    std::vector<std::pair<std::string, std::string>> tagDict;  /* tag name and type, first-seen order */
    std::map<std::string, int> tagId;                          /* name -> slot in tagDict */
    std::vector<std::vector<uint8_t>> recIds;                  /* one line's tag slots, in order */
    std::vector<std::vector<std::string>> tagVals;             /* one column of values per tag */
    uint32_t srcLen = 0;
};

/* A table indexed by data-line number, the shape of this actuator's per-line tracking tables.
   Backed by a deque plus a presence bitmap, with std::map's interface (operator[], find/end,
   count, ascending range-for, clear) so call sites are unchanged. */
template <typename V>
class LineTable {
public:
    struct Node {
        uint32_t first = 0;
        V        second{};
    };

    class Iterator {
    public:
        Iterator() = default;
        Iterator(LineTable* table, uint32_t at) : table_(table), at_(at) {}

        Node& operator*() const { return table_->nodes_[at_]; }
        Node* operator->() const { return &table_->nodes_[at_]; }
        Iterator& operator++()
        {
            ++at_;
            skipAbsent();
            return *this;
        }
        bool operator!=(const Iterator& other) const { return at_ != other.at_; }
        bool operator==(const Iterator& other) const { return at_ == other.at_; }

        /* Advances past the absent slots; public because begin() uses it. */
        void skipAbsent()
        {
            while (at_ < table_->nodes_.size() && !table_->isPresent(at_)) {
                ++at_;
            }
        }

    private:
        LineTable* table_ = nullptr;
        uint32_t   at_ = 0;
    };

    Iterator begin()
    {
        Iterator it(this, 0);
        it.skipAbsent();
        return it;
    }
    Iterator end() { return Iterator(this, (uint32_t)nodes_.size()); }

    Iterator find(uint32_t lineIdx)
    {
        return isPresent(lineIdx) ? Iterator(this, lineIdx) : end();
    }
    size_t count(uint32_t lineIdx) const { return isPresent(lineIdx) ? 1u : 0u; }
    size_t size() const { return count_; }
    bool empty() const { return count_ == 0; }

    V& operator[](uint32_t lineIdx)
    {
        if (lineIdx >= nodes_.size()) {
            nodes_.resize(lineIdx + 1);
        }
        const size_t word = lineIdx >> 6;
        if (word >= present_.size()) {
            present_.resize(word + 1, 0);
        }
        if (!isPresent(lineIdx)) {
            present_[word] |= (uint64_t)1 << (lineIdx & 63);
            nodes_[lineIdx].first = lineIdx;
            ++count_;
        }
        return nodes_[lineIdx].second;
    }

    void clear()
    {
        nodes_.clear();
        present_.clear();
        count_ = 0;
    }

private:
    bool isPresent(uint32_t lineIdx) const
    {
        const size_t word = lineIdx >> 6;
        return word < present_.size() && ((present_[word] >> (lineIdx & 63)) & 1u) != 0;
    }

    std::deque<Node>      nodes_;
    std::vector<uint64_t> present_;
    size_t                count_ = 0;
};

class SamCodecActuator : public CodecActuator {
public:
    SamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr, Reference* pRefeGene = nullptr);
    virtual ~SamCodecActuator() override;

    /* preAnalysis' two passes over a block: header lines, then record lines walked as raw (pointer, length) views. */
    int32_t parseHeaderLine(const std::string& line);
    int32_t scanDataLine(const uint8_t* line, uint32_t lineLen, uint32_t idx,
                         std::pair<uint8_t, uint32_t>* qualityFrequnce);
    int32_t preAnalysis();

    int32_t compress() override;
    int32_t decompress() override;

    int32_t decompressHeader(RoughIOBlock* outputBlock);

    // Field-by-field decompression
    /* One line's fields, in order, appended to the output block (see decodeSamFields). */
    int32_t decodeSamFields(uint32_t lineNo, uint32_t fieldCount, Json::Value& streams,
                            RoughIOBlock* outputBlock, uint8_t* pBaseOut, uint32_t& totalBaseLen,
                            SamLineDecodeState& st);

    int32_t decompressSamByFields(RoughIOBlock* outputBlock);

    /* initDecoder's per-field setup: ID/QNAME, SEQ, and the lazily-decoded offset-only columns. */
    int32_t initIdFieldDecoders(Json::Value& idMeta);
    int32_t initSeqFieldDecoders(Json::Value& baseMeta, uint32_t idx, RoughIOBlock* outputBlock);
    bool recordDeferredFieldOffset(uint32_t idx, Json::Value& fieldMeta);
    int32_t initDecoder(RoughIOBlock* outputBlock);

    /* Pre-decode the whole block's POS/CIGAR/PNEXT so TLEN can be reconstructed in the main loop. */
    int32_t preDecodeForTLEN();

    int32_t decompressRegularField(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    /* The OPTION field (all tags from column 12 on) is compressed/decompressed by CRAM-style tag columnization. Currently disabled (OPTION goes through affix); code kept for future use. */
    /* The OPTION column's parse into OptionParseState (see collectOptionTags). */
    int32_t collectOptionTags(OptionParseState& st);

    int32_t compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    int32_t decompressOptionField(uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock,
                                  const Json::Value& fieldMeta);

    int32_t decompressPNextFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    int32_t decompressPosFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock);

    /* One split segment of a line's ID, whatever layout wrote it (see reconstructIdSegment). */
    int32_t reconstructIdSegment(uint32_t splitIdx, Json::Value& splitMeta, RoughIOBlock* outputBlock,
                                 uint32_t& idLength);

    int32_t decompressIdField(uint32_t fieldIdx, Json::Value& fieldMeta, RoughIOBlock* outputBlock);

    int32_t decompressChrName(uint32_t fieldIdx, uint32_t lineNo, RoughIOBlock* outputBlock);

    int32_t decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& totalBaseLen, RoughIOBlock* outputBlock);

    /* Read one record's bytes from the SEQ match stream; sliced from matchBlockBuffer when
       matchBlockDecode, else decoded per line. Returns bytes, or -1. */
    int32_t readMatchLine(uint32_t fieldIdx, uint8_t* dst, uint32_t len, uint32_t lineNo);

    int32_t decompressQuality(uint8_t* basePtr, uint32_t actualBaseLen, RoughIOBlock* outputBlock);

    int32_t decompressTLen(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock, const Json::Value& fieldMeta);

    int32_t computeTLEN(uint32_t lineIdx, bool minusOne);

    /* Fills `in` with one line's TLEN inputs; false when the rule makes TLEN 0 (not a proper pair on one reference). */
    bool gatherTLEN(uint32_t lineIdx, TlenInput& in);

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
    /* The block's closing metadata: the SAM meta, the block hash and the meta encoding. */
    int32_t finalizeSamBlockMeta(Json::Value& samMeta, const Json::Value& streamMeta,
                                 uint32_t lineNumber, uint32_t fieldCount, uint32_t totalSrcLen,
                                 uint32_t totalDstLen);

    /* The QNAME column: affix segmentation vs coder_qname, decided by a trial (see encodeQnameColumn). */
    int32_t encodeQnameColumn(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    int32_t compressSamByFields();

    /* The QNAME column's lines as the layout encoders want them (see sam_qname_column.h). */
    QnameColumnInput qnameColumnInput() const;

    /* The two QNAME layout encoders, writing the whole block (see sam_qname_column.h). */
    int32_t compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    /* QNAME-specific: cross-line deduplication + position-based modeling (see coder_qname.h). */
    int32_t compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // ID field whole compression
    int32_t compressIdFieldInAll(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // Regular field compression
    int32_t compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* The OPTION column with repeated tag names/types left out; false when the data does not suit it. */
    bool compressOptionStrip(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                             uint32_t& dstLen);

    /* Reads back what compressOptionStrip wrote: the values and the sequence index per record. */
    int32_t decodeOptionStripColumn(const Json::Value& fieldMeta);

    int32_t compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* The coder type chosen by preprocessing; returns fallback when the engine provides none or the decision is not yet made. */
    CoderType pickedCoderFor(uint32_t fieldIdx, CoderType fallback) const;

    /* The SEQ match stream coder chosen once per file by the codec pre-selection (see seqMatchCoderFor). */
    CoderType seqMatchCoderFor();

    /* The file-level preprocessing result, or nullptr when there is no engine to hold one. */
    PreprocessInfo* preprocessInfoMut() const;


    /* Parses a CIGAR string into an ordered list of (op, len) operations. */

    /* Only counts CIGAR operations that consume reference sequence (M/D/N/=/X); used for TLEN reconstruction. */

    /* PNEXT's exception stream: (contentIdx, delta) pairs of records that cannot be rebuilt from their mate. */
    uint32_t writePnextExceptions(const std::vector<std::pair<uint32_t, int64_t>>& exc,
                                  Json::Value& fieldMeta, Json::Value& metaStreams,
                                  uint32_t& totalDstLen);

    /* PNEXT is compressed as (PNEXT - POS) delta text; the deltas of consecutive lines are far smaller than the raw values, so bwt_cm compresses them better. */
    int32_t compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* POS is compressed as unsigned varint (LEB128) deltas against the previous line's POS; the baseline resets at each chromosome switch and the decoder reads byte-by-byte. */
    int32_t compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* TLEN is not stored verbatim; it is reconstructed from POS/PNEXT/CIGAR on decompression, storing exceptions only for lines that cannot be reconstructed. */
    int32_t compressTLen(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    template<typename T>
    int32_t compressNumber(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
        std::vector<size_t>& npos = inBlockPtr->getNpos();
        uint32_t lineNum = npos.size();
        uint8_t* buffer = inBlockPtr->getBuffer();

        // Create encoder for regular field compression
        std::shared_ptr<coder_io> numberIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        std::shared_ptr<coder> numberCoder =
            CoderFactory::makeProfileEncoder(fieldIdx, compressMode(), numberIo.get(),
                                             CoderType::BWT_CM, true, engineCompressLevel());

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
            /* fieldSrcLen is the raw text size (with trailing tab); meta["srclen"] is the fixed-width binary size. */
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
        /* Numeric fields: fixed-width binary here, or textual (compressRegularField,
           mode="string") when affix is selected; mode in meta picks the decode path. */
        fieldMeta["mode"] = "number";

        LOG_INFO("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldIdx, fieldSrcLen, numberIo->data_len, (double)(numberIo->data_len * 100)/(double)fieldSrcLen);

        return numberIo->data_len;
    }

    int32_t compressChrName(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressBaseWithoutRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* A reference the block builds from its own alignments, for when no usable file reference exists (see buildBlockConsensus). */
    struct BlockConsensus {
        int64_t windowStart = 0;
        int64_t windowBases = 0;
        std::vector<uint8_t> squash;
    };

    bool buildBlockConsensus(uint32_t fieldIdx, BlockConsensus& out);

    /* SEQ coded against a block-built consensus, or the plain path when there is none (see compressBaseSelfRef). */
    int32_t compressBaseSelfRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* compressBaseWithRef's base-length auxiliary stream (see the definition). */
    int32_t writeSeqBaseLengthStream(Json::Value& metaSubs, Json::Value& metaStreams,
                                     uint32_t& totalSrcLen, uint32_t& totalDstLen);
    /* The SEQ match stream: sub-stream "m" and, under RLE, sub-stream "mval". */
    int32_t writeSeqMatchStreams(const SeqRleSplit& rle, const uint8_t* matchBuffer, uint32_t matchLen,
                                 uint32_t srcLen, coder_io* matchIo, Json::Value& metaSubs,
                                 Json::Value& metaStreams, uint32_t& totalDstLen);

    int32_t compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    /* One SEQ exception sub-stream: encode its positions and describe it in `streamMeta`. */
    int32_t writeSeqExceptionStream(SeqExceptionClass& exc, uint8_t ch, const char* sname,
                                    uint32_t srcLen, const uint8_t* src, uint32_t runCount,
                                    Json::Value& streamMeta, uint32_t& totalSrcLen,
                                    uint32_t& totalDstLen);
    /* All SEQ exception streams, one per character that occurred, in ascending order. */
    int32_t writeSeqExceptionStreams(std::vector<SeqExceptionClass>& excClasses,
                                     std::vector<uint8_t>& excPresent, Json::Value& streamMeta,
                                     uint32_t& totalSrcLen, uint32_t& totalDstLen);
    /* The SEQ reference path (sam_seq_payload.h): shared with the codec pre-selection. */

    /* The SEQ column's decoder setup, by layout: whole column as one stream, or the per-base match stream. */
    int32_t initSeqWholeBlockDecoder(const Json::Value& baseMeta, uint32_t idx,
                                     RoughIOBlock* outputBlock);
    int32_t predecodeSeqMatchStream(const Json::Value& streams, uint32_t& streamId,
                                    uint32_t firstDstLength, coder_io* matchIo, uint32_t idx,
                                    uint32_t& extraOffset);
    int32_t predecodeSeqWholeMatchStream(const Json::Value& streams, uint32_t streamId,
                                         coder_io* matchIo, uint32_t idx);
    int32_t predecodeSeqBaseLengthStream(const Json::Value& streams, uint32_t streamId);

    /* The block's own consensus decoded from the block (see predecodeSeqSelfRefStream). */
    int32_t predecodeSeqSelfRefStream(const Json::Value& streams, uint32_t streamId);
    int32_t predecodeSeqExceptionStreams(const Json::Value& baseMeta, const Json::Value& streams,
                                         uint32_t& streamId);
    /* The QUAL column's stream: its frequency table, then the decoder the magic names. */
    int32_t initQualFieldDecoders(const Json::Value& qualMeta);
    /* Any ordinary column: one stream, decoder built from its magic (see CoderFactory::makeFieldDecoder). */
    int32_t initFieldDecoder(uint32_t fieldIdx, const Json::Value& fieldMeta);
    /* The TLEN column: one exception stream, else an ordinary column. */
    int32_t initTlenFieldDecoders(const Json::Value& tlenMeta, uint32_t fieldIdx);
    int32_t decodeTlenExceptionStream(const Json::Value& stream, const Json::Value& tlenMeta);
    /* The QNAME column's own decoders, rebuilt so the main loop's idDecoders are not consumed. */
    int32_t predecodeIdStreams(IdPredecodeState& state);
    /* FLAG (decides whether PNEXT holds a delta) and RNAME (POS-delta chromosome switches). */
    int32_t predecodeFlagColumn();
    int32_t predecodeRnameColumn();
    /* POS(3)/CIGAR(5)/PNEXT(7): build the column's pre-decoder, then rebuild it line by line. */
    bool buildTlenColumnPredecoder(uint32_t fieldIdx, TlenColumnPredecoder& out);
    int32_t rebuildTlenColumn(uint32_t fieldIdx, const TlenColumnPredecoder& pre);
    /* QNAME (field 0) into decodedQnames - the pairing key the PNEXT mode below needs. */
    int32_t rebuildQnameColumn();
    /* PNEXT (7) in pnext_qname_rebuild mode: exception pairs plus QNAME-grouped mates. */
    int32_t rebuildPnextByQname();
    /* The (pos, pnext) -> lineNo index the TLEN walk reads. */
    void buildTlenMateIndex();

    // Helper methods
    void setReference(Reference* ref) { pRefeGene = ref; }

    /* The QUAL column's record loop; returns the column's source length (without the frequency-table stream). */
    /* The quality column's frequency-table auxiliary stream (see writeQualFreqStream). */
    int32_t writeQualFreqStream(Json::Value& subMeta, Json::Value& streamMeta,
                                uint32_t& totalSrcLength, uint32_t& totalDstLength);

    uint32_t encodeQualRecords(qual_record_encoder* encoder, uint32_t lineNum, uint32_t fieldIdx,
                               bool needStrand);
    int32_t compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t buildSamIndex();

    /* Whether the SEQ column may be coded against the loaded reference (see referenceUsableForSeq); cached per instance. */
    bool referenceUsableForSeq();

 private:
    int64_t headEndLine;
    uint32_t samLine;
    /* Every content line's tab positions, indexed by line number minus headEndLine (see LineTabs). */
    std::vector<LineTabs> contentPos;

    Reference* pRefeGene;

    /* referenceUsableForSeq's cached verdict; false until it has been taken. */
    bool refeUsableChecked = false;
    bool refeUsableForSeq = false;

    /* The block consensus the SEQ column is currently coded against, or null (see compressBaseSelfRef). */
    const BlockConsensus* seqSelfRef = nullptr;

    /* The same consensus on the decoding side, read back out of the block (see
       predecodeSeqSelfRefStream) before any record is decoded. */
    std::vector<uint8_t> baseSelfRefSquash;
    int64_t baseSelfRefWindowStart = 0;
    int64_t baseSelfRefWindowBases = 0;
    bool baseSelfRefActive = false;

    /* What the ID analysis found in this block's QNAME column (see sam_qname_column.h). */
    IdSplitAnalysis idAnalysis;
    uint16_t maxFieldSize = 0;
    std::vector<std::pair<int64_t, uint16_t>> lineFiledCount;
    LineTable<int64_t>  mappedPos;
    LineTable<uint16_t> mappedChr;
    LineTable<int64_t>  nextMappedPos;
    LineTable<uint16_t> nextMappedChr;
    LineTable<uint16_t> mappedFlag;
    /* Reference sequence length consumed per record (sum of CIGAR M/D/N/=/X); used for TLEN reconstruction. */
    LineTable<uint32_t> cigarReadLen;
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

    /* Per-line CIGAR op list, parsed once and reused by the SEQ reference rebuild and TLEN. */
    std::vector<std::vector<CigarOp>> cigarOpList;

private:
    uint32_t baseNCount;
    /* Set while compressBaseWithoutRef drives the shared reference-path compressor: no reference is
       consulted, so every record stores its own bases (the "streams" layout still applies). */
    bool seqNoReference = false;

    /* A decoded SEQ exception stream for the current block: its character and its positions. */
    struct SeqExcStream {
        uint8_t byte;
        std::vector<uint32_t> pos;
        uint32_t off = 0;
    };
    std::vector<SeqExcStream> seqExc;

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
    /* Per-sub-stream offset/length of the QNAME stream, for preDecodeForTLEN to rebuild decoders. */
    std::vector<uint32_t> idStreamOffsets;
    std::vector<uint32_t> idStreamDstLens;
    /* Coder type name per QNAME sub-stream, for preDecodeForTLEN to rebuild decoders. */
    std::vector<std::string> idStreamCoders;
    /* True when QNAME is compressed/decompressed with coder_qname (see decompressIdField). */
    bool idUsesQnameCoder = false;
    /* State for numeric QNAME sub-streams ("mode":"numeric"): whole varint stream decoded once, one value per line. */
    std::vector<uint8_t*> idNumericBufs;
    std::vector<uint32_t> idNumericLens;
    std::vector<uint32_t> idNumericPos;
    std::vector<uint64_t> idNumericAcc;
    /* State for value-domain QNAME sub-streams ("mode":"intd"): a range decoder walks one value per line. */
    std::vector<RangeCoder> idIntCoders;
    std::vector<id_int::Model> idIntModes;
    std::vector<uint8_t> idIntActive;
    /* Text of a constant QNAME sub-stream ("mode":"const"), stored once. */
    std::vector<std::string> idConstTexts;
    /* A "dict" sub-stream's texts, indexed by its per-line index. */
    std::vector<std::vector<std::string>> idDictEntries;
    /* A "hexd" sub-stream's alphabet, one code per character. */
    struct IdHexAlphabet {
        std::string chars;
        uint32_t minLen = 0;
        uint32_t maxLen = 0;
    };
    std::vector<IdHexAlphabet> idHexAlphabets;
    std::vector<RangeCoder> idHexCoders;
    /* Releases the buffers above; called before (re)populating them per block. */
    void clearIdNumericState();
    std::map<uint32_t, std::shared_ptr<coder>> fieldDecoders;
    /* Decoder for the QUAL column (see CoderFactory::makeQualDecoder). */
    std::shared_ptr<qual_record_decoder> qualDecoder;

    uint8_t* baseSquashBuffer;
    uint8_t* baseDiffSquashBuffer;
    uint8_t* refeStrecchBuffer;
    /* SEQ match stream decoded as a whole block (coder_fc); records are sliced from matchBlockBuffer. */
    bool matchBlockDecode = false;
    uint8_t* matchBlockBuffer = nullptr;
    uint32_t matchBlockLength = 0;
    uint32_t matchBlockOffset = 0;
    /* Byte length of the "mval" sub-stream under RLE, to advance readOffset past both. */
    uint32_t matchExtraOffset = 0;
    /* Reads with missing QUAL (*), expanded to seqLen '*' on compression. */
    std::vector<uint8_t> qualMissingBuf;

    /* OPTION tag-split columnization: decodes the whole block once and serves lines from cache. */
    int32_t decodeOptionColumn(const Json::Value& fieldMeta);
    /* Recompute a derived tag (MD:Z / NM:i) for the record currently being decoded. */
    bool deriveOptionTagValue(uint32_t tagIdx, std::string& out);
    /* The decoded OPTION text of every record in the block, with optionRecOff giving each record's extent. */
    std::vector<uint8_t> optionRecText;
    std::vector<uint32_t> optionRecOff;
    /* Per record, whether it has an OPTION column at all (see decodeOptionStripColumn). */
    std::vector<uint8_t> optionRecPresent;
    bool optionCacheEmpty = true;

    /* tag_split ("mode":"tag_split") state, decoded once per block and reassembled per line. */
    struct OptionTagCol {
        std::string enc;                    /* "text" / "varint" / "dvarint" / "derived" */
        std::vector<std::string> vals;      /* normal columns, one entry per tag-bearing record */
        std::vector<uint8_t> flags;         /* derived: per value, 0 = derive, 1 = stored */
        std::vector<std::string> excs;      /* derived: stored values, in record order */
        size_t cursor = 0;
        size_t excCursor = 0;
    };
    std::vector<std::vector<uint8_t>> optionRecIds;     /* per line: tag ids present */
    std::vector<std::string> optionTagNames;
    std::vector<std::string> optionTagTypes;
    std::vector<OptionTagCol> optionTagCols;
    /* Current record's alignment, set by the decode loop for MD/NM regeneration. */
    uint16_t optCtxChrId = 0xFFFF;
    uint64_t optCtxPos = 0;
    uint16_t optCtxFlag = 0;
    const uint8_t* optCtxSeq = nullptr;
    uint32_t optCtxSeqLen = 0;
    int32_t optCtxContentIdx = -1;
    std::vector<uint8_t> optMdScratch;

    const uint8_t atcg4[4] = {'A', 'C', 'T', 'G'};

    uint16_t refPosChrIndex;
    uint32_t refPosBegin;
    uint32_t refPosEnd;
};
