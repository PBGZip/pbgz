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
#include "coder_qcm.h"
#include "coder_factory.h"
/* The reference-coded SEQ payload, and the CIGAR operation it walks: shared with the codec
   pre-selection, which trials that payload on a sample of the first block (see
   sam_seq_payload.h for why it is its own module). */
#include "sam_seq_payload.h"
/* The QNAME column's layouts and the analysis they need: shared with the codec pre-selection,
   which runs both encoders on a sample to decide which one the file's names suit. */
#include "sam_qname_column.h"

/*
 * The QNAME column's pre-decode state (see preDecodeForTLEN): the sub-stream views it reads
 * from - kept alive here, since the decoders point into them - the decoders themselves, one
 * per sub-stream and null where the layout carries no coder, and the numeric sub-streams the
 * rebuild loop serves one value per line from. Defined in sam_actuator.cpp, like the rest of
 * the SEQ/QNAME column's machinery.
 */
struct IdPredecodeState {
    std::vector<std::shared_ptr<coder_io>> streams;
    std::vector<std::shared_ptr<coder>> decoders;
    std::vector<uint8_t*> numericBufs;
    std::vector<uint32_t> numericLens;
    std::vector<uint32_t> numericPos;
    std::vector<uint64_t> numericAcc;
};

/*
 * One of the columns TLEN reconstruction needs (POS, CIGAR, PNEXT) and the decoder built for
 * it: the stream view is kept alive here, since the decoder reads through a raw pointer into
 * it, and the magic travels along because a coder's buffer trait is looked up by it (see
 * CoderFactory::decoderHoldsCallerBuffer).
 */
struct TlenColumnPredecoder {
    std::shared_ptr<coder_io> stream;
    std::shared_ptr<coder> decoder;
    std::string coderName;
    std::string mode;
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

/* What decoding one line's fields leaves behind for the fields after it (see decodeSamFields):
   where this line's SEQ landed and how long it is, which is what QUAL is decoded against. */
struct SamLineDecodeState {
    uint8_t* basePtr = nullptr;
    uint32_t actualBaseLen = 0;
};


/* The OPTION column's parse result: the tag dictionary, what each line carried, and the byte
   count of the OPTION text (see collectOptionTags). */
struct OptionParseState {
    std::vector<std::pair<std::string, std::string>> tagDict;  /* tag name and type, first-seen order */
    std::map<std::string, int> tagId;                          /* name -> slot in tagDict */
    std::vector<std::vector<uint8_t>> recIds;                  /* one line's tag slots, in order */
    std::vector<std::vector<std::string>> tagVals;             /* one column of values per tag */
    uint32_t srcLen = 0;
};

class SamCodecActuator : public CodecActuator {
public:
    SamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr, Reference* pRefeGene = nullptr);
    virtual ~SamCodecActuator() override;

    /* preAnalysis' two passes over a block: the SAM header lines (via the caller's loop, which
       hands them over one by one) and the record lines. */
    int32_t parseHeaderLine(const std::string& line);
    int32_t scanDataLine(const std::string& line, uint32_t idx,
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

    /* initDecoder's per-field steps: the ID/QNAME column (one decoder per split sub-stream,
       plus the payload state each segment layout reads back), the SEQ column (with or without
       a reference), and the fields whose columns are decoded lazily and keep only an offset. */
    int32_t initIdFieldDecoders(Json::Value& idMeta);
    int32_t initSeqFieldDecoders(Json::Value& baseMeta, uint32_t idx, RoughIOBlock* outputBlock);
    bool recordDeferredFieldOffset(uint32_t idx, Json::Value& fieldMeta);
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
    /* The OPTION column's parse: the tag dictionary and what every line carried, into
       OptionParseState (see collectOptionTags). */
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
    /* The block's closing metadata: the SAM meta, the block hash and the meta encoding. */
    int32_t finalizeSamBlockMeta(Json::Value& samMeta, const Json::Value& streamMeta,
                                 uint32_t lineNumber, uint32_t fieldCount, uint32_t totalSrcLen,
                                 uint32_t totalDstLen);

    /* The QNAME column: affix segmentation against coder_qname, decided by a trial and reused
       across blocks (see encodeQnameColumn). */
    int32_t encodeQnameColumn(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    int32_t compressSamByFields();

    /* The QNAME column's lines as the layout encoders want them (see sam_qname_column.h): the
       block being compressed and what its parse recorded. */
    QnameColumnInput qnameColumnInput() const;

    /*
     * The two QNAME layouts, encoding the whole block (see sam_qname_column.h). The module's
     * encoders also take a line bound, which is the trial's business and not this one's: the
     * pre-selection passes it to measure the layout on a sample, while a block being written is
     * always written in full. Which layout runs here is the file-level verdict's
     * (PreprocessInfo::qnameUseAffix, see encodeQnameColumn).
     */
    int32_t compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    /* QNAME-specific: cross-line deduplication + position-based modeling (see coder_qname.h). */
    int32_t compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // ID field whole compression
    int32_t compressIdFieldInAll(uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    // Regular field compression
    int32_t compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    /* The coder type chosen by preprocessing; returns fallback when the engine provides none or the decision is not yet made. */
    CoderType pickedCoderFor(uint32_t fieldIdx, CoderType fallback) const;

    /*
     * The coder for the SEQ match stream, as decided once per file by the codec pre-selection
     * (CodecSelector::selectSeqReferenceCoder): BWT_CM or FC, whichever coded the block's real
     * match payload smaller. A verdict pinned in PreprocessInfo is honored as-is - that is how
     * the tests exercise each candidate - and without one the pre-trial default, coder_bwt_cm,
     * is returned. See the definition.
     */
    CoderType seqMatchCoderFor();

    /* The file-level preprocessing result, or nullptr when there is no engine to hold one. */
    PreprocessInfo* preprocessInfoMut() const;

    uint32_t parseCigar(uint8_t* cigarString, uint32_t cigarLength);

    /* Parses a CIGAR string into an ordered list of (op, len) operations. */
    void parseCigarOps(uint8_t* cigarString, uint32_t cigarLength, std::vector<CigarOp>& ops);

    /* Only counts CIGAR operations that consume reference sequence (M/D/N/=/X); used for TLEN reconstruction. */
    uint32_t parseCigarRefConsumed(uint8_t* cigarString, uint32_t cigarLength);

    /* PNEXT's exception stream: the (contentIdx, delta) pairs of the records that cannot be
       rebuilt from their mate. Returns the stream's encoded size and records in fieldMeta how
       the stream was written (see writePnextExceptions). */
    template<typename CoderType>
    uint32_t writePnextExceptions(const std::vector<std::pair<uint32_t, int64_t>>& exc,
                                  Json::Value& fieldMeta, Json::Value& metaStreams,
                                  uint32_t& totalDstLen);

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

    /* compressBaseWithRef's base-length auxiliary stream, for the layout where base lengths vary
       (see the definition for what it holds and why it is shaped that way). */
    int32_t writeSeqBaseLengthStream(Json::Value& metaSubs, Json::Value& metaStreams,
                                     uint32_t& totalSrcLen, uint32_t& totalDstLen);
    /* The SEQ match stream: sub-stream "m" and, under RLE, sub-stream "mval". */
    int32_t writeSeqMatchStreams(const SeqRleSplit& rle, const uint8_t* matchBuffer, uint32_t matchLen,
                                 uint32_t srcLen, coder_io* matchIo, Json::Value& metaSubs,
                                 Json::Value& metaStreams, uint32_t& totalDstLen);

    int32_t compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);
    /* One SEQ exception sub-stream (see SeqExceptionClass): encode its positions, describe
       the stream in `streamMeta`, and add its sizes to the block's totals. */
    int32_t writeSeqExceptionStream(SeqExceptionClass& exc, uint8_t ch, const char* sname,
                                    uint32_t srcLen, const uint8_t* src, uint32_t runCount,
                                    Json::Value& streamMeta, uint32_t& totalSrcLen,
                                    uint32_t& totalDstLen);
    /* All of them, one per character that occurred, in ascending character order - including
       the per-file choice of form for the 'N' class (see PreprocessInfo::nposForm). */
    int32_t writeSeqExceptionStreams(std::vector<SeqExceptionClass>& excClasses,
                                     std::vector<uint8_t>& excPresent, Json::Value& streamMeta,
                                     uint32_t& totalSrcLen, uint32_t& totalDstLen);
    /* The SEQ column's reference path - whether a record can be coded against the reference, and
       the CIGAR walk that produces its per-base 2-bit payload - is in sam_seq_payload.h: the codec
       pre-selection builds the same payload to trial the match coder on it. */

    /*
     * The SEQ column's decoder setup, split by the two layouts it can be written in (see
     * compressBaseWithRef): the whole column as one stream, or the match stream that carries
     * the per-base 2-bit codes. Both are steps of initDecoder's SEQ branch.
     */
    int32_t initSeqWholeBlockDecoder(const Json::Value& baseMeta, uint32_t idx,
                                     RoughIOBlock* outputBlock);
    int32_t predecodeSeqMatchStream(const Json::Value& streams, uint32_t& streamId,
                                    uint32_t firstDstLength, coder_io* matchIo, uint32_t idx,
                                    uint32_t& extraOffset);
    int32_t predecodeSeqWholeMatchStream(const Json::Value& streams, uint32_t streamId,
                                         coder_io* matchIo, uint32_t idx);
    int32_t predecodeSeqBaseLengthStream(const Json::Value& streams, uint32_t streamId);
    int32_t predecodeSeqExceptionStreams(const Json::Value& baseMeta, const Json::Value& streams,
                                         uint32_t& streamId);
    /* The QUAL column's stream: its frequency table, then the decoder the magic names. */
    int32_t initQualFieldDecoders(const Json::Value& qualMeta);
    /* Any ordinary column: one stream, its decoder built by the magic it carries (see
       CoderFactory::makeFieldDecoder), and its position recorded for preDecodeForTLEN. */
    int32_t initFieldDecoder(uint32_t fieldIdx, const Json::Value& fieldMeta);
    /* The TLEN column: its reconstruction layout carries one exception stream, everything else
       is an ordinary column. */
    int32_t initTlenFieldDecoders(const Json::Value& tlenMeta, uint32_t fieldIdx);
    int32_t decodeTlenExceptionStream(const Json::Value& stream, const Json::Value& tlenMeta);
    /* The QNAME column's own decoders, rebuilt from the recorded sub-stream offsets so the
       main loop's idDecoders are not consumed (see IdPredecodeState). */
    int32_t predecodeIdStreams(IdPredecodeState& state);
    /* The two columns preDecodeForTLEN rebuilds before the POS chain and the TLEN walk: FLAG
       (which decides whether PNEXT holds a delta) and RNAME (which the POS delta chain uses to
       spot chromosome switches). */
    int32_t predecodeFlagColumn();
    int32_t predecodeRnameColumn();
    /* POS(3) / CIGAR(5) / PNEXT(7): build the column's pre-decoder, then rebuild the column
       from it line by line (see TlenColumnPredecoder). */
    bool buildTlenColumnPredecoder(uint32_t fieldIdx, TlenColumnPredecoder& out);
    int32_t rebuildTlenColumn(uint32_t fieldIdx, const TlenColumnPredecoder& pre);
    /* QNAME (field 0) into decodedQnames - the pairing key the PNEXT mode below needs. */
    int32_t rebuildQnameColumn();
    /* PNEXT (field 7) in pnext_qname_rebuild mode: the exception pairs plus the records that
       belong to a QNAME group of two mapped mates. */
    int32_t rebuildPnextByQname();
    /* The (pos, pnext) -> lineNo index the TLEN walk reads. */
    void buildTlenMateIndex();

    // Helper methods
    void setReference(Reference* ref) { pRefeGene = ref; }

    /* The QUAL column's record loop: what the trial picked is fed record by record, so this step
       knows nothing about which coder won (see qual_record_encoder). Returns the column's source
       length (the quality text itself, without the frequency-table auxiliary stream). */
    /* The quality column's frequency-table auxiliary stream (see writeQualFreqStream). */
    int32_t writeQualFreqStream(Json::Value& subMeta, Json::Value& streamMeta,
                                uint32_t& totalSrcLength, uint32_t& totalDstLength);

    uint32_t encodeQualRecords(qual_record_encoder* encoder, uint32_t lineNum, uint32_t fieldIdx,
                               bool needStrand);
    int32_t compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta);

    int32_t buildSamIndex();

 private:
    int64_t headEndLine;
    uint32_t samLine;
    std::vector<std::vector<int64_t>> contentPos;

    Reference* pRefeGene;

    /*
     * What the ID analysis found in this block's QNAME column (see sam_qname_column.h). The
     * encoder reads it to lay the column out, and the reader reads it to rebuild the column
     * segment by segment - which is why the analysis, not the encoder, holds it.
     */
    IdSplitAnalysis idAnalysis;
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

private:
    uint32_t baseNCount;

    /*
     * A decoded SEQ exception stream for the current block: the character it carries, and
     * its positions, which increase strictly. One stream per character that occurs, so the
     * set of streams follows the data instead of a fixed list of classes - a character that
     * never occurs has no stream and no meta at all (see the layout note on
     * SeqExceptionClass in sam_actuator.cpp).
     *
     * The cursor lives here rather than in a local because it is read per record: a record's
     * refill picks up where the previous record stopped.
     */
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
    /* State for value-domain QNAME sub-streams (compressIdFieldSplit "mode":"intd"): the payload is
       the coded values themselves, so no outer coder is involved and a range decoder walks them one
       value per line. Indexes run parallel to idNumericBufs; other sub-streams leave their entry
       untouched, so the two layouts cannot be confused at reconstruction time. */
    std::vector<RangeCoder> idIntCoders;
    std::vector<id_int::Model> idIntModes;
    std::vector<uint8_t> idIntActive;
    /* The text of a constant QNAME sub-stream (compressIdFieldSplit "mode":"const"): such a segment
       repeats the same bytes in every line, so the encoder stores it once and the reconstruction
       appends it without touching a coder. Empty for every other kind of sub-stream. */
    std::vector<std::string> idConstTexts;
    /* A "dict" sub-stream's texts, indexed by the per-line index it carries (see the layout in
       compressIdFieldSplit). Empty for every other kind of sub-stream. */
    std::vector<std::vector<std::string>> idDictEntries;
    /* A "hexd" sub-stream's per-line payloads are drawn from this alphabet, each character at
       log2(size) bits (see the layout in compressIdFieldSplit). Empty for every other layout. */
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
    /* Decoder for the QUAL column, whichever coder the stream says wrote it; built by
       CoderFactory::makeQualDecoder (see qual_coder_factory.cpp). Null until the QUAL stream
       has been set up. */
    std::shared_ptr<qual_record_decoder> qualDecoder;

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
