/*
 * fastq_actuator.h - Header file for pbgz project
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

#include <string>

#include "io_block.h"
#include "reference.h"
#include "codec_actuator.h"
#include "coder.h"
#include "coder_io.h"

// Forward declaration
class CompressEngine;

#ifdef __SSE4_2__
#include <emmintrin.h>
#endif

#ifdef __AVX2__
#include <immintrin.h>
#endif


enum class CommentType {
    PLUS_ONLY,     // Comment line contains only a '+' sign
    SAME_AS_ID,    // Comment line is the same as ID line
    OTHER,         // Comment line is neither of the above, could be any other string including empty string
    UNKNOWN        // Unknown type
};

typedef struct Mapping {
    void set(uint8_t *_squashBuffer, uint32_t _squashBufferLen, uint8_t *_squashBuffCmplt, uint32_t __squashBuffCmpltLen, uint32_t _offset) {
        squashBuffer[0] = _squashBuffer;
        squashBufferLen[0] = _squashBufferLen;
        squashBuffer[1] = _squashBuffCmplt;
        squashBufferLen[1] = __squashBuffCmpltLen;
        offset = _offset;
    }
	
    inline void incOffset() { ++offset; }

    uint8_t* getSquash(bool pair) const {
        return ((pair) ? (squashBuffer[1] - offset) : (squashBuffer[0] + offset));
    }

    /*  squash buffer, 0 for forward strand, 1 for complementary strand */
    uint8_t *squashBuffer[2];      /* Pointer to squash buffer */
    uint32_t squashBufferLen[2];  /* Data length corresponding to squash buffer */
    uint32_t leftUnalignLen[2]; /* Length of unaligned bases on the left */
    uint8_t leftUnalign[2][3]; /* Unaligned bases on the left, storing squashed data */
    uint32_t rightUnalignLen[2]; /* Length of unaligned bases on the right */
    uint8_t rightUnalign[2][3]; /* Unaligned bases on the right, storing squashed data */
    uint32_t offset;  /* Current offset */
} Mapping;


class FastqCodecActuator : public CodecActuator {
public:
    FastqCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr, Reference* pRef = nullptr): CodecActuator(inPtr, outPtr, engine) {
        idPosLength = 0;
        minBaseLength = INT32_MAX;
        maxBaseLength = 0;
        pReference = pRef;
        baseNCount = 0;
        commentType = CommentType::UNKNOWN;
        baseDecoder = nullptr;
        commentDecoder = nullptr;
        qualityDecoder = nullptr;

        isGen2 = false;
        mappingBuffer = nullptr;
        basePairBuffer = nullptr;
        std::fill_n(basePairSquashBuffer, 4, nullptr);
        std::fill_n(baseSquashBuffer, 4, nullptr);
        baseMappedBuffer = nullptr;
        baseMappedLength = 0;
        baseMappedPosBuffer = nullptr;
        baseMappedPairBuffer = nullptr;
        baseStripNBuffer = nullptr;
        baseNPosBuffer = nullptr;
        refeStretchBuffer = nullptr;
        baseMatchRle = false;
        baseMatchBlock = nullptr;
        baseMatchBlockLen = 0;
        baseMatchBlockOffset = 0;
        baseLitBases = false;
    }

    virtual ~FastqCodecActuator() {
        qualityDecoder.reset();
        baseDecoder.reset();
        commentDecoder.reset();
        idDecoders.clear();

        ioVecters.clear();
        MemoryUtil::safeFree(mappingBuffer);
        MemoryUtil::safeFree(baseMatchBlock);
    }

    int32_t decompress() override ;

    int32_t compress() override ;

    int32_t preAnalysis();

    int32_t initEncoder();

private:

    int32_t preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLe );

    int32_t preAnalysisId(uint8_t* pBuffer, uint32_t bufLen);

    int32_t preAnalysisBase(uint8_t* pBuffer, uint32_t bufLen);

    int32_t preAnalysisComment(uint8_t* pBuffer, uint32_t bufLen, uint32_t lineNo);

    int32_t compressId();

    int32_t compressIdInAll();

    int32_t compressIdInSplit();

    int32_t compressIdStream(coder_io* idIo, coder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx);

    int32_t compressBase();

    int32_t compressBaseWithRef();

    int32_t compressBaseWithoutRef();

    int32_t compressComment();

    int32_t compressQuality();

    int32_t initDecoder(RoughIOBlock* outputBlock);

    /*
     * The reference layout's match stream when it was split into runs and values: rebuild the
     * whole stream from the "m"/"mval" sub-streams into matchBlock, so the per-read pulls below
     * are a straight copy. Answers 0, or -1 with the reason already logged.
     */
    int32_t initRleMatchStream(const Json::Value& metaStreams, uint32_t& id, uint32_t& readOffset);

    /*
     * One read's match bytes: out of the rebuilt stream when this block's was split, or straight
     * from its coder otherwise. Returns the number of bytes (== len), or a negative coder error.
     */
    int32_t pullMatchBytes(uint8_t* dst, uint32_t len);

    /*
     * Decode one auxiliary sub-stream of the reference layout whole, into a caller-supplied buffer
     * of exactly the length the stream declares on its source side. The coder is built from the
     * magic the stream carries and the level is replayed from the same meta, so a stream written
     * by a different coder decodes without this side naming it.
     */
    int32_t decodeWholeStream(const Json::Value& streamMeta, const uint8_t* src, uint32_t srcLen,
                              uint8_t* dst, uint32_t dstLen, const char* tag);


    void mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir);


    // SIMD-optimized N character counting methods
    static inline uint32_t countN_SSE2(const uint8_t* data, size_t length);
    static inline uint32_t countN_Unrolled(const uint8_t* data, size_t length);
    static inline uint32_t countN_Optimized(const uint8_t* data, size_t length);

private:
    std::vector<uint8_t> idSplitSymbols;
    uint32_t idPosLength;
    std::vector<uint16_t> idPositions; // Position of each separator in ID field of each record, starting from 0
    std::vector<uint32_t> idSplitMinLen; // Minimum length of each separator
    std::vector<uint32_t> idSplitMaxLen; // Maximum length of each separator
    const std::string idSplitDefault = "/:= _.,-#\r\t\n";
    uint32_t minBaseLength;
    uint32_t maxBaseLength;
    uint64_t baseNCount;
    CommentType commentType;
    std::vector<std::pair<uint16_t, uint16_t>> qualityFreqTable; // Quality value frequency statistics
    Reference* pReference;

    std::vector<std::shared_ptr<coder>> idDecoders;
    std::vector<std::shared_ptr<coder_io>> ioVecters;
    std::shared_ptr<coder> baseDecoder;
    std::shared_ptr<coder> commentDecoder;
    std::shared_ptr<qual_record_decoder> qualityDecoder;

    bool isGen2;

    /// Variables used in reference genome scenario
    uint8_t* mappingBuffer;
    uint8_t* basePairBuffer;
    uint8_t* basePairSquashBuffer[4];
    uint8_t* baseSquashBuffer[4];
    uint8_t* baseMappedBuffer;
    uint32_t baseMappedLength;
    uint64_t* baseMappedPosBuffer;
    uint8_t* baseMappedPairBuffer;
    uint8_t* baseStripNBuffer;
    uint32_t* baseNPosBuffer;
    uint8_t* refeStretchBuffer;

    /*
     * The auxiliary columns of the reference layout, expanded per block.
     *
     * All three are read one entry per read (the length, the reference position) or per N (the N
     * position), so they are expanded once here rather than parsed inside the per-read loop; the
     * streams themselves are varints, which is what the compressed size depends on. The vectors
     * own the storage the pointers above point into - a fixed uint16/uint32/uint64 per read is no
     * longer written or read (see compressBaseWithRef), and the length vector is in bases, not in
     * "bases above the shortest read", so nothing clamps a read to 65535 any more.
     */
    std::vector<uint32_t> baseNPosVec;
    std::vector<uint32_t> baseLengthVec;
    std::vector<uint64_t> baseMappedPosVec;

    /* Whether this block's match stream arrived split into runs and values (see
       splitSeqMatchStream): when it did, matchBlock holds the rebuilt stream and the per-read
       pulls copy out of it instead of asking the coder. */
    bool baseMatchRle;
    uint8_t* baseMatchBlock;
    uint32_t baseMatchBlockLen;
    uint32_t baseMatchBlockOffset;

    /* Whether this block's "no match" records carry their own characters rather than two-bit codes
       (see the literal form in compressBaseWithRef, and the "litbases" member of the block meta). */
    bool baseLitBases;
};
