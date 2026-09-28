/*
 * fastq_actuator.cpp - Source file for pbgz project
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

#include <cstring>
#include <algorithm>
#include <atomic>
#include <memory>

#include "fastq_actuator.h"
#include "compress_engine.h"
#include "reference.h"
#include "log/logger.h"
#include "utils/memory_util.h"
#include "coder_io.h"
#include "coder_bwt_cm.h"
#include "coder_affix_match.h"
#include "coder_fc.h"
#include "coder_qual.h"
#include "sam_field_layout.h"   /* kSeqExc*Name: the names the reader and both writers must agree on */
#include "seq_stream_util.h"
#include "utils/md5_util.h"
#include "coder_json.h"
#include "actg.h"
#include "city.h"
#include "pbgz_stat.h"
#include "profile_stats.h"
#include "config_manager.h"

namespace {
    void recordFastqFieldStats(PbgzEngine* engine, uint16_t objectId, uint32_t srcLen, uint32_t dstLen) {
        if (!engine) return;

        auto compressEngine = dynamic_cast<CompressEngine*>(engine);
        if (!compressEngine || !compressEngine->getStats()) return;

        auto fastqStat = dynamic_cast<FastqStat*>(compressEngine->getStats());
        if (!fastqStat) return;

        if (objectId != 0 && srcLen > 0) {
            fastqStat->addMetricValue(StatUnitIds::COMPRESSION_RATIO, objectId, StatMetricIds::ORIGINAL_SIZE, srcLen);
            fastqStat->addMetricValue(StatUnitIds::COMPRESSION_RATIO, objectId, StatMetricIds::COMPRESSED_SIZE, dstLen);
        }
    }
}

const uint32_t MAPPED_THRESHOLD_GEN2 = 2;

/*
 * Whether any byte of this read is outside A/C/G/T.
 *
 * Such a read - soft-masked lower case, an IUPAC ambiguity code - cannot be held in the two-bit
 * alphabet at all, so it can be neither coded against the reference nor two-bit coded, and it is
 * written with its own characters instead (see the literal form in compressBaseWithRef). Nothing
 * else on the reference path has to ask about a read's characters.
 *
 * allowN says whether an upper-case N is to be taken as ordinary: it is on a read's own line, where
 * the N's have not been cut out yet, and off once they have - a lower-case n is never cut out and
 * so has to be reported on either side (see compressBaseWithRef).
 */
static bool seqHasNonAcgt(const uint8_t* p, uint32_t len, bool allowN)
{
    uint32_t i = 0;
#ifdef __SSE4_2__
    const __m128i acgt[4] = {_mm_set1_epi8('A'), _mm_set1_epi8('C'),
                             _mm_set1_epi8('G'), _mm_set1_epi8('T')};
    for (; i + 16 <= len; i += 16) {
        const __m128i x = _mm_loadu_si128((const __m128i*)(p + i));
        const __m128i ac = _mm_or_si128(_mm_cmpeq_epi8(x, acgt[0]), _mm_cmpeq_epi8(x, acgt[1]));
        __m128i gt = _mm_or_si128(_mm_cmpeq_epi8(x, acgt[2]), _mm_cmpeq_epi8(x, acgt[3]));
        if (allowN) {
            gt = _mm_or_si128(gt, _mm_cmpeq_epi8(x, _mm_set1_epi8('N')));
        }
        if ((uint32_t)_mm_movemask_epi8(_mm_or_si128(ac, gt)) != 0xFFFFu) {
            return true;
        }
    }
#endif
    for (; i < len; ++i) {
        const uint8_t ch = p[i];
        if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T' && !(allowN && ch == 'N')) {
            return true;
        }
    }
    return false;
}

// SIMD-optimized N character counting using SSE4.2
uint32_t FastqCodecActuator::countN_SSE2(const uint8_t* data, size_t length) {
    uint32_t count = 0;
    size_t i = 0;

#ifdef __SSE4_2__
    const __m128i target_upper = _mm_set1_epi8('N');

    for (; i + 16 <= length; i += 16) {
        __m128i chunk = _mm_loadu_si128((__m128i*)(data + i));
        __m128i cmp_upper = _mm_cmpeq_epi8(chunk, target_upper);
        uint32_t mask = _mm_movemask_epi8(cmp_upper);
        count += __builtin_popcount(mask);
    }
#endif

    for (; i < length; i++) {
        count += (data[i] == 'N');
    }

    return count;
}

// Loop-unrolled version for better branch prediction and reduced loop overhead
uint32_t FastqCodecActuator::countN_Unrolled(const uint8_t* data, size_t length) {
    uint32_t count = 0;
    size_t i = 0;

    for (; i + 8 <= length; i += 8) {
        uint8_t b0 = data[i];
        uint8_t b1 = data[i + 1];
        uint8_t b2 = data[i + 2];
        uint8_t b3 = data[i + 3];
        uint8_t b4 = data[i + 4];
        uint8_t b5 = data[i + 5];
        uint8_t b6 = data[i + 6];
        uint8_t b7 = data[i + 7];

        count += (b0 == 'N');
        count += (b1 == 'N');
        count += (b2 == 'N');
        count += (b3 == 'N');
        count += (b4 == 'N');
        count += (b5 == 'N');
        count += (b6 == 'N');
        count += (b7 == 'N');
    }

    for (; i < length; i++) {
        count += (data[i] == 'N');
    }

    return count;
}

// Auto-selecting wrapper that chooses the best implementation
uint32_t FastqCodecActuator::countN_Optimized(const uint8_t* data, size_t length) {
    if (length < 16) {
        return countN_Unrolled(data, length);
    }

#ifdef __SSE4_2__
    return countN_SSE2(data, length);
#else
    return countN_Unrolled(data, length);
#endif
}

int32_t FastqCodecActuator::compress() {
    /*
     * The same guard as preAnalysis() (which fails a block with no line position), for the
     * callers that compress without running the pre-analysis first. A FASTQ record is four
     * lines, all of them found through the reader's line positions, so a block with bytes but
     * no line position cannot be coded as FASTQ: writing it would produce a block that decodes
     * to nothing while its meta claims the md5 of the real data. Fail instead, which lets the
     * pipeline store the block as binary rather than lose it (CompressEngine::actuatorPreProc).
     */
    if (inBlockPtr->getDataLen() > 0 && inBlockPtr->getNpos().empty()) {
        LOG_ERROR("FASTQ block(%ld) carries %ld bytes but no line positions, not FASTQ text.",
                  (long)inBlockPtr->getBlockId(), (long)inBlockPtr->getDataLen());
        return -1;
    }

    if (0 != initEncoder()) {
        LOG_ERROR("Init encoder failed");
        return -1;
    }

    if (0 != compressId()){
        LOG_ERROR("Compress id failed");
        return -1;
    }

    if (0 != compressBase()) {
        LOG_ERROR("Compress base failed");
        return -1;
    }

    if (0 != compressComment()) {
        LOG_ERROR("Compress comment failed");
        return -1;
    }

    if (0 != compressQuality()) {
        LOG_ERROR("Compress quality failed");
        return -1;
    }

    // Calculate MD5 of data block
    std::string md5;
    calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
    meta["md5"] = md5;
    meta["idlines"] = static_cast<uint32_t>(inBlockPtr->getNpos().size() >> 2);

    // Compress meta information
    coder_json metaCoder;
    int32_t metaLen = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
    if (metaLen <= 0) {
        LOG_ERROR("Failed to encode meta information");
        return -1;
    }
    outBlockPtr->setMetaLen(metaLen);

    return 0;
}

/// @brief Parse separators and their positions in the first line
/// @param pBuffer    // Buffer for ID line
/// @param bufLen     // Length of ID line
/// @param idSplitSymbols    // List of separators
/// @return 0 for success, -1 for failure
int32_t FastqCodecActuator::preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLen) {
    if (pBuffer == nullptr || bufLen == 0) {
        return -1;
    }

    std::vector<uint32_t> idSplitPos;
    for (uint32_t i = 1 ; i < bufLen; ++i) {    // First character is @, skip it
        char ch = pBuffer[i];
        if (idSplitDefault.find(ch) != std::string::npos) {
            idSplitSymbols.push_back(ch);
            /* The positions travel as bytes in the split layouts, so they are kept to what a
               byte holds; a QNAME longer than that is not something those layouts describe. */
            idSplitPos.push_back(static_cast<uint8_t>(i));
        }
    }

    // Initialize max and min length for each separator
    for (size_t i = 0; i < idSplitSymbols.size(); ++i) {
        idSplitMinLen.push_back(UINT32_MAX);
        idSplitMaxLen.push_back(0);
    }

    // Copy first line ID analysis information to idPositions
    uint32_t lastPos = 0;
    for (uint32_t idx = 0;  idx < idSplitPos.size(); ++idx) {
        uint32_t pos = idSplitPos[idx];
        uint32_t curLen = pos - lastPos - (0 == idx ? 0 : 1);   // First line no needs offset
        if (curLen < idSplitMinLen[idx]) {
            idSplitMinLen[idx] = curLen;
        }
        if (curLen > idSplitMaxLen[idx]) {
            idSplitMaxLen[idx] = curLen;
        }
        idPositions.push_back(pos);
        idPosLength++;
        lastPos = pos;
    }
    return 0;
}

int32_t FastqCodecActuator::preAnalysisId(uint8_t* pBuffer, uint32_t bufferLen) {
    uint32_t lastPos = 0;
    uint32_t lastFindPos = 0;
    for (uint32_t idx = 0; idx < idSplitSymbols.size(); ++idx) {
        uint8_t symbol = idSplitSymbols[idx];
        void* found = memchr(pBuffer + lastPos, symbol, bufferLen - lastPos);
        if (found == nullptr) {
            LOG_DEBUG("preAnalysisId meet exception: block(%d) will compress id in all mode", inBlockPtr->getBlockId());
            idPosLength = UINT32_MAX;
            break;
        }

        uint32_t pos = (uint8_t*)found - pBuffer;
        uint32_t curLen = pos - lastFindPos - (0 == idx ? 0 : 1);
        if (curLen < idSplitMinLen[idx]) {
            idSplitMinLen[idx] = curLen;
        }
        if (curLen > idSplitMaxLen[idx]) {
            idSplitMaxLen[idx] = curLen;
        }
        idPositions.push_back(static_cast<uint8_t>(pos));   /* same byte-wide bound as above */
        idPosLength++;
        lastPos = pos + 1;
        lastFindPos = pos;
    }
    return 0;
}

int32_t FastqCodecActuator::preAnalysisBase(uint8_t* pBuffer, uint32_t bufLen) {
    uint32_t baseLength = bufLen - 1;    // Remove newline character
    if (baseLength > maxBaseLength) {
        maxBaseLength = baseLength;
    }
    if (baseLength < minBaseLength) {
        minBaseLength = baseLength;
    }

    // Use SIMD-optimized N counting for better performance
    uint32_t segmentLength = bufLen - 1;
    baseNCount += countN_Optimized(pBuffer, segmentLength);

    return 0;
}

int32_t FastqCodecActuator::preAnalysisComment(uint8_t* pBuffer, uint32_t bufLen, uint32_t lineNo) {
    if (commentType != CommentType::OTHER) {
        if (commentType == CommentType::UNKNOWN) {  // First line
            if (*pBuffer == '+' && bufLen == 2) {
                commentType = CommentType::PLUS_ONLY;
            }
            else {
                // Find ID line content
                uint8_t* idStart = nullptr;
                uint32_t idLineNo = lineNo - 2;
                if (idLineNo == 0) {
                    idStart = inBlockPtr->getBuffer();
                } else {
                    idStart = inBlockPtr->getBuffer() + inBlockPtr->getNpos()[idLineNo - 1] + 1;
                }

                uint32_t idLength = inBlockPtr->getNpos()[idLineNo];
                if ((idLength -1) == bufLen - 1 && 0 == memcmp(idStart + 1, pBuffer, bufLen)) {
                    commentType = CommentType::SAME_AS_ID;
                }else {
                    commentType = CommentType::OTHER;
                }
            }
        } else if (commentType == CommentType::PLUS_ONLY) {
            if (*pBuffer != '+' || bufLen != 2) {
                commentType = CommentType::OTHER;
            }
        } else if (commentType == CommentType::SAME_AS_ID) {
            // Find ID line content
            uint8_t* idStart = nullptr;
            uint32_t idLineNo = lineNo - 2;
            if (idLineNo == 0) {
                idStart = inBlockPtr->getBuffer();
            } else {
                idStart = inBlockPtr->getBuffer() + inBlockPtr->getNpos()[idLineNo - 1] + 1;
            }

            uint32_t idLength = inBlockPtr->getNpos()[idLineNo];
            if ((idLength -1) != bufLen - 1 || 0 != memcmp(idStart + 1, pBuffer, bufLen)) {
                commentType = CommentType::OTHER;
            }
        }
    }
    return 0;
}

int32_t FastqCodecActuator::preAnalysis() {
    uint64_t lineNum = inBlockPtr->getNpos().size();
    if (lineNum == 0) {
        LOG_ERROR("line number is zero");
        return -1;
    }

    uint64_t startPos = 0;
    std::pair<uint8_t, uint32_t> qualityFrequnce[256];
    for (int i = 0; i < 256; ++i) {
        qualityFrequnce[i].first = i;
        qualityFrequnce[i].second = 0;
    }
    for (uint32_t lineNo = 0; lineNo < lineNum; ++lineNo) {
        uint64_t endPos = inBlockPtr->getNpos()[lineNo];
        uint32_t lineLength  = endPos - startPos + 1;   // Length needs to include newline character
        switch (lineNo & 0x3)
        {
        case 0: {  // ID line
            if (idPosLength != UINT32_MAX) {
                if (lineNo == 0) {  // First line
                    if (preAnalysisIdFirstLine(inBlockPtr->getBuffer() + startPos, lineLength) != 0) {
                        return -1;
                    }
                } else {
                    if (preAnalysisId(inBlockPtr->getBuffer() + startPos, lineLength) != 0) {
                        return -1;
                    }
                }
            }
            break;
        }
        case 1: { //  Base line
            if (preAnalysisBase(inBlockPtr->getBuffer() + startPos, lineLength) != 0) {
                return -1;
            }
            break;
        }
        case 2: { // Comment line
            if (preAnalysisComment(inBlockPtr->getBuffer() + startPos, lineLength, lineNo) != 0) {
                return -1;
            }
            break;
        }
        case 3: { // Quality line
            for (uint32_t idx = startPos; idx < endPos; ++idx) {
                qualityFrequnce[*(inBlockPtr->getBuffer() + idx)].second++;
            }
            break;
        }
        default:
            return -1;
        }
        startPos = endPos + 1;
    }

    std::sort(qualityFrequnce, qualityFrequnce + 256,
        [](const std::pair<uint8_t, uint32_t> &a, const std::pair<uint8_t, uint32_t> &b){ return a.second > b.second; });
    for (int i = 0; i < 256; ++i) {
        if (qualityFrequnce[i].second == 0) {
            continue;
        }
        qualityFreqTable.push_back(std::make_pair(qualityFrequnce[i].first - '!', 1));
    }
    LOG_DEBUG("minBaseLen=%d, maxBaseLen=%d", minBaseLength, maxBaseLength);
    return 0;
}

int32_t FastqCodecActuator::initEncoder() {
    const uint32_t line4 = (inBlockPtr->getNpos().size() >> 2);
    const uint32_t lmax = inBlockPtr->getMaxLineLen() + 4; /* 4 reserved for base key not 4-aligned */
    const uint32_t lsquash = (lmax >> 2) + !!(lmax & 0x3); /* squash length */

    isGen2 = (inBlockPtr->getBlockType() == FASTQ_GEN2) || (inBlockPtr->getBlockType() == FASTQ_GEN2_GZIP);
    baseMappedLength = (isGen2) ? lmax : (lmax << 1);

    if (pReference) {
       /* base pair + 4 base pair squash + 4 base squash + base mapped + base N pos in block +
        * base delete N + mapped pos + mapped pair each line
        *
        * There is no fixed-width base-length area any more: the lengths go out as varints, built
        * while the reads are walked (see compressBaseWithRef), so this layout only holds what is
        * written per read rather than per read and per length.
        */
        uint32_t n = lmax + (lsquash << 3) + baseMappedLength + (baseNCount << 2);
        n += lmax + (line4 << 3) + line4;
        mappingBuffer = MemoryUtil::safeAlloc<uint8_t>(n);
        uint8_t *p = mappingBuffer;
        basePairBuffer = p;
        p += lmax;
        for (n = 0; n < 4; n++) {
            basePairSquashBuffer[n] = p;
            p += lsquash;
            baseSquashBuffer[n] = p;
            p += lsquash;
        }

        baseMappedBuffer = p;
        p += baseMappedLength;

        baseNPosBuffer = (uint32_t *)p;
        p += (baseNCount << 2);

        baseStripNBuffer = p;
        p += lmax;

        baseMappedPosBuffer = (uint64_t *)p;
        p += (line4 << 3);

        baseMappedPairBuffer = p;
        p += line4;

        /* Gen3 has no mapper of its own: the base path routes it to compression without a
           reference (see compressBase), so there is nothing here for it to point at. */
    }

    return 0;
}

void FastqCodecActuator::mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir) {
    Mapping mappingTable[4];
    uint8_t ch;
    const uint32_t baseGroupLen = pReference->getBaseGroupLength();
    const uint32_t bgMid = baseGroupLen >> 1;
    const uint32_t mapThresh = MAPPED_THRESHOLD_GEN2;
    const uint8_t* pSeq[2] = {base, basePairBuffer};
    const uint8_t bgIsUnalign4 = !!(baseGroupLen & 0x3);
    const uint32_t lenBgs = (baseGroupLen >> 2) + bgIsUnalign4;
    const uint32_t baseSquashAlign4 = (baseLength >> 2) + !!(baseLength & 0x3);
    uint8_t* prefSquash = (uint8_t*)(pReference->getSquash());
    const int64_t refeSquashLen = pReference->getSquashLength();
    uint32_t bestPosInRefe = UINT32_MAX;
    uint32_t bestIsPair = 0;
    uint32_t bestUnmatches = UINT32_MAX;
    uint32_t bestAlign4 = 0;
    uint32_t bestPos = UINT32_MAX;
    const uint64_t xSquashTab[2] = {0xFCFFFFFFFFFFFFFF, 0xFCFFFFFFFFFFFFFF};

    /* case 1: base length is not greater than reference index corresponding base length */
    if (baseLength <= (baseGroupLen + 4)) {
        actgEncode(base, out, baseLength);
        outLength = baseLength;
        mappingPos = 0;
        mappingDir = 2; /* During decompression, check mdir first; if it's 2, it means no match */
        return;
    }

    /* case 2: base length is greater than reference index corresponding base length */
    const uint32_t edge = baseLength - baseGroupLen;
    const auto mapSquash0 = pbgzprof::nowOrZero();
    actgPair(basePairBuffer, base, baseLength);

    /*
     * Calculate align4 squash buffer and pair squash buffer. Unlike the other timers around it,
     * this one is restarted as each round begins (see the two assignments below), so it is the one
     * marker here that is not const.
     */
    auto mapPrep0 = pbgzprof::nowOrZero();
    for (uint32_t n = 0; n < 4; n++) {
        uint32_t squashLength[2];
        squashLength[0] = (baseLength - n) >> 2;
        uint32_t total = squashLength[0] << 2;
        mappingTable[n].leftUnalignLen[0] = n;
        for (uint32_t m = 0; m < n; m++) {
            mappingTable[n].leftUnalign[0][m] = ((*(pSeq[0] + m)) >> 1) & 0x3; /* squash value of bases not 4-aligned on the left */
        }
        mappingTable[n].rightUnalignLen[0] = baseLength - n - total;
        for (uint32_t m = 0; m < mappingTable[n].rightUnalignLen[0]; m++) {
            mappingTable[n].rightUnalign[0][m] = ((*(pSeq[0] + baseLength - mappingTable[n].rightUnalignLen[0] + m)) >> 1) & 0x3;
        }
        actgSquash(pSeq[0] + n, total, baseSquashBuffer[n]);

        squashLength[1] = (baseLength - n + 1) >> 2; /* Add a character to the right for 32-byte alignment; need to handle the last character when matching to mappingTable[0]'s pair and mappingTable[0] offset is 0 */
        total = squashLength[1] << 2;
        mappingTable[n].leftUnalignLen[1] = baseLength + 1 - n - total;
        for (uint32_t m = 0; m < mappingTable[n].leftUnalignLen[1]; m++) {
            mappingTable[n].leftUnalign[1][m] = ((*(pSeq[1] + m)) >> 1) & 0x3;
        }
        mappingTable[n].rightUnalignLen[1] = (n == 0) ? 0 : (n - 1); /* Subtract 1 because one character is added to the right for key alignment */
        for (uint32_t m = 0; m < mappingTable[n].rightUnalignLen[1]; m++) {
            mappingTable[n].rightUnalign[1][m] = ((*(pSeq[1] + baseLength - mappingTable[n].rightUnalignLen[1] + m)) >> 1) & 0x3;
        }
        actgSquash(pSeq[1] + mappingTable[n].leftUnalignLen[1], total, basePairSquashBuffer[n]);

        /* Establish the relationship between base squash and corresponding pair base squash */
        mappingTable[n].set(baseSquashBuffer[n], squashLength[0], basePairSquashBuffer[n] + squashLength[1] - lenBgs, squashLength[1], 0);

        /* do mapping */
        uint32_t align4Curr = n & 0x3;
        uint32_t matchPairOrigin = (*(pSeq[0] + n + bgMid) < *(pSeq[1] + baseLength - baseGroupLen - n + bgMid));

        uint8_t* pSquash = mappingTable[align4Curr].getSquash(matchPairOrigin);
        uint64_t xSquash = *((uint64_t *)(pSquash));
        xSquash &= xSquashTab[matchPairOrigin];
        uint32_t hash32 = (uint32_t)CityHash64((const char *)(&xSquash), lenBgs);
        uint32_t posCnts;
        uint32_t* posVals = (uint32_t *)(pReference->queryPosition(hash32, posCnts));
        pbgzprof::addSince(pbgzprof::MAP_PREP, mapPrep0);
        pbgzprof::add(pbgzprof::MAP_ATTEMPTS, 0);
        const auto mapCand0 = pbgzprof::nowOrZero();

        for (uint32_t o = 0; o < posCnts; o++) {
            uint32_t matchPos = *posVals++;
            uint32_t matchPair = ((matchPos & 0x80000000) >> 31) ^ matchPairOrigin;
            matchPos = (matchPos & 0x7FFFFFFF) << 3; /* to squash reference pos */

            /* check left and right boundary simply */
            if (matchPos + baseSquashAlign4 >= refeSquashLen || matchPos < baseSquashAlign4) {
                continue;
            }

            uint32_t lOffset = (matchPair) ? (mappingTable[align4Curr].squashBufferLen[1] - mappingTable[align4Curr].offset - lenBgs) : (mappingTable[align4Curr].offset);

            /* caculate unmatch count */
            pSquash = mappingTable[align4Curr].getSquash(matchPair) - lOffset;
            uint8_t* pSquashRefe = prefSquash + matchPos - lOffset;

            uint64_t xSquashMatch = ((*((uint64_t *)(mappingTable[align4Curr].getSquash(matchPair)))) & xSquashTab[matchPair]);
            uint64_t xSquashMatchRefe = ((*((uint64_t *)(prefSquash + matchPos))) & xSquashTab[matchPair]);
            if (xSquashMatch != xSquashMatchRefe)  {
                /* key is not same, skip */
                continue;
            }
            /* align 4 */
            uint32_t unmatches = actgSquashDiffCnt(pSquash, pSquashRefe, mappingTable[align4Curr].squashBufferLen[matchPair]);
            if (unmatches >= bestUnmatches) {
                continue;
            }

            /*  Because this case adds a character to the right for 32-byte alignment: need to handle the last character when matching to mappingTable[0]'s pair and mappingTable[0] offset is 0 */
            unmatches -= (matchPair && (mappingTable[align4Curr].offset == 0)) ? ((xSquashMatch & 0x80000000000000) != (xSquashMatchRefe & 0x80000000000000)) : 0;

            /* left unalign */
            uint32_t l = mappingTable[align4Curr].leftUnalignLen[matchPair];
            uint32_t m;
            for (ch = *(pSquashRefe - 1), m = 0; m < mappingTable[align4Curr].leftUnalignLen[matchPair]; m++) {
                unmatches += ((ch >> (m << 1)) & 0x3) != (mappingTable[align4Curr].leftUnalign[matchPair][l - 1]);
                l--;
            }
            /* right unalign */
            l = mappingTable[align4Curr].squashBufferLen[matchPair];
            for (ch = *(pSquashRefe + l), m = 0; m < mappingTable[align4Curr].rightUnalignLen[matchPair]; m++) {
                unmatches += ((ch >> (6 - (m << 1)) & 0x3) != (mappingTable[align4Curr].rightUnalign[matchPair][m]));
            }
            if (unmatches < bestUnmatches) {
                bestPos = (matchPos << 2) - (lOffset << 2) - mappingTable[align4Curr].leftUnalignLen[matchPair];
                bestPosInRefe = matchPos - lOffset;
                bestIsPair = matchPair;
                bestAlign4 = align4Curr;
                bestUnmatches = unmatches;
            }
            if (bestUnmatches <= mapThresh) {
                break;
            }
        }
        pbgzprof::addSince(pbgzprof::MAP_CAND, mapCand0);
        if (bestUnmatches <= mapThresh) {
            break;
        }
        mappingTable[align4Curr].incOffset();
        mapPrep0 = pbgzprof::nowOrZero();
    }
    pbgzprof::addSince(pbgzprof::MAP_SQUASH, mapSquash0);

    const auto mapWindow0 = pbgzprof::nowOrZero();
    if (bestUnmatches > mapThresh) { /* continue mapping */
        for (uint32_t n = 4; n <= edge; n++) {
            /* do mapping */
            uint32_t align4Curr = n & 0x3;
            uint32_t matchPairOrigin = (*(pSeq[0] + n + bgMid) < *(pSeq[1] + baseLength - baseGroupLen - n + bgMid));

            uint8_t* pSquash = mappingTable[align4Curr].getSquash(matchPairOrigin);
            uint64_t xSquash = *((uint64_t *)(pSquash));
            xSquash &= xSquashTab[matchPairOrigin];
            uint32_t hash32 = (uint32_t)CityHash64((const char *)(&xSquash), lenBgs);
            uint32_t posCnts;
            uint32_t* posVals = (uint32_t *)(pReference->queryPosition(hash32, posCnts));
            pbgzprof::addSince(pbgzprof::MAP_PREP, mapPrep0);
            pbgzprof::add(pbgzprof::MAP_ATTEMPTS, 0);
            const auto mapCand1 = pbgzprof::nowOrZero();

            for (uint32_t o = 0; o < posCnts; o++) {
                uint32_t matchPos = *posVals++;
                uint32_t matchPair = ((matchPos & 0x80000000) >> 31) ^ matchPairOrigin;
                matchPos = (matchPos & 0x7FFFFFFF) << 3; /* to squash reference pos */

                /* check left and right boundary simply */
                if (matchPos + baseSquashAlign4 >= refeSquashLen || matchPos < baseSquashAlign4) {
                    continue;
                }

                uint32_t lOffset = (matchPair) ? (mappingTable[align4Curr].squashBufferLen[1] - mappingTable[align4Curr].offset - lenBgs) : (mappingTable[align4Curr].offset);
                /* caculate unmatch count */
                pSquash = mappingTable[align4Curr].getSquash(matchPair) - lOffset;
                uint8_t* pSquashRefe = prefSquash + matchPos - lOffset;

                uint64_t xSquashMatch = ((*((uint64_t *)(mappingTable[align4Curr].getSquash(matchPair)))) & xSquashTab[matchPair]);
                uint64_t xSquashMatchRefe = ((*((uint64_t *)(prefSquash + matchPos))) & xSquashTab[matchPair]);
                if (xSquashMatch != xSquashMatchRefe) { /* key is not same, skip */
                    continue;
                }

                /* align 4 */
                uint32_t unmatches = actgSquashDiffCnt(pSquash, pSquashRefe, mappingTable[align4Curr].squashBufferLen[matchPair]);
                if (unmatches >= bestUnmatches) {
                    continue;
                }

                /*  Because this case adds a character to the right for 32-byte alignment: need to handle the last character when matching to mappingTable[0]'s pair and mappingTable[0] offset is 0 */
                unmatches -= (matchPair && (mappingTable[align4Curr].offset == 0)) ? ((xSquashMatch & 0x80000000000000) != (xSquashMatchRefe & 0x80000000000000)) : 0;

                /* left unalign */
                uint32_t l = mappingTable[align4Curr].leftUnalignLen[matchPair];
                uint32_t m;
                for (ch = *(pSquashRefe - 1), m = 0; m < mappingTable[align4Curr].leftUnalignLen[matchPair]; m++) {
                    unmatches += ((ch >> (m << 1)) & 0x3) != (mappingTable[align4Curr].leftUnalign[matchPair][l - 1]);
                    l--;
                }
                /* right unalign */
                l = mappingTable[align4Curr].squashBufferLen[matchPair];
                for (ch = *(pSquashRefe + l), m = 0; m < mappingTable[align4Curr].rightUnalignLen[matchPair]; m++) {
                    unmatches += ((ch >> (6 - (m << 1)) & 0x3) != (mappingTable[align4Curr].rightUnalign[matchPair][m]));
                }
                if (unmatches < bestUnmatches) {
                    bestPos = (matchPos << 2) - (lOffset << 2) - mappingTable[align4Curr].leftUnalignLen[matchPair];
                    bestPosInRefe = matchPos - lOffset;
                    bestIsPair = matchPair;
                    bestAlign4 = align4Curr;
                    bestUnmatches = unmatches;
                }
                if (bestUnmatches <= mapThresh) {
                    break;
                }
            }

            pbgzprof::addSince(pbgzprof::MAP_CAND, mapCand1);
            if (bestUnmatches <= mapThresh) {
                break;
            }
            mappingTable[align4Curr].incOffset();
            mapPrep0 = pbgzprof::nowOrZero();
        }
    }
    pbgzprof::addSince(pbgzprof::MAP_WINDOW, mapWindow0);

    /* calc the result of base or base pair mapping with the match pos reference */
    const auto mapPayload0 = pbgzprof::nowOrZero();
    outLength = 0;
    if (bestUnmatches != UINT32_MAX)  {
        /* get pos in reference table */
        uint8_t* pSquash = (bestIsPair) ? (mappingTable[bestAlign4].squashBuffer[1] - (mappingTable[bestAlign4].squashBufferLen[1] - lenBgs)) : (mappingTable[bestAlign4].squashBuffer[0]);
        uint8_t* pSquashRefe = prefSquash + bestPosInRefe;

        uint32_t n, o, l;
        o = l = mappingTable[bestAlign4].leftUnalignLen[bestIsPair];
        for (ch = *(pSquashRefe - 1), n = 0; n < mappingTable[bestAlign4].leftUnalignLen[bestIsPair]; n++) {
            out[--o] = (((ch >> (n << 1)) & 0x3) ^ (mappingTable[bestAlign4].leftUnalign[bestIsPair][l - 1]));
            l--;
        }
        outLength += n; /* Note byte order */
        outLength += actgStretchMappingXor(pSquash, pSquashRefe, mappingTable[bestAlign4].squashBufferLen[bestIsPair], out + outLength);

        outLength -= (bestIsPair && bestAlign4 == 0);
        l = mappingTable[bestAlign4].squashBufferLen[bestIsPair];
        for (ch = *(pSquashRefe + l), n = 0; n < mappingTable[bestAlign4].rightUnalignLen[bestIsPair]; n++) {
            out[outLength++] = ((ch >> (6 - (n << 1)) & 0x3) ^ (mappingTable[bestAlign4].rightUnalign[bestIsPair][n]));
        }
    } else {
        /* not match valid pos */
        actgEncode(base, out, baseLength);
        outLength = baseLength;
        bestPos = 0;
        bestIsPair = 2; ///*  During decompression, check mdir first; if it's 2, it means no match */
    }

    pbgzprof::addSince(pbgzprof::MAP_PAYLOAD, mapPayload0);
    mappingPos = bestPos;
    mappingDir = bestIsPair;
    return;
}


/*
 * One ID sub-stream, fed line by line. The coder is the base class now: which one it is, is
 * decided by the field (see compressIdInSplit), and every candidate the FASTQ ID row lists is
 * fed the same way here.
 */
int32_t FastqCodecActuator::compressIdStream(coder_io* idIo, coder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx) {
    int32_t currLineOffset = 0;    // Line offset, start of each line
    uint8_t * data = nullptr;
    int32_t currLen = 0;
    srcDataLen = 0;   // Total length of source content
    for (uint32_t idx = 0; idx < inBlockPtr->getNpos().size(); idx += 4) {
        uint32_t splitStep = (idx / 4) * idSplitSymbols.size();
        if (splitSymIdx == 0) {
            data = inBlockPtr->getBuffer() + currLineOffset;
            currLen =  idPositions[splitSymIdx + splitStep] + 1;  //  currIdPos;
        } else {
            data = inBlockPtr->getBuffer() + currLineOffset + idPositions[splitSymIdx + splitStep - 1] + 1;
            currLen = idPositions[splitSymIdx + splitStep] - idPositions[splitSymIdx + splitStep - 1];
        }

        idCoder->encode_line(data, currLen);
        srcDataLen += currLen;
        currLineOffset = inBlockPtr->getNpos()[idx + 3] + 1;
    }
    idCoder->encode_flush();
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idIo->data_len);

    Json::Value tmpMeta;
    tmpMeta["srclen"] = srcDataLen;
    tmpMeta["dstlen"] = idIo->data_len;
    tmpMeta["coder"] = idIo->meta;

    streamMeta.append(tmpMeta);
    return 0;
}

int32_t FastqCodecActuator::compressIdInAll() {
    Json::Value idMeta;
    Json::Value streamMeta;
    std::shared_ptr<coder_io> idIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "ID");
    std::shared_ptr<coder> idCoder = makeFastqFieldEncoder(FQ_ID, CoderType::BWT_CM, idIo.get(), true);

    uint32_t srcDataLen = 0;
    int32_t startPos = 0;
    for (uint32_t lineId = 0; lineId < inBlockPtr->getNpos().size(); lineId = lineId + 4) {
        srcDataLen += inBlockPtr->getNpos()[lineId] - startPos + 1;
        idCoder->encode_line(inBlockPtr->getBuffer() + startPos, inBlockPtr->getNpos()[lineId] - startPos + 1);
        startPos = inBlockPtr->getNpos()[lineId + 3] + 1;
    }
    idCoder->encode_flush();

    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idIo->data_len);
    Json::Value tmpMeta;
    tmpMeta["srclen"] = srcDataLen;
    tmpMeta["dstlen"] = idIo->data_len;
    tmpMeta["coder"] = idIo->meta;
    streamMeta.append(tmpMeta);

    idMeta["totalsrclen"] = srcDataLen;
    idMeta["totaldstlen"] = idIo->data_len;
    idMeta["splitsym"] = std::string("\n");
    idMeta["streams"] = streamMeta;
    meta["id"] = idMeta;
    LOG_INFO("Compress Id in all: from %d to %d, compress ratio: %.2f%%", srcDataLen, idIo->data_len, ((float)(idIo->data_len * 100)) / srcDataLen);

    // Record statistics for ID compression
    recordFastqFieldStats(pbgzEngine, StatObjectId::FASTQ_ID, srcDataLen, idIo->data_len);

    return (size_t)outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
}

int32_t FastqCodecActuator::compressIdInSplit() {
    Json::Value idMeta;
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    for (uint32_t i = 0; i < idSplitSymbols.size();++i) {
        // Variable length with all digits
        LOG_DEBUG("split symbol(%d, %c), minlen = %d, maxlen = %d.", i, idSplitSymbols[i], idSplitMinLen[i], idSplitMaxLen[i]);
        if (idSplitMaxLen[i] != idSplitMinLen[i]) {
            uint32_t currIdPos = idPositions[i];
            uint32_t currLineOffset = 0;
            // Check if all digits based on first line
            uint8_t * data = nullptr;
            uint32_t currLen = 0;
            if (i == 0) {
                data = inBlockPtr->getBuffer() + currLineOffset + 1;
                currLen = currIdPos - 1;
            } else {
                data = inBlockPtr->getBuffer() + currLineOffset + idPositions[i - 1] + 1;
                currLen = currIdPos - idPositions[i - 1] - 1;
            }

            bool idDigit = true;
            // Optimized digit check using bit operation: digits 0x30-0x39 all have high nibble 0x03
            for (uint32_t j = 0; j < currLen; ++j) {
                if ((data[j] & 0xF0) != 0x30) {
                    idDigit = false;
                    break;
                }
            }

            if (idDigit) {
                LOG_DEBUG("Is all digit(%d), sub-stream starts from coder_bwt_cm.", i);
                std::shared_ptr<coder_io> idIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "ID sub-stream");
                std::shared_ptr<coder> idCoder = makeFastqFieldEncoder(FQ_ID, CoderType::BWT_CM, idIo.get(), true);
                uint32_t srcLength = 0;
                int32_t ret= compressIdStream(idIo.get(), idCoder.get(), streamMeta, srcLength, i);
                if (ret != 0) {
                    LOG_ERROR("Failed to compress ID stream, splitid = %d", i);
                    return -1;
                }
                totalSrcLength += srcLength;
                totalDstLength += idIo->data_len;
                continue;
            }
        }

        // Fixed length or not all digits scenario
        LOG_DEBUG("Not all digit or fix length(%d), sub-stream starts from coder_affix_match.", i);
        std::shared_ptr<coder_io> idIoAM = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "ID sub-stream");
        std::shared_ptr<coder> idCoderAm = makeFastqFieldEncoder(FQ_ID, CoderType::AFFIX_MATCH, idIoAM.get(), true);
        uint32_t srcLength = 0;
        int32_t ret= compressIdStream(idIoAM.get(), idCoderAm.get(), streamMeta, srcLength, i);
        if (ret != 0) {
            LOG_ERROR("Failed to compress ID stream, splitid = %d", i);
            return -1;
        }
        totalSrcLength += srcLength;
        totalDstLength += idIoAM->data_len;
    }

    idMeta["totalsrclen"] = totalSrcLength;
    idMeta["totaldstlen"] = totalDstLength;
    idMeta["splitsym"] = std::string((char*)idSplitSymbols.data(), idSplitSymbols.size());
    idMeta["streams"] = streamMeta;
    meta["id"] = idMeta;
    LOG_INFO("Compress Id in split: from %d to %d, compress ratio: %.2f%%.", totalSrcLength, totalDstLength, ((float)(totalDstLength * 100)) / totalSrcLength);

    // Record statistics for ID compression
    recordFastqFieldStats(pbgzEngine, StatObjectId::FASTQ_ID, totalSrcLength, totalDstLength);

    return (size_t)outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
}

int32_t FastqCodecActuator::compressId() {
    PBGZ_PROF_SCOPE(pbgzprof::FIELD_BASE);
    if (idPosLength != UINT32_MAX) {
        for (uint32_t idx = 0; idx < idSplitSymbols.size(); ++idx) {
            if (idSplitMinLen[idx] == 0) {
                idPosLength = UINT32_MAX;
                break;
            }
        }
    }

    // ID check invalid, compress as whole block
    if (idPosLength == UINT32_MAX) {
        return compressIdInAll();
    }else {
        // ID check valid, compress in split blocks
        return compressIdInSplit();
    }
}

int32_t FastqCodecActuator::compressBase() {
    if (pReference != nullptr && isGen2) {
        return compressBaseWithRef();
    } else {
        /* Once per run, not once per block - and atomically, since the blocks that reach this line
           run on several threads at once. */
        static std::atomic<bool> noticePrinted{false};
        if (!isGen2 && pReference != nullptr && !noticePrinted.exchange(true)) {
            fprintf(stderr ,"This file is Gen3 FASTQ, compression will be performed without a reference genome.\n");
        }
        return compressBaseWithoutRef();
    }
}

int32_t FastqCodecActuator::compressBaseWithRef() {
    bool encBaseLen = (minBaseLength != maxBaseLength);
    uint32_t line = inBlockPtr->getNpos().size();
    uint8_t* ptr = inBlockPtr->getBuffer();
    int64_t nOffset = 0;
    int64_t currPos = 0;
    uint32_t offset = 0;

    // Prepare SIMD targets for N detection
#ifdef __SSE4_2__
    const __m128i target_upper = _mm_set1_epi8('N');
#endif
    Json::Value metaSubs;
    Json::Value metaStreams;
    Json::Value metaBase;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;
    /*
     * Whether any read in this block carries a character the 2-bit alphabet cannot hold - anything
     * but A/C/G/T once the N's are taken out, which in practice means soft-masked (lower case)
     * reads or IUPAC ambiguity codes.
     *
     * The whole block is asked once, and only a block that answers yes pays for a per-read look:
     * such a read cannot be coded against the reference and cannot be 2-bit coded either, so it is
     * written with its own characters (see the literal form below and the "litbases" member the
     * decoder dispatches on).
     */
    bool anyLiteral = false;
    for (uint32_t i = 1; i < line; i += 4) {
        const uint32_t endPos = inBlockPtr->getNpos()[i];
        const uint32_t startPos = inBlockPtr->getNpos()[i - 1] + 1;
        if (seqHasNonAcgt(ptr + startPos, endPos - startPos, true)) {
            anyLiteral = true;
            break;
        }
    }
    bool blockHasLiteral = false;
    uint32_t literalReads = 0;
    const auto seqPrep0 = pbgzprof::nowOrZero();
    /*
     * The match stream is gathered whole before it is coded: whether it pays to replace it with
     * zero runs and the values that follow them is a property of the completed stream, which
     * cannot be known while the reads are still being mapped (see splitSeqMatchStream). The
     * buffer is sized by the block's own byte count, which bounds the stripped base count above.
     */
    std::vector<uint8_t> matchBuf;
    matchBuf.reserve(inBlockPtr->getDataLen());
    /* The per-read lengths, as varints, written only when they vary (see below). */
    std::vector<uint8_t> baseLenBuf;
    for (uint32_t i = 1; i < line; i += 4) {
        uint32_t endPos = inBlockPtr->getNpos()[i];
        uint32_t startPos = inBlockPtr->getNpos()[i - 1] + 1;
        uint8_t* pBuff = baseStripNBuffer;
        uint32_t outLen = 0;

        const auto strip0 = pbgzprof::nowOrZero();
        // Use SIMD-optimized N detection when processing base sequences
        size_t segmentLength = endPos - startPos;
        const uint8_t* segmentStart = ptr + startPos;
        size_t simdIdx = 0;

#ifdef __SSE4_2__
        // Process 16 bytes at a time with SIMD
        for (; simdIdx + 16 <= segmentLength; simdIdx += 16) {
            __m128i chunk = _mm_loadu_si128((__m128i*)(segmentStart + simdIdx));
            __m128i cmp = _mm_cmpeq_epi8(chunk, target_upper);
            uint32_t mask = _mm_movemask_epi8(cmp);

            // Process each byte in the 16-byte chunk
            for (int j = 0; j < 16; j++) {
                if (mask & (1 << j)) {
                    // Found N at position simdIdx + j
                    *(baseNPosBuffer + nOffset) = currPos + simdIdx + j;
                    nOffset++;
                } else {
                    // Non-N character, copy to output buffer
                    *pBuff = segmentStart[simdIdx + j];
                    pBuff++;
                }
            }
        }
#endif

        // Process remaining bytes with scalar code
        for (; simdIdx < segmentLength; simdIdx++) {
            uint8_t ch = segmentStart[simdIdx];
            if (ch == 'N') {
                *(baseNPosBuffer + nOffset) = currPos + simdIdx;
                nOffset++;
            } else {
                *pBuff = ch;
                pBuff++;
            }
        }

        const uint32_t strippedLen = (uint32_t)(pBuff - baseStripNBuffer);
        pbgzprof::addSince(pbgzprof::SEQ_STRIP, strip0);
        if (anyLiteral && seqHasNonAcgt(baseStripNBuffer, strippedLen, false)) {
            /*
             * This read's characters cannot survive the 2-bit alphabet, so it is written as it
             * stands and marked as not coming from the reference (direction 2 below).
             */
            memcpy(baseMappedBuffer, baseStripNBuffer, strippedLen);
            outLen = strippedLen;
            baseMappedPosBuffer[offset] = 0;
            baseMappedPairBuffer[offset] = 2;   /* 2 = no match, the decoder's own marker */
            blockHasLiteral = true;
            literalReads++;
        } else {
            const auto map0 = pbgzprof::nowOrZero();
            mappingFastqGen2(baseStripNBuffer, strippedLen,
                             baseMappedBuffer, outLen, baseMappedPosBuffer[offset], baseMappedPairBuffer[offset]);
            pbgzprof::addSince(pbgzprof::SEQ_MAP, map0);
            if (anyLiteral && baseMappedPairBuffer[offset] == 2) {
                /*
                 * An ordinary read the mapper could not place. It is written literally as well,
                 * not 2-bit: once a block holds one literal read the flag is on for all of them,
                 * and the decoder has only that flag to go by (see initDecoder). The read's
                 * characters are the same either way, so nothing is lost by the switch.
                 */
                memcpy(baseMappedBuffer, baseStripNBuffer, strippedLen);
                outLen = strippedLen;
                blockHasLiteral = true;
                literalReads++;
            }
        }
        matchBuf.insert(matchBuf.end(), baseMappedBuffer, baseMappedBuffer + outLen);
        pReference->updateMatchedGene(baseMappedPosBuffer[offset], outLen);
        if (encBaseLen) {
            /*
             * The length stream goes out as one varint per read, of its base count above the
             * shortest read of the file. The fixed uint16 per read this replaces cost two bytes
             * whatever the spread, and silently truncated a read longer than minBaseLength +
             * 65535; a varint is one byte for the usual spread and has no ceiling.
             */
            appendVarint(baseLenBuf, endPos - startPos - minBaseLength);
        }
        offset++;
        currPos += endPos - startPos;
    }

    pbgzprof::addSince(pbgzprof::SEQ_PREP, seqPrep0);
    const auto seqCodec0 = pbgzprof::nowOrZero();
    /*
     * Sub-stream "m": the run-length half of the split stream, or the whole match stream when the
     * split does not pay. Under the split the surviving values follow as sub-stream "mval" (see
     * seq_stream_util.h for what makes the split worth taking); the decoder tells the two layouts
     * apart by the "rle" member here, so a block written without the split stays byte-for-byte
     * what it was.
     */
    const SeqRleSplit rle = splitSeqMatchStream(matchBuf.data(), (uint32_t)matchBuf.size());
    const uint32_t matchSrcLen = (uint32_t)matchBuf.size();
    {
        const uint32_t payLen = rle.useRle ? rle.runLength : matchSrcLen;
        const uint8_t* payBuf = rle.useRle ? rle.run.get() : matchBuf.data();
        std::shared_ptr<coder_io> matchIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match");
        CoderFactory::applyLevel(matchIo.get(), CoderType::BWT_CM, engineCompressLevel());
        std::shared_ptr<coder_bwt_cm> matchCm = std::make_shared<coder_bwt_cm>(matchIo.get());
        matchCm->encode_line(const_cast<uint8_t*>(payBuf), payLen);
        matchCm->encode_flush();
        if (matchIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ match stream overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)payLen;
        metaSubs["dstlen"] = (Json::Value::Int)(matchIo->data_len);
        metaSubs["coder"] = matchIo->meta;
        metaSubs["sname"] = "m";
        if (rle.useRle) {
            /* orgrawlen = the match stream before the split, which is what the decoder rebuilds;
               rlelen is this sub-stream's decoded length, i.e. the "srclen" above. */
            metaSubs["rle"] = (Json::Value::UInt)1;
            metaSubs["orgrawlen"] = (Json::Value::UInt)matchSrcLen;
            metaSubs["rlelen"] = (Json::Value::UInt)rle.runLength;
        }
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + matchIo->data_len);
        totalDstLen += matchIo->data_len;
        totalSrcLen += matchSrcLen;
        LOG_DEBUG("SEQ match stream: matchLen=%u nonZero=%u (%.2f%%) rle=%d runLen=%u valLen=%u src=%u dst=%u",
                  matchSrcLen, rle.nonZero,
                  (double)rle.nonZero * 100.0 / (double)(matchSrcLen ? matchSrcLen : 1),
                  (int)rle.useRle, rle.runLength, rle.valLength, payLen,
                  (uint32_t)matchIo->data_len);
    }

    /* Sub-stream "mval": the values that survived the split, one byte each. */
    if (rle.useRle && rle.valLength > 0) {
        std::shared_ptr<coder_io> valIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match val");
        CoderFactory::applyLevel(valIo.get(), CoderType::BWT_CM, engineCompressLevel());
        std::shared_ptr<coder_bwt_cm> valCm = std::make_shared<coder_bwt_cm>(valIo.get());
        valCm->encode_line(rle.val.get(), rle.valLength);
        valCm->encode_flush();
        if (valIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ match value stream overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)rle.valLength;
        metaSubs["dstlen"] = (Json::Value::Int)valIo->data_len;
        metaSubs["coder"] = valIo->meta;
        metaSubs["sname"] = "mval";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + valIo->data_len);
        totalDstLen += valIo->data_len;
    }

    /* Second sub-stream: position stream of matches between reads and reference */
    std::shared_ptr<coder_io> posIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ mapped pos");
    uint32_t line4 = inBlockPtr->getNpos().size() >> 2;
    int64_t srcLen = 0;
    /*
     * Where each read was matched in the reference, as forward deltas in varint, zigzagged because
     * a FASTQ is not sorted by position and a read may well sit before the previous one (the sign
     * goes into the low bit, see seq_stream_util.h).
     *
     * The deltas pay only where they carry less than the positions do. The layout they replace is
     * one uint64 per read, eight bytes of which the three or four high ones are always zero -
     * nearly free for the coder, since the zeros repeat at a fixed period - so a column of
     * unrelated deltas can come out *larger* than the positions themselves once both are coded.
     * Measured: with positions running in order the deltas are worth 3.3x on the column and 2.9x on
     * the whole archive, while with the mapper finding no match for most reads (positions mostly 0
     * and occasionally far away) they cost about 1% of the archive instead.
     *
     * Neither the order nor the size of the delta stream separates the two - 7% of the reads go
     * backwards in the first case against 10% in the second, and the delta column comes out at 25%
     * of the absolute one against 22% - so the choice is made on what the two would weigh in
     * varints, byte for byte, which is what the coders are handed: the deltas only when the column
     * of them is the smaller of the two.
     */
    std::vector<uint8_t> posBuf;
    posBuf.reserve((size_t)line4 * 4 + 16);
    uint64_t posDeltaVarBytes = 0;
    uint64_t posAbsVarBytes = 0;
    {
        uint64_t prevPos = 0;
        for (uint32_t i = 0; i < line4; ++i) {
            const uint64_t pos = baseMappedPosBuffer[i];
            const uint64_t zz = zigzag64((int64_t)(pos - prevPos));
            appendVarint64(posBuf, zz);
            posDeltaVarBytes += varintLen64(zz);
            posAbsVarBytes += varintLen64(pos);
            prevPos = pos;
        }
    }
    const uint32_t posFixedLen = (line4 << 3);
    const bool posUseDelta = (line4 > 0) && (posDeltaVarBytes < posAbsVarBytes);
    uint8_t* posPayBuf = posUseDelta ? posBuf.data() : (uint8_t*)baseMappedPosBuffer;
    srcLen = posUseDelta ? (int64_t)posBuf.size() : (int64_t)posFixedLen;
    if (srcLen > FC_MIN_LEN && srcLen < FC_MAX_LEN) {
        std::shared_ptr<coder_fc> posCm = std::make_shared<coder_fc>(posIo.get());
        posCm->encode_line(posPayBuf, (uint32_t)srcLen);
        posCm->encode_flush();
    } else {
        std::shared_ptr<coder_bwt_cm> posCm = std::make_shared<coder_bwt_cm>(posIo.get());
        CoderFactory::applyLevel(posIo.get(), CoderType::BWT_CM, engineCompressLevel());
        posCm->encode_line(posPayBuf, (uint32_t)srcLen);
        posCm->encode_flush();
    }
    metaSubs.clear();
    metaSubs["srclen"] =  (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] =  (Json::Value::Int)posIo->data_len;
    metaSubs["coder"] = posIo->meta;
    metaSubs["sname"] = "mpos";
    if (posUseDelta) {
        /* The delta layout, and how many reads it covers. A stream without these holds the
           absolute uint64 per read the layout above writes, which is also what an archive written
           before the deltas carries. */
        metaSubs["delta"] = (Json::Value::UInt)1;
        metaSubs["count"] = (Json::Value::UInt)line4;
    }
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + posIo->data_len);
    totalSrcLen += (uint32_t)srcLen;
    totalDstLen += posIo->data_len;
    LOG_DEBUG("SEQ mapped pos stream: reads=%u delta=%d src=%u dst=%u (delta varints %llu, absolute varints %llu)",
              line4, (int)posUseDelta, (uint32_t)srcLen, (uint32_t)posIo->data_len,
              (unsigned long long)posDeltaVarBytes, (unsigned long long)posAbsVarBytes);

    /* Third sub-stream: pair identifier stream of matches between reads and reference */
    auto pairIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ mapped pair");
    srcLen = line4;
    if (srcLen > FC_MIN_LEN && srcLen < FC_MAX_LEN) {
        std::shared_ptr<coder_fc> subCoder = std::make_shared<coder_fc>(pairIo.get());
        subCoder->encode_line((uint8_t *)baseMappedPairBuffer, srcLen);
        subCoder->encode_flush();
    } else {
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(pairIo.get());
        CoderFactory::applyLevel(pairIo.get(), CoderType::BWT_CM, engineCompressLevel());
        subCoder->encode_line((uint8_t *)baseMappedPairBuffer, srcLen);
        subCoder->encode_flush();
    }
    metaSubs.clear();
    metaSubs["srclen"] = (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] = (Json::Value::Int)pairIo->data_len;
    metaSubs["coder"] = pairIo->meta;
    metaSubs["sname"] = "mpair";
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + pairIo->data_len);
    totalSrcLen += srcLen;
    totalDstLen += pairIo->data_len;
    LOG_DEBUG("SEQ mapped pair stream: reads=%u src=%u dst=%u", line4, (uint32_t)srcLen,
              (uint32_t)pairIo->data_len);

    /* Fourth sub-stream: positions of all N's in reads */
    if (baseNCount > 0) {
        std::shared_ptr<coder_io> nposIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ npos");
        /*
         * Three ways to describe the positions, which only ever grow within a block, and the
         * smallest source wins - the stream's own name says which one was taken (the decoder
         * dispatches on it, so no marker is needed):
         *   "nposd" forward deltas in varint, one entry per N - the general form, one byte per N
         *           as long as consecutive N's are near each other;
         *   "nposr" runs: one gap and one length per run of consecutive positions, which is what
         *           poly-N stretches come out as (a run of 35 costs two varints instead of 35);
         *   "npos"  the absolute four-byte offsets of the first layout, kept for the case where
         *           neither of the others is smaller (it never wins by much, being four bytes per
         *           N, but it is the form the old archives carry and stays readable).
         * All three are built here rather than measured: they are a few bytes per N, and the
         * choice is not close.
         */
        std::vector<uint8_t> deltaBuf;
        std::vector<uint8_t> gapsBuf;
        std::vector<uint8_t> lensBuf;
        deltaBuf.reserve((size_t)baseNCount + 16);
        gapsBuf.reserve((size_t)baseNCount / 4 + 16);
        lensBuf.reserve((size_t)baseNCount / 4 + 16);
        uint32_t prevPos = 0;
        uint32_t runStart = 0;
        uint32_t runLen = 0;
        uint32_t prevRunEnd = 0;
        for (uint64_t n = 0; n < baseNCount; ++n) {
            const uint32_t pos = baseNPosBuffer[n];
            appendVarint(deltaBuf, pos - prevPos);
            if (runLen != 0 && pos == prevPos + 1) {
                runLen++;
            } else {
                if (runLen != 0) {
                    appendVarint(gapsBuf, runStart - prevRunEnd);
                    appendVarint(lensBuf, runLen);
                    prevRunEnd = prevPos + 1;
                }
                runStart = pos;
                runLen = 1;
            }
            prevPos = pos;
        }
        if (runLen != 0) {
            appendVarint(gapsBuf, runStart - prevRunEnd);
            appendVarint(lensBuf, runLen);
        }
        const uint32_t absBytes = (uint32_t)(baseNCount << 2);
        const uint32_t deltaBytes = (uint32_t)deltaBuf.size();
        /* The run form is one list of gaps followed by one list of lengths, as SAM's is. */
        const uint32_t runBytes = (uint32_t)(gapsBuf.size() + lensBuf.size());

        const char* nposName = kSeqExcDeltaName;
        uint8_t* nposPay = deltaBuf.data();
        uint32_t nposPayLen = deltaBytes;
        if (runBytes < nposPayLen) {
            gapsBuf.insert(gapsBuf.end(), lensBuf.begin(), lensBuf.end());
            nposName = kSeqExcRunName;
            nposPay = gapsBuf.data();
            nposPayLen = runBytes;
        }
        if (absBytes < nposPayLen) {
            nposName = kSeqExcAbsName;
            nposPay = (uint8_t*)baseNPosBuffer;
            nposPayLen = absBytes;
        }
        srcLen = (int64_t)nposPayLen;
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(nposIo.get());
        CoderFactory::applyLevel(nposIo.get(), CoderType::BWT_CM, engineCompressLevel());
        subCoder->encode_line(nposPay, nposPayLen);
        subCoder->encode_flush();
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)srcLen;
        metaSubs["dstlen"] = (Json::Value::Int)nposIo->data_len;
        metaSubs["coder"] = nposIo->meta;
        metaSubs["sname"] = nposName;
        /* How many entries the stream holds: positions for the first two forms, runs for "nposr".
           A "npos" stream without it is the same absolute form, written before this member
           existed; a "npos" with "delta" is the delta form written before it had its own name. */
        metaSubs["count"] = (Json::Value::UInt)(strcmp(nposName, "nposr") == 0
                                                    ? lensBuf.size()
                                                    : (size_t)baseNCount);
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + nposIo->data_len);
        totalSrcLen += (uint32_t)srcLen;
        totalDstLen += nposIo->data_len;
        LOG_DEBUG("SEQ npos stream: %llu N -> %s src=%u dst=%u (delta %u, run %u, absolute %u)",
                  (unsigned long long)baseNCount, nposName, nposPayLen,
                  (uint32_t)nposIo->data_len, deltaBytes, runBytes, absBytes);
    } else {
        LOG_INFO("NCount is 0, base info not contains npos stream");
    }
    if (anyLiteral) {
        LOG_DEBUG("SEQ literal form: %u of %u reads in this block carry characters outside A/C/G/T",
                  literalReads, line >> 2);
    }
    metaBase["ncount"] = (Json::Value::UInt)baseNCount;
    /*
     * The literal form: set when this block holds at least one record written with its own
     * characters rather than two-bit codes, and then *every* record whose direction says "no
     * match" carries its own characters (see the walk above). One member for the whole block is
     * all the decoder needs, since that is exactly the set of records it applies to.
     */
    if (blockHasLiteral) {
        metaBase["litbases"] = (Json::Value::UInt)1;
    }

    /* Fifth stream: length of each base line, one varint per read (see the walk above). */
    if (encBaseLen) {
        std::shared_ptr<coder_io> lenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ length");
        srcLen = (int64_t)baseLenBuf.size();
        std::shared_ptr<coder_bwt_cm> sub_coder = std::make_shared<coder_bwt_cm>(lenIo.get());
        CoderFactory::applyLevel(lenIo.get(), CoderType::BWT_CM, engineCompressLevel());
        sub_coder->encode_line(baseLenBuf.data(), (uint32_t)srcLen);
        sub_coder->encode_flush();
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)srcLen;
        metaSubs["dstlen"] = (Json::Value::Int)lenIo->data_len;
        metaSubs["coder"] = lenIo->meta;
        metaSubs["sname"] = "baselen";
        /* The varint layout, and how many reads it covers. A stream without these holds the fixed
           uint16 per read an older archive holds. */
        metaSubs["delta"] = (Json::Value::UInt)1;
        metaSubs["count"] = (Json::Value::UInt)line4;
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + lenIo->data_len);
        totalSrcLen += (uint32_t)srcLen;
        totalDstLen += lenIo->data_len;
    } else {
        LOG_INFO("Base length is same, base info not contains baselen stream");
    }
    metaBase["minlen"] = (Json::Value::UInt)minBaseLength;
    metaBase["maxlen"] = (Json::Value::UInt)maxBaseLength;
    metaBase["totalsrclen"] = (Json::Value::UInt)totalSrcLen;
    metaBase["totaldstlen"] = (Json::Value::UInt)totalDstLen;
    metaBase["streams"] = metaStreams;
    meta["base"] = metaBase;

    pbgzprof::addSince(pbgzprof::SEQ_CODEC, seqCodec0);
    LOG_INFO("Compress base with reference: from %d to %d, compress ratio: %.2f%%.", totalSrcLen, totalDstLen, ((float)(totalDstLen * 100)) / totalSrcLen);

    // Record statistics for base compression
    recordFastqFieldStats(pbgzEngine, StatObjectId::FASTQ_BASE, totalSrcLen, totalDstLen);

    return 0;
}

int32_t FastqCodecActuator::compressBaseWithoutRef() {
    int32_t addLength = !!(minBaseLength != maxBaseLength);
    uint32_t baseSrcLength = 0;
    for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx +=4) {
        int64_t end = inBlockPtr->getNpos()[idx];
        int64_t start = inBlockPtr->getNpos()[idx-1] + 1;
        baseSrcLength += end - start + addLength;
    }

    /*
     * Two ways of feeding the column, and the coder decides which: coder_fc wants the whole
     * column as one block - it is gathered into scratch space at the end of the output buffer
     * first - while the others are fed line by line straight from the input block.
     *
     * coder_fc also only takes a column inside its size window. A column outside it falls back to
     * coder_bwt_cm even when that is what preprocessing picked, because the size window is a
     * property of the coder and not something the verdict can override - the verdict was reached
     * on a sample that was inside it.
     */
    CoderType baseCoderType = pickedFastqCoder(FQ_SEQ, CoderType::BWT_CM);
    const bool fcSizeOk = (baseSrcLength > FC_MIN_LEN && baseSrcLength < FC_MAX_LEN);
    std::shared_ptr<coder_io> baseIo;

    if (baseCoderType == CoderType::FC && fcSizeOk) {
        uint8_t* tmpBase = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize() - baseSrcLength;
        uint32_t srcLength = 0;
        for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx += 4) {
            int64_t end = inBlockPtr->getNpos()[idx];
            int64_t start = inBlockPtr->getNpos()[idx-1] + 1;
            int32_t lineLength = end - start + addLength;
            memcpy(tmpBase + srcLength, inBlockPtr->getBuffer() + start, lineLength);
            srcLength += lineLength;
        }

        baseIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain() - baseSrcLength, "SEQ");
        std::shared_ptr<coder> baseCoder = makeEncoderOfType(baseCoderType, baseIo.get(), false,
                                                             CoderType::BWT_CM);
        baseCoder->encode_line(tmpBase, srcLength);
        baseCoder->encode_flush();
    } else {
        if (baseCoderType == CoderType::FC) {
            baseCoderType = CoderType::BWT_CM;
        }
        baseIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ");
        std::shared_ptr<coder> baseCoder = makeEncoderOfType(baseCoderType, baseIo.get(), true,
                                                             CoderType::BWT_CM);

        for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx += 4) {
            int64_t end = inBlockPtr->getNpos()[idx];
            int64_t start = inBlockPtr->getNpos()[idx - 1] + 1;
            int32_t lineLength = end - start + addLength;
            baseCoder->encode_line(inBlockPtr->getBuffer() + start, lineLength);
        }
        baseCoder->encode_flush();
    }

    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + baseIo->data_len);
    Json::Value baseMeta;
    baseMeta["minlen"] = minBaseLength;
    baseMeta["maxlen"] = maxBaseLength;
    baseMeta["coder"] = baseIo->meta;
    baseMeta["totalsrclen"] = baseSrcLength;
    baseMeta["totaldstlen"] = baseIo->data_len;
    meta["base"] = baseMeta;

    LOG_INFO("Compress base: from %d to %d, compress ratio: %.2f%%.", baseSrcLength, baseIo->data_len, ((float)(baseIo->data_len * 100)) / baseSrcLength);

    // Record statistics for base compression
    recordFastqFieldStats(pbgzEngine, StatObjectId::FASTQ_BASE, baseSrcLength, baseIo->data_len);

    return 0;
}

int32_t FastqCodecActuator::compressComment() {
    std::shared_ptr<coder_io> commentIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "comment");
    std::shared_ptr<coder> commentCoder = makeFastqFieldEncoder(FQ_COMMENT, CoderType::BWT_CM, commentIo.get(), true);

    Json::Value commentMeta;
    switch (commentType) {
    case CommentType::PLUS_ONLY:
        commentMeta["type"] = "plusonly";
        LOG_INFO("Compress comment: plus only.");
        break;
    case CommentType::SAME_AS_ID:
        commentMeta["type"] = "sameasid";
        LOG_INFO("Compress comment: same as id.");
        break;
    case CommentType::OTHER:{
        commentMeta["type"] = "other";
        int32_t commmentSrcLength = 0;
        for (uint32_t idx = 2; idx < inBlockPtr->getNpos().size(); idx += 4) {
            int32_t end = inBlockPtr->getNpos()[idx];
            int32_t start = inBlockPtr->getNpos()[idx - 1] + 1;
            int32_t lineLength = end - start;
            commentCoder->encode_line(inBlockPtr->getBuffer() + start, lineLength + 1);
            commmentSrcLength += lineLength;
        }
        commentCoder->encode_flush();
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + commentIo->data_len);

        Json::Value subMeta;
        commentMeta["srclen"] = commmentSrcLength;
        commentMeta["dstlen"] = commentIo->data_len;
        commentMeta["coder"] = commentIo->meta;

        LOG_INFO("Compress comment: from %d to %d, compress ratio: %.2f%%.", commmentSrcLength, commentIo->data_len, ((float)(commentIo->data_len * 100)) / commmentSrcLength);
        break;
    }
    default:
        break;
    }
    meta["comment"] = commentMeta;

    return 0;
}

/*
 * The QUAL column, through the same dedicated path SAM's QUAL column takes.
 *
 * A FASTQ record's quality is not just a byte stream: coder_qual codes it against the read's
 * bases, and coder_fcv2 against its sequencing cycle, so both have to be handed a record at a
 * time and both have to be trialled on records - which is what the FASTQ QUAL row of
 * kFastqFieldCoderConfig and QualSelector do. Which of them won is read off the verdict here and
 * turned into an encoder by the factory; all this side then does is hand over records (see
 * qual_record_encoder).
 *
 * When the mode's row runs no trial - fast leaves its list empty - there is no verdict to read,
 * and the column keeps the row's default, coder_qual: still a record-level coder, just not one
 * chosen by comparison.
 *
 * What a FASTQ record does not have is a FLAG, so the strand is always "forward": the cycle
 * fcv2 restores comes from the record's length, which a FASTQ gives it just as an aligned SAM
 * does.
 */
int32_t FastqCodecActuator::compressQuality() {
    std::shared_ptr<coder_io> qualityIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL");

    const uint8_t mode = (pbgzEngine != nullptr) ? pbgzEngine->getParameter().mode
                                                : (uint8_t)PBGZ_MODE_ARCHIVE;
    CoderType pickedQualCoder = fastqFieldDefaultCoder(FQ_QUAL, mode, CoderType::QUAL);
    const PreprocessInfo* preInfo = (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
    if (preInfo != nullptr) {
        pickedQualCoder = preInfo->coderFor(FQ_QUAL, pickedQualCoder);
    }

    /*
     * The prior is trained from the first blocks and written once per file; its absolute
     * address travels in the stream's meta so the decoding side can retrieve the same snapshot.
     * It is registered only when a snapshot was actually loaded, so the decoder is never
     * promised one it cannot get.
     */
    AuxPayloadPtr qualPriorBlob = (pbgzEngine != nullptr) ? pbgzEngine->getQualPrior(0) : AuxPayloadPtr();
    bool qualPriorLoaded = false;
    int64_t qualPriorAddress = (pbgzEngine != nullptr) ? pbgzEngine->getQualPriorAddress() : -1;

    const FieldCodecSelection* qualSel = (preInfo != nullptr) ? preInfo->getField(FQ_QUAL) : nullptr;
    QualCoderArgs qualCoderArgs;
    qualCoderArgs.freqTable = qualityFreqTable.empty() ? nullptr : &qualityFreqTable;
    qualCoderArgs.fcv2Params =
        (qualSel != nullptr && qualSel->selectedCoder == CoderType::FCV2) ? &qualSel->fcv2Params
                                                                          : nullptr;
    qualCoderArgs.priorBlob = qualPriorBlob.get();
    qualCoderArgs.priorLoaded = &qualPriorLoaded;

    const auto qualCtor0 = pbgzprof::nowOrZero();
    std::shared_ptr<qual_record_encoder> qualityCoder =
        CoderFactory::makeQualEncoder(pickedQualCoder, qualityIo.get(), qualCoderArgs,
                                      engineCompressLevel());
    pbgzprof::addSince(pbgzprof::QUAL_CTOR, qualCtor0);

    /*
     * The prior carries the alphabet of the blocks it was trained on, and a stream is coded with
     * one alphabet from beginning to end. A later block whose quality values left that alphabet
     * cannot be coded at all (see coder_io::IO_UNCODABLE), so rather than fail the block the prior
     * is dropped for it: the coder then builds its alphabet from this block's own values, which is
     * what the decoder does for a block with no prior. The prior's address is reported per block,
     * so the decoder follows either way.
     */
    if (pickedQualCoder == CoderType::FCV2 && qualPriorBlob.get() != nullptr && qualityCoder != nullptr) {
        bool covered = true;
        for (size_t i = 0; i < qualityFreqTable.size(); ++i) {
            if (!qualityCoder->coversByte((uint8_t)(qualityFreqTable[i].first + '!'))) {
                covered = false;
                break;
            }
        }
        if (!covered) {
            LOG_INFO("QUAL prior does not cover block %lld's quality values, coding it without the prior.",
                     (long long)inBlockPtr->getBlockId());
            qualCoderArgs.priorBlob = nullptr;
            std::shared_ptr<coder_io> plainIo =
                makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL");
            qualityCoder = CoderFactory::makeQualEncoder(pickedQualCoder, plainIo.get(), qualCoderArgs,
                                                         engineCompressLevel());
            qualityIo = plainIo;
            qualPriorLoaded = false;
        }
    }
    if (!qualPriorLoaded) {
        qualPriorAddress = -1;
    }
    const auto qualCodec0 = pbgzprof::nowOrZero();

    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;
    Json::Value streamMeta;

    // Encode quality data
    uint32_t streamSrcLen = 0;
    for (uint32_t idx = 3; idx < inBlockPtr->getNpos().size(); idx += 4) {
        uint32_t end = inBlockPtr->getNpos()[idx];
        uint32_t start = inBlockPtr->getNpos()[idx - 1] + 1;
        uint32_t lineLength = end - start;
        /* The record's own base line is the context; quality and bases are the same length in
         * FASTQ, so the two lengths are one here. */
        qualityCoder->encode_record(inBlockPtr->getBuffer() + start, lineLength,
                                    inBlockPtr->getBuffer() + inBlockPtr->getNpos()[idx - 3] + 1,
                                    lineLength, false);
        streamSrcLen += lineLength;
    }
    qualityCoder->flush();
    pbgzprof::addSince(pbgzprof::QUAL_CODEC, qualCodec0);
    /*
     * A stream that did not fit the space left in the block is truncated, and a stream that met a
     * quality value its alphabet does not hold (coder_io::IO_UNCODABLE) is one symbol short. Both
     * leave the decoder reading a stream whose tail it cannot reconstruct - it either runs off the
     * end or decodes the rest of the block from the wrong offset - so neither is written out here.
     */
    if (qualityIo->err != coder_io::IO_OK) {
        LOG_ERROR("QUAL stream is not complete (err=%d, wrote %u bytes, block %lld); "
                  "refusing to write a stream that cannot be decoded",
                  (int)qualityIo->err, (unsigned)qualityIo->data_len,
                  (long long)inBlockPtr->getBlockId());
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityIo->data_len);

    Json::Value subMeta;
    subMeta["srclen"] = streamSrcLen;
    subMeta["dstlen"] = qualityIo->data_len;
    subMeta["coder"] = qualityIo->meta;
    subMeta["streamname"] = "qual";
    if (qualPriorAddress >= 0) {
        subMeta["prior"] = (Json::Value::Int64)qualPriorAddress;
    }
    streamMeta.append(subMeta);

    totalSrcLength += streamSrcLen;
    totalDstLength += qualityIo->data_len;

    // Encode quality frequency table
    std::shared_ptr<coder_io> qualityFreqIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL freq table");
    std::shared_ptr<coder_bwt_cm> qualityFreqCoder = std::make_shared<coder_bwt_cm>(qualityFreqIo.get());
    CoderFactory::applyLevel(qualityFreqIo.get(), CoderType::BWT_CM, engineCompressLevel());
    std::shared_ptr<uint16_t[]> qualityFreqArray = std::make_unique<uint16_t[]>(qualityFreqTable.size()<< 1);
    for (uint32_t i = 0; i < qualityFreqTable.size(); ++i) {
        int idx = i << 1;
        qualityFreqArray[idx] = qualityFreqTable[i].first;
        qualityFreqArray[idx + 1] = qualityFreqTable[i].second;
    }

    uint32_t freqSrcLen = (qualityFreqTable.size() << 1) * sizeof(uint16_t);
    qualityFreqCoder->encode_line((uint8_t*)qualityFreqArray.get(), freqSrcLen);

    qualityFreqCoder->encode_flush();
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityFreqIo->data_len);

    subMeta.clear();
    subMeta["srclen"] = freqSrcLen;
    subMeta["dstlen"] = qualityFreqIo->data_len;
    subMeta["coder"] = qualityFreqIo->meta;
    subMeta["streamname"] = "qualityfreq";
    streamMeta.append(subMeta);

    totalSrcLength += freqSrcLen;
    totalDstLength += qualityFreqIo->data_len;

    Json::Value qualityMeta;
    qualityMeta["totalsrclen"] = totalSrcLength;
    qualityMeta["totaldstlen"] = totalDstLength;
    qualityMeta["streams"] = streamMeta;
    meta["quality"] = qualityMeta;

    LOG_INFO("Compress quality: from %d to %d, compress ratio: %.2f%%.", totalSrcLength, totalDstLength, ((float)(totalDstLength * 100)) / totalSrcLength);

    // Record statistics for quality compression
    recordFastqFieldStats(pbgzEngine, StatObjectId::FASTQ_QUALITY, totalSrcLength, totalDstLength);

    return 0;
}

int32_t FastqCodecActuator::initDecoder(RoughIOBlock* outputBlock) {
    /* The actuator is reused across blocks, so the match stream's layout is a per-block fact and
       has to be cleared here: only the layout that says so is served from a rebuilt buffer. */
    baseMatchRle = false;
    baseLitBases = false;
    baseNPosVec.clear();
    baseLengthVec.clear();
    Json::Value& idStreamMeta = meta["id"]["streams"];
    if (idStreamMeta.size() != idSplitSymbols.size()) {
        LOG_ERROR("id stream meta size not match id split symbols size(%d, %d)", idStreamMeta.size(), idSplitSymbols.size());
        return -1;
    }

    uint32_t readOffset = 0;
    /*
     * Initialize ID decoder. The coder is built from the magic the stream carries rather than
     * named here: which coder won the trial is a per-file decision, and a stream says which one
     * it is (see CoderFactory::makeFieldDecoder). The level is replayed from the same meta, so
     * the decoder is on the model the encoder used.
     */
    for (uint32_t idx = 0; idx < idStreamMeta.size(); ++idx) {
        std::string coderName = idStreamMeta[idx]["coder"]["magic"].asString();
        uint32_t dstLen = idStreamMeta[idx]["dstlen"].asUInt();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "ID sub-stream");
        ioVecters.push_back(io);

        FieldDecoderArgs args;
        if (idStreamMeta[idx]["coder"].isMember("level")) {
            args.level = idStreamMeta[idx]["coder"]["level"].asInt();
        }
        std::shared_ptr<coder> decoder = CoderFactory::makeFieldDecoder(coderName, io.get(), args);
        if (decoder == nullptr) {
            LOG_ERROR("unsupported coder name: %s", coderName.c_str());
            return -1;
        }
        idDecoders.push_back(decoder);
        readOffset += dstLen;
    }

    // Initialize Base decoder
    Json::Value& baseMeta = meta["base"];
    minBaseLength = baseMeta["minlen"].asUInt();
    maxBaseLength = baseMeta["maxlen"].asUInt();
    baseNCount = baseMeta["ncount"].asUInt();
    bool isUseReference = (pReference != nullptr) && baseMeta.isMember("streams");
    if (!isUseReference) {
        std::string baseCoderName = baseMeta["coder"]["magic"].asString();
        uint32_t srcBaseLen = baseMeta["totalsrclen"].asUInt();
        uint32_t dstBaseLen = baseMeta["totaldstlen"].asUInt();
        if (baseCoderName == "coder_fc") {
            coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstBaseLen, &ioErrSink, "SEQ");
            baseIo.meta = baseMeta;
            baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
            coder_fc baseDecoderFc(&baseIo);
            baseDecoderFc.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcBaseLen, srcBaseLen,  UINT8_MAX, false);
        } else {
            /*
             * Every other coder is built from the magic, as the ID column builds its
             * decoders: the stream says which coder wrote it, and the level comes from
             * the same meta. The column is then decoded line by line below.
             */
            std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstBaseLen, "SEQ");
            ioVecters.push_back(io);

            FieldDecoderArgs args;
            if (baseMeta["coder"].isMember("level")) {
                args.level = baseMeta["coder"]["level"].asInt();
            }
            baseDecoder = CoderFactory::makeFieldDecoder(baseCoderName, io.get(), args);
            if (baseDecoder == nullptr) {
                LOG_ERROR("unsupported coder name: %s", baseCoderName.c_str());
                return -1;
            }
        }
        readOffset += dstBaseLen;
    } else {
        Json::Value metaStreams = baseMeta["streams"];
        /* The literal form: when the block carries this member, a record whose direction says "no
           match" holds its own characters rather than two-bit codes (see compressBaseWithRef). */
        baseLitBases = baseMeta.isMember("litbases");
        uint32_t maxBaseLen = baseMeta["maxlen"].asUInt();
        uint32_t n = maxBaseLen;  /* sub stream 1, strip n  */
        n += maxBaseLen; /* for somebuffer_refe_stretch */
        for (uint32_t id = 1; id < metaStreams.size(); id++) {
            n += metaStreams[id]["srclen"].asUInt();
        }
        mappingBuffer = MemoryUtil::safeAlloc<uint8_t>(n);
        uint8_t* ps = mappingBuffer;
        baseStripNBuffer = ps;
        ps += maxBaseLen;
        refeStretchBuffer = ps;
        ps += maxBaseLen;

        /* check sub streams 1 */
        uint32_t id = 0;
        if (metaStreams[id]["sname"].asString() != "m") {
            LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
            return -1;
        }
        uint32_t baseDstLen = metaStreams[id]["dstlen"].asUInt();
        /*
         * Two layouts, told apart by the "rle" member the split one carries (see
         * splitSeqMatchStream): a match stream split into runs and values, which is rebuilt whole
         * here and then sliced per read, or the stream itself, which the coder hands out read by
         * read. Blocks written before the split existed carry no "rle" and decode as they always
         * did.
         */
        if (metaStreams[id].isMember("rle")) {
            if (initRleMatchStream(metaStreams, id, readOffset) != 0) {
                return -1;
            }
        } else {
            if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, baseDstLen, "SEQ match");
                ioVecters.push_back(io);
                baseDecoder = std::make_shared<coder_bwt_cm>(io.get());
            } else {
                LOG_ERROR("check sub stream failed: coder name unmatch");
                return -1;
            }
            readOffset += baseDstLen;
        }

        /* check sub streams 2 */
        id++;
        baseMappedPosBuffer = nullptr;
        uint8_t* temBuffer = nullptr;
        uint32_t srcLen = 0;
        uint32_t dstLen = 0;
        if (metaStreams[id]["sname"].asString() != "mpos") {
            LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
            return -1;
        }
        /*
         * The positions are expanded here, once per block, like the other two auxiliary columns:
         * the per-read walk asks for one at a time, and the varint delta stream is only smaller
         * while it stays in that form. Two layouts are read - the deltas this build writes (the
         * "delta" member), or the absolute uint64 per read of an archive written before them.
         */
        baseMappedPosVec.clear();
        temBuffer = inBlockPtr->getBuffer() + readOffset;
        srcLen = metaStreams[id]["srclen"].asUInt();
        dstLen = metaStreams[id]["dstlen"].asUInt();
        {
            std::vector<uint8_t> raw(srcLen + 1);
            if (decodeWholeStream(metaStreams[id], temBuffer, dstLen, raw.data(), srcLen,
                                  "SEQ mapped pos") != 0) {
                return -1;
            }
            if (metaStreams[id].isMember("delta")) {
                uint32_t off = 0;
                uint64_t prev = 0;
                while (off < srcLen) {
                    uint64_t zz = 0;
                    if (!readVarint64(raw.data(), srcLen, off, zz)) {
                        LOG_ERROR("mpos stream truncated in block %llu", meta["block_id"].asInt64());
                        return -1;
                    }
                    prev = (uint64_t)((int64_t)prev + unzigzag64(zz));
                    baseMappedPosVec.push_back(prev);
                }
            } else {
                if (srcLen & 7) {
                    LOG_ERROR("mpos stream is %u bytes, not a whole number of uint64", srcLen);
                    return -1;
                }
                const uint32_t reads = srcLen >> 3;
                baseMappedPosVec.resize(reads);
                for (uint32_t i = 0; i < reads; ++i) {
                    uint64_t v = 0;
                    for (uint32_t b = 0; b < 8; ++b) {
                        v |= (uint64_t)raw[i * 8 + b] << (b * 8);
                    }
                    baseMappedPosVec[i] = v;
                }
            }
        }
        baseMappedPosBuffer = baseMappedPosVec.data();

        readOffset += metaStreams[id]["dstlen"].asUInt();

        /* check sub stream 3 */
        id++;
        baseMappedPairBuffer = nullptr;
        if (metaStreams[id]["sname"].asString() != "mpair") {
            LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
            return -1;
        }
        if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPairBuffer = ps;
            ps += srcLen;
            coder_io pairIo(temBuffer, dstLen, &ioErrSink, "SEQ mapped pair");
            auto pairCm = std::make_unique<coder_bwt_cm>(&pairIo);
            pairCm->decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else if (metaStreams[id]["coder"]["magic"].asString() == "coder_fc") {
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPairBuffer = ps;
            ps += srcLen;
            coder_io pairIo(temBuffer, dstLen, &ioErrSink, "SEQ mapped pair");
            pairIo.meta = metaStreams[id];
            coder_fc pairFc(&pairIo);
            pairFc.decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else{
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }
        readOffset += metaStreams[id]["dstlen"].asUInt();

        /* check sub stream 4 */
        baseNPosBuffer = nullptr;
        baseNPosVec.clear();
        if (baseMeta["ncount"].asUInt()) {
            id++;
            /* Three names for the same column, one per form it is written in (see
               compressBaseWithRef): the deltas, the runs, and the absolute list of the first
               layout. "npos" also covers the delta form written before it had a name of its own,
               which is recognised by its "delta" member below. */
            const std::string nposStreamName = metaStreams[id]["sname"].asString();
            if (nposStreamName != kSeqExcAbsName && nposStreamName != kSeqExcDeltaName && nposStreamName != kSeqExcRunName) {
                LOG_ERROR("check sub stream failed: %s", nposStreamName.c_str());
                return -1;
            }
            const uint32_t ncount = baseMeta["ncount"].asUInt();
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            std::vector<uint8_t> raw(srcLen + 1);
            if (decodeWholeStream(metaStreams[id], temBuffer, dstLen, raw.data(), srcLen, "SEQ npos") != 0) {
                return -1;
            }
            /*
             * The positions are expanded here, once per block, for the same reason the lengths
             * are: the per-read walk asks for them one at a time, and a stream that only carries
             * what changed is smaller than one the walk must parse. Which form it is comes from
             * the stream's own name, and an archive written before the names carries the absolute
             * uint32 per N - that is what a "npos" without "delta" is.
             */
            const std::string& nposName = metaStreams[id]["sname"].asString();
            const bool nposDelta = (nposName == kSeqExcDeltaName) || metaStreams[id].isMember("delta");
            const bool nposRuns = (nposName == kSeqExcRunName);
            baseNPosVec.reserve(ncount);
            if (nposRuns) {
                /* One list of gaps, then one of lengths, both varints, both with the number of
                   runs as their count; each run covers the positions [end of the previous run +
                   gap, that + length). The lengths are read after the gaps, so the walk over the
                   gaps has to finish before they start - which is why the two lists are read in
                   one walk rather than by an offset. */
                const uint32_t runs = metaStreams[id]["count"].asUInt();
                std::vector<uint32_t> gaps;
                std::vector<uint32_t> lens;
                gaps.reserve(runs);
                lens.reserve(runs);
                uint32_t off = 0;
                for (uint32_t r = 0; r < runs; ++r) {
                    uint32_t v = 0;
                    if (!readVarint(raw.data(), srcLen, off, v)) {
                        LOG_ERROR("nposr stream truncated in block %llu", meta["block_id"].asInt64());
                        return -1;
                    }
                    gaps.push_back(v);
                }
                for (uint32_t r = 0; r < runs; ++r) {
                    uint32_t v = 0;
                    if (!readVarint(raw.data(), srcLen, off, v)) {
                        LOG_ERROR("nposr stream truncated in block %llu", meta["block_id"].asInt64());
                        return -1;
                    }
                    lens.push_back(v);
                }
                uint32_t next = 0;
                for (uint32_t r = 0; r < runs; ++r) {
                    next += gaps[r];
                    for (uint32_t k = 0; k < lens[r]; ++k) {
                        baseNPosVec.push_back(next);
                        next++;
                    }
                }
                if (baseNPosVec.size() != ncount) {
                    LOG_ERROR("nposr stream covers %zu positions, expected %u", baseNPosVec.size(), ncount);
                    return -1;
                }
            } else if (nposDelta) {
                uint32_t off = 0;
                uint32_t prev = 0;
                while (off < srcLen) {
                    uint32_t delta = 0;
                    if (!readVarint(raw.data(), srcLen, off, delta)) {
                        LOG_ERROR("npos stream truncated in block %llu", meta["block_id"].asInt64());
                        return -1;
                    }
                    prev += delta;
                    baseNPosVec.push_back(prev);
                }
                if (baseNPosVec.size() != ncount) {
                    LOG_ERROR("npos stream holds %zu positions, expected %u", baseNPosVec.size(), ncount);
                    return -1;
                }
            } else {
                if (srcLen != (ncount << 2)) {
                    LOG_ERROR("npos stream is %u bytes, expected %u", srcLen, ncount << 2);
                    return -1;
                }
                baseNPosVec.resize(ncount);
                for (uint32_t i = 0; i < ncount; ++i) {
                    baseNPosVec[i] = (uint32_t)raw[i * 4] | ((uint32_t)raw[i * 4 + 1] << 8) |
                                     ((uint32_t)raw[i * 4 + 2] << 16) | ((uint32_t)raw[i * 4 + 3] << 24);
                }
            }
            baseNPosBuffer = baseNPosVec.data();
            readOffset += metaStreams[id]["dstlen"].asUInt();
        } else {
            LOG_INFO("NCount is %d, not contain ncount stream", baseMeta["ncount"].asUInt());
        }

        /* check sub stream 5 */
        baseLengthVec.clear();
        if (baseMeta["minlen"].asUInt() != baseMeta["maxlen"].asUInt()) {
            id++;
            if (metaStreams[id]["sname"].asString() != "baselen") {
                LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
                return -1;
            }
            const uint32_t minLen = baseMeta["minlen"].asUInt();
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            std::vector<uint8_t> raw(srcLen + 1);
            if (decodeWholeStream(metaStreams[id], temBuffer, dstLen, raw.data(), srcLen, "SEQ length") != 0) {
                return -1;
            }
            /* Absolute lengths here, not "length above the shortest read": the stores that read
               this below ask for a length, and doing the addition once is doing it once. */
            if (metaStreams[id].isMember("delta")) {
                uint32_t off = 0;
                while (off < srcLen) {
                    uint32_t aboveMin = 0;
                    if (!readVarint(raw.data(), srcLen, off, aboveMin)) {
                        LOG_ERROR("baselen stream truncated in block %llu", meta["block_id"].asInt64());
                        return -1;
                    }
                    baseLengthVec.push_back(minLen + aboveMin);
                }
            } else {
                if (srcLen & 1) {
                    LOG_ERROR("baselen stream is %u bytes, not a whole number of uint16", srcLen);
                    return -1;
                }
                const uint32_t reads = srcLen >> 1;
                baseLengthVec.resize(reads);
                for (uint32_t i = 0; i < reads; ++i) {
                    baseLengthVec[i] = minLen + ((uint32_t)raw[i * 2] | ((uint32_t)raw[i * 2 + 1] << 8));
                }
            }
            readOffset += metaStreams[id]["dstlen"].asUInt();
        } else {
            LOG_INFO("Base length is same(%d, %d),  not contain baselen stream", baseMeta["minlen"].asUInt(), baseMeta["maxlen"].asUInt());
        }
    }

    // Initialize Comment decoder
    Json::Value& commentMeta = meta["comment"];
    if (commentMeta["type"].asString() == "plusonly") {
        commentType = CommentType::PLUS_ONLY;
    } else if (commentMeta["type"].asString() == "sameasid") {
        commentType = CommentType::SAME_AS_ID;
    } else if (commentMeta["type"].asString() == "other") {
        commentType = CommentType::OTHER;
        uint32_t commentDstLen = commentMeta["dstlen"].asUInt();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, commentDstLen, "comment");
        ioVecters.push_back(io);

        FieldDecoderArgs args;
        if (commentMeta["coder"].isMember("level")) {
            args.level = commentMeta["coder"]["level"].asInt();
        }
        commentDecoder = CoderFactory::makeFieldDecoder(commentMeta["coder"]["magic"].asString(), io.get(), args);
        if (commentDecoder == nullptr) {
            LOG_ERROR("Not support coder name = %s.", commentMeta["coder"]["magic"].asCString());
            return -1;
        }
        readOffset += commentDstLen;
    }

    // Initialize Quality decoder
    Json::Value& qualityMeta = meta["quality"];
    Json::Value& streamsMeta = qualityMeta["streams"];
    if (streamsMeta.size() != 2) {
        LOG_ERROR("Invalid quality streams. size = %d", streamsMeta.size());
        return -1;
    }

    if (streamsMeta[1]["coder"]["magic"].asString() == "coder_bwt_cm") {
        uint32_t qualityDstLen = streamsMeta[0]["dstlen"].asUInt();
        uint32_t freqDstLen = streamsMeta[1]["dstlen"].asUInt();
        coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualityDstLen, freqDstLen, &ioErrSink, "QUAL freq table");
        auto qualFreqCoder = std::make_unique<coder_bwt_cm>(&qualFreqIo);

        uint32_t qualFreqSrcLen = streamsMeta[1]["srclen"].asUInt();
        /*
         * The element count must be uint32_t: the array is allocated by this count,
         * while decode_line writes the untruncated qualFreqSrcLen bytes. A narrower
         * count would wrap once the quality-value alphabet exceeds 127 symbols
         * (qualFreqSrcLen > 510), shrinking the allocation while the write does not,
         * which overflows the heap.
         */
        uint32_t qualFreqArrLen = qualFreqSrcLen / sizeof(uint16_t);
        /* A vector rather than new[]/delete[]: the decoder below has a failure return of its own, and
           an early return between an allocation and its delete is a leak. */
        std::vector<uint16_t> qualFreqArray(qualFreqArrLen);
        uint32_t qualFreq = qualFreqCoder->decode_line((uint8_t*)qualFreqArray.data(), qualFreqSrcLen, UINT8_MAX, false);
        if (qualFreq != qualFreqSrcLen) {
            LOG_ERROR("Decode quality frequncy failed.");
            return -1;
        }

        for (uint32_t idx = 0; idx < qualFreqArrLen ; idx += 2) {
            qualityFreqTable.push_back(std::make_pair(qualFreqArray[idx], qualFreqArray[idx + 1]));
        }

        /*
         * The value stream's decoder is the one its magic names - the QUAL column is decided
         * per file, so this side cannot name a coder class (see CoderFactory::makeQualDecoder).
         * What it supplies is the alphabet just decoded above and, when the stream says the
         * encoder started from one, the prior the engine holds.
         */
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualityDstLen, "QUAL");
        ioVecters.push_back(io);

        QualDecoderArgs qualArgs;
        qualArgs.freqTable = qualityFreqTable.empty() ? nullptr : &qualityFreqTable;
        qualArgs.priorRequired = streamsMeta[0].isMember("prior");
        qualArgs.priorAddress = qualArgs.priorRequired ? streamsMeta[0]["prior"].asInt64() : -1;
        AuxPayloadPtr qualPriorBlob =
            (qualArgs.priorRequired && pbgzEngine != nullptr)
                ? pbgzEngine->getQualPrior(inBlockPtr->getPackageIndex())
                : AuxPayloadPtr();
        qualArgs.priorBlob = qualPriorBlob.get();

        qualityDecoder = CoderFactory::makeQualDecoder(streamsMeta[0]["coder"]["magic"].asString(),
                                                       io.get(), qualArgs);
        if (qualityDecoder == nullptr) {
            LOG_ERROR("Unsupport coder type: %s", streamsMeta[0]["coder"]["magic"].asString().c_str());
            return -1;
        }
    } else {
        LOG_ERROR("Unsupport coder type: %s", streamsMeta[1]["coder"]["magic"].asString().c_str());
        return -1;
    }

    return 0;
}

int32_t FastqCodecActuator::decodeWholeStream(const Json::Value& streamMeta, const uint8_t* src,
                                              uint32_t srcLen, uint8_t* dst, uint32_t dstLen,
                                              const char* tag)
{
    const std::string magic = streamMeta["coder"]["magic"].asString();
    std::shared_ptr<coder_io> io = makeCoderIo(src, srcLen, tag);
    io->meta = streamMeta;
    ioVecters.push_back(io);

    FieldDecoderArgs args;
    if (streamMeta["coder"].isMember("level")) {
        args.level = streamMeta["coder"]["level"].asInt();
    }
    std::shared_ptr<coder> decoder = CoderFactory::makeFieldDecoder(magic, io.get(), args);
    if (decoder == nullptr) {
        LOG_ERROR("Unsupported coder in %s stream: %s", tag, magic.c_str());
        return -1;
    }

    const int32_t got = decoder->decode_line(dst, dstLen, UINT8_MAX, false);
    if (got < 0 || (uint32_t)got != dstLen) {
        LOG_ERROR("%s stream decoded %d bytes, expected %u", tag, got, dstLen);
        return -1;
    }
    return 0;
}

int32_t FastqCodecActuator::initRleMatchStream(const Json::Value& metaStreams, uint32_t& id,
                                               uint32_t& readOffset)
{
    /*
     * Sub-stream "m" is the zero runs, sub-stream "mval" the values that ended them. Both carry
     * their own coder in their meta - the two halves are written independently - so each is
     * decoded with the one it names, and an archive written under a different choice still reads.
     */
    {
        const uint32_t runDstLen = metaStreams[id]["dstlen"].asUInt();
        const uint32_t runLen = metaStreams[id]["srclen"].asUInt();
        const uint32_t matchLen = metaStreams[id]["orgrawlen"].asUInt();

        std::vector<uint8_t> runBuf(runLen + 1);
        if (decodeWholeStream(metaStreams[id], inBlockPtr->getBuffer() + readOffset, runDstLen,
                              runBuf.data(), runLen, "SEQ match run") != 0) {
            return -1;
        }
        readOffset += runDstLen;

        std::vector<uint8_t> valBuf;
        uint32_t valLen = 0;
        if (id + 1 < metaStreams.size() && metaStreams[id + 1]["sname"].asString() == "mval") {
            id++;
            const uint32_t valDstLen = metaStreams[id]["dstlen"].asUInt();
            valLen = metaStreams[id]["srclen"].asUInt();
            valBuf.resize(valLen + 1);
            if (decodeWholeStream(metaStreams[id], inBlockPtr->getBuffer() + readOffset, valDstLen,
                                  valBuf.data(), valLen, "SEQ match val") != 0) {
                return -1;
            }
            readOffset += valDstLen;
        } else {
            /*
             * No surviving values: every byte of the match stream was zero. A block whose reads
             * all come from the reference is written this way, without the sub-stream at all, so
             * the streams that follow are addressed as if it had never existed.
             */
        }

        /*
         * Rebuild the stream the reads are served from: each run is followed by one surviving
         * value. safeAlloc is a calloc, so the runs themselves - the zeros - are already there,
         * and the trailing run after the last value needs no action either.
         */
        MemoryUtil::safeFree(baseMatchBlock);
        baseMatchBlock = MemoryUtil::safeAlloc<uint8_t>(matchLen + 1);
        if (baseMatchBlock == nullptr) {
            LOG_ERROR("Alloc SEQ match block buffer failed");
            return -1;
        }

        uint32_t rp = 0;
        uint32_t out = 0;
        for (uint32_t vp = 0; vp < valLen; ++vp) {
            uint32_t run = 0;
            uint32_t shift = 0;
            while (rp < runLen) {
                const uint8_t b = runBuf[rp++];
                run |= (uint32_t)(b & 0x7F) << shift;
                shift += 7;
                if ((b & 0x80) == 0) {
                    break;
                }
            }
            out += run;
            if (out < matchLen) {
                baseMatchBlock[out++] = valBuf[vp];
            }
        }
        baseMatchBlockLen = matchLen;
        baseMatchBlockOffset = 0;
        baseMatchRle = true;
    }
    return 0;
}

int32_t FastqCodecActuator::pullMatchBytes(uint8_t* dst, uint32_t len)
{
    if (!baseMatchRle) {
        return baseDecoder->decode_line(dst, len, UINT8_MAX, false);
    }
    if (baseMatchBlockOffset + len > baseMatchBlockLen) {
        LOG_ERROR("SEQ match stream ran out in block %llu: asked for %u at %u of %u",
                  meta["block_id"].asInt64(), len, baseMatchBlockOffset, baseMatchBlockLen);
        return -1;
    }
    std::memcpy(dst, baseMatchBlock + baseMatchBlockOffset, len);
    baseMatchBlockOffset += len;
    return (int32_t)len;
}

int32_t FastqCodecActuator::decompress() {
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());

    // Parse meta
    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);

    std::string idSplit = meta["id"]["splitsym"].asString();
    for (uint32_t i = 0; i < idSplit.length(); ++i) {
        idSplitSymbols.push_back(idSplit.c_str()[i]);
    }

    /*
     * Pre-allocation at the block entry: block_size from the file header (the upper
     * bound determined during compression) x 2 — a definite value, not an estimate.
     * Same scheme as on the SAM side (see SamCodecActuator::decompress): header output
     * <= block_size, and coder_fc's trailing SEQ staging <= block_size, exactly twice.
     * When block_size is 0 (not written by old files), fall back to the default block
     * size for the compression level.
     * ensureCapacity must run before initDecoder: inside initDecoder, coder_fc writes
     * the SEQ data at the end of the buffer, which would dangle after a realloc.
     */
    {
        size_t bs = pbgzEngine->getFileBlockSize();
        if (bs == 0) {
            bs = ConfigManager::getInstance().getBlockSizeByCompressLevel(pbgzEngine->getParameter().compressLevel);
        }
        if (outBlockPtr->ensureCapacity(bs * 2) != 0) {
            LOG_ERROR("preallocate output buffer failed, block_size=%zu", bs);
            return -1;
        }
    }

    if (0 != initDecoder(outBlockPtr)) {
        LOG_ERROR("Initital decoder failed.");
        return -1;
    }

    uint32_t lineNum = meta["idlines"].asUInt();
    uint8_t* pBaseOut = nullptr;
    uint32_t baseLines = 0;
    uint8_t* pBaseEnd = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize();
    if (meta["base"]["coder"]["magic"].asString() == "coder_fc") {
        pBaseOut = pBaseEnd - meta["base"]["totalsrclen"].asUInt();
    }
    uint32_t totalBaseLen = 0;
    uint64_t nposOffset = 0;
    bool isUseReference = (pReference != nullptr) && meta["base"].isMember("streams");
    for (uint32_t line = 0; line < lineNum; ++line) {
        uint8_t* idPtr = outBlockPtr->getCurrent();
        uint32_t idLen = 0;
        for (uint32_t splitIdx = 0; splitIdx < idDecoders.size(); ++splitIdx) {
            uint32_t idSplitLen = idDecoders[splitIdx]->decode_line(idPtr, outBlockPtr->getRemain(), idSplitSymbols[splitIdx]);
            idPtr += idSplitLen;
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idSplitLen);
            idLen += idSplitLen;
        }

        uint8_t* basePtr = outBlockPtr->getCurrent();
        uint32_t actualBaseLen = 0;
        if (!isUseReference) {
            if (minBaseLength == maxBaseLength) {
                actualBaseLen = maxBaseLength;
                if (meta["base"]["coder"]["magic"].asString() == "coder_fc") {
                    memcpy(basePtr, pBaseOut, actualBaseLen);
                    pBaseOut += actualBaseLen;
                } else {
                    baseDecoder->decode_line(basePtr, actualBaseLen, UINT8_MAX, false);
                }
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
                *outBlockPtr->getCurrent() = '\n';
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
            } else {
                if (meta["base"]["coder"]["magic"].asString() == "coder_fc") {
                    uint8_t* pBaseTmp = basePtr;
                    uint8_t* ptr = pBaseOut;
                    for (; ptr < pBaseEnd; ++ptr) {
                        *pBaseTmp++ = *ptr;
                        if (*ptr == '\n') {
                            break;
                        }
                    }
                    actualBaseLen = ptr - pBaseOut + 1;
                    pBaseOut = pBaseOut + actualBaseLen;
                } else {
                    /*
                     * Read lengths vary, so the line has to be delimited: the encoding side
                     * wrote each line with its trailing newline (see compressBaseWithoutRef),
                     * and decode_line returns the length it read including that delimiter. The
                     * newline is stripped below, exactly as the coder_fc branch above does.
                     *
                     * The capacity is the longest line plus its delimiter: decode_line fills
                     * the buffer up to out_len and only then looks for the separator, so a
                     * maxBaseLength-long line would hit the capacity check right before its
                     * newline and be reported as a short buffer.
                     */
                    int32_t decodedLen = baseDecoder->decode_line(basePtr, maxBaseLength + 1, '\n', false);
                    if (decodedLen <= 0) {
                        LOG_ERROR("base decode failed in block %llu, decoded %d",
                                  meta["block_id"].asInt64(), decodedLen);
                        return -1;
                    }
                    actualBaseLen = (uint32_t)decodedLen;
                }
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
                actualBaseLen -= 1;  // Remove newline character
            }
        } else {
            /* decode base mapping stream with strip N */;
            /* The length stream is only written when they vary (see compressBaseWithRef), so an
               empty column of them means every read is the longest one. */
            actualBaseLen = baseLengthVec.empty() ? maxBaseLength : baseLengthVec[baseLines];
            uint8_t* pout = outBlockPtr->getCurrent();

            const uint8_t actg4[4] = {'A', 'C', 'T', 'G'};
            if (baseNCount && nposOffset < baseNCount) { /* This block has N, and not all N's have been processed */
                /* Calculate the number of N's in the current base line */
                uint32_t n = totalBaseLen + actualBaseLen;
                uint64_t o = nposOffset;
                for (; o < baseNCount; o++) {
                    if (baseNPosBuffer[o] + 1 > n) {
                        break;
                    }
                }
                uint32_t ncntCurrLine = o - nposOffset;

                /* Decompress the mapping stream of the current base line */
                uint32_t stripNLength = actualBaseLen - ncntCurrLine;
                if (stripNLength > 0) {
                    const int32_t decoderLen = pullMatchBytes(baseStripNBuffer, stripNLength);
                    if (decoderLen < 0 || (uint32_t)decoderLen != stripNLength) {
                        LOG_ERROR("base decode failed in block %llu, expect len %u, actual %d", meta["block_id"].asInt64(), stripNLength, decoderLen);
                        return -1;
                    }
                }
                uint8_t* pdata;
                if (2 == baseMappedPairBuffer[baseLines]) {
                    if (baseLitBases) {
                        /* The read's own characters, written literally because it does not use the
                           reference (see compressBaseWithRef); they are already decrypted. */
                        pdata = baseStripNBuffer;
                    } else {
                        for (o = 0; o < stripNLength; o++) {
                            baseStripNBuffer[o] = actg4[baseStripNBuffer[o]];
                        }
                        pdata = baseStripNBuffer;
                    }
                } else {
                    /* Get the reference at the corresponding position */
                    pReference->getStretch2Bits1Char(refeStretchBuffer, stripNLength, baseMappedPosBuffer[baseLines]);

                    /* Restore the base squash stream after mapping */
                    actgXor(baseStripNBuffer, refeStretchBuffer, baseStripNBuffer, stripNLength);

                    pReference->getActgFrom2Bits(baseStripNBuffer, stripNLength, refeStretchBuffer);
                    if (baseMappedPairBuffer[baseLines]) {
                        actgPair(baseStripNBuffer, refeStretchBuffer, stripNLength);
                        pdata = baseStripNBuffer;
                    } else {
                        pdata = refeStretchBuffer;
                    }
                }

                /* Copy the restored data to the output buffer */
                for (o = 0, n = 0; n < actualBaseLen; n++) {
                    if (nposOffset < baseNCount && (totalBaseLen + n) == baseNPosBuffer[nposOffset]) { /* Current position is N */
                        *pout++ = 'N' ;
                        nposOffset++;
                    } else {
                        /* Current position is not N */
                        *pout++ = (pdata[o++]);
                    }
                }
            } else { /* No N in the block */
                /* Decompress the mapping stream of the current base line */
                uint32_t stripNLength = actualBaseLen;
                const int32_t decoderLen = pullMatchBytes(baseStripNBuffer, stripNLength);
                if (decoderLen < 0 || (uint32_t)decoderLen != stripNLength) {
                    LOG_ERROR("base decode failed in block %llu, expect len %u, actual %d", meta["block_id"].asInt64(), stripNLength, decoderLen);
                    return -1;
                }
                if (2 == baseMappedPairBuffer[baseLines]) {
                    if (baseLitBases) {
                        /* The read's own characters, written literally (see compressBaseWithRef). */
                        memcpy(pout, baseStripNBuffer, stripNLength);
                    } else {
                        for (uint32_t o = 0; o < stripNLength; o++) {
                            pout[o] = actg4[baseStripNBuffer[o]];
                        }
                    }
                } else {
                    /* Get the reference at the corresponding position */
                    pReference->getStretch2Bits1Char(refeStretchBuffer, stripNLength, baseMappedPosBuffer[baseLines]);
                    /* Restore the base squash stream after mapping */
                    actgXor(baseStripNBuffer, refeStretchBuffer, baseStripNBuffer, stripNLength);
                    if (baseMappedPairBuffer[baseLines]) {
                        pReference->getActgFrom2Bits(baseStripNBuffer, stripNLength, refeStretchBuffer);
                        actgPair(pout, refeStretchBuffer, stripNLength);
                    } else {
                        pReference->getActgFrom2Bits(baseStripNBuffer, stripNLength, pout);
                    }
                }
            }
            totalBaseLen += actualBaseLen;
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
            *(outBlockPtr->getCurrent()) = '\n';
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
        }
        ++baseLines;

        // Decode comment
        uint8_t* commentPtr = outBlockPtr->getCurrent();
        switch (commentType) {
            case CommentType::PLUS_ONLY: {
                *commentPtr++ = '+';
                *commentPtr = '\n';
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 2);
                break;
            }
            case CommentType::SAME_AS_ID: {
                memcpy(commentPtr, idPtr, idLen);
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idLen);
                break;
            }
            case CommentType::OTHER: {
                int32_t commentLen = commentDecoder->decode_line(commentPtr, outBlockPtr->getRemain(), '\n', false);
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + commentLen);
                break;
            }
            default:
                break;
        }

        // Decode quality values, one record at a time: the record's base line is the context,
        // whichever coder the stream's magic named.
        qualityDecoder->decode_record(outBlockPtr->getCurrent(), actualBaseLen, basePtr, actualBaseLen);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
        *outBlockPtr->getCurrent() = '\n';
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
    }

    // Check file MD5
    std::string md5;
    calcMd5sum(md5, outBlockPtr->getBuffer(), outBlockPtr->getDataLen());
    if (md5 != meta["md5"].asString()) {
        LOG_ERROR("Md5 check failed for block(%d), expect %s, got %s.", inBlockPtr->getBlockId(),
                meta["md5"].asCString(), md5.c_str());
        return -1;
    }

    return 0;
}
