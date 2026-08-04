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
#include "utils/md5_util.h"
#include "coder_json.h"
#include "io_wrapper.h"
#include "actg.h"
#include "city.h"
#include "pbgz_manager.h"
#include "pbgz_stat.h"

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

// SIMD-optimized N character counting using SSE4.2
uint32_t FastqCodecActuator::countN_SSE2(const uint8_t* data, size_t length) {
    uint32_t count = 0;
    size_t i = 0;

#ifdef __SSE4_2__
    const __m128i target_upper = _mm_set1_epi8('N');
    const __m128i target_lower = _mm_set1_epi8('n');

    for (; i + 16 <= length; i += 16) {
        __m128i chunk = _mm_loadu_si128((__m128i*)(data + i));
        __m128i cmp_upper = _mm_cmpeq_epi8(chunk, target_upper);
        __m128i cmp_lower = _mm_cmpeq_epi8(chunk, target_lower);
        __m128i cmp = _mm_or_si128(cmp_upper, cmp_lower);
        uint32_t mask = _mm_movemask_epi8(cmp);
        count += __builtin_popcount(mask);
    }
#endif

    for (; i < length; i++) {
        count += (data[i] == 'N' || data[i] == 'n');
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

        count += (b0 == 'N' || b0 == 'n');
        count += (b1 == 'N' || b1 == 'n');
        count += (b2 == 'N' || b2 == 'n');
        count += (b3 == 'N' || b3 == 'n');
        count += (b4 == 'N' || b4 == 'n');
        count += (b5 == 'N' || b5 == 'n');
        count += (b6 == 'N' || b6 == 'n');
        count += (b7 == 'N' || b7 == 'n');
    }

    for (; i < length; i++) {
        count += (data[i] == 'N' || data[i] == 'n');
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
        idPositions.push_back(static_cast<uint8_t>(pos));
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
    uint32_t baselenLen = ((minBaseLength == maxBaseLength) ? 0 : (isGen2 ? (line4 << 1) : (line4 << 2)));

    if (pReference) {
       /* base pair + 4 base pair squash + 4 base squash + base mapped + base N pos in block +
        * base delete N + mapped pos + mapped pair + baselen each line
        */ 
        uint32_t n = lmax + (lsquash << 3) + baseMappedLength + (baseNCount << 2);
        n += lmax + (line4 << 3) + line4 + baselenLen;
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

        if (baselenLen) {
            if (isGen2) {
                baseLengthGen2Buffer = (uint16_t *)p;
            } else {
                baseLengthGen3Buffer = (uint32_t *)p;
            }
            p += baselenLen;
        }
        mapping = (isGen2) ? (&FastqCodecActuator::mappingFastqGen2) : (&FastqCodecActuator::mappingFastQGen3);
    }

    return 0;
}

void FastqCodecActuator::mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir) {
    Mapping mappingTable[4];
    uint8_t ch;
    const uint32_t baseGroupLen = pReference->getBaseGroupLength();
    const uint32_t bgMid = baseGroupLen >> 1;
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
    actgPair(basePairBuffer, base, baseLength);

    /* Calculate align4 squash buffer and pair squash buffer */
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
            if (bestUnmatches <= MAPPED_THRESHOLD_GEN2) {
                break;
            }
        }
        if (bestUnmatches <= MAPPED_THRESHOLD_GEN2) {
            break;
        }
        mappingTable[align4Curr].incOffset();
    }

    if (bestUnmatches > MAPPED_THRESHOLD_GEN2) { /* continue mapping */
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
                if (bestUnmatches <= MAPPED_THRESHOLD_GEN2) {
                    break;
                }
            }

            if (bestUnmatches <= MAPPED_THRESHOLD_GEN2) {
                break;
            }
            mappingTable[align4Curr].incOffset();
        }
    }

    /* calc the result of base or base pair mapping with the match pos reference */
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

    mappingPos = bestPos;
    mappingDir = bestIsPair;
    return;
}

void FastqCodecActuator::mappingFastQGen3(const uint8_t*, uint32_t, uint8_t*&, uint32_t&, uint64_t&, uint8_t&) {
    LOG_ERROR("Not support FASTQ Gen3");
    return;
}

template <typename TCoder>
int32_t FastqCodecActuator::compressIdStream(coder_io* idIo, TCoder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx) {
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
    std::shared_ptr<coder_io> idIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> idCoder = std::make_shared<coder_bwt_cm>(idIo.get());

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
    
    return outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
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
                LOG_DEBUG("Is all digit(%d), use coder_bwt_cm.", i);
                std::shared_ptr<coder_io> idIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
                std::shared_ptr<coder_bwt_cm> idCoder = std::make_shared<coder_bwt_cm>(idIo.get());
                uint32_t srcLength = 0;
                int32_t ret= compressIdStream<coder_bwt_cm>(idIo.get(), idCoder.get(), streamMeta, srcLength, i);
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
        LOG_DEBUG("Not all digit or fix length(%d), use coder_affix_match.", i);
        std::shared_ptr<coder_io> idIoAM = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        std::shared_ptr<coder_affix_match> idCoderAm = std::make_shared<coder_affix_match>(idIoAM.get());
        uint32_t srcLength = 0;
        int32_t ret= compressIdStream<coder_affix_match>(idIoAM.get(), idCoderAm.get(), streamMeta, srcLength, i);
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
    
    return outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
}

int32_t FastqCodecActuator::compressId() {
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
        static bool isPrint = false;
        if (!isGen2 && pReference != nullptr && !isPrint) {
            fprintf(stderr ,"This file is Gen3 FASTQ, compression will be performed without a reference genome.\n");
            isPrint = true;
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
    const __m128i target_lower = _mm_set1_epi8('n');
#endif
    Json::Value metaSubs;
    Json::Value metaStreams;
    Json::Value metaBase;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;
    std::shared_ptr<coder_io> matchIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> matchCm = std::make_shared<coder_bwt_cm>(matchIo.get());
    int64_t srcLen = 0;
    for (uint32_t i = 1; i < line; i += 4) {
        uint32_t endPos = inBlockPtr->getNpos()[i];
        uint32_t startPos = inBlockPtr->getNpos()[i - 1] + 1;
        uint8_t* pBuff = baseStripNBuffer;
        uint32_t outLen = 0;

        // Use SIMD-optimized N detection when processing base sequences
        size_t segmentLength = endPos - startPos;
        const uint8_t* segmentStart = ptr + startPos;
        size_t simdIdx = 0;

#ifdef __SSE4_2__
        // Process 16 bytes at a time with SIMD
        for (; simdIdx + 16 <= segmentLength; simdIdx += 16) {
            __m128i chunk = _mm_loadu_si128((__m128i*)(segmentStart + simdIdx));
            __m128i cmp_upper = _mm_cmpeq_epi8(chunk, target_upper);
            __m128i cmp_lower = _mm_cmpeq_epi8(chunk, target_lower);
            __m128i cmp = _mm_or_si128(cmp_upper, cmp_lower);
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
            if (ch == 'N' || ch == 'n') {
                *(baseNPosBuffer + nOffset) = currPos + simdIdx;
                nOffset++;
            } else {
                *pBuff = ch;
                pBuff++;
            }
        }

        (this->*mapping)(baseStripNBuffer, pBuff - baseStripNBuffer,
                         baseMappedBuffer, outLen, baseMappedPosBuffer[offset], baseMappedPairBuffer[offset]);
        srcLen += pBuff - baseStripNBuffer;
        matchCm->encode_line(baseMappedBuffer, outLen);
        pReference->updateMatchedGene(baseMappedPosBuffer[offset], outLen);
        if (encBaseLen) {
            baseLengthGen2Buffer[offset] = endPos - startPos - minBaseLength;
        }
        offset++;
        currPos += endPos - startPos;
    }
    /* First sub-stream: match stream between reads and reference */
    matchCm->encode_flush();
    metaSubs.clear();
    metaSubs["srclen"] =  (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] = (Json::Value::Int)(matchIo->data_len);
    metaSubs["coder"] = matchIo->meta;
    metaSubs["sname"] = "m";
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + matchIo->data_len);
    totalDstLen += matchIo->data_len;
    totalSrcLen += srcLen;

    /* Second sub-stream: position stream of matches between reads and reference */
    std::shared_ptr<coder_io> posIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    uint32_t line4 = inBlockPtr->getNpos().size() >> 2;
    srcLen = (line4 << 3);
    if (srcLen > FC_MIN_LEN && srcLen < FC_MAX_LEN) {
        std::shared_ptr<coder_fc> posCm = std::make_shared<coder_fc>(posIo.get());
        posCm->encode_line((uint8_t *)baseMappedPosBuffer, srcLen);
        posCm->encode_flush();
    } else {
        std::shared_ptr<coder_bwt_cm> posCm = std::make_shared<coder_bwt_cm>(posIo.get());
        posCm->encode_line((uint8_t *)baseMappedPosBuffer, srcLen);
        posCm->encode_flush();
    }
    metaSubs.clear();
    metaSubs["srclen"] =  (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] =  (Json::Value::Int)posIo->data_len;
    metaSubs["coder"] = posIo->meta;
    metaSubs["sname"] = "mpos";
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + posIo->data_len);
    totalSrcLen += srcLen;
    totalDstLen += posIo->data_len;

    /* Third sub-stream: pair identifier stream of matches between reads and reference */
    auto pairIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    srcLen = line4;
    if (srcLen > FC_MIN_LEN && srcLen < FC_MAX_LEN) {
        std::shared_ptr<coder_fc> subCoder = std::make_shared<coder_fc>(pairIo.get());
        subCoder->encode_line((uint8_t *)baseMappedPairBuffer, srcLen);
        subCoder->encode_flush();
    } else {
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(pairIo.get());
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

    /* Fourth sub-stream: positions of all N's in reads */
    if (baseNCount > 0) {
        std::shared_ptr<coder_io> nposIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        srcLen = (baseNCount << 2);
        if (false) {
            std::shared_ptr<coder_fc> subCoder = std::make_shared<coder_fc>(nposIo.get());
            subCoder->encode_line((uint8_t *)baseNPosBuffer, srcLen);
            subCoder->encode_flush();
        } else {
            std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(nposIo.get());
            subCoder->encode_line((uint8_t *)baseNPosBuffer, srcLen);
            subCoder->encode_flush();
        }
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)srcLen;
        metaSubs["dstlen"] = (Json::Value::Int)nposIo->data_len;
        metaSubs["coder"] = nposIo->meta;
        metaSubs["sname"] = "npos";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + nposIo->data_len);
        totalSrcLen += srcLen;
        totalDstLen += nposIo->data_len;
    } else {
        LOG_INFO("NCount is 0, base info not contains npos stream");
    }
    metaBase["ncount"] = (Json::Value::UInt)baseNCount;

    /* Fifth stream: length of each base line */
    if (encBaseLen) {
        std::shared_ptr<coder_io> lenIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        srcLen = (line4 << 1);
        std::shared_ptr<coder_bwt_cm> sub_coder = std::make_shared<coder_bwt_cm>(lenIo.get());
        sub_coder->encode_line((uint8_t *)baseLengthGen2Buffer, srcLen);
        sub_coder->encode_flush();
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)srcLen;
        metaSubs["dstlen"] = (Json::Value::Int)lenIo->data_len;
        metaSubs["coder"] = lenIo->meta;
        metaSubs["sname"] = "baselen";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + lenIo->data_len);
        totalSrcLen += srcLen;
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

    std::shared_ptr<coder_io> baseIo;
    if (baseSrcLength <= FC_MIN_LEN || baseSrcLength >= FC_MAX_LEN) {
        baseIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        std::shared_ptr<coder_bwt_cm> baseCoder = std::make_shared<coder_bwt_cm>(baseIo.get());

        for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx += 4) {
            int64_t end = inBlockPtr->getNpos()[idx];
            int64_t start = inBlockPtr->getNpos()[idx - 1] + 1;
            int32_t lineLength = end - start + addLength;
            baseCoder->encode_line(inBlockPtr->getBuffer() + start, lineLength);
        }
        baseCoder->encode_flush();
    } else {
        uint8_t* tmpBase = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize() - baseSrcLength;
        uint32_t srcLength = 0;
        for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx += 4) {
            int64_t end = inBlockPtr->getNpos()[idx];
            int64_t start = inBlockPtr->getNpos()[idx-1] + 1;
            int32_t lineLength = end - start + addLength;
            memcpy(tmpBase + srcLength, inBlockPtr->getBuffer() + start, lineLength);
            srcLength += lineLength;
        }

        baseIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain() - baseSrcLength);
        std::shared_ptr<coder_fc> baseCoder = std::make_shared<coder_fc>(baseIo.get());
        baseCoder->encode_line(tmpBase, srcLength);
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
    std::shared_ptr<coder_io> commentIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> commentCoder = std::make_shared<coder_bwt_cm>(commentIo.get());
    
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

int32_t FastqCodecActuator::compressQuality() {
    std::shared_ptr<coder_io> qualityIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_qual> qualityCoder = std::make_shared<coder_qual>(qualityIo.get(), true, qualityFreqTable);

    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;
    Json::Value streamMeta;

    // Encode quality data
    uint32_t streamSrcLen = 0;
    for (uint32_t idx = 3; idx < inBlockPtr->getNpos().size(); idx += 4) {
        uint32_t end = inBlockPtr->getNpos()[idx];
        uint32_t start = inBlockPtr->getNpos()[idx - 1] + 1;
        uint32_t lineLength = end - start;
        qualityCoder->encode_qual_gen2(inBlockPtr->getBuffer() + inBlockPtr->getNpos()[idx - 3] + 1,
                                       inBlockPtr->getBuffer() + start, lineLength);
        streamSrcLen += lineLength;
    }
    qualityCoder->encode_flush();
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityIo->data_len);
    
    Json::Value subMeta;
    subMeta["srclen"] = streamSrcLen;
    subMeta["dstlen"] = qualityIo->data_len;
    subMeta["coder"] = qualityIo->meta;
    subMeta["streamname"] = "qual";
    streamMeta.append(subMeta);

    totalSrcLength += streamSrcLen;
    totalDstLength += qualityIo->data_len;

    // Encode quality frequency table
    std::shared_ptr<coder_io> qualityFreqIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> qualityFreqCoder = std::make_shared<coder_bwt_cm>(qualityFreqIo.get());
    std::shared_ptr<uint16_t[]> qualiltyFreqArray = std::make_unique<uint16_t[]>(qualityFreqTable.size()<< 1);
    for (uint32_t i = 0; i < qualityFreqTable.size(); ++i) {
        int idx = i << 1;
        qualiltyFreqArray[idx] = qualityFreqTable[i].first;
        qualiltyFreqArray[idx + 1] = qualityFreqTable[i].second;
    }

    uint32_t freqSrcLen = (qualityFreqTable.size() << 1) * sizeof(uint16_t);
    qualityFreqCoder->encode_line((uint8_t*)qualiltyFreqArray.get(), freqSrcLen);

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
    Json::Value& idStreamMeta = meta["id"]["streams"];
    if (idStreamMeta.size() != idSplitSymbols.size()) {
        LOG_ERROR("id stream meta size not match id split symbols size(%d, %d)", idStreamMeta.size(), idSplitSymbols.size());
        return -1;
    }

    uint32_t readOffset = 0;
    // Initialize ID decoder
    for (uint32_t idx = 0; idx < idStreamMeta.size(); ++idx) {
        std::string coderName = idStreamMeta[idx]["coder"]["magic"].asString();
        uint32_t dstLen = idStreamMeta[idx]["dstlen"].asUInt();
        if (coderName == "coder_affix_match") {
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLen);
            ioVecters.push_back(io);
            idDecoders.push_back(std::make_shared<coder_affix_match>(io.get()));
            idDecoders.back()->set_level(idStreamMeta[idx]["coder"]["level"].asInt());
        } else if (coderName == "coder_bwt_cm") {
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLen);
            ioVecters.push_back(io);
            idDecoders.push_back(std::make_shared<coder_bwt_cm>(io.get()));
            idDecoders.back()->set_level(idStreamMeta[idx]["coder"]["level"].asInt());
        } else {
            LOG_ERROR("unsupported coder name: %s", coderName.c_str());
            return -1;
        }
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
            coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstBaseLen);
            baseIo.meta = baseMeta;
            baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
            coder_fc baseDecoderFc(&baseIo);
            baseDecoderFc.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcBaseLen, srcBaseLen,  UINT8_MAX, false);
        } else if (baseCoderName == "coder_bwt_cm") {
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstBaseLen);
            ioVecters.push_back(io);
            baseDecoder = std::make_shared<coder_bwt_cm>(io.get());
        } else {
            LOG_ERROR("unsupported coder name: %s", baseCoderName.c_str());
            return -1;
        }
        readOffset += dstBaseLen;
    } else {
        Json::Value metaStreams = baseMeta["streams"];
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
        if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {   
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, baseDstLen);
            ioVecters.push_back(io);
            baseDecoder = std::make_shared<coder_bwt_cm>(io.get());
        } else {
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }
        readOffset += baseDstLen;

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
        if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asInt();
            baseMappedPosBuffer = (uint64_t *)ps;
            ps += srcLen;
            coder_io posIo(temBuffer, dstLen);
            auto posCm = std::make_unique<coder_bwt_cm>(&posIo);
            posCm->decode_line((uint8_t *)baseMappedPosBuffer, srcLen, UINT8_MAX, false);
        } else if (metaStreams[id]["coder"]["magic"].asString() == "coder_fc") {
            temBuffer =  inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPosBuffer = (uint64_t *)ps;
            ps += srcLen;
            coder_io posIo(temBuffer, dstLen);
            posIo.meta = metaStreams[id];
            // posIo.meta["dstlen"] = metaStreams[id]["dstlen"]; // Compatible field, can be unified
            coder_fc posFc(&posIo);
            posFc.decode_line((uint8_t*)baseMappedPosBuffer, srcLen, UINT8_MAX, false);
        } else {
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }
            
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
            coder_io pairIo(temBuffer, dstLen);
            auto pairCm = std::make_unique<coder_bwt_cm>(&pairIo);
            pairCm->decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else if (metaStreams[id]["coder"]["magic"].asString() == "coder_fc") {
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPairBuffer = ps;
            ps += srcLen;
            coder_io pairIo(temBuffer, dstLen);
            pairIo.meta = metaStreams[id];
            // pairIo.meta["tot_dstlen"] = metaStreams[id]["dstlen"]; // Compatible field, can be unified
            coder_fc pairFc(&pairIo);
            pairFc.decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else{
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }        
        readOffset += metaStreams[id]["dstlen"].asUInt();

        /* check sub stream 4 */
        baseNPosBuffer = nullptr;
        if (baseMeta["ncount"].asUInt()) {
            id++;
            if (metaStreams[id]["sname"].asString() != "npos") {
                LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
                return -1;
            }
            if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
                temBuffer = inBlockPtr->getBuffer() + readOffset;
                srcLen = metaStreams[id]["srclen"].asUInt();
                dstLen = metaStreams[id]["dstlen"].asUInt();
                baseNPosBuffer = (uint32_t *)ps;
                ps += srcLen;
                coder_io nposIo(temBuffer, dstLen);
                auto nposCm = std::make_unique<coder_bwt_cm>(&nposIo);
                nposCm->decode_line((uint8_t *)baseNPosBuffer, srcLen, UINT8_MAX, false);
            } else {
                LOG_ERROR("check sub stream failed: coder name unmatch");
                return -1;
            }
            readOffset += metaStreams[id]["dstlen"].asUInt();
        } else {
            LOG_INFO("NCount is %d, not contain ncount stream", baseMeta["ncount"].asUInt());
        }

        /* check sub stream 5 */
        baseLengthGen2Buffer = nullptr;
        if (baseMeta["minlen"].asUInt() != baseMeta["maxlen"].asUInt()) {
            id++;
            if (metaStreams[id]["sname"].asString() != "baselen") {
                LOG_ERROR("check sub stream failed: %s", metaStreams[id]["sname"].asString().c_str());
                return -1;
            }        
            if (metaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm")  {
                temBuffer = inBlockPtr->getBuffer() + readOffset;
                srcLen = metaStreams[id]["srclen"].asUInt();
                dstLen = metaStreams[id]["dstlen"].asUInt();
                baseLengthGen2Buffer = (uint16_t *)ps;
                ps += srcLen;
                coder_io lenIo(temBuffer, dstLen);
                auto lenCm = std::make_unique<coder_bwt_cm>(&lenIo);
                lenCm->decode_line((uint8_t *)baseLengthGen2Buffer, srcLen, UINT8_MAX, false);
            } else  {
                LOG_ERROR("check sub stream failed: coder name unmatch");
                return -1;
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
        if (commentMeta["coder"]["magic"].asString() == "coder_bwt_cm") {
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, commentDstLen);
            ioVecters.push_back(io);
            commentDecoder = std::make_shared<coder_bwt_cm>(io.get());
        } else {
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
        coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualityDstLen, freqDstLen);
        auto qualFreqCoder = std::make_unique<coder_bwt_cm>(&qualFreqIo);

        uint32_t qualFreqSrcLen = streamsMeta[1]["srclen"].asUInt();
        /*
         * 元素个数必须用 uint32_t: 数组按这个数分配, 而 decode_line 按未截断的
         * qualFreqSrcLen 字节写入。原来是 uint8_t, 一旦质量值字母表超过 127 个符号
         * (qualFreqSrcLen > 510) 计数就回绕, 分配变小而写入不变, 直接堆越界。
         */
        uint32_t qualFreqArrLen = qualFreqSrcLen / sizeof(uint16_t);
        uint16_t* qualFreqArray = new uint16_t[qualFreqArrLen];
        uint32_t qualFreq = qualFreqCoder->decode_line((uint8_t*)qualFreqArray, qualFreqSrcLen, UINT8_MAX, false);
        if (qualFreq != qualFreqSrcLen) {
            LOG_ERROR("Decode quality frequncy failed.");
            return -1;
        }

        for (uint32_t idx = 0; idx < qualFreqArrLen ; idx += 2) {
            qualityFreqTable.push_back(std::make_pair(qualFreqArray[idx], qualFreqArray[idx + 1]));
        }

        delete [] qualFreqArray;

        if (streamsMeta[0]["coder"]["magic"].asString() == "coder_qual") {
            std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, qualityDstLen);
            ioVecters.push_back(io);
            qualityDecoder = std::make_shared<coder_qual>(io.get(), true, qualityFreqTable);
        } else {
            LOG_ERROR("Unsupport coder type: %s", streamsMeta[0]["coder"]["magic"].asString().c_str());
            return -1;
        }
    } else {
        LOG_ERROR("Unsupport coder type: %s", streamsMeta[1]["coder"]["magic"].asString().c_str());
        return -1;
    }   

    return 0;
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
                } else if (meta["base"]["coder"]["magic"].asString() == "coder_bwt_cm") {
                    baseDecoder->decode_line(basePtr, actualBaseLen, UINT8_MAX, false);
                } else {
                    LOG_ERROR("Unsupport coder type: %s", meta["base"]["coder"]["magic"].asString().c_str());
                    return -1;
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
                } else if (meta["base"]["coder"]["magic"].asString() == "coder_bwt_cm") {
                    baseDecoder->decode_line(basePtr, maxBaseLength, '\n', false);
                } else {
                    LOG_ERROR("Unsupport coder type: %s", meta["base"]["coder"]["magic"].asString().c_str());
                    return -1;
                }
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
                actualBaseLen -= 1;  // Remove newline character
            }
        } else {
            /* decode base mapping stream with strip N */;
            actualBaseLen = (baseLengthGen2Buffer) ? (baseLengthGen2Buffer[baseLines] + minBaseLength) : maxBaseLength;
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
                    uint32_t decoderLen = baseDecoder->decode_line(baseStripNBuffer, stripNLength, UINT8_MAX, false);
                    if (stripNLength != decoderLen) {
                        LOG_ERROR("base decode failed in block %llu, expect len %u, actual %u", meta["block_id"].asInt64(), stripNLength, decoderLen);
                        return -1;
                    }
                }
                uint8_t* pdata;
                if (2 == baseMappedPairBuffer[baseLines]) {
                    for (o = 0; o < stripNLength; o++) {
                        baseStripNBuffer[o] = actg4[baseStripNBuffer[o]];
                    }
                    pdata = baseStripNBuffer;
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
                uint32_t decoderLen = baseDecoder->decode_line(baseStripNBuffer, stripNLength, UINT8_MAX, false);
                if (stripNLength != decoderLen) {
                    LOG_ERROR("base decode failed in block %llu, expect len %u, actual %u", meta["block_id"].asInt64(), stripNLength, decoderLen);
                    return -1;
                }
                if (2 == baseMappedPairBuffer[baseLines]) {
                    for (uint32_t o = 0; o < stripNLength; o++) {
                        pout[o] = actg4[baseStripNBuffer[o]];
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

        // Decode quality values
        qualityDecoder->decode_qual_gen2(basePtr, outBlockPtr->getCurrent(), actualBaseLen);
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
