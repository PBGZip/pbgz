#include <cstring>
#include <algorithm>
#include <memory>

#include "fastq_actuator.h"
#include "reference.h"
#include "log/logger.h"
#include "utils/memory_util.h"
#include "coder_io.h"
#include "coder_bwt_cm.h"
#include "coder_affix_match.h"
#include "coder_fc.h"
#include "coder_qual.h"
#include "utils/md5util.h"
#include "coder_json.h"
#include "io_wrapper.h"
#include "actg.h"
#include "city.h"
#include "pbgz_manager.h"

const uint32_t MAPPED_THRESHOLD_GEN2 = 2;

int32_t FastqActuator::compress() {
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
    meta["idlines"] = inBlockPtr->getNpos().size() >> 2;

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
int32_t FastqActuator::preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLen) {
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

int32_t FastqActuator::preAnalysisId(uint8_t* pBuffer, uint32_t bufferLen) {
    uint32_t lastPos = 0;
    uint32_t lastFindPos = 0;
    for (uint32_t idx = 0; idx < idSplitSymbols.size(); ++idx) {
        uint8_t symbol = idSplitSymbols[idx];
        uint32_t pos = lastPos;
        for (; pos < bufferLen; ++pos) {
            if (*(pBuffer + pos) == symbol) {
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
                break;
            }
            lastPos = pos + 1;
        }

        if (pos > bufferLen) {
            idPosLength = UINT32_MAX; // Mark as unavailable
            break;
        }
    }
    return 0;
}

int32_t FastqActuator::preAnalysisBase(uint8_t* pBuffer, uint32_t bufLen) {
    uint32_t baseLength = bufLen - 1;    // Remove newline character
    if (baseLength > maxBaseLength) {
        maxBaseLength = baseLength;
    }
    if (baseLength < minBaseLength) {
        minBaseLength = baseLength;
    }

    for (uint32_t idx = 0; idx < bufLen - 1; ++idx) {
        uint8_t ch = *(pBuffer + idx); 
        switch (ch)
        {
        case 'a':
        case 'c':
        case 't':
        case 'g':
        case 'A':
        case 'C':
        case 'T':
        case 'G':   
            break;
        case 'N':
        case 'n':
            baseNCount++;
            break;
        
        default:
            break;
        }  
    }

    return 0;
}
   
int32_t FastqActuator::preAnalysisComment(uint8_t* pBuffer, uint32_t bufLen, uint32_t lineNo) {
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

int32_t FastqActuator::preAnalysis() {
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
    return 0;
}

int32_t FastqActuator::initEncoder() {
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
        mapping = (isGen2) ? (&FastqActuator::mappingFastqGen2) : (&FastqActuator::mappingFastQGen3);
    }

    return 0;
}

void FastqActuator::mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir) {
    Mapping mt[4];
    uint32_t n, m, o, l, squashLen[2], loffset;
    uint32_t align4_curr;
    uint8_t *psquash, *psquash_refe, ch;
    const uint32_t baseGroupLen = pReference->getBaseGroupLength();
    const uint32_t bg_mid = baseGroupLen >> 1;
    const uint8_t *pseq[2] = {base, basePairBuffer};
    uint32_t hash32, pos_cnts, *pos_vals, total;
    const uint8_t bg_is_unalign4 = !!(baseGroupLen & 0x3);
    const uint32_t len_bgs = (baseGroupLen >> 2) + bg_is_unalign4;

    uint32_t base_squash_align4 = (baseLength >> 2) + !!(baseLength & 0x3);
    uint32_t match_pair, match_pair_origin;
    uint32_t match_pos, unmatchs = baseLength;
    uint8_t *prefe_squash = (uint8_t*)(pReference->getSquash());
    int64_t refe_squashlen = pReference->getSquashLength();
    uint32_t best_pos_inrefe = UINT32_MAX, best_is_pair = 0; /* Record the position and strand direction in reference for current best match */
    uint32_t best_unmatchs = UINT32_MAX, best_align4;
    uint32_t best_pos = UINT32_MAX; /* Position converted to reference original base */
    uint64_t xsquash, xsquash_match, xsquash_macth_refe;
    const uint64_t xsquash_tab[2] = {0xFCFFFFFFFFFFFFFF, 0xFCFFFFFFFFFFFFFF};

    /* case 1: base length is not greater than reference index corresponding base length */
    if (baseLength <= (baseGroupLen + 4)) {
        actgEncode(base, out, baseLength);
        outLength = baseLength;
        mappingPos = 0;
        mappingDir = 2; /* During decompression, check mdir first; if it's 2, it means no match */
        return;
    }

    /* case 2: base length is greater than reference index corresponding base length */
    uint32_t e = baseLength - baseGroupLen;
    actgPair(basePairBuffer, base, baseLength);

    /* Calculate align4 squash buffer and pair squash buffer */
    for (n = 0; n < 4; n++) {
        squashLen[0] = (baseLength - n) >> 2;
        total = squashLen[0] << 2;
        mt[n].leftUnalignLen[0] = n;
        for (m = 0; m < n; m++) {
            mt[n].leftUnalign[0][m] = ((*(pseq[0] + m)) >> 1) & 0x3; /* squash value of bases not 4-aligned on the left */
        }
        mt[n].rightUnalignLen[0] = baseLength - n - total;
        for (m = 0; m < mt[n].rightUnalignLen[0]; m++) {
            mt[n].rightUnalign[0][m] = ((*(pseq[0] + baseLength - mt[n].rightUnalignLen[0] + m)) >> 1) & 0x3;
        }
        actgSquash(pseq[0] + n, total, baseSquashBuffer[n]);

        squashLen[1] = (baseLength - n + 1) >> 2; /* Add a character to the right for 32-byte alignment; need to handle the last character when matching to mt[0]'s pair and mt[0] offset is 0 */
        total = squashLen[1] << 2;
        mt[n].leftUnalignLen[1] = baseLength + 1 - n - total;
        for (m = 0; m < mt[n].leftUnalignLen[1]; m++) {
            mt[n].leftUnalign[1][m] = ((*(pseq[1] + m)) >> 1) & 0x3;
        }
        mt[n].rightUnalignLen[1] = (n == 0) ? 0 : (n - 1); /* Subtract 1 because one character is added to the right for key alignment */
        for (m = 0; m < mt[n].rightUnalignLen[1]; m++) {
            mt[n].rightUnalign[1][m] = ((*(pseq[1] + baseLength - mt[n].rightUnalignLen[1] + m) >> 1) & 0x3);
        }
        actgSquash(pseq[1] + mt[n].leftUnalignLen[1], total, basePairSquashBuffer[n]);

        /* Establish the relationship between base squash and corresponding pair base squash */
        mt[n].set(baseSquashBuffer[n], squashLen[0], basePairSquashBuffer[n] + squashLen[1] - len_bgs, squashLen[1], 0);

        /* do mapping */
        align4_curr = n & 0x3;
        match_pair_origin = (*(pseq[0] + n + bg_mid) < *(pseq[1] + baseLength - baseGroupLen - n + bg_mid));

        psquash = mt[align4_curr].getSquash(match_pair_origin);
        xsquash = *((uint64_t *)(psquash));
        xsquash &= xsquash_tab[match_pair_origin];
        hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
        pos_vals = (uint32_t *)(pReference->queryPosition(hash32, pos_cnts));

        for (o = 0; o < pos_cnts; o++) {
            match_pos = *pos_vals++;
            match_pair = ((match_pos & 0x80000000) >> 31) ^ match_pair_origin;
            match_pos = (match_pos & 0x7FFFFFFF) << 3; /* to squash reference pos */

            /* check left and right boundary simply */
            if (match_pos + base_squash_align4 >= refe_squashlen || match_pos < base_squash_align4) {
                continue;
            }

            loffset = (match_pair) ? (mt[align4_curr].squashBufferLen[1] - mt[align4_curr].offset - len_bgs) : (mt[align4_curr].offset);

            /* caculate unmatch count */
            psquash = mt[align4_curr].getSquash(match_pair) - loffset;
            psquash_refe = prefe_squash + match_pos - loffset;

            xsquash_match = ((*((uint64_t *)(mt[align4_curr].getSquash(match_pair)))) & xsquash_tab[match_pair]);
            xsquash_macth_refe = ((*((uint64_t *)(prefe_squash + match_pos))) & xsquash_tab[match_pair]);
            if (xsquash_match != xsquash_macth_refe)  {
                /* key is not same, skip */
                continue;
            } 
            /* align 4 */
            unmatchs = actgSquashDiffCnt(psquash, psquash_refe, mt[align4_curr].squashBufferLen[match_pair]);
            if (unmatchs >= best_unmatchs) {
                continue;
            }

            /*  Because this case adds a character to the right for 32-byte alignment: need to handle the last character when matching to mt[0]'s pair and mt[0] offset is 0 */
            unmatchs -= (match_pair && (mt[align4_curr].offset == 0)) ? ((xsquash_match & 0x80000000000000) != (xsquash_macth_refe & 0x80000000000000)) : 0;

            /* left unalign */
            l = mt[align4_curr].leftUnalignLen[match_pair];
            for (ch = *(psquash_refe - 1), m = 0; m < mt[align4_curr].leftUnalignLen[match_pair]; m++) {
                unmatchs += ((ch >> (m << 1)) & 0x3) != (mt[align4_curr].leftUnalign[match_pair][l - 1]);
                l--;
            }
            /* right unalign */
            l = mt[align4_curr].squashBufferLen[match_pair];
            for (ch = *(psquash_refe + l), m = 0; m < mt[align4_curr].rightUnalignLen[match_pair]; m++) {
                unmatchs += ((ch >> (6 - (m << 1)) & 0x3) != (mt[align4_curr].rightUnalign[match_pair][m]));
            }
            if (unmatchs < best_unmatchs) {
                best_pos = (match_pos << 2) - (loffset << 2) - mt[align4_curr].leftUnalignLen[match_pair];
                best_pos_inrefe = match_pos - loffset;
                best_is_pair = match_pair;
                best_align4 = align4_curr;
                best_unmatchs = unmatchs;
            }
            if (best_unmatchs <= MAPPED_THRESHOLD_GEN2) {
                break;
            }
        }
        if (best_unmatchs <= MAPPED_THRESHOLD_GEN2) {
            break;
        }
        mt[align4_curr].incOffset();
    }

    if (best_unmatchs > MAPPED_THRESHOLD_GEN2) { /* continue mapping */
        for (n = 4; n <= e; n++) {
            /* do mapping */
            align4_curr = n & 0x3;
            match_pair_origin = (*(pseq[0] + n + bg_mid) < *(pseq[1] + baseLength - baseGroupLen - n + bg_mid));

            psquash = mt[align4_curr].getSquash(match_pair_origin);
            xsquash = *((uint64_t *)(psquash));
            xsquash &= xsquash_tab[match_pair_origin];
            hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
            pos_vals = (uint32_t *)(pReference->queryPosition(hash32, pos_cnts));

            for (o = 0; o < pos_cnts; o++) {
                match_pos = *pos_vals++;
                match_pair = ((match_pos & 0x80000000) >> 31) ^ match_pair_origin;
                match_pos = (match_pos & 0x7FFFFFFF) << 3; /* to squash reference pos */

                /* check left and right boundary simply */
                if (match_pos + base_squash_align4 >= refe_squashlen || match_pos < base_squash_align4) {
                    continue;
                }

                loffset = (match_pair) ? (mt[align4_curr].squashBufferLen[1] - mt[align4_curr].offset - len_bgs) : (mt[align4_curr].offset);
                /* caculate unmatch count */
                psquash = mt[align4_curr].getSquash(match_pair) - loffset;
                psquash_refe = prefe_squash + match_pos - loffset;

                xsquash_match = ((*((uint64_t *)(mt[align4_curr].getSquash(match_pair)))) & xsquash_tab[match_pair]);
                xsquash_macth_refe = ((*((uint64_t *)(prefe_squash + match_pos))) & xsquash_tab[match_pair]);
                if (xsquash_match != xsquash_macth_refe) { /* key is not same, skip */
                    continue;
                }

                /* align 4 */
                unmatchs = actgSquashDiffCnt(psquash, psquash_refe, mt[align4_curr].squashBufferLen[match_pair]);
                if (unmatchs >= best_unmatchs) {
                    continue;
                }

                /*  Because this case adds a character to the right for 32-byte alignment: need to handle the last character when matching to mt[0]'s pair and mt[0] offset is 0 */
                unmatchs -= (match_pair && (mt[align4_curr].offset == 0)) ? ((xsquash_match & 0x80000000000000) != (xsquash_macth_refe & 0x80000000000000)) : 0;

                /* left unalign */
                l = mt[align4_curr].leftUnalignLen[match_pair];
                for (ch = *(psquash_refe - 1), m = 0; m < mt[align4_curr].leftUnalignLen[match_pair]; m++) {
                    unmatchs += ((ch >> (m << 1)) & 0x3) != (mt[align4_curr].leftUnalign[match_pair][l - 1]);
                    l--;
                }
                /* right unalign */
                l = mt[align4_curr].squashBufferLen[match_pair];
                for (ch = *(psquash_refe + l), m = 0; m < mt[align4_curr].rightUnalignLen[match_pair]; m++) {
                    unmatchs += ((ch >> (6 - (m << 1)) & 0x3) != (mt[align4_curr].rightUnalign[match_pair][m]));
                }
                if (unmatchs < best_unmatchs) {
                    best_pos = (match_pos << 2) - (loffset << 2) - mt[align4_curr].leftUnalignLen[match_pair];
                    best_pos_inrefe = match_pos - loffset;
                    best_is_pair = match_pair;
                    best_align4 = align4_curr;
                    best_unmatchs = unmatchs;
                }
                if (best_unmatchs <= MAPPED_THRESHOLD_GEN2) {
                    break;
                }
            }

            if (best_unmatchs <= MAPPED_THRESHOLD_GEN2) {
                break;
            }
            mt[align4_curr].incOffset();
        }
    }

    /* calc the result of base or base pair mapping with the match pos reference */
    outLength = 0;
    if (best_unmatchs != UINT32_MAX)  {
        /* get pos in reference table */
        psquash = (best_is_pair) ? (mt[best_align4].squashBuffer[1] - (mt[best_align4].squashBufferLen[1] - len_bgs)) : (mt[best_align4].squashBuffer[0]);
        psquash_refe = prefe_squash + best_pos_inrefe;

        n = o = l = mt[best_align4].leftUnalignLen[best_is_pair];
        for (ch = *(psquash_refe - 1), m = 0; m < mt[best_align4].leftUnalignLen[best_is_pair]; m++) {
            out[--o] = (((ch >> (m << 1)) & 0x3) ^ (mt[best_align4].leftUnalign[best_is_pair][l - 1]));
            l--;
        }
        outLength += n; /* Note byte order */
        outLength += actgStretchMappingXor(psquash, psquash_refe, mt[best_align4].squashBufferLen[best_is_pair], out + outLength);

        outLength -= (best_is_pair && best_align4 == 0);
        l = mt[best_align4].squashBufferLen[best_is_pair];
        for (ch = *(psquash_refe + l), m = 0; m < mt[best_align4].rightUnalignLen[best_is_pair]; m++) {
            out[outLength++] = ((ch >> (6 - (m << 1)) & 0x3) ^ (mt[best_align4].rightUnalign[best_is_pair][m]));
        }
    } else {
        /* not match valid pos */
        actgEncode(base, out, baseLength);
        outLength = baseLength;
        best_pos = 0;
        best_is_pair = 2; ///*  During decompression, check mdir first; if it's 2, it means no match */
    }

    mappingPos = best_pos;
    mappingDir = best_is_pair;
    return;
}

void FastqActuator::mappingFastQGen3(const uint8_t*, uint32_t, uint8_t*&, uint32_t&, uint64_t&, uint8_t&) {
    LOG_ERROR("Not support FASTQ Gen3");
    return;
}

template <typename TCoder>
int32_t FastqActuator::compressIdStream(coder_io* idIo, TCoder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx) {
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

int32_t FastqActuator::compressIdInAll() {
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
    return outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
}

int32_t FastqActuator::compressIdInSplit() {
    Json::Value idMeta;
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    for (uint32_t i = 0; i < idSplitSymbols.size();++i) {
        // Variable length with all digits
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
            for (uint32_t j = 0; j < currLen; ++j) {
                if (*(data + j) < '0' || *(data + j) < '9') {
                    idDigit = false;
                    break;
                }
            }

            if (idDigit) {
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
    return outBlockPtr->getDataLen() > outBlockPtr->getBufferSize() ? -1 : 0;
}

int32_t FastqActuator::compressId() {
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

int32_t FastqActuator::compressBase() {
    if (pReference != nullptr) {
        return compressBaseWithRef();
    } else {
        return compressBaseWithoutRef();
    }
}

int32_t FastqActuator::compressBaseWithRef() {
    bool encBaseLen = (minBaseLength != maxBaseLength);
    uint32_t line = inBlockPtr->getNpos().size();
    uint8_t* ptr = inBlockPtr->getBuffer();
    int64_t nOffset = 0;
    int64_t currPos = 0;
    uint32_t offset = 0;
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

        for (uint32_t n = startPos; n < endPos; n++) {
            char ch = *(ptr + n);
            if (ch == 'n' || ch == 'N') {
                *(baseNPosBuffer + nOffset) = currPos + n - startPos;
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
    }
    metaBase["minlen"] = (Json::Value::UInt)minBaseLength;
    metaBase["maxlen"] = (Json::Value::UInt)maxBaseLength;
    metaBase["totalsrclen"] = (Json::Value::UInt)totalSrcLen;
    metaBase["totaldstlen"] = (Json::Value::UInt)totalDstLen;
    metaBase["streams"] = metaStreams;
    meta["base"] = metaBase;
    return 0;
}

int32_t FastqActuator::compressBaseWithoutRef() {
    int32_t addLength = !!(minBaseLength != maxBaseLength);
    uint32_t baseSrcLength = 0;
    for (uint32_t idx = 1; idx < inBlockPtr->getNpos().size(); idx +=4) {
        int64_t end = inBlockPtr->getNpos()[idx];
        int64_t start = inBlockPtr->getNpos()[idx-1] + 1;
        baseSrcLength += end - start + addLength;
    }

    std::shared_ptr<coder_io> baseIo;
    if (baseSrcLength <= FC_MIN_LEN || baseSrcLength > FC_MAX_LEN) {
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
    return 0;
}

int32_t FastqActuator::compressComment() {
    std::shared_ptr<coder_io> commentIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_fc> commentCoder = std::make_shared<coder_fc>(commentIo.get());
    
    Json::Value commentMeta;
    switch (commentType) {
    case CommentType::PLUS_ONLY:
        commentMeta["type"] = "plusonly";
        break;
    case CommentType::SAME_AS_ID:
        commentMeta["type"] = "sameasid";
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
        break;
    }
    default:
        break;
    }
    meta["comment"] = commentMeta;
    return 0;
}

int32_t FastqActuator::compressQuality() {
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

    return 0;
}

int32_t FastqActuator::initDecoder(RoughIOBlock* outputBlock) {
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
            idDecoders.push_back(new coder_affix_match(new coder_io(inBlockPtr->getBuffer() + readOffset, dstLen)));
            idDecoders.back()->set_level(idStreamMeta[idx]["coder"]["level"].asInt());
        } else if (coderName == "coder_bwt_cm") {
            idDecoders.push_back(new coder_bwt_cm(new coder_io(inBlockPtr->getBuffer() + readOffset, dstLen)));
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
    if (pReference == nullptr) {
        std::string baseCoderName = baseMeta["coder"]["magic"].asString();
        uint32_t srcBaseLen = baseMeta["totalsrclen"].asUInt();
        uint32_t dstBaseLen = baseMeta["totaldstlen"].asUInt();
        if (baseCoderName == "coder_fc") {
            coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstBaseLen);
            baseIo.meta = baseMeta;
            baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
            coder_fc baseDecoderFc = coder_fc(&baseIo);
            // Put decompression result at the end
            baseDecoderFc.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcBaseLen, srcBaseLen,  UINT8_MAX, false);
        } else if (baseCoderName == "coder_bwt_cm") {
            baseDecoder = new coder_bwt_cm(new coder_io(inBlockPtr->getBuffer() + readOffset, dstBaseLen));
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
            baseDecoder = new coder_bwt_cm(new coder_io(inBlockPtr->getBuffer() + readOffset, baseDstLen));
        } else {
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }
        readOffset += baseDstLen;

        /* check sub streams 2 */
        id = 1;
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
            coder_bwt_cm posCm(&posIo);
            posCm.decode_line((uint8_t *)baseMappedPosBuffer, srcLen, UINT8_MAX, false);
        } else if (metaStreams[id]["coder"]["magic"].asString() == "coder_fc") {
            temBuffer =  inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPosBuffer = (uint64_t *)ps;
            ps += srcLen;
            coder_io posIo(temBuffer, dstLen);
            posIo.meta = metaStreams[id];
            posIo.meta["dstlen"] = metaStreams[id]["dstlen"]; // Compatible field, can be unified
            coder_fc posCm(&posIo);
            posCm.decode_line((uint8_t*)baseMappedPosBuffer, srcLen, UINT8_MAX, false);
        } else {
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }
            
        readOffset += metaStreams[id]["dstlen"].asUInt();

        /* check sub stream 3 */
        id = 2;
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
            coder_bwt_cm pairCm(&pairIo);
            pairCm.decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else if (metaStreams[id]["coder"]["magic"].asString() == "coder_fc") {
            temBuffer = inBlockPtr->getBuffer() + readOffset;
            srcLen = metaStreams[id]["srclen"].asUInt();
            dstLen = metaStreams[id]["dstlen"].asUInt();
            baseMappedPairBuffer = ps;
            ps += srcLen;
            coder_io pairIo(temBuffer, dstLen);
            pairIo.meta = metaStreams[id];
            pairIo.meta["tot_dstlen"] = metaStreams[id]["dstlen"]; // Compatible field, can be unified
            coder_fc pairCm(&pairIo);
            pairCm.decode_line(baseMappedPairBuffer, srcLen, UINT8_MAX, false);
        } else{
            LOG_ERROR("check sub stream failed: coder name unmatch");
            return -1;
        }        
        readOffset += metaStreams[id]["dstlen"].asUInt();

        /* check sub stream 4 */
        baseNPosBuffer = nullptr;
        if (baseMeta["ncount"].asUInt()) {
            id = 3;
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
                coder_bwt_cm nposCm(&nposIo);
                nposCm.decode_line((uint8_t *)baseNPosBuffer, srcLen, UINT8_MAX, false);
            } else {
                LOG_ERROR("check sub stream failed: coder name unmatch");
                return -1;
            }
            readOffset += metaStreams[id]["dstlen"].asUInt();
        }

        /* check sub stream 5 */
        baseLengthGen2Buffer = nullptr;
        if (baseMeta["minlen"].asUInt() != baseMeta["maxlen"].asUInt()) {
            id = 4;
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
                coder_bwt_cm lenCm(&lenIo);
                lenCm.decode_line((uint8_t *)baseLengthGen2Buffer, srcLen, UINT8_MAX, false);
            } else  {
                LOG_ERROR("check sub stream failed: coder name unmatch");
                return -1;
            }
            readOffset += metaStreams[id]["dstlen"].asUInt();
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
            commentDecoder = new coder_bwt_cm(new coder_io(inBlockPtr->getBuffer() + readOffset, commentDstLen));
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
        coder_bwt_cm qualFreqCoder(&qualFreqIo);

        uint32_t qualFreqSrcLen = streamsMeta[1]["srclen"].asUInt();
        uint8_t qualFreqArrLen = qualFreqSrcLen / sizeof(uint16_t);
        uint16_t* qualFreqArray = new uint16_t[qualFreqArrLen];
        uint32_t qualFreq = qualFreqCoder.decode_line((uint8_t*)qualFreqArray, qualFreqSrcLen, UINT8_MAX, false);
        if (qualFreq != qualFreqSrcLen) {
            LOG_ERROR("Decode quality frequncy failed.");
            return -1;
        }

        for (uint32_t idx = 0; idx < qualFreqArrLen ; idx += 2) {
            qualityFreqTable.push_back(std::make_pair(qualFreqArray[idx], qualFreqArray[idx + 1]));
        }

        if (streamsMeta[0]["coder"]["magic"].asString() == "coder_qual") {
            qualityDecoder = new coder_qual(new coder_io(inBlockPtr->getBuffer() + readOffset, qualityDstLen),
                                       true, qualityFreqTable);
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

int32_t FastqActuator::decompress() {
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
        if (pReference == nullptr) {
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
