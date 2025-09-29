#include <cstring>
#include <algorithm>
#include <memory>

#include "fastq_actuator.h"
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


int32_t FastqActuator::compress() {
    if (0 != initEncoder()) {
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

    //计算数据块的Md5
    std::string md5;
    calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
    meta["md5"] = md5;
    meta["idlines"] = inBlockPtr->getNpos().size() >> 2;

    // 压缩Meta信息
    coder_json metaCoder;
    int32_t metaLen = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
    if (metaLen <= 0) {
        LOG_ERROR("Failed to encode meta information");
        return -1;
    }
    outBlockPtr->setMetaLen(metaLen);

    return 0;
}

/// @brief 解析首行存在的分隔符以及分隔符所在的位置
/// @param pBuffer    // ID行的buffer
/// @param bufLen     // ID行的长度
/// @param idSplitSymbols    // 分隔符列表
/// @return 0 成功，-1失败
int32_t FastqActuator::preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLen) {
    if (pBuffer == nullptr || bufLen == 0) {
        return -1;
    }

    std::vector<uint32_t> idSplitPos;  
    for (uint32_t i = 1 ; i < bufLen; ++i) {    // 首个字符为@，跳过
        char ch = pBuffer[i];
        if (idSplitDefault.find(ch) != std::string::npos) {
            idSplitSymbols.push_back(ch);
            idSplitPos.push_back(static_cast<uint8_t>(i));
        }
    }

    // 初始化每个分隔符的最大和最小长度
    for (size_t i = 0; i < idSplitSymbols.size(); ++i) {
        idSplitMinLen.push_back(UINT32_MAX);
        idSplitMaxLen.push_back(0);
    }

    // 将第一行的ID行分析信息拷贝到idPositions中
    uint32_t lastPos = 0;
    for (uint32_t idx = 0;  idx < idSplitPos.size(); ++idx) {
        uint32_t pos = idSplitPos[idx]; 
        uint32_t curLen = pos - lastPos - (0 == idx ? 0 : 1);   // 首行不需要偏移
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
            idPosLength = UINT32_MAX; // 标记为不可用
            break;
        }
    }
    return 0;
}

int32_t FastqActuator::preAnalysisBase(uint8_t* pBuffer, uint32_t bufLen) {
    uint32_t baseLength = bufLen - 1;    // 去掉换行符号
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
        if (commentType == CommentType::UNKNOWN) {  // 首行
            if (*pBuffer == '+' && bufLen == 2) {
                commentType = CommentType::PLUS_ONLY;
            }
            else {
                // 找到ID行的内容
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
            // 找到ID行的内容
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
        uint32_t lineLength  = endPos - startPos + 1;   // 长度需要带换行符
        switch (lineNo & 0x3)
        {
        case 0: {  // ID行
            if (idPosLength != UINT32_MAX) {
                if (lineNo == 0) {  // 首行
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
        case 1: { //  碱基行
            if (preAnalysisBase(inBlockPtr->getBuffer() + startPos, lineLength) != 0) {
                return -1;
            }
            break;
        }
        case 2: { // 注释行
            if (preAnalysisComment(inBlockPtr->getBuffer() + startPos, lineLength, lineNo) != 0) { 
                return -1;
            }
            break;
        }
        case 3: { // 质量值行
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
    return 0;
}

template <typename TCoder>
int32_t FastqActuator::compressIdStream(coder_io* idIo, TCoder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx) {
    int32_t currLineOffset = 0;    // 行偏移，即每行的开头
    uint8_t * data = nullptr;
    int32_t currLen = 0;
    srcDataLen = 0;   // 源内容的总长度
    for (uint32_t idx = 0; idx < inBlockPtr->getNpos().size(); idx += 4) {
        if (splitSymIdx == 0) {
            data = inBlockPtr->getBuffer() + currLineOffset;
            currLen =  idPositions[splitSymIdx + idx] + 1;  //  currIdPos;  
        } else {
            data = inBlockPtr->getBuffer() + currLineOffset + idPositions[splitSymIdx + idx - 1] + 1;
            currLen = idPositions[splitSymIdx + idx] - idPositions[splitSymIdx + idx - 1];
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
        // 全是数字的变长
        if (idSplitMaxLen[i] != idSplitMinLen[i]) {
            uint32_t currIdPos = idPositions[i];
            uint32_t currLineOffset = 0;
            // 根据首行判断是否全是数字 
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

        // 定长或者不是全是数字场景
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

    // ID检查不合法，整块压缩
    if (idPosLength == UINT32_MAX) {
        return compressIdInAll();
    }else {
        // ID检查合法，分块压缩
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

    // 编码质量数据
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

    // 编码质量符号表
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

int32_t FastqActuator::initDecoder() {
    Json::Value& idStreamMeta = meta["id"]["streams"];
    if (idStreamMeta.size() != idSplitSymbols.size()) {
        LOG_ERROR("id stream meta size not match id split symbols size(%d, %d)", idStreamMeta.size(), idSplitSymbols.size());
        return -1;
    }

    uint32_t readOffset = 0;
    // 初始化id解码器
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

    // 初始化Base解码器
    Json::Value& baseMeta = meta["base"];
    minBaseLength = baseMeta["minlen"].asUInt();
    maxBaseLength = baseMeta["maxlen"].asUInt();
    if (pReference == nullptr) {
        std::string baseCoderName = baseMeta["coder"]["magic"].asString();
        uint32_t srcBaseLen = baseMeta["totalsrclen"].asUInt();
        uint32_t dstBaseLen = baseMeta["totaldstlen"].asUInt();
        if (baseCoderName == "coder_fc") {
            coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstBaseLen);
            baseIo.meta = baseMeta;
            baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
            coder_fc baseDecoderFc = coder_fc(&baseIo);
            // 解压结果放在最后
            baseDecoderFc.decode_line(outBlockPtr->getBuffer() + outBlockPtr->getBufferSize() - srcBaseLen, srcBaseLen,  UINT8_MAX, false);
        } else if (baseCoderName == "coder_bwt_cm") {
            baseDecoder = new coder_bwt_cm(new coder_io(inBlockPtr->getBuffer() + readOffset, dstBaseLen));
        } else {
            LOG_ERROR("unsupported coder name: %s", baseCoderName.c_str());
            return -1;
        }

        readOffset += dstBaseLen;
    } 

    // 初始化Comment解码器
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

    // 初始化Quality解码器
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
        coder_bwt_cm qualFreaCoder(&qualFreqIo);

        uint32_t qualFreqSrcLen = streamsMeta[1]["srclen"].asUInt();
        uint8_t qualFreqArrLen = qualFreqSrcLen / sizeof(uint16_t);
        uint16_t* qualFreqArray = new uint16_t[qualFreqArrLen];
        uint32_t qualFreq = qualFreaCoder.decode_line((uint8_t*)qualFreqArray, qualFreqSrcLen, UINT8_MAX, false);
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

    // 解析meta
    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);

    std::string idSplit = meta["id"]["splitsym"].asString();
    for (uint32_t i = 0; i < idSplit.length(); ++i) {
        idSplitSymbols.push_back(idSplit.c_str()[i]);
    }

    if (0 != initDecoder()) {
        LOG_ERROR("Initital decoder failed.");
        return -1;
    }

    uint32_t lineNum = meta["idlines"].asUInt();
    uint8_t* pBaseOut = nullptr;
    uint8_t* pBaseEnd = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize();
    if (meta["base"]["coder"]["magic"].asString() == "coder_fc") {
        pBaseOut = pBaseEnd - meta["base"]["totalsrclen"].asUInt();
    }

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
                actualBaseLen -= 1;  // 去掉换行符
            }
        } else {
            // 补充支持参考基因场景
        }

        // 解码comment
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

        // 解码质量值
        qualityDecoder->decode_qual_gen2(basePtr, outBlockPtr->getCurrent(), actualBaseLen);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
        *outBlockPtr->getCurrent() = '\n';
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
    }

    // uint32_t outBaseLen = 0;
    // if (pReference == nullptr) {
    //     outBaseLen = meta["base"]["totalsrclen"].asUInt() - (minBaseLength == maxBaseLength ? 0 : lineNum);
    // } else {
    //     outBaseLen = meta["base"]["streams"][0]["srclen"].asUInt() + baseNCount;
    // }

    // uint32_t outQualityLen = meta["quality"]["streams"][0]["srclen"].asUInt();
    // if (outBaseLen == outQualityLen) {
    //     outBlockPtr->setDataLen(outBlockPtr->getDataLen() - 1);
    // } else {
    //     LOG_ERROR("check quality length faild, base length = %d, quality length = %d", outBaseLen, outQualityLen);
    //     return -1;
    // }

    // 检查文件的MD5
    std::string md5;
    calcMd5sum(md5, outBlockPtr->getBuffer(), outBlockPtr->getDataLen());
    if (md5 != meta["md5"].asString()) {
        LOG_ERROR("Md5 check failed, expect %s, got %s.", meta["md5"].asCString(), md5.c_str());
        return -1;
    }

    return 0;
}
