#pragma once

#include <string>

#include "io_block.h"
#include "reference.h"
#include "actuator.h"
#include "coder.h"
#include "coder_io.h"
#include "coder_qual.h"


enum class CommentType {
    PLUS_ONLY,     // Comment行只有一个+号
    SAME_AS_ID,    // Comment行和ID行相同
    OTHER,         // Comment行非上述两种情况，即可能为其他任意字符串，包括空字符串
    UNKNOWN        // 未知的类型
};


class FastqActuator : public Actuator {
public:
    FastqActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, Reference* pRef = nullptr): Actuator(inPtr, outPtr) {
        idPosLength = 0;
        minBaseLength = INT32_MAX;
        maxBaseLength = 0;
        pReference = pRef;
        baseNCount = 0;
        commentType = CommentType::UNKNOWN;
        baseDecoder = nullptr;
        commentDecoder = nullptr;
        qualityDecoder = nullptr;
    }

    virtual ~FastqActuator() {
        for (auto ptr : idDecoders) {
            delete ptr;
        }

        idDecoders.clear();
        if (baseDecoder != nullptr) {
            delete baseDecoder; 
            baseDecoder = nullptr;
        }

        if (commentDecoder != nullptr) {
            delete commentDecoder;
            commentDecoder = nullptr;
        }

        if (qualityDecoder != nullptr) {
            delete qualityDecoder;
            qualityDecoder = nullptr;
        }

    }

    int32_t decompress() override ;

    int32_t compress() override ;

    int32_t preAnalysis();

private:
    int32_t initEncoder();

    int32_t preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLe );

    int32_t preAnalysisId(uint8_t* pBuffer, uint32_t bufLen);

    int32_t preAnalysisBase(uint8_t* pBuffer, uint32_t bufLen);

    int32_t preAnalysisComment(uint8_t* pBuffer, uint32_t bufLen, uint32_t lineNo);

    int32_t compressId();

    int32_t compressIdInAll();

    int32_t compressIdInSplit();

    template <typename TCoder>
    int32_t compressIdStream(coder_io* idIo, TCoder* idCoder, Json::Value& streamMeta, uint32_t& srcDataLen, int32_t splitSymIdx);

    int32_t compressBase();

    int32_t compressBaseWithRef();

    int32_t compressBaseWithoutRef();

    int32_t compressComment();

    int32_t compressQuality();

    int32_t initDecoder();

private:
    std::vector<uint8_t> idSplitSymbols;
    uint32_t idPosLength;
    std::vector<uint16_t> idPositions; // 每条记录ID字段中各个分隔符的位置，从0开始
    std::vector<uint32_t> idSplitMinLen; // 每个分隔符的最小长度
    std::vector<uint32_t> idSplitMaxLen; // 每个分隔符的最大长度
    const std::string idSplitDefault = "/:= _.,-#\r\t\n";
    uint32_t minBaseLength;
    uint32_t maxBaseLength;
    uint64_t baseNCount;
    CommentType commentType;
    std::vector<std::pair<uint16_t, uint16_t>> qualityFreqTable; // 质量值频率统计
    Reference* pReference;

    std::vector<coder*> idDecoders;
    coder* baseDecoder;
    coder* commentDecoder;
    coder_qual* qualityDecoder;
};


