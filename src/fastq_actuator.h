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

typedef struct {
    void set(uint8_t *_squashBuffer, uint32_t _squashBufferLen, uint8_t *_squashBuffCmplt, uint32_t __squashBuffCmpltLen, uint32_t _offset) {
        squashBuffer[0] = _squashBuffer;
        squashBufferLen[0] = _squashBufferLen;
        squashBuffer[1] = _squashBuffCmplt;
        squashBufferLen[1] = __squashBuffCmpltLen;
        offset = _offset;
    }
    inline void incOffset() { ++offset; }
    uint8_t *getSquash(bool pair) const {return ((pair) ? (squashBuffer[1] - offset) : (squashBuffer[0] + offset));}

    /*  squash buffer , 0表示正向，1表示互补链 */
    uint8_t *squashBuffer[2];      /* 指向squash buffer */
    uint32_t squashBufferLen[2];  /* 对应squash buffer的数据长度 */
    uint32_t leftUnalignLen[2]; /* 左边没有对齐的碱基长度 */
    uint8_t leftUnalign[2][3]; /* 左边没有对齐的碱基，这里存的是squash之后的数据 */
    uint32_t rightUnalignLen[2]; /* 右边没有对齐的碱基长度 */
    uint8_t rightUnalign[2][3]; /* 右边没有对齐的碱基，这里存的是squash之后的数据 */
    uint32_t offset;  /* 当前偏移 */
} Mapping;


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
        baseLengthGen2Buffer = nullptr;
        baseLengthGen3Buffer = nullptr;
        refeStretchBuffer = nullptr;
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

    int32_t initEncoder();

private:

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

    void (FastqActuator::*mapping)(const uint8_t*, uint32_t, uint8_t*&, uint32_t&, uint64_t&, uint8_t&);

    void mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir);

    void mappingFastQGen3(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir);

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

    bool isGen2;

    /// 带参考基因场景使用的变量
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
    uint16_t* baseLengthGen2Buffer;
    uint32_t* baseLengthGen3Buffer;
    uint8_t* refeStretchBuffer;
};
