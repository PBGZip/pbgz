#pragma once

#include <string>

#include "io_block.h"
#include "reference.h"
#include "actuator.h"
#include "coder.h"
#include "coder_io.h"
#include "coder_qual.h"


enum class CommentType {
    PLUS_ONLY,     // Comment line contains only a '+' sign
    SAME_AS_ID,    // Comment line is the same as ID line
    OTHER,         // Comment line is neither of the above, could be any other string including empty string
    UNKNOWN        // Unknown type
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
        refeBeginPos = 0;
        refeEndPos = 0;
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

    std::vector<int64_t>& getRefeMappedPositons() {
        return refeMappedPositons;
    }

    void setRefeBeginEndPos(int64_t beginPos, int64_t endPos) {
        refeBeginPos = beginPos;
        refeEndPos = endPos; 
    }

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

    int32_t initDecoder(RoughIOBlock* outputBlock);

    void (FastqActuator::*mapping)(const uint8_t*, uint32_t, uint8_t*&, uint32_t&, uint64_t&, uint8_t&);

    void mappingFastqGen2(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir);

    void mappingFastQGen3(const uint8_t* base, uint32_t baseLength, uint8_t*& out, uint32_t& outLength, uint64_t& mappingPos, uint8_t& mappingDir);

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

    std::vector<coder*> idDecoders;
    coder* baseDecoder;
    coder* commentDecoder;
    coder_qual* qualityDecoder;

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
    uint16_t* baseLengthGen2Buffer;
    uint32_t* baseLengthGen3Buffer;
    uint8_t* refeStretchBuffer;

    std::vector<int64_t> refeMappedPositons;
    uint64_t refeBeginPos;
    uint64_t refeEndPos;
};
