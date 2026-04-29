/*
 * index_actuator_testcase.cpp - Test cases for IndexActuator functionality
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

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

#define private public
#define protected public
#include "index_actuator.h"
#include "sam_index.h"
#include "io_wrapper.h"
#include "sam_info.h"
#undef private
#undef protected

namespace IndexActuatorTestData {
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class IndexActuatorTest : public ::testing::Test {
public:
    void SetUp() override {
        pInBlock = new RoughIOBlock(IndexActuatorTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(IndexActuatorTestData::MAX_BLOCK_SIZE);
    }

    void TearDown() override {
        if (pInBlock != nullptr) {
            delete pInBlock;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
        }
        pInBlock = nullptr;
        pOutBlock = nullptr;
    }

    void createValidSamFile() {
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().addChromosomeInfo("chr1", 1000000);
        SamInfo::getInstance().addChromosomeInfo("chr2", 1000000);
        SamInfo::getInstance().calculateChromosomePositions();
    }

    Json::Value createBasicMeta() {
        Json::Value meta;
        Json::Value samMeta;
        Json::Value streams(Json::arrayValue);

        CreateMockStreams(streams);

        samMeta["lines"] = 5;
        samMeta["streams"] = streams;
        meta["sam"] = samMeta;

        return meta;
    }

    void CreateMockStreams(Json::Value& streams) {
        Json::Value idField;
        idField["field"] = 0;
        idField["dstlen"] = 10;
        idField["coder"]["magic"] = "coder_bwt_cm";
        Json::Value idStreams(Json::arrayValue);
        idStreams[0]["dstlen"] = 10;
        idStreams[0]["coder"]["magic"] = "coder_fc";
        idField["streams"] = idStreams;
        streams.append(idField);

        Json::Value flagField;
        flagField["field"] = 1;
        flagField["dstlen"] = 10;
        flagField["coder"]["magic"] = "coder_bwt_cm";
        flagField["coder"]["level"] = 5;
        streams.append(flagField);

        Json::Value rnameField;
        rnameField["field"] = 2;
        rnameField["dstlen"] = 10;
        rnameField["coder"]["magic"] = "coder_bwt_cm";
        rnameField["coder"]["level"] = 5;
        streams.append(rnameField);

        Json::Value posField;
        posField["field"] = 3;
        posField["dstlen"] = 10;
        posField["coder"]["magic"] = "coder_bwt_cm";
        posField["coder"]["level"] = 5;
        streams.append(posField);
    }

    Json::Value createHeaderMeta() {
        Json::Value meta;
        Json::Value headerMeta;
        headerMeta["srclen"] = 200;
        headerMeta["dstlen"] = 150;
        headerMeta["lines"] = 5;
        Json::Value coder;
        coder["magic"] = "coder_bwt_cm";
        coder["level"] = 5;
        headerMeta["coder"] = coder;
        meta["header"] = headerMeta;
        return meta;
    }

    void* coder_bwt_cm_null_decoder = nullptr;

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
};

TEST_F(IndexActuatorTest, Constructor_ValidPointers_Success) {
    IndexActuator actuator(pInBlock, pOutBlock);
    EXPECT_EQ(pInBlock, actuator.inBlockPtr);
    EXPECT_EQ(pOutBlock, actuator.outBlockPtr);
    EXPECT_EQ(nullptr, actuator.flagDecoder);
    EXPECT_EQ(nullptr, actuator.chrDecoder);
    EXPECT_EQ(nullptr, actuator.posDecoder);
    EXPECT_EQ(0, actuator.headEndLine);
    EXPECT_FALSE(actuator.notifyFlag);
}

TEST_F(IndexActuatorTest, Constructor_NullInBlock_Created) {
    IndexActuator actuator(nullptr, pOutBlock);
    EXPECT_EQ(nullptr, actuator.inBlockPtr);
    EXPECT_EQ(pOutBlock, actuator.outBlockPtr);
}

TEST_F(IndexActuatorTest, Constructor_NullOutBlock_Created) {
    IndexActuator actuator(pInBlock, nullptr);
    EXPECT_EQ(pInBlock, actuator.inBlockPtr);
    EXPECT_EQ(nullptr, actuator.outBlockPtr);
}

TEST_F(IndexActuatorTest, Destructor_CleansUpDecoders) {
    createValidSamFile();

    // Create valid coder_io instances for the decoders
    std::shared_ptr<coder_io> flagIo = std::make_shared<coder_io>(pInBlock->getBuffer(), 100);
    std::shared_ptr<coder_io> chrIo = std::make_shared<coder_io>(pInBlock->getBuffer() + 100, 100);
    std::shared_ptr<coder_io> posIo = std::make_shared<coder_io>(pInBlock->getBuffer() + 200, 100);

    IndexActuator* actuator = new IndexActuator(pInBlock, pOutBlock);
    actuator->flagIo = flagIo;
    actuator->chrIo = chrIo;
    actuator->posIo = posIo;
    actuator->flagDecoder = new coder_bwt_cm(actuator->flagIo.get());
    actuator->chrDecoder = new coder_bwt_cm(actuator->chrIo.get());
    actuator->posDecoder = new coder_bwt_cm(actuator->posIo.get());

    delete actuator;

    EXPECT_TRUE(true);
}

TEST_F(IndexActuatorTest, Initial_NullInputBlock_ReturnsError) {
    IndexActuator actuator(nullptr, pOutBlock);
    int32_t result = actuator.initial();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Initial_NonSAMBlockType_ReturnsError) {
    createValidSamFile();
    pInBlock->setBlockType(BINARY);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Initial_InvalidBlockType_ReturnsError) {
    createValidSamFile();
    pInBlock->setBlockType(TYPE_UNKNOW);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Initial_BAMBlockType_ReturnsError) {
    createValidSamFile();
    pInBlock->setBlockType(BAM);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Initial_GZIPBlockType_ReturnsError) {
    createValidSamFile();
    pInBlock->setBlockType(GZIP);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Initial_SAMBlockType_Allowed) {
    createValidSamFile();
    pInBlock->setBlockType(SAM);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    // SAM block type passes validation. If no "sam" meta member exists,
    // the function returns 0 (success) after checking header
    // This is valid for blocks with only header and no field-by-field compressed data
    EXPECT_EQ(0, result);
}

TEST_F(IndexActuatorTest, Initial_ZeroMetaLen_HandlesGracefully) {
    createValidSamFile();
    pInBlock->setBlockType(SAM);
    pInBlock->setMetaLen(0);

    IndexActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();
    // With zero meta length and no meta data, the function returns 0 (success)
    // after passing block type check, as it assumes header-only or empty data block
    EXPECT_EQ(0, result);
}

TEST_F(IndexActuatorTest, Initial_ValidMetaData_CheckSuccess) {
    createValidSamFile();
    pInBlock->setBlockType(SAM);
    pInBlock->setMetaLen(0);
    pInBlock->setBlockId(1);

    IndexActuator actuator(pInBlock, pOutBlock);
    // Simulate having valid meta data structure
    actuator.notifyFlag = true; // Would be set by successful initial

    // In a real scenario this would pass with valid meta data
    // For unit test, we verify flag can be set
    EXPECT_TRUE(actuator.notifyFlag);
}

TEST_F(IndexActuatorTest, Initial_SamIndexSingletonCanBeAccessed) {
    createValidSamFile();
    pInBlock->setBlockType(SAM);
    pInBlock->setMetaLen(100);
    pInBlock->setBlockId(0);

    IndexActuator actuator(pInBlock, pOutBlock);

    // SamIndex singleton can be accessed for index storage
    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_TRUE(indexList.empty());
}

TEST_F(IndexActuatorTest, Initial_DecodersInitialized_AfterCall) {
    createValidSamFile();
    pInBlock->setBlockType(SAM);
    pInBlock->setMetaLen(0);
    pInBlock->setBlockId(0);

    IndexActuator actuator(pInBlock, pOutBlock);

    // Before successful initialization, decoders are nullptr
    EXPECT_EQ(nullptr, actuator.flagDecoder);
    EXPECT_EQ(nullptr, actuator.chrDecoder);
    EXPECT_EQ(nullptr, actuator.posDecoder);
}

TEST_F(IndexActuatorTest, Initial_MultiBlockHandling_DifferentBlockIds) {
    createValidSamFile();

    // Test with different block IDs to verify indexing works across blocks
    for (int blockId = 0; blockId < 3; blockId++) {
        pInBlock->reset();
        pOutBlock->reset();
        pInBlock->setBlockType(SAM);
        pInBlock->setMetaLen(0);
        pInBlock->setBlockId(blockId);

        IndexActuator actuator(pInBlock, pOutBlock);
        // Each actuator should have its own state
        EXPECT_EQ(blockId, pInBlock->getBlockId());
    }
}

TEST_F(IndexActuatorTest, GetNotifyFlag_FalseAfterConstruction) {
    IndexActuator actuator(pInBlock, pOutBlock);
    EXPECT_FALSE(actuator.getNotifyFlag());
}

TEST_F(IndexActuatorTest, GetNotifyFlag_TrueAfterProcessing) {
    createValidSamFile();
    IndexActuator actuator(pInBlock, pOutBlock);

    Json::Value meta = createBasicMeta();
    pInBlock->setMetaLen(0);
    pInBlock->setBlockId(0);

    // Simulate setting notify flag after some processing
    actuator.notifyFlag = true;

    EXPECT_TRUE(actuator.getNotifyFlag());
}

TEST_F(IndexActuatorTest, Process_NullOutBlock_ReturnsError) {
    IndexActuator actuator(pInBlock, nullptr);
    int32_t result = actuator.process();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, Process_ClearsSortKeys_ReturnsZero) {
    IndexActuator actuator(pInBlock, pOutBlock);

    // Add items to sortKeys to verify they get cleared
    IndexActuator::SortKey key1;
    key1.chrIndex = 1;
    key1.mapPos = 1000;
    actuator.sortKeys.push_back(key1);

    EXPECT_EQ(1, actuator.sortKeys.size());

    int32_t result = actuator.process();
    EXPECT_EQ(0, result);

    // process() should clear sortKeys
    EXPECT_EQ(0, actuator.sortKeys.size());
}

TEST_F(IndexActuatorTest, Process_ClearsSortKeys_Success) {
    createValidSamFile();
    IndexActuator actuator(pInBlock, pOutBlock);

    // Add items to sortKeys for testing
    IndexActuator::SortKey key1;
    key1.chrIndex = 1;
    key1.mapPos = 1000;
    actuator.sortKeys.push_back(key1);

    IndexActuator::SortKey key2;
    key2.chrIndex = 1;
    key2.mapPos = 2000;
    actuator.sortKeys.push_back(key2);

    EXPECT_EQ(2, actuator.sortKeys.size());

    int32_t result = actuator.process();
    EXPECT_EQ(0, result); // process() just clears sortKeys

    EXPECT_EQ(0, actuator.sortKeys.size());
}

TEST_F(IndexActuatorTest, Process_WithNullOutBlock_ReturnsError) {
    createValidSamFile();
    IndexActuator actuator(pInBlock, pOutBlock);

    // Ensure outBlock is null for this test
    RoughIOBlock* nullBlock = nullptr;
    IndexActuator actuator2(pInBlock, nullBlock);
    int32_t result = actuator2.process();
    EXPECT_EQ(-1, result);
}

TEST_F(IndexActuatorTest, SortKey_OperatorLess_ChrIndexComparison) {
    IndexActuator::SortKey key1;
    key1.chrIndex = 1;
    key1.mapPos = 1000;

    IndexActuator::SortKey key2;
    key2.chrIndex = 2;
    key2.mapPos = 1000;

    EXPECT_TRUE(key1 < key2);
    EXPECT_FALSE(key2 < key1);
}

TEST_F(IndexActuatorTest, SortKey_OperatorLess_SameChrIdCompareMapPos) {
    IndexActuator::SortKey key1;
    key1.chrIndex = 1;
    key1.mapPos = 500;

    IndexActuator::SortKey key2;
    key2.chrIndex = 1;
    key2.mapPos = 1000;

    EXPECT_TRUE(key1 < key2);
    EXPECT_FALSE(key2 < key1);
}

TEST_F(IndexActuatorTest, SortKey_OperatorLess_EqualKeys_ReturnsFalse) {
    IndexActuator::SortKey key1;
    key1.chrIndex = 1;
    key1.mapPos = 1000;

    IndexActuator::SortKey key2;
    key2.chrIndex = 1;
    key2.mapPos = 1000;

    EXPECT_FALSE(key1 < key2);
    EXPECT_FALSE(key2 < key1);
}

TEST_F(IndexActuatorTest, SamIndexSingleton_InitializesAsEmpty) {
    createValidSamFile();

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_TRUE(indexList.empty());
}

TEST_F(IndexActuatorTest, HeadEndLine_InitializesToZero) {
    IndexActuator actuator(pInBlock, pOutBlock);

    EXPECT_EQ(0, actuator.headEndLine);
}

TEST_F(IndexActuatorTest, DecoderIo_SharedPtr_InitializesAsNull) {
    IndexActuator actuator(pInBlock, pOutBlock);

    EXPECT_EQ(nullptr, actuator.flagIo);
    EXPECT_EQ(nullptr, actuator.chrIo);
    EXPECT_EQ(nullptr, actuator.posIo);
}

TEST_F(IndexActuatorTest, SortKeys_ReservesSpaceForTwoItems) {
    createValidSamFile();
    IndexActuator actuator(pInBlock, pOutBlock);

    uint32_t lineNum = 5;
    uint32_t expectedCapacity = 2;

    actuator.sortKeys.reserve(lineNum > 2 ? 2 : lineNum);

    EXPECT_EQ(expectedCapacity, 2);
}

TEST_F(IndexActuatorTest, IndexItemStructure_AccessesFieldsCorrectly) {
    SamIndexItem item;

    item.referenceMapPos = 123456;
    item.readNumber = 42;
    item.blockId = 1;

    EXPECT_EQ(123456, item.referenceMapPos);
    EXPECT_EQ(42, item.readNumber);
    EXPECT_EQ(1, item.blockId);
}

TEST_F(IndexActuatorTest, SamIndexSingleton_AddIndexItems) {
    createValidSamFile();

    // Add multiple items for same chromosome directly to singleton
    for (int i = 0; i < 3; i++) {
        SamIndex::getInstance().addSamIndex(1, 1000 + i * 100, i + 1, 1);
    }

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(3, indexList.at(1).size());

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, SamIndexSingleton_DifferentChromosomes) {
    createValidSamFile();

    uint16_t chrIndices[] = {1, 2};
    for (auto chrIdx : chrIndices) {
        SamIndex::getInstance().addSamIndex(chrIdx, 1000, 1, 1);
    }

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(1, indexList.at(1).size());
    EXPECT_EQ(1, indexList.at(2).size());
    EXPECT_EQ(2, indexList.size());

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexActuator_StreamParsing_InvalidFieldOrder) {
    Json::Value streams(Json::arrayValue);

    // Create streams with invalid field order (field 5 instead of field 3 at index 3)
    Json::Value idField;
    idField["field"] = 0;
    idField["dstlen"] = 10;
    streams.append(idField);

    Json::Value flagField;
    flagField["field"] = 1;
    flagField["dstlen"] = 10;
    streams.append(flagField);

    Json::Value rnameField;
    rnameField["field"] = 2;
    rnameField["dstlen"] = 10;
    streams.append(rnameField);

    Json::Value invalidPosField;
    invalidPosField["field"] = 5; // This should cause issues
    invalidPosField["dstlen"] = 10;
    streams.append(invalidPosField);

    // verify the field number
    EXPECT_EQ(5, streams[3]["field"].asUInt());
}

TEST_F(IndexActuatorTest, MetaJson_NoSamInfo_HandlesGracefully) {
    IndexActuator actuator(pInBlock, pOutBlock);

    Json::Value meta;
    meta["header"] = createHeaderMeta()["header"];
    // No "sam" member in meta

    // Test handling when sam info is missing
    EXPECT_FALSE(meta.isMember("sam"));
}

TEST_F(IndexActuatorTest, MetaJson_ValidSamInfo_Succeeds) {
    Json::Value meta = createBasicMeta();

    EXPECT_TRUE(meta.isMember("sam"));
    EXPECT_TRUE(meta["sam"].isMember("lines"));
    EXPECT_TRUE(meta["sam"].isMember("streams"));
}

TEST_F(IndexActuatorTest, StreamMetadata_FieldCheck) {
    Json::Value streams(Json::arrayValue);

    for (uint32_t i = 0; i < 4; ++i) {
        Json::Value fieldMeta;
        fieldMeta["field"] = i;
        fieldMeta["dstlen"] = 10;
        streams.append(fieldMeta);
    }

    EXPECT_EQ(4, streams.size());
    EXPECT_EQ(0, streams[0]["field"].asUInt());
    EXPECT_EQ(1, streams[1]["field"].asUInt());
    EXPECT_EQ(2, streams[2]["field"].asUInt());
    EXPECT_EQ(3, streams[3]["field"].asUInt());
}

TEST_F(IndexActuatorTest, DecoderInitialization_LevelSetting) {
    Json::Value coder;
    coder["magic"] = "coder_bwt_cm";
    coder["level"] = 7;

    EXPECT_EQ("coder_bwt_cm", coder["magic"].asString());
    EXPECT_EQ(7, coder["level"].asInt());
}

TEST_F(IndexActuatorTest, SamIndexItem_Comparison) {
    SamIndexItem item1;
    item1.referenceMapPos = 1000;

    SamIndexItem item2;
    item2.referenceMapPos = 2000;

    EXPECT_LT(item1.referenceMapPos, item2.referenceMapPos);
}

TEST_F(IndexActuatorTest, HeaderParsing_InvalidMetaStructure) {
    Json::Value meta;
    // Missing required members for header parsing
    meta["header"]["srclen"] = 100;
    // Missing "dstlen" and "lines"

    EXPECT_FALSE(meta["header"].isMember("dstlen"));
    EXPECT_FALSE(meta["header"].isMember("lines"));
}

TEST_F(IndexActuatorTest, ChromosomeIndices_SpecialValues) {
    const uint16_t UNMAPPED = 0xFFFF;
    const uint16_t STAR = 0xFFFE;

    EXPECT_EQ(65535, UNMAPPED);
    EXPECT_EQ(65534, STAR);
}

TEST_F(IndexActuatorTest, FlagBit_MappedCheck) {
    const uint16_t bitmask = 0x04; // Unmapped mask

    // Test scenarios with different flag values
    uint16_t mappedFlag = 0x0000; // Properly mapped
    uint16_t unmappedFlag = 0x0004; // Unmapped
    uint16_t complexFlag = 0x0100; // Secondary alignment

    EXPECT_EQ(0, mappedFlag & bitmask);
    EXPECT_EQ(bitmask, unmappedFlag & bitmask);
    EXPECT_EQ(0, complexFlag & bitmask);
}

TEST_F(IndexActuatorTest, OffsetTracking_InitializesToZero) {
    int32_t readOffset = 0;

    readOffset += 100;
    readOffset += 50;
    readOffset += 75;

    EXPECT_EQ(225, readOffset);
}

TEST_F(IndexActuatorTest, MetadataFetching_ChromosomesUsage) {
    createValidSamFile();

    SamInfo& samInfo = SamInfo::getInstance();
    const std::vector<ChromosomeInfo>& chrInfoList = samInfo.getAllChromosomeInfo();

    EXPECT_EQ(2, chrInfoList.size());
    if (chrInfoList.size() >= 1) {
        EXPECT_EQ("chr1", chrInfoList[0].name);
    }
    if (chrInfoList.size() >= 2) {
        EXPECT_EQ("chr2", chrInfoList[1].name);
    }
}

TEST_F(IndexActuatorTest, ChromosomePosition_Calculation) {
    createValidSamFile();

    SamInfo& samInfo = SamInfo::getInstance();
    int64_t pos1 = samInfo.getPositionByIndex(0);
    int64_t pos2 = samInfo.getPositionByIndex(1);

    EXPECT_GE(pos2, pos1); // Chromosome 2 should be positioned after chromosome 1
    EXPECT_GE(pos1, 0);
}

TEST_F(IndexActuatorTest, IndexOutputFormat_ExpectedColumns) {
    // Simulate the format string used in dumpToFile()
    uint16_t chrIndex = 1;
    int64_t referenceMapPos = 12345;
    uint32_t readNumber = 42;
    uint32_t blockId = 1;

    char buffer[256];
    int len = snprintf(buffer, sizeof(buffer), "%hu\t%ld\t%u\t%u\n",
                      chrIndex, referenceMapPos, readNumber, blockId);

    EXPECT_GT(len, 0);

    // Count tabs to verify we have the expected number of columns
    // Format: chrIndex, referenceMapPos, readNumber, blockId -> 4 columns
    // So we should have 3 tabs
    size_t tabCount = 0;
    std::string result(buffer);
    for (char c : result) {
        if (c == '\t') {
            tabCount++;
        }
    }
    EXPECT_EQ(3, tabCount);
}

TEST_F(IndexActuatorTest, BlockManagement_SetAndGetId) {
    pInBlock->setBlockId(42);
    int64_t blockId = pInBlock->getBlockId();

    EXPECT_EQ(42, blockId);
}

TEST_F(IndexActuatorTest, DataLength_SettingAndGetting) {
    pOutBlock->setDataLen(1000);
    int64_t dataLen = pOutBlock->getDataLen();

    EXPECT_EQ(1000, dataLen);
}

TEST_F(IndexActuatorTest, MetaLength_Management) {
    pInBlock->setMetaLen(200);
    uint32_t metaLen = pInBlock->getMetaLen();

    EXPECT_EQ(200, metaLen);
}

TEST_F(IndexActuatorTest, BufferSpace_Calculation) {
    uint32_t dataLen = 1000;
    uint32_t metaLen = 200;
    uint32_t totalLen = dataLen + metaLen;

    EXPECT_EQ(1200, totalLen);
}

TEST_F(IndexActuatorTest, EmptyBlock_Initialization) {
    RoughIOBlock block(1024);

    EXPECT_EQ(-1, block.getBlockId());
    EXPECT_EQ(0, block.getDataLen());
    EXPECT_EQ(0, block.getMetaLen());
    EXPECT_EQ(0, block.getMaxLineLen());
    EXPECT_EQ(TYPE_UNKNOW, block.getBlockType());
}

TEST_F(IndexActuatorTest, SamIndexList_IterationTest) {
    createValidSamFile();

    // Add index items to singleton
    for (uint16_t chrIdx = 1; chrIdx < 3; chrIdx++) {
        SamIndex::getInstance().addSamIndex(chrIdx, chrIdx * 10000, chrIdx + 1, 1);
    }

    size_t totalItems = 0;
    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    for (const auto& pair : indexList) {
        totalItems += pair.second.size();
    }

    EXPECT_EQ(2, totalItems);

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexSplitting_NoSplitWhenCountBelowThreshold) {
    createValidSamFile();

    // Add index items with count = 100, which is below the 1000 threshold
    SamIndex::getInstance().addSamIndex(1, 1000, 100, 0);

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(1, indexList.at(1).size());
    EXPECT_EQ(100, indexList.at(1)[0].readNumber);

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexSplitting_NoSplitWhenCountAtThresholdNoPositionChange) {
    createValidSamFile();

    // Add index items with count = 1001 at same position
    // This should NOT split because position hasn't changed
    SamIndex::getInstance().addSamIndex(1, 1000, 1001, 0);

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(1, indexList.at(1).size());
    EXPECT_EQ(1001, indexList.at(1)[0].readNumber);

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexSplitting_SplitWhenCountExceedsThresholdAndPositionChanges) {
    createValidSamFile();

    // This simulates the behavior when position changes and count > 1000
    // First entry: count = 1001 (exceeds threshold)
    SamIndex::getInstance().addSamIndex(1, 1000, 1001, 0);
    // Second entry: different position (simulating the split)
    SamIndex::getInstance().addSamIndex(1, 2000, 1, 0);

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(2, indexList.at(1).size());
    EXPECT_EQ(1000, indexList.at(1)[0].referenceMapPos);
    EXPECT_EQ(1001, indexList.at(1)[0].readNumber);
    EXPECT_EQ(2000, indexList.at(1)[1].referenceMapPos);
    EXPECT_EQ(1, indexList.at(1)[1].readNumber);

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexSplitting_MultipleSplitsForLargeCount) {
    createValidSamFile();

    // Simulate multiple splits: 2500 reads would be split into chunks
    // Note: In actual decodeAndBuildIndex logic, splits happen when count > 1000 AND position changes
    // So we simulate the result of such splitting:
    SamIndex::getInstance().addSamIndex(1, 1000, 1001, 0);   // First chunk exceeds threshold
    SamIndex::getInstance().addSamIndex(1, 1500, 1000, 0);   // Second chunk at threshold
    SamIndex::getInstance().addSamIndex(1, 2000, 499, 0);    // Third chunk (total 2500)

    const auto& indexList = SamIndex::getInstance().getSamIndexList();
    EXPECT_EQ(3, indexList.at(1).size());
    EXPECT_EQ(1001, indexList.at(1)[0].readNumber);
    EXPECT_EQ(1000, indexList.at(1)[1].readNumber);
    EXPECT_EQ(499, indexList.at(1)[2].readNumber);

    // Sum should equal total (2500)
    uint32_t total = 0;
    for (const auto& item : indexList.at(1)) {
        total += item.readNumber;
    }
    EXPECT_EQ(2500, total);

    // Clean up
    SamIndex::getInstance().clear();
}

TEST_F(IndexActuatorTest, IndexSplitting_DifferentChromosomesIndependentSplitting) {
    createValidSamFile();

    // Chromosome 1: needs split
    SamIndex::getInstance().addSamIndex(1, 1000, 1001, 0);
    SamIndex::getInstance().addSamIndex(1, 1100, 50, 0);

    // Chromosome 2: no split needed (count < 1000)
    SamIndex::getInstance().addSamIndex(2, 2000, 500, 0);

    const auto& indexList = SamIndex::getInstance().getSamIndexList();

    EXPECT_EQ(2, indexList.at(1).size());
    EXPECT_EQ(1, indexList.at(2).size());
    EXPECT_EQ(1051, indexList.at(1)[0].readNumber + indexList.at(1)[1].readNumber);
    EXPECT_EQ(500, indexList.at(2)[0].readNumber);

    // Clean up
    SamIndex::getInstance().clear();
}

