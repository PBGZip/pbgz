
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>

#define private public
#include "fastq_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#undef private

namespace FastqTestData {
    const std::string samllFastQFile = "small.fastq";
    const std::string bigFastQFile = "test.fastq";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};



class FastQActuatorTest : public ::testing::Test {
public:
    // 为测试准备数据对象
	void SetUp() override {
		pInBlock = new RoughIOBlock(FastqTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(FastqTestData::MAX_BLOCK_SIZE);
	}

	// 清除资源
	void TearDown() override {
		if (pInBlock != nullptr) {
			delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }
	}

    void loadFastQData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        IOReader* pIoReader = new FileReader(filename);
        pIoReader->openIO(); 
        BlockReader*  pBlockReader = new BlockReader(pIoReader); 
        pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
};


TEST_F(FastQActuatorTest, testInitialize) {
    FastqActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initEncoder();
    EXPECT_EQ(result, 0);
}

TEST_F(FastQActuatorTest, testPreAnalysisFirstLine) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqActuator actuator(pInBlock, pOutBlock);
    actuator.preAnalysisIdFirstLine(pInBlock->getBuffer(), pInBlock->getNpos()[0] + 1);
    EXPECT_EQ(actuator.idSplitSymbols.size(), 4); // 

    EXPECT_EQ(actuator.idSplitSymbols[0], '.');
    EXPECT_EQ(actuator.idSplitSymbols[1], ' ');
    EXPECT_EQ(actuator.idSplitSymbols[2], '/');
    EXPECT_EQ(actuator.idSplitSymbols[3], '\n');

    
    EXPECT_EQ(actuator.idPositions[0], 12);
    EXPECT_EQ(actuator.idPositions[1], 14);
    EXPECT_EQ(actuator.idPositions[2], 16);
    EXPECT_EQ(actuator.idPositions[3], 18);

    EXPECT_EQ(actuator.idSplitMinLen[0], 12);
    EXPECT_EQ(actuator.idSplitMinLen[1], 1);
    EXPECT_EQ(actuator.idSplitMinLen[2], 1);
    EXPECT_EQ(actuator.idSplitMinLen[3], 1);

    EXPECT_EQ(actuator.idSplitMaxLen[0], 12);    
    EXPECT_EQ(actuator.idSplitMaxLen[1], 1);
    EXPECT_EQ(actuator.idSplitMaxLen[2], 1);
    EXPECT_EQ(actuator.idSplitMaxLen[3], 1);
}

TEST_F(FastQActuatorTest, testPreAnalysisSecondLine) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqActuator actuator(pInBlock, pOutBlock);
    actuator.preAnalysisIdFirstLine(pInBlock->getBuffer(), pInBlock->getNpos()[0] + 1);
    actuator.preAnalysisId(pInBlock->getBuffer() + pInBlock->getNpos()[3] + 1 , pInBlock->getNpos()[4] + 1);

    EXPECT_EQ(actuator.idPositions[4], 12);
    EXPECT_EQ(actuator.idPositions[5], 14);
    EXPECT_EQ(actuator.idPositions[6], 16);
    EXPECT_EQ(actuator.idPositions[7], 18);

    EXPECT_EQ(actuator.idSplitMinLen[0], 12);
    EXPECT_EQ(actuator.idSplitMinLen[1], 1);
    EXPECT_EQ(actuator.idSplitMinLen[2], 1);
    EXPECT_EQ(actuator.idSplitMinLen[3], 1);

    EXPECT_EQ(actuator.idSplitMaxLen[0], 12);    
    EXPECT_EQ(actuator.idSplitMaxLen[1], 1);
    EXPECT_EQ(actuator.idSplitMaxLen[2], 1);
    EXPECT_EQ(actuator.idSplitMaxLen[3], 1);
}

TEST_F(FastQActuatorTest, testPreAnalysisBase) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.preAnalysisBase(pInBlock->getBuffer() + pInBlock->getNpos()[0] + 1 , pInBlock->getNpos()[1] - pInBlock->getNpos()[0]);
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.maxBaseLength, 150);
    EXPECT_EQ(actuator.minBaseLength, 150);
    EXPECT_EQ(actuator.baseNCount, 0);
}

TEST_F(FastQActuatorTest, testPreAnalysisComment) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.preAnalysisComment(pInBlock->getBuffer() + pInBlock->getNpos()[1] + 1,
                                                    pInBlock->getNpos()[2] - pInBlock->getNpos()[1], 2);
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.commentType, CommentType::PLUS_ONLY);
}

TEST_F(FastQActuatorTest, testPreAnalysis) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.preAnalysis();

    EXPECT_EQ(result, 0);

    EXPECT_EQ(actuator.idSplitSymbols[0], '.');
    EXPECT_EQ(actuator.idSplitSymbols[1], ' ');
    EXPECT_EQ(actuator.idSplitSymbols[2], '/');
    EXPECT_EQ(actuator.idSplitSymbols[3], '\n');

    
    EXPECT_EQ(actuator.idPositions[0], 12);
    EXPECT_EQ(actuator.idPositions[1], 14);
    EXPECT_EQ(actuator.idPositions[2], 16);
    EXPECT_EQ(actuator.idPositions[3], 18);

    EXPECT_EQ(actuator.idSplitMinLen[0], 12);
    EXPECT_EQ(actuator.idSplitMinLen[1], 1);
    EXPECT_EQ(actuator.idSplitMinLen[2], 1);
    EXPECT_EQ(actuator.idSplitMinLen[3], 1);

    
    EXPECT_EQ(actuator.idSplitMaxLen[0], 12);    
    EXPECT_EQ(actuator.idSplitMaxLen[1], 5);
    EXPECT_EQ(actuator.idSplitMaxLen[2], 5);
    EXPECT_EQ(actuator.idSplitMaxLen[3], 1);

    EXPECT_EQ(actuator.idPositions[4], 12);
    EXPECT_EQ(actuator.idPositions[5], 14);
    EXPECT_EQ(actuator.idPositions[6], 16);
    EXPECT_EQ(actuator.idPositions[7], 18);

    EXPECT_EQ(actuator.maxBaseLength, 150);
    EXPECT_EQ(actuator.minBaseLength, 150);
    EXPECT_EQ(actuator.baseNCount, 73);

    EXPECT_EQ(actuator.commentType, CommentType::PLUS_ONLY);

    EXPECT_EQ(actuator.qualityFreqTable[0].first, 'F' - '!');
    EXPECT_EQ(actuator.qualityFreqTable[1].first, ':' - '!');
    EXPECT_EQ(actuator.qualityFreqTable[2].first, ',' - '!');
}

