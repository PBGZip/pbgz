/*
 * fastq_actuator_testcase.cpp - Test cases for FASTQ actuator
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

#define private public
#include "fastq_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include <random>
#undef private

namespace FastqTestData {
    const std::string samllFastQFile = "small.fastq";
    const std::string bigFastQFile = "test.fastq";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class FastQActuatorTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
		pInBlock = new RoughIOBlock(FastqTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(FastqTestData::MAX_BLOCK_SIZE);

        generateFastqFile(FastqTestData::samllFastQFile, 25);
        generateFastqFile(FastqTestData::bigFastQFile, 250000);
	}

	// Clean up resources
	void TearDown() override {
		if (pInBlock != nullptr) {
			delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }

        std::remove(FastqTestData::samllFastQFile.c_str());
        std::remove(FastqTestData::bigFastQFile.c_str());        
	}

    void loadFastQData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        IOReader* pIoReader = new FileReader(filename);
        pIoReader->openIO(); 
        BlockReader*  pBlockReader = new BlockReader(pIoReader); 
        pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);
        
        // Clean up allocated resources
        delete pBlockReader;
        delete pIoReader;
    }

    void generateFastqFile(const std::string& filename, int numRecords) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            // If unable to open in current directory, try creating in test directory
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // If still fails, try creating in build directory
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }
        
        if (!file.is_open()) {
            return; // Unable to create file
        }
        
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for (int i = 1; i <= numRecords; ++i) {
            // Header line: @SRR12922210.i 1/1
            file << "@SRR12922210." << i << " " << i << "/1\n";
            
            // Sequence line: random DNA sequence
            int seqLen = 150;
            for (int j = 0; j < seqLen; ++j) {
                char base = "ACGTN"[gen() % 5];
                file << base;
            }
            file << "\n";
            
            // Plus line
            file << "+\n";
            
            // Quality line: weighted random quality scores
            // F (highest probability), : (second probability), , (lowest probability)
            std::discrete_distribution<int> quality_dist({70, 20, 10}); // F:70%, :20%, ,10%
            for (int j = 0; j < seqLen; ++j) {
                int quality_choice = quality_dist(gen);
                char quality_char;
                switch (quality_choice) {
                    case 0: quality_char = 'F'; break;  // F - highest probability
                    case 1: quality_char = ':'; break;  // : - second probability  
                    case 2: quality_char = ','; break;  // , - lowest probability
                    default: quality_char = 'F'; break; // fallback
                }
                file << quality_char;
            }
            file << "\n";
        }
        
        file.close();
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
    PbgzParameter para;
};


TEST_F(FastQActuatorTest, testInitialize) {
    FastqCodecActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initEncoder();
    EXPECT_EQ(result, 0);
}

TEST_F(FastQActuatorTest, testPreAnalysisFirstLine) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, nullptr);
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
    FastqCodecActuator actuator(pInBlock, pOutBlock, nullptr);
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
    FastqCodecActuator actuator(pInBlock, pOutBlock, nullptr);
    int32_t result = actuator.preAnalysisBase(pInBlock->getBuffer() + pInBlock->getNpos()[0] + 1 , pInBlock->getNpos()[1] - pInBlock->getNpos()[0]);
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.maxBaseLength, 150);
    EXPECT_EQ(actuator.minBaseLength, 150);
    EXPECT_GE(actuator.baseNCount, 0);
}

TEST_F(FastQActuatorTest, testPreAnalysisComment) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, nullptr);
    int32_t result = actuator.preAnalysisComment(pInBlock->getBuffer() + pInBlock->getNpos()[1] + 1,
                                                    pInBlock->getNpos()[2] - pInBlock->getNpos()[1], 2);
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.commentType, CommentType::PLUS_ONLY);
}

TEST_F(FastQActuatorTest, testPreAnalysis) {
    loadFastQData(FastqTestData::bigFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, nullptr);
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
    EXPECT_GT(actuator.baseNCount, 0);

    EXPECT_EQ(actuator.commentType, CommentType::PLUS_ONLY);

    EXPECT_EQ(actuator.qualityFreqTable[0].first, 'F' - '!');
    EXPECT_EQ(actuator.qualityFreqTable[1].first, ':' - '!');
    EXPECT_EQ(actuator.qualityFreqTable[2].first, ',' - '!');
}
