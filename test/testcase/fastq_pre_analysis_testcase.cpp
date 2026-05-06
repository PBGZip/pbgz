/*
 * fastq_pre_analysis_testcase.cpp - Test cases for FASTQ pre-analysis
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

namespace FastqPreAnalysisTestData {
    const std::string testFastQFile = "test.fq";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class FastqPreAnalysisTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
		pInBlock = new RoughIOBlock(FastqPreAnalysisTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(FastqPreAnalysisTestData::MAX_BLOCK_SIZE);

        generateTestFastqFile(FastqPreAnalysisTestData::testFastQFile);
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

        std::remove(FastqPreAnalysisTestData::testFastQFile.c_str());
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

    void generateTestFastqFile(const std::string& filename) {
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
        
        // Create test file based on FASTQ data in current file
        file << "@chr1_41586121_41586549_15:5:0_10:3:2_0/1\n";
        file << "TGTCCAGTCAGATCCCTGGGTATTCACTGGAAATATGGGGGTCATTTCCCTCACTTGGAGGAGCAAAGGTGTACTGCATATAGCCCCACGTGTGTGCCTCCTCCTCGTGCCCGCTTCAAGCAGGGCAGTCTCAGCTTGCAGGGCTCAGGT\n";
        file << "+\n";
        file << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_152370116_152370409_20:3:3_17:3:3_1/1\n";
        file << "TATAATGGCTGCATACAACTCATCCTTTTTTATGGCTGGATAGTATTCCATGGTGACTATGTGTCACAATTTCATAATCCAGTCTGATCATTTTTGGACATTTGGGTTGTGTTACAACGAGCCAATCCGCGAATAGACCCCAATTGTCCT\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_17297153_17297533_17:5:3_19:4:0_2/1\n";
        file << "TCAGTAAAACTTGGGTCAAGTCCTTATCTATTTTACTGAGCGATTTGTGATCATGGATGGGACATTCCAGTTCTCTGATCCTATATCGTTGCCCTACAAAATGGGGTCAATAATGCCTTCCTCATGGGGTTTAGATAAGAATTTATGAGC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_55329955_55330247_22:2:3_15:2:1_3/1\n";
        file << "TTTTCCTGCAAGCTCCCGTCTTGCATTGTCAACAACTCATCAATACATAACTTTGTTTTATGTATTTCTGATTCAACACACCTTATTAAACAAATACTGTTGATTTATTAAACATAGAACTCACATCCAATAGGCGTGCAGCATATACAC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_11312507_11312801_22:3:1_15:2:2_4/1\n";
        file << "GACTCTATTAGATGACCGACCGCTGCCAGCTCTGGCGGTGAGTAGGGCTGAGCCGCTCACCGCAGCCAGGCGGCCCTGAGAACGCAGCATAGTGATGACAAAAAGACCTGAAGGTCGTGGGCGATGTTTTTCGGTGAGTTGGGAGGTGCC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_205902928_205903303_11:1:1_18:6:0_5/1\n";
        file << "AAGCAGGAAGAAGTCGCATCACTGACTAGACCAATAACAGGTTCTGAAATTGATGCAATAATTGCCTACCAACCAACAAAAGTCCAGGACCAGATGGATTCCCAGGTTAATTGTACGAGAGGTACAAAGAGGAGCTGGTACCCTTCCTTC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_67278157_67278448_15:1:2_17:4:1_6/1\n";
        file << "TACTAGTCTTCCACTGGGACCCGTTTGATTCATAAAGAAGTCCTCAAAAAAGATCTTCAAAGAAGTCACATGAAAATGGGTCCCTTCCAACAAAAAAATCCGTGACTACCTCATGTGTGTTACAGAATGTGCAGCCAAACTGAAAAGGAC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_225077116_225077417_15:4:1_9:3:0_7/1\n";
        file << "TAGTTTATATTATCAATGAATTATCTCAAATGGTTATTCAATACATCTATTGCGAAAATCTTGCAGTTGTTTTCCTTTAATCTTTTACTGAGTGGAGATATCTTAATAGACTTACTAAGGATAATCTCACATCCACATTGATCCGATTTA\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_221304894_221305124_17:4:1_20:2:0_8/1\n";
        file << "ACACAAGTAGTCTCAGGTCGGAGCATATATTTATATTCAAGCTGGGAAAACCATGCTATGAATCCTTTAAAAACTAAAAATGCGTATACATAGATGTTCATTTAAAAATGGAAGAATGGTTTATCAGTGCCAATCCATCATCTTGTCTTT\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_172301339_172301575_15:1:1_17:2:1_9/1\n";
        file << "TATTAATTTCAGAAGAAGTATTCCTCTAAATCGAGTTACAGGAACAAAAGGTATTTTGTTTTTGATTATTGATAGAGACATGATGTCTCCTGATGCCTGAGCTGGAGTTCAGTGGCACAGTGATAGCTCACTTTAGGCTCCAACTCCTGG\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_158527486_158527796_13:8:3_16:3:1_a/1\n";
        file << "CTACTCCAGATACTACAGTTGCACCACTTCGAGAAGACACGGAAAGGAAAACAATAAATCAGGAAAGAAAAGAAAATCCTATCCACCCAAAAATAATTACAAAAATTAGAAGAGCCATGCACTCCAGAAGACCCGGAGCCGTATAAGTAT\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_157144833_157145062_18:4:1_15:6:2_b/1\n";
        file << "TACATGGAGTAAATGAGAAGCAGGGGGTGCAACAAACTCAGAGTAACTCCAAGTACCTGGTACCCAGCAGCAGTGTACTTGGGCACCGTCGGACAACCTGTGTTAACTGATAGATGATTCTTCTCTGGGCTAAACTCAGACCCAAGGGTC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file << "@chr1_107966086_107966417_14:2:1_12:6:0_c/1\n";
        file << "GCGATTCTCCAACCTCAGCCTCCCGATCCGCTGGGACTACAGGTGCGTGCCACCACGGCCGGCCAATCTTTTGTAATATTAGTAGAGTCGGGATTTCACCGTCTTAGCCCGGATTCTCTCGATCTCCTGAGCTCGAGATCCGCCCGCGTC\n";
        file << "+\n";
        file << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        
        file.close();
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
    PbgzParameter para;
};

// Test preAnalysisIdFirstLine method
TEST_F(FastqPreAnalysisTest, testPreAnalysisIdFirstLine) {
    loadFastQData(FastqPreAnalysisTestData::testFastQFile);
    PbgzParameter para;
    FastqCodecActuator actuator(pInBlock, pOutBlock, para);
    
    // Test analysis of first line ID
    int32_t result = actuator.preAnalysisIdFirstLine(pInBlock->getBuffer(), pInBlock->getNpos()[0] + 1);
    EXPECT_EQ(result, 0);
    
    // Verify separator recognition
    EXPECT_EQ(actuator.idSplitSymbols.size(), 11); //
    EXPECT_EQ(actuator.idSplitSymbols[0], '_');  // First underscore
    EXPECT_EQ(actuator.idSplitSymbols[1], '_');  // Second underscore
    EXPECT_EQ(actuator.idSplitSymbols[2], '_');  // Third underscore
    EXPECT_EQ(actuator.idSplitSymbols[3], ':');  // 
    EXPECT_EQ(actuator.idSplitSymbols[4], ':'); // 
    EXPECT_EQ(actuator.idSplitSymbols[5], '_'); // 
    EXPECT_EQ(actuator.idSplitSymbols[6], ':'); // 
    EXPECT_EQ(actuator.idSplitSymbols[7], ':'); // 
    EXPECT_EQ(actuator.idSplitSymbols[8], '_'); // 
    EXPECT_EQ(actuator.idSplitSymbols[9], '/'); // 
    EXPECT_EQ(actuator.idSplitSymbols[10], '\n'); // 
    
    // "@chr1_41586121_41586549_15:5:0_10:3:2_0/1\n";
    // Verify position information
    EXPECT_EQ(actuator.idPositions[0], 5);  // First underscore position
    EXPECT_EQ(actuator.idPositions[1], 14);  // 
    EXPECT_EQ(actuator.idPositions[2], 23); // 
    EXPECT_EQ(actuator.idPositions[3], 26); //
    EXPECT_EQ(actuator.idPositions[4], 28); // 
    EXPECT_EQ(actuator.idPositions[5], 30); // 
    EXPECT_EQ(actuator.idPositions[6], 33); // 
    EXPECT_EQ(actuator.idPositions[7], 35); // 
    EXPECT_EQ(actuator.idPositions[8], 37); // 
    EXPECT_EQ(actuator.idPositions[9], 39); // 
    EXPECT_EQ(actuator.idPositions[10], 41); // 
}

// Test preAnalysisId method - handle multiple ID records
TEST_F(FastqPreAnalysisTest, testPreAnalysisIdMultipleRecords) {
    loadFastQData(FastqPreAnalysisTestData::testFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, para);
    
    // Analyze first line first
    actuator.preAnalysisIdFirstLine(pInBlock->getBuffer(), pInBlock->getNpos()[0] + 1);
    
    // Analyze second line ID (fourth line)
    uint8_t* secondIdStart = pInBlock->getBuffer() + pInBlock->getNpos()[3] + 1;
    uint32_t secondIdEnd = pInBlock->getNpos()[4] + 1;
    int32_t result = actuator.preAnalysisId(secondIdStart, secondIdEnd);
    EXPECT_EQ(result, 0);
    
    // Verify second ID record is also correctly analyzed
    // "@chr1_152370116_152370409_20:3:3_17:3:3_1/1\n";
    EXPECT_EQ(actuator.idPositions[11], 5);  // First underscore position
    EXPECT_EQ(actuator.idPositions[12], 15); // Second underscore position  
    EXPECT_EQ(actuator.idPositions[13], 25); // Third underscore position
    EXPECT_EQ(actuator.idPositions[14], 28); // First colon position
    EXPECT_EQ(actuator.idPositions[15], 30); // Second colon position
    EXPECT_EQ(actuator.idPositions[16], 32); // Fourth underscore position
    EXPECT_EQ(actuator.idPositions[17], 35); // Third colon position
    EXPECT_EQ(actuator.idPositions[18], 37); // Fourth colon position
    EXPECT_EQ(actuator.idPositions[19], 39); // Fifth underscore position
    EXPECT_EQ(actuator.idPositions[20], 41); // Slash position
    EXPECT_EQ(actuator.idPositions[21], 43); // Newline position
}

// Test ID format min and max length
TEST_F(FastqPreAnalysisTest, testPreAnalysisIdLengthAnalysis) {
    loadFastQData(FastqPreAnalysisTestData::testFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, para);
    
    // Analyze all records
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Verify ID field length analysis
    EXPECT_EQ(actuator.idSplitMinLen.size(), 11);
    EXPECT_EQ(actuator.idSplitMaxLen.size(), 11);
    
    // Verify chromosome field length (part before first underscore)
    EXPECT_EQ(actuator.idSplitMinLen[0], 5); // "chr1"
    EXPECT_EQ(actuator.idSplitMaxLen[0], 5); // "chr1"
    
    // Verify position field length analysis
    EXPECT_EQ(actuator.idSplitMinLen[1], 8); // Position field
    EXPECT_EQ(actuator.idSplitMaxLen[1], 9); // Position field
}

// Test complete preAnalysis method
TEST_F(FastqPreAnalysisTest, testPreAnalysisComplete) {
    loadFastQData(FastqPreAnalysisTestData::testFastQFile);
    FastqCodecActuator actuator(pInBlock, pOutBlock, para);
    
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Verify ID analysis results
    EXPECT_EQ(actuator.idSplitSymbols[0], '_');
    EXPECT_EQ(actuator.idSplitSymbols[1], '_');
    EXPECT_EQ(actuator.idSplitSymbols[2], '_');
    EXPECT_EQ(actuator.idSplitSymbols[3], ':');
    EXPECT_EQ(actuator.idSplitSymbols[4], ':');
    EXPECT_EQ(actuator.idSplitSymbols[5], '_');
    EXPECT_EQ(actuator.idSplitSymbols[6], ':');
    EXPECT_EQ(actuator.idSplitSymbols[7], ':');
    EXPECT_EQ(actuator.idSplitSymbols[8], '_');
    EXPECT_EQ(actuator.idSplitSymbols[9], '/');
    EXPECT_EQ(actuator.idSplitSymbols[10], '\n');
    
    // Verify Base sequence analysis
    EXPECT_EQ(actuator.maxBaseLength, 150); // Based on actual test data length
    EXPECT_EQ(actuator.minBaseLength, 150); // All sequences same length
    EXPECT_EQ(actuator.baseNCount, 0);      // Test data should have no N bases
    
    // Verify comment line analysis
    EXPECT_EQ(actuator.commentType, CommentType::PLUS_ONLY);
    
    // Verify quality score analysis
    EXPECT_GT(actuator.qualityFreqTable.size(), 0);
}