/*
 * sam_reader_testcase.cpp - Test cases for SAM reader functionality
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
#include <stdio.h>
#include <stdlib.h>

#define private public
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include <actg.h>
#undef private

namespace SamReaderTestData {
    const std::string smallSamFile = "small_test.sam";
    const std::string largeSamFile = "large_test.sam";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamReaderTest : public ::testing::Test {
public:
    void SetUp() override {
        // Generate small SAM file
        generateSmallSamFile(SamReaderTestData::smallSamFile);
        // Generate large SAM file
        generateLargeSamFile(SamReaderTestData::largeSamFile);
    }

    void TearDown() override {
        std::remove(SamReaderTestData::smallSamFile.c_str());
        std::remove(SamReaderTestData::largeSamFile.c_str());
    }

    void generateSmallSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }
        
        if (!file.is_open()) {
            return;
        }
        
        // Write SAM header
        file << "@HD    VN:1.6  SO:coordinate\n";
        file << "@SQ	SN:chr1	LN:248956422\n";
        file << "@SQ	SN:chr2	LN:242193529\n";
        file << "@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem hg38.p14.fa ERR016060_1.fastq\n";
        
        // Write several SAM records
        file << "ERR016060.1	4	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	AS:i:0	XS:i:0\n";
        file << "ERR016060.2	4	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	AS:i:0	XS:i:0\n";
        file << "ERR016060.3	4	*	0	0	*	*	0	0	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	AS:i:0	XS:i:0\n";
        
        file.close();
    }

    void generateLargeSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }
        
        if (!file.is_open()) {
            return;
        }
        
        // Write SAM header
        file << "@HD    VN:1.6  SO:coordinate\n";
        file << "@SQ	SN:chr1	LN:248956422\n";
        file << "@SQ	SN:chr2	LN:242193529\n";
        file << "@SQ	SN:chr3	LN:198295559\n";
        file << "@SQ	SN:chr4	LN:190214555\n";
        file << "@SQ	SN:chr5	LN:181538259\n";
        file << "@SQ	SN:chr6	LN:170805979\n";
        file << "@SQ	SN:chr7	LN:159345973\n";
        file << "@SQ	SN:chr8	LN:145138636\n";
        file << "@SQ	SN:chr9	LN:138394717\n";
        file << "@SQ	SN:chr10	LN:133797422\n";
        file << "@SQ	SN:chr11	LN:135086622\n";
        file << "@SQ	SN:chr12	LN:133275309\n";
        file << "@SQ	SN:chr13	LN:114364328\n";
        file << "@SQ	SN:chr14	LN:107043718\n";
        file << "@SQ	SN:chr15	LN:101991189\n";
        file << "@SQ	SN:chr16	LN:90338345\n";
        file << "@SQ	SN:chr17	LN:83257441\n";
        file << "@SQ	SN:chr18	LN:80373285\n";
        file << "@SQ	SN:chr19	LN:58617616\n";
        file << "@SQ	SN:chr20	LN:64444167\n";
        file << "@SQ	SN:chr21	LN:46709983\n";
        file << "@SQ	SN:chr22	LN:50818468\n";
        file << "@SQ	SN:chrX	LN:156040895\n";
        file << "@SQ	SN:chrY	LN:57227415\n";
        file << "@SQ	SN:chrM	LN:16569\n";
        file << "@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem hg38.p14.fa ERR016060_1.fastq\n";
        
        // Write many SAM records (100 records)
        for (int i = 1; i <= 100; i++) {
            file << "ERR016060." << i << "\t4\t*\t0\t0\t*\t*\t0\t0\t";
            file << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\t";
            file << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\t";
            file << "AS:i:0\tXS:i:0\n";
        }
        
        file.close();
    }

    bool isSamFile(const std::string& filename) {
        RoughIOBlock* pBlock = new RoughIOBlock(SamReaderTestData::MAX_BLOCK_SIZE);
        
        IOReader* pIoReader = new FileReader(filename);
        pIoReader->openIO();
        BlockReader* pBlockReader = new BlockReader(pIoReader);
        
        pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
        
        bool isSam = (pBlock->getBlockType() == SAM);
        
        delete pBlockReader;
        delete pIoReader;
        delete pBlock;
        
        return isSam;
    }
};

TEST_F(SamReaderTest, TestSmallSamFile) {
    // Test small SAM file
    bool isSam = isSamFile(SamReaderTestData::smallSamFile);
    EXPECT_TRUE(isSam) << "Small SAM file should be recognized as SAM format";
}

TEST_F(SamReaderTest, TestLargeSamFile) {
    // Test large SAM file
    bool isSam = isSamFile(SamReaderTestData::largeSamFile);
    EXPECT_TRUE(isSam) << "Large SAM file should be recognized as SAM format";
}

TEST_F(SamReaderTest, TestSamFileReading) {
    // Test reading SAM file and verify block type
    RoughIOBlock* pBlock = new RoughIOBlock(SamReaderTestData::MAX_BLOCK_SIZE);
    
    IOReader* pIoReader = new FileReader(SamReaderTestData::smallSamFile);
    pIoReader->openIO();
    BlockReader* pBlockReader = new BlockReader(pIoReader);
    
    int32_t result = pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
    
    // Verify successful reading
    EXPECT_GT(result, 0) << "Reading SAM file should succeed";
    
    // Verify block type is SAM
    EXPECT_EQ(pBlock->getBlockType(), SAM) << "Read block type should be SAM";
    
    // Verify block has data
    EXPECT_GT(pBlock->getDataLen(), 0) << "Block should have data";
    
    delete pBlockReader;
    delete pIoReader;
    delete pBlock;
}

TEST_F(SamReaderTest, TestLargeSamFileReading) {
    // Test reading large SAM file and verify block type
    RoughIOBlock* pBlock = new RoughIOBlock(SamReaderTestData::MAX_BLOCK_SIZE);
    
    IOReader* pIoReader = new FileReader(SamReaderTestData::largeSamFile);
    pIoReader->openIO();
    BlockReader* pBlockReader = new BlockReader(pIoReader);
    
    int32_t result = pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
    
    // Verify successful reading
    EXPECT_GT(result, 0) << "Reading large SAM file should succeed";
    
    // Verify block type is SAM
    EXPECT_EQ(pBlock->getBlockType(), SAM) << "Read block type should be SAM";
    
    // Verify block has data
    EXPECT_GT(pBlock->getDataLen(), 0) << "Block should have data";
    
    // Large file should have more data
    EXPECT_GT(pBlock->getDataLen(), 1000) << "Large SAM file should contain more data";
    
    delete pBlockReader;
    delete pIoReader;
    delete pBlock;
}
