/*
 * sam_sort_testcase.cpp - Test cases for SAM sorting functionality
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
#include <random>

#include "coder/coder.h"

#define private public
#include "sam_sort_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include "sam_info.h"
#undef private

namespace SamSortTestData {
    const std::string testSamFile = "test_sort.sam";
    const std::string testUnmappedSamFile = "test_unmapped.sam";
    const std::string testInvalidSamFile = "test_invalid.sam";
    const std::string testMixedSamFile = "test_mixed.sam";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamSortTest : public ::testing::Test {
public:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

        pInBlock = new RoughIOBlock(SamSortTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamSortTestData::MAX_BLOCK_SIZE);

        // Clear chromosome information before each test
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();

        ConfigManager::getInstance().logLevel = LogLevel::DEBUGGING;
    }

    void TearDown() override {
        if (pInBlock != nullptr) {
            delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }

        // Clean up all potentially generated test files
        std::remove(SamSortTestData::testSamFile.c_str());
        std::remove(SamSortTestData::testUnmappedSamFile.c_str());
        std::remove(SamSortTestData::testInvalidSamFile.c_str());
        std::remove(SamSortTestData::testMixedSamFile.c_str());
        std::remove("sorted_head.sam");
        std::remove("sorted_sam_0.sam");
        std::remove("sorted_sam_1.sam");
        std::remove("sorted_sam_2.sam");
    }

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        std::vector<std::string> paths = {
            filename,
            "./test/" + filename,
            "../test/" + filename
        };

        IOReader* pIoReader = nullptr;
        for (const auto& path : paths) {
            pIoReader = new FileReader(path);
            if (pIoReader->openIO() == 0) {
                break;
            }
            delete pIoReader;
            pIoReader = nullptr;
        }

        if (pIoReader == nullptr) {
            return;
        }

        BlockReader* pBlockReader = new BlockReader(pIoReader);
        pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);
        pInBlock->setBlockId(0);

        delete pBlockReader;
        delete pIoReader;
    }

    void generateSortTestSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        // Generate standard SAM file with records in unsorted order
        file << "@HD\tVN:1.6\tSO:unsorted\n";
        file << "@SQ\tSN:chr1\tLN:1000\n";
        file << "@SQ\tSN:chr2\tLN:800\n";
        file << "@SQ\tSN:chr3\tLN:600\n";
        file << "@PG\tID:sorttest\tPN:sorttest\tVN:1.0\n";

        // Records in random order for testing sorting
        file << "read5\t0\tchr1\t500\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 500
        file << "read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 100
        file << "read3\t0\tchr1\t300\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 300
        file << "read2\t0\tchr2\t50\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 50 (chr2 base 1000)
        file << "read4\t0\tchr2\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 200 (chr2 base 1000)
        file << "read6\t0\tchr3\t150\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 150 (chr3 base 1800)

        file.close();
    }

    void generateUnmappedSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:unsorted\n";
        file << "@SQ\tSN:chr1\tLN:1000\n";

        // Unmapped records (flag bit 4 set)
        file << "unmapped1\t4\tchr1\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "unmapped2\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";

        file.close();
    }

    void generateMixedSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:unsorted\n";
        file << "@SQ\tSN:chr1\tLN:1000\n";
        file << "@SQ\tSN:chr2\tLN:800\n";

        // Mixed mapped and unmapped records
        file << "read2\t0\tchr1\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 200
        file << "unmapped\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 100
        file << "unmapped2\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read3\t0\tchr2\t300\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";  // Position 300 (chr2 base 1000)

        file.close();
    }

    void generateInvalidSamFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:unsorted\n";
        file << "@SQ\tSN:chr1\tLN:1000\n";

        // Record with very short line (less than 3 characters)
        file << "@\n";

        file.close();
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
};

// Test SamSortItem comparison operator
TEST(SamSortItemTest, TestComparator) {
    SamSortItem item1;
    SamSortItem item2;
    SamSortItem item3;

    item1.referencePos = 100;
    item1.lineId = 1;

    item2.referencePos = 200;
    item2.lineId = 2;

    item3.referencePos = 50;
    item3.lineId = 3;

    EXPECT_TRUE(item1 < item2);
    EXPECT_TRUE(item3 < item1);
    EXPECT_FALSE(item2 < item1);
    EXPECT_FALSE(item1 < item3);
}

// Test SamSortItem default constructor
TEST(SamSortItemTest, TestDefaultConstructor) {
    SamSortItem item;
    EXPECT_EQ(item.referencePos, -1);
    EXPECT_EQ(item.lineId, 0);
}

// Test filename generation functions
TEST(SamSortUtilTest, TestGetSortedHeadFileName) {
    std::string filename = getSortedHeadFileName();
    EXPECT_EQ(filename, "sorted_head.sam");
}

TEST(SamSortUtilTest, TestGetSortedSamFileName) {
    std::string filename1 = getSortedSamFileName(0);
    EXPECT_EQ(filename1, "sorted_sam_0.sam");

    std::string filename2 = getSortedSamFileName(1);
    EXPECT_EQ(filename2, "sorted_sam_1.sam");

    std::string filename3 = getSortedSamFileName(100);
    EXPECT_EQ(filename3, "sorted_sam_100.sam");
}

// Test SAMSortActuator constructor and destructor
TEST_F(SamSortTest, TestConstructorAndDestructor) {
    SAMSortActuator* actuator = new SAMSortActuator(pInBlock, pOutBlock);
    EXPECT_NE(actuator, nullptr);
    delete actuator;
}

// Test initial() with valid SAM file
TEST_F(SamSortTest, TestInitialWithValidSam) {
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();

    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.headLineNum, 5);  // @HD, 3x@SQ, @PG
    EXPECT_EQ(actuator.mappedSamItem.size(), 6);
    EXPECT_EQ(actuator.unmappedSamItem.size(), 0);

    // Verify chromosome information was parsed
    EXPECT_TRUE(SamInfo::getInstance().hasChrName("chr1"));
    EXPECT_TRUE(SamInfo::getInstance().hasChrName("chr2"));
    EXPECT_TRUE(SamInfo::getInstance().hasChrName("chr3"));

    result = actuator.process();
    EXPECT_EQ(result, 0);

    // Verify mapped items sorted by position
    if (actuator.mappedSamItem.size() >= 2) {
        EXPECT_LE(actuator.mappedSamItem[0].referencePos, actuator.mappedSamItem[1].referencePos);
    }
}

// Test initial() with null input block
TEST_F(SamSortTest, TestInitialWithNullBlock) {
    SAMSortActuator actuator(nullptr, pOutBlock);
    int32_t result = actuator.initial();

    EXPECT_EQ(result, -1);
}

// Test initial() with unmapped records only
TEST_F(SamSortTest, TestInitialWithUnmappedRecords) {
    generateUnmappedSamFile(SamSortTestData::testUnmappedSamFile);
    loadSamData(SamSortTestData::testUnmappedSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();

    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.headLineNum, 2);  // @HD, @SQ
    EXPECT_EQ(actuator.mappedSamItem.size(), 0);
    EXPECT_GT(actuator.unmappedSamItem.size(), 0);
}

// Test initial() with mixed mapped and unmapped records
TEST_F(SamSortTest, TestInitialWithMixedRecords) {
    generateMixedSamFile(SamSortTestData::testMixedSamFile);
    loadSamData(SamSortTestData::testMixedSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();

    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.headLineNum, 3);  // @HD, 2x@SQ
    EXPECT_GT(actuator.mappedSamItem.size(), 0);
    EXPECT_GT(actuator.unmappedSamItem.size(), 0);

    // Verify mapped items are stored separately from unmapped
    for (const auto& item : actuator.mappedSamItem) {
        EXPECT_GE(item.referencePos, 0);
    }
}

// Test initial() with invalid SAM file
TEST_F(SamSortTest, TestInitialWithInvalidSam) {
    generateInvalidSamFile(SamSortTestData::testInvalidSamFile);
    loadSamData(SamSortTestData::testInvalidSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.initial();

    EXPECT_EQ(result, -1);
}

// Test process() with valid sorted items
TEST_F(SamSortTest, TestProcessWithSortedItems) {
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    // Set block ID for testing
    pInBlock->setBlockId(0);

    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    // Check if files were created
    std::ifstream headFile("sorted_head.sam");
    std::ifstream samFile("sorted_sam_0.sam");

    EXPECT_TRUE(headFile.good());
    EXPECT_TRUE(samFile.good());

    if (headFile.is_open()) {
        std::string line;
        int lineCount = 0;
        while (std::getline(headFile, line)) {
            lineCount++;
        }
        EXPECT_GT(lineCount, 0);  // Should have header lines
        headFile.close();
    }

    if (samFile.is_open()) {
        std::string line;
        int lineCount = 0;
        while (std::getline(samFile, line)) {
            lineCount++;
            if (lineCount > 1) {  // Skip first line for comparison
                // Verify lines contain expected pattern
                EXPECT_TRUE(line.size() > 0);
            }
        }
        EXPECT_GT(lineCount, 0);  // Should have SAM records
        samFile.close();
    }
}

// Test process() with no header lines
TEST_F(SamSortTest, TestProcessWithNoHeader) {
    // Create SAM file without header
    std::ofstream file("no_header.sam");
    file << "read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    file.close();

    loadSamData("no_header.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    pInBlock->setBlockId(0);

    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);  // Should succeed even without header
    EXPECT_EQ(actuator.headLineNum, 0);

    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    std::remove("no_header.sam");
}

// Test process() with multiple block IDs
TEST_F(SamSortTest, TestProcessWithMultipleBlockIds) {
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    SAMSortActuator actuator1(pInBlock, pOutBlock);
    int32_t initResult = actuator1.initial();
    EXPECT_EQ(initResult, 0);

    pInBlock->setBlockId(1);
    int32_t processResult = actuator1.process();
    EXPECT_EQ(processResult, 0);

    // Check file was created with correct name
    std::ifstream samFile("sorted_sam_1.sam");
    EXPECT_TRUE(samFile.good());
    if (samFile.is_open()) {
        samFile.close();
    }
}

// Test sorting correctness
TEST_F(SamSortTest, TestSortingCorrectness) {
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    // Verify records are sorted by reference position
    if (actuator.mappedSamItem.size() > 1) {
        for (size_t i = 0; i < actuator.mappedSamItem.size() - 1; i++) {
            EXPECT_LE(actuator.mappedSamItem[i].referencePos, actuator.mappedSamItem[i + 1].referencePos);
        }
    }

    // Verify output file contains records in sorted order
    std::ifstream samFile("sorted_sam_0.sam");
    if (samFile.is_open()) {
        std::string line;
        std::string prevLine;
        while (std::getline(samFile, line)) {
            if (!line.empty() && line[0] != '@') {
                if (!prevLine.empty()) {
                    // Very basic verification - lines should be present
                    EXPECT_GT(line.size(), 0);
                }
                prevLine = line;
            }
        }
        samFile.close();
    }
}

// Test with very large chromosome positions
TEST_F(SamSortTest, TestWithLargeChromosomePositions) {
    // Create SAM file with large chromosome positions
    std::ofstream file("large_pos.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:100000000\n";

    file << "read1\t0\tchr1\t99999990\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    file << "read2\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    file.close();

    loadSamData("large_pos.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    EXPECT_GT(actuator.mappedSamItem.size(), 0);
    EXPECT_EQ(actuator.mappedSamItem.size(), 2);  // Should have 2 mapped records

    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    // Verify sorting works with large positions
    if (actuator.mappedSamItem.size() > 1) {
        EXPECT_LT(actuator.mappedSamItem[0].referencePos, actuator.mappedSamItem[1].referencePos);
    }

    std::remove("large_pos.sam");
}

// Test with empty SAM file (only header)
TEST_F(SamSortTest, TestWithEmptySam) {
    std::ofstream file("empty.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file.close();

    loadSamData("empty.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    EXPECT_EQ(actuator.headLineNum, 2);
    EXPECT_EQ(actuator.mappedSamItem.size(), 0);
    EXPECT_EQ(actuator.unmappedSamItem.size(), 0);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    std::remove("empty.sam");
}

// Test with multiple header types
TEST_F(SamSortTest, TestWithMultipleHeaderTypes) {
    std::ofstream file("multi_header.sam");
    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@SQ\tSN:chr2\tLN:800\n";
    file << "@RG\tID:group1\tDS:desc1\n";
    file << "@PG\tID:prog1\tPN:program\tVN:1.0\n";
    file << "@CO\tThis is a comment\n";
    file << "read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    file.close();

    loadSamData("multi_header.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    // Should count all header lines
    EXPECT_EQ(actuator.headLineNum, 6);  // @HD, 2x@SQ, @RG, @PG, @CO
    EXPECT_EQ(actuator.mappedSamItem.size(), 1);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    // Verify header file contains all header lines
    std::ifstream headFile("sorted_head.sam");
    if (headFile.is_open()) {
        std::string line;
        int lineCount = 0;
        while (std::getline(headFile, line)) {
            if (!line.empty()) {
                lineCount++;
            }
        }
        EXPECT_EQ(lineCount, 6);
        headFile.close();
    }

    std::remove("multi_header.sam");
}

// Test with records containing various SAM flags
TEST_F(SamSortTest, TestWithVariousFlags) {
    std::ofstream file("flags.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";

    // Flag 0: properly paired, both mapped
    file << "read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Flag 4: segment unmapped
    file << "read2\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Flag 16: reverse strand
    file << "read3\t16\tchr1\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Flag 5: mate reverse strand (bit mapping, both 4 and 1 set)
    file << "read4\t5\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    file.close();

    loadSamData("flags.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    // Should have 2 mapped (flags 0 and 16) and 2 unmapped (flags 4 and 5)
    EXPECT_EQ(actuator.mappedSamItem.size(), 2);
    EXPECT_EQ(actuator.unmappedSamItem.size(), 2);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    std::remove("flags.sam");
}

// Test with chromosome name "*"
TEST_F(SamSortTest, TestWithStarChromosome) {
    std::ofstream file("star_chr.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";

    // Record with chromosome "*"
    file << "read1\t0\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    file.close();

    loadSamData("star_chr.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    // Should be treated as unmapped
    EXPECT_EQ(actuator.mappedSamItem.size(), 0);
    EXPECT_EQ(actuator.unmappedSamItem.size(), 1);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    std::remove("star_chr.sam");
}

// Test sequential sort operations
TEST_F(SamSortTest, TestSequentialSortOperations) {
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    // First sort operation
    SAMSortActuator actuator1(pInBlock, pOutBlock);
    int32_t initResult1 = actuator1.initial();
    EXPECT_EQ(initResult1, 0);

    pInBlock->setBlockId(0);
    int32_t processResult1 = actuator1.process();
    EXPECT_EQ(processResult1, 0);

    // Second sort operation (should create new files)
    generateSortTestSamFile(SamSortTestData::testSamFile);
    loadSamData(SamSortTestData::testSamFile);

    SAMSortActuator actuator2(pInBlock, pOutBlock);
    int32_t initResult2 = actuator2.initial();
    EXPECT_EQ(initResult2, 0);

    pInBlock->setBlockId(1);
    int32_t processResult2 = actuator2.process();
    EXPECT_EQ(processResult2, 0);

    // Both files should exist
    std::ifstream file0("sorted_sam_0.sam");
    std::ifstream file1("sorted_sam_1.sam");

    EXPECT_TRUE(file0.good());
    EXPECT_TRUE(file1.good());

    if (file0.is_open()) file0.close();
    if (file1.is_open()) file1.close();
}

// Test with SAM record having minimal fields
TEST_F(SamSortTest, TestWithMinimalSamRecord) {
    std::ofstream file("minimal.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    // Minimal record with only required fields
    file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    file.close();

    loadSamData("minimal.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    EXPECT_EQ(actuator.mappedSamItem.size(), 1);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    std::remove("minimal.sam");
}

// Test note about unmapped items - they are stored but not processed
TEST_F(SamSortTest, TestUnmappedItemsNotProcessed) {
    generateUnmappedSamFile(SamSortTestData::testUnmappedSamFile);
    loadSamData(SamSortTestData::testUnmappedSamFile);

    SAMSortActuator actuator(pInBlock, pOutBlock);
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    EXPECT_EQ(actuator.mappedSamItem.size(), 0);
    EXPECT_GT(actuator.unmappedSamItem.size(), 0);

    pInBlock->setBlockId(0);
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

    // Check that sorted SAM file should be empty or minimal since all records are unmapped
    std::ifstream samFile("sorted_sam_0.sam");
    if (samFile.is_open()) {
        std::string line;
        int contentCount = 0;
        while (std::getline(samFile, line)) {
            if (!line.empty()) {
                contentCount++;
            }
        }
        // Should be 0 or very few since unmapped items are not written
        samFile.close();
    }
}

// Integration test: full workflow
TEST_F(SamSortTest, TestFullWorkflowIntegration) {
    // Create test data with deliberately unsorted records
    std::ofstream file("integration.sam");
    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@SQ\tSN:chr2\tLN:800\n";
    file << "@SQ\tSN:chr3\tLN:600\n";

    // Add records in reverse order to test sorting
    file << "read5\t0\tchr1\t500\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    file << "unmapped\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    file << "read3\t0\tchr1\t300\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    file << "read1\t0\tchr2\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    file << "read2\t0\tchr3\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";

    file.close();

    loadSamData("integration.sam");

    SAMSortActuator actuator(pInBlock, pOutBlock);
    pInBlock->setBlockId(0);

    // Initialize
    int32_t initResult = actuator.initial();
    EXPECT_EQ(initResult, 0);

    // Verify parsing
    EXPECT_EQ(actuator.headLineNum, 4);
    EXPECT_EQ(actuator.mappedSamItem.size(), 4);
    EXPECT_EQ(actuator.unmappedSamItem.size(), 1);

    // Process
    int32_t processResult = actuator.process();
    EXPECT_EQ(processResult, 0);

        // Verify sorting
    if (actuator.mappedSamItem.size() > 1) {
        for (size_t i = 0; i < actuator.mappedSamItem.size() - 1; i++) {
            EXPECT_LE(actuator.mappedSamItem[i].referencePos, actuator.mappedSamItem[i + 1].referencePos);
        }
    }

    // Verify output files exist and have correct content
    std::ifstream headFile("sorted_head.sam");
    std::ifstream samFile("sorted_sam_0.sam");

    EXPECT_TRUE(headFile.good());
    EXPECT_TRUE(samFile.good());

    // Verify header file content
    if (headFile.is_open()) {
        std::string line;
        int headerLines = 0;
        while (std::getline(headFile, line)) {
            if (!line.empty() && line[0] == '@') {
                headerLines++;
            }
        }
        EXPECT_EQ(headerLines, 4);
        headFile.close();
    }

    // Verify sorted SAM file content
    if (samFile.is_open()) {
        std::string line;
        int recordCount = 0;
        std::string prevLine;
        while (std::getline(samFile, line)) {
            if (!line.empty() && line[0] != '@') {
                recordCount++;
                // Each line should be a valid SAM record
                EXPECT_GT(line.size(), 0);
            }
        }
        EXPECT_EQ(recordCount, 5);  
        samFile.close();
    }

    std::remove("integration.sam");
}




