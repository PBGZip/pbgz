/*
 * compress_selector_testcase.cpp - Test cases for compression selector
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
#include <cstring>
#include <string>
#include <vector>
#include <cstdint>

#include "compress_selector.h"
#include "io_block.h"
#include "utils/memory_util.h"
#include "sam_info.h"
#include "coder/coder.h"

class CompressSelectorTest : public ::testing::Test {
protected:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();

        CompressionSelectorManager::getInstance().clearAllSelections();

        pBlock = new RoughIOBlock(8 << 20);
    }

    void TearDown() override {
        CompressionSelectorManager::getInstance().clearAllSelections();
        if (pBlock) {
            delete pBlock;
            pBlock = nullptr;
        }
    }

    void populateBlock(const std::string& samContent) {
        uint32_t len = static_cast<uint32_t>(samContent.length());
        memcpy(pBlock->getBuffer(), samContent.c_str(), len);
        pBlock->setDataLen(len);

        std::vector<uint32_t>& npos = pBlock->getNpos();
        npos.clear();
        for (uint32_t i = 0; i < len; ++i) {
            if (samContent[i] == '\n') {
                npos.push_back(i);
            }
        }
    }

    // SAM field: 0:QNAME 1:FLAG 2:RNAME 3:POS 4:MAPQ 5:CIGAR 6:RNEXT 7:PNEXT 8:TLEN 9:SEQ 10:QUAL 11:OPTION
    static std::string buildSamLine(const std::vector<std::string>& fields) {
        std::string line;
        for (size_t i = 0; i < fields.size(); ++i) {
            if (i > 0) line += '\t';
            line += fields[i];
        }
        line += '\n';
        return line;
    }

    RoughIOBlock* pBlock = nullptr;
};

TEST_F(CompressSelectorTest, TestFlagField) {
    // 2 lines with FLAG field
    std::string content = buildSamLine({"read1", "0",  "chr1", "1",  "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "16", "chr1", "77", "60", "76M", "*", "0", "0", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(1, pBlock);  // FLAG

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(1).size(), 0u);
    mgr.printCompressionStats();
}

TEST_F(CompressSelectorTest, TestPosField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1",    "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "152",  "60", "76M", "*", "0", "0", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(3, pBlock);  // POS

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(3).size(), 0u);
}

TEST_F(CompressSelectorTest, TestMapqField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "1", "30", "76M", "*", "0", "0", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(4, pBlock);  // MAPQ

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(4).size(), 0u);
}

TEST_F(CompressSelectorTest, TestCigarField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M",   "*", "0", "0", "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "1", "60", "35M1I40M", "*", "0", "0", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(5, pBlock);  // CIGAR

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(5).size(), 0u);
}

TEST_F(CompressSelectorTest, TestPnextField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0",   "0", "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "1", "60", "76M", "*", "152", "0", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(7, pBlock);  // PNEXT

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(7).size(), 0u);
}

TEST_F(CompressSelectorTest, TestTlenField) {
    std::string content = buildSamLine({"read1", "0",  "chr1", "1",  "60", "76M", "*", "0", "0",  "ACGT", "HHHH", "NM:i:1"})
                        + buildSamLine({"read2", "16", "chr1", "77", "60", "76M", "*", "0", "151", "TGCA", "IIII", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(8, pBlock);  // TLEN

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(8).size(), 0u);
}

TEST_F(CompressSelectorTest, TestQualField) {
    // SEQ and QUAL must have the same length for QUAL_GEN2
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGTACGTACGTACGT", "HHHHHHHHHHHHHHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "1", "60", "76M", "*", "0", "0", "NNNNNNNNNNNNNNNN", "!!!!!HHHHHHHHHHH", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(10, pBlock);  // QUAL

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(10).size(), 0u);
}

TEST_F(CompressSelectorTest, TestOptionField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1\tMD:Z:75A0"})
                        + buildSamLine({"read2", "0", "chr1", "1", "60", "76M", "*", "0", "0", "TGCA", "IIII", "NM:i:3\tMD:Z:37C37C0"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(11, pBlock);  // OPTION

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(11).size(), 0u);
}

TEST_F(CompressSelectorTest, TestUnconfiguredField) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"});

    populateBlock(content);

    // QNAME (0) is not in compressFieldConfg, should return early without crash
    SamCompressionSlector::testSamRegularFiled(0, pBlock);
    SamCompressionSlector::testSamRegularFiled(2, pBlock);  // RNAME
    SamCompressionSlector::testSamRegularFiled(6, pBlock);  // RNEXT
    SamCompressionSlector::testSamRegularFiled(9, pBlock);  // SEQ

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_EQ(mgr.getStatsForField(0).size(), 0u);
    EXPECT_EQ(mgr.getStatsForField(2).size(), 0u);
    EXPECT_EQ(mgr.getStatsForField(6).size(), 0u);
    EXPECT_EQ(mgr.getStatsForField(9).size(), 0u);
}

TEST_F(CompressSelectorTest, TestGetBestModel) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1",  "60", "76M", "*", "0", "0", "ACGTACGT", "HHHHHHHH", "NM:i:1"})
                        + buildSamLine({"read2", "0", "chr1", "77", "60", "76M", "*", "0", "0", "NNNNNNNN", "!!!!!!!!", "NM:i:2"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(3, pBlock);  // POS

    auto& mgr = CompressionSelectorManager::getInstance();
    CompressionModel best = mgr.getBestModelForField(3);
    // Should return a valid model (BWT_CM or AFFIX_MATCH), not crash
    EXPECT_TRUE(best == CompressionModel::MODEL_BWT_CM || best == CompressionModel::MODEL_AFFIX_MATCH);
}

TEST_F(CompressSelectorTest, TestPrintStatsEmpty) {
    // Should not crash when no data has been collected
    auto& mgr = CompressionSelectorManager::getInstance();
    mgr.printCompressionStats();
    SUCCEED();
}

TEST_F(CompressSelectorTest, TestClearAllSelections) {
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(1, pBlock);  // FLAG
    SamCompressionSlector::testSamRegularFiled(3, pBlock);  // POS

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(1).size(), 0u);
    EXPECT_GT(mgr.getStatsForField(3).size(), 0u);

    mgr.clearAllSelections();
    EXPECT_EQ(mgr.getStatsForField(1).size(), 0u);
    EXPECT_EQ(mgr.getStatsForField(3).size(), 0u);
}

TEST_F(CompressSelectorTest, TestSingleLineInsufficientFields) {
    // SAM line without enough fields (only 3 tabs = 4 fields)
    std::string content = "read1\t0\tchr1\t1\n";
    populateBlock(content);

    // FLAG is field 1, exists in this line (has enough tabs)
    SamCompressionSlector::testSamRegularFiled(1, pBlock);  // FLAG

    auto& mgr = CompressionSelectorManager::getInstance();
    EXPECT_GT(mgr.getStatsForField(1).size(), 0u);

    mgr.clearAllSelections();

    // POS is field 3, also exists
    SamCompressionSlector::testSamRegularFiled(3, pBlock);  // POS
    EXPECT_GT(mgr.getStatsForField(3).size(), 0u);

    mgr.clearAllSelections();

    // CIGAR is field 5, doesn't exist (only 3 tabs)
    SamCompressionSlector::testSamRegularFiled(5, pBlock);
    EXPECT_EQ(mgr.getStatsForField(5).size(), 0u);
}

TEST_F(CompressSelectorTest, TestQualWithMismatchedSeqLength) {
    // QUAL length (4) != SEQ length (3), these lines should be skipped
    std::string content = buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACG", "HHHH", "NM:i:1"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(10, pBlock);  // QUAL

    auto& mgr = CompressionSelectorManager::getInstance();
    // Should be empty because QUAL/SEQ length mismatch causes skip
    EXPECT_EQ(mgr.getStatsForField(10).size(), 0u);
}

TEST_F(CompressSelectorTest, TestHeaderLineSkipped) {
    // Header lines (starting with @) should be skipped
    std::string content =
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        + buildSamLine({"read1", "0", "chr1", "1", "60", "76M", "*", "0", "0", "ACGT", "HHHH", "NM:i:1"});

    populateBlock(content);

    SamCompressionSlector::testSamRegularFiled(1, pBlock);  // FLAG

    auto& mgr = CompressionSelectorManager::getInstance();
    // Should only have stats from the 1 data line (not the 2 header lines)
    EXPECT_GT(mgr.getStatsForField(1).size(), 0u);
}
