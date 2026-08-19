
/*
 * sam_actuator_testcase.cpp - Test cases for SAM actuator functionality
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

#define private public
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include <compress_engine.h>
#undef private

namespace SamTestData {
    const std::string testSamFile = "test.sam";
    const std::string testFastqFile = "test.fa";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamActuatorTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

        // Ensure SamInfo state is clean
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();

		pInBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);

        generateSamFile(SamTestData::testSamFile);

        ConfigManager::getInstance().logLevel = LogLevel::WARNING;
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
        // Clean up all potentially generated test files
        std::remove(SamTestData::testSamFile.c_str());
        std::remove("./test/test.sam");
        std::remove("../test/test.sam");
        std::remove("test_reference.fa");
        std::remove("./test/test_reference.fa");
        std::remove("../test/test_reference.fa");
	}

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        // Try multiple possible paths, starting with current directory and test directory
        std::vector<std::string> paths = {
            filename,           // current directory
            "./test/" + filename, // test directory
            "../test/" + filename // test directory above build directory
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
            // If all paths fail, use test_data/test.sam as final fallback
            pIoReader = new FileReader("test_data/test.sam");
            if (pIoReader->openIO() != 0) {
                delete pIoReader;
                return; // Cannot open file
            }
        }

        /*
         * 直接按原始内容载入整块（含 @SQ 头部），并记录换行符位置。
         * 头部在 reader 里独立成块后，actuator 测试需要一个自含头部的完整块。
         */
        const size_t blockSize = pInBlock->getBlockSize();
        const size_t readLen = pIoReader->readIO(pInBlock->getBuffer(), blockSize);
        pInBlock->setDataLen((int64_t)readLen);
        pInBlock->setBlockType(SAM);
        std::vector<size_t>& npos = pInBlock->getNpos();
        const char* buf = reinterpret_cast<const char*>(pInBlock->getBuffer());
        for (size_t i = 0; i < readLen; ++i) {
            if (buf[i] == '\n') {
                npos.push_back(i);
            }
        }

        delete pIoReader;
    }

    void generateSamFile(const std::string& filename) {
        // Directly write real SAM data from test_data, not reading from file
        std::ofstream file(filename);
        if (!file.is_open()) {
            // If cannot open in current directory, try to create in test directory
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // If still fails, try to create in build directory
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }

        if (!file.is_open()) {
            return; // Cannot create file
        }

        // Generate standard SAM file content, all records are 76bp length
        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";
        file << "@SQ\tSN:chr3\tLN:339\n";
        file << "@SQ\tSN:chr4\tLN:339\n";
        file << "@SQ\tSN:chrX\tLN:340\n";
        file << "@SQ\tSN:chrY\tLN:339\n";
        file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188\tCL:bwa mem test_data/reference.fa test_data/test_reads.fastq\n";

        // All records are 76M, sequence length and quality values are both 76
        file << "read1_11.1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read3_11.2\t0\tchr1\t153\t60\t74M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCCGCAAG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:3\tMD:Z:37C37C0\tAS:i:72\tXS:i:0\n";
        file << "read4_11.3\t0\tchr2\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGACAATCGCCGGAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read5_11.4\t0\tchr2\t77\t60\t76M\t*\t0\t0\tATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAACGCAATCGATCGA\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read6_11.5\t16\tchr3\t1\t60\t76M\t*\t0\t0\tCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCACAATCGGCGGCCGAATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read7_11.6\t0\tchr3\t77\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATACGCAATCGGCGGCCGAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:74A0\tAS:i:75\tXS:i:0\n";
        file << "read8_11.7\t16\tchr4\t1\t60\t76M\t*\t0\t0\tTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGCAACAATCGGCGGCC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read9_11.8\t0\tchrX\t1\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGACGCAATCGCCGGCCTAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read10_11.9\t16\tchrY\t1\t60\t75M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAACGCAATCGCCGGCCTTGTTC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read10_11.10\t16\tchrY\t1\t60\t75M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAACGCAATCGCCGGCCTTGTTC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";

        file.close();
    }

    void generateFasta(const std::string& filename) {
        std::ifstream srcFile("test_data/reference.fa");
        if (!srcFile.is_open()) {
            // If cannot open test_data/reference.fa, generate reference genome with real ATCG sequences
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
                return; // Cannot create file
            }

            // Generate reference genome with real ATCG sequences
            file << ">chr1\n";
            for (int i = 0; i < 17; i++) {
                file << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
            }
            file << "ATCGATCGATCGATCGATCG\n";  // Pad to appropriate length

            file << ">chr2\n";
            for (int i = 0; i < 16; i++) {
                file << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
            }
            file << "GCTAGCTAGCTAGCTAGC\n";

            file << ">chr3\n";
            for (int i = 0; i < 17; i++) {
                file << "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
            }
            file << "TACGTACGTACGTACGTACGT\n";

            file << ">chr4\n";
            for (int i = 0; i < 17; i++) {
                file << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
            }
            file << "CATGCATGCATGCATGCATGC\n";

            file << ">chrX\n";
            for (int i = 0; i < 17; i++) {
                file << "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
            }
            file << "ATGCATGCATGCATGCATGC\n";

            file << ">chrY\n";
            for (int i = 0; i < 17; i++) {
                file << "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\n";
            }
            file << "GCATGCATGCATGCATGCATG\n";

            file.close();
            return;
        }

            // If can open test_data/reference.fa, copy its content
        std::string content((std::istreambuf_iterator<char>(srcFile)),
                           std::istreambuf_iterator<char>());
        srcFile.close();

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
            return; // Cannot create file
        }

        file << content;
        file.close();
    }

    // Create Reference object for testing
    Reference createTestReference() {
        std::string tempFastaFile = "test_reference.fa";
        generateFasta(tempFastaFile);
        Reference ref(tempFastaFile, 1);
        ref.makeIndex();
        return ref;
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;

    // Mapping data for testing
    std::map<uint32_t, uint16_t> mappedFlag;
    std::map<uint32_t, uint64_t> mappedPos;
};

TEST_F(SamActuatorTest, testPreAnalysis) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.headEndLine, 8);  // Header has 8 lines (@HD, 6@SQ, 1@PG, may have other header lines)
    EXPECT_EQ(actuator.contentPos.size(), 10);  // Has 10 alignment records
}

TEST_F(SamActuatorTest, testPreAnalysisIdInvalid) {
    // Generate large SAM file for performance testing
    std::ofstream file("idinvlid_pre_analysis.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:248956422\n";
    file << "@SQ\tSN:chr2\tLN:242193529\n";
    file << "@SQ\tSN:chr3\tLN:198295559\n";
    file << "@SQ\tSN:chr4\tLN:190214555\n";
    file << "@SQ\tSN:chr5\tLN:181538259\n";
    file << "read1_11.SRR001.1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
    file << "read3_11.SRR001.2\t0\tchr1\t153\t60\t74M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCCGCAAG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:3\tMD:Z:37C37C0\tAS:i:72\tXS:i:0\n";
    file << "read4_11.SRR001.3\t0\tchr2\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGACAATCGCCGGAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
    file << "read5_11.SRR001.4\t0\tchr2\t77\t60\t76M\t*\t0\t0\tATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAACGCAATCGATCGA\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
    file << "read6_11.5\t16\tchr3\t1\t60\t76M\t*\t0\t0\tCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCACAATCGGCGGCCGAATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";

    file.close();

    loadSamData("idinvlid_pre_analysis.sam");
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);
    int32_t result = actuator.preAnalysis();
    std::remove("idinvlid_pre_analysis.sam");
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.idPosLength, UINT32_MAX);
}

TEST_F(SamActuatorTest, testCompressQuality) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressQuality method with field index 10 (QUAL field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressQuality(10, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressChrName) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressChrName method with field index 2 (RNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressChrName(2, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressWithoutRef) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressWithoutRef
    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testCompressWithRef) {
    loadSamData(SamTestData::testSamFile);
    Reference refGene = createTestReference();
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine, &refGene);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressWithRef
    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testCompressBaseWithoutRef) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressBaseWithoutRef method with field index 9 (SEQ field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressBaseWithoutRef(9, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressBaseWithRef) {
    loadSamData(SamTestData::testSamFile);
    Reference refGene = createTestReference();
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine, &refGene);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressBaseWithRef method with field index 9 (SEQ field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressBaseWithRef(9, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressIdFieldSplit) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressIdFieldSplit method with field index 0 (QNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressIdFieldSplit(fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressIdFieldQname) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressIdFieldQname method with field index 0 (QNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressIdFieldQname(fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressIdFieldInAll) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressIdFieldInAll method with field index 0 (QNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressIdFieldInAll(fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressRegularField) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compressRegularField method with field index 6 (PNEXT field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressRegularField(6, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testSetReference) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    Reference refGene = createTestReference();

    // Test setReference
    actuator.setReference(&refGene);

    // Verify reference is set by trying to compress
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testDecompress) {
    loadSamData(SamTestData::testSamFile);

    // Create SamActuator object for compression
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    pInBlock->reset();

    // Copy compressed output Block content to new input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    // Reset output Block
    pOutBlock->reset();

    // Create SamActuator object for decompression
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine);

    // Decompression doesn't need preAnalysis, directly call decompress
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // Basic check: ensure decompression produced data
    EXPECT_GT(pOutBlock->getDataLen(), 0);
}

/*
 * 配对读的字节级往返：PNEXT 走差值编码、TLEN 走推算（异常单独存）。逐字节比较
 * 解压结果与原始数据，验证差值还原、伙伴索引推算和异常回填三者叠加仍无损。
 */
TEST_F(SamActuatorTest, testPNextDeltaTlenRoundTripByteExact) {
    /* 用带配对读的 SAM 覆盖 fixture 生成的文件，便于驱动 TLEN 推算路径。 */
    std::string seq = std::string(76, 'A');
    std::string qual = std::string(76, '!');
    std::string seq2 = std::string(50, 'G');
    std::string qual2 = std::string(50, 'I');
    std::ofstream file(SamTestData::testSamFile);
    ASSERT_TRUE(file.is_open());
    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000000\n";
    file << "@SQ\tSN:chr2\tLN:1000000\n";
    /* 正常配对（99/147），RNEXT 用 = 表示同参照，TLEN 与公式一致。 */
    file << "read1\t99\tchr1\t100\t60\t76M\t=\t300\t276\t" << seq << "\t" << qual << "\tNM:i:0\n";
    file << "read2\t147\tchr1\t300\t60\t76M\t=\t100\t-276\t" << seq << "\t" << qual << "\tNM:i:0\n";
    /* 50M 短读配对，CIGAR 参考跨度即 50。 */
    file << "read3\t99\tchr1\t1000\t60\t50M\t=\t1200\t250\t" << seq2 << "\t" << qual2 << "\tNM:i:0\n";
    file << "read4\t147\tchr1\t1200\t60\t50M\t=\t1000\t-250\t" << seq2 << "\t" << qual2 << "\tNM:i:0\n";
    /* 不配对（FLAG=0）且 TLEN 非 0：推算为 0，这一行必须作为异常存下来。 */
    file << "read5\t0\tchr1\t2000\t60\t76M\t*\t0\t123\t" << seq << "\t" << qual << "\tNM:i:0\n";
    file.close();

    loadSamData(SamTestData::testSamFile);
    std::string original((char*)pInBlock->getBuffer(), pInBlock->getDataLen());

    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::string roundtrip((char*)pOutBlock->getBuffer(), pOutBlock->getDataLen());
    EXPECT_EQ(roundtrip, original);
}

TEST_F(SamActuatorTest, testDecompressWithRef) {
    loadSamData(SamTestData::testSamFile);

    Reference reference = createTestReference();

    // Create SamActuator object for compression
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    pInBlock->reset();

    // Copy compressed output Block content to new input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    // Reset output Block
    pOutBlock->reset();

    // Create SamActuator object for decompression
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);

    // Decompression doesn't need preAnalysis, directly call decompress
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // Basic check: ensure decompression produced data
    EXPECT_GT(pOutBlock->getDataLen(), 0);
}


TEST_F(SamActuatorTest, testCigarParse) {
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Helper function: convert string to uint8_t* to pass to parseCigar
    auto testParse = [](SamCodecActuator& a, const std::string& cigar) -> uint32_t {
        return a.parseCigar(const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(cigar.c_str())), cigar.length());
    };

    // Test basic M operations
    EXPECT_EQ(testParse(actuator, "100M"), 100);
    EXPECT_EQ(testParse(actuator, "76M"), 76);
    EXPECT_EQ(testParse(actuator, "1M"), 1);

    // Test I operations
    EXPECT_EQ(testParse(actuator, "10I"), 10);
    EXPECT_EQ(testParse(actuator, "5I"), 5);

    // Test S operations
    EXPECT_EQ(testParse(actuator, "15S"), 15);
    EXPECT_EQ(testParse(actuator, "3S"), 3);

    // Test = operations
    EXPECT_EQ(testParse(actuator, "50="), 50);
    EXPECT_EQ(testParse(actuator, "25="), 25);

    // Test X operations
    EXPECT_EQ(testParse(actuator, "30X"), 30);
    EXPECT_EQ(testParse(actuator, "12X"), 12);

    // Test composite operations - only accumulate M、I、S、=、X
    EXPECT_EQ(testParse(actuator, "10M5I3S2="), 20);  // 10+5+3+2 = 20
    EXPECT_EQ(testParse(actuator, "5M10I15S5X10="), 45);  // 5+10+15+5+10 = 45
    EXPECT_EQ(testParse(actuator, "1M1I1S1=1X"), 5);  // 1+1+1+1+1 = 5

    // Test composite CIGAR with other operations - other operations should be ignored
    EXPECT_EQ(testParse(actuator, "10M5D3I2H10N"), 13);  // Only count 10M+3I = 13, ignore 5D, 2H, 10N
    EXPECT_EQ(testParse(actuator, "5M10D5I10P5N5S"), 15);  // Only count 5M+5I+5S = 15
    EXPECT_EQ(testParse(actuator, "100M50D50N"), 100);  // Only count 100M, ignore 50D and 50N

    // Test real-world CIGAR strings
    EXPECT_EQ(testParse(actuator, "76M"), 76);  // Complete match
    EXPECT_EQ(testParse(actuator, "3S73M"), 76);  // 3 soft clipping + 73 matches
    EXPECT_EQ(testParse(actuator, "10M5I60M5D"), 75);  // 10+5+60 = 75, ignore 5D
    EXPECT_EQ(testParse(actuator, "1S20M1I30M1D10M1S"), 63);  // 1+20+1+30+10+1 = 63, ignore 1D

    // Test edge cases
    EXPECT_EQ(testParse(actuator, ""), 0);  // Empty string

    // Test CIGAR with only non-counting operations
    EXPECT_EQ(testParse(actuator, "100D"), 0);  // Only deletions, should not accumulate
    EXPECT_EQ(testParse(actuator, "50N"), 0);   // Only reference skips, should not accumulate
    EXPECT_EQ(testParse(actuator, "10H5P"), 0);  // Only hard clipping and padding, should not accumulate

    // Test complex real-world CIGAR scenarios
    EXPECT_EQ(testParse(actuator, "1S10M1I10M1D10M1I10M1D10M1D10M1I10M1S1H"), 75);  // Complex alignment scenario
    EXPECT_EQ(testParse(actuator, "35M1I39M"), 75);  // Alignment with insertion in middle
    EXPECT_EQ(testParse(actuator, "2S50M2I20M1D5M3S"), 82);  // Soft clipping on both ends

    // Test large numbers
    EXPECT_EQ(testParse(actuator, "1000M"), 1000);
    EXPECT_EQ(testParse(actuator, "10000M500I200S100="), 10800);

    // Test mixed case (although SAM specification usually uses uppercase)
    EXPECT_EQ(testParse(actuator, "10m5i3s2="), 20);  // Lowercase should also work
    EXPECT_EQ(testParse(actuator, "10M5i3S2x"), 20);  // Mixed case

    // Test invalid CIGAR formats
    EXPECT_EQ(testParse(actuator, "M"), 0);  // Missing number
    EXPECT_EQ(testParse(actuator, "invalid"), 0);  // Completely invalid
    EXPECT_EQ(testParse(actuator, "10M5"), 10);  // Missing operator at end
    EXPECT_EQ(testParse(actuator, "M10I"), 10);  // Missing number at start
}

TEST_F(SamActuatorTest, testBuildSamIndexDisabled) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    para.isMakeIndex = false;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test buildSamIndex when disabled
    result = actuator.buildSamIndex();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testBuildSamIndexSuccess) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    para.isMakeIndex = true;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Test buildSamIndex with valid sorted SAM data
    result = actuator.buildSamIndex();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testBuildSamIndexUnsorted) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    para.isMakeIndex = true;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Manually set up mapping data to simulate unsorted order
    // First read maps to chr2 position 100
    actuator.mappedFlag[8] = 0;
    actuator.mappedChr[8] = 1;  // chr2
    actuator.mappedPos[8] = 100;
    actuator.mappedFlag[9] = 0;
    actuator.mappedChr[9] = 0;  // chr1
    actuator.mappedPos[9] = 1;

    // Test buildSamIndex with unsorted data should fail
    result = actuator.buildSamIndex();
    EXPECT_EQ(result, -1);
}

TEST_F(SamActuatorTest, testBuildSamIndexSkipUnmapped) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    para.isMakeIndex = true;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Set up mapping data - some reads are unmapped (flag & 0x04)
    actuator.mappedFlag[8] = 0x04;  // Unmapped
    actuator.mappedChr[8] = 0;
    actuator.mappedPos[8] = 0;
    actuator.mappedFlag[9] = 0;
    actuator.mappedChr[9] = 1;  // chr2
    actuator.mappedPos[9] = 100;

    // Test buildSamIndex should skip unmapped reads
    result = actuator.buildSamIndex();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testBuildSamIndexSkipInvalidChr) {
    loadSamData(SamTestData::testSamFile);
    PbgzParameter para;
    para.isMakeIndex = true;
    CompressEngine engine(para);
    SamCodecActuator actuator(pInBlock, pOutBlock, &engine);

    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);

    // Set up mapping data with invalid chromosome indices
    actuator.mappedFlag[8] = 0;
    actuator.mappedChr[8] = 0xFFFF;  // Unmapped
    actuator.mappedPos[8] = 0;
    actuator.mappedFlag[9] = 0;
    actuator.mappedChr[9] = 0xFFFE;  // '=' marker
    actuator.mappedPos[9] = 100;

    // Test buildSamIndex should skip reads with invalid chr indices
    result = actuator.buildSamIndex();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testDecompressHeaderWithOutputBlock) {
    loadSamData(SamTestData::testSamFile);

    Reference reference = createTestReference();

    // Create SamActuator object for compression
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    // Create a new block for decompressing header
    RoughIOBlock* headerBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
    ASSERT_NE(headerBlock, nullptr);

    // Reset and copy compressed data to new input block
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    // Create decompressor and init metadata
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);
    decompressor.initMetaInfo();

    // Test decompressHeader with custom output block
    result = decompressor.decompressHeader(headerBlock);
    EXPECT_EQ(result, 0);

    // Verify header was decompressed to headerBlock
    EXPECT_GT(headerBlock->getDataLen(), 0);

    // Verify header contains @HD
    std::string headerContent((char*)headerBlock->getBuffer(), headerBlock->getDataLen());
    EXPECT_TRUE(headerContent.find("@HD") != std::string::npos);

    // Verify header contains @SQ lines
    EXPECT_TRUE(headerContent.find("@SQ") != std::string::npos);

    delete headerBlock;
}

TEST_F(SamActuatorTest, testDecompressSamByFieldsWithOutputBlock) {
    loadSamData(SamTestData::testSamFile);

    Reference reference = createTestReference();

    // Create SamActuator object for compression
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    // Create a new block for decompressing SAM fields
    RoughIOBlock* samBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
    ASSERT_NE(samBlock, nullptr);

    // Reset and copy compressed data to new input block
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    // Create decompressor
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);
    decompressor.initMetaInfo();

    // Decompress header first
    result = decompressor.decompressHeader(samBlock);
    EXPECT_EQ(result, 0);

    // Test decompressSamByFields with custom output block
    result = decompressor.decompressSamByFields(samBlock);
    EXPECT_EQ(result, 0);

    // Verify SAM data was decompressed to samBlock
    EXPECT_GT(samBlock->getDataLen(), 0);

    // Verify decompressed data contains SAM records
    std::string samContent((char*)samBlock->getBuffer(), samBlock->getDataLen());
    EXPECT_TRUE(samContent.find("read") != std::string::npos || samContent.find("\t0\t") != std::string::npos);

    // Verify SAM records have proper tabs (field separators)
    EXPECT_TRUE(samContent.find("\t") != std::string::npos);

    delete samBlock;
}

TEST_F(SamActuatorTest, testDecompressWithRefGenePosScenario) {
    loadSamData(SamTestData::testSamFile);

    Reference reference = createTestReference();

    // Initialize SamIndex with some test data for refGenePos scenario
    SamInfo::getInstance().clearChromosomeInfo();
    SamInfo::getInstance().resetChrIdCounter();

    // Create SamActuator object for compression
    PbgzParameter para;
    para.isMakeIndex = true;  // Enable index creation
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    // Build index (simulating what happens during compression)
    result = compressor.buildSamIndex();
    EXPECT_EQ(result, 0);

    pInBlock->reset();

    // Copy compressed output Block content to new input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    // Create separate blocks for header and data (simulating refGenePos scenario)
    RoughIOBlock* headerBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
    RoughIOBlock* dataBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
    ASSERT_NE(headerBlock, nullptr);
    ASSERT_NE(dataBlock, nullptr);

    pOutBlock->reset();

    // Create SamActuator object for decompression
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);
    decompressor.initMetaInfo();

    // Simulate refGenePos scenario: decompress header to separate block
    result = decompressor.decompressHeader(headerBlock);
    EXPECT_EQ(result, 0);
    EXPECT_GT(headerBlock->getDataLen(), 0);

    // Verify header content
    std::string headerContent((char*)headerBlock->getBuffer(), headerBlock->getDataLen());
    EXPECT_TRUE(headerContent.find("@HD") != std::string::npos);
    EXPECT_TRUE(headerContent.find("@SQ") != std::string::npos);

    // Simulate refGenePos scenario: decompress SAM fields to separate data block
    result = decompressor.decompressSamByFields(dataBlock);
    EXPECT_EQ(result, 0);
    EXPECT_GT(dataBlock->getDataLen(), 0);

    // Verify data content
    std::string dataContent((char*)dataBlock->getBuffer(), dataBlock->getDataLen());
    EXPECT_TRUE(dataContent.find("\t0\t") != std::string::npos);
    EXPECT_TRUE(dataContent.find("\t60\t") != std::string::npos);

    delete headerBlock;
    delete dataBlock;
}

TEST_F(SamActuatorTest, testDecompressWithRefGenePosFiltering) {
    // Generate SAM file with specific positions for testing filtering
    std::ofstream file("filter_test.sam");
    if (!file.is_open()) {
        std::string testPath = "./test/filter_test.sam";
        file.open(testPath);
        if (!file.is_open()) {
            testPath = "../test/filter_test.sam";
            file.open(testPath);
        }
    }

    if (!file.is_open()) {
        return;
    }

    // Generate SAM file with reads at different positions
    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@SQ\tSN:chr2\tLN:1000\n";
    file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n";

    // Read at position 1 (will be included if range is 1-100)
    file << "read_at_pos1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Read at position 50 (will be included if range is 1-100)
    file << "read_at_pos50\t0\tchr1\t50\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Read at position 100 (will be included if range is 1-100)
    file << "read_at_pos100\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Read at position 150 (will NOT be included if range is 1-100)
    file << "read_at_pos150\t0\tchr1\t150\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Read at position 200 (will NOT be included if range is 1-100)
    file << "read_at_pos200\t0\tchr1\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Read on different chromosome (will NOT be included)
    file << "read_on_chr2\t0\tchr2\t50\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";
    // Unmapped read (will NOT be included)
    file << "read_unmapped\t4\t*\t0\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\ttRG:Z:test\n";

    file.close();

    // Load the test file
    loadSamData("filter_test.sam");

    Reference reference = createTestReference();

    // Create SamActuator object for compression with filtering range
    PbgzParameter para;
    para.isMakeIndex = true;
    para.refeGenePos = "chr1:1-100";  // Set filter before compression
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    // Compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    // Build index
    result = compressor.buildSamIndex();
    EXPECT_EQ(result, 0);

    pInBlock->reset();

    // Copy compressed output Block content to new input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    pOutBlock->reset();

    // Create SamActuator object for decompression
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);

    // Decompress with filtering
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // Verify filtered results
    EXPECT_GT(pOutBlock->getDataLen(), 0);
    std::string outputContent((char*)pOutBlock->getBuffer(), pOutBlock->getDataLen());

    // Verify refPosChrIndex, refPosBegin, refPosEnd were parsed correctly
    // chr1 should have index 0 (registered in SamInfo from @SQ header)
    uint16_t chrIndex = SamInfo::getInstance().getChrNameIndex("chr1");
    EXPECT_NE(chrIndex, 65535) << "chr1 should be registered";
    EXPECT_EQ(decompressor.refPosChrIndex, chrIndex);
    EXPECT_EQ(decompressor.refPosBegin, 1);
    EXPECT_EQ(decompressor.refPosEnd, 100);

    // Should contain reads at positions 1, 50, 100
    EXPECT_TRUE(outputContent.find("read_at_pos1") != std::string::npos) << "Should contain read_at_pos1";
    EXPECT_TRUE(outputContent.find("read_at_pos50") != std::string::npos) << "Should contain read_at_pos50";
    EXPECT_TRUE(outputContent.find("read_at_pos100") != std::string::npos) << "Should contain read_at_pos100";

    // Should NOT contain reads at positions 150, 200, chr2, or unmapped
    EXPECT_FALSE(outputContent.find("read_at_pos150") != std::string::npos) << "Should NOT contain read_at_pos150";
    EXPECT_FALSE(outputContent.find("read_at_pos200") != std::string::npos) << "Should NOT contain read_at_pos200";
    EXPECT_FALSE(outputContent.find("read_on_chr2") != std::string::npos) << "Should NOT contain read_on_chr2";
    EXPECT_FALSE(outputContent.find("read_unmapped") != std::string::npos) << "Should NOT contain read_unmapped";

    std::remove("filter_test.sam");
}

TEST_F(SamActuatorTest, testDecompressWithRefGenePosAllIncluded) {
    // Generate SAM file with specific positions for testing filtering
    std::ofstream file("range_test.sam");
    if (!file.is_open()) {
        std::string testPath = "./test/range_test.sam";
        file.open(testPath);
        if (!file.is_open()) {
            testPath = "../test/range_test.sam";
            file.open(testPath);
        }
    }

    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n";

    // Read at position 1
    file << "read_1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Read at position 200
    file << "read_200\t0\tchr1\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Read at position 500
    file << "read_500\t0\tchr1\t500\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Read at position 1000
    file << "read_1000\t0\tchr1\t1000\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    file.close();

    loadSamData("range_test.sam");

    Reference reference = createTestReference();

    PbgzParameter para;
    para.isMakeIndex = true;
    para.refeGenePos = "chr1:1-10000";  // Wide range to include all reads
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    result = compressor.compress();
    EXPECT_EQ(result, 0);

    result = compressor.buildSamIndex();
    EXPECT_EQ(result, 0);

    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // Verify refPosChrIndex, refPosBegin, refPosEnd are set correctly
    uint16_t chrIndex = SamInfo::getInstance().getChrNameIndex("chr1");
    EXPECT_NE(chrIndex, 65535) << "chr1 should be registered";
    EXPECT_EQ(decompressor.refPosChrIndex, chrIndex);
    EXPECT_EQ(decompressor.refPosBegin, 1);
    EXPECT_EQ(decompressor.refPosEnd, 10000);

    EXPECT_GT(pOutBlock->getDataLen(), 0);
    std::string outputContent((char*)pOutBlock->getBuffer(), pOutBlock->getDataLen());

    // All reads should be included
    EXPECT_TRUE(outputContent.find("read_1") != std::string::npos) << "Should contain read_1";
    EXPECT_TRUE(outputContent.find("read_200") != std::string::npos) << "Should contain read_200";
    EXPECT_TRUE(outputContent.find("read_500") != std::string::npos) << "Should contain read_500";
    EXPECT_TRUE(outputContent.find("read_1000") != std::string::npos) << "Should contain read_1000";

    std::remove("range_test.sam");
}

TEST_F(SamActuatorTest, testDecompressWithRefGenePosNoneMatched) {
    // Test scenario where no reads match the filter criteria
    std::ofstream file("no_match_test.sam");
    if (!file.is_open()) {
        std::string testPath = "./test/no_match_test.sam";
        file.open(testPath);
        if (!file.is_open()) {
            testPath = "../test/no_match_test.sam";
            file.open(testPath);
        }
    }

    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@SQ\tSN:chr2\tLN:1000\n";
    file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n";

    // Read at position 500 on chr1
    file << "read_500\t0\tchr1\t500\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // Read at position 200 on chr2
    file << "read_chr2_200\t0\tchr2\t200\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    file.close();

    loadSamData("no_match_test.sam");

    Reference reference = createTestReference();

    PbgzParameter para;
    para.isMakeIndex = true;
    para.refeGenePos = "chr1:1-100";  // No reads match this range
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    result = compressor.compress();
    EXPECT_EQ(result, 0);

    result = compressor.buildSamIndex();
    EXPECT_EQ(result, 0);

    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // Should contain header but no reads
    std::string outputContent((char*)pOutBlock->getBuffer(), pOutBlock->getDataLen());
    EXPECT_FALSE(outputContent.find("@HD") != std::string::npos) << "Should not contain header";
    EXPECT_FALSE(outputContent.find("read_500") != std::string::npos) << "Should NOT contain read_500";

    // Verify refPosChrIndex is not 65535 (chr1 is registered), but no reads match range
    EXPECT_NE(decompressor.refPosChrIndex, 65535);
    EXPECT_EQ(decompressor.refPosBegin, 1);
    EXPECT_EQ(decompressor.refPosEnd, 100);

    std::remove("no_match_test.sam");
}

TEST_F(SamActuatorTest, testDecompressWithRefGenePosInvalidChr) {
    // Test with invalid chromosome name
    std::ofstream file("invalid_chr_test.sam");
    if (!file.is_open()) {
        std::string testPath = "./test/invalid_chr_test.sam";
        file.open(testPath);
        if (!file.is_open()) {
            testPath = "../test/invalid_chr_test.sam";
            file.open(testPath);
        }
    }

    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n";

    file << "read_1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

    file.close();

    loadSamData("invalid_chr_test.sam");

    Reference reference = createTestReference();

    PbgzParameter para;
    para.isMakeIndex = true;
    para.refeGenePos = "chr999:1-100";  // Non-existent chromosome
    CompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &reference);

    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);

    result = compressor.compress();
    EXPECT_EQ(result, 0);

    result = compressor.buildSamIndex();
    EXPECT_EQ(result, 0);

    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());

    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &reference);

    // Test decompress - invalid chromosome should result in no filtering
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);

    // When chromosome name is not found, refPosChrIndex remains 65535,
    // so no filtering should happen and all reads should be output
    std::string outputContent((char*)pOutBlock->getBuffer(), pOutBlock->getDataLen());
    EXPECT_TRUE(outputContent.find("read_1") != std::string::npos) << "Should contain read when chr name is invalid";

    // Verify refPosChrIndex is 65535 because chr999 doesn't exist
    EXPECT_EQ(decompressor.refPosChrIndex, 65535);

    std::remove("invalid_chr_test.sam");
}
