/*
 * sam_decompress_testcase.cpp - Test cases for SAM decompression functionality
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
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <random>
#include <iterator>

#define private public
#include "coder/coder.h"
#include "coder_json.h"
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include "utils/memory_util.h"
#undef private
#include <compress_engine.h>
#include <decompress_engine.h>
#include <sam_info.h>
#include <pbgz_index.h>

namespace SamDecompressData {
    const std::string testSamFile = "test.sam";
    const std::string testFastqFile = "test.fa";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamDecompressTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

		pInBlock = new RoughIOBlock(SamDecompressData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamDecompressData::MAX_BLOCK_SIZE);

        generateSamFile(SamDecompressData::testSamFile);

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
        // Clean up all possible generated test files
        std::remove(SamDecompressData::testSamFile.c_str());
        std::remove("./test/test.sam");
        std::remove("../test/test.sam");
        std::remove("test_reference.fa");
        std::remove("./test/test_reference.fa");
        std::remove("../test/test_reference.fa");
        std::remove("test_id_separators.sam");
        std::remove("test_missing_fields.sam");
        std::remove("test_flag_matching.sam");
        std::remove("test_base_fixed.sam");
        std::remove("test_base_variable.sam");
        std::remove("test_base_fastq3.sam");
        std::remove("test_optional_fields.sam");
        std::remove("test_no_optional_fields.sam");
        std::remove("test_field_count_too_few.sam");
        std::remove("test_field_missing.sam");
        std::remove("test_mixed.sam");
        std::remove("test_large.sam");
        std::remove("test_reference_roundtrip.sam");
        std::remove("test_unused_ref_no_sq.sam");

        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();
        SamIndex::getInstance().clear();
	}

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        // Try multiple possible paths, first try current directory and test directory
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
                return; // Unable to open file
            }
        }

        /*
         * Load the whole block directly from the raw content (including the @SQ header) and
         * record newline positions. Since the reader separates the header into its own block,
         * the actuator tests need a complete, self-contained block with the header included.
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

    /*
     * Load the block directly from the raw content and record newline positions without format
     * detection. Used to intentionally construct input that "looks like SAM but is invalid",
     * bypassing BlockFactory's strict format checks to directly verify the actuator's pre-analysis
     * rejection of invalid content.
     */
    void loadSamDataRaw(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        std::vector<std::string> paths = { filename, "./test/" + filename, "../test/" + filename };
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
        // Directly write real SAM data from test_data, not read from file
        std::ofstream file(filename);
        if (!file.is_open()) {
            // If unable to open in current directory, try to create in test directory
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // If still fails, try to create in build directory
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }

        if (!file.is_open()) {
            return; // Unable to create file
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
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read3\t0\tchr1\t153\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCCGCAAG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:3\tMD:Z:37C37C0\tAS:i:72\tXS:i:0\n";
        file << "read4\t0\tchr2\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGACAATCGCCGGAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read5\t0\tchr2\t77\t60\t76M\t*\t0\t0\tATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAACGCAATCGATCGA\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read6\t16\tchr3\t1\t60\t76M\t*\t0\t0\tCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCACAATCGGCGGCCGAATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read7\t0\tchr3\t77\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATACGCAATCGGCGGCCGAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:74A0\tAS:i:75\tXS:i:0\n";
        file << "read8\t16\tchr4\t1\t60\t76M\t*\t0\t0\tTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGCAACAATCGGCGGCC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read9\t0\tchrX\t1\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGACGCAATCGCCGGCCTAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read10\t0\tchrY\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAACGCAATCGCCGGCCTTGTTC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";

        file.close();
    }

    void generateFasta(const std::string& filename) {
        std::ifstream srcFile("reference.fa");
        if (!srcFile.is_open()) {
            // If unable to open reference.fa, generate reference genome with real ATCG sequences
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
                return; // Unable to create file
            }

            // Generate reference genome with real ATCG sequences
            file << ">chr1\n";
            for (int i = 0; i < 17; i++) {
                file << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
            }
            file << "ATCGATCGATCGATCGATCG\n";  // Fill to appropriate length

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
                file << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
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

        // If able to open reference.fa, copy its content
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
            return; // Unable to create file
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

    // Generate SAM file with various ID separators
    void generateSamFileWithSeparators(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";
        file << "@SQ\tSN:chr3\tLN:339\n";

        // Each read contains all separators, in consistent order: : . # / | - _ $
        file << "read1:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read2:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read3:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read4:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t4\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read5:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t5\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read6:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t6\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read7:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t7\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read8:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t8\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // Generate SAM file with missing fields (for testing segment compression logic)
    void generateSamFileWithMissingFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Normal record
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // Record with missing some fields (only ID and FLAG)
        file << "read2\t0\n";
        // Record with missing middle fields
        file << "read3\t0\tchr1\t2\t60\n";
        // Record with missing more fields
        file << "read4\t0\tchr2\t3\t76M\t*\t0\n";
        // Record with only ID
        file << "read5\n";

        file.close();
    }

    // Generate SAM file with different FLAG matching scenarios
    void generateSamFileWithFlagVariants(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Matching FLAG (65)
        file << "read1\t65\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // Unmatched FLAG (4)
        file << "read2\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:0\n";
        // Forward matching FLAG (0)
        file << "read3\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // Reverse matching FLAG (16)
        file << "read4\t16\tchr2\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // Generate SAM file with fixed-length Base fields
    void generateSamFileWithFixedBase(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // All sequences have same length (fixed length)
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // Generate SAM file with variable-length Base fields
    void generateSamFileWithVariableBase(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Different sequence lengths, ensure quality value count matches base count
        // First line: 50M CIGAR, 50 bases, 50 quality values
        std::string seq1 = "";
        for (int i = 0; i < 50; i++) {
            seq1 += "ATCG"[i % 4];
        }
        std::string qual1 = std::string(50, '!');
        file << "read1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t" << seq1 << "\t" << qual1 << "\tNM:i:1\n";
        
        // Second line: 75M CIGAR, 75 bases, 75 quality values
        std::string seq2 = "";
        for (int i = 0; i < 75; i++) {
            seq2 += "ATCG"[i % 4];
        }
        std::string qual2 = std::string(75, '!');
        file << "read2\t0\tchr1\t2\t60\t75M\t*\t0\t0\t" << seq2 << "\t" << qual2 << "\tNM:i:1\n";
        
        // Third line: 100M CIGAR, 100 bases, 100 quality values
        std::string seq3 = "";
        for (int i = 0; i < 100; i++) {
            seq3 += "ATCG"[i % 4];
        }
        std::string qual3 = std::string(100, '!');
        file << "read3\t0\tchr2\t3\t60\t100M\t*\t0\t0\t" << seq3 << "\t" << qual3 << "\tNM:i:1\n";

        file.close();
    }

    // Generate 3rd generation fastq style SAM file
    void generateSamFileWithFastq3Base(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Include long sequences and special bases, ensure quality value count matches base count
        file << "read1\t0\tchr1\t1\t60\t4000M\t*\t0\t0\t";
        for (int i = 0; i < 1000; i++) {
            file << "ATCG";
        }
        file << "\t";
        for (int i = 0; i < 4000; i++) {
            file << "!";
        }
        file << "\tNM:i:10\n";

        file.close();
    }

    // Generate SAM file with optional fields
    void generateSamFileWithOptionalFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Include multiple optional fields
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\tRG:Z:group1\tBC:Z:AGCTCT\tQT:Z:FFFFFFFFFFFFFF\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74T0\tAS:i:74\tXS:i:0\tRG:Z:group2\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:0\tMD:Z:76\tAS:i:76\tXS:i:0\tXA:Z:chr2,+5000,76M,1;chr3,+10000,76M,1;\n";

        file.close();
    }

    // Generate SAM file without optional fields
    void generateSamFileWithoutOptionalFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // Only required fields
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

        file.close();
    }

    // Compression and decompression helper function (optimized version: direct memory operations without files)
    void compressAndDecompress(const std::string& inputFile) {
        // Load input file
        loadSamData(inputFile);

        // Create compressor
        Reference ref = createTestReference();
        PbgzParameter para;
        DecompressEngine engine(para);
        SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
        ASSERT_EQ(compressor.preAnalysis(), 0);
        ASSERT_EQ(compressor.compress(), 0);

        // After compression, copy output from outBlock to inBlock
        pInBlock->reset();
        memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
        pInBlock->setDataLen(pOutBlock->getDataLen());
        pInBlock->setMetaLen(pOutBlock->getMetaLen());

        // Reset outBlock
        pOutBlock->reset();

        // Use these two blocks for decompression, only need to verify decompression return code
        SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);

        int32_t ret = decompressor.decompress();
        std::remove(inputFile.c_str());
        EXPECT_EQ(ret, 0);
    }

protected:
    RoughIOBlock* pInBlock = nullptr;
    RoughIOBlock* pOutBlock = nullptr;

    PbgzParameter para;
};

// Test 1: ID field covering all separator characters
TEST_F(SamDecompressTest, TestIDFieldSeparators) {
    generateSamFileWithSeparators("test_id_separators.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_id_separators.sam"));
}

// Test 2: Missing fields scenario, triggering segment compression logic
TEST_F(SamDecompressTest, TestMissingFieldsSegmentCompression) {
    generateSamFileWithMissingFields("test_missing_fields.sam");
        // Only the compression and decompression return codes are verified here; the file content is not compared
        // Load input file (invalid content is loaded as raw bytes to verify the actuator's pre-analysis rejection logic)
        loadSamDataRaw("test_missing_fields.sam");

        // Create compressor
        Reference ref = createTestReference();
        PbgzParameter para;
        CompressEngine engine(para);
        SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
        ASSERT_EQ(compressor.preAnalysis(), -1);
}

// Test 3: FLAG field matching, unmatched, forward matching, reverse matching
TEST_F(SamDecompressTest, TestFlagVariants) {
    generateSamFileWithFlagVariants("test_flag_matching.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_flag_matching.sam"));
}

// Test 4.1: Fixed-length Base field scenario
TEST_F(SamDecompressTest, TestBaseFieldFixed) {
    generateSamFileWithFixedBase("test_base_fixed.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_base_fixed.sam"));
}

// Test 4.2: Variable-length Base field scenario
TEST_F(SamDecompressTest, TestBaseFieldVariable) {
    generateSamFileWithVariableBase("test_base_variable.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_base_variable.sam"));
}

// Test 4.3: 3rd generation fastq scenario
TEST_F(SamDecompressTest, TestBaseFieldFastq3) {
    generateSamFileWithFastq3Base("test_base_fastq3.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_base_fastq3.sam"));
}

// Test 5.1: With optional fields
TEST_F(SamDecompressTest, TestOptionalFieldsPresent) {
    generateSamFileWithOptionalFields("test_optional_fields.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_optional_fields.sam"));
}

// Test 5.2: Without optional fields
TEST_F(SamDecompressTest, TestNoOptionalFields) {
    generateSamFileWithoutOptionalFields("test_no_optional_fields.sam");
    // Only the compression/decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_no_optional_fields.sam"));
}

/* Reads with missing quality (a single '*'): the compression side expands them to SEQ length and
 * the decompression side folds them back to '*', keeping the quality-value stream aligned with
 * record lengths. A '*' read must not disturb the rest of the QUAL column. */
TEST_F(SamDecompressTest, TestMissingQuality) {
    std::ofstream file("test_missing_quality.sam");
    ASSERT_TRUE(file.is_open());
    file << "@HD\tVN:1.6\n";
    file << "@SQ\tSN:chr1\tLN:1000\n";
    file << "m1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
    file << "u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t*\n";
    file << "u2\t4\t*\t0\t0\t*\t*\t0\t0\tTTTTTTTTTT\t*\tXA:B:i,1,2,3\n";
    file << "s1\t0\tchr1\t21\t60\t2S2M\t*\t0\t0\tACGT\t*\n";
    file.close();

    EXPECT_NO_THROW(compressAndDecompress("test_missing_quality.sam"));
    std::remove("test_missing_quality.sam");
}


// Comprehensive test: Mix all scenarios
TEST_F(SamDecompressTest, TestMixedScenarios) {
    std::ofstream file("test_mixed.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:340\n";
    file << "@SQ\tSN:chr2\tLN:338\n";

        // Mix various separators
    file << "read1:part1.part2#part3/part4|part5-part6_part7$part8\t65\tchr1\t1\t60\t100M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:5\tMD:Z:95A0\tAS:i:95\tXS:i:0\n";
        // Records with missing fields
    file << "read2\t4\t*\t0\t0\t*\t*\t0\t0\n";
        // Forward matching
    file << "read3\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        // Reverse matching
    file << "read4\t16\tchr2\t3\t60\t200M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG_ATCGATCG_ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:10\n";
    file.close();

    // Only the compression and decompression return codes are verified here; the file content is not compared
    // Load input file (contains invalid records with missing fields; load as raw bytes to verify the actuator's pre-analysis rejection logic)
    loadSamDataRaw("test_mixed.sam");

    // Create compressor
    Reference ref = createTestReference();
    PbgzParameter para;
    CompressEngine engine(para);

    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), -1);
}

// Compression performance test (large file)
TEST_F(SamDecompressTest, TestCompressionPerformance) {
    // Generate large SAM file for performance testing
    std::ofstream file("test_large.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:248956422\n";
    file << "@SQ\tSN:chr2\tLN:242193529\n";
    file << "@SQ\tSN:chr3\tLN:198295559\n";
    file << "@SQ\tSN:chr4\tLN:190214555\n";
    file << "@SQ\tSN:chr5\tLN:181538259\n";

    // Generate 10000 records
    for (int i = 0; i < 10000; i++) {
        file << "read_" << i << "_test:some/complex|id#$" << i % 10 << "\t";
        file << (i % 4 == 0 ? 4 : (i % 2 == 0 ? 0 : 16)) << "\t";
        file << "chr" << (i % 5 + 1) << "\t" << (i * 100 + 1) << "\t60\t";
        file << (50 + i % 100) << "M\t*\t0\t0\t";

        // Generate variable-length sequences
        int seqLen = 50 + i % 100;
        for (int j = 0; j < seqLen; j++) {
            file << "ATCG"[j % 4];
        }
        file << "\t";

        // Generate quality values
        for (int j = 0; j < seqLen; j++) {
            file << "!";
        }

        if (i % 3 == 0) {
            // Some records contain optional fields
            file << "\tNM:i:" << (i % 10) << "\tMD:Z:" << seqLen << "\tAS:i:" << (90 - i % 10);
        }
        file << "\n";
    }

    file.close();

    // Only the compression and decompression return codes are verified here; the file content is not compared
    EXPECT_NO_THROW(compressAndDecompress("test_large.sam"));
}

// Test 6.1: Field count too few scenario
TEST_F(SamDecompressTest, TestFieldCountTooFew) {
    std::ofstream file("test_field_count_too_few.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:340\n";
    file << "@SQ\tSN:chr2\tLN:338\n";

    // Normal record
    file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
    // Record with too few fields (only 8 fields, missing 3 required fields)
    file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\n";
    // Another normal record
    file << "read3\t0\tchr2\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
    file.close();

    // Load input file (contains records with missing fields; load as raw bytes to verify the actuator's pre-analysis rejection logic)
    loadSamDataRaw("test_field_count_too_few.sam");

    // Create compressor
    Reference ref = createTestReference();
    PbgzParameter para;
    CompressEngine engine(para);

    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    int32_t ret = compressor.preAnalysis();
    std::remove("test_field_count_too_few.sam");
    EXPECT_EQ(ret, -1);
}

// Test 6.2: Field missing scenario
TEST_F(SamDecompressTest, TestFieldMissing) {
    std::ofstream file("test_field_missing.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:340\n";
    file << "@SQ\tSN:chr2\tLN:338\n";

    // Normal record
    file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
    // Record with missing SEQ field (sequence)
    file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\t\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
    // Another normal record
    file << "read3\t0\tchr2\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
    file.close();

    // Load input file
    loadSamData("test_field_missing.sam");

    // Create compressor
    Reference ref = createTestReference();
    PbgzParameter para;
    CompressEngine engine(para);

    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    int32_t ret = compressor.preAnalysis();
    std::remove("test_field_missing.sam");
    EXPECT_EQ(ret, -1);
}

/*
 * A roundtrip that keeps the reference in play and compares the text: compressAndDecompress above
 * only asks whether the decompressor answered 0, which a reference applied to half the records - or
 * to none of them - can still do. These records sit at the reference's own coordinates, so their SEQ
 * column travels the reference-coded path (see sam_seq_payload), and what comes back has to be the
 * bytes that went in.
 */
/*
 * A file may carry no @SQ at all and yet be compressed with a reference supplied on the
 * command line. No record then uses it, so the file-level verdict is "reference unused"
 * and decompression skips loading it. Each block must still be self-contained: the SEQ
 * layout has to say "noref" even though a reference object was supplied to the encoder.
 */
TEST_F(SamDecompressTest, TestUnusedReferenceBlockDecodesWithoutReference) {
    const std::string file = "test_unused_ref_no_sq.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:unsorted\n";
        out << "@RG\tID:A\tSM:test\n";
        for (int r = 0; r < 12; ++r) {
            std::string seq;
            for (int i = 0; i < 50; ++i) {
                seq += "ACGT"[(r + i) % 4];
            }
            out << "read" << (r + 1) << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
                << seq << "\t" << std::string(50, 'I') << "\tRG:Z:A\n";
        }
    }

    loadSamData(file);
    /* Supplied, but no record has a position that could address it. */
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);
    ASSERT_FALSE(engine.isReferenceUsed());

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    /* The no-reference decode path taken when the archive says refe.required=false. */
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, nullptr);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

TEST_F(SamDecompressTest, TestReferenceRoundtripPreservesTheText) {
    const std::string file = "test_reference_roundtrip.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n";
        out << "@SQ\tSN:chr1\tLN:340\n";
        for (int r = 0; r < 12; ++r) {
            const int pos = 1 + r * 20;
            std::string seq;
            for (int i = 0; i < 50; ++i) {
                seq += "ATCG"[(pos + i) % 4];
            }
            out << "read" << (r + 1) << "\t0\tchr1\t" << pos << "\t60\t50M\t*\t0\t0\t"
                << seq << "\t" << std::string(50, 'I') << "\n";
        }
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected;
    {
        char buf[4096];
        while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
            expected.append(buf, (size_t)in.gcount());
        }
    }
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

/*
 * SEQ '*' means "unavailable" and may still carry a real CIGAR (PacBio secondary alignments).
 * Its field length is 1, so the CIGAR cannot tell the decoder how many bases to expect; the
 * record's length has to travel in the baselen stream. A following record must stay aligned.
 */
TEST_F(SamDecompressTest, TestMissingSeqWithCigarRoundtrip) {
    const std::string file = "test_missing_seq_cigar.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n";
        out << "@SQ\tSN:chr1\tLN:340\n";
        out << "m1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n";
        out << "star1\t256\tchr1\t21\t0\t10M5I2D8M\t*\t0\t0\t*\t*\n";
        out << "m2\t0\tchr1\t41\t60\t12M\t*\t0\t0\tACGTACGTACGT\tIIIIIIIIIIII\n";
        out << "star2\t256\tchr1\t61\t0\t3S20M\t*\t0\t0\t*\t*\n";
        out << "m3\t16\tchr1\t81\t60\t9M\t*\t0\t0\tTTTTAAAAC\tIIIIIIIII\n";
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected;
    {
        char buf[4096];
        while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
            expected.append(buf, (size_t)in.gcount());
        }
    }
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

/*
 * A QUAL column of two or more '*' is real Phred-9 data, not the single '*' missing marker.
 * The "qmiss" flag stream lets the decoder keep the run instead of folding it to one '*'.
 */
TEST_F(SamDecompressTest, TestAllStarQualityRoundtrip) {
    const std::string file = "test_all_star_quality.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n";
        out << "@SQ\tSN:chr1\tLN:340\n";
        out << "m1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        out << "star1\t0\tchr1\t5\t60\t2M\t*\t0\t0\tAC\t**\n";
        out << "star2\t4\t*\t0\t0\t*\t*\t0\t0\tACG\t***\n";
        out << "m2\t0\tchr1\t9\t60\t4M\t*\t0\t0\tACGT\t****\n";
        out << "miss\t0\tchr1\t13\t60\t3M\t*\t0\t0\tACG\t*\n";
        out << "m3\t16\tchr1\t17\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected;
    {
        char buf[4096];
        while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
            expected.append(buf, (size_t)in.gcount());
        }
    }
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

/*
 * QNAME layouts other than plain text (dictionary, fixed alphabet, numeric/counter) are what
 * the old, separate PNEXT pre-decoder could not reproduce: it dropped their separators and
 * skipped several layouts, so mate pairing by QNAME failed and PNEXT became 0 (TLEN -POS).
 * This block mixes all of them, checks that those layouts are really exercised, and then
 * compares the decoded text byte for byte (PNEXT/TLEN come from the rebuild and exceptions).
 */
TEST_F(SamDecompressTest, TestQnameLayoutsPnextRebuildRoundtrip) {
    const std::string file = "test_qname_pnext.sam";
    const std::string seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const std::string qual(36, 'I');
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n";
        out << "@SQ\tSN:chr1\tLN:2000000\n";
        static const char* inst[3] = {"INST0", "INST1", "INST2"};
        static const char alnum[] = "ABCD";
        for (int i = 0; i < 400; ++i) {
            std::string hex;
            uint32_t rnd = (uint32_t)i * 2654435761u + 12345u;
            for (int j = 0; j < 4; ++j) {
                rnd = rnd * 1664525u + 1013904223u;
                hex.push_back(alnum[(rnd >> 16) & 3u]);
            }
            char num[16];
            snprintf(num, sizeof(num), "%06d", i * 7);
            const std::string qname = std::string(inst[i % 3]) + ":" + hex + ":" + num +
                                      ":read" + std::to_string(i);
            const int p1 = 101 + i * 100;
            const int p2 = p1 + 50;
            out << qname << "\t99\tchr1\t" << p1 << "\t60\t36M\t=\t" << p2
                << "\t86\t" << seq << "\t" << qual << "\n";
            out << qname << "\t147\tchr1\t" << p2 << "\t60\t36M\t=\t" << p1
                << "\t-86\t" << seq << "\t" << qual << "\n";
        }
        /* No mate in the block: PNEXT has to travel in the exception stream. */
        out << "INST3:EXC:999999:orphan\t99\tchr1\t500001\t60\t36M\t=\t600001\t100036\t"
            << seq << "\t" << qual << "\n";
        /* Mate unmapped: PNEXT is 0 and no exception is written. */
        out << "INST3:EXC:999998:nomate\t4\t*\t0\t0\t*\t*\t0\t0\t"
            << seq << "\t" << qual << "\n";
        /* Non-reciprocal pair: both records have to come from the exception stream. */
        out << "INST3:EXC:999997:nonrecip\t99\tchr1\t700001\t60\t36M\t=\t700999\t1000\t"
            << seq << "\t" << qual << "\n";
        out << "INST3:EXC:999997:nonrecip\t147\tchr1\t700999\t60\t36M\t=\t700500\t1000\t"
            << seq << "\t" << qual << "\n";
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    /* The crafted names must really reach the dictionary, fixed-alphabet and value layouts;
       those are the paths the deleted pre-decoder used to skip or mangle. */
    bool sawDict = false;
    bool sawValue = false;
    Json::Value compressedMeta;
    coder_json metaJson;
    metaJson.decoder(pOutBlock->getMetaBuffer(), (int32_t)pOutBlock->getMetaLen(), compressedMeta);
    const Json::Value& qstreams = compressedMeta["sam"]["streams"][0]["streams"];
    std::string modes;
    for (Json::ArrayIndex i = 0; i < qstreams.size(); ++i) {
        const std::string mode = qstreams[i].isMember("mode")
                                     ? qstreams[i]["mode"].asString() : std::string();
        modes += (i ? "," : "") + (mode.empty() ? std::string("text") : mode);
        sawDict = sawDict || mode == "dict";
        sawValue = sawValue || mode == "numeric" || mode == "cnt" || mode == "intd" ||
                   mode == "intu";
    }
    /* Dictionary and value-domain layouts are the ones the old pre-decoder mangled for
       ctvt (missing separators); the same unified decode also covers the fixed-alphabet
       layout, which this encoder profile rarely selects over the text coder. */
    EXPECT_TRUE(sawDict) << "modes=" << modes;
    EXPECT_TRUE(sawValue) << "modes=" << modes;

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected;
    {
        char buf[4096];
        while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
            expected.append(buf, (size_t)in.gcount());
        }
    }
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

/*
 * The qmiss stream is written only when a block needs it (a real all-'*' column or a truly
 * missing QUAL). Blocks without either keep the two-stream QUAL layout, so their bytes and
 * their compatibility with older decoders are unchanged.
 */
TEST_F(SamDecompressTest, TestQualityMissingStreamOnlyWhenNeeded) {
    auto qualStreamCount = [&](const std::string& name, const std::string& qual) -> size_t {
        {
            std::ofstream out(name);
            EXPECT_TRUE(out.is_open());
            out << "@HD\tVN:1.6\n";
            out << "@SQ\tSN:chr1\tLN:1000\n";
            out << "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\t" << qual << "\n";
        }
        loadSamData(name);
        Reference ref = createTestReference();
        PbgzParameter para;
        DecompressEngine engine(para);
        SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
        EXPECT_EQ(compressor.preAnalysis(), 0);
        EXPECT_EQ(compressor.compress(), 0);
        Json::Value compressedMeta;
        coder_json metaJson;
        metaJson.decoder(pOutBlock->getMetaBuffer(), (int32_t)pOutBlock->getMetaLen(),
                         compressedMeta);
        std::remove(name.c_str());
        return compressedMeta["sam"]["streams"][10]["streams"].size();
    };

    EXPECT_EQ(qualStreamCount("q_plain.sam", "IIII"), 2u);   /* no marker needed */
    EXPECT_EQ(qualStreamCount("q_stars.sam", "**"), 3u);     /* real all-'*' column */
    EXPECT_EQ(qualStreamCount("q_missing.sam", "*"), 3u);    /* missing QUAL */
}

/*
 * QNAMEs whose separator shape changes from line to line (real RNA-seq/ONT names do this)
 * cannot use the split layout. The whole-column fallback has to be taken without indexing a
 * short per-line position vector, which used to crash the compressor.
 */
TEST_F(SamDecompressTest, TestVariableShapeQnamesRoundtrip) {
    const std::string file = "test_variable_qnames.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n";
        out << "@SQ\tSN:chr1\tLN:2000000\n";
        for (int i = 0; i < 200; ++i) {
            std::string qname;
            if (i % 2 == 0) {
                qname = "Other_RNA_exon_exon_ENSRNOG00000032748|Abcc1|10|11261142|11261318///"
                        "MN01080:7:000H37TCT:1:21103:19852:";
            } else {
                qname = "Genome///MN01080:7:000H37TCT:1:11102:17240:";
            }
            qname += std::to_string(i);
            out << qname << "\t0\tchr1\t" << (1 + i * 10) << "\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        }
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);
    EXPECT_TRUE(compressor.idAnalysis.unavailable);

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected;
    {
        char buf[4096];
        while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
            expected.append(buf, (size_t)in.gcount());
        }
    }
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}

namespace {

/* The main SEQ payload substream is named "m"; its coder is the one this test pins. */
std::string seqPayloadCoderMagic(RoughIOBlock* block)
{
    Json::Value meta;
    coder_json json;
    json.decoder(block->getMetaBuffer(), (int32_t)block->getMetaLen(), meta);
    const Json::Value& streams = meta["sam"]["streams"][9]["streams"];
    std::string fallback;
    uint64_t best = 0;
    for (Json::ArrayIndex i = 0; i < streams.size(); ++i) {
        const Json::Value& s = streams[i];
        if (!s.isMember("coder") || !s["coder"].isMember("magic")) {
            continue;
        }
        const std::string magic = s["coder"]["magic"].asString();
        const std::string name = s.isMember("sname") ? s["sname"].asString() : std::string();
        if (name == "m") {
            return magic;
        }
        const uint64_t dst = s.isMember("dstlen") ? s["dstlen"].asUInt64() : 0;
        if (dst >= best) {
            best = dst;
            fallback = magic;
        }
    }
    return fallback;
}

void writeAcgtSam(const std::string& file, uint32_t reads, uint32_t len, bool same, uint32_t seed)
{
    std::ofstream out(file);
    ASSERT_TRUE(out.is_open());
    out << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n";
    std::mt19937 rng(seed);
    auto nextBase = [&rng]() { return "ACGT"[rng() & 3u]; };
    std::string first;
    first.reserve(len);
    for (uint32_t i = 0; i < len; ++i) {
        first.push_back(nextBase());
    }
    const std::string qual(len, 'I');
    for (uint32_t i = 0; i < reads; ++i) {
        std::string seq;
        if (same) {
            seq = first;
        } else {
            seq.reserve(len);
            for (uint32_t j = 0; j < len; ++j) {
                seq.push_back(nextBase());
            }
        }
        out << 'r' << i << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << seq << '\t' << qual << '\n';
    }
}

std::string readTextFile(const std::string& file)
{
    std::ifstream in(file, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}

}  // namespace

/*
 * Long all-ACGT reads at -l8 must not take the packed order-1 fast path. It replaced
 * the selected coder with coder_arith_k, which cannot see the cross-read redundancy
 * of long-read data (SRR12223958 was 1.45MB instead of 127KB). The unit path has no
 * engine pre-selection, so either coder_fc or the coder_bwt_cm fallback is valid.
 */
TEST_F(SamDecompressTest, TestLongReadAllAcgtSeqAvoidsArithKAtLevel8)
{
    const std::string file = "test_long_read_fc.sam";
    writeAcgtSam(file, 200, 1000, true, 11);
    const std::string expected = readTextFile(file);

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    para.compressLevel = 8;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);
    const std::string coder = seqPayloadCoderMagic(pOutBlock);
    EXPECT_NE(coder, "coder_arith_k");
    EXPECT_TRUE(coder == "coder_fc" || coder == "coder_bwt_cm") << "coder=" << coder;

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);
    EXPECT_EQ(std::string((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen()),
              expected);
    std::remove(file.c_str());
}

/*
 * Short all-ACGT reads are handled by the data-driven SEQ coder selection at -l8.
 * This synthetic file is 4000 copies of one read, so the redundant-data rule may pick
 * coder_fc (or the engine fallback); only a SEQ coders set and a clean roundtrip are
 * required here. The order-1 speed guard is exercised by the large-file benchmark.
 */
TEST_F(SamDecompressTest, TestShortReadAllAcgtSeqUsesSeqCoderAtLevel8)
{
    const std::string file = "test_short_read_arithk.sam";
    writeAcgtSam(file, 4000, 100, true, 22);
    const std::string expected = readTextFile(file);

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    para.compressLevel = 8;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);
    EXPECT_TRUE(seqPayloadCoderMagic(pOutBlock) == "coder_fc" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_bwt_cm" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_arith_ctx" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_arith_k");

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();
    SamCodecActuator decompressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);
    EXPECT_EQ(std::string((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen()),
              expected);
    std::remove(file.c_str());
}

/*
 * Long random all-ACGT reads are the high-entropy case. The data-driven selection
 * compares fc/arith_ctx/order-1 on the actual match payload; a small unit file may
 * legitimately pick any of the SEQ coders, so only the SEQ coder set and a clean
 * roundtrip are asserted. The large-file speed guard is covered by the benchmark.
 */
TEST_F(SamDecompressTest, TestLongRandomAcgtSeqUsesSeqCoderAtLevel8)
{
    const std::string file = "test_long_random_arithk.sam";
    writeAcgtSam(file, 3000, 1000, false, 33);

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    para.compressLevel = 8;
    DecompressEngine engine(para);
    SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);
    EXPECT_TRUE(seqPayloadCoderMagic(pOutBlock) == "coder_fc" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_bwt_cm" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_arith_ctx" ||
                seqPayloadCoderMagic(pOutBlock) == "coder_arith_k");
    std::remove(file.c_str());
}

/*
 * The archive records whether the supplied reference was actually used. A mapped read
 * with a matching @SQ marks it; an unmapped read does not. Decompression uses the
 * same flag to decide whether -r is mandatory for a non-packed reference.
 */
TEST_F(SamDecompressTest, TestReferenceUsedFlagTracksActualUse) {
    /* A 1000 bp chr1 whose @SQ length matches the reference, so the actuator really
       takes the reference path for the mapped record. */
    const std::string refFile = "test_ref_used.fa";
    {
        std::ofstream out(refFile);
        ASSERT_TRUE(out.is_open());
        out << ">chr1\n";
        for (int i = 0; i < 250; ++i) out << "ACGT";
        out << "\n";
    }
    Reference ref(refFile, 1);
    ref.makeIndex();

    const std::string mapped = "test_ref_used_mapped.sam";
    {
        std::ofstream out(mapped);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
        out << "r1\t0\tchr1\t1\t60\t100M\t*\t0\t0\t"
            << std::string("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"
                           "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")
            << "\t" << std::string(100, 'I') << "\n";
    }
    {
        loadSamData(mapped);
        PbgzParameter para;
        DecompressEngine engine(para);
        SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
        ASSERT_EQ(compressor.preAnalysis(), 0);
        ASSERT_EQ(compressor.compress(), 0);
        EXPECT_TRUE(engine.isReferenceUsed());
    }

    const std::string unmapped = "test_ref_unused_unmapped.sam";
    {
        std::ofstream out(unmapped);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n";
        out << "r1\t4\t*\t0\t0\t*\t*\t0\t0\t"
            << std::string("ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC"
                           "ACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC")
            << "\t" << std::string(100, 'I') << "\n";
    }
    {
        loadSamData(unmapped);
        PbgzParameter para;
        DecompressEngine engine(para);
        SamCodecActuator compressor(pInBlock, pOutBlock, &engine, &ref);
        ASSERT_EQ(compressor.preAnalysis(), 0);
        ASSERT_EQ(compressor.compress(), 0);
        EXPECT_FALSE(engine.isReferenceUsed());
    }

    std::remove(mapped.c_str());
    std::remove(unmapped.c_str());
    std::remove(refFile.c_str());
    std::remove((refFile + ".fai").c_str());
    std::remove((refFile + ".ni").c_str());
}

/*
 * The QNAME dedup layout (affix on first occurrences + a separate repeat stream) is a
 * file-level verdict from the codec pre-selection. Pin it, then roundtrip a block whose
 * names repeat at distances that cross many lines, with mapped and unmapped records
 * interleaved in an unsorted order; the decoded text must be byte-identical.
 */
TEST_F(SamDecompressTest, TestQnameDedupLayoutRoundtrip) {
    const std::string file = "test_qname_dedup.sam";
    {
        std::ofstream out(file);
        ASSERT_TRUE(out.is_open());
        out << "@HD\tVN:1.6\tSO:unsorted\n";
        out << "@SQ\tSN:chr1\tLN:2000000\n";
        std::vector<std::string> pool;
        char name[128];
        for (int i = 0; i < 4000; ++i) {
            snprintf(name, sizeof(name), "HWUSI-EAS492_0029_FC:8:%d:%d:%d#0",
                     1 + i % 120, 1000 + (i * 3) % 19000, 1000 + (i * 7) % 19000);
            pool.emplace_back(name);
        }
        uint32_t rnd = 2463534242u;
        for (int i = 0; i < 20000; ++i) {
            rnd = rnd * 1664525u + 1013904223u;
            const std::string& qname = pool[(rnd >> 9) % pool.size()];
            rnd = rnd * 1664525u + 1013904223u;
            if ((rnd & 1u) != 0) {
                out << qname << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
            } else {
                const uint32_t flag = (rnd & 2u) ? 0u : 16u;
                out << qname << '\t' << flag << "\tchr1\t" << (1 + (rnd >> 2) % 1900000u)
                    << "\t60\t*\t*\t0\t0\t*\t*\n";
            }
        }
    }

    loadSamData(file);
    Reference ref = createTestReference();
    PbgzParameter para;
    CompressEngine ceng(para);
    ASSERT_NE(ceng.getPreprocessInfoMut(), nullptr);
    ceng.getPreprocessInfoMut()->setQnameUseAffix(2);
    SamCodecActuator compressor(pInBlock, pOutBlock, &ceng, &ref);
    ASSERT_EQ(compressor.preAnalysis(), 0);
    ASSERT_EQ(compressor.compress(), 0);

    const int64_t compressedLen = pOutBlock->getDataLen();
    const int64_t metaLen = pOutBlock->getMetaLen();
    pInBlock->reset();
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), (size_t)(compressedLen + metaLen));
    pInBlock->setDataLen(compressedLen);
    pInBlock->setMetaLen(metaLen);
    pOutBlock->reset();

    DecompressEngine deng(para);
    SamCodecActuator decompressor(pInBlock, pOutBlock, &deng, &ref);
    ASSERT_EQ(decompressor.decompress(), 0);

    std::ifstream in(file, std::ios::binary);
    std::string expected((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string got((const char*)pOutBlock->getBuffer(), (size_t)pOutBlock->getDataLen());
    EXPECT_EQ(got.size(), expected.size());
    EXPECT_EQ(got, expected);
    std::remove(file.c_str());
}
