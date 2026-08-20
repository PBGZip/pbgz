/*
 * sam_sort_combine_testcase.cpp - Test cases for SAM sorting and merging functionality
 */

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <thread>
#include <mutex>
#include <gmock/gmock.h>

#include "coder/coder.h"
#include "blocking_queue.h"
#include "io_block.h"
#include "utils/memory_util.h"
#include "io_wrapper.h"

#include "pbgz_testcase_util.h"

namespace SamSortTestData {
    const std::string sortedHeadFile = "sorted_head.sam";
    const std::string sortedSamFile0 = "sorted_sam_0_0.sam";
    const std::string sortedSamFile1 = "sorted_sam_0_1.sam";
    const std::string sortedSamFile2 = "sorted_sam_0_2.sam";
    const std::string mergedOutputFile = "merged_output.sam";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};


#define protected public
#define private public
#include "pbgz_types.h"
#include "sort_engine.h"
#include "sam_sort_actuator.h"
#include "config_manager.h"
#include "sam_info.h"
#include "pbgz_manager.h"
#undef private
#undef protected


// MockSortEngine - Shadow SortEngine's queues with MockBlockingQueue
class MockSortEngine : public SortEngine {
public:
    MockSortEngine(const PbgzParameter& para, uint32_t blockSz = 8 << 20)
        : SortEngine(para) {
        freeOutputPool = std::make_unique<MockBlockingQueue>(blockSz);
        outputDataPool = std::make_unique<MockBlockingQueue>(blockSz);
    }

    int32_t init(const std::string& outputFile) {
        ioWriter = MemoryUtil::safeNewClass<FileWriter>(outputFile);
        ioWriter->openIO();
        ((MockBlockingQueue*)freeOutputPool.get())->setIOWriter(ioWriter);
        ((MockBlockingQueue*)outputDataPool.get())->setIOWriter(ioWriter);
        return 0;
    }

    void finish() {
        ioWriter->closeIO();
    }

    virtual ~MockSortEngine(){
        MemoryUtil::safeDeleteClass(ioWriter);
    }
};

class SamSortCombineTest : public ::testing::Test {
public:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

        // Initialize chromosome info before each test
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();

        ConfigManager::getInstance().logLevel = LogLevel::WARNING;

        // Clean up any existing output files
        std::remove(SamSortTestData::sortedHeadFile.c_str());
        std::remove(SamSortTestData::sortedSamFile0.c_str());
        std::remove(SamSortTestData::sortedSamFile1.c_str());
        std::remove(SamSortTestData::sortedSamFile2.c_str());
        std::remove(SamSortTestData::mergedOutputFile.c_str());

        std::remove("merge_output.sam");
        std::remove("final_output.bam");
        std::remove("sorted_head.sam");
        std::remove("sorted_sam_0_0.sam");
        std::remove("sorted_sam_0_1.sam");
        std::remove("sorted_sam_0_2.sam");
        std::remove("merged_output.sam");
    }

    void TearDown() override {
        // Clean up all generated test files
        std::remove(SamSortTestData::sortedHeadFile.c_str());
        std::remove(SamSortTestData::sortedSamFile0.c_str());
        std::remove(SamSortTestData::sortedSamFile1.c_str());
        std::remove(SamSortTestData::sortedSamFile2.c_str());
        std::remove(SamSortTestData::mergedOutputFile.c_str());
        std::remove("merge_output.sam");
        std::remove("final_output.bam");
        std::remove("sorted_head.sam");
        std::remove("sorted_sam_0_0.sam");
        std::remove("sorted_sam_0_1.sam");
        std::remove("sorted_sam_0_2.sam");
        std::remove("merged_output.sam");
        std::remove("test_merged.sam");

    }

    void createTestSortedHeadFile(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to create test header file");
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:248956422\n";
        file << "@SQ\tSN:chr2\tLN:242193529\n";
        file << "@PG\tID:combine_test\tPN:pbgz_test\tVN:1.0\n";

        file.close();
    }

    void createTestSortedSamFile(const std::string& filename, const std::vector<std::string>& lines) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to create test SAM file");
        }

        for (const auto& line : lines) {
            file << line << "\n";
        }

        file.close();
    }

    void setupSortedSamFiles() {
        createTestSortedHeadFile(SamSortTestData::sortedHeadFile);

        createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
            "100:read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATATCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "300:read2\t0\tchr1\t300\t60\t76M\t*\t0\t0\tACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        });

        createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
            "200:read3\t0\tchr1\t200\t60\t76M\t*\t0\t0\tAGCTAGCTAGCTAGCTAGCTAGCT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "248956472:read4\t0\tchr2\t50\t60\t76M\t*\t0\t0\tTCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        });

        createTestSortedSamFile(SamSortTestData::sortedSamFile2, {
            "248956622:read5\t0\tchr2\t200\t60\t76M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
            "248956822:read6\t0\tchr2\t400\t60\t76M\t*\t0\t0\tCTAGCTAGCTAGCTAGCTAGCTAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        });

        SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
        SamInfo::getInstance().addChromosomeInfo("chr2", 242193529);
        SamInfo::getInstance().calculateChromosomePositions();
    }

    void verifyOutputContent(const std::string& filename, int expectedHeaderLines, int expectedRecords) {
        std::ifstream file(filename);
        ASSERT_TRUE(file.good()) << "Output file does not exist: " << filename;

        std::string line;
        int headerCount = 0;
        int recordCount = 0;

        while (std::getline(file, line)) {
            if (line.find("@HD") == 0 || line.find("@SQ") == 0 || line.find("@PG") == 0 || line.find("@CO") == 0) {
                headerCount++;
            } else if (!line.empty() && line[0] != '@') {
                recordCount++;
            }
        }

        file.close();

        EXPECT_EQ(headerCount, expectedHeaderLines) << "Expected " << expectedHeaderLines << " header lines, got " << headerCount;
        EXPECT_EQ(recordCount, expectedRecords) << "Expected " << expectedRecords << " records, got " << recordCount;
    }
};

TEST_F(SamSortCombineTest, SamSortCombineFile) {
    setupSortedSamFiles();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 3;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 6);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortSingleFile) {
    // Test the merge scenario for a single file
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "100:read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATATCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "200:read2\t0\tchr1\t200\t60\t76M\t*\t0\t0\tACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 1;  // only one file needs merging
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 2);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortWithUnmappedRecords) {
    // Test the scenario with unmapped (referencePos < 0) records
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "100:read1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATATCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:read_unmapped\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDSEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "200:read2\t0\tchr1\t200\t60\t76M\t*\t0\t0\tACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 2;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 3);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortSameChromosomeDifferentPositions) {
    // Test sorting of different positions on the same chromosome
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "500:read_pos500\t0\tchr1\t500\t60\t76M\t*\t0\t0\tPOS500SEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "100:read_pos100\t0\tchr1\t100\t60\t76M\t*\t0\t0\tPOS100SEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile2, {
        "300:read_pos300\t0\tchr1\t300\t60\t76M\t*\t0\t0\tPOS300SEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 3;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 3);

    // Verify the output order is correct (ascending by position)
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    int recordIndex = 0;
    std::vector<std::string> expectedReads = {"read_pos100", "read_pos300", "read_pos500"};
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            EXPECT_TRUE(line.find(expectedReads[recordIndex]) == 0)
                << "Expected line " << recordIndex + 1 << " to start with " << expectedReads[recordIndex];
            recordIndex++;
        }
    }
    file.close();
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortMultipleChromosomes) {
    // Test sorting records across chromosomes
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "248956522:read_chr2\t0\tchr2\t100\t60\t76M\t*\t0\t0\tCHR2POS100SEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "50:read_chr1\t0\tchr1\t50\t60\t76M\t*\t0\t0\tCHR1POS50SEQUENCE\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info (chr2's position should come after chr1)
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().addChromosomeInfo("chr2", 242193529);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 1;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 2);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortNoInputFiles) {
    // Test the scenario with no input files (blockCount = 0)
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().addChromosomeInfo("chr2", 242193529);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 0;  // no input files to process
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 0);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortLargeVolumeRecords) {
    // Test a scenario with many records to verify the chunking behavior
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);

    std::vector<std::string> records0, records1;
    for (int i = 0; i < 20; i++) {
        records0.push_back(std::to_string(i * 10) + ":read_file1_pos" + std::to_string(i * 10) + "\t0\tchr1\t" + std::to_string(i * 10) + "\t60\t76M\t*\t0\t0\tSEQUENCE" + std::to_string(i) + "\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        records1.push_back(std::to_string(i * 10 + 5) + ":read_file2_pos" + std::to_string(i * 10 + 5) + "\t0\tchr1\t" + std::to_string(i * 10 + 5) + "\t60\t76M\t*\t0\t0\tSEQUENCE" + std::to_string(i) + "\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    }
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, records0);
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, records1);

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 2;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 40);
    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamsortMultipleUnmappedInSingleFile) {
    // Test the scenario where a single file contains multiple unmapped records
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "-1:unmapped1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDSEQUENCE1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDSEQUENCE2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 1;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);

    // Verify the output contains 2 unmapped records
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    int unmappedCount = 0;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            if (line.find("unmapped") == 0) {
                unmappedCount++;
            }
        }
    }
    file.close();
    EXPECT_EQ(unmappedCount, 2);

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortMultipleFilesWithMultipleUnmapped) {
    // Test the scenario where multiple files each contain multiple unmapped records
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "-1:unmapped_f1_1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF1SEQ1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:mapped_f1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tMAPPEDF1SEQ\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped_f1_2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF1SEQ2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "-1:unmapped_f2_1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF2SEQ1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "200:mapped_f2\t0\tchr1\t200\t60\t76M\t*\t0\t0\tMAPPEDF2SEQ\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped_f2_2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF2SEQ2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile2, {
        "-1:unmapped_f3_1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF3SEQ1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "300:mapped_f3\t0\tchr1\t300\t60\t76M\t*\t0\t0\tMAPPEDF3SEQ\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped_f3_2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPEDF3SEQ2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().addChromosomeInfo("chr2", 242193529);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 3;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);

    // Verify output: 3 mapped + 6 unmapped = 9 records
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 9);

    // Verify unmapped records are handled correctly
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    int unmappedCount = 0;
    int mappedCount = 0;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            if (line.find("unmapped") == 0) {
                unmappedCount++;
            } else if (line.find("mapped") == 0) {
                mappedCount++;
            }
        }
    }
    file.close();
    EXPECT_EQ(unmappedCount, 6);
    EXPECT_EQ(mappedCount, 3);

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortMixedMappedAndUnmappedInSequence) {
    // Test the scenario where mapped and unmapped records alternate within a single file
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "-1:unmapped_1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:mapped_100\t0\tchr1\t100\t60\t76M\t*\t0\t0\tMAPPED100\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped_2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "200:mapped_200\t0\tchr1\t200\t60\t76M\t*\t0\t0\tMAPPED200\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped_3\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED3\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 1;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);

    // Verify the output
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    int unmappedCount = 0;
    int mappedCount = 0;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            if (line.find("unmapped") == 0) {
                unmappedCount++;
            } else if (line.find("mapped") == 0) {
                mappedCount++;
            }
        }
    }
    file.close();
    EXPECT_EQ(mappedCount, 2);
    EXPECT_EQ(unmappedCount, 3);

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortEmptyHeaderFile) {
    // Test the scenario where the first line read back is an unmapped record
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);

    // File 0's first line is an unmapped record
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "-1:first_line_unmapped\t4\t*\t0\t0\t76M\t*\t0\t0\tFIRSTUNMAPPED\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:mapped_after\t0\tchr1\t100\t60\t76M\t*\t0\t0\tMAPPEDAFTER\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "50:second_file_mapped\t0\tchr1\t50\t60\t76M\t*\t0\t0\tSECONDMAPPED\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 2;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);

    // Verify the output contains 1 unmapped and 2 mapped records
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    int recordCount = 0;
    int unmappedCount = 0;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            recordCount++;
            if (line.find("first_line_unmapped") == 0) {
                unmappedCount++;
            }
        }
    }
    file.close();
    EXPECT_EQ(recordCount, 3);
    EXPECT_EQ(unmappedCount, 1);

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortAllUnmappedRecords) {
    // Test the scenario where all records are unmapped
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "-1:unmapped1\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED1\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped2\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED2\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "-1:unmapped3\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED3\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "-1:unmapped4\t4\t*\t0\t0\t76M\t*\t0\t0\tUNMAPPED4\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 2;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);

    // Verify all 4 records are unmapped
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 4);

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortSamePositionDifferentIdsInSingleFile) {
    // Test records with different IDs at the same position within a single file
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "100:read_id1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATATCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:read_id2\t0\tchr1\t100\t60\t76M\t*\t0\t0\tACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:read_id3\t0\tchr1\t100\t60\t76M\t*\t0\t0\tAGCTAGCTAGCTAGCTAGCTAGCT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "200:read_id4\t0\tchr1\t200\t60\t76M\t*\t0\t0\tTCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 1;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 4);

    // Verify the output contains all four records in their original order
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    std::vector<std::string> foundIds;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            size_t tabPos = line.find('\t');
            if (tabPos != std::string::npos) {
                std::string readName = line.substr(0, tabPos);
                if (readName.find("read_id") == 0) {
                    foundIds.push_back(readName);
                }
            }
        }
    }
    file.close();
    EXPECT_EQ(foundIds.size(), 4);
    // Verify order: the records at position 100 keep their original order, followed by the record at position 200
    // Since they share the same position (100), they follow the file's original order: read_id1, read_id2, read_id3, read_id4
    std::vector<std::string> expectedIds = {"read_id1", "read_id2", "read_id3", "read_id4"};
    for (size_t i = 0; i < expectedIds.size(); ++i) {
        EXPECT_EQ(foundIds[i], expectedIds[i])
            << "Expected " << expectedIds[i] << " at position " << i << ", got " << foundIds[i];
    }

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, SamSortSamePositionDifferentIdsInMultipleFiles) {
    // Test records with different IDs at the same position across multiple files
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);
    createTestSortedSamFile(SamSortTestData::sortedSamFile0, {
        "100:read_f1_id1\t0\tchr1\t100\t60\t76M\t*\t0\t0\tATATCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:read_f1_id2\t0\tchr1\t100\t60\t76M\t*\t0\t0\tACGTACGTACGTACGTACGTACGT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile1, {
        "100:read_f2_id3\t0\tchr1\t100\t60\t76M\t*\t0\t0\tAGCTAGCTAGCTAGCTAGCTAGCT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:read_f2_id4\t0\tchr1\t100\t60\t76M\t*\t0\t0\tTCGATCGATCGATCGATCGATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });
    createTestSortedSamFile(SamSortTestData::sortedSamFile2, {
        "200:read_f3_id5\t0\tchr1\t200\t60\t76M\t*\t0\t0\tGCTAGCTAGCTAGCTAGCTAGCT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
        "100:read_f3_id6\t0\tchr1\t100\t60\t76M\t*\t0\t0\tCTAGCTAGCTAGCTAGCTAGCTAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    });

    // Initialize chromosome info
    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);
    engine->blockCount = 3;
    int32_t result = engine->startEnginePostProc();
    engine->finish();
    EXPECT_EQ(result, 0);
    verifyOutputContent(SamSortTestData::mergedOutputFile, 4, 6);

    // Verify the output contains all six records in their original order
    std::ifstream file(SamSortTestData::mergedOutputFile);
    std::string line;
    std::vector<std::string> foundIds;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@') {
            size_t tabPos = line.find('\t');
            if (tabPos != std::string::npos) {
                std::string readName = line.substr(0, tabPos);
                if (readName.find("read_") == 0) {
                    foundIds.push_back(readName);
                }
            }
        }
    }
    file.close();
    EXPECT_EQ(foundIds.size(), 6);
    // Verify order: the five records at position 100 follow the file read order, then the record at position 200
    // File read order: file0 -> file1 -> file2
    std::vector<std::string> expectedIds = {"read_f1_id1", "read_f1_id2", "read_f2_id3", "read_f2_id4", "read_f3_id5", "read_f3_id6"};
    for (size_t i = 0; i < expectedIds.size(); ++i) {
        EXPECT_EQ(foundIds[i], expectedIds[i])
            << "Expected " << expectedIds[i] << " at position " << i << ", got " << foundIds[i];
    }

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, CombineSamFileWithFileWriter) {
    setupSortedSamFiles();

    std::string outputFile = "test_merged.sam";

    std::remove(outputFile.c_str());

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);

    std::vector<std::string> fileList = {
        SamSortTestData::sortedSamFile0,
        SamSortTestData::sortedSamFile1,
        SamSortTestData::sortedSamFile2
    };

    SamCombineOutputFileWriter fileWriter(outputFile);
    int32_t result = engine->combineSamFile(fileList, &fileWriter);
    fileWriter.close();

    EXPECT_EQ(result, 0);
    verifyOutputContent(outputFile, 0, 6);

    // Verify the merged order is correct
    std::ifstream file(outputFile);
    std::string line;
    std::vector<std::string> sortedLines;
    while (std::getline(file, line)) {
        if (!line.empty() && line[0] != '@' && line[0] != '\n') {
            sortedLines.push_back(line);
            LOG_DEBUG("%s, %d", line.c_str(), line.length());
        }
    }
    file.close();

    std::vector<std::string> expectedOrder = {
        "100:read1", "200:read3", "300:read2",
        "248956472:read4", "248956622:read5", "248956822:read6"
    };
    ASSERT_EQ(sortedLines.size(), expectedOrder.size());
    for (size_t i = 0; i < sortedLines.size(); ++i) {
        EXPECT_TRUE(sortedLines[i].find(expectedOrder[i]) == 0)
            << "Expected line " << i << " to start with " << expectedOrder[i];
    }

    // Clean up generated test file
    std::remove(outputFile.c_str());

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, CombineSamFileWithBlockWriter) {
    setupSortedSamFiles();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);
    engine->init(SamSortTestData::mergedOutputFile);

    SamCombineOutputBlockWriter blockWriter(engine->freeOutputPool.get(), engine->outputDataPool.get());
    blockWriter.initial(0);

    std::vector<std::string> fileList = {
        SamSortTestData::sortedSamFile0,
        SamSortTestData::sortedSamFile1,
        SamSortTestData::sortedSamFile2
    };
    int32_t result = engine->combineSamFile(fileList, &blockWriter);
    EXPECT_EQ(result, 0);
    blockWriter.close();
    engine->finish();

    MemoryUtil::safeDeleteClass(engine);
}

TEST_F(SamSortCombineTest, CombineAllSamFileExceeds128Files) {
    createTestSortedHeadFile(SamSortTestData::sortedHeadFile);

    const uint32_t fileCount = 130;

    for (uint32_t i = 0; i < fileCount; ++i) {
        std::string filename = "sorted_sam_0_" + std::to_string(i) + ".sam";
        std::vector<std::string> lines;
        lines.push_back(std::to_string(i * 100) + ":read_" + std::to_string(i) + "\t0\tchr1\t" + std::to_string(i * 100) + "\t60\t76M\t*\t0\t0\tSEQUENCE" + std::to_string(i) + "\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        createTestSortedSamFile(filename, lines);
    }

    SamInfo::getInstance().addChromosomeInfo("chr1", 248956422);
    SamInfo::getInstance().calculateChromosomePositions();

    PbgzParameter para;
    MockSortEngine* engine = MemoryUtil::safeNewClass<MockSortEngine>(para);

    std::vector<std::string> outputFiles;
    int32_t result = engine->combineAllSamFile(0, fileCount, outputFiles);

    EXPECT_EQ(result, 0);
    EXPECT_EQ(outputFiles.size(), 2);

    PbgzParameter para2;
    MockSortEngine* engine2 = MemoryUtil::safeNewClass<MockSortEngine>(para2);
    engine2->init(SamSortTestData::mergedOutputFile);

    SamCombineOutputBlockWriter blockWriter(engine2->freeOutputPool.get(), engine2->outputDataPool.get());
    blockWriter.initial(outputFiles.size());
    result = engine2->combineSamFile(outputFiles, &blockWriter);
    EXPECT_EQ(result, 0);
    blockWriter.close();
    engine2->finish();

    int recordCount = 0;
    std::vector<std::string> foundPositions;
     std::string line;
    std::ifstream mergeFile(SamSortTestData::mergedOutputFile);
    while (std::getline(mergeFile, line)) {
        if (!line.empty() && line[0] != '@' && line[0] != '\n') {
            recordCount++;
            size_t pos = line.find('\t');
            if (pos != std::string::npos) {
                foundPositions.push_back(line.substr(0, pos));
            }
        }
    }
    mergeFile.close();

    EXPECT_EQ(recordCount, fileCount) << "Should have all " << fileCount << " records in merged output";

    std::vector<std::string> expectedPositions;
    for (uint32_t i = 0; i < fileCount; ++i) {
        expectedPositions.push_back("read_" + std::to_string(i));
    }

    EXPECT_EQ(foundPositions, expectedPositions) << "Should have all positions in correct order";

    // Clean up all generated intermediate files
    for (const auto& file : outputFiles) {
        std::remove(file.c_str());
    }

    // Clean up all temporary SAM files created for this test
    for (uint32_t i = 0; i < fileCount; ++i) {
        std::string filename = "sorted_sam_0_" + std::to_string(i) + ".sam";
        std::remove(filename.c_str());
    }

    MemoryUtil::safeDeleteClass(engine);
    MemoryUtil::safeDeleteClass(engine2);
}