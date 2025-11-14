#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <memory>
#include <chrono>
#include <thread>
#include <vector>
// #define private public 
#include "reference.h"
//#undef private
#include "pbgz_errno.h"
#include "utils/path_util.h"
#include "coder.h"
#include "utils/memory_util.h"

class ReferenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create temporary FASTA file for testing
        testFastaFile = "test_reference.fasta";
        std::ofstream fastaFile(testFastaFile);
        fastaFile << ">test_sequence\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
        fastaFile.close();
        

        quashLen = 31;
        quash = MemoryUtil::safeAlloc<uint8_t>(quashLen); 
        for (int i = 0; i < 16; ++i) {
            quash[i] = 0x27;
        }
         for (int i = 16; i < 31; ++i) {
            quash[i] = 0xD8;
        }

        // Create test directory
        testDir = "test_reference_dir";
        std::filesystem::create_directory(testDir);

        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove(testFastaFile);
        std::filesystem::remove_all(testDir);
    }

    void printBufferHex(const uint8_t* buffer, uint32_t length) {
        for (uint32_t i = 0; i < length; ++i) {
            fprintf(stderr, "%X", buffer[i]);
        }
        fprintf(stderr, "\n");
    }

    std::string testFastaFile;
    std::string testDir;

    uint32_t quashLen;
    uint8_t* quash = nullptr; 
};

// Test constructor
TEST_F(ReferenceTest, Constructor) {
    Reference ref(testFastaFile, 4);
    
    EXPECT_EQ(ref.getFastaFileName(), testFastaFile);
    EXPECT_EQ(ref.getBaseGroupLength(), 31);
    EXPECT_EQ(ref.getBaseGroupStep(), 32);
    EXPECT_EQ(ref.getSquash(), nullptr);
    EXPECT_EQ(ref.getSquashLength(), 0);
}

TEST_F(ReferenceTest, ReferenceMakeIndex) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();

    EXPECT_EQ(refe.getFastaFileName(), testFastaFile);
    std::string nifile;
    nifile = refe.getNiFilePath(); 
    EXPECT_EQ(PathUtil::getFileName(nifile), "test_reference.fasta.8d969b67.ni");

    // printBufferHex(refe.getSquash(), refe.getSquashLength());
    EXPECT_EQ(refe.getSquashLength(), quashLen);
    EXPECT_EQ(memcmp(refe.getSquash(), quash, quashLen), 0);
}

