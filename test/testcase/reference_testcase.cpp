#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <memory>
#include <chrono>
#include <thread>
#include <vector>
#include "reference.h"
#include "pbgz_errno.h"
#include "utils/path_util.h"

class ReferenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建测试用的临时FASTA文件
        testFastaFile = "test_reference.fasta";
        std::ofstream fastaFile(testFastaFile);
        fastaFile << ">test_sequence\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
        fastaFile.close();
        
        // 创建测试目录
        testDir = "test_reference_dir";
        std::filesystem::create_directory(testDir);
    }

    void TearDown() override {
        // 清理测试文件
        std::filesystem::remove(testFastaFile);
        std::filesystem::remove_all(testDir);
    }

    std::string testFastaFile;
    std::string testDir;
};

// 测试构造函数
TEST_F(ReferenceTest, Constructor) {
    Reference ref(testFastaFile, 4);
    
    EXPECT_EQ(ref.getFastaFileName(), testFastaFile);
    EXPECT_EQ(ref.getBaseGroupLength(), 31);
    EXPECT_EQ(ref.getBaseGroupStep(), 32);
    EXPECT_EQ(ref.getSquash(), nullptr);
    EXPECT_EQ(ref.getSquashLength(), 0);
}
