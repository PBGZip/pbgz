#include <cstdint>
#include <string>
#include <fstream>
#include <random>
#include <gtest/gtest.h>
#define  private public
#include "block_wrapper.h"
#undef private

namespace BlockWrapperTestData {
    // Need to construct a 100-line file first
    const std::string samllFastQFile = "small.fastq";
    // Need to construct a 10000-line file first
    const std::string bigFastQFile = "test.fastq";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;

    const char* fisrtLine = "@SRR12922210.1 1/1\n";
};

class BlockWrapperTest : public ::testing::Test {
public:
    void SetUp() override { 
        generateFastqFiles();
    }

    void TearDown() override {
        cleanupFastqFiles();
    }

private:
    void generateFastqFiles() {
        generateFastqFile(BlockWrapperTestData::samllFastQFile, 25);
        generateFastqFile(BlockWrapperTestData::bigFastQFile, 50000);
    }

    void cleanupFastqFiles() {
        std::remove(BlockWrapperTestData::samllFastQFile.c_str());
        std::remove(BlockWrapperTestData::bigFastQFile.c_str());
    }

    void generateFastqFile(const std::string& filename, int numRecords) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            // 如果无法在当前目录打开，尝试在测试目录中创建
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // 如果还是失败，尝试在构建目录中创建
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }
        
        if (!file.is_open()) {
            return; // 无法创建文件
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
};

TEST_F(BlockWrapperTest, TestBlockReaderSmallFile) {
    IOReader* pIoReader = new FileReader(BlockWrapperTestData::samllFastQFile);
    pIoReader->openIO(); 
    BlockReader*  pBlockReader = new BlockReader(pIoReader); 
    RoughIOBlock* pInBlock = new RoughIOBlock(BlockWrapperTestData::MAX_BLOCK_SIZE);
    pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);

    EXPECT_EQ(pInBlock->getBlockType(), FASTQ_GEN2);
    EXPECT_EQ(pInBlock->getNpos().size(), 100);
    EXPECT_EQ(memcmp(pInBlock->getBuffer(), BlockWrapperTestData::fisrtLine, pInBlock->getNpos()[0]),  0);

    delete pIoReader;
    delete pBlockReader;
    delete pInBlock;

}

TEST_F(BlockWrapperTest, TestBlockReaderBigFile) {
    IOReader* pIoReader = new FileReader(BlockWrapperTestData::bigFastQFile);
    pIoReader->openIO(); 
    BlockReader*  pBlockReader = new BlockReader(pIoReader); 
    RoughIOBlock* pInBlock = new RoughIOBlock(BlockWrapperTestData::MAX_BLOCK_SIZE);
    pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);

    EXPECT_EQ(pInBlock->getBlockType(), FASTQ_GEN2);
    EXPECT_EQ(pInBlock->getNpos().size() % 4, 0);
    EXPECT_EQ(pInBlock->getNpos()[0], 18);
    EXPECT_EQ(pInBlock->getNpos()[1], 169);
    EXPECT_EQ(pInBlock->getNpos()[2], 171);
    EXPECT_EQ(pInBlock->getNpos()[3], 322);
    EXPECT_EQ(memcmp(pInBlock->getBuffer(), BlockWrapperTestData::fisrtLine, pInBlock->getNpos()[0]),  0);

    delete pIoReader;
    delete pBlockReader;
    delete pInBlock;
}

TEST_F(BlockWrapperTest, TestBlockReaderTwoBlock) {
    IOReader* pIoReader = new FileReader(BlockWrapperTestData::bigFastQFile);
    pIoReader->openIO(); 
    BlockReader*  pBlockReader = new BlockReader(pIoReader); 
    RoughIOBlock* pInBlock = new RoughIOBlock(BlockWrapperTestData::MAX_BLOCK_SIZE);
    // first block
    pBlockReader->readBlock(pInBlock);
    BlockType fileType = pInBlock->getBlockType();
    EXPECT_EQ(pInBlock->getBlockId(), 0);

    // second block
    pInBlock->reset();
    pBlockReader->readBlock(pInBlock, fileType);
    EXPECT_EQ(pInBlock->getBlockId(), 1);

    EXPECT_EQ(pInBlock->getBlockType(), FASTQ_GEN2);
    EXPECT_EQ(pInBlock->getNpos().size() % 4, 0);

    // printf("%d", pInBlock->getNpos()[0]);
    // char buffer[256];
    // memcpy(buffer, pInBlock->getBuffer(),  32);
    // buffer[32] = 0;
    // printf("%s",buffer);

    // @SRR12922210.25411 25411/1
    EXPECT_EQ(pInBlock->getNpos()[0], 26);
    EXPECT_EQ(pInBlock->getNpos()[1], 177);
    EXPECT_EQ(pInBlock->getNpos()[2], 179);
    EXPECT_EQ(pInBlock->getNpos()[3], 330);

    const char* line = "@SRR12922210.25411 25411/1\n";
    EXPECT_EQ(memcmp(pInBlock->getBuffer(), line, pInBlock->getNpos()[0]),  0);

    delete pIoReader;
    delete pBlockReader;
    delete pInBlock;
}
