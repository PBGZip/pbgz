#include <cstdint>
#include <string>
#include <gtest/gtest.h>
#define  private public
#include "block_wrapper.h"
#undef private

namespace BlockWrapperTestData {
    // 需要先构造一个100行的文件
    const std::string samllFastQFile = "small.fastq";
    // 需要先构造一个10000行的文件
    const std::string bigFastQFile = "test.fastq";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;

    const char* fisrtLine = "@SRR12922210.1 1/1\n";
};

class BlockWrapperTest : public ::testing::Test {
public:
    void SetUp() override { }

    void TearDown() override { }

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



