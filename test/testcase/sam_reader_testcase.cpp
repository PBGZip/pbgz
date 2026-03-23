#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>

#define private public
#include "../src/io_wrapper.h"
#include "../src/block_wrapper.h"
#include "../src/io_block.h"
#undef private

namespace SamReaderTestData {
    const std::string smallSamFile = "small_test.sam";
    const std::string largeSamFile = "large_test.sam";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamReaderTest : public ::testing::Test {
public:
    void SetUp() override {
        // 生成小SAM文件
        generateSmallSamFile(SamReaderTestData::smallSamFile);
        // 生成大SAM文件
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
        
        // 写入SAM头部
        file << "@HD    VN:1.6  SO:coordinate\n";
        file << "@SQ	SN:chr1	LN:248956422\n";
        file << "@SQ	SN:chr2	LN:242193529\n";
        file << "@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem hg38.p14.fa ERR016060_1.fastq\n";
        
        // 写入几条SAM记录
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
        
        // 写入SAM头部
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
        
        // 写入大量SAM记录（100条）
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
        
        int32_t result = pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
        
        bool isSam = (pBlock->getBlockType() == SAM);
        
        delete pBlockReader;
        delete pIoReader;
        delete pBlock;
        
        return isSam;
    }
};

TEST_F(SamReaderTest, TestSmallSamFile) {
    // 测试小SAM文件
    bool isSam = isSamFile(SamReaderTestData::smallSamFile);
    EXPECT_TRUE(isSam) << "小SAM文件应该被识别为SAM格式";
}

TEST_F(SamReaderTest, TestLargeSamFile) {
    // 测试大SAM文件
    bool isSam = isSamFile(SamReaderTestData::largeSamFile);
    EXPECT_TRUE(isSam) << "大SAM文件应该被识别为SAM格式";
}

TEST_F(SamReaderTest, TestSamFileReading) {
    // 测试读取SAM文件并验证block类型
    RoughIOBlock* pBlock = new RoughIOBlock(SamReaderTestData::MAX_BLOCK_SIZE);
    
    IOReader* pIoReader = new FileReader(SamReaderTestData::smallSamFile);
    pIoReader->openIO();
    BlockReader* pBlockReader = new BlockReader(pIoReader);
    
    int32_t result = pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
    
    // 验证读取成功
    EXPECT_GT(result, 0) << "读取SAM文件应该成功";
    
    // 验证block类型为SAM
    EXPECT_EQ(pBlock->getBlockType(), SAM) << "读取的block类型应该是SAM";
    
    // 验证block中有数据
    EXPECT_GT(pBlock->getDataLen(), 0) << "block中应该有数据";
    
    delete pBlockReader;
    delete pIoReader;
    delete pBlock;
}

TEST_F(SamReaderTest, TestLargeSamFileReading) {
    // 测试读取大SAM文件并验证block类型
    RoughIOBlock* pBlock = new RoughIOBlock(SamReaderTestData::MAX_BLOCK_SIZE);
    
    IOReader* pIoReader = new FileReader(SamReaderTestData::largeSamFile);
    pIoReader->openIO();
    BlockReader* pBlockReader = new BlockReader(pIoReader);
    
    int32_t result = pBlockReader->readBlock(pBlock, TYPE_UNKNOW);
    
    // 验证读取成功
    EXPECT_GT(result, 0) << "读取大SAM文件应该成功";
    
    // 验证block类型为SAM
    EXPECT_EQ(pBlock->getBlockType(), SAM) << "读取的block类型应该是SAM";
    
    // 验证block中有数据
    EXPECT_GT(pBlock->getDataLen(), 0) << "block中应该有数据";
    
    // 大文件应该有更多数据
    EXPECT_GT(pBlock->getDataLen(), 1000) << "大SAM文件应该包含更多数据";
    
    delete pBlockReader;
    delete pIoReader;
    delete pBlock;
}
