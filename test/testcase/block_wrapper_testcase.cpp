/*
 * block_wrapper_testcase.cpp - Test cases for block wrapper
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

#include <cstdint>
#include <cstring>
#include <string>
#include <fstream>
#include <random>
#include <gtest/gtest.h>
#include <htslib/bgzf.h>

#define  private public
#define  protected public
#include "block_wrapper.h"
#include "io_wrapper.h"
#undef protected
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
    BlockReader*  pBlockReader = BlockFactory::createBlockReader(pIoReader, 3); 
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
    BlockReader*  pBlockReader = BlockFactory::createBlockReader(pIoReader, 3); 
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
    BlockReader*  pBlockReader = BlockFactory::createBlockReader(pIoReader, 3); 
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

namespace {

    const uint32_t kBlockSize = 8 << 20;

    void writeFile(const std::string& path, const std::string& data) {
        std::ofstream f(path, std::ios::binary);
        f.write(data.data(), static_cast<std::streamsize>(data.size()));
        f.close();
    }

    void writeBgzfFile(const std::string& path, const std::string& data) {
        BGZF* f = bgzf_open(path.c_str(), "w");
        ASSERT_NE(f, nullptr);
        ASSERT_GT(bgzf_write(f, data.data(), data.size()), 0);
        bgzf_close(f);
    }

    /* 普通 gzip（非 BGZF）：zlib gzopen 写出的 .gz 没有 FEXTRA 标志 */
    void writePlainGzipFile(const std::string& path, const std::string& data) {
        gzFile f = gzopen(path.c_str(), "wb");
        ASSERT_NE(f, nullptr);
        ASSERT_GT(gzwrite(f, data.data(), (unsigned)data.size()), 0);
        gzclose(f);
    }

    const std::string kSeq76 = [] {
        std::string s;
        for (int i = 0; i < 19; ++i) s += "ACGT";
        return s;
    }();
    const std::string kQual76 = [] {
        std::string s;
        for (int i = 0; i < 76; ++i) s += 'I';
        return s;
    }();

    std::string makeFastq(int numRecords) {
        std::string s;
        for (int i = 0; i < numRecords; ++i) {
            s += "@read_" + std::to_string(i) + "\n";
            s += kSeq76 + "\n";
            s += "+\n";
            s += kQual76 + "\n";
        }
        return s;
    }

    std::string makeSam(int numReads) {
        std::string s = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n";
        for (int i = 0; i < numReads; ++i) {
            s += "read" + std::to_string(i) + "\t0\tchr1\t" + std::to_string(i * 10 + 1) +
                 "\t60\t76M\t*\t0\t0\t" + kSeq76 + "\t" + kQual76 + "\n";
        }
        return s;
    }

    /* 最小 BAM（头部 + 2 条记录） */
    std::string makeBam() {
        std::string bam = "BAM\x1";
        const std::string header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n";
        int32_t lText = static_cast<int32_t>(header.size());
        bam.append(reinterpret_cast<char*>(&lText), 4);
        bam += header;
        int32_t nRef = 1;
        bam.append(reinterpret_cast<char*>(&nRef), 4);
        int32_t lName = 5;
        bam.append(reinterpret_cast<char*>(&lName), 4);
        bam += "chr1"; bam.push_back('\0');
        int32_t refLen = 100000;
        bam.append(reinterpret_cast<char*>(&refLen), 4);
        for (int i = 0; i < 2; ++i) {
            std::string rec;
            rec.append(4, '\0');                       // block_size 占位
            int32_t refID = 0;
            int32_t pos = i * 2;
            rec.append(reinterpret_cast<char*>(&refID), 4);
            rec.append(reinterpret_cast<char*>(&pos), 4);
            uint8_t lrn = 3;
            uint8_t mapq = 60;
            rec.push_back(static_cast<char>(lrn));
            rec.push_back(static_cast<char>(mapq));
            uint16_t bin = 0, ncig = 1, flag = 0;
            rec.append(reinterpret_cast<char*>(&bin), 2);
            rec.append(reinterpret_cast<char*>(&ncig), 2);
            rec.append(reinterpret_cast<char*>(&flag), 2);
            int32_t lseq = 76, nid = -1, npos = -1, tlen = 0;
            rec.append(reinterpret_cast<char*>(&lseq), 4);
            rec.append(reinterpret_cast<char*>(&nid), 4);
            rec.append(reinterpret_cast<char*>(&npos), 4);
            rec.append(reinterpret_cast<char*>(&tlen), 4);
            rec += "r" + std::to_string(i);
            rec.push_back('\0');
            uint32_t cigar = (76u << 4);
            rec.append(reinterpret_cast<char*>(&cigar), 4);
            for (int j = 0; j < 38; ++j) {
                rec.push_back(static_cast<char>(0x12 + (j % 2) * 0x36));   // ACGT/ACGT
            }
            for (int j = 0; j < 76; ++j) rec.push_back(static_cast<char>('I' - 33));
            rec += "NM";
            rec.push_back('i');
            int32_t nm = 0;
            rec.append(reinterpret_cast<char*>(&nm), 4);
            int32_t bs = static_cast<int32_t>(rec.size() - 4);
            memcpy(&rec[0], &bs, 4);
            bam += rec;
        }
        return bam;
    }

    std::string blockText(RoughIOBlock& block) {
        return std::string(reinterpret_cast<char*>(block.getBuffer()),
                           static_cast<size_t>(block.getDataLen()));
    }

    /* ---- BamWriter 测试辅助 ---- */

    const std::string kSamHeader =
        "@HD\tVN:1.6\n"
        "@SQ\tSN:chr1\tLN:100000\n"
        "@SQ\tSN:chr2\tLN:200000\n";

    std::string makeSamRead(int i, int pos) {
        return "r" + std::to_string(i) + "\t0\tchr1\t" + std::to_string(pos) +
               "\t60\t76M\t=\t" + std::to_string(pos + 10) + "\t0\t" + kSeq76 + "\t" +
               kQual76 + "\tNM:i:0\tAS:i:" + std::to_string(60 + i) + "\n";
    }

    /* 把 SAM 文本填入一个类型为 SAM 的块 */
    void fillSamBlock(RoughIOBlock* block, const std::string& text) {
        block->reset();
        ASSERT_GE(block->getBufferSize(), text.size());
        memcpy(block->getBuffer(), text.data(), text.size());
        block->setDataLen((int64_t)text.size());
        block->setBlockType(SAM);
        block->getNpos().clear();
        for (size_t i = 0; i < text.size(); ++i) {
            if (text[i] == '\n') {
                block->getNpos().push_back(i);
            }
        }
    }

    std::string readFile(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        return s;
    }

}  // namespace

/* ---- BlockFactory 格式探测 ---- */

TEST(BlockFactoryFormatTest, DetectFastq) {
    writeFile("ut_f.fq", makeFastq(4));
    IOReader* io = new FileReader("ut_f.fq");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<FastqBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), FASTQ_GEN2);
    delete reader;
    delete io;
    std::remove("ut_f.fq");
}

TEST(BlockFactoryFormatTest, DetectGzFastq) {
    writeBgzfFile("ut_f.fq.gz", makeFastq(4));
    IOReader* io = new GzFileReader("ut_f.fq.gz", 1);
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<FastqBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    /* GZ 输入由 io 层透明解压，块类型仍按非 GZ 的 FASTQ 设置 */
    EXPECT_EQ(block.getBlockType(), FASTQ_GEN2);
    delete reader;
    delete io;
    std::remove("ut_f.fq.gz");
}

TEST(BlockFactoryFormatTest, DetectSam) {
    writeFile("ut_s.sam", makeSam(3));
    IOReader* io = new FileReader("ut_s.sam");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<SamBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), SAM);
    EXPECT_FALSE(reader->blockHasData(&block));   // 首块是纯头部块
    delete reader;
    delete io;
    std::remove("ut_s.sam");
}

TEST(BlockFactoryFormatTest, DetectGzSam) {
    writeBgzfFile("ut_s.sam.gz", makeSam(3));
    IOReader* io = new GzFileReader("ut_s.sam.gz", 1);
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<SamBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    /* GZ 输入由 io 层透明解压，块类型仍按非 GZ 的 SAM 设置 */
    EXPECT_EQ(block.getBlockType(), SAM);
    delete reader;
    delete io;
    std::remove("ut_s.sam.gz");
}

TEST(BlockFactoryFormatTest, DetectBam) {
    writeFile("ut_b.bam", makeBam());
    IOReader* io = new FileReader("ut_b.bam");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BamBlockReader*>(reader), nullptr);
    EXPECT_EQ(dynamic_cast<BamGzBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);         // BAM 转成 SAM 块，但类型标记为 BAM
    EXPECT_FALSE(reader->blockHasData(&block));   // 首块为头部块
    delete reader;
    delete io;
    std::remove("ut_b.bam");
}

TEST(BlockFactoryFormatTest, DetectGzBam) {
    writeBgzfFile("ut_b.bam.gz", makeBam());
    IOReader* io = new FileReader("ut_b.bam.gz");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BamGzBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);
    delete reader;
    delete io;
    std::remove("ut_b.bam.gz");
}

TEST(BlockFactoryFormatTest, DetectBinary) {
    std::string bin;
    for (int i = 0; i < 1024; ++i) bin += static_cast<char>((i * 37) & 0xFF);
    writeFile("ut_bin.dat", bin);
    IOReader* io = new FileReader("ut_bin.dat");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BinaryBlockReader*>(reader), nullptr);
    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BINARY);
    delete reader;
    delete io;
    std::remove("ut_bin.dat");
}

/* ---- 文件级 BAM 判定（引擎据此决定 ioReader 是否走透明 gz 解压） ---- */

TEST(BlockFactoryFormatTest, IsBamFile) {
    writeFile("ut_i.bam", makeBam());                        // 原始 BAM（非 gz）
    writeBgzfFile("ut_i.bam.bgz", makeBam());                // BGZF 压缩的 BAM
    writeBgzfFile("ut_i.sam.gz", makeSam(3));                // gz 压缩的 SAM
    writeFile("ut_i.sam", makeSam(3));                       // 原始 SAM

    /* BGZF 的 BAM 逐字节是 gzip 流，必须能识别为 BAM，避免走透明 gz 解压丢文件总长 */
    EXPECT_TRUE(BlockUtil::isBamFile("ut_i.bam"));
    EXPECT_TRUE(BlockUtil::isBamFile("ut_i.bam.bgz"));
    EXPECT_FALSE(BlockUtil::isBamFile("ut_i.sam"));
    EXPECT_FALSE(BlockUtil::isBamFile("ut_i.sam.gz"));

    std::remove("ut_i.bam");
    std::remove("ut_i.bam.bgz");
    std::remove("ut_i.sam.gz");
    std::remove("ut_i.sam");
}

/* ---- 按内容探测输入格式（SAM/BAM 决定参考是否只需 squash） ---- */

TEST(BlockFactoryFormatTest, DetectInputFileType) {
    writeFile("ut_t.sam", makeSam(3));
    writeBgzfFile("ut_t.sam.gz", makeSam(3));
    writePlainGzipFile("ut_t.sam.gz2", makeSam(3));
    writeFile("ut_t.fq", makeFastq(4));
    writeBgzfFile("ut_t.fq.gz", makeFastq(4));
    writePlainGzipFile("ut_t.fq.gz2", makeFastq(4));
    writeFile("ut_t.bam", makeBam());
    writeBgzfFile("ut_t.bam.bgz", makeBam());
    std::string bin;
    for (int i = 0; i < 512; ++i) bin += static_cast<char>((i * 13) & 0xFF);
    writeFile("ut_t.bin", bin);

    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.sam"), SAM);
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.sam.gz"), SAM);     // gz 自动解压探测
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.sam.gz2"), SAM);    // 普通 gzip（无 FEXTRA）
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.fq"), FASTQ_GEN2);
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.fq.gz"), FASTQ_GEN2);
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.fq.gz2"), FASTQ_GEN2);
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.bam"), BAM);
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.bam.bgz"), BAM);     // BGZF 解压后是 BAM
    EXPECT_EQ(BlockUtil::detectInputFileType("ut_t.bin"), BINARY);
    EXPECT_EQ(BlockUtil::detectInputFileType("no_such_file.dat"), TYPE_UNKNOW);

    std::remove("ut_t.sam");
    std::remove("ut_t.sam.gz");
    std::remove("ut_t.sam.gz2");
    std::remove("ut_t.fq");
    std::remove("ut_t.fq.gz");
    std::remove("ut_t.fq.gz2");
    std::remove("ut_t.bam");
    std::remove("ut_t.bam.bgz");
    std::remove("ut_t.bin");
}

/* ---- SAM 块读取：头部独立成块 + 按 read 数分块 ---- */

TEST(SamBlockReaderTest, HeaderIndependentBlockAndReadCountSplit) {
    writeFile("ut_sam.sam", makeSam(7));
    IOReader* io = new FileReader("ut_sam.sam");
    ASSERT_EQ(io->openIO(), 0);
    /* 堆上构造：BlockReader 带 256MB cache 成员，栈上会溢出 */
    SamBlockReader* reader = new SamBlockReader(io, nullptr, 0, 3);   // readsPerBlock=3, splitHeader=true

    RoughIOBlock block(kBlockSize);
    // block 0：纯头部块
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), SAM);
    EXPECT_FALSE(reader->blockHasData(&block));
    EXPECT_EQ(block.getNpos().size(), 2);
    EXPECT_EQ(blockText(block), "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n");

    // 3 / 3 / 1 条 read 数据块
    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 3);
    EXPECT_TRUE(reader->blockHasData(&block));

    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 3);

    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 1);

    block.reset();
    EXPECT_EQ(reader->readBlock(&block, TYPE_UNKNOW), 0);   // EOF

    delete reader;
    delete io;
    std::remove("ut_sam.sam");
}

TEST(SamBlockReaderTest, HeaderMergedMode) {
    writeFile("ut_sam.sam", makeSam(5));
    IOReader* io = new FileReader("ut_sam.sam");
    ASSERT_EQ(io->openIO(), 0);
    /* 堆上构造：BlockReader 带 256MB cache 成员，栈上会溢出 */
    SamBlockReader* reader = new SamBlockReader(io, nullptr, 0, 3, false);   // splitHeader=false

    RoughIOBlock block(kBlockSize);
    // block 0：头部（2 行）+ 3 条 read
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 5);
    EXPECT_TRUE(reader->blockHasData(&block));

    // block 1：剩余 2 条 read
    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 2);

    delete reader;
    delete io;
    std::remove("ut_sam.sam");
}

/* 压缩级别决定 readsPerBlock：1-5->10000，6-7->25000，8-9->100000 */
TEST(SamBlockReaderTest, ReadsPerBlockByCompressLevel) {
    writeFile("ut_sam.sam", makeSam(4));
    IOReader* io = new FileReader("ut_sam.sam");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 1);
    ASSERT_NE(reader, nullptr);
    SamBlockReader* samReader = dynamic_cast<SamBlockReader*>(reader);
    ASSERT_NE(samReader, nullptr);
    EXPECT_EQ(samReader->readsPerBlock, 10000u);
    delete reader;
    delete io;
    std::remove("ut_sam.sam");
}

/* ---- BamWriter（SAM -> BAM 转换） ---- */

/* SAM 头 + 数据块写入 BamWriter 后，产出可被 BamBlockReader 读回的原 SAM 内容 */
TEST(BamWriterTest, SamToBamRoundtrip) {
    const std::string data = kSamHeader + makeSamRead(0, 1) + makeSamRead(1, 21) + makeSamRead(2, 41);

    IOWriter* io = new FileWriter("ut_bw.bam");
    ASSERT_EQ(io->openIO(), 0);
    BamWriter* writer = new BamWriter(io);

    RoughIOBlock block(kBlockSize);
    fillSamBlock(&block, kSamHeader);
    ASSERT_EQ(writer->writeBlock(&block), 0);

    block.reset();
    fillSamBlock(&block, makeSamRead(0, 1) + makeSamRead(1, 21) + makeSamRead(2, 41));
    ASSERT_EQ(writer->writeBlock(&block), 0);

    ASSERT_EQ(writer->finish(), 0);
    delete writer;
    io->closeIO();
    delete io;

    /* 用 BAM 读取侧（BGZF 分支）读回 */
    IOReader* rd = new FileReader("ut_bw.bam");
    ASSERT_EQ(rd->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(rd, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BamGzBlockReader*>(reader), nullptr);

    RoughIOBlock out(kBlockSize);
    ASSERT_GT(reader->readBlock(&out, TYPE_UNKNOW), 0);
    EXPECT_EQ(out.getBlockType(), BAM);
    EXPECT_EQ(blockText(out), kSamHeader);   // 头块

    std::string lines;
    while (true) {
        out.reset();
        int64_t ret = reader->readBlock(&out, TYPE_UNKNOW);
        if (ret <= 0) {
            break;
        }
        lines += blockText(out);
    }
    EXPECT_EQ(lines, makeSamRead(0, 1) + makeSamRead(1, 21) + makeSamRead(2, 41));

    delete reader;
    delete rd;
    std::remove("ut_bw.bam");
}

/* 头部未独立成块（头部行 + 数据行同块）：BamWriter 先写头再转数据行 */
TEST(BamWriterTest, MergedHeaderAndDataBlock) {
    const std::string text = kSamHeader + makeSamRead(0, 1) + makeSamRead(1, 21);

    IOWriter* io = new FileWriter("ut_bw.bam");
    ASSERT_EQ(io->openIO(), 0);
    BamWriter* writer = new BamWriter(io);

    RoughIOBlock block(kBlockSize);
    fillSamBlock(&block, text);
    ASSERT_EQ(writer->writeBlock(&block), 0);
    ASSERT_EQ(writer->finish(), 0);
    delete writer;
    io->closeIO();
    delete io;

    IOReader* rd = new FileReader("ut_bw.bam");
    ASSERT_EQ(rd->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(rd, 5);
    ASSERT_NE(reader, nullptr);

    RoughIOBlock out(kBlockSize);
    ASSERT_GT(reader->readBlock(&out, TYPE_UNKNOW), 0);
    EXPECT_EQ(blockText(out), kSamHeader);
    out.reset();
    ASSERT_GT(reader->readBlock(&out, TYPE_UNKNOW), 0);
    EXPECT_EQ(blockText(out), makeSamRead(0, 1) + makeSamRead(1, 21));

    delete reader;
    delete rd;
    std::remove("ut_bw.bam");
}

/* 非 SAM 块透传：-b 误用于 FASTQ/二进制时原样写出，不丢数据 */
TEST(BamWriterTest, NonSamPassThrough) {
    const std::string bin(100, '\xAA');

    IOWriter* io = new FileWriter("ut_bw.bin");
    ASSERT_EQ(io->openIO(), 0);
    BamWriter* writer = new BamWriter(io);

    RoughIOBlock block(kBlockSize);
    block.reset();
    memcpy(block.getBuffer(), bin.data(), bin.size());
    block.setDataLen((int64_t)bin.size());
    block.setBlockType(BINARY);
    ASSERT_EQ(writer->writeBlock(&block), 0);
    ASSERT_EQ(writer->finish(), 0);
    delete writer;
    io->closeIO();
    delete io;

    EXPECT_EQ(readFile("ut_bw.bin"), bin);
    std::remove("ut_bw.bin");
}

/* BGZF 输出是标准 BAM 容器：首字节为 gzip 魔数且带 BC/BS 扩展头 */
TEST(BamWriterTest, BgzfContainer) {
    IOWriter* io = new FileWriter("ut_bw.bam");
    ASSERT_EQ(io->openIO(), 0);
    BamWriter* writer = new BamWriter(io);

    RoughIOBlock block(kBlockSize);
    fillSamBlock(&block, kSamHeader);
    ASSERT_EQ(writer->writeBlock(&block), 0);
    ASSERT_EQ(writer->finish(), 0);
    delete writer;
    io->closeIO();
    delete io;

    const std::string out = readFile("ut_bw.bam");
    ASSERT_GE(out.size(), 28u);
    EXPECT_EQ((uint8_t)out[0], 0x1f);
    EXPECT_EQ((uint8_t)out[1], 0x8b);
    EXPECT_EQ((uint8_t)out[2], 0x08);
    EXPECT_EQ((uint8_t)out[3], 0x04);
    EXPECT_EQ((uint8_t)out[10], 0x06);   // XLEN = 6
    EXPECT_EQ((uint8_t)out[12], 0x42);   // "BC"
    EXPECT_EQ((uint8_t)out[13], 0x43);
    std::remove("ut_bw.bam");
}
