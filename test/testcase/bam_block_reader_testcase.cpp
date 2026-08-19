/*
 * bam_block_reader_testcase.cpp - BAM 块读取（BAM -> SAM 转换）测试
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

#include <cstring>
#include <fstream>
#include <gtest/gtest.h>
#include <htslib/bgzf.h>

#include "block_wrapper.h"
#include "io_wrapper.h"

namespace {

    const std::string kBamHeader = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n";
    const uint32_t kBlockSize = 8 << 20;

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

    void putI32(std::string& s, int32_t v) { s.append(reinterpret_cast<char*>(&v), 4); }
    void putU16(std::string& s, uint16_t v) { s.append(reinterpret_cast<char*>(&v), 2); }
    void putU8(std::string& s, uint8_t v) { s.push_back(static_cast<char>(v)); }

    int nibble(char c) {
        switch (c) {
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 4;
            case 'T': return 8;
            case 'N': return 15;
            default:  return 0;
        }
    }

    /* 构建一条 BAM 记录：单 CIGAR op（M），aux NM:i:0 / AS:i:70，QUAL 存 phred（ascii-33） */
    std::string buildRecord(const std::string& qname, int32_t flag, int32_t refID, int32_t pos,
                            uint8_t mapq, const std::string& seq, const std::string& qual) {
        std::string rec;
        putI32(rec, 0);                      // block_size 占位
        putI32(rec, refID);
        putI32(rec, pos);                    // 0-based
        putU8(rec, static_cast<uint8_t>(qname.size() + 1));
        putU8(rec, mapq);
        putU16(rec, 0);                      // bin
        putU16(rec, 1);                      // n_cigar_op
        putU16(rec, static_cast<uint16_t>(flag));
        putI32(rec, static_cast<int32_t>(seq.size()));
        putI32(rec, -1);                     // next_refID
        putI32(rec, -1);                     // next_pos
        putI32(rec, 0);                      // tlen
        rec += qname;
        rec.push_back('\0');
        uint32_t cigar = (static_cast<uint32_t>(seq.size()) << 4) | 0;   // L<<4 | M
        rec.append(reinterpret_cast<char*>(&cigar), 4);
        for (size_t i = 0; i < seq.size(); i += 2) {
            const uint8_t hi = static_cast<uint8_t>(nibble(seq[i]) & 0xF);
            const uint8_t lo = (i + 1 < seq.size()) ? static_cast<uint8_t>(nibble(seq[i + 1]) & 0xF) : 0;
            rec.push_back(static_cast<char>((hi << 4) | lo));
        }
        for (char c : qual) {
            rec.push_back(static_cast<char>(c - 33));
        }
        rec += "NM";
        rec.push_back('i');
        putI32(rec, 0);
        rec += "AS";
        rec.push_back('i');
        putI32(rec, 70);
        int32_t bs = static_cast<int32_t>(rec.size() - 4);
        memcpy(&rec[0], &bs, 4);
        return rec;
    }

    /* 构建原始 BAM 字节流（未压缩） */
    std::string buildBam(int numReads, const std::string& headerText) {
        std::string bam = "BAM\x1";
        putI32(bam, static_cast<int32_t>(headerText.size()));
        bam += headerText;
        putI32(bam, 1);                      // n_ref
        putI32(bam, 5);                      // "chr1\0"
        bam += "chr1"; bam.push_back('\0');
        putI32(bam, 100000);
        for (int i = 0; i < numReads; ++i) {
            bam += buildRecord("r" + std::to_string(i), 0, 0, i * 2, 60, kSeq76, kQual76);
        }
        return bam;
    }

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

    std::string blockText(RoughIOBlock& block) {
        return std::string(reinterpret_cast<char*>(block.getBuffer()),
                           static_cast<size_t>(block.getDataLen()));
    }

}  // namespace

/* BAM 头部独立成块，数据区解压成 SAM 行组成 SAM 块 */
TEST(BamBlockReaderTest, HeaderIndependentBlockAndSamConversion) {
    const std::string bam = buildBam(5, kBamHeader);
    writeFile("ut_bam.bam", bam);

    IOReader* io = new FileReader("ut_bam.bam");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BamBlockReader*>(reader), nullptr);
    EXPECT_EQ(dynamic_cast<BamGzBlockReader*>(reader), nullptr);

    RoughIOBlock block(kBlockSize);

    // block 0：SAM 头部块（只含 @ 行，不携带数据）
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);
    EXPECT_FALSE(reader->blockHasData(&block));
    EXPECT_EQ(blockText(block), kBamHeader);
    EXPECT_EQ(block.getNpos().size(), 2);   // @HD / @SQ 两行

    // block 1：数据块（5 条 read 转成 SAM 行）
    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);
    EXPECT_TRUE(reader->blockHasData(&block));
    EXPECT_EQ(block.getNpos().size(), 5);

    // 校验第一条转换出的 SAM 行
    const std::string line = std::string(
        reinterpret_cast<char*>(block.getBuffer()), block.getNpos()[0]);
    const std::string expected =
        "r0\t0\tchr1\t1\t60\t76M\t*\t0\t0\t" + kSeq76 + "\t" + kQual76 + "\tNM:i:0\tAS:i:70";
    EXPECT_EQ(line, expected);

    // EOF
    block.reset();
    EXPECT_EQ(reader->readBlock(&block, TYPE_UNKNOW), 0);

    delete reader;
    delete io;
    std::remove("ut_bam.bam");
}

/* BAM 头文本缺失 @SQ 时，由参考序列表生成 @SQ 行 */
TEST(BamBlockReaderTest, GenerateSqFromRefs) {
    const std::string headerOnly = "@HD\tVN:1.6\n";
    const std::string bam = buildBam(2, headerOnly);
    writeFile("ut_bam.bam", bam);

    IOReader* io = new FileReader("ut_bam.bam");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);

    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(blockText(block), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n");

    delete reader;
    delete io;
    std::remove("ut_bam.bam");
}

/* 数据区按 readsPerBlock 条 read 分块 */
TEST(BamBlockReaderTest, ReadCountSplit) {
    const std::string bam = buildBam(7, kBamHeader);
    writeFile("ut_bam.bam", bam);

    IOReader* io = new FileReader("ut_bam.bam");
    ASSERT_EQ(io->openIO(), 0);
    /* 堆上构造：BlockReader 带 256MB cache 成员，栈上会溢出 */
    BamBlockReader* reader = new BamBlockReader(io, nullptr, 0, 3);   // readsPerBlock=3, splitHeader=true

    RoughIOBlock block(kBlockSize);
    // block 0：头部块
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getNpos().size(), 2);
    EXPECT_FALSE(reader->blockHasData(&block));

    // 3 / 3 / 1 条 read 三个数据块
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
    std::remove("ut_bam.bam");
}

/* splitHeader=false：头部与首个数据块合并，排序等下游可自含 @SQ */
TEST(BamBlockReaderTest, HeaderMergedMode) {
    const std::string bam = buildBam(5, kBamHeader);
    writeFile("ut_bam.bam", bam);

    IOReader* io = new FileReader("ut_bam.bam");
    ASSERT_EQ(io->openIO(), 0);
    /* 堆上构造：BlockReader 带 256MB cache 成员，栈上会溢出 */
    BamBlockReader* reader = new BamBlockReader(io, nullptr, 0, 3, false);   // splitHeader=false

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
    std::remove("ut_bam.bam");
}

/* BGZF 压缩的 BAM（.bam.gz）：BamGzBlockReader 先 inflate 再转 SAM */
TEST(BamBlockReaderTest, GzBamInflate) {
    const std::string bam = buildBam(4, kBamHeader);
    writeBgzfFile("ut_bam.bam.gz", bam);

    IOReader* io = new FileReader("ut_bam.bam.gz");
    ASSERT_EQ(io->openIO(), 0);
    BlockReader* reader = BlockFactory::createBlockReader(io, 5);
    ASSERT_NE(reader, nullptr);
    EXPECT_NE(dynamic_cast<BamGzBlockReader*>(reader), nullptr);

    RoughIOBlock block(kBlockSize);
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);
    EXPECT_EQ(block.getNpos().size(), 2);     // 头部块

    block.reset();
    ASSERT_GT(reader->readBlock(&block, TYPE_UNKNOW), 0);
    EXPECT_EQ(block.getBlockType(), BAM);
    EXPECT_EQ(block.getNpos().size(), 4);     // 数据块

    const std::string line = std::string(
        reinterpret_cast<char*>(block.getBuffer()), block.getNpos()[0]);
    EXPECT_EQ(line,
              "r0\t0\tchr1\t1\t60\t76M\t*\t0\t0\t" + kSeq76 + "\t" + kQual76 + "\tNM:i:0\tAS:i:70");

    delete reader;
    delete io;
    std::remove("ut_bam.bam.gz");
}
