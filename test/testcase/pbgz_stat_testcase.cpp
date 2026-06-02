/*
 * pbgz_stat_testcase.cpp - Test cases for pbgz_stat
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
#include <iostream>
#include <sstream>
#include <vector>
#include "../src/pbgz_stat.h"

class PbgzStatTest : public ::testing::Test {
public:
    void SetUp() override {
        // 测试前的准备工作
    }
    
    void TearDown() override {
        // 测试后的清理工作
    }
};

// 测试 FastqStat 初始化
TEST_F(PbgzStatTest, FastqStatInitialization) {
    FastqStat fastqStat;
    int result = fastqStat.init();
    EXPECT_EQ(result, 0);
}

// 测试 SamStat 初始化
TEST_F(PbgzStatTest, SamStatInitialization) {
    SamStat samStat;
    int result = samStat.init();
    EXPECT_EQ(result, 0);
}

// 测试 FastqStat 设置指标值
TEST_F(PbgzStatTest, FastqStatSetValue) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 100ULL);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 50);
    value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 50ULL);
}

// 测试累加指标值
TEST_F(PbgzStatTest, FastqStatAccumulateValue) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 10);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 20);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 15);
    
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 45ULL);
}

// 测试压缩率计算
TEST_F(PbgzStatTest, CompressionRatioCalculation) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    // 压缩率 = (200 / 1000) * 10000 = 2000 (即 20%)
    EXPECT_EQ(ratio, 2000ULL);
}

// 测试原始大小为0时的压缩率计算
TEST_F(PbgzStatTest, CompressionRatioZeroOriginal) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 0);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO);
    // 原始大小为0时，压缩率应为0
    EXPECT_EQ(ratio, 0ULL);
}

// 测试 SamStat 设置指标值
TEST_F(PbgzStatTest, SamStatSetValue) {
    SamStat samStat;
    samStat.init();
    
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 200);
    uint64_t value = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 200ULL);
    
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSED_SIZE, 80);
    value = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 80ULL);
}

// 测试 SamStat 累加指标值
TEST_F(PbgzStatTest, SamStatAccumulateValue) {
    SamStat samStat;
    samStat.init();
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE, 5);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE, 10);
    
    uint64_t value = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 15ULL);
}

// 测试获取不存在的指标
TEST_F(PbgzStatTest, GetNonExistentMetric) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 获取一个未设置的指标
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 0ULL);
}

// 测试全局指标设置
TEST_F(PbgzStatTest, GlobalMetrics) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 注意：先调用 calculateAllMetrics 来确保压缩率计算
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    
    // 现在测试压缩率是否正确计算
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 2000ULL);  // (200/1000)*10000 = 2000
}

// 测试字段名称存储
TEST_F(PbgzStatTest, FieldNamesStorage) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 由于printStats()现在使用fprintf输出到stderr，我们直接调用即可验证功能
    // 这个测试主要验证printStats()不会崩溃，并且能正常输出
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    // printStats()会计算并输出统计数据
    fastqStat.printStats();
    
    // 验证统计数据计算正确
    uint64_t originalSize = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t compressedSize = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(originalSize, 100ULL);
    EXPECT_EQ(compressedSize, 20ULL);
}

// 测试命名空间访问
TEST_F(PbgzStatTest, NamespaceAccess) {
    // 测试 StatUnitId 命名空间
    EXPECT_EQ(StatObjectId::FASTQ_ID, 1);
    EXPECT_EQ(StatObjectId::FASTQ_BASE, 2);
    EXPECT_EQ(StatObjectId::FASTQ_COMMENT, 3);
    EXPECT_EQ(StatObjectId::FASTQ_QUALITY, 4);
    
    EXPECT_EQ(StatObjectId::SAM_QNAME, 5);
    EXPECT_EQ(StatObjectId::SAM_FLAG, 6);
    EXPECT_EQ(StatObjectId::SAM_RNAME, 7);
    EXPECT_EQ(StatObjectId::SAM_POS, 8);
    EXPECT_EQ(StatObjectId::SAM_SEQ, 14);
    EXPECT_EQ(StatObjectId::SAM_QUAL, 15);
    
    // 测试 StatMetricIds 命名空间
    EXPECT_EQ(StatMetricIds::ORIGINAL_SIZE, 1);
    EXPECT_EQ(StatMetricIds::COMPRESSED_SIZE, 2);
    EXPECT_EQ(StatMetricIds::COMPRESSION_RATIO, 3);
    
    // 测试 StatUnitNames 字符串
    EXPECT_EQ(StatObjectNames::Fastq::ID, "ID");
    EXPECT_EQ(StatObjectNames::Fastq::BASE, "Base");
    EXPECT_EQ(StatObjectNames::Fastq::COMMENT, "Comment");
    EXPECT_EQ(StatObjectNames::Fastq::QUALITY, "Quality");
    
    EXPECT_EQ(StatObjectNames::Sam::QNAME, "QNAME");
    EXPECT_EQ(StatObjectNames::Sam::FLAG, "FLAG");
    EXPECT_EQ(StatObjectNames::Sam::SEQ, "SEQ");
    EXPECT_EQ(StatObjectNames::Sam::QUAL, "QUAL");
}

// 测试多个字段的统计
TEST_F(PbgzStatTest, MultipleFields) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 为多个字段设置数据
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 500);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE, 10);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE, 40);
    
    // 计算压缩率
    
    // 验证每个字段的压缩率都正确计算
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
}

// 测试不同压缩率的计算
TEST_F(PbgzStatTest, DifferentCompressionRatios) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 设置不同的压缩率
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);  // 20% -> 2000
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 2000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 400); // 20% -> 2000
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 500);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE, 150); // 30% -> 3000
    
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSION_RATIO), 3000ULL);
}

// 测试重新计算压缩率
TEST_F(PbgzStatTest, RecalculateCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 初始设置
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    
    // 更新压缩后大小，重新计算
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 100);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 3000ULL);  // (300/1000)*10000
}

// 测试多种指标类型的累加
TEST_F(PbgzStatTest, AccumulationAndSetting) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 先 set，再add
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 50);
    
    // 原始大小类型为ACCUMULATED，所以 add 会累加
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 150ULL);  // 100 + 50 = 150
    
    // COMPRESSED_SIZE 类型是 ACCUMULATED，add 会有作用
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 30);
    
    value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 180ULL);
}

// 测试计算指标的更新顺序
TEST_F(PbgzStatTest, CalculatedMetricUpdateOrder) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 只设置原始大小
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    
    // 压缩率因为原大非0但压缩为0，应该为0
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 0ULL);
    
    // 设置压缩后大小后重新计算
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    // 现在应该有正确的压缩率
    ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 2000ULL);
}

// 测试独立的统计对象
TEST_F(PbgzStatTest, IndependentStatObjects) {
    FastqStat fastqStat1;
    FastqStat fastqStat2;
    
    fastqStat1.init();
    fastqStat2.init();
    
    // 不同的对象，数据互不影响
    fastqStat1.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat1.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat2.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat2.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 40);
    
    
    // 验证两个对象的数据独立
    EXPECT_EQ(fastqStat1.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat2.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 200ULL);
}

// 测试重复初始化
TEST_F(PbgzStatTest, ReInitialization) {
    FastqStat fastqStat;
    
    fastqStat.init();
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    
    // 重新初始化应该清除之前的设置
    fastqStat.init();
    
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 0ULL);
}

// 测试所有字段统计
TEST_F(PbgzStatTest, AllFastqFields) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 测试所有4个FASTQ字段都正确注册
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 100);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE), 100ULL);
}

// 测试所有SAM字段统计
TEST_F(PbgzStatTest, AllSamFields) {
    SamStat samStat;
    samStat.init();
    
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatMetricIds::ORIGINAL_SIZE, 100);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::ORIGINAL_SIZE), 100ULL);
}

// 测试FastqStat多次累加数据
TEST_F(PbgzStatTest, FastqStatMultipleAccumulations) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 第一次累加ID字段
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 50);
    
    // 第二次累加ID字段
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 150);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 30);
    
    // 验证累加结果：100+150=250, 50+30=80
    uint64_t totalOriginal = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 250ULL);
    EXPECT_EQ(totalCompressed, 80ULL);
}

// 测试SamStat多次累加数据
TEST_F(PbgzStatTest, SamStatMultipleAccumulations) {
    SamStat samStat;
    samStat.init();
    
    // 第一次累加QNAME字段
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 25);
    
    // 第二次累加QNAME字段
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 200);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 45);
    
    // 第三次累加QNAME字段
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 150);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 60);
    
    // 验证累加结果：100+200+150=450, 25+45+60=130
    uint64_t totalOriginal = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 450ULL);
    EXPECT_EQ(totalCompressed, 130ULL);
}

// 测试累加和设置混合使用
TEST_F(PbgzStatTest, MixedSetAndAccumulate) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 先设置初始值
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 80);
    
    // 然后累加
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    // 最后再设置，但通常累加用于动态累加，set用于初始化
    // 这里我们验证结果是set + add = 150, 100
    uint64_t totalOriginal = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 150ULL);
    EXPECT_EQ(totalCompressed, 100ULL);
}

// 测试空字段的压缩率
TEST_F(PbgzStatTest, EmptyFieldCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 0);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 0);
    
    
    // 当原始大小为0时，压缩率应该为0
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 0ULL);
}

// 测试所有FASTQ字段的累加
TEST_F(PbgzStatTest, FastqStatAllFieldsAccumulate) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 对所有4个字段进行多次累加
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 500);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 10);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE, 5);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 300);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE, 150);
    
    // 第二次累加ID字段
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 40);
    
    
    // 验证累加结果
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 300ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE), 60ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE), 500ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE), 100ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE), 10ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE), 5ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE), 300ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE), 150ULL);
    
    // 验证压缩率
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL); // 60/300*10000=2000
}

// 测试所有SAM字段的累加
TEST_F(PbgzStatTest, SamStatAllFieldsAccumulate) {
    SamStat samStat;
    samStat.init();
    
    // 对所有11个字段进行累加
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 25);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::ORIGINAL_SIZE, 50);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE, 10);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::ORIGINAL_SIZE, 75);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::COMPRESSED_SIZE, 15);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::ORIGINAL_SIZE, 80);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::COMPRESSED_SIZE, 16);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::ORIGINAL_SIZE, 40);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::COMPRESSED_SIZE, 8);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::ORIGINAL_SIZE, 120);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::COMPRESSED_SIZE, 60);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::ORIGINAL_SIZE, 75);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::COMPRESSED_SIZE, 15);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::ORIGINAL_SIZE, 80);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::COMPRESSED_SIZE, 32);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::ORIGINAL_SIZE, 60);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::COMPRESSED_SIZE, 48);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::ORIGINAL_SIZE, 500);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSED_SIZE, 250);
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatMetricIds::ORIGINAL_SIZE, 300);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatMetricIds::COMPRESSED_SIZE, 150);
    
    
    // 验证累加结果
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE), 25ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::ORIGINAL_SIZE), 50ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE), 10ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::ORIGINAL_SIZE), 75ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatMetricIds::COMPRESSED_SIZE), 15ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::ORIGINAL_SIZE), 80ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatMetricIds::COMPRESSED_SIZE), 16ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::ORIGINAL_SIZE), 40ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatMetricIds::COMPRESSED_SIZE), 8ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::ORIGINAL_SIZE), 120ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatMetricIds::COMPRESSED_SIZE), 60ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::ORIGINAL_SIZE), 75ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatMetricIds::COMPRESSED_SIZE), 15ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::ORIGINAL_SIZE), 80ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatMetricIds::COMPRESSED_SIZE), 32ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::ORIGINAL_SIZE), 60ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatMetricIds::COMPRESSED_SIZE), 48ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::ORIGINAL_SIZE), 500ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSED_SIZE), 250ULL);
    
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatMetricIds::ORIGINAL_SIZE), 300ULL);
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatMetricIds::COMPRESSED_SIZE), 150ULL);
    
    // 验证压缩率计算：QNAME像 25/100=25%, SEQ像 250/500=50%
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSION_RATIO), 2500ULL); // 25/100*10000
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSION_RATIO), 5000ULL);  // 250/500*10000
}

// 测试高压缩率的场景
TEST_F(PbgzStatTest, HighCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 设置高压缩率场景：1000字节压缩到50字节（5%压缩率）
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 50);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 500ULL); // 50/1000*10000=500
}

// 测试低压缩率的场景  
TEST_F(PbgzStatTest, LowCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 设置低压缩率场景：100字节压缩到90字节（90%压缩率）
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 90);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 9000ULL); // 90/100*10000=9000
}

// 测试大数据量的累加
TEST_F(PbgzStatTest, LargeValueAccumulation) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // 测试大数值的累加
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000000ULL);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 100000ULL);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 2000000ULL);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 150000ULL);
    
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 3000000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE), 250000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 833ULL); // 250000/3000000*10000≈833
}