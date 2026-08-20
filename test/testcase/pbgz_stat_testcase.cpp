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
        // Setup before the test
    }
    
    void TearDown() override {
        // Cleanup after the test
    }
};

// Test FastqStat initialization
TEST_F(PbgzStatTest, FastqStatInitialization) {
    FastqStat fastqStat;
    int result = fastqStat.init();
    EXPECT_EQ(result, 0);
}

// Test SamStat initialization
TEST_F(PbgzStatTest, SamStatInitialization) {
    SamStat samStat;
    int result = samStat.init();
    EXPECT_EQ(result, 0);
}

// Test FastqStat setting metric values
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

// Test accumulating metric values
TEST_F(PbgzStatTest, FastqStatAccumulateValue) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 10);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 20);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 15);
    
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 45ULL);
}

// Test compression ratio calculation
TEST_F(PbgzStatTest, CompressionRatioCalculation) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    // ratio = (200 / 1000) * 10000 = 2000 (i.e., 20%)
    EXPECT_EQ(ratio, 2000ULL);
}

// Test compression ratio calculation when the original size is 0
TEST_F(PbgzStatTest, CompressionRatioZeroOriginal) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 0);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO);
    // With an original size of 0, the compression ratio should be 0
    EXPECT_EQ(ratio, 0ULL);
}

// Test SamStat setting metric values
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

// Test SamStat accumulating metric values
TEST_F(PbgzStatTest, SamStatAccumulateValue) {
    SamStat samStat;
    samStat.init();
    
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE, 5);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE, 10);
    
    uint64_t value = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 15ULL);
}

// Test fetching a non-existent metric
TEST_F(PbgzStatTest, GetNonExistentMetric) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Fetch a metric that was never set
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 0ULL);
}

// Test global metric setting
TEST_F(PbgzStatTest, GlobalMetrics) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Note: call calculateAllMetrics first to ensure the compression ratio is computed
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    
    // Now test that the compression ratio is computed correctly
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 2000ULL);  // (200/1000)*10000 = 2000
}

// Test field name storage
TEST_F(PbgzStatTest, FieldNamesStorage) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Since printStats() now outputs to stderr via fprintf, we can call it directly to verify functionality
    // This test mainly verifies that printStats() does not crash and produces output normally
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    // printStats() computes and outputs the statistics
    fastqStat.printStats();
    
    // Verify the statistics are computed correctly
    uint64_t originalSize = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t compressedSize = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(originalSize, 100ULL);
    EXPECT_EQ(compressedSize, 20ULL);
}

// Test namespace access
TEST_F(PbgzStatTest, NamespaceAccess) {
    // Test the StatUnitId namespace
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
    
    // Test the StatMetricIds namespace
    EXPECT_EQ(StatMetricIds::ORIGINAL_SIZE, 1);
    EXPECT_EQ(StatMetricIds::COMPRESSED_SIZE, 2);
    EXPECT_EQ(StatMetricIds::COMPRESSION_RATIO, 3);
    
    // Test the StatUnitNames strings
    EXPECT_EQ(StatObjectNames::Fastq::ID, "ID");
    EXPECT_EQ(StatObjectNames::Fastq::BASE, "Base");
    EXPECT_EQ(StatObjectNames::Fastq::COMMENT, "Comment");
    EXPECT_EQ(StatObjectNames::Fastq::QUALITY, "Quality");
    
    EXPECT_EQ(StatObjectNames::Sam::QNAME, "QNAME");
    EXPECT_EQ(StatObjectNames::Sam::FLAG, "FLAG");
    EXPECT_EQ(StatObjectNames::Sam::SEQ, "SEQ");
    EXPECT_EQ(StatObjectNames::Sam::QUAL, "QUAL");
}

// Test statistics across multiple fields
TEST_F(PbgzStatTest, MultipleFields) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Set data for multiple fields
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 500);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE, 10);
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE, 40);
    
    // Compute the compression ratio
    
    // Verify each field's compression ratio is computed correctly
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
}

// Test computation of different compression ratios
TEST_F(PbgzStatTest, DifferentCompressionRatios) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Set different compression ratios
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

// Test recomputing the compression ratio
TEST_F(PbgzStatTest, RecalculateCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Initial setup
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL);
    
    // Update the compressed size and recompute
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 100);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 3000ULL);  // (300/1000)*10000
}

// Test accumulation across multiple metric types
TEST_F(PbgzStatTest, AccumulationAndSetting) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // First set, then add
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 50);
    
    // ORIGINAL_SIZE is of type ACCUMULATED, so add accumulates
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 150ULL);  // 100 + 50 = 150
    
    // COMPRESSED_SIZE is of type ACCUMULATED, so add takes effect
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 30);
    
    value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE);
    EXPECT_EQ(value, 180ULL);
}

// Test the update order of computed metrics
TEST_F(PbgzStatTest, CalculatedMetricUpdateOrder) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Set only the original size
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    
    // The ratio should be 0 because the original size is non-zero but the compressed size is 0
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 0ULL);
    
    // Set the compressed size, then recompute
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 200);
    
    // Now the compression ratio should be correct
    ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 2000ULL);
}

// Test independent statistic objects
TEST_F(PbgzStatTest, IndependentStatObjects) {
    FastqStat fastqStat1;
    FastqStat fastqStat2;
    
    fastqStat1.init();
    fastqStat2.init();
    
    // Different objects do not affect each other's data
    fastqStat1.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat1.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat2.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat2.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 40);
    
    
    // Verify the two objects' data are independent
    EXPECT_EQ(fastqStat1.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat2.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 200ULL);
}

// Test repeated initialization
TEST_F(PbgzStatTest, ReInitialization) {
    FastqStat fastqStat;
    
    fastqStat.init();
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    
    // Re-initialization should clear the previous settings
    fastqStat.init();
    
    uint64_t value = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    EXPECT_EQ(value, 0ULL);
}

// Test statistics for all fields
TEST_F(PbgzStatTest, AllFastqFields) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Verify all 4 FASTQ fields are registered correctly
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 100);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE), 100ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE), 100ULL);
}

// Test statistics for all SAM fields
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

// Test FastqStat accumulating data multiple times
TEST_F(PbgzStatTest, FastqStatMultipleAccumulations) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // First accumulation on the ID field
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 50);
    
    // Second accumulation on the ID field
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 150);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 30);
    
    // Verify the accumulated results: 100+150=250, 50+30=80
    uint64_t totalOriginal = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 250ULL);
    EXPECT_EQ(totalCompressed, 80ULL);
}

// Test SamStat accumulating data multiple times
TEST_F(PbgzStatTest, SamStatMultipleAccumulations) {
    SamStat samStat;
    samStat.init();
    
    // First accumulation on the QNAME field
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 100);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 25);
    
    // Second accumulation on the QNAME field
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 200);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 45);
    
    // Third accumulation on the QNAME field
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE, 150);
    samStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE, 60);
    
    // Verify the accumulated results: 100+200+150=450, 25+45+60=130
    uint64_t totalOriginal = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 450ULL);
    EXPECT_EQ(totalCompressed, 130ULL);
}

// Test mixing set and add
TEST_F(PbgzStatTest, MixedSetAndAccumulate) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // First set the initial values
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 80);
    
    // Then accumulate
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 50);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    // Finally set again, but typically add is used for dynamic accumulation and set for initialization
    // Here we verify the result is set + add = 150, 100
    uint64_t totalOriginal = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE);
    uint64_t totalCompressed = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE);
    
    EXPECT_EQ(totalOriginal, 150ULL);
    EXPECT_EQ(totalCompressed, 100ULL);
}

// Test the compression ratio for empty fields
TEST_F(PbgzStatTest, EmptyFieldCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 0);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 0);
    
    
    // When the original size is 0, the compression ratio should be 0
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 0ULL);
}

// Test accumulation across all FASTQ fields
TEST_F(PbgzStatTest, FastqStatAllFieldsAccumulate) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Accumulate multiple times across all 4 fields
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 20);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 500);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 100);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE, 10);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE, 5);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE, 300);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE, 150);
    
    // Second accumulation on the ID field
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 200);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 40);
    
    
    // Verify the accumulated results
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 300ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE), 60ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE), 500ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE), 100ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::ORIGINAL_SIZE), 10ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatMetricIds::COMPRESSED_SIZE), 5ULL);
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::ORIGINAL_SIZE), 300ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatMetricIds::COMPRESSED_SIZE), 150ULL);
    
    // Verify the compression ratios
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 2000ULL); // 60/300*10000=2000
}

// Test accumulation across all SAM fields
TEST_F(PbgzStatTest, SamStatAllFieldsAccumulate) {
    SamStat samStat;
    samStat.init();
    
    // Accumulate across all 11 fields
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
    
    
    // Verify the accumulated results
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
    
    // Verify the ratio computations: QNAME like 25/100=25%, SEQ like 250/500=50%
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatMetricIds::COMPRESSION_RATIO), 2500ULL); // 25/100*10000
    EXPECT_EQ(samStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatMetricIds::COMPRESSION_RATIO), 5000ULL);  // 250/500*10000
}

// Test a high compression ratio scenario
TEST_F(PbgzStatTest, HighCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Set up a high compression ratio scenario: 1000 bytes compressed to 50 bytes (5% ratio)
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 50);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 500ULL); // 50/1000*10000=500
}

// Test a low compression ratio scenario  
TEST_F(PbgzStatTest, LowCompressionRatio) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Set up a low compression ratio scenario: 100 bytes compressed to 90 bytes (90% ratio)
    fastqStat.setMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::ORIGINAL_SIZE, 100);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSED_SIZE, 90);
    
    
    uint64_t ratio = fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatMetricIds::COMPRESSION_RATIO);
    EXPECT_EQ(ratio, 9000ULL); // 90/100*10000=9000
}

// Test accumulation with large data volumes
TEST_F(PbgzStatTest, LargeValueAccumulation) {
    FastqStat fastqStat;
    fastqStat.init();
    
    // Test accumulation of large values
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 1000000ULL);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 100000ULL);
    
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE, 2000000ULL);
    fastqStat.addMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE, 150000ULL);
    
    
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::ORIGINAL_SIZE), 3000000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSED_SIZE), 250000ULL);
    EXPECT_EQ(fastqStat.getMetricValue(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatMetricIds::COMPRESSION_RATIO), 833ULL); // 250000/3000000*10000≈833
}