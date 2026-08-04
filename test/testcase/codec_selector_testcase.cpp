/*
 * codec_selector_testcase.cpp - Tests for the file preprocessing codec selector
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

#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "codec_selector.h"
#include "preprocess_info.h"
#include "io_block.h"
#include "coder.h"
#include "utils/memory_util.h"

namespace {

/*
 * Fill n quality bytes mimicking real Illumina QUAL: values sit in a narrow
 * high-quality range and tend to REPEAT, forming runs of identical bytes.
 * Those runs are what make the BWT->MTF output heavily skewed toward symbol 0,
 * which is exactly the distribution coder_simple_rc compresses best.
 */
void fillQuality(uint8_t* data, uint32_t n, std::mt19937& rng)
{
    std::uniform_int_distribution<int> valDist(28, 40);
    std::uniform_int_distribution<int> repeatDist(0, 99);
    int q = valDist(rng);
    for (uint32_t i = 0; i < n; ++i) {
        if (repeatDist(rng) >= 18) {   /* ~82% chance: extend the current run */
            q = valDist(rng);
        }
        data[i] = (uint8_t)(q + 33);
    }
}

std::vector<uint8_t> makeQualityData(uint32_t n, uint32_t seed)
{
    std::vector<uint8_t> data(n);
    std::mt19937 rng(seed);
    fillQuality(data.data(), n, rng);
    return data;
}

void appendQuality(std::string& out, uint32_t n, std::mt19937& rng)
{
    std::vector<uint8_t> tmp(n);
    fillQuality(tmp.data(), n, rng);
    out.append((const char*)tmp.data(), n);
}

/* Locate the real SAM fixture relative to the test working directory. */
std::string findRealSamFile()
{
    const char* candidates[] = {
        "test/con_sorted.sam",
        "../test/con_sorted.sam",
        "../../test/con_sorted.sam",
    };
    for (const char* path : candidates) {
        std::ifstream probe(path);
        if (probe.good()) {
            return path;
        }
    }
    return "";
}

/* Concatenate the QUAL column (11th SAM field) from the first data lines. */
std::vector<uint8_t> extractRealQuality(const std::string& samPath, uint32_t targetBytes)
{
    std::vector<uint8_t> qual;
    std::ifstream in(samPath);
    if (!in) {
        return qual;
    }
    std::string line;
    while (qual.size() < targetBytes && std::getline(in, line)) {
        if (line.empty() || line[0] == '@') {
            continue;
        }
        uint32_t fieldIdx = 0;
        size_t fieldStart = 0;
        for (size_t pos = 0; pos <= line.size() && fieldIdx <= SAM_QUAL; ++pos) {
            if (pos == line.size() || line[pos] == '\t') {
                if (fieldIdx == SAM_QUAL && pos > fieldStart) {
                    qual.insert(qual.end(), line.begin() + fieldStart, line.begin() + pos);
                }
                ++fieldIdx;
                fieldStart = pos + 1;
            }
        }
    }
    return qual;
}

/* Build a block of synthetic SAM alignment lines (11 mandatory fields). */
void buildSamBlock(RoughIOBlock* block, uint32_t numLines, uint32_t readLen)
{
    std::vector<size_t>& npos = block->getNpos();
    uint8_t* buffer = block->getBuffer();
    uint32_t offset = 0;
    std::mt19937 rng(1234);

    for (uint32_t i = 0; i < numLines; ++i) {
        std::string line;
        line += "read" + std::to_string(i);
        line += "\t99\tchrI\t";
        line += std::to_string(i * 100 + 1);
        line += "\t60\t";
        line += std::to_string(readLen) + "M\t";
        line += "=\t";
        line += std::to_string(i * 100 + 100);
        line += "\t200\t";
        for (uint32_t j = 0; j < readLen; ++j) {
            line += "ACGT"[rng() % 4];
        }
        line += "\t";
        appendQuality(line, readLen, rng);
        line += "\n";

        memcpy(buffer + offset, line.data(), line.size());
        offset += (uint32_t)line.size();
        npos.push_back(offset - 1);
    }
    block->setDataLen(offset);
    block->setBlockType(SAM);
}

/* Build a block of synthetic FASTQ records (4 lines each). */
void buildFastqBlock(RoughIOBlock* block, uint32_t numRecords, uint32_t readLen)
{
    std::vector<size_t>& npos = block->getNpos();
    uint8_t* buffer = block->getBuffer();
    uint32_t offset = 0;
    std::mt19937 rng(5678);

    auto appendLine = [&](const std::string& s) {
        memcpy(buffer + offset, s.data(), s.size());
        offset += (uint32_t)s.size();
        buffer[offset++] = '\n';
        npos.push_back(offset - 1);
    };

    for (uint32_t i = 0; i < numRecords; ++i) {
        appendLine("@read" + std::to_string(i));

        std::string seq;
        for (uint32_t j = 0; j < readLen; ++j) seq += "ACGT"[rng() % 4];
        appendLine(seq);

        appendLine("+");

        std::string qual;
        appendQuality(qual, readLen, rng);
        appendLine(qual);
    }
    block->setDataLen(offset);
    block->setBlockType(FASTQ_GEN2);
}

} // namespace

class CodecSelectorTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();
    }
};

TEST_F(CodecSelectorTest, SmallSampleIsSkipped)
{
    std::vector<uint8_t> data(1024, 'A');
    FieldCodecSelection sel = CodecSelector::selectCoder(data.data(), (uint32_t)data.size());
    EXPECT_EQ(sel.status, FieldStatus::SKIPPED);
}

TEST_F(CodecSelectorTest, NullDataIsSkipped)
{
    FieldCodecSelection sel = CodecSelector::selectCoder(nullptr, 1 << 20);
    EXPECT_EQ(sel.status, FieldStatus::SKIPPED);
}

TEST_F(CodecSelectorTest, QualitySampleIsSelectedAndCompresses)
{
    std::vector<uint8_t> data = makeQualityData(256 << 10, 42);
    FieldCodecSelection sel = CodecSelector::selectCoder(data.data(), (uint32_t)data.size());
    EXPECT_EQ(sel.status, FieldStatus::SELECTED);
    EXPECT_GT(sel.bestCompLen, 0u);
    EXPECT_LT(sel.bestCompLen, sel.sampleLen);
}

/*
 * Measured on 1MB of real QUAL extracted from con_sorted.sam:
 *   bwt_cm 26.60% < fc 27.23% < simple_rc 32.46%.
 * So the selector must pick bwt_cm here; that ratio also beats CRAM FQZ
 * (~28.7%), which is the compression target for the QUAL field.
 */
TEST_F(CodecSelectorTest, RealQualityPrefersBwtCm)
{
    std::string samPath = findRealSamFile();
    if (samPath.empty()) {
        GTEST_SKIP() << "con_sorted.sam fixture not found";
    }
    std::vector<uint8_t> qual = extractRealQuality(samPath, 1 << 20);
    ASSERT_GE(qual.size(), (size_t)(64u << 10)) << "not enough QUAL data extracted";

    FieldCodecSelection sel = CodecSelector::selectCoder(qual.data(), (uint32_t)qual.size());
    ASSERT_EQ(sel.status, FieldStatus::SELECTED);
    EXPECT_EQ(sel.selectedCoder, CoderType::BWT_CM);
    EXPECT_LT(sel.ratio(), 0.29) << "selected coder should beat CRAM FQZ on QUAL";
}

TEST_F(CodecSelectorTest, AnalyzeNullBlockFails)
{
    PreprocessInfo info;
    EXPECT_EQ(CodecSelector::analyze(nullptr, 0, info), -1);
}

TEST_F(CodecSelectorTest, AnalyzeUnsupportedTypeIsNoOp)
{
    RoughIOBlock block(1 << 20);
    block.setBlockType(BINARY);
    PreprocessInfo info;
    EXPECT_EQ(CodecSelector::analyze(&block, 0, info), 0);
    EXPECT_FALSE(info.isDone());
    EXPECT_TRUE(info.fields.empty());
}

TEST_F(CodecSelectorTest, AnalyzeSamBlockSelectsLargeFields)
{
    RoughIOBlock block(1 << 20);
    buildSamBlock(&block, 1000, 90);

    PreprocessInfo info;
    ASSERT_EQ(CodecSelector::analyze(&block, 0, info), 0);
    info.markDone();
    EXPECT_TRUE(info.isDone());
    EXPECT_EQ(info.fileType, SAM);
    ASSERT_EQ(info.fields.size(), (size_t)SAM_FIELD_COUNT);

    EXPECT_EQ(info.fields[SAM_QUAL].status, FieldStatus::SELECTED);
    EXPECT_EQ(info.fields[SAM_SEQ].status, FieldStatus::SELECTED);
    EXPECT_LT(info.fields[SAM_QUAL].bestCompLen, info.fields[SAM_QUAL].sampleLen);

    EXPECT_EQ(info.fields[SAM_FLAG].status, FieldStatus::SKIPPED);
    EXPECT_EQ(info.fields[SAM_MAPQ].status, FieldStatus::SKIPPED);
}

TEST_F(CodecSelectorTest, AnalyzeFastqBlockSelectsSeqAndQual)
{
    RoughIOBlock block(1 << 20);
    buildFastqBlock(&block, 1000, 90);

    PreprocessInfo info;
    ASSERT_EQ(CodecSelector::analyze(&block, 0, info), 0);
    info.markDone();
    EXPECT_TRUE(info.isDone());
    ASSERT_EQ(info.fields.size(), (size_t)FQ_FIELD_COUNT);

    EXPECT_EQ(info.fields[FQ_SEQ].status, FieldStatus::SELECTED);
    EXPECT_EQ(info.fields[FQ_QUAL].status, FieldStatus::SELECTED);
    EXPECT_LT(info.fields[FQ_QUAL].bestCompLen, info.fields[FQ_QUAL].sampleLen);
}

/*
 * 先验是一次固定支出，只有 QUAL 总量足够大时才摊得回来。小输入必须不训练，
 * 否则每个小文件都白白多背一个辅助块。
 */
TEST_F(CodecSelectorTest, QualPriorSkippedForSmallInput)
{
    RoughIOBlock block(1 << 20);
    buildSamBlock(&block, 1000, 90);

    PreprocessInfo info;
    ASSERT_EQ(CodecSelector::analyze(&block, 1ull << 20, info), 0);
    EXPECT_TRUE(info.qualPrior().empty());
    EXPECT_EQ(info.qualPriorTrainedBytes(), 0u);
}

/* 输入长度不可知（管道）时按"写"处理：漏写的损失没有上界，误写的成本有。 */
TEST_F(CodecSelectorTest, QualPriorKeptWhenInputSizeUnknown)
{
    RoughIOBlock blockSmall(1 << 20);
    buildSamBlock(&blockSmall, 1000, 90);
    PreprocessInfo gated;
    ASSERT_EQ(CodecSelector::analyze(&blockSmall, 1ull << 20, gated), 0);

    RoughIOBlock block(1 << 20);
    buildSamBlock(&block, 1000, 90);
    PreprocessInfo info;
    ASSERT_EQ(CodecSelector::analyze(&block, 0, info), 0);

    if (info.fields[SAM_QUAL].selectedCoder == CoderType::FCV2) {
        EXPECT_FALSE(info.qualPrior().empty());
        EXPECT_TRUE(gated.qualPrior().empty());
    }
}

TEST_F(CodecSelectorTest, CoderForFallsBackWhenNotSelected)
{
    PreprocessInfo info;
    info.fields.resize(SAM_FIELD_COUNT);
    info.fields[SAM_QUAL].status = FieldStatus::SELECTED;
    info.fields[SAM_QUAL].selectedCoder = CoderType::SIMPLE_RC;
    info.fields[SAM_SEQ].status = FieldStatus::SKIPPED;
    info.markDone();

    EXPECT_EQ(info.coderFor(SAM_QUAL, CoderType::BWT_CM), CoderType::SIMPLE_RC);
    EXPECT_EQ(info.coderFor(SAM_SEQ, CoderType::BWT_CM), CoderType::BWT_CM);
    EXPECT_EQ(info.coderFor(999, CoderType::FC), CoderType::FC);
}
