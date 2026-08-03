/*
 * codec_selector.cpp - Codec pre-selection for the file preprocessing module
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

#include "codec_selector.h"

#include <chrono>

#include <climits>
#include <string>
#include <vector>

#include "io_block.h"
#include "block_wrapper.h"
#include "log/logger.h"
#include "safe_line_reader.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"

namespace {

/*
 * Amount of raw line data scanned from the start of the first block.
 * 4MB is enough that the dominant fields (SEQ/QUAL ~ 37% each for typical
 * short reads) each yield a ~1.5MB sample, which is comfortably large enough
 * to separate the candidate coders.
 */
const uint32_t SAMPLE_TARGET = 4u << 20;   /* 4 MB */

/*
 * Minimum per-field sample required before we trust a codec comparison.
 * Each coder carries a small fixed overhead (BWT index, range-coder flush,
 * EOF marker, on the order of tens of bytes). At 64KB a real 1.7pp coder
 * difference is ~1.1KB, i.e. an order of magnitude above that overhead, so
 * the comparison is meaningful. Fields with a smaller sample are SKIPPED and
 * keep their default coder.
 */
const uint32_t MIN_SELECT_SAMPLE = 64u << 10;   /* 64 KB */

/* coder_bwt_cm internal block sizes per level, mirroring coder_bwt_cm. */
const uint32_t BWT_LEVEL_SIZE[10] = {
    0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
    1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
};

/*
 * Trial-compress a byte stream with one coder type and report the compressed
 * size. All coders follow the same encode_line()/encode_flush() contract and
 * write into coder_io::data, so a single template covers them.
 *
 * The output buffer is sized generously (2x input + margin): arithmetic coding
 * never expands much, but coder_fc writes up to 2x input internally, and the
 * coder_io sink is not bounds-checked.
 */
template <typename CoderT>
bool trialEncode(const uint8_t* data, uint32_t len, int bwtLevel, uint32_t& outLen, uint32_t& usec)
{
    const auto t0 = std::chrono::steady_clock::now();
    uint32_t cap = (len << 1) + 65536;
    std::vector<uint8_t> outBuf(cap, 0);
    coder_io io(outBuf.data(), (int32_t)cap);
    if (bwtLevel > 0) {
        io.set_level(bwtLevel);   /* only coder_bwt_cm consumes this */
    }
    {
        CoderT coder(&io);
        coder.encode_line(data, len);
        coder.encode_flush();
    }
    outLen = (uint32_t)io.data_len;
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return (outLen > 0 && outLen < cap);
}

} // namespace

int CodecSelector::pickBwtLevel(uint32_t sampleLen)
{
    /* Smallest level whose block size fits the whole sample in one block. */
    for (int lv = 1; lv <= 9; ++lv) {
        if (BWT_LEVEL_SIZE[lv] >= sampleLen) {
            return lv;
        }
    }
    return 9;
}

FieldCodecSelection CodecSelector::selectCoder(const uint8_t* data, uint32_t len)
{
    FieldCodecSelection sel;
    sel.sampleLen = len;
    if (data == nullptr || len < MIN_SELECT_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    uint32_t bestLen = UINT32_MAX;
    CoderType bestCoder = CoderType::BWT_CM;
    bool anyOk = false;
    uint32_t outLen = 0;

    /* coder_bwt_cm: use a block size that fits the sample in one block. */
    /*
     * coder_simple_rc 已从候选中移除：实测它是有损的，往返校验全部失败。
     * 试压只比较压缩后的大小、不验证能否原样解回来，留着它迟早会选出一个解不出
     * 原始数据的编码器。simpleRcLen 保留为 0，只是为了不改动 -v 的打印结构。
     */
    uint32_t bwtCmLen = 0, fcLen = 0;
    uint32_t bwtCmUs = 0, fcUs = 0, outUs = 0;

    if (trialEncode<coder_bwt_cm>(data, len, pickBwtLevel(len), outLen, outUs)) {
        bwtCmLen = outLen; bwtCmUs = outUs;
        if (outLen < bestLen) { bestLen = outLen; bestCoder = CoderType::BWT_CM; }
        anyOk = true;
    }

    if (trialEncode<coder_fc>(data, len, 0, outLen, outUs)) {
        fcLen = outLen; fcUs = outUs;
        if (outLen < bestLen) { bestLen = outLen; bestCoder = CoderType::FC; }
        anyOk = true;
    }

    sel.trialBwtCmLen = bwtCmLen;
    sel.trialFcLen = fcLen;
    sel.trialBwtCmUs = bwtCmUs;
    sel.trialFcUs = fcUs;

    if (!anyOk) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    sel.status = FieldStatus::SELECTED;
    sel.selectedCoder = bestCoder;
    sel.bestCompLen = bestLen;
    return sel;
}

uint32_t CodecSelector::extractSamFieldSamples(RoughIOBlock* block,
                                           std::vector<std::string>& fieldBufs,
                                           uint32_t sampleBudget)
{
    fieldBufs.assign(SAM_FIELD_COUNT, std::string());
    SafeLineReader reader(block);

    const uint8_t* line = nullptr;
    uint32_t lineLen = 0;
    while (reader.nextLine(line, lineLen)) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }
        if (lineLen == 0 || line[0] == '@') {
            continue;
        }

        uint32_t fieldIdx = 0;
        uint32_t pos = 0;
        while (pos <= lineLen && fieldIdx < SAM_FIELD_COUNT) {
            uint32_t tabPos = pos;
            while (tabPos < lineLen && line[tabPos] != '\t') {
                ++tabPos;
            }
            fieldBufs[fieldIdx].append((const char*)(line + pos), (size_t)(tabPos - pos));
            ++fieldIdx;
            pos = tabPos + 1;
        }
    }
    return reader.scannedBytes();
}

uint32_t CodecSelector::extractFastqFieldSamples(RoughIOBlock* block,
                                             std::vector<std::string>& fieldBufs,
                                             uint32_t sampleBudget)
{
    fieldBufs.assign(FQ_FIELD_COUNT, std::string());
    SafeLineReader reader(block);

    const uint8_t* lines[4];
    uint32_t lens[4];
    uint32_t slot = 0;

    while (reader.nextLine(lines[slot], lens[slot])) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }

        ++slot;
        if (slot < 4) {
            continue;
        }

        /* A complete 4-line FASTQ record is now buffered in lines[0..3]. */
        const uint8_t* id = lines[0];
        const uint8_t* comment = lines[2];
        if (lens[0] > 0 && lens[2] > 0 && id[0] == '@' && comment[0] == '+') {
            if (lens[0] > 1) fieldBufs[FQ_ID].append((const char*)(id + 1), (size_t)(lens[0] - 1));
            if (lens[1] > 0) fieldBufs[FQ_SEQ].append((const char*)lines[1], (size_t)lens[1]);
            if (lens[3] > 0) fieldBufs[FQ_QUAL].append((const char*)lines[3], (size_t)lens[3]);
            if (lens[2] > 1) fieldBufs[FQ_COMMENT].append((const char*)(comment + 1), (size_t)(lens[2] - 1));
        }

        slot = 0;
    }
    return reader.scannedBytes();
}

int32_t CodecSelector::analyzeSam(RoughIOBlock* block, PreprocessInfo& info)
{
    std::vector<std::string> fieldBufs;
    info.scannedBytes = extractSamFieldSamples(block, fieldBufs, SAMPLE_TARGET);

    info.fields.resize(SAM_FIELD_COUNT);
    uint64_t totalSample = 0;
    for (uint32_t f = 0; f < SAM_FIELD_COUNT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();
        if (buf.size() < MIN_SELECT_SAMPLE) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }
        info.fields[f] = selectCoder((const uint8_t*)buf.data(), (uint32_t)buf.size());
        LOG_DEBUG("Preprocess SAM field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }
    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyzeFastq(RoughIOBlock* block, PreprocessInfo& info)
{
    std::vector<std::string> fieldBufs;
    info.scannedBytes = extractFastqFieldSamples(block, fieldBufs, SAMPLE_TARGET);

    info.fields.resize(FQ_FIELD_COUNT);
    uint64_t totalSample = 0;
    for (uint32_t f = 0; f < FQ_FIELD_COUNT; ++f) {
        const std::string& buf = fieldBufs[f];
        totalSample += buf.size();
        if (buf.size() < MIN_SELECT_SAMPLE) {
            info.fields[f].status = FieldStatus::SKIPPED;
            info.fields[f].sampleLen = (uint32_t)buf.size();
            continue;
        }
        info.fields[f] = selectCoder((const uint8_t*)buf.data(), (uint32_t)buf.size());
        LOG_DEBUG("Preprocess FASTQ field %u: sample=%u, coder=%s, comp=%u (%.2f%%)",
                  f, info.fields[f].sampleLen, coderTypeToMagic(info.fields[f].selectedCoder),
                  info.fields[f].bestCompLen, info.fields[f].ratio() * 100.0);
    }
    info.sampleBytes = (uint32_t)totalSample;
    return 0;
}

int32_t CodecSelector::analyze(RoughIOBlock* block, PreprocessInfo& info)
{
    if (block == nullptr) {
        return -1;
    }
    BlockType type = block->getBlockType();
    info.reset(type);

    if (BlockUtil::isSAMBlock(type)) {
        return analyzeSam(block, info);
    }
    if (BlockUtil::isFastqBlock(type)) {
        return analyzeFastq(block, info);
    }

    /* Unsupported type: leave info empty; actuators use their defaults. */
    return 0;
}
