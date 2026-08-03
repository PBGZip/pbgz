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
#include "coder/coder_fcv2.h"

namespace {

/*
 * Amount of raw line data scanned from the start of the first block.
 * 4MB is enough that the dominant fields (SEQ/QUAL ~ 37% each for typical
 * short reads) each yield a ~1.5MB sample, which is comfortably large enough
 * to separate the candidate coders.
 */
const uint32_t SAMPLE_TARGET = 4u << 20;   /* 4 MB */

/*
 * 先验训练最多使用 45 MB QUAL。实测到此处收益约为 1.84 个百分点，继续加到 90 MB
 * 只多约 0.10 个百分点；多出的时间却发生在切换到并行压缩之前的串行阶段，会直接
 * 推迟后续工作线程启动，所以必须在收益曲线变平前截断。
 */
const uint32_t QUAL_PRIOR_TRAIN_MAX = 45u * 1024u * 1024u;

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

/*
 * 领先者需要领先多少才算胜负已分。
 *
 * 取压缩后大小的相对差：领先者比次优者小 3% 以上就定案。差距小于这个数时，
 * 换更多样本重测有可能反转，值得再跑一轮；大于这个数基本不会反转。
 */
const double SETTLE_MARGIN = 0.03;

/* 评估的起步样本量。太小容易被局部数据带偏，取和最小可信样本一致。 */
const uint32_t PROBE_START = MIN_SELECT_SAMPLE;

FieldCodecSelection CodecSelector::selectCoder(const uint8_t* data, uint32_t len)
{
    FieldCodecSelection sel;
    sel.sampleLen = len;
    if (data == nullptr || len < MIN_SELECT_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    /*
     * 从小样本起步，每轮加倍，一旦领先者拉开足够差距就收手。
     *
     * 每一轮都是拿新实例重压 [0, probe)，而不是在上一轮的基础上增量喂。这样做有两个
     * 原因：coder_fc 只支持整块压缩，第二次调用 encode_line 会直接报错；而且重压
     * 恰好模拟了实际压缩时"一个数据块从头压到尾"的情形，量出来的数字更贴近真实表现。
     *
     * 重压带来的额外开销是有界的：样本量按 2 倍递增，所有轮次加起来不超过最后一轮的
     * 两倍。而多数字段在头一两轮就定案了，实际反而比一次性压满整个采样更省。
     */
    uint32_t bwtCmLen = 0, fcLen = 0;
    uint32_t bwtCmUs = 0, fcUs = 0, outUs = 0, outLen = 0;
    uint32_t bestLen = UINT32_MAX;
    CoderType bestCoder = CoderType::BWT_CM;
    bool anyOk = false;
    uint32_t probe = (PROBE_START < len) ? PROBE_START : len;

    while (true) {
        bwtCmLen = fcLen = 0;
        bestLen = UINT32_MAX;
        bestCoder = CoderType::BWT_CM;
        anyOk = false;
        uint32_t runnerUp = UINT32_MAX;
        sel.rounds++;

        if (trialEncode<coder_bwt_cm>(data, probe, pickBwtLevel(probe), outLen, outUs)) {
            bwtCmLen = outLen; bwtCmUs = outUs;
            bestLen = outLen; bestCoder = CoderType::BWT_CM;
            anyOk = true;
        }
        if (trialEncode<coder_fc>(data, probe, 0, outLen, outUs)) {
            fcLen = outLen; fcUs = outUs;
            if (outLen < bestLen) { runnerUp = bestLen; bestLen = outLen; bestCoder = CoderType::FC; }
            else { runnerUp = outLen; }
            anyOk = true;
        }

        if (!anyOk) {
            break;
        }
        /* 样本已经用尽，没有更多数据可加，只能就此定案。 */
        if (probe >= len) {
            break;
        }
        /* 只有一个候选能跑通，再加数据也没有比较对象。 */
        if (runnerUp == UINT32_MAX) {
            break;
        }
        /* 领先者已经拉开足够差距，再加数据不会反转。 */
        if ((double)(runnerUp - bestLen) / (double)runnerUp >= SETTLE_MARGIN) {
            break;
        }
        probe = (probe > len / 2) ? len : (probe * 2);
    }

    sel.decidedLen = probe;
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

void CodecSelector::extractQualSamples(RoughIOBlock* block,
                                       std::vector<QualSampleRecord>& records,
                                       std::vector<uint32_t>& freqByByte,
                                       uint32_t sampleBudget,
                                       uint64_t qualBudget)
{
    records.clear();
    freqByByte.assign(256, 0);
    SafeLineReader reader(block);
    uint64_t collectedQualBytes = 0;

    const uint8_t* line = nullptr;
    uint32_t lineLen = 0;
    while (reader.nextLine(line, lineLen)) {
        if (reader.scannedBytes() >= sampleBudget) {
            break;
        }
        if (lineLen == 0 || line[0] == '@') {
            continue;
        }

        /* 逐个 tab 切出前 11 个必需字段的起止位置。 */
        uint32_t beg[SAM_FIELD_COUNT] = {0};
        uint32_t end[SAM_FIELD_COUNT] = {0};
        uint32_t fieldIdx = 0;
        uint32_t pos = 0;
        while (pos <= lineLen && fieldIdx < SAM_FIELD_COUNT) {
            uint32_t tabPos = pos;
            while (tabPos < lineLen && line[tabPos] != '\t') {
                ++tabPos;
            }
            beg[fieldIdx] = pos;
            end[fieldIdx] = tabPos;
            ++fieldIdx;
            pos = tabPos + 1;
        }
        if (fieldIdx < SAM_FIELD_COUNT) {
            continue;   /* 字段不全的行跳过，不参与评估 */
        }

        uint32_t qBeg = beg[SAM_QUAL], qEnd = end[SAM_QUAL];
        if (qEnd <= qBeg) {
            continue;
        }
        /* QUAL 为单个 '*' 表示缺失，这类记录没有质量值可压。 */
        if (qEnd - qBeg == 1 && line[qBeg] == '*') {
            continue;
        }

        uint32_t qualLen = qEnd - qBeg;
        /*
         * 只接收完整记录，不能为了凑满上限截断 QUAL：fcv2 依赖真实读长来还原循环位置，
         * 截断会悄悄训练出错误的上下文。超过上限后不再保留记录，但仍扫描到块尾，保证
         * 训练路径处理的是完整首块而非选择阶段的 4 MB 原始样本。
         */
        if (collectedQualBytes + qualLen > qualBudget) {
            continue;
        }

        QualSampleRecord rec;
        rec.qual.assign((const char*)(line + qBeg), qualLen);
        rec.seq.assign((const char*)(line + beg[SAM_SEQ]), end[SAM_SEQ] - beg[SAM_SEQ]);

        /* 手工解析 FLAG，取 0x10 位。字段没有以 \0 结尾，不能直接用 strtol。 */
        long flagVal = 0;
        for (uint32_t p = beg[SAM_FLAG]; p < end[SAM_FLAG]; ++p) {
            uint8_t ch = line[p];
            if (ch < '0' || ch > '9') {
                break;
            }
            flagVal = flagVal * 10 + (ch - '0');
        }
        rec.rev = (flagVal & 16) != 0;

        for (uint32_t p = qBeg; p < qEnd; ++p) {
            freqByByte[line[p]]++;
        }
        records.push_back(rec);
        collectedQualBytes += qualLen;
    }
}

/*
 * 在首块的完整 QUAL 上训练 fcv2 先验。编码器选择仍只用小样本；小样本足以稳定比较
 * 候选，而模型学习需要更多质量值，复用它会把约 1.18 MB 的选择样本误当成训练集，损失
 * 大部分先验收益。这里故意只修改 PreprocessInfo：它是流水线与编码器决策的唯一接口，
 * 让 coder/ 层继续不知道块和文件偏移，维持既有分层。
 */
void CodecSelector::trainQualPrior(RoughIOBlock* block, PreprocessInfo& info)
{
    if (block == nullptr || block->getDataLen() <= 0) {
        return;
    }

    try {
        std::vector<QualSampleRecord> records;
        std::vector<uint32_t> freqByByte;
        extractQualSamples(block, records, freqByByte, UINT32_MAX, QUAL_PRIOR_TRAIN_MAX);

        uint64_t trainedBytes = 0;
        for (size_t i = 0; i < records.size(); ++i) {
            trainedBytes += records[i].qual.size();
        }
        if (trainedBytes == 0) {
            return;
        }

        std::vector<uint8_t> scratch((size_t)trainedBytes * 2 + (1u << 16), 0);
        coder_io io(scratch.data(), (int32_t)scratch.size());
        coder_fcv2 coder(&io, freqByByte);
        for (size_t i = 0; i < records.size(); ++i) {
            const QualSampleRecord& record = records[i];
            coder.encode_record((const uint8_t*)record.qual.data(),
                                (uint32_t)record.qual.size(), record.rev);
        }
        if (coder.encode_flush() <= 0) {
            return;
        }

        std::vector<uint8_t> snapshot;
        if (!coder.export_model(snapshot) || snapshot.empty()) {
            return;
        }
        info.setQualPrior(std::move(snapshot), trainedBytes);
    } catch (...) {
        /* 训练只是可选优化；任意分配或未知数据失败都必须回到无先验的既有压缩路径。 */
    }
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
        if (f == (uint32_t)SAM_QUAL) {
            /*
             * 质量值列走专用评估：候选是 coder_qual 和 fcv2，而不是通用的字节流
             * 压缩器。此前这里比较的是 coder_bwt_cm 和 coder_fc，而压缩时用的是
             * 另外两个，选出来的结果实际上没有被使用。
             */
            std::vector<QualSampleRecord> qualRecs;
            std::vector<uint32_t> qualFreq;
            extractQualSamples(block, qualRecs, qualFreq, SAMPLE_TARGET);
            info.fields[f] = QualSelector::select(qualRecs, qualFreq);
            if (info.fields[f].status == FieldStatus::SELECTED &&
                info.fields[f].selectedCoder == CoderType::FCV2) {
                trainQualPrior(block, info);
            }
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
