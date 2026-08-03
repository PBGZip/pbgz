/*
 * qual_selector.cpp - 质量值列的编码器评估
 *
 * 本文件刻意不包含 coder_fc.h 或任何会引入 fc/rangecoder.h 的头文件：那里的
 * RangeCoder 与 coder_qual.h 间接引入的 clr.h 中的同名类冲突。质量值的评估之所以
 * 从 codec_selector.cpp 里分出来，这是原因之一，详见 qual_selector.h 的说明。
 */

#include "qual_selector.h"

#include <string.h>
#include <chrono>
#include <memory>

#include "coder/coder_io.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_bwt_cm.h"
#include "log/logger.h"

namespace {

/* 样本少于这个字节数就不做评估，直接沿用默认编码器。太少的样本判断不可靠。 */
const uint32_t MIN_QUAL_SAMPLE = 64u << 10;

/* 试压输出缓冲的富余量。质量值几乎不可能压不动，留一倍余量足够安全。 */
inline size_t trialCapacity(size_t srcLen)
{
    return (srcLen << 1) + (1u << 16);
}

/*
 * 用 coder_qual 试压。
 *
 * 它以对应的碱基序列作为上下文，所以必须逐条记录喂入，同时给出 seq 和 qual。
 * 频率表的格式沿用 sam_actuator 的既有约定：按出现次数降序排列的字母表，
 * 第二个元素固定为 1（实际计数没有被保留下来）。
 */
bool trialQual(const std::vector<QualSampleRecord>& records,
               const std::vector<uint32_t>& freqByByte,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < records.size(); i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<std::pair<uint16_t, uint16_t>> freqTable;
    for (uint32_t b = 0; b < 256; b++) {
        if (b < freqByByte.size() && freqByByte[b] > 0) {
            freqTable.push_back(std::make_pair((uint16_t)(b - '!'), (uint16_t)1));
        }
    }
    if (freqTable.empty()) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_qual coder(&io, true, freqTable);
        for (size_t i = 0; i < records.size(); i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_qual_gen2((uint8_t*)r.seq.data(),
                                   (uint8_t*)r.qual.data(),
                                   (uint32_t)r.qual.size());
        }
        coder.encode_flush();
        outLen = (uint32_t)io.data_len;
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/*
 * 用 fcv2 试压。
 *
 * 它不需要碱基序列，但需要每条记录的长度和链方向。链方向由它自己写进码流，
 * 所以这里传进去之后解码端不必再提供一次。
 */
bool trialFcv2(const std::vector<QualSampleRecord>& records,
               const std::vector<uint32_t>& freqByByte,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < records.size(); i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_fcv2 coder(&io, freqByByte);
        for (size_t i = 0; i < records.size(); i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_record((const uint8_t*)r.qual.data(),
                                (uint32_t)r.qual.size(), r.rev);
        }
        outLen = (uint32_t)coder.encode_flush();
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

/*
 * 用 coder_bwt_cm 试压。
 *
 * 这个候选原先不在质量值列的候选集里——通用字段走 CodecSelector 才会试 bwt_cm，
 * 质量值列只在 coder_qual 和 fcv2 之间二选一。实测发现这是个真实缺口：4 MB 真实
 * 质量值上 fcv2 25.81%、bwt_cm 26.31%、coder_qual 33.74%，fcv2 依然胜出，
 * 但 bwt_cm 比 coder_qual 好 7.4 个百分点，而 coder_qual 恰恰是当前的兜底选择。
 *
 * fcv2 有明确的不适用场景：它需要每条记录的长度和链方向，只有比对后 SAM 的 QUAL
 * 列能提供（见 CoderFactory::coderSupports）。那些场景下回退到 bwt_cm 而不是
 * coder_qual，代价从 7.4 个百分点降到 0.49 个百分点。
 *
 * 按记录逐条 encode_line，与 sam_actuator 里的实际调用方式保持一致。
 */
bool trialBwtCm(const std::vector<QualSampleRecord>& records,
                uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < records.size(); i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_bwt_cm coder(&io);
        for (size_t i = 0; i < records.size(); i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_line((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size());
        }
        coder.encode_flush();
        outLen = (io.data_len > 0) ? (uint32_t)io.data_len : 0;
    }
    usec = (uint32_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return outLen > 0 && outLen < buf.size();
}

} /* namespace */

FieldCodecSelection QualSelector::select(const std::vector<QualSampleRecord>& records,
                                         const std::vector<uint32_t>& freqByByte)
{
    FieldCodecSelection sel;

    uint32_t sampleLen = 0;
    for (size_t i = 0; i < records.size(); i++) {
        sampleLen += (uint32_t)records[i].qual.size();
    }
    sel.sampleLen = sampleLen;
    sel.decidedLen = sampleLen;
    sel.rounds = 1;

    if (records.empty() || sampleLen < MIN_QUAL_SAMPLE) {
        sel.status = FieldStatus::SKIPPED;
        return sel;
    }

    uint32_t qualLen = 0, qualUs = 0;
    uint32_t fcv2Len = 0, fcv2Us = 0;
    uint32_t cmLen = 0, cmUs = 0;
    bool qualOk = trialQual(records, freqByByte, qualLen, qualUs);
    bool fcv2Ok = trialFcv2(records, freqByByte, fcv2Len, fcv2Us);
    bool cmOk = trialBwtCm(records, cmLen, cmUs);

    sel.trialCount = 0;
    sel.addTrial(CoderType::QUAL, qualLen, qualUs);
    sel.addTrial(CoderType::FCV2, fcv2Len, fcv2Us);
    sel.addTrial(CoderType::BWT_CM, cmLen, cmUs);

    if (!qualOk && !fcv2Ok && !cmOk) {
        /* 候选全都跑不通，退回默认编码器，绝不因为评估失败而让压缩无法进行。 */
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    /*
     * 取压缩后最小的那个。三个候选的速度差异不足以改变结论：实测 fcv2 与
     * coder_qual 都在 30 MB/s 上下，bwt_cm 约 5.7 MB/s，而 bwt_cm 只在 fcv2
     * 跑不通时才可能被选中，那种场景下慢一些也比多付 7.4 个百分点划算。
     */
    sel.selectedCoder = CoderType::QUAL;
    sel.bestCompLen = 0;
    bool picked = false;
    for (uint32_t i = 0; i < sel.trialCount; ++i) {
        const bool ok = (sel.trialCoder[i] == CoderType::QUAL)   ? qualOk :
                        (sel.trialCoder[i] == CoderType::FCV2)   ? fcv2Ok : cmOk;
        if (!ok) {
            continue;
        }
        if (!picked || sel.trialLen[i] < sel.bestCompLen) {
            sel.selectedCoder = sel.trialCoder[i];
            sel.bestCompLen = sel.trialLen[i];
            picked = true;
        }
    }
    sel.status = FieldStatus::SELECTED;

    LOG_DEBUG("Qual codec trial: coder_qual=%u (%u us), fcv2=%u (%u us), bwt_cm=%u (%u us), picked=%s",
              qualLen, qualUs, fcv2Len, fcv2Us, cmLen, cmUs, coderTypeToMagic(sel.selectedCoder));
    return sel;
}
