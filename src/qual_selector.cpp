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
    bool qualOk = trialQual(records, freqByByte, qualLen, qualUs);
    bool fcv2Ok = trialFcv2(records, freqByByte, fcv2Len, fcv2Us);

    /*
     * 复用既有的两个试压结果字段来携带这两个候选的数字：这里的 trialBwtCmLen 存
     * coder_qual 的结果，trialFcLen 存 fcv2 的结果。字段名沿用通用路径的命名，
     * 是为了不改动 -v 的打印结构；质量值这一行的含义以本函数为准。
     */
    sel.trialBwtCmLen = qualLen;
    sel.trialBwtCmUs = qualUs;
    sel.trialFcLen = fcv2Len;
    sel.trialFcUs = fcv2Us;

    if (!qualOk && !fcv2Ok) {
        /* 两个候选都跑不通，退回默认编码器，绝不因为评估失败而让压缩无法进行。 */
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    if (fcv2Ok && (!qualOk || fcv2Len < qualLen)) {
        sel.selectedCoder = CoderType::FCV2;
        sel.bestCompLen = fcv2Len;
    } else {
        sel.selectedCoder = CoderType::QUAL;
        sel.bestCompLen = qualLen;
    }
    sel.status = FieldStatus::SELECTED;

    LOG_DEBUG("Qual codec trial: coder_qual=%u (%u us), fcv2=%u (%u us), picked=%s",
              qualLen, qualUs, fcv2Len, fcv2Us, coderTypeToMagic(sel.selectedCoder));
    return sel;
}
