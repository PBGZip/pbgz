/*
 * qual_selector.cpp - 质量值列的编码器评估
 *
 * 本文件刻意不包含 coder_fc.h 或任何会引入 fc/rangecoder.h 的头文件：那里的
 * RangeCoder 与 coder_qual.h 间接引入的 clr.h 中的同名类冲突。质量值的评估之所以
 * 从 codec_selector.cpp 里分出来，这是原因之一，详见 qual_selector.h 的说明。
 *
 * 与通用字段的 selectCoder 对应，质量值列也走"多轮收敛"（策略 7）：从 64KB 起步，
 * 每轮翻倍样本量，一旦领先者拉开 3% 以上的差距就定案，避免在小样本上被自适应编码器
 * 尚未收敛的假象带偏。此外还按数据特征给 fcv2 选上下文参数档位（策略 2，见
 * fqzcomp 的 fqz_pick_parameters 思路），把档位也作为候选参与试压，胜出的档位随
 * FieldCodecSelection.fcv2Params 交给压缩端与先验训练。
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

/* 领先者需要领先多少才算胜负已分，与通用字段 selectCoder 同一阈值。 */
const double SETTLE_MARGIN = 0.03;

/* 返回样本里前 n 条记录的累计 QUAL 字节数不超过 budget 的最大 n（只收完整记录）。 */
size_t recordsForBudget(const std::vector<QualSampleRecord>& records, size_t budget)
{
    size_t used = 0;
    size_t n = 0;
    for (; n < records.size(); n++) {
        if (used + records[n].qual.size() > budget) {
            break;
        }
        used += records[n].qual.size();
    }
    return n;
}

/* coder_bwt_cm 内部块大小，与 codec_selector.cpp 的 BWT_LEVEL_SIZE 保持一致。 */
const uint32_t kBwtLevelSize[10] = {
    0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
    1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
};

/* 取最小的、内部块能装下本轮样本的档位。 */
int bwtLevelFor(uint32_t sampleLen)
{
    for (int lv = 1; lv <= 9; ++lv) {
        if (kBwtLevelSize[lv] >= sampleLen) {
            return lv;
        }
    }
    return 9;
}

/*
 * 用 coder_qual 试压。
 *
 * 它以对应的碱基序列作为上下文，所以必须逐条记录喂入，同时给出 seq 和 qual。
 * 频率表的格式沿用 sam_actuator 的既有约定：按出现次数降序排列的字母表，
 * 第二个元素固定为 1（实际计数没有被保留下来）。
 */
bool trialQual(const std::vector<QualSampleRecord>& records, size_t recordCount,
               const std::vector<uint32_t>& freqByByte,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
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
        for (size_t i = 0; i < recordCount; i++) {
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
 * 所以这里传进去之后解码端不必再提供一次。cfg 决定上下文参数档位，见策略 2。
 */
bool trialFcv2(const std::vector<QualSampleRecord>& records, size_t recordCount,
               const std::vector<uint32_t>& freqByByte, const Fcv2Cfg& cfg,
               uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_fcv2 coder(&io, freqByByte, cfg);
        for (size_t i = 0; i < recordCount; i++) {
            const QualSampleRecord& r = records[i];
            if (r.qual.empty()) {
                continue;
            }
            coder.encode_record((const uint8_t*)r.qual.data(),
                                (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
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
 * 按记录逐条 encode_line，与 sam_actuator 里的实际调用方式保持一致。内部块大小按
 * 本轮样本量选档（与通用字段 selectCoder 的 pickBwtLevel 同理），小样本用不到大块。
 */
bool trialBwtCm(const std::vector<QualSampleRecord>& records, size_t recordCount,
                int bwtLevel, uint32_t& outLen, uint32_t& usec)
{
    size_t total = 0;
    for (size_t i = 0; i < recordCount; i++) {
        total += records[i].qual.size();
    }
    if (total == 0) {
        return false;
    }

    std::vector<uint8_t> buf(trialCapacity(total), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        io.set_level(bwtLevel);
        coder_bwt_cm coder(&io);
        for (size_t i = 0; i < recordCount; i++) {
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

/* 把 coder 层的档位翻译成 PreprocessInfo 携带的参数（字段一一对应）。 */
QualFcv2Params toQualParams(const Fcv2Cfg& cfg)
{
    QualFcv2Params p;
    p.cycleMax = cfg.cycleMax;
    p.cycleBucket = cfg.cycleBucket;
    p.deltaMax = cfg.deltaMax;
    p.deltaBucket = cfg.deltaBucket;
    p.prevShift = cfg.prevShift;
    p.useDelta = cfg.useDelta;
    p.useDedup = cfg.useDedup;
    p.useQa = cfg.useQa;
    return p;
}

/*
 * 按数据特征给 fcv2 选候选档位（策略 2）。
 *
 * 思路取自 fqzcomp 的 fqz_pick_parameters：质量值字母表小（NovaSeq/HiSeqX 类，符号数
 * 少）时上下文不稀疏，可以用更细的位置分档、不对前驱质量值做量化；字母表大或样本
 * 小时上下文稀疏，分档变粗、前驱右移更多。但统计规则只是先验，实测（合成 30MB 质量值
 * 上，细档比粗档好约 0.3 个百分点）显示哪个档位真正胜出与数据有关，所以把多档都放进
 * 候选，交由 select() 一并试压取压缩后更小的——参数选择以实测为准。
 *
 * 唯一的数据特征判断是样本量：细档需要足够样本收敛，小样本下放进候选只会平白增加
 * 试压时间，甚至被未收敛的假象误选，所以小样本时只保留默认档与粗档。
 */
std::vector<Fcv2Cfg> candidateFcv2Cfgs(size_t sampleLen)
{
    std::vector<Fcv2Cfg> cfgs;
    cfgs.push_back(Fcv2Cfg());   /* 默认档 */

    Fcv2Cfg dedup = Fcv2Cfg();   /* 默认档 + 相邻重复 read 去重（策略 3） */
    dedup.useDedup = true;
    cfgs.push_back(dedup);

    Fcv2Cfg qa = Fcv2Cfg();      /* 默认档 + read 平均质量档（策略 4） */
    qa.useQa = true;
    cfgs.push_back(qa);

    Fcv2Cfg fine;                /* 小字母表风格：细档、前驱不量化 */
    fine.cycleBucket = 24;
    fine.deltaBucket = 12;
    fine.prevShift = 0;

    Fcv2Cfg coarse;              /* 大字母表/小样本风格：粗档、前驱右移两位 */
    coarse.cycleBucket = 8;
    coarse.deltaBucket = 4;
    coarse.prevShift = 2;

    Fcv2Cfg coarseQa = coarse;   /* 粗档 + read 平均质量档（策略 4，实测最优组合） */
    coarseQa.useQa = true;
    cfgs.push_back(coarseQa);

    if (sampleLen >= (1u << 20)) {
        cfgs.push_back(fine);
    }
    cfgs.push_back(coarse);
    return cfgs;
}

/* 一次"按字节预算压缩全部候选"的回合。返回各候选的压缩大小与是否可用。 */
struct QualRoundResult {
    bool     qualOk = false;
    uint32_t qualLen = 0;
    bool     fcv2Ok = false;
    uint32_t fcv2Len = 0;
    Fcv2Cfg  fcv2Cfg;          /* 多个 fcv2 档位里胜出的那个 */
    bool     cmOk = false;
    uint32_t cmLen = 0;
    uint32_t bestLen = 0;
    CoderType bestCoder = CoderType::QUAL;
    bool     anyOk = false;
    uint32_t fcv2Us = 0;
    uint32_t qualUs = 0;
    uint32_t cmUs = 0;
};

QualRoundResult runRound(const std::vector<QualSampleRecord>& records,
                         const std::vector<uint32_t>& freqByByte,
                         const std::vector<Fcv2Cfg>& cfgs,
                         uint32_t probe)
{
    QualRoundResult r;
    size_t count = recordsForBudget(records, probe);

    if (trialQual(records, count, freqByByte, r.qualLen, r.qualUs)) {
        r.qualOk = true;
        r.anyOk = true;
        r.bestLen = r.qualLen;
        r.bestCoder = CoderType::QUAL;
    }

    /* 多个 fcv2 档位里取最小的那个，与 coder_qual / bwt_cm 比较时用这个代表 fcv2。 */
    bool fcv2Picked = false;
    for (size_t i = 0; i < cfgs.size(); i++) {
        uint32_t len = 0, us = 0;
        if (trialFcv2(records, count, freqByByte, cfgs[i], len, us)) {
            if (!fcv2Picked || len < r.fcv2Len) {
                r.fcv2Len = len;
                r.fcv2Us = us;
                r.fcv2Cfg = cfgs[i];
                fcv2Picked = true;
            }
        }
    }
    if (fcv2Picked) {
        r.fcv2Ok = true;
        r.anyOk = true;
        if (r.fcv2Len < r.bestLen || !r.qualOk) {
            r.bestLen = r.fcv2Len;
            r.bestCoder = CoderType::FCV2;
        }
    }

    if (trialBwtCm(records, count, bwtLevelFor(probe), r.cmLen, r.cmUs)) {
        r.cmOk = true;
        r.anyOk = true;
        if (r.cmLen < r.bestLen || (!r.qualOk && !r.fcv2Ok)) {
            r.bestLen = r.cmLen;
            r.bestCoder = CoderType::BWT_CM;
        }
    }
    return r;
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

    const std::vector<Fcv2Cfg> cfgs = candidateFcv2Cfgs(sampleLen);

    /*
     * 多轮收敛（策略 7）：从小样本起步，每轮翻倍，一旦领先者拉开足够差距就收手。
     * 每轮都是拿新样本量重压全部候选，而不是增量喂——重压恰好模拟了实际压缩时
     * "一个数据块从头压到尾"的情形，量出来的数字更贴近真实表现。
     * 与通用字段 selectCoder 的理由一致，见 codec_selector.cpp 的 selectCoder。
     *
     * 提前定案有一个前提：样本必须大到自适应编码器已经收敛。实测本文件的质量值上，
     * bwt_cm 在小样本（64KB~512KB）时领先 fcv2 达 2%~5%（它在这个区间压得异常好），
     * 但随样本增大优势逐步收窄，约 2~4MB 处被 fcv2 反超。若一上来就按 3% 领先阈值
     * 定案，会在第一轮就误选 bwt_cm。因此设一个最低定案样本量，低于它只翻倍不定案。
     */
    const uint32_t MIN_SETTLE_PROBE = 1u << 20;   /* 1 MB */

    uint32_t probe = (MIN_QUAL_SAMPLE < sampleLen) ? MIN_QUAL_SAMPLE : sampleLen;
    QualRoundResult final;
    bool finalSet = false;

    while (true) {
        QualRoundResult r = runRound(records, freqByByte, cfgs, probe);
        final = r;
        finalSet = true;
        sel.rounds++;

        if (!r.anyOk) {
            break;
        }
        if (probe >= sampleLen) {
            break;
        }
        if (probe < MIN_SETTLE_PROBE) {
            probe = (probe > sampleLen / 2) ? sampleLen : (probe * 2);
            continue;
        }
        /* 只有一个候选能跑通，再加数据也没有比较对象。 */
        {
            uint32_t runnerUp = UINT32_MAX;
            if (r.bestCoder == CoderType::QUAL) {
                if (r.fcv2Ok) runnerUp = r.fcv2Len;
                if (r.cmOk && r.cmLen < runnerUp) runnerUp = r.cmLen;
            } else if (r.bestCoder == CoderType::FCV2) {
                runnerUp = r.qualOk ? r.qualLen : UINT32_MAX;
                if (r.cmOk && r.cmLen < runnerUp) runnerUp = r.cmLen;
            } else { /* BWT_CM */
                runnerUp = r.qualOk ? r.qualLen : UINT32_MAX;
                if (r.fcv2Ok && r.fcv2Len < runnerUp) runnerUp = r.fcv2Len;
            }
            if (runnerUp == UINT32_MAX) {
                break;
            }
            /* 领先者已经拉开足够差距，再加数据不会反转。 */
            if ((double)(runnerUp - r.bestLen) / (double)runnerUp >= SETTLE_MARGIN) {
                break;
            }
        }
        probe = (probe > sampleLen / 2) ? sampleLen : (probe * 2);
    }

    sel.decidedLen = finalSet ? probe : sampleLen;
    sel.trialCount = 0;
    sel.addTrial(CoderType::QUAL, final.qualOk ? final.qualLen : 0, final.qualUs);
    sel.addTrial(CoderType::FCV2, final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us);
    sel.addTrial(CoderType::BWT_CM, final.cmOk ? final.cmLen : 0, final.cmUs);

    if (!final.anyOk) {
        sel.status = FieldStatus::FAILED;
        return sel;
    }

    sel.selectedCoder = final.bestCoder;
    sel.bestCompLen = final.bestLen;
    if (final.bestCoder == CoderType::FCV2) {
        /* 胜出的是 fcv2，把最后一轮里胜出的档位交给上层，编码与先验训练据此保持一致。 */
        sel.fcv2Params = toQualParams(final.fcv2Cfg);
    }
    sel.status = FieldStatus::SELECTED;

    LOG_DEBUG("Qual codec trial: coder_qual=%u (%u us), fcv2=%u (%u us), bwt_cm=%u (%u us), picked=%s",
              final.qualOk ? final.qualLen : 0, final.qualUs,
              final.fcv2Ok ? final.fcv2Len : 0, final.fcv2Us,
              final.cmOk ? final.cmLen : 0, final.cmUs,
              coderTypeToMagic(sel.selectedCoder));
    return sel;
}
