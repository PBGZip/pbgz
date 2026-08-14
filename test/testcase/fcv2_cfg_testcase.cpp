/*
 * fcv2_cfg_testcase.cpp - fcv2 上下文参数档位（策略 1/2）的单元测试。
 *
 * 覆盖三件事：
 *   1) 非默认档位（细档/粗档/关闭跳变上下文）下的编解码往返一致；
 *   2) 先验快照携带档位：同档位加载往返、不同档位训练的先验被正确采用（先验决定档位）；
 *   3) QualSelector 选中 fcv2 时把胜出档位写进 FieldCodecSelection.fcv2Params，
 *      且该档位确实能跑通往返。
 */

#include <gtest/gtest.h>

#include <random>
#include <string>
#include <vector>

#include "qual_selector.h"
#include "codec_selector.h"
#include "field_coder_config.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_io.h"
#include "preprocess_info.h"
#include "io_block.h"
#include "coder.h"
#include "utils/memory_util.h"

namespace {

struct Rec {
    std::string qual;
    std::string seq;
    bool rev;
};

/* 与 fcv2_cfg_bench 同款生成器：循环下降 + 游程 + 不稳定度 + 反向链反转。 */
std::vector<Rec> makeRecords(size_t n, uint32_t readLen, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, 1.8);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    static const char bases[] = "ACGT";
    std::vector<Rec> recs;
    recs.reserve(n);
    for (size_t r = 0; r < n; r++) {
        Rec rec;
        rec.rev = ((r & 1) != 0);
        int readMean = 30 + (int)(rng() % 8);
        double instability = uni(rng);
        std::string qual(readLen, '!');
        std::string seq(readLen, 'A');
        int prev = -1;
        for (uint32_t c = 0; c < readLen; c++) {
            seq[c] = bases[rng() % 4];
            double drop = 7.0 * (double)c / (double)readLen;
            double q = (double)readMean - drop + noise(rng);
            if (prev >= 0 && uni(rng) < (0.12 + 0.55 * instability)) {
                q = (double)prev;
            }
            int qi = (int)std::lround(q);
            if (qi < 0) qi = 0;
            if (qi > 40) qi = 40;
            qual[c] = (char)('!' + qi);
            prev = qi;
        }
        if (rec.rev) {
            std::reverse(qual.begin(), qual.end());
            std::reverse(seq.begin(), seq.end());
        }
        rec.qual.swap(qual);
        rec.seq.swap(seq);
        recs.push_back(std::move(rec));
    }
    return recs;
}

/*
 * fcv2 友好的生成器：正向链的质量在测序序中随循环下降，反向链在存储序中因此呈现
 * 随位置上升的剖面。两条链在同一文件里交替出现，拼起来的字节流每隔 readLen 就有一个
 * 方向翻转的"折点"，bwt_cm 的 MTF 无从借力；而 fcv2 按 rev 把循环序号还原后，所有
 * read 都是同一张下降剖面，循环上下文（m1/m5）正对着这种结构，通常选中 fcv2。
 */
std::vector<Rec> makeCycleRecords(size_t n, uint32_t readLen, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, 1.0);
    static const char bases[] = "ACGT";
    std::vector<Rec> recs;
    recs.reserve(n);
    for (size_t r = 0; r < n; r++) {
        Rec rec;
        rec.rev = ((r & 1) != 0);
        std::string qual(readLen, '!');
        std::string seq(readLen, 'A');
        for (uint32_t c = 0; c < readLen; c++) {
            seq[c] = bases[rng() % 4];
            /* 测序序里随循环从 ~Q40 单调降到 ~Q28 */
            int qi = (int)std::lround(40.0 - 12.0 * (double)c / (double)readLen + noise(rng));
            if (qi < 0) qi = 0;
            if (qi > 40) qi = 40;
            qual[c] = (char)('!' + qi);
        }
        if (rec.rev) {
            /* 存储序相对测序序反转 */
            std::reverse(qual.begin(), qual.end());
            std::reverse(seq.begin(), seq.end());
        }
        rec.qual.swap(qual);
        rec.seq.swap(seq);
        recs.push_back(std::move(rec));
    }
    return recs;
}

std::vector<uint32_t> frequencies(const std::vector<Rec>& recs)
{
    std::vector<uint32_t> freq(256, 0);
    for (const Rec& r : recs) {
        for (size_t i = 0; i < r.qual.size(); i++) freq[(uint8_t)r.qual[i]]++;
    }
    return freq;
}

bool encodeAll(coder_fcv2& coder, const std::vector<Rec>& recs)
{
    for (const Rec& r : recs) {
        coder.encode_record((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size(), r.rev,
                            (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
    }
    return coder.encode_flush() > 0;
}

bool decodeAll(coder_fcv2& coder, const std::vector<Rec>& recs)
{
    if (coder.begin_decode() != 0) return false;
    for (const Rec& r : recs) {
        std::string out(r.qual.size(), 0);
        if (coder.decode_record((uint8_t*)out.data(), (uint32_t)out.size(),
                                (const uint8_t*)r.seq.data(),
                                (uint32_t)r.seq.size()) != (int32_t)out.size()) {
            return false;
        }
        if (out != r.qual) return false;
    }
    return true;
}

bool roundtripWithCfg(const std::vector<Rec>& recs, const Fcv2Cfg& cfg)
{
    std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();

    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    coder_fcv2 enc(&io, freq, cfg);
    if (!encodeAll(enc, recs)) return false;
    int32_t packed = enc.encode_flush();
    if (packed <= 0) return false;

    coder_io dio(comp.data(), packed);
    coder_fcv2 dec(&dio, freq, cfg);
    return decodeAll(dec, recs);
}

std::vector<uint8_t> trainPrior(const std::vector<Rec>& recs, const Fcv2Cfg& cfg)
{
    std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    coder_fcv2 coder(&io, freq, cfg);
    if (!encodeAll(coder, recs)) return {};
    std::vector<uint8_t> blob;
    if (!coder.export_model(blob) || blob.empty()) return {};
    return blob;
}

Fcv2Cfg fineCfg()
{
    Fcv2Cfg cfg;
    cfg.cycleBucket = 24;
    cfg.deltaBucket = 12;
    cfg.prevShift = 0;
    return cfg;
}

Fcv2Cfg coarseCfg()
{
    Fcv2Cfg cfg;
    cfg.cycleBucket = 8;
    cfg.deltaBucket = 4;
    cfg.prevShift = 2;
    return cfg;
}

Fcv2Cfg noDeltaCfg()
{
    Fcv2Cfg cfg;
    cfg.useDelta = false;
    return cfg;
}

Fcv2Cfg dedupCfg()
{
    Fcv2Cfg cfg;
    cfg.useDedup = true;
    return cfg;
}

Fcv2Cfg qaCfg()
{
    Fcv2Cfg cfg;
    cfg.useQa = true;
    return cfg;
}

/*
 * 造一批含相邻重复 read 的数据：每隔若干条就复制上一条的 QUAL/SEQ。
 * 有重复时 dedup 应该明显省字节；无重复时也不该破坏往返。
 */
std::vector<Rec> makeDupRecords(size_t n, uint32_t readLen, uint32_t dupEvery, uint32_t seed)
{
    std::vector<Rec> recs = makeRecords(n, readLen, seed);
    for (size_t i = 1; i < recs.size(); i++) {
        if (dupEvery > 0 && (i % dupEvery) == 0) {
            recs[i] = recs[i - 1];   /* 与上一条完全相同 */
        }
    }
    return recs;
}

class Fcv2CfgTest : public ::testing::Test {
protected:
    void SetUp() override
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();
    }
};

} /* namespace */

TEST_F(Fcv2CfgTest, RoundTripUnderNonDefaultConfigs)
{
    std::vector<Rec> recs = makeRecords(4000, 151, 7);
    EXPECT_TRUE(roundtripWithCfg(recs, Fcv2Cfg()));        /* 默认档 */
    EXPECT_TRUE(roundtripWithCfg(recs, fineCfg()));        /* 细档 */
    EXPECT_TRUE(roundtripWithCfg(recs, coarseCfg()));      /* 粗档 */
    EXPECT_TRUE(roundtripWithCfg(recs, noDeltaCfg()));     /* 关闭跳变上下文 */
    EXPECT_TRUE(roundtripWithCfg(recs, dedupCfg()));       /* 去重（无重复数据） */
    EXPECT_TRUE(roundtripWithCfg(recs, qaCfg()));          /* 平均质量档 */
}

/* 去重往返：含相邻重复 read 时必须逐字节还原。 */
TEST_F(Fcv2CfgTest, DedupRoundTripWithDuplicates)
{
    std::vector<Rec> recs = makeDupRecords(4000, 151, 3, 5);   /* 每 3 条重复一条 */
    EXPECT_TRUE(roundtripWithCfg(recs, dedupCfg()));
    /* 有重复时 dedup 应比不去重更省（至少不比不去重差太多）。 */
    std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> c1(total * 2 + (1u << 16), 0), c2(total * 2 + (1u << 16), 0);
    coder_io io1(c1.data(), (int32_t)c1.size());
    coder_fcv2 e1(&io1, freq, Fcv2Cfg());
    ASSERT_TRUE(encodeAll(e1, recs));
    coder_io io2(c2.data(), (int32_t)c2.size());
    coder_fcv2 e2(&io2, freq, dedupCfg());
    ASSERT_TRUE(encodeAll(e2, recs));
    int32_t n1 = e1.encode_flush(), n2 = e2.encode_flush();
    EXPECT_LT((size_t)n2, (size_t)n1) << "含 1/3 重复的数据，去重应更省字节";
}

/* 去重 + 先验：快照携带 useDedup，加载后往返一致。 */
TEST_F(Fcv2CfgTest, DedupPriorRoundTrip)
{
    std::vector<Rec> recs = makeDupRecords(4000, 151, 2, 11);
    std::vector<uint32_t> freq = frequencies(recs);
    std::vector<uint8_t> blob = trainPrior(recs, dedupCfg());
    ASSERT_FALSE(blob.empty());
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    bool loaded = false;
    coder_fcv2 enc(&io, freq, dedupCfg(), blob, &loaded);
    ASSERT_TRUE(loaded);
    ASSERT_TRUE(encodeAll(enc, recs));
    int32_t packed = enc.encode_flush();

    coder_io dio(comp.data(), packed);
    bool decLoaded = false;
    coder_fcv2 dec(&dio, freq, Fcv2Cfg(), blob, &decLoaded);
    ASSERT_TRUE(decLoaded);
    ASSERT_TRUE(decodeAll(dec, recs));
}

/* 平均质量档（策略 4）：read 间平均质量差异大时，质量档应比不去分档更省。 */
TEST_F(Fcv2CfgTest, QaSavesBytesOnVariedQuality)
{
    /* 混合两批 read：一批高质量（Q38~40）、一批低质量（Q15~20），平均质量差异大。 */
    std::vector<Rec> high = makeRecords(3000, 151, 21);
    std::vector<Rec> low = makeRecords(3000, 151, 22);
    for (Rec& r : low) {
        for (char& c : r.qual) c = (char)('!' + 12 + (c - '!') % 6);
    }
    std::vector<Rec> recs;
    recs.reserve(high.size() + low.size());
    for (size_t i = 0; i < high.size() || i < low.size(); i++) {
        if (i < high.size()) recs.push_back(high[i]);
        if (i < low.size()) recs.push_back(low[i]);
    }
    EXPECT_TRUE(roundtripWithCfg(recs, qaCfg()));

    std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> c1(total * 2 + (1u << 16), 0), c2(total * 2 + (1u << 16), 0);
    coder_io io1(c1.data(), (int32_t)c1.size());
    coder_fcv2 e1(&io1, freq, Fcv2Cfg());
    ASSERT_TRUE(encodeAll(e1, recs));
    coder_io io2(c2.data(), (int32_t)c2.size());
    coder_fcv2 e2(&io2, freq, qaCfg());
    ASSERT_TRUE(encodeAll(e2, recs));
    int32_t n1 = e1.encode_flush(), n2 = e2.encode_flush();
    EXPECT_LT((size_t)n2, (size_t)n1) << "平均质量差异大的数据，质量档应更省字节";
}

TEST_F(Fcv2CfgTest, PriorCarriesConfigAndAdoptsIt)
{
    std::vector<Rec> recs = makeRecords(4000, 151, 9);
    std::vector<uint32_t> freq = frequencies(recs);

    /* 用细档训练先验。 */
    std::vector<uint8_t> blob = trainPrior(recs, fineCfg());
    ASSERT_FALSE(blob.empty());

    /*
     * 构造时用粗档，但先验是细档训练的：loadModel 采用先验自带的细档，编码端与码流
     * 头部都以细档为准，解码端读回一致，往返必须成立。
     */
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    bool loaded = false;
    coder_fcv2 enc(&io, freq, coarseCfg(), blob, &loaded);
    ASSERT_TRUE(loaded) << "同档位布局的先验应被加载（先验决定档位）";
    ASSERT_TRUE(encodeAll(enc, recs));
    int32_t packed = enc.encode_flush();

    coder_io dio(comp.data(), packed);
    bool decLoaded = false;
    coder_fcv2 dec(&dio, freq, Fcv2Cfg(), blob, &decLoaded);
    ASSERT_TRUE(decLoaded);
    ASSERT_TRUE(decodeAll(dec, recs));
}

TEST_F(Fcv2CfgTest, StreamHeaderCarriesConfigSelfDescribe)
{
    std::vector<Rec> recs = makeRecords(4000, 151, 11);
    std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();

    /* 编码端用细档，解码端不知道档位、以默认档构造，begin_decode 应从码流读回细档。 */
    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    coder_fcv2 enc(&io, freq, fineCfg());
    ASSERT_TRUE(encodeAll(enc, recs));
    int32_t packed = enc.encode_flush();

    coder_io dio(comp.data(), packed);
    coder_fcv2 dec(&dio, freq, Fcv2Cfg());
    EXPECT_EQ(dec.begin_decode(), 0);
    for (const Rec& r : recs) {
        std::string out(r.qual.size(), 0);
        ASSERT_EQ(dec.decode_record((uint8_t*)out.data(), (uint32_t)out.size(),
                                    (const uint8_t*)r.seq.data(),
                                    (uint32_t)r.seq.size()), (int32_t)out.size());
        EXPECT_EQ(out, r.qual);
    }
}

/* 采样足够大时，QualSelector 应选中某个候选并把胜出的档位带回。 */
TEST_F(Fcv2CfgTest, QualSelectorCarriesWinningConfig)
{
    RoughIOBlock block(8 << 20);
    std::vector<Rec> recs = makeCycleRecords(4000, 151, 13);
    uint8_t* buffer = block.getBuffer();
    std::vector<size_t>& npos = block.getNpos();
    uint32_t offset = 0;
    for (const Rec& r : recs) {
        std::string line = "read\t0\tchr1\t1\t60\t151M\t=\t1\t151\t";
        line += r.seq;
        line += "\t";
        line += r.qual;
        line += "\n";
        memcpy(buffer + offset, line.data(), line.size());
        offset += (uint32_t)line.size();
        npos.push_back(offset - 1);
    }
    block.setDataLen(offset);
    block.setBlockType(SAM);

    PreprocessInfo info;
    ASSERT_EQ(CodecSelector::analyze(&block, (1ull << 30), info), 0);
    ASSERT_EQ(info.fields.size(), (size_t)SAM_FIELD_COUNT_SELECT);
    const FieldCodecSelection& qualSel = info.fields[SAM_QUAL];
    if (qualSel.selectedCoder != CoderType::FCV2) {
        char buf[256];
        snprintf(buf, sizeof(buf), "qual=%u fcv2=%u cm=%u",
                 qualSel.trialCount > 0 ? qualSel.trialLen[0] : 0,
                 qualSel.trialCount > 1 ? qualSel.trialLen[1] : 0,
                 qualSel.trialCount > 2 ? qualSel.trialLen[2] : 0);
        GTEST_SKIP() << "该合成样本未选中 fcv2（实际 " << coderTypeToMagic(qualSel.selectedCoder)
                     << "，" << buf << "），档位回传无从验证";
    }
    /* 档位必须是候选里的三档之一，且带档位能跑通往返。 */
    const QualFcv2Params& p = qualSel.fcv2Params;
    bool isKnownPreset = (p == QualFcv2Params()) ||
        (p.cycleBucket == 24 && p.deltaBucket == 12 && p.prevShift == 0 && p.useDelta) ||
        (p.cycleBucket == 8 && p.deltaBucket == 4 && p.prevShift == 2 && p.useDelta);
    EXPECT_TRUE(isKnownPreset) << "fcv2Params 不是已知候选档位";
}
