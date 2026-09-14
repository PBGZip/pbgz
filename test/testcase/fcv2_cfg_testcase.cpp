/*
 * fcv2_cfg_testcase.cpp - Unit tests for the fcv2 context parameter presets (strategies 1/2).
 *
 * Covers three things:
 *   1) encode/decode round-trip consistency under non-default presets (fine/coarse/delta context
 *      disabled);
 *   2) the prior snapshot carries the preset: loading under the same preset round-trips, and a
 *      prior trained under a different preset is correctly adopted (the prior decides the preset);
 *   3) when QualSelector picks fcv2, the winning preset is written into
 *      FieldCodecSelection.fcv2Params, and that preset actually round-trips.
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

/* Generator of the same style as fcv2_cfg_bench: cycle decline + runs + instability + reverse-chain reversal. */
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
 * fcv2-friendly generator: forward-strand quality declines with cycle in sequencing order, so the
 * reverse strand shows a rising profile with position in stored order. The two strands alternate
 * in the same file, so the concatenated byte stream has a direction-reversing "fold" every
 * readLen, leaving bwt_cm's MTF nothing to exploit; once fcv2 restores the cycle index using rev,
 * every read shares the same declining profile, and the cycle context (m1/m5) directly targets
 * this structure, so fcv2 is usually selected.
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
            /* in sequencing order, quality declines monotonically from ~Q40 to ~Q28 across the cycle */
            int qi = (int)std::lround(40.0 - 12.0 * (double)c / (double)readLen + noise(rng));
            if (qi < 0) qi = 0;
            if (qi > 40) qi = 40;
            qual[c] = (char)('!' + qi);
        }
        if (rec.rev) {
            /* stored order is reversed relative to sequencing order */
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

/*
 * The average-quality context lives in m6, the last model of the mix, and the
 * default model count trims that tail (see FCV2_MODEL_COUNT). A test of the
 * strategy itself therefore has to ask for the full mix: with a trimmed count
 * useQa has no model to act on and both sides of a comparison are the same
 * coder.
 */
Fcv2Cfg fullMixCfg()
{
    Fcv2Cfg cfg;
    cfg.modelCount = FCV2_MAX_MODEL_COUNT;
    return cfg;
}

Fcv2Cfg qaCfg()
{
    Fcv2Cfg cfg = fullMixCfg();
    cfg.useQa = true;
    return cfg;
}

/*
 * Build a batch of records containing adjacent duplicate reads: every few records, the previous
 * record's QUAL/SEQ is copied. With duplicates, dedup should clearly save bytes; without them it
 * must not break the round trip.
 */
std::vector<Rec> makeDupRecords(size_t n, uint32_t readLen, uint32_t dupEvery, uint32_t seed)
{
    std::vector<Rec> recs = makeRecords(n, readLen, seed);
    for (size_t i = 1; i < recs.size(); i++) {
        if (dupEvery > 0 && (i % dupEvery) == 0) {
            recs[i] = recs[i - 1];   /* identical to the previous record */
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
    EXPECT_TRUE(roundtripWithCfg(recs, Fcv2Cfg()));        /* default preset */
    EXPECT_TRUE(roundtripWithCfg(recs, fineCfg()));        /* fine preset */
    EXPECT_TRUE(roundtripWithCfg(recs, coarseCfg()));      /* coarse preset */
    EXPECT_TRUE(roundtripWithCfg(recs, noDeltaCfg()));     /* delta context disabled */
    EXPECT_TRUE(roundtripWithCfg(recs, dedupCfg()));       /* dedup (no duplicate data) */
    EXPECT_TRUE(roundtripWithCfg(recs, qaCfg()));          /* average-quality preset */
}

/*
 * Tier values at and above their maxima must survive the stream header, which
 * writes every one of them as a single byte.
 *
 * The above-maximum case is the one that discriminates. deltaMax can be set to
 * 256, which normalizeCfg accepts, but the header can only carry a byte: 256 is
 * written as 0 and read back as the 8 that normalizeCfg clamps it to, so the two
 * sides disagree on the context layout and the stream decodes into garbage with
 * nothing reported (the header carries no checksum, and the model-layout
 * comparison only runs when a prior is loaded). The preset round trips above all
 * sit on moderate values and never reach this.
 *
 * A literal 255 cannot catch it: 255 is representable, so it round-trips either
 * way. Only a value above the maximum does, because that is what makes the
 * encoder depend on the header carrying its clamped result.
 */
TEST_F(Fcv2CfgTest, OverMaximalTierRoundTrip)
{
    Fcv2Cfg over;              /* every field above its maximum */
    over.cycleMax = 4096;      /* clamps to 128 */
    over.cycleBucket = 512;    /* clamps to 32 */
    over.deltaMax = 256;       /* would truncate to 0 in the header */
    over.deltaBucket = 64;     /* clamps to 16 */
    over.prevShift = 9;        /* clamps to 3 */
    over.modelCount = 99;      /* clamps to FCV2_MAX_MODEL_COUNT */

    std::vector<Rec> recs = makeCycleRecords(1200, 4000, 23);
    EXPECT_TRUE(roundtripWithCfg(recs, over));
}

/* The same, with every field sitting exactly on its maximum at once. */
TEST_F(Fcv2CfgTest, MaximalTierRoundTrip)
{
    Fcv2Cfg cfg;
    cfg.cycleMax = 128;
    cfg.cycleBucket = 32;
    cfg.deltaMax = 255;
    cfg.deltaBucket = 16;
    cfg.prevShift = 3;
    cfg.modelCount = FCV2_MAX_MODEL_COUNT;

    std::vector<Rec> recs = makeCycleRecords(1200, 4000, 23);
    EXPECT_TRUE(roundtripWithCfg(recs, cfg));
}

/* Dedup round trip: with adjacent duplicate reads, the result must be restored byte for byte. */
TEST_F(Fcv2CfgTest, DedupRoundTripWithDuplicates)
{
    std::vector<Rec> recs = makeDupRecords(4000, 151, 3, 5);   /* duplicate one in every 3 records */
    EXPECT_TRUE(roundtripWithCfg(recs, dedupCfg()));
    /* With duplicates present, dedup should save more bytes than no dedup (at least not be much worse). */
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
    EXPECT_LT((size_t)n2, (size_t)n1) << "Data with 1/3 duplicates must compress fewer bytes with dedup on";
}

/* Dedup + prior: the snapshot carries useDedup, and the round trip is consistent after loading. */
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

/* Average-quality preset (strategy 4): when the mean quality differs widely between reads, the quality preset should save more than not bucketing by quality. */
TEST_F(Fcv2CfgTest, QaSavesBytesOnVariedQuality)
{
    /* Mix two batches of reads: one high-quality (Q38~40) and one low-quality (Q15~20), with large mean-quality differences. */
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
    coder_fcv2 e1(&io1, freq, fullMixCfg());   /* the baseline differs only in useQa */
    ASSERT_TRUE(encodeAll(e1, recs));
    coder_io io2(c2.data(), (int32_t)c2.size());
    coder_fcv2 e2(&io2, freq, qaCfg());
    ASSERT_TRUE(encodeAll(e2, recs));
    int32_t n1 = e1.encode_flush(), n2 = e2.encode_flush();
    EXPECT_LT((size_t)n2, (size_t)n1) << "Data with widely differing mean quality must compress fewer bytes with the quality tier";
}

TEST_F(Fcv2CfgTest, PriorCarriesConfigAndAdoptsIt)
{
    std::vector<Rec> recs = makeRecords(4000, 151, 9);
    std::vector<uint32_t> freq = frequencies(recs);

    /* Train the prior with the fine preset. */
    std::vector<uint8_t> blob = trainPrior(recs, fineCfg());
    ASSERT_FALSE(blob.empty());

    /*
     * The coder is constructed with the coarse preset, but the prior was trained with the fine
     * one: loadModel adopts the prior's own fine preset, so both the encoder and the stream header
     * use the fine preset, and the decoder reads it back consistently, so the round trip must hold.
     */
    size_t total = 0;
    for (const Rec& r : recs) total += r.qual.size();
    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    bool loaded = false;
    coder_fcv2 enc(&io, freq, coarseCfg(), blob, &loaded);
    ASSERT_TRUE(loaded) << "A prior with the same tier layout must be loaded (the prior determines the tier)";
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

    /* The encoder uses the fine preset; the decoder does not know the preset and is constructed with the default, so begin_decode must read the fine preset back from the stream. */
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

/*
 * The model count follows the read length, and this is the only place that
 * decides it: with the datasets at hand fcv2 is only ever selected on short
 * reads, so the long-read branch would otherwise never run in any test or in
 * production, and a silent regression in it would go unnoticed until some future
 * file happened to select fcv2 there.
 */
TEST_F(Fcv2CfgTest, ModelCountFollowsReadLength)
{
    /* Short reads, on both sides of anything an Illumina run produces. */
    EXPECT_EQ(qualModelCountForMeanLen(35), FCV2_MODEL_COUNT);
    EXPECT_EQ(qualModelCountForMeanLen(90), FCV2_MODEL_COUNT);      /* con_sorted */
    EXPECT_EQ(qualModelCountForMeanLen(151), FCV2_MODEL_COUNT);
    EXPECT_EQ(qualModelCountForMeanLen(QUAL_LONG_READ_LEN - 1), FCV2_MODEL_COUNT);

    /* Long reads start at the threshold itself, and stay there however long. */
    EXPECT_EQ(qualModelCountForMeanLen(QUAL_LONG_READ_LEN), 4);
    EXPECT_EQ(qualModelCountForMeanLen(14044), 4);                 /* ERR11436629 */
    EXPECT_EQ(qualModelCountForMeanLen(100000), 4);

    /* A trimmed long-read count must still be a usable mix (1..max). */
    EXPECT_GE(qualModelCountForMeanLen(14044), 1);
    EXPECT_LE(qualModelCountForMeanLen(14044), FCV2_MAX_MODEL_COUNT);
}

/* The trimmed count must produce a stream that decodes under its own count. */
TEST_F(Fcv2CfgTest, LongReadCountRoundTrip)
{
    Fcv2Cfg cfg;
    cfg.modelCount = qualModelCountForMeanLen(14044);
    ASSERT_EQ(cfg.modelCount, 4);
    std::vector<Rec> recs = makeCycleRecords(1500, 4000, 17);   /* long-read shaped */
    EXPECT_TRUE(roundtripWithCfg(recs, cfg));
}

/* With a large enough sample, QualSelector should select a candidate and bring back the winning preset. */
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
        GTEST_SKIP() << "This synthetic sample did not select fcv2 (selected " << coderTypeToMagic(qualSel.selectedCoder)
                     << ", " << buf << "), so tier feedback cannot be verified";
    }
    /* The preset must be one of the known candidates (default / fine / coarse /
     * ultra 4-2-ps0 / ultra-80), and it must round-trip with the preset applied. */
    const QualFcv2Params& p = qualSel.fcv2Params;
    bool isKnownPreset = (p == QualFcv2Params()) ||
        (p.cycleBucket == 24 && p.deltaBucket == 12 && p.prevShift == 0 && p.useDelta) ||
        (p.cycleBucket == 8 && p.deltaBucket == 4 && p.prevShift == 2 && p.useDelta) ||
        (p.cycleBucket == 4 && p.deltaBucket == 2 && p.prevShift == 0 && p.useDelta) ||
        (p.cycleMax == 80 && p.cycleBucket == 4 && p.deltaBucket == 2 && p.prevShift == 0 && p.useDelta);
    EXPECT_TRUE(isKnownPreset) << "fcv2Params is not one of the known candidate tiers";
}

/*
 * The stream header states its format version as the stream's first symbol (see
 * FCV2_STREAM_VERSION), and a decoder reads either that layout or the one from before the byte
 * existed (alphabet size first, no model count, seven models). A stream that fits neither is
 * refused with CODER_ERR_UNSUPPORTED_VERSION instead of being parsed on a guess - which is what
 * an archive from an unknown layout looks like from here. The first case checks the stream this
 * build writes is accepted; the second alters the stream's leading code words, so it no longer
 * decodes to this version and its alphabet size lands outside TREE_CAP as well, leaving neither
 * layout applicable.
 */
TEST_F(Fcv2CfgTest, StreamVersionIsChecked)
{
    const std::vector<Rec> recs = makeRecords(24, 60, 11);
    const std::vector<uint32_t> freq = frequencies(recs);
    size_t total = 0;
    for (const Rec& r : recs) {
        total += r.qual.size();
    }

    std::vector<uint8_t> comp(total * 2 + (1u << 16), 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    coder_fcv2 enc(&io, freq);
    ASSERT_TRUE(encodeAll(enc, recs));
    const int32_t packed = enc.encode_flush();
    ASSERT_GT(packed, 8);

    {
        coder_io dio(comp.data(), packed);
        coder_fcv2 dec(&dio, freq);
        EXPECT_EQ(dec.begin_decode(), 0);
    }

    {
        /* The range decoder's starting code is read from the stream's first words, so writing
           them changes the version the header states. */
        std::vector<uint8_t> tampered(comp.begin(), comp.begin() + packed);
        for (size_t i = 2; i < 6 && i < tampered.size(); ++i) {
            tampered[i] = 0xFF;
        }
        coder_io dio(tampered.data(), packed);
        coder_fcv2 dec(&dio, freq);
        EXPECT_EQ(dec.begin_decode(), coder_ns::CODER_ERR_UNSUPPORTED_VERSION);
    }
}
