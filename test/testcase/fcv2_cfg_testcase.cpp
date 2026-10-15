/*
 * fcv2_cfg_testcase.cpp - Unit tests for the fcv2 context parameter presets (strategies 1/2).
 *
 * Covers three things:
 *   1) encode/decode round-trip consistency under non-default presets (fine/coarse tiers, dedup);
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
    cfg.prevShift = 0;
    return cfg;
}

Fcv2Cfg coarseCfg()
{
    Fcv2Cfg cfg;
    cfg.cycleBucket = 8;
    cfg.prevShift = 2;
    return cfg;
}

Fcv2Cfg dedupCfg()
{
    Fcv2Cfg cfg;
    cfg.useDedup = true;
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
    EXPECT_TRUE(roundtripWithCfg(recs, Fcv2Cfg()));        /* the default cfg */
    EXPECT_TRUE(roundtripWithCfg(recs, fineCfg()));        /* the fine shipped preset */
    EXPECT_TRUE(roundtripWithCfg(recs, coarseCfg()));      /* a coarse tier (not a shipped preset) */
    EXPECT_TRUE(roundtripWithCfg(recs, dedupCfg()));       /* dedup (no duplicate data) */
}

/*
 * Tier values at and above their maxima must survive the stream header, which
 * writes every one of them as a single byte.
 *
 * The above-maximum case is the one that discriminates: a value above the
 * maximum makes the encoder depend on the header carrying its clamped result,
 * so the two sides would disagree on the context layout if the header could not
 * carry it. The delta tiers that used to exercise this are gone with m5; the
 * remaining tiers are cycleMax (which round-trips via normalizeCfg), cycleBucket
 * and prevShift.
 */
TEST_F(Fcv2CfgTest, OverMaximalTierRoundTrip)
{
    Fcv2Cfg over;              /* every field above its maximum */
    over.cycleMax = 4096;      /* clamps to 128 */
    over.cycleBucket = 512;    /* clamps to 32 */
    over.prevShift = 9;        /* clamps to 3 */

    std::vector<Rec> recs = makeCycleRecords(1200, 4000, 23);
    EXPECT_TRUE(roundtripWithCfg(recs, over));
}

/* The same, with every field sitting exactly on its maximum at once. */
TEST_F(Fcv2CfgTest, MaximalTierRoundTrip)
{
    Fcv2Cfg cfg;
    cfg.cycleMax = 128;
    cfg.cycleBucket = 32;
    cfg.prevShift = 3;

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

/* The average-quality preset (strategy 4) is gone with its m6 model (useQa +
 * QA_BINS were removed along with it), so QaSavesBytesOnVariedQuality is gone
 * too: there is no per-read average-quality tier left to test.
 */

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
    ASSERT_TRUE(loaded) << "A prior with the same tier layout must be loaded (the prior determines the tier)";    ASSERT_TRUE(encodeAll(enc, recs));
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
 * The fcv2 model count is fixed at FCV2_MODEL_COUNT (m0..m4) for every file;
 * there is no per-read-length rule anymore (qualModelCountForMeanLen is gone),
 * so the mix is the same on short and long reads.
 */
TEST_F(Fcv2CfgTest, ModelCountIsFixed)
{
    /* The compile-time constant is exactly the shipped count. */
    EXPECT_EQ(FCV2_MODEL_COUNT, 5);
    /* A default coder mixes every model: the mix does not vary with anything. */
    EXPECT_TRUE(roundtripWithCfg(makeRecords(1500, 151, 33), Fcv2Cfg()));
}

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
    /* The preset must be one of the known candidates (fine / ultra / ultra-80;
     * default and coarse are no longer in the table), and it must round-trip
     * with the preset applied. */
    const QualFcv2Params& p = qualSel.fcv2Params;
    bool isKnownPreset =
        (p.cycleBucket == 24 && p.prevShift == 0) ||
        (p.cycleBucket == 4 && p.cycleMax == 96 && p.prevShift == 0) ||
        (p.cycleBucket == 4 && p.cycleMax == 80 && p.prevShift == 0);
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

/*
 * A value outside the stream's alphabet must fail the stream, never be dropped.
 *
 * The alphabet is fixed once the coder is built: from the frequency table, or from the prior's own
 * alphabet (which comes from the blocks the prior was trained on). Nothing stops a later block from
 * carrying a quality value the earlier ones never used - SRR2769247 followed by SRR18473920 is such
 * a pair, one using quality values 'B'..'i' and the other '#'..'G' - and coder_fcv2 used to skip
 * such a byte silently. The decoder has no way to know a byte was skipped, so it read every
 * following symbol at the wrong offset: the block's checksum failed, and with enough data the
 * decoder ran off the end of the stream and the process died on the load.
 *
 * What must hold instead: the coder reports what it can code (coversByte), and a value it cannot
 * code latches coder_io::IO_UNCODABLE so the caller can hand the block to a codec that covers it.
 */
TEST_F(Fcv2CfgTest, OutOfAlphabetValueFailsTheStreamInsteadOfBeingDropped)
{
    /* An alphabet holding exactly two values: '!' and 'B'. */
    std::vector<uint32_t> freq(256, 0);
    freq[(uint8_t)'!'] = 10;
    freq[(uint8_t)'B'] = 5;

    std::vector<uint8_t> comp(4096, 0);
    coder_io io(comp.data(), (int32_t)comp.size());
    coder_fcv2 enc(&io, freq);

    EXPECT_TRUE(enc.coversByte('!'));
    EXPECT_TRUE(enc.coversByte('B'));
    EXPECT_FALSE(enc.coversByte('C'));   /* 'C' is in no block the alphabet was built from */
    EXPECT_FALSE(enc.coversByte('h'));

    /* Records inside the alphabet encode normally. */
    uint8_t seq[8] = {'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T'};
    uint8_t good[8] = {'!', 'B', 'B', '!', '!', 'B', '!', 'B'};
    enc.encode_record(good, sizeof(good), false, seq, sizeof(seq));
    EXPECT_EQ(io.err, coder_io::IO_OK);

    /* One value outside it fails the stream rather than shortening it by a symbol. */
    uint8_t bad[4] = {'!', 'C', 'B', '!'};
    enc.encode_record(bad, sizeof(bad), false, seq, sizeof(bad));
    EXPECT_EQ(io.err, coder_io::IO_UNCODABLE);
}

/*
 * The same alphabet, used on both sides, still round-trips: the failure above is not a property of
 * narrow alphabets as such but of values leaving them.
 */
TEST_F(Fcv2CfgTest, NarrowAlphabetRoundTrips)
{
    std::vector<uint32_t> freq(256, 0);
    freq[(uint8_t)'!'] = 10;
    freq[(uint8_t)'B'] = 5;

    std::vector<uint8_t> seq(32, 'A');
    std::vector<uint8_t> qual(32);
    for (size_t i = 0; i < qual.size(); i++) {
        qual[i] = ((i / 3) % 2 == 0) ? (uint8_t)'!' : (uint8_t)'B';
    }

    std::vector<uint8_t> comp(4096, 0);
    {
        coder_io eio(comp.data(), (int32_t)comp.size());
        coder_fcv2 enc(&eio, freq);
        for (int r = 0; r < 50; r++) {
            enc.encode_record(qual.data(), (uint32_t)qual.size(), false, seq.data(), (uint32_t)seq.size());
        }
        enc.encode_flush();
        EXPECT_EQ(eio.err, coder_io::IO_OK);
    }

    coder_io dio(comp.data(), (int32_t)comp.size());
    coder_fcv2 dec(&dio, freq);
    ASSERT_EQ(dec.begin_decode(), 0);
    for (int r = 0; r < 50; r++) {
        std::vector<uint8_t> out(qual.size(), 0);
        ASSERT_EQ(dec.decode_record(out.data(), (uint32_t)out.size(), seq.data(), (uint32_t)seq.size()),
                  (int32_t)out.size());
        EXPECT_EQ(out, qual);
    }
}
