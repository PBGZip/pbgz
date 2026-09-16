/*
 * qual_tier_scan.cpp - Grid-scan fcv2 context parameters on real SAM QUAL data.
 *
 * For a fixed read count (simulating the "one block / one file of that size",
 * cold start, no prior) every parameter combination is compressed and
 * round-trip verified, and the result is reported against the tier the project
 * currently ships (Fcv2Cfg's defaults). Run it per volume of interest (e.g.
 * 10000 / 25000 / 100000 reads).
 *
 * The axis is selected on the command line because the interesting axes are not
 * the same on different data:
 *
 *   axis 0  cycleMax sweep. cycleMax is the hard bound of the positional context
 *           - cycleOf() clamps every cycle at or beyond it to cycleMax-1, and
 *           m1/m2/m4/m5/m6 all index by it - so a file whose reads are longer
 *           than cycleMax collapses all of its tail positions into one context
 *           slot. 151 bp reads therefore have 56 positions sharing a slot, and
 *           14 kb reads are almost entirely collapsed.
 *
 *   axis 1  cycleMax x cycleBucket. The two are coupled: cycleMax is how far
 *           the position context reaches, cycleBucket is how finely the reach
 *           is divided (it also sets the weight-row count, and m3's position
 *           resolution). Widening one without the other moves the sparsity
 *           somewhere else.
 *
 *   axis 2  deltaMax x deltaBucket, the same question for m5's transition-count
 *           context, whose reach is the number of transitions inside one read
 *           (32 by default).
 *
 *   axis 5  the shipped preset table (FCV2_PRESETS), unpruned.
 *
 *   axis 3  the full grid: cycleMax x (cycleBucket, deltaBucket) x prevShift.
 *           Expensive; use it once the narrower sweeps have picked a region.
 *
 * The coder clamps cycleMax to [32, 128] (CFG_CYCLE_MAX_MIN/MAX in
 * coder_fcv2.cpp) because the model arrays are sized by it, so this tool cannot
 * test beyond 128. If the winner sits exactly at 128, the cap is what binds and
 * raising it is a separate change (array sizes, prior blob size).
 *
 * Each run also compresses with coder_bwt_cm at the level the selector would
 * use, so "could fcv2 win here at all" is answerable on the spot - on long reads
 * bwt_cm beats fcv2 in the real selection.
 *
 * Usage: qual_tier_scan <input.sam> <max_lines> [axis 0..6]
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "coder/coder_io.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_bwt_cm.h"

namespace {

struct Rec {
    std::string qual;
    std::string seq;
    bool rev;
};

void loadRecords(const std::string& path, uint64_t maxLines, std::vector<Rec>& recs,
                 std::vector<uint32_t>& freqByByte)
{
    std::ifstream in(path);
    if (!in.is_open()) {
        fprintf(stderr, "cannot open %s\n", path.c_str());
        return;
    }
    freqByByte.assign(256, 0);
    std::string line;
    uint64_t lineNo = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '@')
            continue;
        ++lineNo;
        if (maxLines > 0 && lineNo > maxLines)
            break;

        std::string f[11];
        size_t start = 0, field = 0;
        for (size_t i = 0; i <= line.size() && field < 11; ++i) {
            if (i == line.size() || line[i] == '\t') {
                f[field] = line.substr(start, i - start);
                start = i + 1;
                ++field;
            }
        }
        if (field < 11)
            continue;

        Rec r;
        r.seq = f[9];
        r.qual = f[10];
        if (r.qual == "*" || r.qual.empty())
            continue;
        long flag = strtol(f[1].c_str(), nullptr, 10);
        r.rev = (flag & 0x10) != 0;
        for (char c : r.qual)
            freqByByte[(uint8_t)c]++;
        recs.push_back(std::move(r));
    }
    fprintf(stderr, "loaded %zu records (max_lines=%llu)\n", recs.size(),
            (unsigned long long)maxLines);
}

uint64_t qualBytesOf(const std::vector<Rec>& recs)
{
    uint64_t n = 0;
    for (const Rec& r : recs)
        n += r.qual.size();
    return n;
}

bool encodeFcv2(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq,
                const Fcv2Cfg& cfg, uint32_t& dstLen)
{
    if (recs.empty())
        return false;
    uint64_t qualBytes = qualBytesOf(recs);
    std::vector<uint8_t> comp((size_t)qualBytes * 2 + (1u << 20), 0);
    {
        coder_io io(comp.data(), (int32_t)comp.size());
        coder_fcv2 coder(&io, freq, cfg);
        for (const Rec& r : recs) {
            coder.encode_record((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        int32_t n = coder.encode_flush();
        if (io.err != coder_io::IO_OK || n <= 0)
            return false;
        dstLen = (uint32_t)n;
    }
    /* Round-trip verify: the tier travels in the stream header, so a decoder
       built from the defaults must still reproduce the input. */
    {
        coder_io io(comp.data(), dstLen);
        coder_fcv2 dec(&io, freq, Fcv2Cfg());
        if (dec.begin_decode() < 0)
            return false;
        for (const Rec& r : recs) {
            std::vector<uint8_t> dst(r.qual.size(), 0);
            int32_t got = dec.decode_record(dst.data(), (uint32_t)r.qual.size(),
                                            (const uint8_t*)r.seq.data(),
                                            (uint32_t)r.seq.size());
            if (got != (int32_t)r.qual.size() ||
                memcmp(dst.data(), r.qual.data(), r.qual.size()) != 0) {
                return false;
            }
        }
    }
    return true;
}

/* Internal block sizes of coder_bwt_cm, same table as BWT_LEVEL_SIZE in
 * codec_selector.cpp; the selector picks the smallest block that holds the
 * sample, and the reference row has to match that choice to be comparable. */
int bwtLevelFor(uint32_t sampleLen)
{
    static const uint32_t kSize[10] = {
        0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
        1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
    };
    for (int lv = 1; lv <= 9; ++lv) {
        if (kSize[lv] >= sampleLen)
            return lv;
    }
    return 9;
}

bool encodeBwtCm(const std::vector<Rec>& recs, uint32_t& dstLen)
{
    uint64_t qualBytes = qualBytesOf(recs);
    if (recs.empty() || qualBytes == 0)
        return false;
    std::vector<uint8_t> buf((size_t)qualBytes * 2 + (1u << 20), 0);
    coder_io io(buf.data(), (int32_t)buf.size());
    io.set_level(bwtLevelFor((uint32_t)qualBytes));
    coder_bwt_cm coder(&io);
    for (const Rec& r : recs) {
        if (r.qual.empty())
            continue;
        coder.encode_line((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size());
    }
    coder.encode_flush();
    int32_t n = io.data_len;
    if (io.err != coder_io::IO_OK || n <= 0)
        return false;
    dstLen = (uint32_t)n;
    return true;
}

std::string nameOf(const Fcv2Cfg& c)
{
    char buf[96];
    snprintf(buf, sizeof(buf), "cm%3d cb%2d dm%3d db%2d ps%d",
             c.cycleMax, c.cycleBucket, c.deltaMax, c.deltaBucket, c.prevShift);
    return buf;
}

/* The preset a configuration corresponds to, if any, so a row is recognizable by
   the name the selector and the logs use. */
const char* labelFor(const Fcv2Cfg& c)
{
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        const Fcv2Cfg& p = FCV2_PRESETS[i].cfg;
        if (c.cycleMax == p.cycleMax && c.cycleBucket == p.cycleBucket &&
            c.deltaMax == p.deltaMax && c.deltaBucket == p.deltaBucket &&
            c.prevShift == p.prevShift) {
            return FCV2_PRESETS[i].name;
        }
    }
    return nullptr;
}

/*
 * The QUAL feature vector a selection rule would have to work from: the read
 * length and the shape of the value distribution. Printed with every scan so the
 * measurements can be read against the features they came from.
 */
struct QualFeatures {
    double   meanLen;
    double   meanQual;      /* mean raw byte value */
    int      alphabet;
    int      minByte;
    int      maxByte;
    int      phredBase;     /* 33 or 64, inferred from the lowest value seen */
    double   fracHigh;      /* share of values at or above the 90th percentile byte */
};

QualFeatures featuresOf(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq)
{
    QualFeatures f = {};
    uint64_t bytes = 0;
    for (size_t b = 0; b < 256; ++b) {
        if (freq[b] == 0) {
            continue;
        }
        if (f.alphabet == 0) {
            f.minByte = (int)b;
        }
        f.maxByte = (int)b;
        ++f.alphabet;
        bytes += freq[b];
        f.meanQual += (double)b * (double)freq[b];
    }
    f.meanLen = bytes ? (double)bytes / (double)recs.size() : 0.0;
    f.meanQual = bytes ? f.meanQual / (double)bytes : 0.0;
    f.phredBase = (f.minByte >= 64) ? 64 : 33;
    /* Any of the highest values, e.g. whether the platform emits a high tail. */
    const uint64_t special = (f.maxByte > f.phredBase) ? freq[f.maxByte] : 0;
    f.fracHigh = bytes ? (double)special / (double)bytes : 0.0;
    return f;
}

struct Result {
    Fcv2Cfg cfg;
    uint32_t len;
    bool ok;
    bool isBwt;
    const char* label;   /* preset name, or nullptr to use the parameter string */
};

/* The tiers the project ships today, as anchors in every run. */
void pushShippedTiers(std::vector<Fcv2Cfg>& grid)
{
    grid.push_back(Fcv2Cfg());                        /* default */
    Fcv2Cfg ultra;                                    /* ultra 4/2/ps0 */
    ultra.cycleBucket = 4;
    ultra.deltaBucket = 2;
    ultra.prevShift = 0;
    grid.push_back(ultra);
    Fcv2Cfg coarse;                                   /* coarse 8/4/ps2 */
    coarse.cycleBucket = 8;
    coarse.deltaBucket = 4;
    coarse.prevShift = 2;
    grid.push_back(coarse);
    Fcv2Cfg fine;                                     /* fine 24/12/ps0 */
    fine.cycleBucket = 24;
    fine.deltaBucket = 12;
    fine.prevShift = 0;
    grid.push_back(fine);
}

std::vector<Fcv2Cfg> buildGrid(int axis)
{
    std::vector<Fcv2Cfg> grid;
    const int kCycleMax[] = {32, 48, 64, 80, 96, 112, 128};

    if (axis == 0) {
        for (int cm : kCycleMax) {
            Fcv2Cfg c;
            c.cycleMax = cm;
            grid.push_back(c);
        }
    } else if (axis == 1) {
        const int kBuckets[] = {4, 8, 16, 24, 32};
        for (int cm : kCycleMax) {
            for (int cb : kBuckets) {
                Fcv2Cfg c;
                c.cycleMax = cm;
                c.cycleBucket = cb;
                grid.push_back(c);
            }
        }
    } else if (axis == 2) {
        /* 255 is the largest representable deltaMax: the field is one byte in the
         * stream header, so 256 would be truncated on write. */
        const int kDeltaMax[] = {8, 16, 32, 64, 128, 255};
        const int kDeltaBucket[] = {4, 8, 16};
        for (int dm : kDeltaMax) {
            for (int db : kDeltaBucket) {
                Fcv2Cfg c;
                c.deltaMax = dm;
                c.deltaBucket = db;
                grid.push_back(c);
            }
        }
    } else if (axis == 6) {
        /*
         * Long-read hypotheses, on the observation that on long-read QUAL the
         * best cycleMax is the *smallest* tried (32) and the whole cycleMax plane
         * only moves 0.05%: the positional context is not what limits fcv2 there,
         * so refining it is not the lever. Two other contexts degenerate on a
         * 14 kb record and are what this axis probes:
         *
         *   m5 (delta = transitions so far). Its bucket is
         *   delta * deltaBucket / deltaMax, so with the default deltaMax=32 a
         *   read with thousands of transitions sits in the last bucket from very
         *   early on and the model carries no information. Raising deltaMax (up
         *   to the 256 cap) is the only way to give it range.
         *
         *   m6 (read mean quality, QA_BINS=4). It is degenerate while useQa is
         *   off - its context equals m1's - and no shipped candidate ever turns
         *   it on. On long reads the per-read mean quality varies much more than
         *   on a short-read flowcell, so it is worth measuring rather than
         *   assuming, which is what the qa variant here does.
         *
         * Crossed with the small cycleMax values that measured best there.
         */
        const int kLongCycleMax[] = {16, 24, 32, 48, 64};
        const int kLongDeltaMax[] = {32, 64, 128, 255};
        const bool kQa[] = {false, true};
        for (int cm : kLongCycleMax) {
            for (int dm : kLongDeltaMax) {
                for (bool qa : kQa) {
                    Fcv2Cfg c;
                    c.cycleMax = cm;
                    c.cycleBucket = 8;
                    c.deltaMax = dm;
                    c.deltaBucket = 16;
                    c.prevShift = 0;
                    c.useQa = qa;
                    grid.push_back(c);
                }
            }
        }
    } else if (axis == 7) {
        /*
         * The order-2 question. A second-order context is already in the mix (m3 uses the previous
         * two quality values), but it is quantized by prevShift and crossed with the cycle bucket,
         * so the axes above measure it diluted. The conditional entropy of this data says the
         * second and third predecessors still carry ~0.02 and ~0.06 bit/symbol - about 0.5% to
         * 1.5% of the QUAL column - so the question is whether the models that exist can reach it
         * when nothing is quantized away: prevShift 0 with a coarse cycle context, a combination no
         * other axis tries (axis 1 sweeps cycleBucket at prevShift 1, axis 6 fixes it at 8).
         */
        const int kOrder2CycleBucket[] = {2, 3, 4, 6, 8};
        const int kOrder2DeltaBucket[] = {4, 8, 16};
        const int kOrder2CycleMax[] = {32, 48};
        for (int cm : kOrder2CycleMax) {
            for (int cb : kOrder2CycleBucket) {
                for (int db : kOrder2DeltaBucket) {
                    Fcv2Cfg c;
                    c.cycleMax = cm;
                    c.cycleBucket = cb;
                    c.deltaMax = 32;
                    c.deltaBucket = db;
                    c.prevShift = 0;
                    c.useDelta = true;
                    c.useDedup = false;
                    c.useQa = false;
                    grid.push_back(c);
                }
            }
        }
    } else if (axis == 5) {
        /*
         * Exactly the shipped preset table (FCV2_PRESETS), without the
         * volume-tier pruning the selector applies: the point is to rank the
         * presets themselves on this data, including the ones the selector would
         * not bother trying at this volume.
         */
        for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
            grid.push_back(FCV2_PRESETS[i].cfg);
        }
    } else if (axis == 4) {
        /*
         * Finalists: the winners of axes 0-3 plus the shipped ultra tier, at
         * several prevShift values. Cheap enough to re-run at every block volume,
         * which is what decides whether a winner found on one sample is the
         * actual optimum or an artifact of that sample size.
         */
        Fcv2Cfg ultra;
        ultra.cycleBucket = 4;
        ultra.deltaBucket = 2;
        ultra.prevShift = 0;
        grid.push_back(ultra);                       /* shipped ultra, cm96 */
        for (int ps = 0; ps <= 2; ++ps) {
            Fcv2Cfg a;                               /* cm80 cb4 db2 */
            a.cycleMax = 80;
            a.cycleBucket = 4;
            a.deltaBucket = 2;
            a.prevShift = ps;
            grid.push_back(a);
            Fcv2Cfg b;                               /* cm80 cb8 db4 */
            b.cycleMax = 80;
            b.cycleBucket = 8;
            b.deltaBucket = 4;
            b.prevShift = ps;
            grid.push_back(b);
            Fcv2Cfg c;                               /* long-read region: small cycleMax */
            c.cycleMax = 32;
            c.cycleBucket = 8;
            c.deltaBucket = 8;
            c.prevShift = ps;
            grid.push_back(c);
        }
    } else {
        const int kBucketPairs[][2] = {{4, 2}, {8, 4}, {16, 8}, {24, 12}};
        for (int cm : kCycleMax) {
            for (const int* p : kBucketPairs) {
                for (int ps = 0; ps <= 2; ++ps) {
                    Fcv2Cfg c;
                    c.cycleMax = cm;
                    c.cycleBucket = p[0];
                    c.deltaBucket = p[1];
                    c.prevShift = ps;
                    grid.push_back(c);
                }
            }
        }
    }
    return grid;
}

const char* axisName(int axis)
{
    switch (axis) {
    case 0: return "cycleMax sweep (cycleBucket/deltaBuckets/prevShift at the shipped defaults)";
    case 1: return "cycleMax x cycleBucket";
    case 2: return "deltaMax x deltaBucket";
    case 4: return "finalists from axes 0-3, per prevShift (for every block volume)";
    case 5: return "the shipped preset table (FCV2_PRESETS), unpruned";
    case 6: return "long-read hypotheses: small cycleMax x deltaMax x useQa";
    case 7: return "order-2 focus: prevShift 0 with coarse cycle buckets";
    default: return "full grid: cycleMax x (cycleBucket, deltaBucket) x prevShift";
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: %s <input.sam> <max_lines> [axis 0..6]\n", argv[0]);
        return 2;
    }
    const uint64_t maxLines = strtoull(argv[2], NULL, 10);
    const int axis = (argc > 3) ? atoi(argv[3]) : 0;

    std::vector<Rec> recs;
    std::vector<uint32_t> freq;
    loadRecords(argv[1], maxLines, recs, freq);
    if (recs.empty()) {
        fprintf(stderr, "no quality records\n");
        return 1;
    }
    const uint64_t qualBytes = qualBytesOf(recs);
    const QualFeatures feat = featuresOf(recs, freq);

    std::vector<Fcv2Cfg> grid = buildGrid(axis);
    std::vector<Fcv2Cfg> shipped;
    pushShippedTiers(shipped);

    printf("=== fcv2 tier scan: %s ===\n", argv[1]);
    printf("reads %zu, QUAL %llu bytes\n", recs.size(), (unsigned long long)qualBytes);
    printf("features: mean length %.1f | alphabet %d | byte range %d..%d (Phred+%d) | "
           "mean byte %.1f (phred %.1f) | top-value share %.4f\n",
           feat.meanLen, feat.alphabet, feat.minByte, feat.maxByte, feat.phredBase,
           feat.meanQual, feat.meanQual - feat.phredBase, feat.fracHigh);
    printf("axis %d: %s\n", axis, axisName(axis));

    std::vector<Result> results;
    results.reserve(grid.size() + shipped.size() + 1);

    for (const Fcv2Cfg& cfg : grid) {
        uint32_t len = 0;
        bool ok = encodeFcv2(recs, freq, cfg, len);
        results.push_back({cfg, len, ok, false, labelFor(cfg)});
        const std::string nm = nameOf(cfg);
        fprintf(stderr, "  %-28s %10u %s\n",
                (labelFor(cfg) != nullptr) ? labelFor(cfg) : nm.c_str(), len, ok ? "OK" : "FAIL");
    }

    uint32_t base = 0;
    {
        uint32_t len = 0;
        base = encodeFcv2(recs, freq, Fcv2Cfg(), len) ? len : 0;
    }
    uint32_t bwtLen = 0;
    const bool bwtOk = encodeBwtCm(recs, bwtLen);

    std::stable_sort(results.begin(), results.end(),
                     [](const Result& a, const Result& b) {
                         if (a.ok != b.ok) return a.ok;
                         return a.len < b.len;
                     });

    /* Rows carry their preset name where they correspond to one. */
    struct RowName {
        const Result& r;
        std::string text() const { return r.label ? std::string(r.label) : nameOf(r.cfg); }
    };

    printf("\n%-28s %12s %10s  %s\n", "cfg", "compressed", "bits/byte", "vs shipped default");
    for (const Result& r : results) {
        const RowName nm{r};
        if (!r.ok) {
            printf("%-28s %12s %10s  FAIL\n", nm.text().c_str(), "-", "-");
            continue;
        }
        printf("%-28s %12u %10.4f  %+.3f%%\n", nm.text().c_str(), r.len,
               qualBytes ? (double)r.len * 8.0 / (double)qualBytes : 0.0,
               (base > 0) ? (double)((int64_t)r.len - (int64_t)base) * 100.0 / (double)base : 0.0);
    }

    /* On the preset axis the rows above already are the reference set. */
    if (axis != 5) {
        printf("\nshipped tiers for reference:\n");
        for (const Fcv2Cfg& c : shipped) {
            uint32_t len = 0;
            bool ok = encodeFcv2(recs, freq, c, len);
            if (!ok) {
                printf("  %-26s FAIL\n", nameOf(c).c_str());
                continue;
            }
            printf("  %-26s %12u %10.4f  %+.3f%%\n", nameOf(c).c_str(), len,
                   qualBytes ? (double)len * 8.0 / (double)qualBytes : 0.0,
                   (base > 0) ? (double)((int64_t)len - (int64_t)base) * 100.0 / (double)base : 0.0);
        }
    }
    if (bwtOk) {
        printf("\n  %-26s %12u %10.4f  %+.3f%%\n", "coder_bwt_cm (reference)", bwtLen,
               qualBytes ? (double)bwtLen * 8.0 / (double)qualBytes : 0.0,
               (base > 0) ? (double)((int64_t)bwtLen - (int64_t)base) * 100.0 / (double)base : 0.0);
    } else {
        printf("\n  coder_bwt_cm reference failed\n");
    }
    printf("\nNote: fcv2 clamps cycleMax to [32,128]; the cap bounds this scan.\n");
    return 0;
}
