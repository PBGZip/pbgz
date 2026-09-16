/*
 * fcv2_len_train.cpp - train fcv2's context parameters and mixing-model count per read-length range.
 *
 * Why this exists
 * ---------------
 * fcv2 is the QUAL coder built for short reads; on long reads it currently loses to coder_bwt_cm,
 * which is what makes QUAL the one column where pbgz trails CRAM on ONT data. Two decisions are
 * behind that, and both are rules rather than measurements:
 *
 *   - the mixing-model count follows the mean read length, one threshold at 500 bytes: 4 models
 *     above it, 6 below (qualModelCountForMeanLen in qual_selector.cpp);
 *   - the context parameters come from five hand-picked presets (FCV2_PRESETS in coder_fcv2.h),
 *     trialled as they are.
 *
 * Read length is the variable that separates the regimes, so the honest way to replace those rules
 * is to measure the encoded size of real records *per read-length bucket* and pick the
 * configuration that wins inside each bucket. That is what this tool does; the table it prints is
 * meant to be pasted back as the per-bucket rule.
 *
 * How the search works
 * --------------------
 * The objective is a black box - one full fcv2 encode of the bucket's sample, measured in bytes -
 * so the search is a coordinate descent over the discrete knobs (each coordinate is stepped
 * through a ladder of candidate values, the best value is kept, and the sweep repeats until
 * nothing improves) plus a sweep over the model count, started from several points (every shipped
 * preset and the current per-length rule) so a local minimum of one start does not decide the
 * answer on its own. A time weight can be added to the objective when the point is speed rather
 * than size.
 *
 * Note where gradient descent really happens: fcv2's mixing weights are trained by gradient
 * descent inside the coder, at encode time, on every stream. What is trained here is the
 * configuration around them, which is discrete and therefore descended rather than differentiated.
 *
 * Usage
 * -----
 *   fcv2_len_train <reads.sam> [options]
 *
 *     --max-lines N        read at most N records from the input (0 = all)
 *     --buckets a,b,c      read-length bucket edges in bytes, default 0,100,300,500,1000,3000,10000
 *     --records-per-bucket N  cap the sample per bucket by record count (default 2000)
 *     --bytes-per-bucket N    and by size, whichever binds first (default 64 MB)
 *     --iters K            descent sweeps per bucket (default 3)
 *     --starts N           how many starting configurations to descend from (default all)
 *     --time-us-weight W   bytes added per microsecond of encode time (default 0 = size only)
 *     --no-search          skip the descent: only measure the shipped default, coder_bwt_cm and
 *                          coder_ppmd on each bucket (usable with a large sample)
 *     --quiet              only print the result table
 *
 * On the sample caps: record count matters more than byte count here. fcv2's contexts are trained
 * inside the coder as it encodes, so a bucket sampled to a fixed number of *bytes* gives its long
 * reads only a handful of records (the ONT buckets below 100 KB of QUAL are four records) and the
 * descent then overfits that handful - cycleMax collapses to 32 because nothing has shown the coder
 * a longer read yet. That is also why the product's QUAL trial, which budgets its sample in bytes,
 * keeps judging fcv2 on ONT data it cannot have converged on; see --records-per-bucket.
 *
 * The input is SAM text (samtools view x.bam > x.sam), because the tool needs the QUAL column and
 * the read length, and record boundaries are what the coder is driven by.
 */

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "coder/coder_bwt_cm.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_io.h"
#include "coder/coder_ppmd.h"

namespace {

struct Rec {
    std::string qual;
    std::string seq;
    bool rev;
};

void loadRecords(const std::string& path, uint64_t maxLines, std::vector<Rec>& recs)
{
    std::ifstream in(path);
    if (!in.is_open()) {
        fprintf(stderr, "cannot open %s\n", path.c_str());
        return;
    }
    std::string line;
    uint64_t lineNo = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '@') {
            continue;
        }
        ++lineNo;
        if (maxLines > 0 && lineNo > maxLines) {
            break;
        }

        std::string f[11];
        size_t start = 0, field = 0;
        for (size_t i = 0; i <= line.size() && field < 11; ++i) {
            if (i == line.size() || line[i] == '\t') {
                f[field] = line.substr(start, i - start);
                start = i + 1;
                ++field;
            }
        }
        if (field < 11) {
            continue;
        }

        Rec r;
        r.seq = f[9];
        r.qual = f[10];
        if (r.qual == "*" || r.qual.empty()) {
            continue;
        }
        const long flag = strtol(f[1].c_str(), nullptr, 10);
        r.rev = (flag & 0x10) != 0;
        recs.push_back(std::move(r));
    }
    fprintf(stderr, "loaded %zu records (max_lines=%llu)\n", recs.size(),
            (unsigned long long)maxLines);
}

/*
 * The frequency table fcv2 is given is a ranked alphabet, not raw counts: the product passes
 * freqByByte[byte] = alphabetSize - rank, most frequent byte first (see
 * SamCodecActuator::compressQuality, which builds it from qualFreqTable). The coder grows its
 * starting Huffman tree from it, so passing counts instead is not a harmless scale change - with
 * counts over a multi-hundred-megabyte sample the tree's accumulator overflows and every encode
 * fails. Built per bucket here, because one bucket stands in for one block of a file whose reads
 * are that long, which is the volume the product would build it from.
 */
void buildRankFreq(const std::vector<Rec>& recs, std::vector<uint32_t>& freq)
{
    std::vector<uint64_t> count(256, 0);
    for (const Rec& r : recs) {
        for (char c : r.qual) {
            count[(uint8_t)c]++;
        }
    }
    std::vector<std::pair<uint64_t, uint32_t>> present;
    for (uint32_t b = 0; b < 256; ++b) {
        if (count[b] > 0) {
            present.push_back(std::make_pair(count[b], b));
        }
    }
    std::sort(present.begin(), present.end(),
              [](const std::pair<uint64_t, uint32_t>& a, const std::pair<uint64_t, uint32_t>& b) {
                  if (a.first != b.first) {
                      return a.first > b.first;
                  }
                  return a.second < b.second;
              });
    freq.assign(256, 0);
    for (size_t i = 0; i < present.size(); ++i) {
        freq[present[i].second] = (uint32_t)(present.size() - i);
    }
}

uint64_t qualBytesOf(const std::vector<Rec>& recs)
{
    uint64_t n = 0;
    for (const Rec& r : recs) {
        n += r.qual.size();
    }
    return n;
}

/* Encode the bucket with one configuration and return its compressed size and wall time.
   The decode side is verified on the winner only (see verifyRoundTrip): running it on every
   evaluation would triple the cost of the search. */
bool encodeFcv2(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq, const Fcv2Cfg& cfg,
                uint32_t& dstLen, uint64_t& usec)
{
    if (recs.empty()) {
        return false;
    }
    const uint64_t qualBytes = qualBytesOf(recs);
    std::vector<uint8_t> comp((size_t)qualBytes * 2 + (1u << 20), 0);
    const auto t0 = std::chrono::steady_clock::now();
    {
        coder_io io(comp.data(), (int32_t)comp.size());
        coder_fcv2 coder(&io, freq, cfg);
        for (const Rec& r : recs) {
            coder.encode_record((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        const int32_t n = coder.encode_flush();
        if (io.err != coder_io::IO_OK || n <= 0) {
            return false;
        }
        dstLen = (uint32_t)n;
    }
    usec = (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t0).count();
    return true;
}

/* Decode with a decoder built from the defaults - the parameters travel in the stream header, so
   this fails if the configuration did not make it there. */
bool verifyRoundTrip(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq,
                     const Fcv2Cfg& cfg)
{
    const uint64_t qualBytes = qualBytesOf(recs);
    std::vector<uint8_t> comp((size_t)qualBytes * 2 + (1u << 20), 0);
    uint32_t dstLen = 0;
    {
        coder_io io(comp.data(), (int32_t)comp.size());
        coder_fcv2 coder(&io, freq, cfg);
        for (const Rec& r : recs) {
            coder.encode_record((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        const int32_t n = coder.encode_flush();
        if (io.err != coder_io::IO_OK || n <= 0) {
            return false;
        }
        dstLen = (uint32_t)n;
    }
    coder_io io(comp.data(), (int32_t)dstLen);
    coder_fcv2 dec(&io, freq, Fcv2Cfg());
    if (dec.begin_decode() < 0) {
        return false;
    }
    for (const Rec& r : recs) {
        std::vector<uint8_t> dst(r.qual.size(), 0);
        const int32_t got = dec.decode_record(dst.data(), (uint32_t)r.qual.size(),
                                              (const uint8_t*)r.seq.data(),
                                              (uint32_t)r.seq.size());
        if (got != (int32_t)r.qual.size() ||
            memcmp(dst.data(), r.qual.data(), r.qual.size()) != 0) {
            return false;
        }
    }
    return true;
}

/* The coder_bwt_cm reference row: the coder fcv2 has to beat for the QUAL decision. */
bool encodeBwtCm(const std::vector<Rec>& recs, uint32_t& dstLen)
{
    static const uint32_t kSize[10] = {
        0, 1u << 20, 1u << 22, 1u << 23, 0x00FFFFFFu,
        1u << 25, 1u << 26, 1u << 27, 1u << 28, 0x7FFFFFFFu
    };
    const uint64_t qualBytes = qualBytesOf(recs);
    if (recs.empty() || qualBytes == 0) {
        return false;
    }
    int bwtLevel = 9;
    for (int lv = 1; lv <= 9; ++lv) {
        if (kSize[lv] >= qualBytes) {
            bwtLevel = lv;
            break;
        }
    }
    std::vector<uint8_t> buf((size_t)qualBytes * 2 + (1u << 20), 0);
    coder_io io(buf.data(), (int32_t)buf.size());
    io.set_level(bwtLevel);
    coder_bwt_cm coder(&io);
    for (const Rec& r : recs) {
        if (!r.qual.empty()) {
            coder.encode_line((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size());
        }
    }
    coder.encode_flush();
    const int32_t n = io.data_len;
    if (io.err != coder_io::IO_OK || n <= 0) {
        return false;
    }
    dstLen = (uint32_t)n;
    return true;
}

/*
 * The other candidate worth measuring against fcv2 and bwt_cm on QUAL: coder_ppmd, a PPMd8
 * (order 6, 8 MB) context model that the product only uses for the reference blob today. It is
 * driven per record like bwt_cm - the QUAL stream is length-delimited, each record's length is
 * known on both sides - and it needs no alphabet or frequency table. Its round trip is verified
 * here as well, because that is the property the product would depend on.
 */
bool encodePpmd(const std::vector<Rec>& recs, uint32_t& dstLen, bool verify)
{
    const uint64_t qualBytes = qualBytesOf(recs);
    if (recs.empty() || qualBytes == 0) {
        return false;
    }
    std::vector<uint8_t> buf((size_t)qualBytes * 2 + (1u << 20), 0);
    Json::Value header;   /* info + fnlen, written by the encoder and read back by the decoder */
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_ppmd coder(&io);
        for (const Rec& r : recs) {
            if (!r.qual.empty()) {
                coder.encode((const uint8_t*)r.qual.data(), (int64_t)r.qual.size());
            }
        }
        coder.encode_flush();
        if (io.err != coder_io::IO_OK || io.data_len <= 0) {
            return false;
        }
        dstLen = (uint32_t)io.data_len;
        header = io.meta;
    }
    if (verify) {
        coder_io io(buf.data(), (int32_t)dstLen);
        io.meta = header;
        coder_ppmd dec(&io);
        for (const Rec& r : recs) {
            std::vector<uint8_t> dst(r.qual.size(), 0);
            if (dec.decode(dst.data(), (int64_t)r.qual.size()) != (int64_t)r.qual.size() ||
                memcmp(dst.data(), r.qual.data(), r.qual.size()) != 0) {
                return false;
            }
        }
    }
    return true;
}

std::string nameOf(const Fcv2Cfg& c)
{
    char buf[128];
    snprintf(buf, sizeof(buf), "cm%3d cb%2d dm%3d db%2d ps%d delta%d dedup%d qa%d mc%d",
             c.cycleMax, c.cycleBucket, c.deltaMax, c.deltaBucket, c.prevShift,
             (int)c.useDelta, (int)c.useDedup, (int)c.useQa, c.modelCount);
    return buf;
}

/* The shipped preset a configuration corresponds to, so a trained row stays recognizable. */
std::string presetNameOf(const Fcv2Cfg& c)
{
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        const Fcv2Cfg& p = FCV2_PRESETS[i].cfg;
        if (c.cycleMax == p.cycleMax && c.cycleBucket == p.cycleBucket &&
            c.deltaMax == p.deltaMax && c.deltaBucket == p.deltaBucket &&
            c.prevShift == p.prevShift && c.useDelta == p.useDelta &&
            c.useDedup == p.useDedup && c.useQa == p.useQa && c.modelCount == p.modelCount) {
            return FCV2_PRESETS[i].name;
        }
    }
    return "-";
}

struct Score {
    bool ok = false;
    uint32_t bytes = 0;
    uint64_t usec = 0;
    double value = 0.0;   /* bytes + timeUsWeight * usec: the descent's objective */
};

Score evaluate(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq, const Fcv2Cfg& cfg,
               double timeUsWeight)
{
    Score s;
    s.ok = encodeFcv2(recs, freq, cfg, s.bytes, s.usec);
    if (s.ok) {
        s.value = (double)s.bytes + timeUsWeight * (double)s.usec;
    }
    return s;
}

/* Ladders of the values each coordinate is stepped through. The ranges come from the shipped
   presets (cycleMax 80/96, cycleBucket 4/8/16/24, deltaBucket 2/4/8/12, prevShift 0/1/2) widened
   one step past both ends, so the search can leave the hand-picked set but stays inside what the
   coder's tables are sized for. */
const int kCycleMaxValues[]     = {32, 48, 64, 80, 96, 128, 160, 256};
const int kCycleBucketValues[]  = {2, 3, 4, 6, 8, 12, 16, 24, 32};
const int kDeltaMaxValues[]     = {16, 24, 32, 48, 64};
const int kDeltaBucketValues[]  = {1, 2, 3, 4, 6, 8, 12, 16};
const int kPrevShiftValues[]    = {0, 1, 2, 3};
const int kModelCountValues[]   = {4, 5, 6, 7};

bool stepCoordinate(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq, Fcv2Cfg& cfg,
                    double timeUsWeight, double& bestValue, const char* what,
                    const int* values, size_t valueCount, int Fcv2Cfg::*field, bool quiet)
{
    bool improved = false;
    for (size_t i = 0; i < valueCount; ++i) {
        if (values[i] == cfg.*field) {
            continue;
        }
        Fcv2Cfg trial = cfg;
        trial.*field = values[i];
        const Score s = evaluate(recs, freq, trial, timeUsWeight);
        if (s.ok && s.value < bestValue) {
            bestValue = s.value;
            cfg = trial;
            improved = true;
            if (!quiet) {
                fprintf(stderr, "    %s -> %d : %.0f\n", what, values[i], s.value);
            }
        }
    }
    return improved;
}

bool stepBool(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq, Fcv2Cfg& cfg,
              double timeUsWeight, double& bestValue, const char* what, bool Fcv2Cfg::*field,
              bool quiet)
{
    Fcv2Cfg trial = cfg;
    trial.*field = !(cfg.*field);
    const Score s = evaluate(recs, freq, trial, timeUsWeight);
    if (s.ok && s.value < bestValue) {
        bestValue = s.value;
        cfg = trial;
        if (!quiet) {
            fprintf(stderr, "    %s -> %d : %.0f\n", what, (int)(cfg.*field), s.value);
        }
        return true;
    }
    return false;
}

/* One descent sweep: every coordinate in turn, each one kept only if it lowers the objective.
   Returns true when something moved, so the caller can sweep again. */
bool descend(const std::vector<Rec>& recs, const std::vector<uint32_t>& freq, Fcv2Cfg& cfg,
             double& bestValue, double timeUsWeight, bool quiet)
{
    bool moved = false;
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "cycleMax",
                            kCycleMaxValues, sizeof(kCycleMaxValues) / sizeof(int),
                            &Fcv2Cfg::cycleMax, quiet);
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "cycleBucket",
                            kCycleBucketValues, sizeof(kCycleBucketValues) / sizeof(int),
                            &Fcv2Cfg::cycleBucket, quiet);
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "deltaMax",
                            kDeltaMaxValues, sizeof(kDeltaMaxValues) / sizeof(int),
                            &Fcv2Cfg::deltaMax, quiet);
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "deltaBucket",
                            kDeltaBucketValues, sizeof(kDeltaBucketValues) / sizeof(int),
                            &Fcv2Cfg::deltaBucket, quiet);
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "prevShift",
                            kPrevShiftValues, sizeof(kPrevShiftValues) / sizeof(int),
                            &Fcv2Cfg::prevShift, quiet);
    moved |= stepCoordinate(recs, freq, cfg, timeUsWeight, bestValue, "modelCount",
                            kModelCountValues, sizeof(kModelCountValues) / sizeof(int),
                            &Fcv2Cfg::modelCount, quiet);
    moved |= stepBool(recs, freq, cfg, timeUsWeight, bestValue, "useDelta", &Fcv2Cfg::useDelta,
                      quiet);
    moved |= stepBool(recs, freq, cfg, timeUsWeight, bestValue, "useDedup", &Fcv2Cfg::useDedup,
                      quiet);
    moved |= stepBool(recs, freq, cfg, timeUsWeight, bestValue, "useQa", &Fcv2Cfg::useQa, quiet);
    return moved;
}

struct BucketResult {
    uint32_t loBytes = 0;
    uint32_t hiBytes = 0;          /* 0 = open-ended */
    uint64_t reads = 0;
    uint64_t qualBytes = 0;        /* whole bucket, for context */
    uint64_t sampleQualBytes = 0;  /* what the ratios are measured against */
    uint32_t meanLen = 0;
    Fcv2Cfg cfg;
    uint32_t bytes = 0;
    uint64_t usec = 0;
    uint32_t baseBytes = 0;        /* the shipped default preset, as the starting point */
    uint64_t baseUsec = 0;
    uint32_t bwtBytes = 0;         /* coder_bwt_cm, the coder fcv2 has to beat */
    uint32_t ppmdBytes = 0;        /* coder_ppmd: the other candidate for a QUAL decision */
    bool ppmdOk = false;
    bool roundTripOk = false;
};

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "usage: %s <reads.sam> [--max-lines N] [--buckets a,b,...] "
                "[--bytes-per-bucket N] [--iters K] [--time-us-weight W] [--quiet]\n",
                argv[0]);
        return 2;
    }

    std::string path = argv[1];
    uint64_t maxLines = 0;
    uint64_t bytesPerBucket = 64u << 20;
    uint64_t recordsPerBucket = 2000;
    int iters = 3;
    int startsLimit = 0;   /* 0 = all */
    double timeUsWeight = 0.0;
    bool quiet = false;
    bool noSearch = false;
    std::vector<uint32_t> edges = {0, 100, 300, 500, 1000, 3000, 10000};

    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](void) -> std::string {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s needs a value\n", a.c_str());
                exit(2);
            }
            return argv[++i];
        };
        if (a == "--max-lines") {
            maxLines = strtoull(next().c_str(), nullptr, 10);
        } else if (a == "--bytes-per-bucket") {
            bytesPerBucket = strtoull(next().c_str(), nullptr, 10);
        } else if (a == "--records-per-bucket") {
            recordsPerBucket = strtoull(next().c_str(), nullptr, 10);
        } else if (a == "--starts") {
            startsLimit = atoi(next().c_str());
        } else if (a == "--iters") {
            iters = atoi(next().c_str());
        } else if (a == "--time-us-weight") {
            timeUsWeight = atof(next().c_str());
        } else if (a == "--quiet") {
            quiet = true;
        } else if (a == "--no-search") {
            noSearch = true;
        } else if (a == "--buckets") {
            edges.clear();
            const std::string v = next();
            size_t pos = 0;
            while (pos <= v.size()) {
                const size_t comma = v.find(',', pos);
                const std::string tok = v.substr(pos, (comma == std::string::npos) ? comma : comma - pos);
                if (!tok.empty()) {
                    edges.push_back((uint32_t)strtoul(tok.c_str(), nullptr, 10));
                }
                if (comma == std::string::npos) {
                    break;
                }
                pos = comma + 1;
            }
            std::sort(edges.begin(), edges.end());
        } else {
            fprintf(stderr, "unknown option %s\n", a.c_str());
            return 2;
        }
    }
    if (edges.empty()) {
        edges.push_back(0);
    }

    std::vector<Rec> recs;
    loadRecords(path, maxLines, recs);
    if (recs.empty()) {
        fprintf(stderr, "no usable records\n");
        return 1;
    }

    /* Group by read length. A record's length is its QUAL length: that is what the coder is fed,
       and it is also what the mean-length rule in the product keys off. */
    std::vector<BucketResult> results(edges.size());
    for (size_t b = 0; b < edges.size(); ++b) {
        BucketResult& r = results[b];
        r.loBytes = edges[b];
        r.hiBytes = (b + 1 < edges.size()) ? edges[b + 1] : 0;

        std::vector<Rec> sample;
        uint64_t sampleBytes = 0;
        for (const Rec& rec : recs) {
            const size_t len = rec.qual.size();
            if (len < r.loBytes) {
                continue;
            }
            if (r.hiBytes != 0 && len >= r.hiBytes) {
                continue;
            }
            ++r.reads;
            r.qualBytes += len;
            if (sample.size() >= recordsPerBucket || sampleBytes >= bytesPerBucket) {
                continue;   /* keep counting the bucket, stop growing the sample */
            }
            sampleBytes += len;
            sample.push_back(rec);
        }
        if (sample.empty()) {
            continue;
        }
        r.sampleQualBytes = sampleBytes;
        r.meanLen = (uint32_t)(sampleBytes / sample.size());

        /* One bucket stands in for one block, so the ranked alphabet is built from this bucket's
           own sample, the way the product builds it from the block it is about to encode. */
        std::vector<uint32_t> freq;
        buildRankFreq(sample, freq);

        fprintf(stderr, "\n[bucket %zu: %u..%s) reads=%llu sample=%llu bytes mean=%u]\n", b,
                r.loBytes, r.hiBytes ? std::to_string(r.hiBytes).c_str() : "inf",
                (unsigned long long)r.reads, (unsigned long long)sampleBytes, r.meanLen);

        /* Reference rows: the shipped default preset and coder_bwt_cm. */
        {
            Fcv2Cfg def;
            Score s = evaluate(sample, freq, def, timeUsWeight);
            if (s.ok) {
                r.baseBytes = s.bytes;
                r.baseUsec = s.usec;
            }
            uint32_t bwtLen = 0;
            if (encodeBwtCm(sample, bwtLen)) {
                r.bwtBytes = bwtLen;
            }
            uint32_t ppmdLen = 0;
            r.ppmdOk = encodePpmd(sample, ppmdLen, true);
            if (r.ppmdOk) {
                r.ppmdBytes = ppmdLen;
            } else {
                fprintf(stderr, "  coder_ppmd encode or round trip failed on this bucket\n");
            }
        }

        if (noSearch) {
            /*
             * Reference-only run: measure the shipped default preset and stop. Used when the point
             * is the coder comparison (fcv2 vs bwt_cm vs ppmd) rather than a trained row, where a
             * sample large enough to be meaningful would make a descent far too slow.
             */
            Fcv2Cfg def;
            const Score s = evaluate(sample, freq, def, timeUsWeight);
            if (!s.ok) {
                continue;
            }
            r.cfg = def;
            r.bytes = s.bytes;
            r.usec = s.usec;
            r.roundTripOk = true;   /* the shipped preset's round trip is covered by the test suite */
            continue;
        }

        /* Multi-start descent: every shipped preset (with both model-count regimes) plus the
           current rule's combination, so the sweep cannot be decided by one unlucky start. */
        Fcv2Cfg best;
        double bestValue = 1e30;
        bool haveBest = false;
        std::vector<Fcv2Cfg> starts;
        for (int p = 0; p < FCV2_PRESET_COUNT; ++p) {
            Fcv2Cfg c = FCV2_PRESETS[p].cfg;
            c.modelCount = FCV2_MODEL_COUNT;
            starts.push_back(c);
            c.modelCount = 4;
            starts.push_back(c);
        }
        const size_t nStarts = (startsLimit > 0 && (size_t)startsLimit < starts.size())
                               ? (size_t)startsLimit : starts.size();
        for (size_t si = 0; si < nStarts; ++si) {
            Fcv2Cfg cfg = starts[si];
            Score s = evaluate(sample, freq, cfg, timeUsWeight);
            if (!s.ok) {
                continue;
            }
            double value = s.value;
            if (!quiet) {
                fprintf(stderr, "  start %-28s %.0f\n", nameOf(cfg).c_str(), value);
            }
            for (int it = 0; it < iters; ++it) {
                if (!descend(sample, freq, cfg, value, timeUsWeight, quiet)) {
                    break;
                }
            }
            if (!haveBest || value < bestValue) {
                bestValue = value;
                best = cfg;
                haveBest = true;
            }
        }
        if (!haveBest) {
            continue;
        }
        r.cfg = best;
        Score s = evaluate(sample, freq, r.cfg, timeUsWeight);
        if (s.ok) {
            r.bytes = s.bytes;
            r.usec = s.usec;
        }
        r.roundTripOk = verifyRoundTrip(sample, freq, r.cfg);
    }

    /* Result table. */
    printf("\n=== fcv2 per read-length bucket: trained configuration ===\n");
    printf("input: %s   sample cap/bucket: %llu records / %llu bytes   iters: %d   starts: %d   "
           "time-us-weight: %.3g\n",
           path.c_str(), (unsigned long long)recordsPerBucket, (unsigned long long)bytesPerBucket,
           iters, startsLimit, timeUsWeight);
    printf("\n%-16s %8s %10s %10s %9s %8s %8s %8s %8s %6s  %s\n", "length range", "reads",
           "wholeQual", "sampleQual", "meanLen", "default", "bwt_cm", "ppmd", "trained", "rt",
           "trained configuration");
    for (size_t b = 0; b < results.size(); ++b) {
        const BucketResult& r = results[b];
        if (r.reads == 0 || r.bytes == 0) {
            printf("%-16s %8llu %10llu %10s %9s %8s %8s %8s %8s %6s  %s\n",
                   (std::to_string(r.loBytes) + ".." +
                    (r.hiBytes ? std::to_string(r.hiBytes) : "inf")).c_str(),
                   (unsigned long long)r.reads, (unsigned long long)r.qualBytes, "-", "-", "-", "-",
                   "-", "-", "-", r.reads ? "no sample" : "empty");
            continue;
        }
        char range[48];
        snprintf(range, sizeof(range), "%u..%s", r.loBytes,
                 r.hiBytes ? std::to_string(r.hiBytes).c_str() : "inf");
        const double denom = (double)(r.sampleQualBytes ? r.sampleQualBytes : 1);
        const double ratio = 100.0 * (double)r.bytes / denom;
        char defStr[32], bwtStr[32], outStr[32], ppmdStr[32];
        snprintf(defStr, sizeof(defStr), "%.3f%%", r.baseBytes ? 100.0 * r.baseBytes / denom : 0.0);
        snprintf(bwtStr, sizeof(bwtStr), "%.3f%%", r.bwtBytes ? 100.0 * r.bwtBytes / denom : 0.0);
        snprintf(outStr, sizeof(outStr), "%.3f%%", ratio);
        if (r.ppmdOk) {
            snprintf(ppmdStr, sizeof(ppmdStr), "%.3f%%", 100.0 * r.ppmdBytes / denom);
        } else {
            snprintf(ppmdStr, sizeof(ppmdStr), "FAIL");
        }
        printf("%-16s %8llu %10llu %10llu %9u %8s %8s %8s %8s %6s  %s\n", range,
               (unsigned long long)r.reads, (unsigned long long)r.qualBytes,
               (unsigned long long)r.sampleQualBytes, r.meanLen, defStr, bwtStr, ppmdStr, outStr,
               r.roundTripOk ? "ok" : "FAIL", nameOf(r.cfg).c_str());
    }

    printf("\nTrained table (paste into the per-length rule):\n");
    printf("  /* mean read length -> fcv2 configuration, trained by test/fcv2_len_train.cpp */\n");
    printf("  struct { uint32_t maxMeanLen; Fcv2Cfg cfg; } kQualFcv2ByLen[] = {\n");
    for (size_t b = 0; b < results.size(); ++b) {
        const BucketResult& r = results[b];
        if (r.bytes == 0) {
            continue;
        }
        printf("      { %9u, { %3d, %2d, %3d, %2d, %d, %s, %s, %s, %d } },   /* %s, preset=%s */\n",
               (r.hiBytes ? r.hiBytes - 1 : 0xFFFFFFFFu), r.cfg.cycleMax, r.cfg.cycleBucket,
               r.cfg.deltaMax, r.cfg.deltaBucket, r.cfg.prevShift,
               r.cfg.useDelta ? "true" : "false", r.cfg.useDedup ? "true" : "false",
               r.cfg.useQa ? "true" : "false", r.cfg.modelCount,
               (std::to_string(r.loBytes) + ".." +
                (r.hiBytes ? std::to_string(r.hiBytes) : "inf")).c_str(),
               presetNameOf(r.cfg).c_str());
    }
    printf("  };\n");

    printf("\nPer-bucket caveats: the sample is capped at %llu bytes per bucket, so a bucket whose\n"
           "real QUAL column is much larger can still favour a different configuration; re-run with\n"
           "a larger --bytes-per-bucket to confirm a row before it becomes a rule.\n",
           (unsigned long long)bytesPerBucket);
    return 0;
}
