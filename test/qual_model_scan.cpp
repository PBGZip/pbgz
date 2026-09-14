/*
 * qual_model_scan - how many mixing models should fcv2 use for QUAL?
 *
 * The QUAL column is the one field where the compressor has to choose a
 * different coder to get both speed and ratio: coder_bwt_cm has run out of
 * headroom, and fcv2 (a context mixer over several probability models of
 * different granularity) is the candidate. The mixer costs one counter update
 * per model per coded bit, so the number of models directly sets its speed, and
 * dropping the tail models (m5, m6) makes it both faster and smaller - fewer
 * weights to fit means less gradient noise.
 *
 * This tool measures that tradeoff on real QUAL data: for every model count it
 * encodes the sampled records, reports compressed size and throughput against
 * the raw QUAL volume, and decodes the stream back to prove the count is
 * self-describing (the decoder learns it from the stream header, not from a
 * compile-time constant).
 *
 * Usage:
 *   qual_model_scan <sam-file> [maxRecords] [tier]
 *
 * tier selects the context parameter preset, matching the candidates the QUAL
 * selector trials (see candidateFcv2Cfgs in qual_selector.cpp):
 *   0 default (cycle 96/16, delta 32/8, prevShift 1)
 *   1 ultra   (cycle 96/4,  delta 32/2, prevShift 0)
 *   2 coarse  (cycle 96/8,  delta 32/4, prevShift 2)
 *   3 fine    (cycle 96/24, delta 32/12, prevShift 0)
 */

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <chrono>
#include <fstream>
#include <string>
#include <vector>

#include "coder/coder_fcv2.h"

namespace {

struct Record {
    std::string qual;
    std::string seq;
    bool        rev;
};

/* One whole SAM alignment line, split into the three fields fcv2 needs. */
bool parseSamLine(const std::string& line, Record& out)
{
    if (line.empty() || line[0] == '@') {
        return false;
    }
    size_t field = 0;
    size_t start = 0;
    const char* flag = nullptr;
    size_t flagLen = 0;
    const char* seq = nullptr;
    size_t seqLen = 0;
    const char* qual = nullptr;
    size_t qualLen = 0;
    for (size_t i = 0; i <= line.size(); ++i) {
        if (i == line.size() || line[i] == '\t') {
            const char* f = line.data() + start;
            const size_t n = i - start;
            if (field == 1) { flag = f; flagLen = n; }
            else if (field == 9) { seq = f; seqLen = n; }
            else if (field == 10) { qual = f; qualLen = n; }
            ++field;
            start = i + 1;
            if (field > 10) {
                break;
            }
        }
    }
    if (flag == nullptr || seq == nullptr || qual == nullptr || qualLen == 0) {
        return false;
    }
    /* FLAG bit 0x10: SEQ/QUAL are stored reversed w.r.t. the sequencing cycle. */
    long flagValue = strtol(std::string(flag, flagLen).c_str(), nullptr, 10);
    out.rev = (flagValue & 0x10) != 0;
    /* '*' is the SAM placeholder for "not stored". */
    if (seqLen == 1 && seq[0] == '*') {
        out.seq.clear();
    } else {
        out.seq.assign(seq, seqLen);
    }
    if (qualLen == 1 && qual[0] == '*') {
        out.qual.clear();
    } else {
        out.qual.assign(qual, qualLen);
    }
    return !out.qual.empty();
}

void applyTier(int tier, Fcv2Cfg& cfg)
{
    switch (tier) {
    case 1:
        cfg.cycleBucket = 4;
        cfg.deltaBucket = 2;
        cfg.prevShift = 0;
        break;
    case 2:
        cfg.cycleBucket = 8;
        cfg.deltaBucket = 4;
        cfg.prevShift = 2;
        break;
    case 3:
        cfg.cycleBucket = 24;
        cfg.deltaBucket = 12;
        cfg.prevShift = 0;
        break;
    default:
        break;
    }
}

const char* tierName(int tier)
{
    switch (tier) {
    case 1: return "ultra(96/4,32/2,ps0)";
    case 2: return "coarse(96/8,32/4,ps2)";
    case 3: return "fine(96/24,32/12,ps0)";
    default: return "default(96/16,32/8,ps1)";
    }
}

/*
 * Encode then decode the whole record set with `cfg`, and report both times.
 * The decode side gets the same base sequence it would have after decoding SEQ,
 * which is what the compression path does.
 */
bool runOne(const std::vector<Record>& records, const std::vector<uint32_t>& freq,
            const Fcv2Cfg& cfg, size_t rawBytes,
            size_t& outLen, double& encMBps, double& decMBps, bool& roundTripOk)
{
    std::vector<uint8_t> buf(rawBytes * 2 + (1u << 16), 0);
    int32_t written = 0;
    {
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_fcv2 coder(&io, freq, cfg);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < records.size(); ++i) {
            const Record& r = records[i];
            coder.encode_record((const uint8_t*)r.qual.data(), (uint32_t)r.qual.size(), r.rev,
                                (const uint8_t*)r.seq.data(), (uint32_t)r.seq.size());
        }
        written = coder.encode_flush();
        const auto t1 = std::chrono::steady_clock::now();
        const double sec = std::chrono::duration<double>(t1 - t0).count();
        encMBps = (sec > 0.0) ? (double)rawBytes / 1e6 / sec : 0.0;
    }
    if (written <= 0) {
        return false;
    }
    outLen = (size_t)written;

    coder_io decIo(buf.data(), written);
    coder_fcv2 decoder(&decIo, freq);
    if (decoder.begin_decode() != 0) {
        return false;
    }
    std::vector<uint8_t> got;
    roundTripOk = true;
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < records.size() && roundTripOk; ++i) {
        const Record& r = records[i];
        got.assign(r.qual.size(), 0);
        const int32_t n = decoder.decode_record(got.data(), (uint32_t)got.size(),
                                                (const uint8_t*)r.seq.data(),
                                                (uint32_t)r.seq.size());
        roundTripOk = (n == (int32_t)r.qual.size() && memcmp(got.data(), r.qual.data(), got.size()) == 0);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    decMBps = (sec > 0.0) ? (double)rawBytes / 1e6 / sec : 0.0;
    return true;
}

}  /* namespace */

int main(int argc, char** argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s <sam-file> [maxRecords] [tier 0..3]\n", argv[0]);
        return 2;
    }
    const char* path = argv[1];
    size_t maxRecords = (argc > 2) ? (size_t)strtoull(argv[2], nullptr, 10) : 200000;
    int tier = (argc > 3) ? atoi(argv[3]) : 0;

    std::ifstream in(path);
    if (!in) {
        fprintf(stderr, "cannot open %s\n", path);
        return 2;
    }

    std::vector<Record> records;
    records.reserve(maxRecords);
    std::vector<uint32_t> freq(256, 0);
    std::string line;
    size_t rawBytes = 0;
    while (records.size() < maxRecords && std::getline(in, line)) {
        Record r;
        if (!parseSamLine(line, r)) {
            continue;
        }
        for (size_t i = 0; i < r.qual.size(); ++i) {
            freq[(uint8_t)r.qual[i]]++;
        }
        rawBytes += r.qual.size();
        records.push_back(std::move(r));
    }
    if (records.empty()) {
        fprintf(stderr, "no usable records in %s\n", path);
        return 1;
    }

    Fcv2Cfg base;
    applyTier(tier, base);

    printf("file        : %s\n", path);
    printf("records     : %zu   raw QUAL bytes: %zu (%.1f MB)\n",
           records.size(), rawBytes, (double)rawBytes / 1e6);
    printf("tier        : %s\n\n", tierName(tier));

    struct Row {
        int    n;
        size_t bytes;
        double encMBps;
        double decMBps;
        bool   rt;
        bool   ok;
    };
    std::vector<Row> rows;
    size_t refLen = 0;
    for (int n = 1; n <= FCV2_MAX_MODEL_COUNT; ++n) {
        Fcv2Cfg cfg = base;
        cfg.modelCount = n;
        Row row;
        row.n = n;
        row.bytes = 0;
        row.encMBps = 0.0;
        row.decMBps = 0.0;
        row.rt = false;
        row.ok = runOne(records, freq, cfg, rawBytes, row.bytes, row.encMBps, row.decMBps, row.rt);
        rows.push_back(row);
        if (n == FCV2_MAX_MODEL_COUNT && row.ok) {
            refLen = row.bytes;
        }
    }

    printf("%3s %10s %8s %9s %9s %9s %4s\n",
           "n", "bytes", "ratio%", "encMB/s", "decMB/s", "vs n=7", "rt");
    for (size_t i = 0; i < rows.size(); ++i) {
        const Row& r = rows[i];
        if (!r.ok) {
            printf("%3d %10s %8s %9s %9s %9s %4s\n", r.n, "-", "-", "-", "-", "-", "ERR");
            continue;
        }
        const double ratio = 100.0 * (double)r.bytes / (double)rawBytes;
        char delta[32] = "-";
        if (refLen != 0) {
            snprintf(delta, sizeof(delta), "%+.2f%%",
                     100.0 * ((double)r.bytes - (double)refLen) / (double)refLen);
        }
        printf("%3d %10zu %7.2f%% %9.3f %9.3f %9s %4s\n",
               r.n, r.bytes, ratio, r.encMBps, r.decMBps, delta, r.rt ? "OK" : "BAD");
    }
    printf("\nratio%% is of the raw QUAL volume; 'vs n=7' is the size change against\n"
           "the full mix, so negative means the trimmed mix compresses smaller.\n");
    return 0;
}
