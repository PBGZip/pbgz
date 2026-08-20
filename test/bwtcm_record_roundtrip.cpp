/*
 * Per-record round-trip verification for coder_bwt_cm (used to pinpoint the root cause of
 * decompression failures when QUAL goes through bwt_cm).
 *
 * Question to answer: if QUAL records are fed one at a time into bwt_cm via encode_line(ptr, len)
 * and read back with decode_line(out, len) at the same lengths, can the bytes be restored exactly?
 *
 * The coder itself must be unit-tested first because an end-to-end failure (decoding 64 MB vs an
 * original 19 MB) has two possible causes: the coder's fixed-length stream interface is broken on
 * its own, or the wiring in sam_actuator is wrong. The two fixes are completely different, and
 * reading the code cannot tell them apart, so a minimal case with no pbgz context is needed to
 * pin the coder down in isolation.
 *
 * Key point: when decode_line is called with split_ch = UINT8_MAX, it takes the pure fixed-length
 * stream branch: it ignores separators and only takes as many bytes as the caller's out_len. The
 * encode side feeding data with no separator is its counterpart, so it should hold in theory; this
 * test verifies that "in theory".
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "coder/coder.h"
#include "coder/coder_bwt_cm.h"

namespace {

void registerCoderCallbacks()
{
    coder_ns::register_alloc_proc([](size_t size) -> uint8_t* {
        return (uint8_t*)malloc(size);
    });
    coder_ns::register_realloc_proc(
        [](size_t& size, uint8_t* ptr, size_t newSize) -> uint8_t* {
            uint8_t* p = (uint8_t*)realloc(ptr, newSize);
            if (p != NULL) size = newSize;
            return p;
        });
    coder_ns::register_free_func([](void*& ptr) { free(ptr); ptr = NULL; });
    coder_ns::resister_logger_proc([](int, const char* msg) {
        fprintf(stderr, "coder log: %s\n", msg != NULL ? msg : "");
    });
    coder_ns::initFcCoder();
}

/*
 * Build a batch of pseudo quality-value records of varying lengths. The lengths are deliberately
 * unequal because QUAL length varies with SEQ in real SAM data; fixed-length records would hide
 * problems such as a record being split across a block boundary.
 */
struct Records {
    std::vector<uint8_t> flat;
    std::vector<uint32_t> lens;
};

Records makeRecords(uint32_t count, uint32_t baseLen)
{
    Records r;
    uint32_t seed = 12345;
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t len = baseLen + (i % 37);
        r.lens.push_back(len);
        for (uint32_t j = 0; j < len; ++j) {
            seed = seed * 1103515245u + 12345u;
            /* quality values are concentrated on a few symbols, close to the real distribution; otherwise there is nothing to compress and problems would not show */
            static const char alphabet[] = "!,:FI";
            r.flat.push_back((uint8_t)alphabet[(seed >> 16) % 5]);
        }
    }
    return r;
}

/* Return 0 to indicate success. */
int runCase(const char* name, uint32_t count, uint32_t baseLen, int level)
{
    Records rec = makeRecords(count, baseLen);
    const uint32_t srcLen = (uint32_t)rec.flat.size();

    /* Give the compression buffer plenty of room so overflow does not get misjudged as a coder defect. */
    std::vector<uint8_t> comp(srcLen * 2 + (1u << 20), 0);

    uint32_t dstLen = 0;
    {
        coder_io io(comp.data(), (uint32_t)comp.size());
        io.meta["level"] = (Json::Value::Int)level;
        coder_bwt_cm enc(&io);
        uint32_t off = 0;
        for (uint32_t i = 0; i < rec.lens.size(); ++i) {
            enc.encode_line(rec.flat.data() + off, rec.lens[i]);
            off += rec.lens[i];
        }
        enc.encode_flush();
        dstLen = io.data_len;
    }

    std::vector<uint8_t> out(srcLen + 1024, 0);
    uint32_t outOff = 0;
    int bad = 0;
    {
        coder_io io(comp.data(), dstLen);
        coder_bwt_cm dec(&io);
        for (uint32_t i = 0; i < rec.lens.size(); ++i) {
            int32_t got = dec.decode_line(out.data() + outOff, rec.lens[i]);
            if (got != (int32_t)rec.lens[i]) {
                fprintf(stderr,
                        "  [%s] record %u length mismatch: expected %u, got %d\n",
                        name, i, rec.lens[i], got);
                bad = 1;
                break;
            }
            outOff += rec.lens[i];
        }
    }

    if (!bad && outOff != srcLen) {
        fprintf(stderr, "  [%s] total length mismatch: expected %u, got %u\n", name, srcLen, outOff);
        bad = 1;
    }
    if (!bad && memcmp(out.data(), rec.flat.data(), srcLen) != 0) {
        for (uint32_t i = 0; i < srcLen; ++i) {
            if (out[i] != rec.flat[i]) {
                fprintf(stderr, "  [%s] first mismatch at offset %u: expected 0x%02X, got 0x%02X\n",
                        name, i, rec.flat[i], out[i]);
                break;
            }
        }
        bad = 1;
    }

    printf("%-28s level=%d  %u records / %u bytes -> %u bytes (%.2f%%)  %s\n",
           name, level, count, srcLen, dstLen,
           srcLen == 0 ? 0.0 : (double)dstLen * 100.0 / (double)srcLen,
           bad ? "FAIL" : "PASS");
    return bad;
}

} // namespace

int main()
{
    registerCoderCallbacks();

    int fail = 0;
    /* Small data: everything falls on encode_flush's single-block path. */
    fail |= runCase("single-block-small-data", 2000, 100, 4);
    /* Large data: must span bsize, triggering whole-block commits in encode_line and records being split. */
    fail |= runCase("cross-block-level1(1MB-block)", 40000, 150, 1);
    /* Multiple cross-block spans. */
    fail |= runCase("multi-cross-block-level1", 200000, 150, 1);
    /* Production default level. */
    fail |= runCase("default-tier-level4", 50000, 150, 4);

    printf("\n%s\n", fail ? "result: some cases failed" : "result: all passed");
    return fail;
}
