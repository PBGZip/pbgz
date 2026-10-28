/*
 * checksum.cpp - checksum algorithms and metadata markers for pbgz
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "utils/checksum.h"

#include <cstddef>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

#include "farmhash/src/farmhash.h"
#include "utils/md5_util.h"

namespace {

/*
 * Streaming form of FarmHash Fingerprint128 (CityHash128).  The algorithm and
 * constants are derived from Google FarmHash, which is distributed under the
 * MIT license (see src/farmhash/COPYING).  Fingerprint128 has no incremental
 * API, so this keeps the rolling 128-bit state and at most one partial block.
 */
static const uint64_t k0 = 0xc3a5c85c97cb3127ULL;
static const uint64_t k1 = 0xb492b66fbe98f273ULL;

inline uint64_t fetch64(const uint8_t* p)
{
    uint64_t v = 0;
    memcpy(&v, p, sizeof(v));
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    v = __builtin_bswap64(v);
#endif
    return v;
}

inline uint64_t rotate64(uint64_t v, int shift)
{
    return (v >> shift) | (v << (64 - shift));
}

/* Hash128to64() from FarmHash: Murmur-inspired 128-bit to 64-bit mixing. */
inline uint64_t hashLen16(uint64_t u, uint64_t v)
{
    const uint64_t kMul = 0x9ddfea08eb382d69ULL;
    uint64_t a = (u ^ v) * kMul;
    a ^= (a >> 47);
    uint64_t b = (v ^ a) * kMul;
    b ^= (b >> 47);
    b *= kMul;
    return b;
}

inline std::pair<uint64_t, uint64_t> weakHashLen32WithSeeds(
    uint64_t w, uint64_t x, uint64_t y, uint64_t z, uint64_t a, uint64_t b)
{
    a += w;
    b = rotate64(b + a + z, 21);
    uint64_t c = a;
    a += x;
    a += y;
    b += rotate64(a, 44);
    return std::make_pair(a + z, b + c);
}

inline std::pair<uint64_t, uint64_t> weakHashLen32WithSeeds(
    const uint8_t* s, uint64_t a, uint64_t b)
{
    return weakHashLen32WithSeeds(fetch64(s), fetch64(s + 8),
                                  fetch64(s + 16), fetch64(s + 24), a, b);
}

struct Farm128State {
    bool mainLoop = false;   /* body length >= 128: CityHash128's block path */
    bool started = false;
    bool finished = false;
    uint64_t totalLen = 0;   /* length declared at construction */
    uint64_t accepted = 0;   /* bytes accepted so far (never exceeds totalLen) */
    uint64_t bodyLen = 0;    /* bytes hashed by the block loop plus the tail */
    uint64_t bodyDone = 0;   /* bytes already consumed by the block loop */
    uint64_t x = 0;
    uint64_t y = 0;
    uint64_t z = 0;
    std::pair<uint64_t, uint64_t> v{0, 0};
    std::pair<uint64_t, uint64_t> w{0, 0};
    /* Input before the state is initialised, then the trailing partial block. */
    std::vector<uint8_t> pending;
    /* First 16 bytes, consumed as the seed by begin(). */
    uint8_t seed[16] = {0};
    bool hasSeed = false;
    /* Last 32 processed bytes: the tail loop may read back into them. */
    std::vector<uint8_t> history;

    void init(uint64_t length)
    {
        totalLen = length;
        bodyLen = (totalLen >= 16) ? totalLen - 16 : totalLen;
        mainLoop = bodyLen >= 128;
    }

    /* CityHash128 consumes the first 16 bytes as the seed. */
    void begin()
    {
        const uint8_t* head = pending.data();
        const uint8_t* body = head + 16;
        memcpy(seed, head, sizeof(seed));
        hasSeed = true;
        x = fetch64(head);
        y = fetch64(head + 8) + k0;
        z = bodyLen * k1;
        v.first = rotate64(y ^ k1, 49) * k1 + fetch64(body);
        v.second = rotate64(v.first, 42) * k1 + fetch64(body + 8);
        w.first = rotate64(y + z, 35) * k1 + x;
        w.second = rotate64(x + fetch64(body + 88), 53) * k1;
        pending.erase(pending.begin(), pending.begin() + 16);
        started = true;
    }

    /* The two 64-byte halves of one CityHash128 block, as in the vendored code. */
    void processBlock(const uint8_t* s)
    {
        for (int half = 0; half < 2; ++half) {
            x = rotate64(x + y + v.first + fetch64(s + 8), 37) * k1;
            y = rotate64(y + v.second + fetch64(s + 48), 42) * k1;
            x ^= w.second;
            y += v.first + fetch64(s + 40);
            z = rotate64(z + w.first, 33) * k1;
            v = weakHashLen32WithSeeds(s, v.second * k1, x + w.first);
            w = weakHashLen32WithSeeds(s + 32, z + w.second, y + fetch64(s + 16));
            std::swap(z, x);
            s += 64;
        }
    }

    void consume()
    {
        size_t pos = 0;
        while (bodyLen - bodyDone >= 128 && pending.size() - pos >= 128) {
            processBlock(pending.data() + pos);
            pos += 128;
            bodyDone += 128;
            history.assign(pending.data() + pos - 32, pending.data() + pos);
        }
        if (pos > 0) {
            pending.erase(pending.begin(), pending.begin() + (std::ptrdiff_t)pos);
        }
    }

    void update(const uint8_t* data, size_t len)
    {
        if (finished || len == 0 || accepted >= totalLen) {
            return;
        }
        if ((uint64_t)len > totalLen - accepted) {
            len = (size_t)(totalLen - accepted);
        }
        accepted += len;

        if (!mainLoop) {
            /* The whole input is shorter than 144 bytes; keep it for finish(). */
            pending.insert(pending.end(), data, data + len);
            return;
        }

        pending.insert(pending.end(), data, data + len);
        if (!started) {
            /* Seed (16) plus the 96 bytes the state initialisation reads. */
            if (pending.size() < 112) {
                return;
            }
            begin();
        }
        consume();
    }

    std::pair<uint64_t, uint64_t> finish()
    {
        finished = true;
        if (bodyDone == 0) {
            /* Fewer than 128 body bytes were consumed (short input, or the
               caller declared more than it supplied): hash what actually
               arrived instead of running the block/tail path. */
            std::vector<uint8_t> actual;
            if (hasSeed) {
                actual.insert(actual.end(), seed, seed + sizeof(seed));
            }
            actual.insert(actual.end(), pending.begin(), pending.end());
            auto value = util::Fingerprint128(reinterpret_cast<const char*>(actual.data()),
                                              actual.size());
            return std::make_pair(util::Uint128Low64(value), util::Uint128High64(value));
        }

        /* The tail loop may read up to 31 bytes before the trailing block. */
        uint64_t tailLen = pending.size();
        x += rotate64(v.first + z, 49) * k0;
        y = y * k0 + rotate64(w.second, 37);
        z = z * k0 + rotate64(w.first, 27);
        w.first *= 9;
        v.first *= k0;

        std::vector<uint8_t> tail;
        tail.reserve(history.size() + tailLen);
        tail.insert(tail.end(), history.begin(), history.end());
        tail.insert(tail.end(), pending.begin(), pending.end());
        const uint8_t* s = tail.data() + history.size();
        for (size_t tailDone = 0; tailDone < tailLen; ) {
            tailDone += 32;
            y = rotate64(x + y, 42) * k0 + v.second;
            w.first += fetch64(s + tailLen - tailDone + 16);
            x = x * k0 + w.first;
            z += w.second + fetch64(s + tailLen - tailDone);
            w.second += v.first;
            v = weakHashLen32WithSeeds(s + tailLen - tailDone, v.first + z, v.second);
            v.first *= k0;
        }

        x = hashLen16(x, v.first);
        y = hashLen16(y + z, w.first);
        return std::make_pair(hashLen16(x + v.second, w.second) + y,
                              hashLen16(x + w.second, y + v.second));
    }
};

std::string farm128Hex(uint64_t low, uint64_t high)
{
    char buffer[33];
    snprintf(buffer, sizeof(buffer), "%016llx%016llx",
             (unsigned long long)high, (unsigned long long)low);
    return std::string(buffer);
}

}  // namespace

const char* checksumTypeName(ChecksumType type)
{
    return (type == ChecksumType::FarmHash128) ? "farmhash128" : "md5";
}

const char* checksumTypeTag(ChecksumType type)
{
    return (type == ChecksumType::FarmHash128) ? "f" : "m";
}

ChecksumType checksumTypeFromName(const std::string& name)
{
    return (name == "f") ? ChecksumType::FarmHash128 : ChecksumType::MD5;
}

std::string checksumHex(ChecksumType type, const uint8_t* data, size_t len)
{
    if (type == ChecksumType::FarmHash128) {
        auto value = util::Fingerprint128(reinterpret_cast<const char*>(data), len);
        return farm128Hex(util::Uint128Low64(value), util::Uint128High64(value));
    }

    uint8_t digest[16];
    char hex[33];
    md5_compute(data, len, digest);
    md5_to_hex(digest, hex);
    return std::string(hex);
}

struct ChecksumStream::Impl {
    Impl(ChecksumType t, uint64_t totalLen) : type(t) {
        if (type == ChecksumType::MD5) {
            md5_init(&md5);
        } else {
            farm.init(totalLen);
        }
    }

    ChecksumType type;
    md5_ctx md5;
    Farm128State farm;
    bool finished = false;
};

ChecksumStream::ChecksumStream(ChecksumType t, uint64_t totalLen)
    : type(t), impl(std::make_unique<Impl>(t, totalLen)) {
}

ChecksumStream::~ChecksumStream() = default;

void ChecksumStream::update(const uint8_t* data, size_t len)
{
    if (impl->finished || data == nullptr || len == 0) {
        return;
    }
    if (type == ChecksumType::MD5) {
        md5_update(&impl->md5, data, len);
    } else {
        impl->farm.update(data, len);
    }
}

std::string ChecksumStream::finish()
{
    if (impl->finished) {
        return std::string();
    }
    impl->finished = true;
    if (type == ChecksumType::MD5) {
        uint8_t digest[16];
        char hex[33];
        md5_final(digest, &impl->md5);
        md5_to_hex(digest, hex);
        return std::string(hex);
    }
    std::pair<uint64_t, uint64_t> value = impl->farm.finish();
    return farm128Hex(value.first, value.second);
}

void setChecksumMeta(Json::Value& meta, ChecksumType type, const std::string& value)
{
    meta[kChecksumMethodKey] = checksumTypeTag(type);
    if (type == ChecksumType::FarmHash128) {
        meta["hash"] = value;
        if (meta.isMember("md5")) {
            meta.removeMember("md5");
        }
    } else {
        meta["md5"] = value;
        if (meta.isMember("hash")) {
            meta.removeMember("hash");
        }
    }
}

std::string getChecksumMeta(const Json::Value& meta, ChecksumType& type)
{
    type = checksumTypeFromName(meta.isMember(kChecksumMethodKey)
                                    ? meta[kChecksumMethodKey].asString()
                                    : std::string());
    const char* key = (type == ChecksumType::FarmHash128) ? "hash" : "md5";
    return meta.isMember(key) ? meta[key].asString() : std::string();
}
