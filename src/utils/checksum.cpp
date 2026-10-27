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
 * Streaming form of FarmHash Fingerprint64 (farmhashna::Hash64).
 * The algorithm and constants are derived from Google FarmHash, which is
 * distributed under the MIT license (see src/farmhash/COPYING).  Fingerprint64
 * itself has no incremental API, so this keeps the algorithm's rolling state
 * and at most one partial block in memory.
 */
static const uint64_t k0 = 0xc3a5c85c97cb3127ULL;
static const uint64_t k1 = 0xb492b66fbe98f273ULL;
static const uint64_t k2 = 0x9ae16a3b2f90404fULL;

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

inline uint64_t shiftMix(uint64_t v)
{
    return v ^ (v >> 47);
}

inline uint64_t hashLen16Mul(uint64_t u, uint64_t v, uint64_t mul)
{
    uint64_t a = (u ^ v) * mul;
    a ^= (a >> 47);
    uint64_t b = (v ^ a) * mul;
    b ^= (b >> 47);
    b *= mul;
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

struct FarmState {
    bool started = false;
    uint64_t x = 0;
    uint64_t y = 0;
    uint64_t z = 0;
    std::pair<uint64_t, uint64_t> v{0, 0};
    std::pair<uint64_t, uint64_t> w{0, 0};
    uint64_t totalLen = 0;
    std::vector<uint8_t> pending;
    std::vector<uint8_t> last64;

    void init(const uint8_t* firstBlock)
    {
        x = 81;
        y = 81 * k1 + 113;
        z = shiftMix(y * k2 + 113) * k2;
        v = std::make_pair(0ULL, 0ULL);
        w = std::make_pair(0ULL, 0ULL);
        x = x * k2 + fetch64(firstBlock);
        started = true;
    }

    void processBlock(const uint8_t* s)
    {
        x = rotate64(x + y + v.first + fetch64(s + 8), 37) * k1;
        y = rotate64(y + v.second + fetch64(s + 48), 42) * k1;
        x ^= w.second;
        y += v.first + fetch64(s + 40);
        z = rotate64(z + w.first, 33) * k1;
        v = weakHashLen32WithSeeds(s, v.second * k1, x + w.first);
        w = weakHashLen32WithSeeds(s + 32, z + w.second, y + fetch64(s + 16));
        std::swap(z, x);
    }

    void update(const uint8_t* data, size_t len)
    {
        if (len == 0) {
            return;
        }
        totalLen += len;

        /* Keep the last 64 bytes for the finalization step. */
        if (len >= 64) {
            last64.assign(data + len - 64, data + len);
        } else {
            const size_t excess = last64.size() + len > 64 ? last64.size() + len - 64 : 0;
            if (excess > 0) {
                last64.erase(last64.begin(), last64.begin() + (std::ptrdiff_t)excess);
            }
            last64.insert(last64.end(), data, data + len);
        }

        pending.insert(pending.end(), data, data + len);

        if (!started && pending.size() >= 65) {
            init(pending.data());
        }
        if (!started) {
            return;
        }

        size_t pos = 0;
        while (pending.size() - pos >= 65) {
            processBlock(pending.data() + pos);
            pos += 64;
        }
        if (pos > 0) {
            pending.erase(pending.begin(), pending.begin() + (std::ptrdiff_t)pos);
        }
    }

    uint64_t finish() const
    {
        if (!started) {
            /* The whole input fits in one HashLen0to64 call. */
            return util::Fingerprint64(reinterpret_cast<const char*>(pending.data()),
                                       pending.size());
        }

        const uint8_t* s = last64.data();
        const uint64_t mul = k1 + ((z & 0xff) << 1);
        std::pair<uint64_t, uint64_t> vv = v;
        std::pair<uint64_t, uint64_t> ww = w;
        ww.first += ((totalLen - 1) & 63);
        vv.first += ww.first;
        ww.first += vv.first;

        uint64_t xx = rotate64(x + y + vv.first + fetch64(s + 8), 37) * mul;
        uint64_t yy = rotate64(y + vv.second + fetch64(s + 48), 42) * mul;
        xx ^= ww.second * 9;
        yy += vv.first * 9 + fetch64(s + 40);
        uint64_t zz = rotate64(z + ww.first, 33) * mul;
        vv = weakHashLen32WithSeeds(s, vv.second * mul, xx + ww.first);
        ww = weakHashLen32WithSeeds(s + 32, zz + ww.second, yy + fetch64(s + 16));
        std::swap(zz, xx);

        return hashLen16Mul(hashLen16Mul(vv.first, ww.first, mul) +
                                shiftMix(yy) * k0 + zz,
                            hashLen16Mul(vv.second, ww.second, mul) + xx, mul);
    }
};

std::string farmHashHex(uint64_t value)
{
    char buffer[17];
    snprintf(buffer, sizeof(buffer), "%016llx", (unsigned long long)value);
    return std::string(buffer);
}

}  // namespace

const char* checksumTypeName(ChecksumType type)
{
    return (type == ChecksumType::FarmHash64) ? "farmhash64" : "md5";
}

ChecksumType checksumTypeFromName(const std::string& name)
{
    return (name == "farmhash64") ? ChecksumType::FarmHash64 : ChecksumType::MD5;
}

std::string checksumHex(ChecksumType type, const uint8_t* data, size_t len)
{
    if (type == ChecksumType::FarmHash64) {
        return farmHashHex(util::Fingerprint64(reinterpret_cast<const char*>(data), len));
    }

    uint8_t digest[16];
    char hex[33];
    md5_compute(data, len, digest);
    md5_to_hex(digest, hex);
    return std::string(hex);
}

struct ChecksumStream::Impl {
    explicit Impl(ChecksumType t) : type(t) {
        if (type == ChecksumType::MD5) {
            md5_init(&md5);
        }
    }

    ChecksumType type;
    md5_ctx md5;
    FarmState farm;
    bool finished = false;
};

ChecksumStream::ChecksumStream(ChecksumType t)
    : type(t), impl(std::make_unique<Impl>(t)) {
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
    return farmHashHex(impl->farm.finish());
}

void setChecksumMeta(Json::Value& meta, ChecksumType type, const std::string& value)
{
    meta["hashmethod"] = checksumTypeName(type);
    if (type == ChecksumType::FarmHash64) {
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
    type = checksumTypeFromName(meta.isMember("hashmethod")
                                    ? meta["hashmethod"].asString()
                                    : std::string());
    if (type == ChecksumType::FarmHash64) {
        if (meta.isMember("hash")) {
            return meta["hash"].asString();
        }
        /* Tolerate a marker with the older value key. */
        return meta.isMember("md5") ? meta["md5"].asString() : std::string();
    }
    if (meta.isMember("md5")) {
        return meta["md5"].asString();
    }
    return meta.isMember("hash") ? meta["hash"].asString() : std::string();
}
