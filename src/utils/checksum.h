/*
 * checksum.h - checksum algorithms and metadata markers for pbgz
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

#pragma once

#include <stddef.h>
#include <stdint.h>

#include <memory>
#include <string>

#include <json/json.h>

/*
 * Checksum algorithm recorded in metadata.  Files written before the marker
 * existed (or files marked "md5") are verified with MD5; new files use
 * FarmHash's portable Fingerprint64.
 */
enum class ChecksumType {
    MD5 = 0,
    FarmHash64 = 1,
};

const char* checksumTypeName(ChecksumType type);

/* Unknown/absent names fall back to MD5 for backward compatibility. */
ChecksumType checksumTypeFromName(const std::string& name);

/* Hex string: 32 chars for MD5, 16 chars for FarmHash64. */
std::string checksumHex(ChecksumType type, const uint8_t* data, size_t len);

inline std::string checksumHex(ChecksumType type, const std::string& data) {
    return checksumHex(type, reinterpret_cast<const uint8_t*>(data.data()), data.size());
}

/*
 * Incremental checksum for streaming inputs.  MD5 uses the existing
 * incremental implementation; FarmHash64 is a streaming form of
 * Fingerprint64 (farmhashna::Hash64) so large reference files never need to
 * be resident in memory.
 */
class ChecksumStream {
public:
    explicit ChecksumStream(ChecksumType type);
    ~ChecksumStream();

    void update(const uint8_t* data, size_t len);

    std::string finish();

private:
    ChecksumType type;
    struct Impl;
    std::unique_ptr<Impl> impl;
};

/*
 * Block metadata: new files carry {"hashmethod":"farmhash64","hash":"..."};
 * old files carry only {"md5":"..."} or {"hashmethod":"md5","md5":"..."}.
 */
void setChecksumMeta(Json::Value& meta, ChecksumType type, const std::string& value);

/* Returns the expected value; type is set to the algorithm selected by the marker. */
std::string getChecksumMeta(const Json::Value& meta, ChecksumType& type);
