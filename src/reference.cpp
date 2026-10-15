/*
 * reference.h - Header file for pbgz project
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

#include <algorithm>   /* std::sort */
#include <atomic>
#include <utility>

#include "spinlock/spinlock-pthread.h"
#include "city.h"
#include "utils/md5_util.h"
#include "reference.h"
#include "pbgz_errno.h"
#include "log/logger.h"
#include "io_wrapper.h"
#include "pbgz_file.h"
#include "utils/memory_util.h"
#include "coder_json.h"
#include "utils/path_util.h"
#include "pbgz_manager.h"
#include "actg.h"
#ifdef __SSE4_2__
#include "hardware.h"
#endif

Reference::Reference(const std::string& fastaName, uint32_t threadNum) {
    // Initialize pointers to nullptr
    refGeneSquash = nullptr;
    refGeneSquashlen = 0;
    hashBucketCnt = nullptr;
    hashTableBuffer = nullptr;
    refGeneSquashMatched = nullptr;
    refGeneSquashMatchedlen = 0;

    // Set reference file path and thread count
    refGeneFile = fastaName;
    parallel = threadNum;
    guardBar = nullptr;
}

int32_t Reference::referencCheck() {
    if ((baseGroupStep & 0x3) != 0) {
        LOG_ERROR("Invalid base group step %d.", baseGroupStep);
        return pbgz::PBGZ_ERR_INTERNEL;
    }

    if (baseGroupLen != 31) {
        LOG_ERROR("Invalid base group length %d.", baseGroupLen);
        return pbgz::PBGZ_ERR_INTERNEL;
    }

    if (baseGroupLen > baseGroupStep) {
        LOG_ERROR("Reference base group length %d must be less than base group step.", baseGroupLen, baseGroupStep);
        return pbgz::PBGZ_ERR_INTERNEL;
    }

    return pbgz::PBGZ_ERR_OK;
}

Reference::~Reference() {
    // Free allocated memory
    MemoryUtil::safeFree(refGeneSquash);
    MemoryUtil::safeFree(hashBucketCnt);
    MemoryUtil::safeFree(hashTableBuffer);
    MemoryUtil::safeFree(refGeneSquashMatched);
    MemoryUtil::safeDeleteClass(guardBar);
}

/* Read the whole index, verify its checksum, and detach from disk. */
bool Reference::loadFromNiFile(const std::string& niFile) {
    const int64_t magicLen = (int64_t)PBGZ_FILE_MAGIC.length();
    const int64_t headerLen = magicLen + 3 + 2 + (int64_t)sizeof(int64_t);
    uint8_t* image = nullptr;

    /* An unusable index is not fatal: warn and fall back to the FASTA. */
    auto reject = [&](const std::string& why) {
        MemoryUtil::safeFree(image);
        fprintf(stderr, "warning: index %s is unusable (%s), building the index from %s instead.\n",
                niFile.c_str(), why.c_str(), refGeneFile.c_str());
        LOG_WARNING("ni file %s rejected: %s", niFile.c_str(), why.c_str());
        return false;
    };

    FileReader niReader(niFile);
    if (0 != niReader.openIO()) {
        return reject("cannot be opened");
    }

    const int64_t fileSize = niReader.getFileSize();
    if (fileSize <= headerLen) {
        niReader.closeIO();
        return reject("too short to be an index");
    }

    image = MemoryUtil::safeAlloc<uint8_t>(fileSize);
    if (image == nullptr) {
        niReader.closeIO();
        return reject("cannot be loaded into memory");
    }

    const int64_t readLen = niReader.readIO(image, fileSize);
    niReader.closeIO();
    if (readLen != fileSize) {
        return reject("is being written or was truncated while reading");
    }

    if (0 != memcmp(image, PBGZ_FILE_MAGIC.c_str(), magicLen) ||
        0 != memcmp(image + magicLen + 3, "ni", 2)) {
        return reject("is not a pbgz index");
    }

    int64_t dataLen = 0;
    memcpy(&dataLen, image + magicLen + 3 + 2, sizeof(dataLen));
    if (dataLen <= 0 || dataLen >= fileSize - headerLen) {
        return reject("declares a payload length that does not fit the file");
    }
    const int64_t metaLen = fileSize - headerLen - dataLen;

    Json::Value niMeta;
    try {
        coder_json cmCoder;
        cmCoder.decoder(image + headerLen + dataLen, metaLen, niMeta);
    } catch (const std::exception& e) {
        return reject(std::string("carries unreadable metadata: ") + e.what());
    }

    /* The index must have been built from the reference selected by -r. */
    const std::string wantName = PathUtil::getFileName(refGeneFile);
    const std::string haveName = niMeta["refe_name"].asString();
    if (haveName != wantName) {
        return reject("was built from '" + haveName + "', not from '" + wantName + "'");
    }

    const std::string expectMd5 = niMeta["ni_data_md5"].asString();
    if (expectMd5.empty()) {
        return reject("carries no payload checksum");
    }

    md5_ctx dataMd5;
    uint8_t digest[16];
    md5_init(&dataMd5);
    md5_update(&dataMd5, image + headerLen, dataLen);
    md5_final(digest, &dataMd5);
    if (expectMd5 != std::string(md5_to_string(digest).c_str())) {
        return reject("failed its own checksum, it was modified or written concurrently");
    }

    memmove(image, image + headerLen, dataLen);
    refGeneSquash = image;
    refGeneSquashlen = dataLen;
    fastaChecksum = niMeta["refe_orgfile_md5"].asString();
    fastaLength = niMeta["refe_len"].asInt64();

    fprintf(stderr, "Using index %s for reference %s (verified in memory).\n",
            niFile.c_str(), refGeneFile.c_str());
    return true;
}

bool Reference::initSquashFromFasta() {
    std::string refGeneMd5;
    int64_t refGeneLen;
    calculateReferenceMd5(refGeneFile, refGeneMd5, refGeneLen);
    fastaChecksum = refGeneMd5;
    fastaLength = refGeneLen;

    std::unique_ptr<IOReader> reader = createIOReader(refGeneFile);
    if (reader == nullptr || reader->openIO() != 0) {
        LOG_ERROR("Failed to open reference FASTA: %s", refGeneFile.c_str());
        return false;
    }

    /* Build squash in memory: 4 bases per byte (2 bits each). */
    std::vector<uint8_t> squash;
    squash.reserve(refGeneLen / 4 + 1024);

    uint8_t cacheActg[4];
    uint32_t cacheLen = 0;
    std::string line;

    /* Base offset at which the next record starts. */
    uint64_t baseCount = 0;
    refSequences.clear();

    while (reader->readLine(line) > 0) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            /* Record name: first whitespace-delimited token of the header. */
            const size_t nameEnd = line.find_first_of(" \t", 1);
            refSequences.emplace_back(line.substr(1, nameEnd - 1), baseCount);
            continue;
        }

        uint32_t lineLen = (uint32_t)line.length();
        uint32_t pos = 0;
        baseCount += lineLen;

        if (cacheLen > 0) {
            uint32_t need = 4 - cacheLen;
            uint32_t take = (lineLen < need) ? lineLen : need;
            memcpy(cacheActg + cacheLen, line.c_str(), take);
            cacheLen += take;
            pos += take;
            if (cacheLen == 4) {
                uint8_t ch;
                actgSquash(cacheActg, 4, &ch);
                squash.push_back(ch);
                cacheLen = 0;
            }
        }

        uint32_t remaining = lineLen - pos;
        uint32_t aligned = remaining >> 2 << 2;
        if (aligned > 0) {
            /* Convert the aligned part straight into the squash tail. */
            const size_t at = squash.size();
            squash.resize(at + (aligned >> 2));
            actgSquash((const uint8_t*)(line.c_str() + pos), aligned, squash.data() + at);
            pos += aligned;
        }

        uint32_t tail = remaining - aligned;
        if (tail > 0) {
            memcpy(cacheActg, line.c_str() + pos, tail);
            cacheLen = tail;
        }
    }

    if (cacheLen > 0) {
        uint8_t ch = 0;
        for (uint32_t n = 0; n < cacheLen; n++) {
            ch |= ((cacheActg[n] & 0x06) >> 1) << (6 - (n << 1));
        }
        squash.push_back(ch);
    }

    refGeneSquashlen = (int64_t)squash.size();
    refGeneSquash = MemoryUtil::safeAlloc<uint8_t>(refGeneSquashlen);
    if (refGeneSquash == nullptr) {
        return false;
    }
    memcpy(refGeneSquash, squash.data(), refGeneSquashlen);
    reader->closeIO();
    return true;
}

uint8_t* Reference::initSquashByStream(int64_t squashLength) {
    refGeneSquashlen = squashLength;
    refGeneSquash = MemoryUtil::safeAlloc<uint8_t>(refGeneSquashlen);
    return refGeneSquash;
}

void Reference::calculateReferenceMd5(const std::string& refGeneFile, std::string& refGeneMd5, int64_t& refGeneLen) {
    std::thread calc_md5([&refGeneFile, &refGeneMd5, &refGeneLen]() {
        uint32_t bufferLen = 4096, len;
        refGeneLen = 0;
        md5_ctx md5;
        uint8_t digest[16];
        uint8_t *buffer = MemoryUtil::safeAlloc<uint8_t>(bufferLen);
        md5_init(&md5);

        // Detect gzip input
        bool isGzFile = PathUtil::isGzFile(refGeneFile);
        // Choose a reader based on file type
        if (isGzFile) {
            // Choose gzip reader based on SIMD support
            bool useFastGzReader = false;
#ifdef __SSE4_2__
            useFastGzReader = Hardware::isSupportSimd();
#endif
            if (useFastGzReader) {
                // Intel ISA-L accelerated reader
                FastGzFileReader fastGzReader(refGeneFile);
                if (fastGzReader.openIO() == 0) {
                    for (;;) {
                        len = fastGzReader.readIO(buffer, bufferLen);
                        if (len == 0 && bufferLen > 0) { // Check error return value - 0 means read failed or file ended
                            LOG_ERROR("Failed to read gzip file with FastGzFileReader: %s", refGeneFile.c_str());
                            refGeneMd5 = "";
                            refGeneLen = 0;
                            MemoryUtil::safeFree(buffer);
                            return;
                        }
                        refGeneLen += len;
                        md5_update(&md5, buffer, len);
                        if (len != bufferLen) {
                            break;
                        }
                    }
                    fastGzReader.closeIO();
                } else {
                    LOG_ERROR("Failed to open gzip file with FastGzFileReader: %s", refGeneFile.c_str());
                    refGeneMd5 = "";
                    refGeneLen = 0;
                    MemoryUtil::safeFree(buffer);
                    return;
                }
            } else {
                // Standard gzip reader
                GzFileReader gzReader(refGeneFile, 1); // Use single thread reading
                if (gzReader.openIO() == 0) {
                    for (;;) {
                        len = gzReader.readIO(buffer, bufferLen);
                        refGeneLen += len;
                        md5_update(&md5, buffer, len);
                        if (len != bufferLen) {
                            break;
                        }
                    }
                    gzReader.closeIO();
                } else {
                    LOG_ERROR("Failed to open gzip file: %s", refGeneFile.c_str());
                    refGeneMd5 = "";
                    refGeneLen = 0;
                    MemoryUtil::safeFree(buffer);
                    return;
                }
            }
        } else {
            // Plain text reader
            FileReader fileReader(refGeneFile);
            if (fileReader.openIO() == 0) {
                for (;;) {
                    len = fileReader.readIO(buffer, bufferLen);
                    refGeneLen += len;
                    md5_update(&md5, buffer, len);
                    if (len != bufferLen) {
                        break;
                    }
                }
                fileReader.closeIO();
            } else {
                LOG_ERROR("Failed to open file: %s", refGeneFile.c_str());
                refGeneMd5 = "";
                refGeneLen = 0;
                MemoryUtil::safeFree(buffer);
                return;
            }
        }

        md5_final(digest, &md5);
        refGeneMd5 = md5_to_string(digest);
        MemoryUtil::safeFree(buffer);
    });
    calc_md5.join();
}

bool Reference::writeNiFileHeader(FileWriter& niWriter, int64_t& offset) {
    niWriter.writeIO(PBGZ_FILE_MAGIC.c_str(), PBGZ_FILE_MAGIC.length());
    niWriter.writeIO(PbgzManager::getInstance().getVersionAsArray().data(), 3);
    niWriter.writeIO((uint8_t *)"ni", 2);
    offset = PBGZ_FILE_MAGIC.length() + 3 + 2;
    niWriter.writeIO((uint8_t *)(&offset), sizeof(offset));
    return true;
}


// Create an IOReader based on file type and hardware support
std::unique_ptr<IOReader> Reference::createIOReader(const std::string& fileName) {
    bool isGz = PathUtil::isGzFile(fileName);
    if (isGz) {
        // Choose gzip reader based on SIMD support
        bool useFastGzReader = false;
#ifdef __SSE4_2__
        useFastGzReader = Hardware::isSupportSimd();
#endif
        if (useFastGzReader) {
            // Intel ISA-L accelerated reader
            // FastGzFileReader requires a non-const string reference
            std::string mutableFileName = fileName;
            return std::make_unique<FastGzFileReader>(mutableFileName);
        } else {
            // Standard gzip reader
            return std::make_unique<GzFileReader>(fileName, 1); // Use single thread reading
        }
    } else {
        // Plain text reader
        return std::make_unique<FileReader>(fileName);
    }
}

/* Serialize the already-built squash; parsing lives in initSquashFromFasta. */
bool Reference::writeNiFile(const std::string& niFile) {
    if (refGeneSquash == nullptr || refGeneSquashlen <= 0) {
        LOG_ERROR("Squash is not built yet, nothing to write into %s.", niFile.c_str());
        return false;
    }

    FileWriter niWriter(niFile);
    if (niWriter.openIO() != 0) {
        LOG_ERROR("Failed to open NI file for writing: %s", niFile.c_str());
        return false;
    }

    int64_t offset;
    if (!writeNiFileHeader(niWriter, offset)) {
        niWriter.closeIO();
        return false;
    }

    niWriter.writeIO(refGeneSquash, refGeneSquashlen);
    niWriter.writeIOAt(offset, (uint8_t *)(&refGeneSquashlen), sizeof(refGeneSquashlen));

    md5_ctx md5;
    uint8_t digest[16];
    md5_init(&md5);
    md5_update(&md5, refGeneSquash, refGeneSquashlen);
    md5_final(digest, &md5);
    std::string md5String = md5_to_string(digest);

    if (!writeNiFileMetadata(niWriter, fastaChecksum, fastaLength, md5String)) {
        niWriter.closeIO();
        return false;
    }
    niWriter.closeIO();

    if (0 != niWriter.getWriteError()) {
        LOG_ERROR("Write ni file %s failed: %s", niFile.c_str(), strerror(niWriter.getWriteError()));
        return false;
    }
    return true;
}

bool Reference::writeNiFileMetadata(FileWriter& niWriter, const std::string& refGeneMd5, int64_t refGeneLen, std::string& md5) {
    std::string refename = PathUtil::getFileName(refGeneFile);
    Json::Value meta;
    meta["refe_name"] = refename;
    meta["refe_len"] = Json::Value::UInt64(refGeneLen);
    meta["refe_orgfile_md5"] = refGeneMd5; /* MD5 of the reference file as read. */
    meta["ni_data_md5"] = md5;  /* MD5 of the NI file payload. */
    fastaLength = refGeneLen;
    fastaChecksum = refGeneMd5;

    coder_json cmeta;
    std::string metaStr;
    cmeta.encoder(meta, metaStr);
    niWriter.writeIO((uint8_t *)(metaStr.c_str()), metaStr.length());
    return true;
}

/* Write the index only to the caller's path; no implicit caching. */
bool Reference::makeNiFile(const std::string& niFile) {
    if (!initSquashFromFasta()) {
        LOG_ERROR("Build squash from %s failed.", refGeneFile.c_str());
        return false;
    }
    return writeNiFile(niFile);
}

/* Identify by file header, not suffix. */
bool Reference::isNiFile(const std::string& file) {
    const size_t magicLen = PBGZ_FILE_MAGIC.length();
    const size_t headLen = magicLen + 3 + 2; /* magic + version + "ni" */
    uint8_t head[64];
    if (headLen > sizeof(head)) {
        return false;
    }

    FILE* fp = fopen(file.c_str(), "rb");
    if (fp == nullptr) {
        return false;
    }
    size_t got = fread(head, 1, headLen, fp);
    fclose(fp);
    if (got != headLen) {
        return false;
    }
    return 0 == memcmp(head, PBGZ_FILE_MAGIC.c_str(), magicLen) &&
           0 == memcmp(head + magicLen + 3, "ni", 2);
}

bool Reference::makeIndex() {
    Timer cost_ms(true);
    BaseGroupHash *bgHash;
    HashTable hashTable;
    std::string tips;
    /* Per-bucket write cursor into hash_buff, for concurrent writes. */
    uint32_t *hashBucketCurpos;

    hashBucketCurpos = MemoryUtil::safeAlloc<uint32_t>(hashBuckets);
    /* Fall back to the FASTA if the index cannot be trusted. */
    if (!loadSquashAndMatched()) {
        MemoryUtil::safeFree(hashBucketCurpos);
        return false;
    }

    cost_ms.reset();
    makeIndexFetchBaseGroup(bgHash);
    makeIndexCalcHashTableSize(hashTable);
    makeIndexInitHashTable(hashTable, hashBucketCurpos);
    makeIndexBuildHashTable(bgHash, hashBucketCurpos);
    makeIndexSortHashTable();
    tips = "elapsed ms: ";
    tips += std::to_string((int64_t)(cost_ms.elapsed()));
    if (guardBar != nullptr) {
        guardBar->done(tips);
    }

    MemoryUtil::safeFree(bgHash);
    MemoryUtil::safeFree(hashBucketCurpos);
    return true;
}

/* Load the reference squash and allocate the matched tracking buffer. */
bool Reference::loadSquashAndMatched() {
    const int64_t supportMax = ((int64_t)2 << 30) * baseGroupStep;
    /* Fall back to the FASTA if the index cannot be trusted. */
    bool squashReady = !niIndexFile.empty() && loadFromNiFile(niIndexFile);
    if (!squashReady) {
        squashReady = initSquashFromFasta();
    }
    if (!squashReady) {
        LOG_ERROR("initialize reference failed");
        return false;
    }
    if ((refGeneSquashlen << 2) >= supportMax) {
        LOG_ERROR("Reference max support %lu(M), current size %ld(M)\n",
                   supportMax / 1024 / 1024, (refGeneSquashlen << 2) / 1024 / 1024);
    } else {
        fprintf(stderr, "Reference max support %lu(M), current size %ld(M)\n",
                supportMax / 1024 / 1024, (refGeneSquashlen << 2) / 1024 / 1024);
    }

    refGeneSquashMatchedlen = refGeneSquashlen; /* One byte per squash byte: matched flag. */
    /* safeAlloc zeroes the buffer, i.e. nothing matched yet. */
    refGeneSquashMatched = MemoryUtil::safeAlloc<std::atomic<uint8_t>>(refGeneSquashMatchedlen);
    return true;
}

/* Load only the squash, without building the hash table. */
bool Reference::makeSquashIndex() {
    Timer cost_ms(true);
    if (!loadSquashAndMatched()) {
        return false;
    }
    fprintf(stderr, "Reference squash ready, elapsed ms: %ld\n", (long)cost_ms.elapsed());
    return true;
}

void Reference::makeIndexFetchBaseGroup(BaseGroupHash* &bgHash) {
    int64_t pcnt = this->parallel;
    const uint32_t bg_len = baseGroupLen;
    const uint32_t *actg_stretch_tab = actgStretch;
    const int32_t hmask = hashMask;

    int64_t total = (refGeneSquashlen << 2) / baseGroupStep;
    int64_t each = total / pcnt;
    int64_t remain = (total - (each * pcnt));
    uint8_t* p = this->refGeneSquash;
    bgHash = MemoryUtil::safeAlloc<BaseGroupHash>(total);
    hashBucketCnt = MemoryUtil::safeAlloc<uint32_t>(hashBuckets);
    spinlock* slocks = MemoryUtil::safeAlloc<spinlock>(hashBuckets);
    uint32_t* hb_cnt = hashBucketCnt;
    BaseGroupHash* bhash = bgHash;
    BaseGroupHash* bhash_start = bgHash;

    std::vector<std::thread> tpools;
    int64_t current, offset_start;
    int64_t *pnn = nullptr;
    for (int64_t n = 0; n < pcnt; n++) {
        current = (n + 1 == pcnt) ? (each + remain) : each; /* Keys processed by this thread */
        if ((n + 1) == pcnt) {
            gbTotal = 0;
            gbCurrent = current;
            guardBar = MemoryUtil::safeNewClass<GuardBar>(gbCurrent, gbTotal, "Building index from reference");
            guardBar->start();
            pnn = &gbTotal;
        }

        offset_start = bhash - bhash_start;

        /* bg_len is a compile-time constant here, so no capture needed. */
        tpools.push_back(std::thread([p, current, offset_start, hmask, &bhash_start,
                                      &actg_stretch_tab, &hb_cnt, &slocks, pnn]() {
            int64_t n = 0, m = 0;
            uint8_t *s = p;
            uint32_t kpos = bg_len >> 1;
            int64_t align4 = bg_len >> 2 << 2;
            uint32_t hash32, curr_bucket;
            int64_t offset = offset_start, pos;
            int64_t *pnnn = (pnn) ? pnn : (&n);
            uint64_t xsquash;

            /* Size scratch buffers from the compile-time baseGroupLen to avoid VLAs. */
            static constexpr uint32_t len_bgs = (baseGroupLen >> 2) + ((baseGroupLen & 0x3) ? 1 : 0);
            const char actg4[4] = {'A', 'C', 'T', 'G'};
            char actg_bg[baseGroupLen + 1];
            char actg_bg_pair[baseGroupLen + 1]; /* Complementary base */
            char actg_bgs[len_bgs];    /* squash base group */
            char actg_bgs_pair[len_bgs];

            for (*pnnn = 0; *pnnn < current; *pnnn = (*pnnn) + 1) {
                /* step 1 : get base group */
                for (m = 0; m < align4; m += 4) {
                    *((uint32_t *)(actg_bg + m)) = actg_stretch_tab[*s++];
                }
                switch (bg_len - align4) {
                case 0:
                    break;
                case 1:
                    m = align4;
                    *(actg_bg + m) = actg4[(((*s) >> 6) & 0x3)];
                    s++;
                    break;
                case 2:
                    m = align4;
                    *(actg_bg + m++) = actg4[(((*s) >> 6) & 0x3)];
                    *(actg_bg + m) = actg4[(((*s) >> 4) & 0x3)];
                    s++;
                    break;
                case 3:
                    m = align4;
                    *(actg_bg + m++) = actg4[(((*s) >> 6) & 0x3)];
                    *(actg_bg + m++) = actg4[(((*s) >> 4) & 0x3)];
                    *(actg_bg + m) = actg4[(((*s) >> 2) & 0x3)];
                    s++;
                    break;
                default:
                    break;
                }

                /* step 2: get base group pair */
                actgPair((uint8_t *)actg_bg_pair, (uint8_t *)actg_bg, bg_len);

                /* step 3:  save current base group with direction*/
                if (actg_bg[kpos] < actg_bg_pair[kpos]) {
                    actgSquash((uint8_t *)actg_bg_pair, bg_len, (uint8_t *)actg_bgs_pair);
                    xsquash = *((uint64_t *)(actg_bgs_pair));
                    xsquash &= 0xFCFFFFFFFFFFFFFF;
                    hash32 = (uint32_t) CityHash64((const char *)(&xsquash), len_bgs);
                    pos = offset | 0x80000000;
                } else {
                    actgSquash((uint8_t *)actg_bg, bg_len, (uint8_t *)actg_bgs);
                    xsquash = *((uint64_t *)(actg_bgs));
                    xsquash &= 0xFCFFFFFFFFFFFFFF;
                    hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
                    pos = offset;
                }
                curr_bucket = hash32 & hmask;
                (bhash_start + offset)->hashBucket = curr_bucket;
                (bhash_start + offset)->baseGroupPos = pos;
                offset++;
                {
                    spinlock &sl = slocks[curr_bucket];
                    spin_lock(&sl);
                    hb_cnt[curr_bucket]++;
                    spin_unlock(&sl);
                }
            }
        }));

        p += (current * (baseGroupStep >> 2)); /* Advances past the bases consumed by this thread. */
        bhash += current;
    }

    for (int64_t n = 0; n < pcnt; n++) {
        if (tpools[n].joinable()) {
            tpools[n].join();
        }
    }
    free(slocks);
}

void Reference::makeIndexCalcHashTableSize(HashTable& hashTable) {
    uint32_t pcnt = this->parallel;
    /* First segment: hash buffer content and length; second: bucket length. */
    hashTable.resize(pcnt);

    int64_t total = hashBuckets;
    int64_t each = hashBuckets / pcnt;
    int64_t remain = (total - (each * pcnt));
    uint32_t* p = hashBucketCnt;

    int64_t current;
    std::vector<std::thread> tpools;
    for (uint32_t n = 0; n < pcnt; n++) {
        current = (n + 1 == pcnt) ? (each + remain) : each;

        std::pair<std::pair<uint32_t *, uint32_t>, uint32_t> &hash_buffer = hashTable[n];
        tpools.push_back(std::thread([p, current, &hash_buffer]() {
            int64_t len_hash = 0, len_bucket = 0;
            int64_t n = 0;
            uint32_t current_cnt;

            for (n = 0; n < current; n++) {
                current_cnt = *(p + n);
                len_bucket++;
                /* Multi-value buckets reserve a slot for their buffer offset. */
                len_hash += (current_cnt <= 1) ? 1 : (current_cnt + 1);
            }
            hash_buffer.first.second = len_hash;
            hash_buffer.second = len_bucket;
        }));

        p += current;
    }

    for (uint32_t n = 0; n < pcnt; n++) {
        if (tpools[n].joinable()) {
            tpools[n].join();
        }
    }
}

void Reference::makeIndexInitHashTable(const HashTable& hashTable, uint32_t* &hashBucketCurPos) {
    uint32_t hashBufflen = 0;
    uint32_t hashBucketlen = 0;
    uint32_t n;
    for (n = 0; n < hashTable.size(); n++) {
        hashBufflen += hashTable[n].first.second;
        hashBucketlen += hashTable[n].second;
    }

    hashTableBuffer = MemoryUtil::safeAlloc<uint32_t>(hashBufflen);

    int64_t total = hashBuckets;
    uint32_t pcnt = this->parallel;
    int64_t each = hashBuckets / pcnt;
    int64_t remain = (total - (each * pcnt));
    uint32_t* phbucket = hashBucketCnt;
    uint32_t* phbuff = hashTableBuffer;
    uint32_t* phbuffStart = hashTableBuffer;
    uint32_t* phbuffConflict = hashTableBuffer + hashBucketlen;
    uint32_t* pbucketNext = hashBucketCurPos;

    int64_t current;
    std::vector<std::thread> tpools;
    for (n = 0; n < pcnt; n++) {
        current = (n + 1 == pcnt) ? (each + remain) : each;
        tpools.push_back(std::thread([phbucket, phbuff, phbuffConflict, phbuffStart, pbucketNext, &current]() {
            uint32_t n = 0;
            uint32_t current_cnt;
            uint32_t *p = phbuff;
            uint32_t *pc = phbuffConflict;
            uint32_t *pn = pbucketNext;

            for (n = 0; n < current; n++) {
                current_cnt = *(phbucket + n);
                switch (current_cnt)
                {
                case 0:
                    pn++;
                    p++;
                    break;
                case 1:
                    *pn++ = p - phbuffStart; /* No conflict: bucket points directly into the hash buffer. */
                    p++;
                    break;
                default:
                    *pn++ = pc - phbuffStart;
                    *p = pc - phbuffStart; /* Conflict: bucket stores the offset of its conflict buffer. */
                    pc += current_cnt;      /* Reserve one offset per hash value. */
                    p++;
                    break;
                }
            }
        }));

        pbucketNext += current;
        phbucket += current;
        phbuff += hashTable[n].second;
        phbuffConflict += hashTable[n].first.second - hashTable[n].second;
    }

    for (n = 0; n < pcnt; n++) {
        if (tpools[n].joinable()) {
            tpools[n].join();
        }
    }
}

void Reference::makeIndexBuildHashTable(const BaseGroupHash* bgHash, uint32_t* &hashBucketCurPos) {
    uint32_t *phbuffStart, id;
    std::vector<std::thread> tpools;
    uint32_t *hb_cnt = hashBucketCnt;
    int64_t n, pcnt = this->parallel;
    uint32_t *hb_curpos = hashBucketCurPos;

    int64_t total = (refGeneSquashlen << 2) / baseGroupStep;
    const BaseGroupHash *h = (BaseGroupHash*)bgHash;
    phbuffStart = hashTableBuffer;

    for (n = 0; n < pcnt; n++) {
        id = n;
        tpools.push_back(std::thread([h, total, &phbuffStart, &hb_curpos, &hb_cnt, id, pcnt]() {
            Timer cost_ms(true);
            uint32_t n = 0, bucket, pos;

            for (n = 0; n < total; n++)
            {
                bucket = (h + n)->hashBucket;
                if ((bucket % pcnt) == id) {
                    switch (hb_cnt[bucket])
                    {
                    case 0:
                        break;
                    case 1:
                        pos = hb_curpos[bucket];
                        *(phbuffStart + pos) = (h + n)->baseGroupPos;
                        break;
                    default:
                        pos = hb_curpos[bucket];
                        *(phbuffStart + pos) = (h + n)->baseGroupPos;
                        hb_curpos[bucket]++; /* Advance the bucket cursor. */
                        break;
                    }
                }
            }
        }));
    }

    for (n = 0; n < pcnt; n++) {
        if (tpools[n].joinable()) {
            tpools[n].join();
        }
    }
}

void Reference::makeIndexSortHashTable() {
    uint32_t *phb_cnt, *phb_cnt_start;
    uint32_t *phb_buff, *phb_buff_start, offset = 0;
    int64_t each, current, remain, total;
    int64_t n, pcnt = this->parallel;
    std::vector<std::thread> tpools;

    total = hashBuckets;
    each = total / pcnt;
    remain = (total - (each * pcnt));
    phb_buff_start = phb_buff = hashTableBuffer;
    phb_cnt_start = phb_cnt = hashBucketCnt;

    for (n = 0; n < pcnt; n++) {
        current = (n + 1 == pcnt) ? (each + remain) : each;
        tpools.push_back(std::thread([phb_buff, current, phb_buff_start, &phb_cnt_start, offset]() {
            uint32_t current_cnt;
            uint32_t n, *p;
            uint32_t *pstart = phb_buff_start;

            for (n = 0; n < current; n++) {
                current_cnt = *(phb_cnt_start + offset + n);
                if (current_cnt > 1) {
                    if (current_cnt < 255) {
                        p = pstart + phb_buff[n];
                        std::sort(p, p + current_cnt, [](const uint32_t &p1, const uint32_t &p2) -> bool {
                             return p1 < p2;
                        }); // Stable ordering, avoids out-of-bounds on equal elements
                    }
                    *(phb_cnt_start + offset + n) = std::min((uint32_t)16, current_cnt);
                }
            }
        }));
        offset += current;
        phb_buff += current;
    }

    for (n = 0; n < pcnt; n++) {
        if (tpools[n].joinable()) {
            tpools[n].join();
        }
    }
}

void Reference::dumpHashTable() {
    FILE *fp = fopen("hash_table", "wb");
    if (fp == nullptr) {
        LOG_ERROR("Failed to open hash_table file for writing");
        return;
    }

    for (int32_t n = 0; n < hashBuckets; n++) {
        if (hashBucketCnt[n] <= 0) {
            continue;
        }
        fprintf(fp, "bucket %u : ", n);
        uint32_t *p = (hashBucketCnt[n] > 1) ? (hashTableBuffer + hashTableBuffer[n]) : (hashTableBuffer + n);

        for (uint32_t m = 0; m < hashBucketCnt[n]; m++) {
            fprintf(fp, "%u,", *(p + m));
        }
        fprintf(fp, "\n");
    }
    fflush(fp);
    fclose(fp);
}

void Reference::dumpSquash() {
    FILE *fp = fopen("squash.txt", "wb");
    if (fp == nullptr) {
        LOG_ERROR("Failed to open squash file for writing");
        return;
    }

    for(int64_t idx = 0; idx < refGeneSquashlen; idx++) {
        fprintf(fp, "%X", refGeneSquash[idx]);
        if (idx % 64 == 0 && idx != 0) {
            fprintf(fp,"\n");
        }
    }

    fflush(fp);
    fclose(fp);
}


uint32_t Reference::getBaseGroupLength() const {
    return this->baseGroupLen;
}

uint32_t Reference::getBaseGroupStep() const {
    return this->baseGroupStep;
}

const uint32_t* Reference::queryPosition(const uint32_t &hash, uint32_t &length) {
    uint32_t hbucket = hash & hashMask;
    length = hashBucketCnt[hbucket];
    return (length == 1) ? (hashTableBuffer + hbucket) : (hashTableBuffer + hashTableBuffer[hbucket]);
}

void Reference::updateMatchedGene(uint64_t actgPos, uint32_t matchLength) {
    if (matchLength == 0) {
        return;
    }
    const uint64_t sposStart = actgPos >> 2;
    const uint64_t sposEnd = (actgPos + matchLength - 1) >> 2;
    for (uint64_t n = sposStart; n <= sposEnd; ++n) {
        refGeneSquashMatched[n].store(1, std::memory_order_relaxed);
    }
}

/* Get reference squash buffer */
const uint8_t* Reference::getSquash() const {
    return this->refGeneSquash;
}

/* Get reference squash buffer length */
int64_t Reference::getSquashLength() const {
    return this->refGeneSquashlen;
}

/* Get fasta file name */
const std::string& Reference::getFastaFileName() const {
    return this->refGeneFile;
}

/* Get fasta file content length */
int64_t Reference::getFastaLength() const {
    return this->fastaLength;
}

/* Get fasta file content MD5 */
const std::string& Reference::getFastaChecksum() const {
    return this->fastaChecksum;
}

/* Clear unmatched reference to zero */
void Reference::sanitizeRefSquash(int64_t startSquashPos, int64_t len) {
    uint64_t e = startSquashPos + len;
    for (uint64_t n = startSquashPos; n < e; n++) {
        *(refGeneSquash + n) = refGeneSquashMatched[n].load(std::memory_order_relaxed)
                                   ? (*(refGeneSquash + n)) : 0;
    }
}

/* Get ACTG bases of specified length at specified position */
void Reference::getStretchActg(uint8_t *out, uint32_t outLen, uint64_t actgPos) {
    uint64_t squashPos = actgPos >> 2;
    uint8_t* p = refGeneSquash + squashPos;
    uint8_t ch = *p;

    uint32_t offset  = 0;
    const char actg4[4] = {'A', 'C', 'T', 'G'};
    /* Left unaligned part */
    switch (actgPos & 0x3)
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    case 2:
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    case 3:
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    default:
        break;
    }
    if (offset == outLen) {
        return;
    }

    /* Aligned part */
    uint32_t lenNeed = (outLen - offset) >> 2;
    for (uint32_t n = 0; n < lenNeed; n++) {
        *((uint32_t *)(out + offset)) = actgStretch[*p];
        offset += 4;
        p++;
    }
    if (offset == outLen) {
        return;
    }

    /* Right unaligned part */
    lenNeed = outLen - offset;
    ch = *p;
    switch (lenNeed & 0x3)
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        break;
    case 2:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        break;
    case 3:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        break;
    default:
        break;
    }
}

/* Convert to ACTG based on the last 2 bits of each byte */
void Reference::getActgFrom2Bits(const uint8_t *src2Bits, uint32_t src2BitsLen, uint8_t *dstActg) {
    uint8_t *s = (uint8_t *)src2Bits, *p;
    uint32_t align4 = src2BitsLen >> 2 << 2, offset = 0;
    const uint8_t actg4[4] = {'A', 'C', 'T', 'G'};

    uint32_t n = 0;
    for ( ;n < align4; n += 4) {
        p = s + n;
        *((uint32_t *)(dstActg + offset)) = actgStretch[((*p) << 6) | ((*(p + 1)) << 4) | ((*(p + 2)) << 2) | (*(p + 3))];
        offset += 4;
    }

    if (n == src2BitsLen) {
        return;
    }

    for (; n < src2BitsLen; n++) {
        *(dstActg + offset++) = actg4[*(src2Bits + n)];
    }
}

/* Get 2-bit squash bases for the given range. */
void Reference::getStretch2Bits1Char(uint8_t *out, uint32_t outLen, uint64_t actgPos) {
    uint32_t offset = 0;
    uint64_t squashPos = actgPos >> 2;

    uint8_t *p = refGeneSquash + squashPos;
    uint8_t ch = *p;

    /* Left unaligned part */
    switch (actgPos & 0x3)
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = (ch >> 4) & 0x3;
        *(out + offset++) = (ch >> 2) & 0x3;
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    case 2:
        *(out + offset++) = (ch >> 2) & 0x3;
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    case 3:
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    default:
        break;
    }
    /* The unaligned part may write up to 3 bytes even when outLen is smaller. */
    if (offset >= outLen) {
        return;
    }

    /* Aligned part */
    uint32_t lenNeed = (outLen - offset) >> 2;
    for (uint32_t n = 0; n < lenNeed; n++) {
        *((uint32_t *)(out + offset)) = actgStretch2bits[*p];
        offset += 4;
        p++;
    }
    if (offset >= outLen) {
        return;
    }

    /* Right unaligned part */
    lenNeed = outLen - offset;
    ch = *p;
    switch (lenNeed & 0x3)
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = (ch >> 6) & 0x3;
        break;
    case 2:
        *(out + offset++) = (ch >> 6) & 0x3;
        *(out + offset++) = (ch >> 4) & 0x3;
        break;
    case 3:
        *(out + offset++) = (ch >> 6) & 0x3;
        *(out + offset++) = (ch >> 4) & 0x3;
        *(out + offset++) = (ch >> 2) & 0x3;
        break;
    default:
        break;
    }
}
