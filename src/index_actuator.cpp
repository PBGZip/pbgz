/*
 * index_actuator.cpp - Index actuator implementation file
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
 * FITNESS FOR A PARTICULAR PURPOSE OR NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>

#include "index_actuator.h"
#include "pbgz_engine.h"
#include "pbgz_index.h"   /* SamIndex, used by the index paths below */
#include "sam_info.h"
#include "log/logger.h"
#include "coder/coder_io.h"
#include "coder/coder_json.h"
#include "coder/coder_bwt_cm.h"
#include "coder_factory.h"

IndexActuator::IndexActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine)
    : Actuator(inPtr, outPtr, engine), flagDecoder(nullptr), chrDecoder(nullptr), posDecoder(nullptr) {
    headEndLine = 0;
}

IndexActuator::~IndexActuator() {
    flagDecoder.reset();
    chrDecoder.reset();
    posDecoder.reset();
}

int32_t IndexActuator::initial() {
    if (inBlockPtr == nullptr) {
        return -1;
    }

    if (!supports(inBlockPtr->getBlockType())) {
        LOG_ERROR("Block type %d cannot be indexed.", inBlockPtr->getBlockType());
        return -1;
    }

    // Parse meta information
    coder_json metaCoder;
    Json::Value meta;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);

    // Parse SAM header
    if (meta.isMember("header")) {
        if (parseHeader(meta) != 0) {
            LOG_ERROR("Parse SAM header failed");
            return -1;
        }
    }

    // Decode data fields and build index
    if (!meta.isMember("sam")) {
        LOG_INFO("No SAM info for field-by-field decompression");
        return 0;
    }

    Json::Value& samMeta = meta["sam"];
    Json::Value& streams = samMeta["streams"];
    uint32_t lineNum = samMeta["lines"].asUInt();

    int32_t readOffset = 0;
    if (meta.isMember("header")) {
        readOffset = meta["header"]["dstlen"].asUInt();
    }

    if (initDecoders(streams, readOffset) != 0) {
        LOG_ERROR("Initialize decoders failed");
        return -1;
    }

    if (skipUnneededFields(streams, readOffset) != 0) {
        LOG_ERROR("Skip unneeded fields failed");
        return -1;
    }

    if (decodeAndBuildIndex(lineNum) != 0) {
        LOG_ERROR("Decode and build index failed");
        return -1;
    }

    return 0;
}

int32_t IndexActuator::parseHeader(Json::Value& meta) {
    Json::Value& headerMeta = meta["header"];
    if (!headerMeta.isMember("srclen") || !headerMeta.isMember("dstlen") ||
        !headerMeta.isMember("lines") || !headerMeta.isMember("coder")) {
        LOG_ERROR("Invalid SAM header metadata for decompression");
        return -1;
    }

    headEndLine = headerMeta["lines"].asInt64();
    uint32_t dstLen = headerMeta["dstlen"].asUInt();

    // Create SAM file header decompressor from the recorded magic (BWT in archive, RANS in fast).
    coder_io headerIo(inBlockPtr->getBuffer(), dstLen, &ioErrSink, "SAM header");
    FieldDecoderArgs headerArgs;
    if (headerMeta["coder"].isMember("level")) {
        headerArgs.level = headerMeta["coder"]["level"].asInt();
    }
    std::shared_ptr<coder> headerDecoder =
        CoderFactory::makeFieldDecoder(headerMeta["coder"]["magic"].asString(), &headerIo,
                                       headerArgs);
    if (headerDecoder == nullptr) {
        LOG_ERROR("Unsupported SAM header coder: %s",
                  headerMeta["coder"]["magic"].asString().c_str());
        return -1;
    }

    // Decompress SAM file header data and parse chromosome info
    uint32_t lineCount = 0;
    uint8_t tempBuffer[1024];
    while (lineCount < headEndLine) {
        uint32_t decodedLen = headerDecoder->decode_line(tempBuffer, sizeof(tempBuffer), '\n', false);
        if (decodedLen == 0) {
            break;
        }
        std::string headStr = std::string((char*)tempBuffer, decodedLen);
        if (headStr.length() >= 3 && headStr.substr(0, 3) == "@SQ") {
            SamUtil::parseChromosomeInfo(headStr);
        }
        lineCount++;
    }

    return 0;
}

int32_t IndexActuator::initDecoders(const Json::Value& streams, int32_t& readOffset) {
    auto openField = [&](const Json::Value& fieldMeta, const char* name,
                         std::shared_ptr<coder_io>& io) -> std::shared_ptr<coder> {
        if (!fieldMeta.isMember("coder")) {
            LOG_ERROR("Index: field %s has no coder meta", name);
            return nullptr;
        }
        const uint32_t dstLen = fieldMeta["dstlen"].asUInt();
        io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, name);
        FieldDecoderArgs args;
        if (fieldMeta["coder"].isMember("level")) {
            args.level = fieldMeta["coder"]["level"].asInt();
        }
        AuxPayloadPtr posPrior;
        if (fieldMeta["coder"]["magic"].asString() == "coder_arith" && pbgzEngine != nullptr) {
            posPrior = pbgzEngine->getPosPrior();
            args.posPrior = posPrior.get();
        }
        std::shared_ptr<coder> dec =
            CoderFactory::makeFieldDecoder(fieldMeta["coder"]["magic"].asString(), io.get(), args);
        if (dec == nullptr) {
            LOG_ERROR("Index: unsupported coder %s for field %s",
                      fieldMeta["coder"]["magic"].asString().c_str(), name);
            return nullptr;
        }
        readOffset += dstLen;
        return dec;
    };

    // Initialize ID field decoders if needed (field 0)
    if (streams.size() > 0 && streams[0]["field"].asUInt() == 0) {
        const Json::Value& idMeta = streams[0];
        const Json::Value& idStreamMeta = idMeta["streams"];
        for (uint32_t i = 0; i < idStreamMeta.size(); ++i) {
            readOffset += idStreamMeta[i]["dstlen"].asUInt();
        }
    }

    // FLAG (field 1): fixed-width number or text lines.
    if (streams.size() > 1 && streams[1]["field"].asUInt() == 1) {
        flagStringMode = streams[1].isMember("mode") && streams[1]["mode"].asString() == "string";
        flagDecoder = openField(streams[1], "FLAG", flagIo);
        if (flagDecoder == nullptr) {
            return -1;
        }
    }

    // RNAME (field 2): chromosome index, one uint16 per record.
    if (streams.size() > 2 && streams[2]["field"].asUInt() == 2) {
        chrDecoder = openField(streams[2], "RNAME", chrIo);
        if (chrDecoder == nullptr) {
            return -1;
        }
    }

    // POS (field 3): delta-varint, text, or fixed uint32 depending on the writer's mode.
    if (streams.size() > 3 && streams[3]["field"].asUInt() == 3) {
        const std::string mode = streams[3].isMember("mode")
                                     ? streams[3]["mode"].asString() : std::string();
        posStringMode = (mode == "string");
        posDeltaMode = (mode == "pos_delta");
        posDecoder = openField(streams[3], "POS", posIo);
        if (posDecoder == nullptr) {
            return -1;
        }
    }

    return 0;
}

int32_t IndexActuator::skipUnneededFields(const Json::Value& streams, int32_t& readOffset) {
    for (uint32_t idx = 4; idx < streams.size(); ++idx) {
        const Json::Value& fieldMeta = streams[idx];
        uint32_t fieldIdx = fieldMeta["field"].asUInt();

        if (fieldIdx == 9) {
            // Base field - skip it and any sub-streams
            if (fieldMeta.isMember("streams")) {
                const Json::Value& baseStreams = fieldMeta["streams"];
                for (uint32_t i = 0; i < baseStreams.size(); ++i) {
                    readOffset += baseStreams[i]["dstlen"].asUInt();
                }
            } else {
                readOffset += fieldMeta["totaldstlen"].asUInt();
            }
        } else if (fieldIdx == 10) {
            // Quality field - skip it and its sub-streams
            if (fieldMeta.isMember("streams")) {
                const Json::Value& qualStreams = fieldMeta["streams"];
                for (uint32_t i = 0; i < qualStreams.size(); ++i) {
                    readOffset += qualStreams[i]["dstlen"].asUInt();
                }
            } else {
                readOffset += fieldMeta["totaldstlen"].asUInt();
            }
        } else if (fieldIdx != 1 && fieldIdx != 2 && fieldIdx != 3) {
            readOffset += fieldMeta["dstlen"].asUInt();
        }
    }

    return 0;
}

int32_t IndexActuator::decodeAndBuildIndex(uint32_t lineNum) {
    SortKey lastKey = {0, 0};
    // Reserve space for validation
    sortKeys.reserve(lineNum > 2 ? 2 : lineNum);

    // Track for each chromosome: vector of items
    // Each item: first position, last position, count
    // Split when count >= 1000 and next position != last position
    std::map<uint16_t, std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>> chrBlockStats;
    const uint32_t MAX_SPLITS_PER_CHR = 1000;

    for (uint32_t lineNo = 0; lineNo < lineNum; ++lineNo) {
        uint16_t chrIndex = 0xFFFF;
        uint16_t flag = 0;
        uint32_t mapPos = 0;

        // Decode FLAG (field 1): text lines or fixed uint16.
        if (flagDecoder != nullptr) {
            if (flagStringMode) {
                uint8_t tempBuffer[64];
                int32_t decodedLen = flagDecoder->decode_line(tempBuffer, sizeof(tempBuffer), '\t', false);
                if (decodedLen > 0) {
                    size_t textLen = (size_t)decodedLen;
                    if (tempBuffer[textLen - 1] == '\t' || tempBuffer[textLen - 1] == '\n') {
                        --textLen;
                    }
                    if (textLen > 0) {
                        flag = (uint16_t)strtoul(std::string((char*)tempBuffer, textLen).c_str(),
                                                 nullptr, 10);
                    }
                }
            } else {
                uint8_t tempBuffer[2];
                int32_t decodedLen = flagDecoder->decode_line(tempBuffer, sizeof(tempBuffer), UINT8_MAX, false);
                if (decodedLen == sizeof(uint16_t)) {
                    memcpy(&flag, tempBuffer, sizeof(uint16_t));
                }
            }
        }

        // Decode RNAME (field 2) - chromosome index
        if (chrDecoder != nullptr) {
            uint8_t tempBuffer[2];
            int32_t decodedLen = chrDecoder->decode_line(tempBuffer, sizeof(tempBuffer), UINT8_MAX, false);
            if (decodedLen == sizeof(uint16_t)) {
                memcpy(&chrIndex, tempBuffer, sizeof(uint16_t));
            }
        }

        // Decode POS (field 3): delta-varint (current layout), text, or fixed uint32.
        if (posDecoder != nullptr) {
            if (posDeltaMode) {
                if (chrIndex != posLastChr) {
                    posDeltaPrev = 0;
                    posLastChr = chrIndex;
                }
                uint64_t delta = 0;
                int shift = 0;
                bool ok = true;
                for (;;) {
                    uint8_t b = 0;
                    if (posDecoder->decode_line(&b, 1, UINT8_MAX, false) != 1) {
                        ok = false;
                        break;
                    }
                    delta |= (uint64_t)(b & 0x7f) << shift;
                    if ((b & 0x80) == 0) {
                        break;
                    }
                    shift += 7;
                    if (shift >= 64) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) {
                    LOG_ERROR("Index: POS delta stream truncated at line %u", lineNo);
                    return -1;
                }
                posDeltaPrev += (int64_t)delta;
                mapPos = (uint32_t)posDeltaPrev;
            } else if (posStringMode) {
                uint8_t tempBuffer[32];
                int32_t decodedLen = posDecoder->decode_line(tempBuffer, sizeof(tempBuffer), '\t', false);
                if (decodedLen > 0) {
                    size_t textLen = (size_t)decodedLen;
                    if (tempBuffer[textLen - 1] == '\t' || tempBuffer[textLen - 1] == '\n') {
                        --textLen;
                    }
                    if (textLen > 0) {
                        mapPos = (uint32_t)strtoul(std::string((char*)tempBuffer, textLen).c_str(),
                                                   nullptr, 10);
                    }
                }
            } else {
                uint8_t tempBuffer[4];
                int32_t decodedLen = posDecoder->decode_line(tempBuffer, sizeof(tempBuffer), UINT8_MAX, false);
                if (decodedLen == sizeof(uint32_t)) {
                    memcpy(&mapPos, tempBuffer, sizeof(uint32_t));
                }
            }
        }

        // Validate ordering for mapped reads
        if ((flag & 0x04) == 0 && chrIndex != 0xFFFF && chrIndex != 0xFFFE) {
            SortKey currentKey = {chrIndex, mapPos};

            // Check if the current line is in ascending order
            if (currentKey < lastKey) {
                LOG_ERROR("SAM block %d is not sorted by chrIndex and mapPos: line %u (chr=%u,pos=%u) < line %u-1 (chr=%u,pos=%u)",
                    inBlockPtr->getBlockId(), lineNo, chrIndex, mapPos, lineNo,
                    lastKey.chrIndex, lastKey.mapPos);
                return -1;
            }

            if (lineNo < 2) {
                sortKeys.push_back(currentKey);
            }
            lastKey = currentKey;

            // Track chromosome statistics for this block
            auto& items = chrBlockStats[chrIndex];
            if (items.empty()) {
                // New chromosome, add first item: firstPos, lastPos, count = 1
                items.emplace_back(mapPos, mapPos, 1);
            } else {
                // Get current item's data
                uint32_t& lastPos = std::get<1>(items.back());
                uint32_t& count = std::get<2>(items.back());

                // Split when count reaches threshold and next position is different from last position
                if (count >= MAX_SPLITS_PER_CHR && mapPos != lastPos) {
                    // Start a new item with current position
                    items.emplace_back(mapPos, mapPos, 1);
                } else {
                    // Add to current item: update lastPos and increment count
                    lastPos = mapPos;
                    count++;
                }
            }
        }
    }

    // Add block-level index to SamIndex singleton
    for (const auto& chrPair : chrBlockStats) {
        uint16_t chrIndex = chrPair.first;
        const auto& items = chrPair.second;

        int64_t refPos = SamInfo::getInstance().getPositionByIndex(chrIndex);
        if (refPos == -1) {
            continue;
        }

        for (const auto& item : items) {
            uint32_t firstPos = std::get<0>(item);
            uint32_t count = std::get<2>(item);
            SamIndex::getInstance().addSamIndex(chrIndex, firstPos, count, inBlockPtr->getBlockId());
        }
    }

    return 0;
}

int32_t IndexActuator::process() {
    if (outBlockPtr == nullptr) {
        return -1;
    }

    sortKeys.clear();
    return 0;
}
