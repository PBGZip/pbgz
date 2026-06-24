/*
 * binary_actuator.cpp - Source file for pbgz project
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

#include "binary_actuator.h"
#include "coder_io.h"
#include "coder_fc.h"
#include "coder_bwt_cm.h"
#include "utils/md5_util.h"
#include "coder_json.h"
#include "log/logger.h"


int32_t BinaryCodecActuator::compress() {
    coder_io io (outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    uint32_t srcLength = inBlockPtr->getDataLen();
    if (srcLength <= FC_MIN_LEN || srcLength > FC_MAX_LEN) {
        coder_bwt_cm coder(&io);
        coder.encode_line(inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
        coder.encode_flush();
    } else {
        coder_fc coder(&io);
        coder.encode_line(inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
        coder.encode_flush();
    }

    outBlockPtr->setDataLen(io.data_len);
    if (inBlockPtr->getBlockType() != BINARY) {
        LOG_DEBUG("Reset block type to %d", BINARY);
        inBlockPtr->setBlockType(BINARY);
    }

    Json::Value subMeta;
    Json::Value meta;
    subMeta["srclen"] = srcLength;
    subMeta["dstlen"] = io.data_len;
    subMeta["coder"] = io.meta;
    meta["streams"] = subMeta;

    std::string md5;
    calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
    meta["md5"] = md5;

    // Compress meta information
    coder_json metaCoder;
    int64_t metaLength = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
    if (metaLength < 0) {
        LOG_ERROR("Failed to compress meta information");
        return -1;
    }
    outBlockPtr->setMetaLen((uint32_t)metaLength);
    LOG_DEBUG("Compress binary from %d to %d.", srcLength, io.data_len);
    return 0;
}

int32_t BinaryCodecActuator::decompress() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter , inBlockPtr or outBlockPtr is nullptr");
        return -1;
    }

    // First parse out meta information
    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);
    uint32_t decoderLen = 0;
    uint32_t decSrcLen = meta["streams"]["srclen"].asUInt();
    uint32_t decDstLen = meta["streams"]["dstlen"].asUInt();
    if (inBlockPtr->getDataLen() != decDstLen) {
        LOG_ERROR("Dst length not match.");
        return -1;
    }
    coder_io io(inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
    if (meta["streams"]["coder"]["magic"] == "coder_bwt_cm") {
        coder_bwt_cm coder(&io);
        decoderLen = coder.decode_line(outBlockPtr->getBuffer(), decSrcLen, UINT8_MAX, false);
    } else if (meta["streams"]["coder"]["magic"] == "coder_fc") {
        io.meta = meta["streams"];
        coder_fc coder(&io);
        decoderLen = coder.decode_line(outBlockPtr->getBuffer(), decSrcLen, UINT8_MAX, false);
    } else {
        LOG_ERROR("Not support yet for coder %s", meta["streams"]["coder"]["magic"].asString().c_str());
        return -1;
    }

    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + decoderLen);
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());
    
    // Check checksum of source content
    std::string md5;
    calcMd5sum(md5, outBlockPtr->getBuffer(), outBlockPtr->getDataLen());
    if (md5 != meta["md5"].asString()) {
        LOG_ERROR("check md5 failed, blockid = %d", outBlockPtr->getBlockId());
        return -1;
    }

    return 0;
}
