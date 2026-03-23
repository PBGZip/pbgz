/*
 * pbgz_file.cpp - cpp file for pbgz project
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

#include <iostream>

#include "pbgz_file.h"
#include "log/logger.h"
#include "coder_json.h"
#include "farmhash/src/farmhash.h"


std::string PbgzFileHeader::getVersionStr() {
    // Convert version array to string format "x.x.x"
    // Assuming version is in the format [major, minor, patch]
    // where each element is a single character representing a digit.
    // Adjust the conversion if the version format is different.
    if (version[0] == '\0') {
        return "Unknown Version";
    }
    // Ensure that version elements are valid characters    
    char buffer[12];
    snprintf(buffer, sizeof(buffer), "%d.%d.%d", version[0], version[1], version[2]);
    return std::string(buffer);
}

int32_t PbgzDataBlock::setBlockData(uint8_t* data, uint32_t length) {
    if (data == nullptr) {
        return -1;
    }
    pBlockData = data;
    blockDataLength = length;
    return 0;
}

void PbgzDataBlock::calcChecksum() {
    dataBlockChecksum = util::Hash64((char*)pBlockData, (size_t)blockDataLength);
    Json::StreamWriterBuilder writer;
    std::string jsonStr = Json::writeString(writer, dataMetaInfo);
    metaChecksum = util::Hash64(jsonStr);
}

int32_t PbgzDataBlock::verifyCheckSum() {
    if (dataBlockChecksum == util::Hash64((char*)pBlockData, (size_t)blockDataLength)) {
        Json::StreamWriterBuilder writer;
        std::string jsonStr = Json::writeString(writer, dataMetaInfo);
        if (metaChecksum == util::Hash64(jsonStr)) {
            return 0;
        }
        return -1;
    }
    return -1;
}
