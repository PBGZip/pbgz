/*
 * pbgz_file_handler.cpp - CPP file for pbgz_v2 project
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

#include "pbgz_file_handler.h"
#include "log/logger.h"

int32_t PbgzFileReader::open() {
    // Implement the read logic for PBGZ file format
    // This is a placeholder implementation
    LOG_STDOUT(LOG_INFO, "Reading PBGZ file...");
    FILE* pFile = fopen(fileName.c_str(), "rb");
    if (!pFile) {
        LOG_STDOUT(LOG_ERROR, "Failed to open file: %s", fileName.c_str());
        return -1; // File open error
    }

    // Read the file header
    char buffer[7];
    size_t readLen = fread(buffer, sizeof(buffer), 1, pFile);
    if (readLen !=  sizeof(buffer)) {
        LOG_STDOUT(LOG_ERROR, "Failed to read file header from: %s", fileName.c_str());
        fclose(pFile);
        return -1; // File read error
    }

    if (strncmp(buffer, PBGZ_FILE_MAGIC.c_str(), 4) != 0) {
        // 安照16进制打印，自行解析
        LOG_STDOUT(LOG_ERROR, "%s is not a valid pbgz file, magic no is %X", fileName.c_str(), (uint32_t)buffer);
        fclose(pFile);
        return -1; // Invalid magic
    }

    LOG_STDOUT(LOG_INFO, "PBGZ file opened successfully: %s", fileName.c_str());

    PbgzFileHeader header;
    header.unserialize(reinterpret_cast<uint8_t*>(buffer), sizeof(buffer));
    fileHeaderMap[currentFileIndex++] = header;

    // Read file meta information
    PbgzFileMeta fileMeta;

    return 0; // Return success status
}