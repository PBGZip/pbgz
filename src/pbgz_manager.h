/*
 * pbgz_manager.h - Header file for pbgz project
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

#include <vector>
#include <string>
#include <string.h>

#include "coder/coder.h"
#include "log/logger.h"
#include "io_block.h"
#include "pbgz_types.h"
#include "utils/timer.h"

class PbgzManager {
public:
    PbgzManager() {
        totalReadLen = 0;
        totalWriteLen = 0;
    }

    PbgzManager& operator=(const PbgzManager& item) = delete;
    PbgzManager(const PbgzManager& item) = delete;

    ~PbgzManager() {

    }

    void exitProc(int errorCode, const char* errorMessage);

    static PbgzManager& getInstance();

    std::string getVersion();

    std::vector<char> getVersionAsArray();

    void updateReadDataLen(RoughIOBlock* blockPtr);
    
    void updateWriteDataLen(RoughIOBlock* blockPtr);

    void printHeadInfo(PbgzParameter& para);

    void printFileType(BlockType blockType);

    void printTailInfo(Timer costTime, bool isPrintRatio);

    // Add output file to manage list, will be cleaned up on failure
    void addOutputFile(const std::string& fileName);
    
    // Clean up all output files that were marked for cleanup
    void cleanupOutputFiles();

    void printBufferContent(uint8_t* buffer, uint32_t bufferLen) {
        char temp[2048 + 1] = {0};
        bufferLen = std::min<uint32_t>(bufferLen, 2048);
        memcpy(temp, (char*)buffer, bufferLen);
        temp[2048] = 0;
        LOG_DEBUG("%s", temp);
    }

    void printBufferContentHex(uint8_t* buffer, uint32_t bufferLen) {
        bufferLen = std::min<uint32_t>(bufferLen, 2048);
        for (uint32_t i = 0; i < bufferLen; ++i) {
            fprintf(stderr, "%02X", buffer[i]);
        }
        fprintf(stderr, "\n");
    }

    void printBufferContentBinary(uint8_t* buffer, uint32_t bufferLen) {
        for (uint32_t j = 0; j < bufferLen; j++) {
            // Output src[i] in binary format
            for (int n = 7; n >= 0; n--) {
                fprintf(stderr, "%d", (buffer[j] >> n) & 1);
            }
            fprintf(stderr, "\t");

            if ((j+1) % 8 == 0) {
                fprintf(stderr, "\n");
            }
        }
        fprintf(stderr, "\n");
    }
private:
    void updateDataInfo();

    std::vector<std::pair<std::string, bool>> outfiles; 

    uint64_t totalReadLen;

    uint64_t totalWriteLen;
};

void pbgzExitProc(int errorCode, const char* errorMessage);

void coderLog(int logLevel, const char* logMessage);

int32_t powerof2Proximal(int32_t i);
