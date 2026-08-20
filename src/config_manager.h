/*
 * config_manager.h - Header file for pbgz project
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

#include <map>

#include "log/logger.h"
#include "pbgz_types.h"


class ConfigManager{

public:
    LogLevel getLogLevel();

    LogAppender getLogAppender();

    void init(PbgzParameter& para);

    static ConfigManager& getInstance();

    std::string& getLogFileName();

    uint32_t getBlockSizeByCompressLevel(uint8_t compressLevel) {
        const uint32_t defaultBlockSize = 32 << 20;
        if (compressLevel < 1 || compressLevel > 9) {
            return defaultBlockSize;
        }
        auto item = compressLevelBlockMap.find(compressLevel);
        if (item == compressLevelBlockMap.end()) {
            return defaultBlockSize;
        }
        return item->second;
    }

private:
    ConfigManager();
    ConfigManager(ConfigManager&) = delete;
    const ConfigManager& operator=(const ConfigManager&) = delete;

private:
    LogLevel logLevel;
    LogAppender logAppender;
    std::string logFile;
    std::map<uint8_t, uint32_t> compressLevelBlockMap = {
        {1, 256 << 10},    // 256K
        {2, 512 << 10},      // 512K
        {3, 1 << 20},      // 1M 
        {4, 2 << 20},     // 2M
        {5, 4 << 20},     // 4M
        {6, 8 << 20},     // 8M
        {7, 16 << 20},    // 16M
        {8, 32 << 20},    // 32M
        {9, 64 << 20},    // 64M
    };
};


