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
        {1, 512 << 10},    // 512K
        {2, 2 << 20},      // 2M
        {3, 8 << 20},      // 8M 
        {4, 16 << 20},     // 16M
        {5, 64 << 20},     // 64M
        {6, 128 << 20},    // 128M
        {7, 256 << 20},    // 256M
        {8, 512 << 20},    // 512M
        {9, 1024 << 20},   // 1024M
    };
};


