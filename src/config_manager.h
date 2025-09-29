#pragma once

#include "log/logger.h"

class ConfigManager{

public:
    LogLevel getLogLevel();

    LogAppender getLogAppender();

    static ConfigManager& getInstance();

private:
    ConfigManager();
    ConfigManager(ConfigManager&) = delete;
    const ConfigManager& operator=(const ConfigManager&) = delete;

private:
    LogLevel logLevel;
    LogAppender logAppender;
    std::string logFile;
    uint32_t maxLogFileSize;
    std::string logIpAddress;
    uint16_t logPort;
    std::string logUrl;
    int runMode;   // 运行模式，Debug/Release
};


