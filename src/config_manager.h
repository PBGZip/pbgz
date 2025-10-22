#pragma once

#include "log/logger.h"
#include "pbgz_types.h"

class ConfigManager{

public:
    LogLevel getLogLevel();

    LogAppender getLogAppender();

    void init(PbgzParameter& para);

    static ConfigManager& getInstance();

    std::string& getLogFileName();

private:
    ConfigManager();
    ConfigManager(ConfigManager&) = delete;
    const ConfigManager& operator=(const ConfigManager&) = delete;

private:
    LogLevel logLevel;
    LogAppender logAppender;
    std::string logFile;
};


