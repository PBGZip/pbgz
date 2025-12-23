#include <fstream>
#include <json/json.h>
#include <iostream>
#include <unistd.h>
#include <pwd.h>
#include <filesystem>

#include "config_manager.h"


ConfigManager& ConfigManager::getInstance() {
    static ConfigManager instance;
    return instance;
}

LogLevel ConfigManager::getLogLevel() {
    return logLevel;
}

LogAppender ConfigManager::getLogAppender(){
    return logAppender;
}

std::string& ConfigManager::getLogFileName() {
    return logFile;
}

ConfigManager::ConfigManager() {
    logLevel = LogLevel::OFF;
    logAppender = LogAppender::CONSOLE;
}

void ConfigManager::init(PbgzParameter& para) {
    switch (para.logLevel) {
    case 1:
        logLevel = LogLevel::DEBUGGING;
        break;
    case 2:
        logLevel = LogLevel::INFO;
        break;
    case 3:
        logLevel = LogLevel::WARNING;
        break;
    case 4:
        logLevel = LogLevel::ERROR;
        break;
    case 5:
        logLevel = LogLevel::FATAL;
        break;
    case 6:
    default:
        logLevel = LogLevel::OFF;
        break;
    }
    
    if (logFile.empty()) {
        logAppender = LogAppender::CONSOLE;
    } else {
        logAppender = LogAppender::FILE;
        logFile = para.logFile;
    }
}



