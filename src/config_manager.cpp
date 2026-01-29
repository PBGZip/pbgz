/*
 * config_manager.cpp - Source file for pbgz project
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



