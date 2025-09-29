/*
 * logger.cpp - CPP file for pbgz_v2 project
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

#include <cstdarg>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <time.h>

#include "logger.h"
#include "config_manager.h"


Logger& Logger::getInstance() {
    static Logger logInstance;
    return logInstance;
}

std::string Logger::getLogLevelString(LogLevel logLevel) {
    switch (logLevel) {
        case LogLevel::TRACE: return "TRACE";
        case LogLevel::DEBUGGING: return "DEBUG";
        case LogLevel::INFO: return "INFO";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::FATAL: return "FATAL";
        case LogLevel::OFF: return "OFF";
        default: return "UNKNOWN";
    }
}

int32_t Logger::buildLogString(char* logBuffer, uint32_t bufferLen, LogLevel logLevel, int line, const char* function, const char* fileName, 
                const char* logContent) {
    time_t now = time(nullptr);
    struct tm *ltm = localtime(&now);
    char timeBuffer[32];
    strftime(timeBuffer, sizeof(timeBuffer), "%Y-%m-%d %H:%M:%S", ltm); 

    std::string logLevelStr = getLogLevelString(logLevel);

    /// 不做长度判断，如果空间不够，则会截断
    snprintf(logBuffer, bufferLen, "[%s][%s][%s]%s[%s:%d]\n",
         getLogLevelString(logLevel).c_str(), function, timeBuffer, logContent, fileName, line);

    return 0;
}

void Logger::log(LogLevel logLevel, int line, const char* function, const char* fileName, 
                const char* logFormat, ...) {
    LogLevel configLevel = ConfigManager::getInstance().getLogLevel();
    if (configLevel == LogLevel::OFF) {
        return;
    }
    if (configLevel > logLevel) {
        return;
    }

    char logContent[512];
    va_list args;
    va_start(args, logFormat);
    vsnprintf(logContent, sizeof(logContent), logFormat, args);
    va_end(args);

    char logBuffer[1024];
    buildLogString(logBuffer, sizeof(logBuffer), logLevel, line, function, fileName, logContent);
    
    return log(configLevel, logBuffer);
}

void Logger::log(LogLevel logLevel, const char* logMessage) {
    LogLevel configLevel = ConfigManager::getInstance().getLogLevel();
    if (configLevel == LogLevel::OFF) {
        return;
    } 

    if (configLevel > logLevel) {
        return;
    }

    // 先简单实现，后续可优化为异步写文件方式，防止IO阻塞主进程
    LogAppender appender = ConfigManager::getInstance().getLogAppender();
    switch (appender) {
        case LogAppender::CONSOLE:{
            if (logLevel > LogLevel::ERROR) {
                std::cerr << logMessage;
            } else {
                std::cout << logMessage;
            }
            break;
        }
        case LogAppender::FILE:{
            FILE* pLogFile = fopen(logFileName.c_str(), "a");
            if (pLogFile != nullptr) {
                fwrite(logMessage, 1, strlen(logMessage), pLogFile);
                fflush(pLogFile);
                fclose(pLogFile);
                pLogFile = nullptr;
            }
            break;
        }
        case LogAppender::NETWORK:  // 暂不支持输出到网络
        default:
            break;
    }   
}


