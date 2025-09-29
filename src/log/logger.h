/*
 * logger.h - Header file for pbgz_v2 project
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

#ifndef LOGGER_H
#define LOGGER_H

#include <stdint.h>
#include <string>

enum class LogLevel{
    TRACE = 0,     // 跟踪日志
    DEBUGGING = 1, // 调试日志
    INFO = 2,      // 重要的信息
    WARNING = 3,   // 告警日志
    ERROR = 4,     // 错误日志
    FATAL = 5,     // 致命错误
    OFF = 6        // 关闭所有日志
 };

 enum class LogAppender{
    CONSOLE = 0,  // 输出到控制台
    FILE = 1,     // 输出到文件
    NETWORK = 2,  // 输出到网络
 };


class Logger {
public:
    void log(LogLevel logLevel, int line, const char* function, const char* fileName, 
                const char* logFormat, ...);

    static Logger& getInstance();

    void log(LogLevel logLevel, const char* logMessage);

private:
    std::string getLogLevelString(LogLevel logLevel);

    int32_t buildLogString(char* logBuffer, uint32_t bufferLen, LogLevel logLevel, int line, const char* function, const char* fileName, 
                const char* logContent);
    
    Logger():logFileName("pbgz.log"){}

private:
    std::string logFileName;
};




#define LOG_TRACE(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::TRACE, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 

#define LOG_DEBUG(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::DEBUGGING, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 

#define LOG_INFO(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::INFO, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 

#define LOG_WARNING(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::WARNING, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 

#define LOG_ERROR(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::ERROR, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 

#define LOG_FATAL(logFormat, ... )  \
    Logger::getInstance().log(LogLevel::FATAL, __LINE__, __FUNCTION__, __FILE__, logFormat, ##__VA_ARGS__) 


#endif

