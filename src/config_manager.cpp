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

ConfigManager::ConfigManager() {
    // 读取配置文件 pbgz_conf.json
    std::ifstream configFile("pbgz_conf.json");
    if (!configFile.is_open()) {
        // 当前目录没有配置文件则到根目录获取
        struct passwd *pw = getpwuid(getuid());
        if (pw != NULL) {
            char separator = std::filesystem::path::preferred_separator;
            std::string configName =  std::string(pw->pw_dir) + separator + "pbgz_conf.json";
            configFile.open(configName);
        }
    }

    if (configFile.is_open()) {
        Json::Value root;
        Json::CharReaderBuilder builder;
        std::string errs;
        if (Json::parseFromStream(builder, configFile, &root, &errs)) {
            if (root.isMember("log")) {
                Json::Value logConfig = root["log"];
                // 读取日志等级
                if (logConfig.isMember("level")) {
                    std::string levelStr = logConfig["level"].asString();
                    if (levelStr == "TRACE") logLevel = LogLevel::TRACE;
                    else if (levelStr == "DEBUG") logLevel = LogLevel::DEBUGGING;
                    else if (levelStr == "INFO") logLevel = LogLevel::INFO;
                    else if (levelStr == "WARNING") logLevel = LogLevel::WARNING;
                    else if (levelStr == "ERROR") logLevel = LogLevel::ERROR;
                    else if (levelStr == "FATAL") logLevel = LogLevel::FATAL;
                    else logLevel = LogLevel::OFF;
                }
                // 读取日志输出方式
                if (logConfig.isMember("appender-ref")) {
                    std::string appenderName = logConfig["appender-ref"].asString();
                    if (logConfig.isMember("appender")) {
                        Json::Value appenderList = logConfig["appender"];
                        if (appenderList.isMember(appenderName)){
                            Json::Value appender = appenderList[appenderName];
                            std::string type = appender["type"].asString();
                            if (type == "CONSOLE") {
                                logAppender = LogAppender::CONSOLE;
                            }
                            else if (type == "FILE") {
                                logAppender = LogAppender::FILE;
                                logFile = appender["file_name"].asString();
                                maxLogFileSize =  appender["max_size"].asUInt();
                            } else if (type == "NETWORK") {
                                logIpAddress = appender["ip"].asString();
                                logPort = appender["port"].asUInt();
                                logUrl = appender["url"].asString();
                            }else{
                                // 默认输出到控制台
                                logAppender = LogAppender::CONSOLE;
                            }
                        }
                    }
                }
            }
        } else {
            fprintf(stderr, "load config file err: %s \n", errs.c_str());
        }
        configFile.close();
    } else {
       // 没有配置文件使用默认的配置
       logLevel = LogLevel::OFF;
       logAppender = LogAppender::CONSOLE;
    }
}


