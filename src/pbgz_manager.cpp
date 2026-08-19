/*
 * pbgz_manager.cpp - Source file for pbgz project
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


#include <unistd.h>
#include <cstring>

#include "pbgz_manager.h"
#include "pbgz_errno.h"
#include "utils/path_util.h"
#include "block_wrapper.h"


PbgzManager& PbgzManager::getInstance() {
    static PbgzManager instance;
    return instance;
}

void PbgzManager::exitProc(int errorCode, const char* errorMessage){
    if (errorCode < pbgz::PBGZ_ERR_OK) {
        for (auto &currFile : outfiles) {
            const std::string& fileName = currFile.first;
            if (!currFile.second && PathUtil::fileExists(fileName)) {
                unlink(fileName.c_str());
            }
        }
    } 

    if (errorMessage != nullptr && strlen(errorMessage) > 0) {
        fprintf(stdout, "pbgz exit : %s \n", errorMessage);
    }

    _Exit(errorCode);
}

void PbgzManager::addOutputFile(const std::string& fileName) {
    /* second 置 false：异常退出（errorCode < 0）时删除该输出文件。 */
    outfiles.push_back(std::make_pair(fileName, false));
}

void pbgzExitProc(int errorCode, const char* errorMessage) {
    return PbgzManager::getInstance().exitProc(errorCode, errorMessage);
}

void coderLog(int logLevel, const char* logMessage) {
    LogLevel pbgzLogLevel = LogLevel::OFF;
    switch (logLevel)
    {
    case coder_ns::DEBUGGING:
        pbgzLogLevel = LogLevel::DEBUGGING;
        break;
    case coder_ns::INFO:
        pbgzLogLevel = LogLevel::INFO;
        break;
    case coder_ns::WARNING:
        pbgzLogLevel = LogLevel::WARNING;
        break;
    case coder_ns::ERROR:
        pbgzLogLevel = LogLevel::ERROR;
        break;
    case coder_ns::FATAL:
        pbgzLogLevel = LogLevel::FATAL;
        break;    
    default:
        break;
    }

    return Logger::getInstance().log(pbgzLogLevel, logMessage);
}

std::string PbgzManager::getVersion() {
    char buffer[12];
    snprintf(buffer, sizeof(buffer), "%d.%d.%d", PBGZ_VERSION_MAJOR, PBGZ_VERSION_MINOR, PBGZ_VERSION_PATCH);
    return std::string(buffer);
}

std::vector<char> PbgzManager::getVersionAsArray() {
    return std::vector<char>{PBGZ_VERSION_MAJOR, PBGZ_VERSION_MINOR, PBGZ_VERSION_PATCH};
}

void PbgzManager::updateReadDataLen(RoughIOBlock* blockPtr) {
    totalReadLen += blockPtr->getDataLen();
    updateDataInfo();
}
    
void PbgzManager::updateWriteDataLen(RoughIOBlock* blockPtr) {
    totalWriteLen += blockPtr->getDataLen() + blockPtr->getMetaLen();
    updateDataInfo();
}   

void PbgzManager::updateDataInfo() {
    /*
     * 进度条原地刷新：先清行（\033[K）再写内容，末尾 \r 回到行首。
     * 不用清行的话，较长的上一次内容会残留，且与其他 stderr 输出（如 File type）
     * 交错时会把对方末尾盖在行尾，出现 "File type: SAM40/0]---" 这类乱码。
     */
    fprintf(stderr, "\r\033[K\033[37mfrom/to ---[%ld/%ld]---\r", totalReadLen, totalWriteLen);
}

void PbgzManager::printFileType(BlockType blockType) {
    /* 可能插在进度条之后打印：先回行首并清行，避免覆盖进度条残留 */
    fprintf(stderr, "\r\033[KFile type: %s\n", BlockUtil::getBlockTypeName(blockType).c_str());
}

void PbgzManager::printHeadInfo(PbgzParameter& para) {
    fprintf(stderr, "\033[37mpbgz version => %s\033[0m\n\n", getVersion().c_str());
    fprintf(stderr, "Parallel set: %d\n", para.threadNum);
}

void PbgzManager::printTailInfo(Timer costTime, bool isPrintRatio) {
    fprintf(stderr, "\r\033[K\033[37mfrom/to ---[%ld/%ld]---\033[0m\n", totalReadLen, totalWriteLen);
    if (isPrintRatio) {
        fprintf(stderr, "\nCompress finish, cost %um%us.\n", (costTime.elapsed() / 1000) / 60, (costTime.elapsed() / 1000) % 60);
        fprintf(stderr, "Total size_dest size %lu bytes, compressed to %lu bytes, ratio %0.2f%%\n", totalReadLen, totalWriteLen, (totalWriteLen * 1.0) * 100 / totalReadLen);
    } else {
        fprintf(stderr, "\nDecompress finish, cost %um%us.\n", (costTime.elapsed() / 1000) / 60, (costTime.elapsed() / 1000) % 60);
    }
}

/* Convert n to the nearest integer power of 2 */
int32_t powerof2Proximal(int32_t i)
{
    i |= (i >> 1);
    i |= (i >> 2);
    i |= (i >> 4);
    i |= (i >> 8);
    i |= (i >> 16);
    return i - (i >> 1);
}
