/*
 * pbgz_file_handler.h - Header file for pbgz_v2 project
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

#ifndef PBGZ_FILE_HANDLER_H
#define PBGZ_FILE_HANDLER_H 

#include <stdint.h>
#include <string>

#include "pbgz_file.h"

uint8_t* getFileReadBuffer();

/// @brief pbgz格式文件的读取器
class PbgzFileReader {
public:
    int32_t open();

    const std::map<int32_t, PbgzFileHeader>& getFileHeader() const;

    const std::map<int32_t, PbgzFileMeta>& getFileMeta() const;

    int32_t readDataBlock(PbgzDataBlock& dataBlock);

    PbgzFileReader(const std::string& fileName) : fileName(fileName) {
        currentFileIndex = -1;  // 初始化文件索引为-1，表示未读取任何文件
    }

    void close();

    virtual ~PbgzFileReader();

private:
    int32_t initFileHeadAndMeta(bool isCheckMagic = false);

    static uint8_t* getFileReadBuffer();

private:
    std::string fileName;
    std::map<int32_t, PbgzFileHeader> fileHeaderMap;  // 文件头, 一个文件可能多个压缩包拼接而成
    std::map<int32_t, PbgzFileMeta> fileMetaMap;      // 文件元信息
    int32_t currentFileIndex; // 当前文件序号，表示当前读到哪个文件了
    FILE *pFile;
};

/// @brief  pbgz格式文件的写入器
/// Note: 
class PbgzFileWriter {
public:
    uint32_t write();

    PbgzFileWriter(const std::string& fileName) : fileName(fileName), pFile(nullptr) {}

private:
    static uint8_t* getFileWriteBuffer();


private:
    std::string fileName;
    FILE *pFile; // 文件指针
    
};

#endif