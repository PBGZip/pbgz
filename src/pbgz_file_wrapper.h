/*
 * pbgz_file_handler.h - Header file for pbgz project
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

#pragma once

#include <stdint.h>
#include <string>

#include "pbgz_file.h"
#include "io_wrapper.h"

class RoughIOBlock;


/// @brief Reader for pbgz format files
class PbgzFileReader {
public:
    int32_t open();

    PbgzFileHeader& getFileHeader();

    PbgzFileMeta& getBaseFileMeta();

    PbgzFileMeta& getDynamicFileMeta();

    /*
     * 当前正在读的那个 pbgz 包在整个输入流里的起始字节位置。
     *
     * 存在的理由是 cat 拼接：多个独立压缩包首尾相接成一个文件之后，从第二个包开始，
     * 包内记录的所有偏移都不再等于它在拼接文件里的实际位置。压缩时输出总是从 0 开始
     * 写，所以包里存的偏移本质上是"相对本包头部的距离"；单包读取时包起点恰好是 0，
     * 相对值与绝对值相等，问题被掩盖，只有拼接才会暴露。
     *
     * 读侧拿到这个起点，用「包起点 + 包内偏移」还原真实位置即可，写出的格式一个字节
     * 都不用改，也不影响已有文件的兼容性。
     */
    uint64_t getCurrentFileStart() const;

    /*
     * 最近一次 readDataBlock 所读那个块的绝对起点（魔数所在字节）。
     *
     * 辅助块的身份就是它的绝对地址：包起点只在"每个包恰好一份"时够用，先验一旦分片
     * 就不成立。缓存按这个地址做键，重复遇到同一个块时可以直接命中而不必重读。
     */
    uint64_t getCurrentBlockStart() const { return currentBlockStart; }

    /*
     * 当前正在读的是第几个 pbgz 包（cat 拼接后依次递增）。
     * 与 getCurrentFileStart 不同，这个序号在管道输入下同样有效——
     * 它只依赖已解析过的包头个数，不依赖任何文件位置。
     */
    int32_t getCurrentFileIndex() const { return currentFileIndex; }

    /*
     * 读一个数据块。dst 非空时，块数据先按 dataLength 扩容到 dst（输入块缓冲按
     * 参数级 block_size 预分配，而文件的块可能更大，如 -l 9 的 512MB 块压缩后仍达
     * ~72MB），再把数据读进 dst 的缓冲。
     */
    int32_t readDataBlock(PbgzDataBlock& dataBlock, RoughIOBlock* dst = nullptr);

    PbgzFileReader(IOReader* pReader) : ioReader(pReader) {
        currentFileIndex = -1;  // Initialize file index to -1, indicating no file has been read yet
        currentBlockStart = 0;
    }

    void close();

    virtual ~PbgzFileReader() { }

private:
    uint64_t currentBlockStart;

    int32_t initFileHeadAndMeta(bool isCheckMagic = false);

    static uint8_t* getFileReadBuffer();

    int32_t readFileMeta(PbgzFileMeta& fileMeta, bool isCheckMagic = true);

private:
    std::map<int32_t, PbgzFileHeader> fileHeaderMap;  // File headers, a file may consist of multiple compressed packages concatenated together
    std::map<int32_t, PbgzFileMeta> baseFileMetaMap;      // Base file metadata
    std::map<int32_t, PbgzFileMeta> dynamicFileMetaMap;   // Dynamic file meta data, only exits in file reader
    std::map<int32_t, uint64_t> fileStartMap;             // 每个拼接包在输入流里的起始位置，用于把包内相对偏移还原成真实位置
    int32_t currentFileIndex; // Current file sequence number, indicating which file is currently being read
    IOReader* ioReader;
};

/// @brief Writer for pbgz format files
/// @note: 
class PbgzFileWriter {
public:
    int32_t open();

    int32_t close();

    int32_t writeFileMeta(PbgzFileMeta& fileMeta);

    int32_t writeBaseFileMeta() {
        return writeFileMeta(baseFileMeta);
    }

    int32_t writeDynamicFileMeta() {
        return writeFileMeta(dynamicFileMeta);
    }

    int32_t writeBlockData(PbgzDataBlock& dataBlock);

    int32_t setBaseFileMeta(const PbgzFileMeta& metaInfo) {
        baseFileMeta = metaInfo;
        return 0;
    }

     int32_t setDynamicFileMeta(const PbgzFileMeta& metaInfo) {
        dynamicFileMeta = metaInfo;
        return 0;
    }

    PbgzFileWriter(IOWriter* pWriter ) : ioWriter(pWriter) {}

    virtual ~PbgzFileWriter() { 
        ioWriter = nullptr;
    }

    PbgzFileHeader& getFileHeader() {
        return fileHeader;
    }

    PbgzFileMeta& getBaseFileMeta() {
        return baseFileMeta;
    }

    PbgzFileMeta& getDynamicFileMeta() {
        return dynamicFileMeta;
    }

    void updateMetaOffset(uint64_t dynamicMetaOffset);

private:
    static uint8_t* getFileWriteBuffer();

    int32_t initFileHead();

private:
    PbgzFileHeader fileHeader; // File header
    PbgzFileMeta baseFileMeta; // Base file meta
    PbgzFileMeta dynamicFileMeta; // Dynamic file meta 
    IOWriter* ioWriter;
};
