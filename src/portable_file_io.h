/*
 * portable_file_io.h - Portable streaming file I/O for pbgz
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

#include <cstdio>
#include <string>

#include "io_interface.h"

/*
 * Streaming file reader based on the C standard library.  It keeps only the
 * stdio buffer in memory, so resident memory does not grow with the input file
 * size (unlike the mmap FileReader).  Uses 64-bit offsets on every platform.
 */
class BufferedFileReader : public IOReader, public SeekableReader {
public:
    explicit BufferedFileReader(const std::string& fileName);
    ~BufferedFileReader() override;

    int32_t openIO() override;
    void closeIO() override;
    size_t readIO(void* pBuffer, size_t readSize) override;
    size_t readLine(std::string& line) override;
    size_t readLine(std::string& line, ChecksumStream* sink) override;

    int32_t seekIO(uint64_t seekOffset,
                   IOWhence whence = IOWhence::IO_WHENCE_SET) override;
    size_t getFileSize() override { return fileSize; }
    size_t getCurrentPos() override { return position; }

private:
    std::string fileName;
    FILE* fp = nullptr;
    size_t fileSize = 0;
    size_t position = 0;
};

/*
 * Streaming file writer based on the C standard library.  Sequential writes
 * append at the logical end of the file; writeIOAt() patches an earlier offset.
 * Memory usage is bounded by the stdio buffer.
 */
class BufferedFileWriter : public IOWriter, public SeekableWriter {
public:
    explicit BufferedFileWriter(const std::string& fileName);
    ~BufferedFileWriter() override;

    int32_t openIO() override;
    void closeIO() override;
    size_t writeIO(const void* pBuffer, size_t writeLen) override;

    int32_t seekIO(uint64_t seekOffset,
                   IOWhence whence = IOWhence::IO_WHENCE_SET) override;
    int32_t seekToEnd() override;
    int32_t writeIOAt(size_t seekOffset, const void* pBuffer, size_t writeLen) override;
    void flushIO() override;
    size_t getCurrentPos() override { return position; }
    size_t getFileSize() override { return fileSize; }

private:
    std::string fileName;
    FILE* fp = nullptr;
    size_t fileSize = 0;
    size_t position = 0;
};
