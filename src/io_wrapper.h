/*
 * io_wrapper.h - Header file for pbgz project
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

#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdint.h>
#include <htslib/bgzf.h>
#include <isa-l/igzip_lib.h>
#include <zlib.h>

enum class IOWhence {
    IO_WHENCE_SET = 0,
    IO_WHENCE_CUR = 1,
    IO_WHENCE_END = 2, 
};

class IOReader {
public:
    virtual int32_t openIO() = 0;

    virtual void closeIO() = 0;

    virtual size_t readIO(void* pBuffer, size_t readSize) = 0;

    virtual ~IOReader() {};
};

class IOWriter {
public:
    virtual int32_t openIO() = 0;

    virtual void closeIO() = 0;

    virtual size_t writeIO(const void* pBuffer, size_t writeLen) = 0;

    virtual ~IOWriter() {}
};

class FileOperator {
public:
    friend class FileReader;
    friend class FileWriter;
public:
    int32_t openIO();

    void closeIO();

    FileOperator(const std::string& inFileName) : fileName(inFileName) {
        mappedAddress = nullptr;
        fileSize = 0;
        fd = -1;
        position = 0;
    }

    virtual ~FileOperator() {
        closeIO();
    } 

    int32_t seekIO(size_t seekOffset, IOWhence whence);

protected:
    int32_t mmapFile(size_t mapFileSize);

protected:
    std::string fileName;
    int fd;
    uint8_t* mappedAddress;
    size_t fileSize;
    size_t position;
    int32_t mapMode;
    int32_t openMode;
};

class FileReader : public IOReader {

public:
    int32_t openIO();

    void closeIO() {return fo.closeIO();}

    size_t readIO(void* pBuffer, size_t readSize);

    FileReader(const std::string& inFileName): fo(inFileName) {
        fo.openMode = O_RDONLY;
        fo.mapMode = PROT_READ;
    }

    void* getAt(size_t pos);

    int32_t seekIO(uint64_t seekOffset, IOWhence whence = IOWhence::IO_WHENCE_SET) {
        return fo.seekIO(seekOffset, whence);
    }

    size_t getFileSize() {
        return fo.fileSize;
    }

    size_t getCurrentPos() {
        return fo.position;
    }

private:
    FileOperator fo;
};

class FileWriter : public IOWriter {
public:
    int32_t openIO();

    void closeIO();

    size_t writeIO(const void* pBuffer, size_t writeLen);

    FileWriter(const std::string& inFileName): fo(inFileName)  {
        fo.openMode = O_RDWR | O_CREAT;
        fo.mapMode = PROT_READ|PROT_WRITE;
        mapSize = 0;
    }

    int32_t seekIO(size_t seekOffset, IOWhence whence = IOWhence::IO_WHENCE_SET) {
        return fo.seekIO( seekOffset, whence);
    }

    int32_t writeIOAt(size_t seekOffset, const void* pBuffer, size_t writeLen);

    void flushIO() {
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
    }

    ~FileWriter() {
    }

    size_t getMappedSize(size_t fileSize);

    size_t getCurrentPos() {
        return fo.position;
    }

protected:
    FileOperator fo;
    size_t mapSize;
};

class PipeReader : public IOReader {
public:
    PipeReader() { }

    // Pipe reading does not require open and close operations
    int openIO() { return 0; }
    void closeIO() { return; }

    size_t readIO(void* pBuffer, size_t readSize);
};


class PipeWriter : public IOWriter {
public:
    PipeWriter() { }

    // Pipe writing does not require open and close operations
    int openIO() { return 0; }
    void closeIO() { return; }

    size_t writeIO(const void* pBuffer, size_t writeLen);
};

class GzFileReader : public IOReader {
public:
    GzFileReader(const std::string& fileName, uint32_t parall) {
        gzFileName = fileName;
        parallel = parall;
    }

    int openIO() override;

    size_t readIO(void* pBuffer, size_t readSize) override;

    void closeIO() override;

    ~GzFileReader() {
        closeIO();
    }
private:
    std::string gzFileName;
    BGZF* fpGz;
    uint32_t parallel;
};

class GzFileWriter : public IOWriter {
public:
    GzFileWriter(std::string& fileName, uint32_t parall) {
        gzFileName = fileName;
        parallel = parall;
    }

    int openIO() override;

    size_t writeIO(const void* pBuffer, size_t writeLen) override;

    void closeIO() override;

    ~GzFileWriter() {
        closeIO();
    }

private:
    std::string gzFileName;
    BGZF* fpGz;
    uint32_t parallel;
};

class FastGzFileReader : public IOReader {
public:
    FastGzFileReader(std::string& fileName) : fileReader(fileName) {
        gzFileName = fileName;
    }

    int openIO() override;

    size_t readIO(void* pBuffer, size_t readSize) override;

    void closeIO() override;

    ~FastGzFileReader() {
        closeIO();
    }

private:
    std::string gzFileName;
    FileReader fileReader;

    struct inflate_state isalState;
    struct isal_gzip_header isalHeader;
    size_t isalExtraEach;
    uint8_t* isalExtraBuffer;
    uint32_t isalExtraLen;
    uint8_t* isalInputBuffer;
    int64_t isalInputLength;
    int64_t firstReadLen;

    // Helper function
    bool readMoreDataAndReset(uint8_t* isalIn, uint32_t isalInLen, uint8_t* isalOut, uint32_t isalOutLen, bool isMultipleHeaders);
    int64_t processMultipleGzipHeaders();
};

class GzPipeReader : public IOReader {
public:
    GzPipeReader(uint32_t paral);

    int openIO() override;

    size_t readIO(void* pBuffer, size_t readSize) override;

    void closeIO() override;
    
    ~GzPipeReader() {
        closeIO();
    }
private:
    BGZF* fpGz;
    uint32_t parallel;
};

class GzPipeWriter : public IOWriter {
public:
    GzPipeWriter(uint32_t paral);

    int openIO() override;

    size_t writeIO(const void* pBuffer, size_t writeLen) override;

    void closeIO() override;
    
    ~GzPipeWriter() {
        closeIO();
    }
private:
    BGZF* fpGz;
    uint32_t parallel;
};

class FastGzPipeReader : public IOReader {
public:
    FastGzPipeReader();

    int openIO() override;
    size_t readIO(void* pBuffer, size_t readSize) override;
    void closeIO() override;

    ~FastGzPipeReader() {
        closeIO();
    }
private:
    PipeReader pipeReader;
    
    // ISAL decompression state
    struct inflate_state state;
    bool stateInitialized;
    
    // ISAL gzip header
    struct isal_gzip_header gzipHeader;
    
    // Buffers
    uint8_t* inputBuffer;
    uint8_t* outputBuffer;
    uint8_t* remainingBuffer;
    size_t inputBufferSize;
    size_t outputBufferSize;
    size_t remainingBufferSize;
    size_t remainingSize;
    
    // Status flags
    bool streamEnded;
};
