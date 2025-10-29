
#pragma once

#include <string>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdint.h>

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

    int32_t seekIO(int32_t seekOffset, IOWhence whence);

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

    int32_t seekIO(int32_t seekOffset, IOWhence whence = IOWhence::IO_WHENCE_SET) {
        return fo.seekIO(seekOffset, whence);
    }

    size_t getFileSize() {
        return fo.fileSize;
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

    int32_t seekIO(int32_t seekOffset, IOWhence whence = IOWhence::IO_WHENCE_SET) {
        return fo.seekIO( seekOffset, whence);
    }

    int32_t writeIOAt(int32_t seekOffset, const void* pBuffer, size_t writeLen);

    void flushIO() {
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
    }

    ~FileWriter() {
    }

    size_t getMappedSize(size_t fileSize);

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