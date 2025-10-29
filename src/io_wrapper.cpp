
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>

#include "io_wrapper.h"
#include "log/logger.h"



int FileOperator::openIO() {
    if (fileName.empty()) {
        LOG_ERROR("IO Open failed, No file name gived");
        return -1;
    }

    if (fd != -1) {
        close(fd);
        fd = -1;
    }

    fd = open(fileName.c_str(), openMode, 0644);
    if (fd == -1) {
        LOG_ERROR("IO Open failed, file name = %s", fileName.c_str());
        return -1;
    }

    struct stat fileStat;
    if (fstat(fd, &fileStat) == -1) {
        LOG_ERROR("Faile to get file stat");
        close(fd);
    }

    fileSize = fileStat.st_size;

    return 0;
}

int32_t FileOperator::mmapFile(size_t mapFileSize) {
    mappedAddress = static_cast<uint8_t*>(mmap(nullptr, mapFileSize, mapMode, MAP_SHARED, fd, 0));
    if (mappedAddress == MAP_FAILED) {
        LOG_ERROR("file mmaped failed.file name = %s, errno = %d", fileName.c_str(), errno);
        close(fd);
        return -1;
    }

    return 0;
}

void FileOperator::closeIO() {
    if (mappedAddress != nullptr) {
        munmap(mappedAddress, fileSize);
        mappedAddress = nullptr;
    }

    if (fd != -1) {
        close(fd);
        fd = -1;
    }

    return;
}

int32_t FileOperator::seekIO(int32_t seekOffset, IOWhence whence) {
    switch (whence) {
    case IOWhence::IO_WHENCE_SET: {
        if (seekOffset < 0) {
            return -1;
        }

        position = seekOffset;
        break;
    }
    case IOWhence::IO_WHENCE_CUR: {
        size_t newPos = position + seekOffset;
        if (newPos > fileSize) {
            return -1;
        }

        position = newPos;
        break;
    }

    case IOWhence::IO_WHENCE_END: {
        if (seekOffset > 0) {
            return -1;
        }
        position = fileSize + seekOffset;
        break;
    }

    default:
        return -1;
    }

    return 0;
}

int32_t FileReader::openIO() {
    if (0 != fo.openIO()) {
        LOG_ERROR("IOOperator failed");
        return -1;
    }

    return fo.mmapFile(fo.fileSize);
}

size_t FileReader::readIO(void* pBuffer, size_t readSize) {
    if (pBuffer == nullptr) {
        return 0;
    }

    if (fo.mappedAddress == nullptr) {
        LOG_ERROR("File not mapped");
        return 0;
    }

    size_t left = fo.fileSize - fo.position;
    size_t realRead = readSize > left ? left : readSize; 
    memcpy(pBuffer, fo.mappedAddress + fo.position, realRead);
    
    fo.position += realRead;
    return realRead;
}

void* FileReader::getAt(size_t pos) {
    return fo.mappedAddress + pos;
}

int32_t FileWriter::openIO() {
    if (0 != fo.openIO()) {
        LOG_ERROR("File writer open failed.");
        return -1;
    }
    mapSize = getMappedSize(fo.fileSize);
    ftruncate(fo.fd, mapSize);
    if (0 != fo.mmapFile(mapSize)) {
        LOG_ERROR("File writer mmap failed.");
        return -1;
    }
    return 0;
}

void FileWriter::closeIO() {
        // Force flush data before exit
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
        ftruncate(fo.fd, fo.fileSize);
        return fo.closeIO();
    }

size_t FileWriter::getMappedSize(size_t fileSize) {
    const size_t mapSizeUint = 16 * 1024 * 1024;
    return ((fileSize / mapSizeUint) + 1) * mapSizeUint;
}

size_t FileWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (pBuffer == nullptr || writeLen == 0) {
        LOG_ERROR("Input invalid.");
        return 0;
    }

    if (fo.mappedAddress == nullptr) {
        LOG_ERROR("mappedAddress is nullptr.");
        return 0;
    }

    size_t newSize = fo.fileSize + writeLen;
    if (newSize > mapSize) {
        // If file size exceeds mapping size, need to remap
        // Flush content to disk first
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
        // Unmap
        munmap(fo.mappedAddress, mapSize);
        // Calculate new mapping size
        mapSize = getMappedSize(newSize);
        // Extend file
        ftruncate(fo.fd, mapSize);
        // Remap file
        if (0 != fo.mmapFile(mapSize)) {
            LOG_ERROR("File not mapped");
            return 0;
        }
    }
   
    (void)memcpy(fo.mappedAddress + fo.fileSize, pBuffer, writeLen);
    fo.fileSize = newSize;
    return newSize;
}

int32_t FileWriter::writeIOAt(int32_t seekOffset, const void* pBuffer, size_t writeLen) {
    if (pBuffer == nullptr || writeLen == 0) {
        LOG_ERROR("Input invalid.");
        return 0;
    }

    size_t newFileSize = fo.fileSize;
    if (seekOffset + writeLen > mapSize) {
         // If file size exceeds mapping size, need to remap
        // Flush content to disk first
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
        // Unmap
        munmap(fo.mappedAddress, mapSize);
        // New file size
        newFileSize = seekOffset + writeLen;
        mapSize = getMappedSize(newFileSize);
        ftruncate(fo.fd, mapSize);
        if (0 != fo.mmapFile(mapSize)) {
            LOG_ERROR("File not mapped");
            return 0;
        }
    }

    (void)memcpy(fo.mappedAddress + seekOffset, pBuffer, writeLen);
    fo.fileSize = newFileSize;
    return newFileSize;
}

size_t PipeReader::readIO(void* pBuffer, size_t readSize) {
    if (!pBuffer) {
        return -1;
    }

    return read(STDIN_FILENO, pBuffer, readSize);
}


size_t PipeWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (!pBuffer) {
        return -1;
    }

    return write(STDOUT_FILENO, pBuffer, writeLen);
}
