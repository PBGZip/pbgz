
#include <unistd.h>
#include <sys/stat.h>
#include <string.h>

#include "io_wrapper.h"
#include "log/logger.h"
#include "utils/memory_util.h"


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

int32_t FileOperator::seekIO(size_t seekOffset, IOWhence whence) {
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
    if (fo.mappedAddress != nullptr) {
        msync(fo.mappedAddress, fo.fileSize, MS_SYNC);
        // Truncate file to actual size before unmapping
        if (fo.fd != -1) {
            ftruncate(fo.fd, fo.fileSize);
        }
        // Unmap with the actual mapped size, not fileSize
        munmap(fo.mappedAddress, mapSize);
        fo.mappedAddress = nullptr;
    }
    
    // Close file descriptor
    if (fo.fd != -1) {
        close(fo.fd);
        fo.fd = -1;
    }
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
    fo.position = fo.fileSize;
    return newSize;
}

int32_t FileWriter::writeIOAt(size_t seekOffset, const void* pBuffer, size_t writeLen) {
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
    fo.position = seekOffset + writeLen;
    return newFileSize;
}

size_t PipeReader::readIO(void* pBuffer, size_t readSize) {
    if (!pBuffer) {
        return -1;
    }

    size_t totalRead = 0;
    uint8_t* buffer = static_cast<uint8_t*>(pBuffer);
    
    // 循环读取直到达到请求的大小或没有更多数据
    while (totalRead < readSize) {
        ssize_t bytesRead = read(STDIN_FILENO, buffer + totalRead, readSize - totalRead);
        if (bytesRead < 0) {
            if (errno == EINTR) {
                continue; // 被信号中断，重试
            } else if (errno == EAGAIN || errno == EWOULDBLOCK) {
                break; // 非阻塞模式下无数据可读
            } else {
                return -1; // 其他错误
            }
        } else if (bytesRead == 0) {
            break; // EOF，写入端关闭
        }
        
        totalRead += bytesRead;
    }
    
    return totalRead;
}

size_t PipeWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (!pBuffer) {
        return -1;
    }

    return write(STDOUT_FILENO, pBuffer, writeLen);
}

int GzFileReader::openIO() {
    fpGz = bgzf_open(gzFileName.c_str(), "rb");
    if (fpGz == nullptr) {
        LOG_ERROR("bgzf open failed.");
        return -1;
    }

    if (fpGz->is_gzip) {
        bgzf_mt(fpGz, parallel, 256);
    }

    return 0;
}

size_t GzFileReader::readIO(void* pBuffer, size_t readSize) {
    ssize_t readLen = bgzf_read(fpGz, pBuffer, readSize);
    if (readLen < 0) {
        LOG_ERROR("Read file from %s failed", gzFileName.c_str());
        return 0;
    }

    return (size_t)readLen;
}

void GzFileReader::closeIO() {
    if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    }
}

int GzFileWriter::openIO() {
    fpGz = bgzf_open(gzFileName.c_str(), "wb");
    if (fpGz == nullptr) {
        LOG_ERROR("bgzf open failed.");
        return -1;
    }

    fpGz->compress_level = 6;
    bgzf_mt(fpGz, parallel, 256);
    return 0;
}

size_t GzFileWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (pBuffer == nullptr || writeLen == 0) {
        return 0;
    }

    size_t len = bgzf_write(fpGz, pBuffer, writeLen);
    if (len != writeLen) {
        LOG_ERROR("Write gz file error, write length = %d, expect %d", len, writeLen);
        return -1;
    }
    return len;
}

void GzFileWriter::closeIO() {
    if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    }
}

int FastGzFileReader::openIO() {
    fileReader.openIO();
    isalExtraEach = 2048;
    isalExtraLen = 2 << 20;
    isalExtraBuffer = MemoryUtil::safeAlloc<uint8_t>(isalExtraLen);
    if (isalExtraBuffer == nullptr) {
        return -1;
    }

    isalInputLength = 128 << 20;
    isalInputBuffer = MemoryUtil::safeAlloc<uint8_t>(isalInputLength);
    if (isalInputBuffer == nullptr) {
        MemoryUtil::safeFree(isalExtraBuffer);
        return -1;
    }

    isal_gzip_header_init(&isalHeader);
    isal_inflate_init(&isalState);
    isalState.crc_flag = ISAL_GZIP_NO_HDR_VER;
    isalState.next_in = isalInputBuffer;

    int64_t readLen = std::max<int64_t>(64<<20, isalInputLength);
    isalState.avail_in = fileReader.readIO(isalState.next_in, readLen);
    firstReadLen = isalState.avail_in;
    if (0 != isal_read_gzip_header(&isalState, &isalHeader)) {
        LOG_ERROR("found in valid gzip header in file %s.", gzFileName.c_str());
        MemoryUtil::safeFree(isalExtraBuffer);
        MemoryUtil::safeFree(isalInputBuffer);
        return -1;
    }
    return 0;
}

size_t FastGzFileReader::readIO(void* pBuffer, size_t readSize) {
    uint64_t unzippedLen = 0;
    uint64_t isalRes = ISAL_DECOMP_OK;
    isalState.next_out = static_cast<uint8_t*>(pBuffer);
    isalState.avail_out = readSize;
    
    for (;;) {
        if (unzippedLen >= readSize || isalState.next_in == nullptr) {
            break;
        }  

        for (;;) {
            uint8_t* isalIn = isalState.next_in;
            uint32_t isalInLen = isalState.avail_in;
            uint8_t* isalOut = isalState.next_out;
            uint32_t isalOutLen = isalState.avail_out;

            isalRes = isal_inflate(&isalState);
            if (isalRes == ISAL_DECOMP_OK) {
                break;
            }

            if (isalRes != ISAL_END_INPUT) {
                LOG_ERROR("Invalid gz file. errno = %d.", isalRes);
                return -1;
            }

            // 读取更多数据并重置状态
            if (!readMoreDataAndReset(isalIn, isalInLen, isalOut, isalOutLen, false)) {
                return -1;
            }
        }

        unzippedLen = isalState.next_out - static_cast<uint8_t*>(pBuffer);
        if (unzippedLen != readSize) {
            bool fileDone = true;
            if (isalState.avail_in == 0) {
                isalState.next_in = isalInputBuffer;
                int64_t readLen = isalInputLength;
                isalState.avail_in = fileReader.readIO(isalState.next_in, readLen);
                fileDone = isalState.avail_in == 0;
            } else {
                // 处理多gzip头的情况
                isalRes = processMultipleGzipHeaders();
                if (isalRes != ISAL_DECOMP_OK) {
                    return -1;
                }
                fileDone = isalRes != ISAL_DECOMP_OK;
            }

            if (fileDone) {
                isalState.next_in = nullptr;
                return unzippedLen;
            }
        }
    }
    return unzippedLen;
}

void FastGzFileReader::closeIO() {
    fileReader.closeIO();
    MemoryUtil::safeFree(isalExtraBuffer);
    MemoryUtil::safeFree(isalInputBuffer);
}

bool FastGzFileReader::readMoreDataAndReset(uint8_t* isalIn, uint32_t isalInLen, uint8_t* isalOut, uint32_t isalOutLen, bool isMultipleHeaders) {
    memset(&isalHeader, 0, sizeof(isal_gzip_header));
    isal_gzip_header_init(&isalHeader);
    memset(&isalState, 0, sizeof(inflate_state));
    isal_inflate_init(&isalState);
    isalState.crc_flag = isMultipleHeaders ? ISAL_GZIP_NO_HDR_VER : ISAL_GZIP_NO_HDR;

    int64_t readLen = isalExtraEach;
    size_t in2CacheLen = isalInLen + readLen;
    isalExtraBuffer = MemoryUtil::safeRealloc<uint8_t>((size_t&)isalExtraLen, isalExtraBuffer, in2CacheLen);
    memcpy(isalExtraBuffer, isalIn, isalInLen);
    size_t actualReadLen = fileReader.readIO(isalExtraBuffer + isalInLen, readLen);

    isalInputBuffer = MemoryUtil::safeRealloc<uint8_t>((size_t&)isalInputLength, isalInputBuffer, in2CacheLen);
    memcpy(isalInputBuffer, isalExtraBuffer, in2CacheLen);

    isalState.avail_in = isalInLen + actualReadLen;
    isalState.next_in = isalInputBuffer;
    isalState.avail_out = isalOutLen;
    isalState.next_out = isalOut;

    // 如果读取的数据为0，说明文件已结束
    if (actualReadLen == 0) {
        LOG_ERROR("Invalid gz file. errno = %d.", ISAL_END_INPUT);
        return false;
    }
    return true;
}

int64_t FastGzFileReader::processMultipleGzipHeaders() {
    for (;;) {
        /* 记录当前解压位置信息，包括原始数据的位置的长度，输出数据的位置和长度 */
        uint8_t* isalIn = isalState.next_in;
        int64_t isalInLen = isalState.avail_in;
        uint8_t* isalOut = isalState.next_out;
        int64_t isalOutLen = isalState.avail_out;

        /* 当遇到有些特别特殊的，有多个头的gz文件时，intel gz需要重置header */
        isal_inflate_reset(&isalState);
        int64_t isalRes = isal_read_gzip_header(&isalState, &isalHeader);
        if (ISAL_DECOMP_OK == isalRes) {/* 解压gz ok，退出*/
            break;
        }
        if (isalRes != ISAL_END_INPUT) {
            LOG_ERROR("invalid gz header, errno %d", isalRes);
            return -1;
        }

        if (!readMoreDataAndReset(isalIn, isalInLen, isalOut, isalOutLen, true)) {
            return -1;
        }
    }
    return ISAL_DECOMP_OK;
}

GzPipeReader::GzPipeReader(uint32_t paral) {
    parallel = paral;
}

int GzPipeReader::openIO() {
    fpGz = bgzf_dopen(fileno(stdin), "rb");
    if (fpGz == nullptr) {
        LOG_ERROR("bgzf open failed.");
        return -1;
    }

    if (fpGz->is_gzip) {
        bgzf_mt(fpGz, parallel, 256);
    }

    return 0;
}

size_t GzPipeReader::readIO(void* pBuffer, size_t readSize) {
    ssize_t readLen = bgzf_read(fpGz, pBuffer, readSize);
    if (readLen < 0) {
        LOG_ERROR("Read file from pipe failed");
        return 0;
    }

    return (size_t)readLen;
}

void GzPipeReader::closeIO() {
    if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    }
}

GzPipeWriter::GzPipeWriter(uint32_t paral) {
    parallel = paral;
}

int GzPipeWriter::openIO() {
    fpGz = bgzf_dopen(fileno(stdout), "wb");
    if (fpGz == nullptr) {
        LOG_ERROR("bgzf open failed.");
        return -1;
    }

    fpGz->compress_level = 6;
    bgzf_mt(fpGz, parallel, 256);
    return 0;
}

size_t GzPipeWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (pBuffer == nullptr || writeLen == 0) {
        return 0;
    }

    size_t len = bgzf_write(fpGz, pBuffer, writeLen);
    if (len != writeLen) {
        LOG_ERROR("Write gz file error, write length = %d, expect %d", len, writeLen);
        return -1;
    }
    return len;
}

void GzPipeWriter::closeIO() {
   if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    } 
}

FastGzPipeReader::FastGzPipeReader() : stateInitialized(false), inputBuffer(nullptr), 
    outputBuffer(nullptr), remainingBuffer(nullptr), inputBufferSize(0), 
    outputBufferSize(0), remainingBufferSize(0), remainingSize(0), streamEnded(false) {
}

int FastGzPipeReader::openIO() {
    // 初始化管道读取器
    int ret = pipeReader.openIO();
    if (ret != 0) {
        return ret;
    }

    // 初始化ISAL解压状态
    isal_inflate_init(&state);
    state.crc_flag = ISAL_GZIP_NO_HDR_VER;
    stateInitialized = true;

    // 初始化gzip头
    isal_gzip_header_init(&gzipHeader);

    // 分配缓冲区
    inputBufferSize = 64 * 1024; // 64KB输入缓冲区
    outputBufferSize = 128 * 1024; // 128KB输出缓冲区
    remainingBufferSize = outputBufferSize;
    
    inputBuffer = MemoryUtil::safeAlloc<uint8_t>(inputBufferSize);
    outputBuffer = MemoryUtil::safeAlloc<uint8_t>(outputBufferSize);
    remainingBuffer = MemoryUtil::safeAlloc<uint8_t>(remainingBufferSize);
    
    if (inputBuffer == nullptr || outputBuffer == nullptr || remainingBuffer == nullptr) {
        LOG_ERROR("Failed to allocate buffers for FastPileGzReader");
        closeIO();
        return -1;
    }

    remainingSize = 0;
    streamEnded = false;
    
    return 0;
}

size_t FastGzPipeReader::readIO(void* pBuffer, size_t readSize) {
    if (pBuffer == nullptr || readSize == 0) {
        return 0;
    }

    if (!stateInitialized) {
        LOG_ERROR("FastPileGzReader not initialized");
        return 0;
    }

    // 直接使用ISAL状态，参考FastGzFileReader的模式
    state.next_out = static_cast<uint8_t*>(pBuffer);
    state.avail_out = readSize;
    
    size_t totalCopied = 0;
    
    while (totalCopied < readSize) {
        // 如果没有更多输入数据，需要从管道读取
        if (state.avail_in == 0) {
            size_t bytesRead = pipeReader.readIO(inputBuffer, inputBufferSize);
            if (bytesRead == 0) {
                // 没有更多数据可读，结束
                streamEnded = true;
                break;
            }
            
            state.next_in = inputBuffer;
            state.avail_in = bytesRead;
            
            // 如果是第一次读取或有新的gzip流，需要读取gzip头
            if (state.avail_in > 0) {
                int ret = isal_read_gzip_header(&state, &gzipHeader);
                if (ret != ISAL_DECOMP_OK) {
                    if (ret == ISAL_END_INPUT) {
                        // 需要更多数据才能解析头
                        continue;
                    } else {
                        LOG_ERROR("Invalid gzip header, error: %d", ret);
                        return totalCopied;
                    }
                }
            }
        }
        
        // 执行解压
        uint32_t ret = isal_inflate(&state);
        
        // 计算本次解压的数据量
        size_t decompressedSize = readSize - state.avail_out;
        totalCopied += decompressedSize;
        
        if (ret == ISAL_DECOMP_OK) {
            // 解压完成，继续下一轮
            continue;
        } else if (ret == ISAL_END_INPUT) {
            // 当前gzip流结束，准备下一个流
            isal_inflate_reset(&state);
            state.crc_flag = ISAL_GZIP_NO_HDR_VER;
            isal_gzip_header_init(&gzipHeader);
            
            // 如果还有输入数据，尝试解析下一个gzip头
            if (state.avail_in > 0) {
                ret = isal_read_gzip_header(&state, &gzipHeader);
                if (ret != ISAL_DECOMP_OK) {
                    if (ret == ISAL_END_INPUT) {
                        // 需要更多数据
                        continue;
                    } else {
                        LOG_ERROR("Invalid gzip header in next stream, error: %d", ret);
                        return totalCopied;
                    }
                }
            }
        } else {
            LOG_ERROR("Decompression failed with error: %d", ret);
            return totalCopied;
        }
        
        // 如果已经满足请求的大小，退出
        if (totalCopied >= readSize) {
            break;
        }
    }
    
    return totalCopied;
}

void FastGzPipeReader::closeIO() {
    pipeReader.closeIO();
    
    // 清理缓冲区
    if (inputBuffer) {
        MemoryUtil::safeFree(inputBuffer);
        inputBuffer = nullptr;
    }
    if (outputBuffer) {
        MemoryUtil::safeFree(outputBuffer);
        outputBuffer = nullptr;
    }
    if (remainingBuffer) {
        MemoryUtil::safeFree(remainingBuffer);
        remainingBuffer = nullptr;
    }
    
    // 重置状态
    remainingSize = 0;
    streamEnded = false;
    inputBufferSize = 0;
    outputBufferSize = 0;
    remainingBufferSize = 0;
    stateInitialized = false;
}
