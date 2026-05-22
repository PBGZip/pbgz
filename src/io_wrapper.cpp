/*
 * io_wrapper.cpp - Source file for pbgz project
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
#include <sys/stat.h>
#include <string.h>

#include "io_wrapper.h"
#include "log/logger.h"
#include "utils/memory_util.h"


int FileOperator::openIO() {
    if (fileName.empty()) {
        LOG_ERROR("IO Open failed, No file name given");
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
        LOG_ERROR("Failed to get file stat");
        close(fd);
    }

    fileSize = fileStat.st_size;

    return 0;
}

int32_t FileOperator::mmapFile(size_t mapFileSize) {
    if (mapFileSize == 0) {
        return 0;
    }
    
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
        if (seekOffset > fileSize) {
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
    
    const uint8_t* src = fo.mappedAddress + fo.position;
    uint8_t* dst = (uint8_t*)pBuffer;
    
    // Conservative SIMD optimization: only for large data (>= 1MB)
    // Use aligned SIMD copy for better performance on large reads
    if (realRead >= 1024 * 1024) {
        // Enhanced prefetch for large data
        const size_t prefetchDistance = 2048;
        
        // Prefetch source regions
        for (size_t i = 0; i < realRead; i += 32) {
            if (i + prefetchDistance < realRead) {
#if defined(__x86_64__)
                _mm_prefetch((char*)(src + i + prefetchDistance), _MM_HINT_NTA);
                _mm_prefetch((char*)(src + i + prefetchDistance + 64), _MM_HINT_NTA);
#elif defined(__aarch64__)
                __builtin_prefetch(src + i + prefetchDistance, 0, 0);
                __builtin_prefetch(src + i + prefetchDistance + 64, 0, 0);
#endif
            }
        }
        
        // SIMD copy with boundary handling
        size_t i = 0;
        
        // Handle leading unaligned bytes (<16 bytes)
        size_t misalignment = ((uintptr_t)(src + i)) & 15;
        if (misalignment > 0 && i < realRead) {
            size_t prefix = (misalignment > (realRead - i)) ? (realRead - i) : misalignment;
            memcpy(dst + i, src + i, prefix);
            i += prefix;
        }
        
        // Main SIMD copy loop (16-byte aligned)
        if (i + 16 <= realRead) {
#if defined(__x86_64__) && defined(__AVX2__)
            while (i + 16 <= realRead) {
                __m128i data = _mm_loadu_si128((__m128i*)(src + i));
                _mm_storeu_si128((__m128i*)(dst + i), data);
                i += 16;
            }
#elif defined(__aarch64__)
            while (i + 16 <= realRead) {
                uint8x16_t data = vld1q_u8(src + i);
                vst1q_u8(dst + i, data);
                i += 16;
            }
#endif
        }
        
        // Handle trailing bytes
        if (i < realRead) {
            memcpy(dst + i, src + i, realRead - i);
        }
    } else {
        // Small data: use standard memcpy and light prefetch
        if (realRead >= 64 * 1024) {
            const size_t prefetchDistance = 2048;
            for (size_t i = 0; i < realRead; i += 32) {
                if (i + prefetchDistance < realRead) {
#if defined(__x86_64__)
                    _mm_prefetch((char*)(src + i + prefetchDistance), _MM_HINT_NTA);
                    _mm_prefetch((char*)(src + i + prefetchDistance + 64), _MM_HINT_NTA);
#elif defined(__aarch64__)
                    __builtin_prefetch(src + i + prefetchDistance, 0, 0);
                    __builtin_prefetch(src + i + prefetchDistance + 64, 0, 0);
#endif
                }
            }
        }
        memcpy(dst, src, realRead);
    }
    
    fo.position += realRead;
    if (fo.position >= fo.fileSize) {
        eofFlag = true;
    }
    return realRead;
}

void* FileReader::getAt(size_t pos) {
    return fo.mappedAddress + pos;
}

size_t FileReader::readLine(std::string& line) {
    line.clear();
    
    if (fo.mappedAddress == nullptr || fo.position >= fo.fileSize) {
        eofFlag = true;
        return 0;
    }
    
    size_t startPos = fo.position;
    const size_t remaining = fo.fileSize - fo.position;
    const char* current = (const char*)(fo.mappedAddress + fo.position);
    
    // Conservative SIMD optimization: only for long lines (>500 bytes)
    // Avoid overhead for common short lines, use standard method for consistency
    if (remaining > 500) {
#if defined(__x86_64__) && defined(__AVX2__)
        const __m256i newline_vec = _mm256_set1_epi8('\n');
        size_t processed = 0;
        
        while (processed + 32 <= remaining) {
            __m256i data = _mm256_loadu_si256((__m256i*)(current + processed));
            __m256i cmp = _mm256_cmpeq_epi8(data, newline_vec);
            unsigned mask = _mm256_movemask_epi8(cmp);
            
            if (mask != 0) {
                // Found newline: use standard method for consistent behavior
                const void* newlinePos = memchr(current + processed, '\n', remaining - processed);
                if (newlinePos != nullptr) {
                    const char* foundPos = static_cast<const char*>(newlinePos);
                    size_t lineLength = foundPos - current;
                    
                    // Handle carriage return
                    if (lineLength > 0 && current[lineLength - 1] == '\r') {
                        lineLength--;
                    }
                    
                    if (lineLength > 0) {
                        line.append(current, lineLength);
                    }
                    
                    fo.position = foundPos - (const char*)fo.mappedAddress + 1; // Skip newline
                    
                    if (fo.position >= fo.fileSize) {
                        eofFlag = true;
                    }
                    
                    return line.length();
                }
                break;
            }
            processed += 32;
        }
#endif
    }
    
    // Standard fallback method for all cases
    const void* newlinePos = memchr(current, '\n', remaining);
    
    if (newlinePos != nullptr) {
        const char* foundPos = static_cast<const char*>(newlinePos);
        fo.position = foundPos - (const char*)fo.mappedAddress;
    } else {
        // No newline found, go to end of file
        fo.position = fo.fileSize;
    }
    
    size_t lineLength = fo.position - startPos;
    
    // If newline found, skip it
    if (fo.position < fo.fileSize && fo.mappedAddress[fo.position] == '\n') {
        fo.position++;
    }
    
    // Handle Windows-style line breaks \r\n
    if (lineLength > 0 && fo.mappedAddress[startPos + lineLength - 1] == '\r') {
        lineLength--;
    }
    
    // Build string
    if (lineLength > 0) {
        line.append(reinterpret_cast<char*>(fo.mappedAddress + startPos), lineLength);
    }
    
    if (fo.position >= fo.fileSize) {
        eofFlag = true;
    }
    
    return lineLength;
}

int32_t FileWriter::openIO() {
    if (0 != fo.openIO()) {
        LOG_ERROR("File writer open failed.");
        return -1;
    }
    mapSize = calculateNextMapSize();
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
    return ((fileSize / initialMapSize) + 1) * initialMapSize;
}

size_t FileWriter::calculateNextMapSize() {
    const size_t INITIAL_SIZE = 64 * 1024 * 1024; // 64MB initial size
    
    if (initialMapSize == 0) {
        initialMapSize = INITIAL_SIZE;
        return INITIAL_SIZE;
    }
    
    return mapSize * 2;
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
        // Flush content to disk first - use ASYNC for better performance
        msync(fo.mappedAddress, fo.fileSize, MS_ASYNC);
        // Unmap
        munmap(fo.mappedAddress, mapSize);
        // Calculate new mapping size using exponential growth strategy
        while (mapSize < newSize) {
            mapSize = calculateNextMapSize();
        }
        // Extend file
        ftruncate(fo.fd, mapSize);
        // Remap file
        if (0 != fo.mmapFile(mapSize)) {
            LOG_ERROR("File not mapped");
            return 0;
        }
    }
    
    const uint8_t* src = (const uint8_t*)pBuffer;
    uint8_t* dst = fo.mappedAddress + fo.fileSize;
    
    // Conservative SIMD optimization: only for large data (>= 1MB)
    // Use aligned SIMD copy for better performance on large writes
    if (writeLen >= 1024 * 1024) {
        // Prefetch destination regions for write-combining optimization
        const size_t prefetchDistance = 2048;
        for (size_t i = 0; i < writeLen; i += 32) {
            if (i + prefetchDistance < writeLen) {
#if defined(__x86_64__)
                _mm_prefetch((char*)(dst + i + prefetchDistance), _MM_HINT_NTA);
                _mm_prefetch((char*)(dst + i + prefetchDistance + 64), _MM_HINT_NTA);
#elif defined(__aarch64__)
                __builtin_prefetch(dst + i + prefetchDistance, 0, 0);
                __builtin_prefetch(dst + i + prefetchDistance + 64, 0, 0);
#endif
            }
        }
        
        // SIMD copy with boundary handling
        size_t i = 0;
        
        // Handle leading unaligned bytes (<16 bytes)
        size_t misalignment = (16 - (uintptr_t)(dst + i)) & 15;
        if (misalignment > 0 && i < writeLen) {
            size_t prefix = (misalignment > (writeLen - i)) ? (writeLen - i) : misalignment;
            memcpy(dst + i, src + i, prefix);
            i += prefix;
        }
        
        // Main SIMD copy loop (16-byte aligned)
        if (i + 16 <= writeLen) {
#if defined(__x86_64__) && defined(__AVX2__)
            while (i + 16 <= writeLen) {
                __m128i data = _mm_loadu_si128((__m128i*)(src + i));
                _mm_storeu_si128((__m128i*)(dst + i), data);
                i += 16;
            }
#elif defined(__aarch64__)
            while (i + 16 <= writeLen) {
                uint8x16_t data = vld1q_u8(src + i);
                vst1q_u8(dst + i, data);
                i += 16;
            }
#endif
        }
        
        // Handle trailing bytes
        if (i < writeLen) {
            memcpy(dst + i, src + i, writeLen - i);
        }
    } else {
        // Small data: use standard memcpy and light prefetch
        if (writeLen >= 64 * 1024) {
            const size_t prefetchDistance = 2048;
            for (size_t i = 0; i < writeLen; i += 32) {
                if (i + prefetchDistance < writeLen) {
#if defined(__x86_64__)
                    _mm_prefetch((char*)(dst + i + prefetchDistance), _MM_HINT_NTA);
                    _mm_prefetch((char*)(dst + i + prefetchDistance + 64), _MM_HINT_NTA);
#elif defined(__aarch64__)
                    __builtin_prefetch(dst + i + prefetchDistance, 0, 0);
                    __builtin_prefetch(dst + i + prefetchDistance + 64, 0, 0);
#endif
                }
            }
        }
        memcpy(dst, src, writeLen);
    }
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
    size_t requiredSize = seekOffset + writeLen;
    
    if (requiredSize > mapSize) {
         // If file size exceeds mapping size, need to remap
        // Flush content to disk first - use ASYNC for better performance
        msync(fo.mappedAddress, fo.fileSize, MS_ASYNC);
        // Unmap
        munmap(fo.mappedAddress, mapSize);
        // New file size
        newFileSize = requiredSize;
        // Calculate new mapping size using exponential growth strategy
        while (mapSize < requiredSize) {
            mapSize = calculateNextMapSize();
        }
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
    
    // Loop read until reaching requested size or no more data
    while (totalRead < readSize) {
        ssize_t bytesRead = read(STDIN_FILENO, buffer + totalRead, readSize - totalRead);
        if (bytesRead < 0) {
            if (errno == EINTR) {
                continue; // Interrupted by signal, retry
            } else if (errno == EAGAIN || errno == EWOULDBLOCK) {
                break; // No data available in non-blocking mode
            } else {
                return -1; // Other errors
            }
        } else if (bytesRead == 0) {
            eofFlag = true;
            break; // EOF, write end closed
        }
        
        totalRead += bytesRead;
    }
    
    return totalRead;
}

size_t PipeReader::readLine(std::string& line) {
    line.clear();
    
    // Initialize buffer (if not already allocated)
    if (lineBuffer == nullptr) {
        bufferSize = 8192;  // 8KB buffer
        lineBuffer = new char[bufferSize];
        bufferPos = 0;
        bytesRead = 0;
    }
    
    size_t totalRead = 0;
    
    while (true) {
        // If buffer data has been processed, read new data
        if (bufferPos >= static_cast<size_t>(bytesRead)) {
            bytesRead = read(STDIN_FILENO, lineBuffer, bufferSize);
            if (bytesRead < 0) {
                if (errno == EINTR) {
                continue; // Interrupted by signal, retry
                } else if (errno == EAGAIN || errno == EWOULDBLOCK) {
                    break; // No data available in non-blocking mode
                } else {
                    return -1; // Other errors
                }
            } else if (bytesRead == 0) {
                // EOF or error, return number of characters read
                eofFlag = true;
                return totalRead;
            }
            bufferPos = 0;
        }
        
        // Process characters in buffer
        while (bufferPos < static_cast<size_t>(bytesRead)) {
            char c = lineBuffer[bufferPos++];
            totalRead++;
            
            // If newline encountered, stop reading
            if (c == '\n') {
                return totalRead;
            }
            
            // If carriage return encountered, skip it (Windows-style line breaks)
            if (c == '\r') {
                continue;
            }
            
            // Add character to line
            line += c;
        }
    }
    
    return totalRead;
}

PipeReader::~PipeReader() {
    // Clean up buffer
    if (lineBuffer != nullptr) {
        delete[] lineBuffer;
        lineBuffer = nullptr;
    }
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
        return -1;
    } else if (readLen == 0) {
        eofFlag = true;
        return 0;
    }

    return (size_t)readLen;
}

void GzFileReader::closeIO() {
    if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    }
    
    // Clean up buffer
    if (lineBuffer != nullptr) {
        delete[] lineBuffer;
        lineBuffer = nullptr;
    }
}

size_t GzFileReader::readLine(std::string& line) {
    line.clear();
    
    if (fpGz == nullptr) {
        return 0;
    }
    
    // Initialize buffer (if not already allocated)
    if (lineBuffer == nullptr) {
        bufferSize = 8192;  // 8KB buffer
        lineBuffer = new char[bufferSize];
        bufferPos = 0;
        bytesRead = 0;
    }
    
    size_t totalRead = 0;
    
    while (true) {
        // If buffer data has been processed, read new data
        if (bufferPos >= static_cast<size_t>(bytesRead)) {
            bytesRead = bgzf_read(fpGz, lineBuffer, bufferSize);
            if (bytesRead <= 0) {
                // EOF or error, return number of characters read
                if (bytesRead == 0) {
                    eofFlag = true;
                }
                return totalRead;
            }
            bufferPos = 0;
        }
        
        // Process characters in buffer
        while (bufferPos < static_cast<size_t>(bytesRead)) {
            char c = lineBuffer[bufferPos++];
            totalRead++;
            
            // If newline encountered, stop reading
            if (c == '\n') {
                return totalRead;
            }
            
            // If carriage return encountered, skip it (Windows-style line breaks)
            if (c == '\r') {
                continue;
            }
            
            // Add character to line
            line += c;
        }
    }
    
    return totalRead;
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

            // Read more data and reset state
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
                // Handle multiple gzip headers case
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
    
    // Clean up buffer
    if (lineBuffer != nullptr) {
        delete[] lineBuffer;
        lineBuffer = nullptr;
    }
}

size_t FastGzFileReader::readLine(std::string& line) {
    line.clear();
    
    // Initialize buffer (if not already allocated)
    if (lineBuffer == nullptr) {
        bufferSize = 8192;  // 8KB buffer
        lineBuffer = new char[bufferSize];
        bufferPos = 0;
        bytesRead = 0;
    }
    
    size_t totalRead = 0;
    
    while (true) {
        // If buffer data has been processed, read new data
        if (bufferPos >= bytesRead) {
            bytesRead = readIO(lineBuffer, bufferSize);
            if (bytesRead == 0) {
                // EOF, return number of characters read
                return totalRead;
            }
            bufferPos = 0;
        }
        
        // Process characters in buffer
        while (bufferPos < bytesRead) {
            char c = lineBuffer[bufferPos++];
            totalRead++;
            
            // If newline encountered, stop reading
            if (c == '\n') {
                return totalRead;
            }
            
            // If carriage return encountered, skip it (Windows-style line breaks)
            if (c == '\r') {
                continue;
            }
            
            // Add character to line
            line += c;
        }
    }
    
    return totalRead;
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

    // If read data is 0, it means file has ended
    if (actualReadLen == 0) {
        LOG_ERROR("Invalid gz file. errno = %d.", ISAL_END_INPUT);
        return false;
    }
    return true;
}

int64_t FastGzFileReader::processMultipleGzipHeaders() {
    for (;;) {
        /* Record current decompression position info, including original data position and length, output data position and length */
        uint8_t* isalIn = isalState.next_in;
        int64_t isalInLen = isalState.avail_in;
        uint8_t* isalOut = isalState.next_out;
        int64_t isalOutLen = isalState.avail_out;

        /* When encountering some special gz files with multiple headers, intel gz needs to reset header */
        isal_inflate_reset(&isalState);
        int64_t isalRes = isal_read_gzip_header(&isalState, &isalHeader);
        if (ISAL_DECOMP_OK == isalRes) {/* Gz decompression ok, exit*/
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
    lineBuffer = nullptr;
    bufferSize = 0;
    bufferPos = 0;
    bytesRead = 0;
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
    if (readSize == 0 || pBuffer == nullptr) {
        return 0;
    }
    ssize_t readLen = bgzf_read(fpGz, pBuffer, readSize);
    if (readLen < 0) {
        LOG_ERROR("Read file from pipe failed");
        return -1;
    } else if (readLen == 0) {
        eofFlag = true;
        return 0;
    }

    return (size_t)readLen;
}

void GzPipeReader::closeIO() {
    if (fpGz != nullptr) {
        bgzf_close(fpGz);
        fpGz = nullptr;
    }
    
    // Clean up buffer
    if (lineBuffer != nullptr) {
        delete[] lineBuffer;
        lineBuffer = nullptr;
    }
}

size_t GzPipeReader::readLine(std::string& line) {
    line.clear();
    
    if (fpGz == nullptr) {
        return 0;
    }
    
    // Initialize buffer (if not already allocated)
    if (lineBuffer == nullptr) {
        bufferSize = 8192;  // 8KB buffer
        lineBuffer = new char[bufferSize];
        bufferPos = 0;
        bytesRead = 0;
    }
    
    size_t totalRead = 0;
    
    while (true) {
        // If buffer data has been processed, read new data
        if (bufferPos >= static_cast<size_t>(bytesRead)) {
            bytesRead = bgzf_read(fpGz, lineBuffer, bufferSize);
            if (bytesRead <= 0) {
                // EOF or error, return number of characters read
                if (bytesRead == 0) {
                    eofFlag = true;
                }
                return totalRead;
            }
            bufferPos = 0;
        }
        
        // Process characters in buffer
        while (bufferPos < static_cast<size_t>(bytesRead)) {
            char c = lineBuffer[bufferPos++];
            totalRead++;
            
            // If newline encountered, stop reading
            if (c == '\n') {
                return totalRead;
            }
            
            // If carriage return encountered, skip it (Windows-style line breaks)
            if (c == '\r') {
                continue;
            }
            
            // Add character to line
            line += c;
        }
    }
    
    return totalRead;
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
    // Initialize pipe reader
    int ret = pipeReader.openIO();
    if (ret != 0) {
        return ret;
    }

    // Initialize ISAL decompression state
    isal_inflate_init(&state);
    state.crc_flag = ISAL_GZIP_NO_HDR_VER;
    stateInitialized = true;

    // Initialize gzip header
    isal_gzip_header_init(&gzipHeader);

    // Allocate buffers
    inputBufferSize = 64 * 1024; // 64KB input buffer
    outputBufferSize = 128 * 1024; // 128KB output buffer
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

    // Directly use ISAL state, referencing FastGzFileReader pattern
    state.next_out = static_cast<uint8_t*>(pBuffer);
    state.avail_out = readSize;
    
    size_t totalCopied = 0;
    
    while (totalCopied < readSize) {
        // If no more input data available, read from pipe
        if (state.avail_in == 0) {
            size_t bytesRead = pipeReader.readIO(inputBuffer, inputBufferSize);
            if (bytesRead == 0) {
                // No more data to read, end
                streamEnded = true;
                break;
            }
            
            state.next_in = inputBuffer;
            state.avail_in = bytesRead;
            
            // If first read or new gzip stream, need to read gzip header
            if (state.avail_in > 0) {
                int ret = isal_read_gzip_header(&state, &gzipHeader);
                if (ret != ISAL_DECOMP_OK) {
                    if (ret == ISAL_END_INPUT) {
                        // Need more data to parse header
                        continue;
                    } else {
                        LOG_ERROR("Invalid gzip header, error: %d", ret);
                        return totalCopied;
                    }
                }
            }
        }
        
        // Perform decompression
        uint32_t ret = isal_inflate(&state);
        
        // Calculate amount of data decompressed this time
        size_t decompressedSize = readSize - state.avail_out;
        totalCopied += decompressedSize;
        
        if (ret == ISAL_DECOMP_OK) {
            // Decompression complete, continue to next round
            continue;
        } else if (ret == ISAL_END_INPUT) {
            // Current gzip stream ended, prepare for next stream
            isal_inflate_reset(&state);
            state.crc_flag = ISAL_GZIP_NO_HDR_VER;
            isal_gzip_header_init(&gzipHeader);
            
            // If still have input data, try to parse next gzip header
            if (state.avail_in > 0) {
                ret = isal_read_gzip_header(&state, &gzipHeader);
                if (ret != ISAL_DECOMP_OK) {
                    if (ret == ISAL_END_INPUT) {
                        // Need more data
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
        
        // If requested size is satisfied, exit
        if (totalCopied >= readSize) {
            break;
        }
    }
    
    return totalCopied;
}

void FastGzPipeReader::closeIO() {
    pipeReader.closeIO();
    
    // Clean up buffers
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
    
    // Reset state
    remainingSize = 0;
    streamEnded = false;
    inputBufferSize = 0;
    outputBufferSize = 0;
    remainingBufferSize = 0;
    stateInitialized = false;
}

size_t FastGzPipeReader::readLine(std::string& line) {
    line.clear();
    
    // Use larger buffer for batch reading to improve efficiency
    const size_t BUFFER_SIZE = 8192;  // 8KB buffer
    char buffer[BUFFER_SIZE];
    size_t totalRead = 0;
    size_t bufferPos = 0;
    size_t bytesRead = 0;
    
    while (true) {
        // If buffer data has been processed, read new data
        if (bufferPos >= bytesRead) {
            bytesRead = readIO(buffer, BUFFER_SIZE);
            if (bytesRead == 0) {
                // EOF, return number of characters read
                eofFlag = true;
                return totalRead;
            }
            bufferPos = 0;
        }
        
        // Process characters in buffer
        while (bufferPos < bytesRead) {
            char c = buffer[bufferPos++];
            totalRead++;
            
            // If newline encountered, stop reading
            if (c == '\n') {
                return totalRead;
            }
            
            // If carriage return encountered, skip it (Windows-style line breaks)
            if (c == '\r') {
                continue;
            }
            
            // Add character to line
            line += c;
        }
    }
    
    return totalRead;
}
