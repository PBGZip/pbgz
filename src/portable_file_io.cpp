/*
 * portable_file_io.cpp - Portable streaming file I/O for pbgz
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

#include "portable_file_io.h"

#include <cerrno>
#include <cstring>

#include "log/logger.h"

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {

/*
 * 64-bit stdio positioning.  fseeko/ftello are POSIX; Windows provides the
 * _fseeki64/_ftelli64 pair.  Both accept offsets beyond 2 GiB.
 */
int seekFile(FILE* fp, int64_t offset, int whence)
{
#if defined(_WIN32)
    return _fseeki64(fp, (__int64)offset, whence);
#else
    return fseeko(fp, (off_t)offset, whence);
#endif
}

int64_t tellFile(FILE* fp)
{
#if defined(_WIN32)
    return (int64_t)_ftelli64(fp);
#else
    return (int64_t)ftello(fp);
#endif
}

/* Flush the C library buffer and ask the OS to persist the file. */
int syncFile(FILE* fp)
{
    if (fflush(fp) != 0) {
        return -1;
    }
#if defined(_WIN32)
    return _commit(_fileno(fp));
#else
    return fsync(fileno(fp));
#endif
}

}  // namespace

BufferedFileReader::BufferedFileReader(const std::string& name)
    : fileName(name) {
}

BufferedFileReader::~BufferedFileReader() {
    closeIO();
}

int32_t BufferedFileReader::openIO() {
    closeIO();
    if (fileName.empty()) {
        LOG_ERROR("IO Open failed, No file name given");
        return -1;
    }
    fp = fopen(fileName.c_str(), "rb");
    if (fp == nullptr) {
        LOG_ERROR("IO Open failed, file name = %s, errno = %d", fileName.c_str(), errno);
        return -1;
    }

    if (0 != seekFile(fp, 0, SEEK_END)) {
        LOG_ERROR("Failed to seek %s to end: %s", fileName.c_str(), strerror(errno));
        fclose(fp);
        fp = nullptr;
        return -1;
    }
    const int64_t size = tellFile(fp);
    if (size < 0 || 0 != seekFile(fp, 0, SEEK_SET)) {
        LOG_ERROR("Failed to get file size for %s: %s", fileName.c_str(), strerror(errno));
        fclose(fp);
        fp = nullptr;
        return -1;
    }

    fileSize = (size_t)size;
    position = 0;
    eofFlag = false;
    return 0;
}

void BufferedFileReader::closeIO() {
    if (fp != nullptr) {
        fclose(fp);
        fp = nullptr;
    }
    eofFlag = false;
}

size_t BufferedFileReader::readIO(void* pBuffer, size_t readSize) {
    if (fp == nullptr || pBuffer == nullptr || readSize == 0) {
        return 0;
    }

    const size_t n = fread(pBuffer, 1, readSize, fp);
    position += n;
    if (n < readSize) {
        if (ferror(fp)) {
            LOG_ERROR("Read file %s failed: %s", fileName.c_str(), strerror(errno));
        }
        eofFlag = true;
    }
    if (position >= fileSize) {
        eofFlag = true;
    }
    return n;
}

size_t BufferedFileReader::readLine(std::string& line) {
    line.clear();
    if (fp == nullptr) {
        eofFlag = true;
        return 0;
    }

    char buffer[8192];
    while (fgets(buffer, (int)sizeof(buffer), fp) != nullptr) {
        const size_t n = strlen(buffer);
        position += n;

        if (n > 0 && buffer[n - 1] == '\n') {
            size_t lineLen = n - 1;
            if (lineLen > 0 && buffer[lineLen - 1] == '\r') {
                --lineLen;
            }
            line.append(buffer, lineLen);
            if (position >= fileSize) {
                eofFlag = true;
            }
            return line.size();
        }

        line.append(buffer, n);
        if (position >= fileSize) {
            break;
        }
    }

    eofFlag = true;
    return line.size();
}

int32_t BufferedFileReader::seekIO(uint64_t seekOffset, IOWhence whence) {
    if (fp == nullptr) {
        return -1;
    }

    size_t newPos = 0;
    switch (whence) {
    case IOWhence::IO_WHENCE_SET:
        if (seekOffset > fileSize) {
            return -1;
        }
        newPos = (size_t)seekOffset;
        break;
    case IOWhence::IO_WHENCE_CUR:
        if (seekOffset > fileSize - position) {
            return -1;
        }
        newPos = position + (size_t)seekOffset;
        break;
    case IOWhence::IO_WHENCE_END:
        if (seekOffset > fileSize) {
            return -1;
        }
        newPos = fileSize - (size_t)seekOffset;
        break;
    default:
        return -1;
    }

    if (0 != seekFile(fp, (int64_t)newPos, SEEK_SET)) {
        LOG_ERROR("Seek file %s failed: %s", fileName.c_str(), strerror(errno));
        return -1;
    }
    position = newPos;
    clearerr(fp);
    eofFlag = false;
    return 0;
}

BufferedFileWriter::BufferedFileWriter(const std::string& name)
    : fileName(name) {
}

BufferedFileWriter::~BufferedFileWriter() {
    closeIO();
}

int32_t BufferedFileWriter::openIO() {
    closeIO();
    if (fileName.empty()) {
        LOG_ERROR("IO Open failed, No file name given");
        latchWriteError(EINVAL);
        return -1;
    }
    fp = fopen(fileName.c_str(), "w+b");
    if (fp == nullptr) {
        LOG_ERROR("IO Open failed, file name = %s, errno = %d", fileName.c_str(), errno);
        latchWriteError(errno);
        return -1;
    }
    fileSize = 0;
    position = 0;
    return 0;
}

void BufferedFileWriter::closeIO() {
    if (fp == nullptr) {
        return;
    }

    if (0 != fflush(fp)) {
        LOG_ERROR("Flush file %s failed: %s", fileName.c_str(), strerror(errno));
        latchWriteError(errno);
    }
#if defined(_WIN32)
    if (0 != _commit(_fileno(fp))) {
#else
    if (0 != fsync(fileno(fp))) {
#endif
        LOG_ERROR("fsync file %s failed: %s", fileName.c_str(), strerror(errno));
        latchWriteError(errno);
    }
    if (0 != fclose(fp)) {
        LOG_ERROR("Close file %s failed: %s", fileName.c_str(), strerror(errno));
        latchWriteError(errno);
    }
    fp = nullptr;
}

size_t BufferedFileWriter::writeIO(const void* pBuffer, size_t writeLen) {
    if (fp == nullptr || pBuffer == nullptr || writeLen == 0) {
        LOG_ERROR("Input invalid.");
        return 0;
    }

    /* Sequential writes always append at the logical file end, matching the
       mmap writer (writeIOAt may have moved the stdio position earlier). */
    if (position != fileSize) {
        if (0 != seekFile(fp, (int64_t)fileSize, SEEK_SET)) {
            LOG_ERROR("Seek file %s to end failed: %s", fileName.c_str(), strerror(errno));
            latchWriteError(errno);
            return 0;
        }
        position = fileSize;
    }

    const uint8_t* src = (const uint8_t*)pBuffer;
    size_t written = 0;
    while (written < writeLen) {
        const size_t n = fwrite(src + written, 1, writeLen - written, fp);
        if (n == 0) {
            LOG_ERROR("Write file %s failed: %s", fileName.c_str(), strerror(errno));
            latchWriteError(errno != 0 ? errno : EIO);
            return 0;
        }
        written += n;
    }

    position += written;
    fileSize = position;
    return written;
}

int32_t BufferedFileWriter::seekIO(uint64_t seekOffset, IOWhence whence) {
    if (fp == nullptr) {
        return -1;
    }

    size_t newPos = 0;
    switch (whence) {
    case IOWhence::IO_WHENCE_SET:
        if (seekOffset > fileSize) {
            return -1;
        }
        newPos = (size_t)seekOffset;
        break;
    case IOWhence::IO_WHENCE_CUR:
        if (seekOffset > fileSize - position) {
            return -1;
        }
        newPos = position + (size_t)seekOffset;
        break;
    case IOWhence::IO_WHENCE_END:
        if (seekOffset > fileSize) {
            return -1;
        }
        newPos = fileSize - (size_t)seekOffset;
        break;
    default:
        return -1;
    }

    if (0 != seekFile(fp, (int64_t)newPos, SEEK_SET)) {
        LOG_ERROR("Seek file %s failed: %s", fileName.c_str(), strerror(errno));
        return -1;
    }
    position = newPos;
    return 0;
}

int32_t BufferedFileWriter::seekToEnd() {
    if (fp == nullptr) {
        return -1;
    }
    if (0 != seekFile(fp, (int64_t)fileSize, SEEK_SET)) {
        LOG_ERROR("Seek file %s to end failed: %s", fileName.c_str(), strerror(errno));
        return -1;
    }
    position = fileSize;
    return 0;
}

int32_t BufferedFileWriter::writeIOAt(size_t seekOffset, const void* pBuffer, size_t writeLen) {
    if (fp == nullptr || pBuffer == nullptr || writeLen == 0) {
        LOG_ERROR("Input invalid.");
        return 0;
    }

    if (0 != seekFile(fp, (int64_t)seekOffset, SEEK_SET)) {
        LOG_ERROR("Seek file %s to %zu failed: %s", fileName.c_str(), seekOffset, strerror(errno));
        latchWriteError(errno);
        return 0;
    }

    const uint8_t* src = (const uint8_t*)pBuffer;
    size_t written = 0;
    while (written < writeLen) {
        const size_t n = fwrite(src + written, 1, writeLen - written, fp);
        if (n == 0) {
            LOG_ERROR("Write file %s failed: %s", fileName.c_str(), strerror(errno));
            latchWriteError(errno != 0 ? errno : EIO);
            return 0;
        }
        written += n;
    }

    position = seekOffset + written;
    if (position > fileSize) {
        fileSize = position;
    }
    /* Match the mmap writer's historical return value (the new logical size). */
    return (int32_t)fileSize;
}

void BufferedFileWriter::flushIO() {
    if (fp == nullptr) {
        return;
    }
    if (0 != syncFile(fp)) {
        LOG_ERROR("Flush file %s failed: %s", fileName.c_str(), strerror(errno));
        latchWriteError(errno);
    }
}
