/*
 * io_interface.h - Platform-independent I/O interfaces for pbgz
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
#include <stddef.h>
#include <string>

class ChecksumStream;

enum class IOWhence {
    IO_WHENCE_SET = 0,
    IO_WHENCE_CUR = 1,
    IO_WHENCE_END = 2,
};

class IOReader {
public:
    IOReader() {
        eofFlag = false;
    }
    virtual int32_t openIO() = 0;

    virtual void closeIO() = 0;

    virtual size_t readIO(void* pBuffer, size_t readSize) = 0;

    virtual ~IOReader() {};

    virtual size_t readLine(std::string& line) = 0;

    /* Line readers that can expose the raw bytes (including delimiters) feed
       them to the sink while parsing; used to hash a reference FASTA in the
       same pass that builds its squash. */
    virtual size_t readLine(std::string& line, ChecksumStream* /*sink*/) {
        return readLine(line);
    }

    virtual bool isEOF() { return eofFlag; }

protected:
    bool eofFlag;
};

class IOWriter {
public:
    virtual int32_t openIO() = 0;

    virtual void closeIO() = 0;

    virtual size_t writeIO(const void* pBuffer, size_t writeLen) = 0;

    virtual ~IOWriter() {}

    /*
     * Sticky write error: writes run on a separate thread and closeIO() returns
     * void, so the first error is latched here for the writer thread to pick up.
     */
    void latchWriteError(int32_t e) { if (writeErr == 0) writeErr = e; }

    int32_t getWriteError() const { return writeErr; }

protected:
    int32_t writeErr = 0;
};

/*
 * Capabilities shared by the mmap file reader and the portable streaming file
 * reader.  Callers use these interfaces instead of the concrete mmap classes so
 * either implementation can be selected without changing the data path.
 */
class SeekableReader {
public:
    virtual ~SeekableReader() {}
    virtual int32_t seekIO(uint64_t seekOffset,
                           IOWhence whence = IOWhence::IO_WHENCE_SET) = 0;
    virtual size_t getFileSize() = 0;
    virtual size_t getCurrentPos() = 0;
};

class SeekableWriter {
public:
    virtual ~SeekableWriter() {}
    virtual int32_t seekIO(uint64_t seekOffset,
                           IOWhence whence = IOWhence::IO_WHENCE_SET) = 0;
    virtual int32_t seekToEnd() = 0;
    virtual int32_t writeIOAt(size_t seekOffset, const void* pBuffer, size_t writeLen) = 0;
    virtual void flushIO() = 0;
    virtual size_t getCurrentPos() = 0;
    virtual size_t getFileSize() = 0;
};
