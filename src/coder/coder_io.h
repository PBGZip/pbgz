/*
 * coder_io.h - Header file for pbgz project
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

#ifndef _CODER_IO_H_
#define _CODER_IO_H_

#include <stdint.h>
#include <json/json.h>

struct coder_err_sink;

struct coder_io
{
    enum mode
    {
        MENC,
        MDEC,
        MUNSET
    };

    /*
     * Stream error flag, set by write overflow and read-side exhaustion so the
     * higher-level coder can propagate the error through its return value.
     */
    enum io_err {
        IO_OK = 0,
        IO_BUF_FULL = -1,    // Write side: data_len has reached data_capacity; writing further overflows
        IO_READ_EMPTY = -2,  // Read side: data_len has reached data_capacity; no data left to read
        /* Write side: a symbol outside the coder's alphabet; the stream is failed and the caller chooses another codec. */
        IO_UNCODABLE = -3,
    };

    /* Constructor for standalone use; overflow is a selection criterion, not a
       failure, so no aggregation point is attached. */
    coder_io(const uint8_t *buff, int32_t buff_len)
    {
        init(buff, buff_len, nullptr, "");
    }

    /* Constructor for a block-processing pass; overflow is reported to the
       aggregation point, with name identifying the stream. */
    coder_io(const uint8_t *buff, int32_t buff_len, coder_err_sink *sink, const char *name)
    {
        init(buff, buff_len, sink, name);
    }

    /* Single error-setting entry point; err is sticky and the aggregation point
       is notified at once. */
    void set_err(int32_t e)
    {
        if (err != IO_OK) {
            return;
        }
        err = e;
        report(e);
    }

    /* Append coder identifier */
    void appen_magic(const std::string magic)
    {
        meta["magic"] = magic;
    }

    std::string get_magic() const
    {
        return meta["magic"].asString();
    }

    /* Set level */
    void set_level(int32_t level)
    {
        meta["level"] = level;
    }

    int get_level() const
    {
        return meta["level"].asInt();
    }

    /* Write one byte; on a full buffer nothing is written and IO_BUF_FULL is set. */
    void putc(uint8_t c)
    {
        if (data_len >= data_capacity) {
            set_err(IO_BUF_FULL);
            return;
        }
        *(data + data_len++) = c;
    }

    /* Read one byte; on exhaustion returns '\0' but sets IO_READ_EMPTY so a real
       zero can be told apart from "nothing left to read". */
    uint8_t getc()
    {
        if (data_len >= data_capacity) {
            set_err(IO_READ_EMPTY);
            return '\0';
        }
        return *(data + data_len++);
    }

    /* IO mode */
    mode m;

    uint8_t *data;
    /* Total length of data */
    int32_t data_capacity;
    /* Currently processed length */
    int32_t data_len;
    /* Stream error flag; see io_err. */
    int32_t err;
    /* Encoder parameter input and output metadata interaction through meta */
    Json::Value meta;

private:
    void init(const uint8_t *buff, int32_t buff_len, coder_err_sink *sink, const char *name)
    {
        data = (uint8_t *)buff;
        data_capacity = buff_len;
        data_len = 0;
        err = IO_OK;
        err_sink = sink;
        stream_name = (name != nullptr) ? name : "";
        meta.clear();
        m = MUNSET;
    }

    /* Defined after coder_err_sink */
    void report(int32_t e);

    coder_err_sink *err_sink;
    const char *stream_name;
};

/*
 * Aggregation point for overflow errors: the moment a coder_io view sets its
 * error it is reported here, and the caller asks once at the exit of the block
 * pass. One per block, thread-exclusive.
 */
struct coder_err_sink
{
    /* First error, sticky. */
    int32_t err = coder_io::IO_OK;
    /* Name of the stream that errored */
    const char *what = "";

    bool ok() const { return err == coder_io::IO_OK; }

    void latch(int32_t e, const char *name)
    {
        if (err != coder_io::IO_OK) {
            return;
        }
        err = e;
        what = (name != nullptr) ? name : "";
    }
};

inline void coder_io::report(int32_t e)
{
    if (err_sink != nullptr) {
        err_sink->latch(e, stream_name);
    }
}

#endif
