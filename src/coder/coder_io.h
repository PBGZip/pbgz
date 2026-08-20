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
     * Stream error flag. Both write overflow (buffer too small) and read-side
     * exhaustion set it, so the higher-level coder can inspect it during decode
     * and propagate the error to the actuator through the return value—instead of
     * silently writing past the end and corrupting the heap, or silently returning
     * '\0' and producing a "successful" file with wrong content.
     *
     * Read-side exhaustion used to return '\0', but '\0' can also be real data
     * (quality value '!'-'!'=0), so the return value alone cannot distinguish a
     * real zero from "nothing left to read"; a flag is used instead.
     */
    enum io_err {
        IO_OK = 0,
        IO_BUF_FULL = -1,    // Write side: data_len has reached data_capacity; writing further overflows
        IO_READ_EMPTY = -2,  // Read side: data_len has reached data_capacity; no data left to read
    };

    /*
     * Constructor for standalone use: in contexts such as trial compression and
     * unit tests, overflow only means "this codec cannot fit the data"; it is a
     * selection criterion rather than a failure, so no aggregation point is
     * attached.
     */
    coder_io(const uint8_t *buff, int32_t buff_len)
    {
        init(buff, buff_len, nullptr, "");
    }

    /*
     * Constructor for ownership by a block-processing pass: overflow here means
     * the block pass itself failed and must be reported to the aggregation point.
     * name is the stream name, used to identify which stream overflowed on error.
     */
    coder_io(const uint8_t *buff, int32_t buff_len, coder_err_sink *sink, const char *name)
    {
        init(buff, buff_len, sink, name);
    }

    /*
     * The single error-setting entry point. err itself remains sticky (only the
     * first error is recorded), and the aggregation point is notified at once—
     * reporting happens at the moment the error is set rather than in the
     * destructor, so whether coder_io is a local or a member, and regardless of
     * destruction order, the error is always visible.
     */
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

    /*
     * Write one byte. When buf is full, nothing is written and IO_BUF_FULL is
     * set—this is the root defense against heap overflow: the original
     * implementation *(data+data_len++)=c performed no check and overflowed
     * directly into the heap. The higher-level coder checks err at the end of
     * decode and returns a negative value on error, prompting the actuator to
     * grow the buffer and retry.
     */
    void putc(uint8_t c)
    {
        if (data_len >= data_capacity) {
            set_err(IO_BUF_FULL);
            return;
        }
        *(data + data_len++) = c;
    }

    /*
     * Read one byte. On read-side exhaustion, returns '\0' (preserving the
     * original behavior so the normal path is unaffected), but sets
     * IO_READ_EMPTY so the higher level can tell a real zero from "nothing left
     * to read". The old silent '\0' made the decompressor think data remained,
     * producing wrong content that was still reported as "success".
     */
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
    /* Stream error flag; see io_err. The coder checks it at the end of decode and returns the error code. */
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
 * Aggregation point for overflow errors.
 *
 * coder_io is a bounded view over the block buffer; processing one block opens
 * over a dozen views. "Did this block pass overflow?" is a property of the block
 * pass as a whole, not of any individual view; it used to be split across a
 * dozen unrelated local err values, only reachable if the caller remembered to
 * check—as a result, SAM checked 12 of them, FASTQ checked none, and the index
 * checked none. Patching in the missing 30 checks would just repeat that
 * forgetfulness; the next newly added stream would still leak.
 *
 * So the answer is collected in one place: the moment a view sets its error it
 * is reported here, and the caller asks once at the exit of the block pass. The
 * aggregation point hangs off the Actuator, and since the executor is created
 * per block and destroyed after use, it is naturally "one per block,
 * thread-exclusive": no cleanup needed and never shared across threads.
 */
struct coder_err_sink
{
    /* First error, sticky; later errors are usually its knock-on effects, so keeping the first one is most useful for diagnosis */
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
