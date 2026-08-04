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

struct coder_io
{
    enum mode
    {
        MENC,
        MDEC,
        MUNSET
    };

    /*
     * 流错误标志。写越界（buf 不够）和读端耗尽都置位，让上层 coder 在 decode 过程中
     * 能查到、通过返回值把错误传给 actuator——不再静默越界写坏堆、或静默返回 '\0'
     * 产出错误却"成功"的文件。
     *
     * 读端耗尽原本返回 '\0'，但 '\0' 本身也可能是真实数据（质量值 '!'-'!'=0），
     * 所以不能靠返回值区分"真读到 0"和"没得读"；用标志区分。
     */
    enum io_err {
        IO_OK = 0,
        IO_BUF_FULL = -1,    // 写端：data_len 已到 data_capacity，再写越界
        IO_READ_EMPTY = -2,  // 读端：data_len 已到 data_capacity，无数据可读
    };

    coder_io(const uint8_t *buff, int32_t buff_len)
    {
        data = (uint8_t *)buff;
        data_capacity = buff_len;
        data_len = 0;
        err = IO_OK;
        meta.clear();
        m = MUNSET;
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
     * 写一个字节。buf 满时不再写、置 IO_BUF_FULL——这是堆溢出的根源防线：
     * 原实现 *(data+data_len++)=c 不检查，越界直接写坏堆。上层 coder 在 decode
     * 结束时查 err，有错返回负数，actuator 据此扩容重试。
     */
    void putc(uint8_t c)
    {
        if (data_len >= data_capacity) {
            err = IO_BUF_FULL;
            return;
        }
        *(data + data_len++) = c;
    }

    /*
     * 读一个字节。读端耗尽返回 '\0'（保持原行为，不破坏正常路径），但置
     * IO_READ_EMPTY 让上层可区分"真读到 0"与"没得读"。原来静默返回 '\0' 会让
     * 解压器以为还有数据，产出错误内容却"成功"。
     */
    uint8_t getc()
    {
        if (data_len >= data_capacity) {
            err = IO_READ_EMPTY;
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
    /* 流错误标志，见 io_err；coder 在 decode 结束时检查并返回错误码 */
    int32_t err;
    /* Encoder parameter input and output metadata interaction through meta */
    Json::Value meta;
};

#endif
