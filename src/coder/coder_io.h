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

    /*
     * 独立使用的构造：试压、单元测试这类场合，越界只说明"这个编码器装不下"，
     * 是选择依据而不是故障，所以不接汇聚点。
     */
    coder_io(const uint8_t *buff, int32_t buff_len)
    {
        init(buff, buff_len, nullptr, "");
    }

    /*
     * 归属某次块处理的构造：越界就是这次块处理失败，必须报到汇聚点。
     * name 是流名，出错时用来说清是哪一路流越的界。
     */
    coder_io(const uint8_t *buff, int32_t buff_len, coder_err_sink *sink, const char *name)
    {
        init(buff, buff_len, sink, name);
    }

    /*
     * 唯一置错入口。err 本身仍是粘性的（只记第一个），同时立刻上报汇聚点——
     * 上报放在置错的瞬间而不是析构时，这样 coder_io 是局部量还是成员、
     * 谁先析构，都不影响错误能否被看到。
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
     * 写一个字节。buf 满时不再写、置 IO_BUF_FULL——这是堆溢出的根源防线：
     * 原实现 *(data+data_len++)=c 不检查，越界直接写坏堆。上层 coder 在 decode
     * 结束时查 err，有错返回负数，actuator 据此扩容重试。
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
     * 读一个字节。读端耗尽返回 '\0'（保持原行为，不破坏正常路径），但置
     * IO_READ_EMPTY 让上层可区分"真读到 0"与"没得读"。原来静默返回 '\0' 会让
     * 解压器以为还有数据，产出错误内容却"成功"。
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
    /* 流错误标志，见 io_err；coder 在 decode 结束时检查并返回错误码 */
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

    /* 定义在 coder_err_sink 之后 */
    void report(int32_t e);

    coder_err_sink *err_sink;
    const char *stream_name;
};

/*
 * 越界错误的汇聚点。
 *
 * coder_io 是块缓冲上的一个有界视图，处理一个块要开出十几个视图。"这次块处理
 * 有没有越界"是块处理这一整件事的性质，不是每个视图各自的性质；原来它被拆成
 * 十几个互不相干的局部 err，能不能被问到全凭调用方自觉——结果 SAM 问了 12 个、
 * FASTQ 一个都没问、索引一个都没问。补上缺的那 30 处只是把自觉重复一遍，下一个
 * 新增的流照样会漏。
 *
 * 所以把答案收到一处：视图置错的瞬间就报到这里，调用方在块处理的出口问一次。
 * 汇聚点挂在 Actuator 上，而执行器是每块新建、用完即删的，因此它天然是
 * "每块一个、线程独享"，不需要清理，也不存在跨线程共享。
 */
struct coder_err_sink
{
    /* 首个错误，粘性；后续错误多半是它的连锁反应，留第一个最有定位价值 */
    int32_t err = coder_io::IO_OK;
    /* 出错的流名 */
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
