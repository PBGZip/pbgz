/*
 * coder.h - Header file for pbgz project
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

#ifndef _CODER_H_
#define _CODER_H_


#include <stdint.h>
#include <cstdarg>
#include <functional>

#include "coder_io.h"

using coder_exit_func = std::function<void(int, const char*)>;

using coder_alloc_func = std::function<uint8_t*(size_t)> ;

using coder_realloc_func = std::function<uint8_t*(size_t&, uint8_t*, size_t)>;

using coder_free_func = std::function<void(void*&)>;

using coder_logger_func = std::function<void(int, const char*)>;


namespace coder_ns {
    void register_alloc_proc(coder_alloc_func proc); 

    void register_realloc_proc(coder_realloc_func proc);

    void register_exit_proc(coder_exit_func proc);

    void resister_logger_proc(coder_logger_func proc);

    void register_free_func(coder_free_func proc);

    void initFcCoder();

    enum coder_errcode {
        CODER_ERR_BAD_ARGS = (-1),
        CODER_ERR_MEM_ALLOC_FAIL = (CODER_ERR_BAD_ARGS) - 1,
        CODER_ERR_INNER = (CODER_ERR_MEM_ALLOC_FAIL) -1,
        CODER_ERR_INVALID_STATE = CODER_ERR_INNER - 1,

    };

    enum coder_log_level {
        DEBUGGING,
        INFO,
        WARNING,
        ERROR,
        FATAL,
    };
}

void* safe_realloc_init(uint32_t& size, uint8_t* ptr, size_t new_size, char ch); 

void* safe_realloc(uint32_t& size, uint8_t* ptr, size_t new_size); 

void* safe_alloc(size_t size);

void* safe_alloc_init(size_t size, char ch);

void* safe_calloc(size_t num, size_t size) ;

void safe_free(void** ptr) ;

void coder_logger(coder_ns::coder_log_level level, const char* log_format, ...);

// The following functions will cause process exit, not recommended for use within coder library
void coder_exit(int16_t exit_code, const char* exit_msg_format, ...);

void check_exit(bool condition, int16_t exit_code, const char* exit_msg_format, ...);


class coder
{
public:
    virtual ~coder() {}

    /*
     * 编码一段数据。
     *
     * 统一取三参数形式：coder_bwt_cm 原本只有两个参数，为了能通过基类指针调用而补齐，
     * 它内部并不使用 need2hold。第三个参数有默认值，两参数的调用方式不受影响。
     *
     * 默认空实现：只做解码、不做编码的子类不必覆盖。
     */
    virtual void encode_line([[maybe_unused]] const uint8_t *in,
                             [[maybe_unused]] const uint32_t in_len,
                             [[maybe_unused]] bool need2hold = false) {}

    /* 编码收尾，把残留在内部缓冲里的数据全部吐出去。默认空实现。 */
    virtual void encode_flush() {}

    /*
     * 本编码器能否以"逐行累积"的方式使用，即对同一个 coder_io 反复调用 encode_line，
     * 最后统一 encode_flush。
     *
     * 这不是可有可无的偏好，而是硬约束。有的编码器（如 coder_fc）在内部一次性完成
     * 预处理、变换和熵编码，第二次调用 encode_line 会直接报错退出，输入太短也会报错。
     * 执行器里各字段的用法并不统一：QNAME、RNAME、CIGAR 这些是按行循环喂进去的，
     * 而 SEQ 是先拼进一块临时缓冲再整块喂一次。
     *
     * 预处理阶段试压时统一按整块方式测量，所以试压结果不能无条件套用到逐行字段上，
     * 必须先用本函数过滤掉不支持逐行的编码器。
     */
    virtual bool supportsLineMode() const { return true; }

    /*
     * 声明本编码器适用于哪种文件类型的哪个字段。
     *
     * 绝大多数编码器是通用的字节流压缩器，放到哪个字段上都能跑，所以默认返回 true。
     * 少数编码器有额外的前置依赖——比如质量值上下文混合编码器需要每条记录的读长和
     * 链方向，只有 SAM 的 QUAL 列能提供——这类编码器覆盖本函数，把自己的适用范围
     * 声明清楚，预处理阶段试压时会先问一句，不适用就直接跳过，不会被误选。
     *
     * fileType 取 BlockType 的值，这里用 uint32_t 是为了避免 coder 层反向依赖
     * io_block.h。fieldIdx 的含义随 fileType 而定：SAM 用 SamField，FASTQ 用
     * FastqField，两套编号各自独立。
     */
    virtual bool supports([[maybe_unused]] uint32_t fileType,
                          [[maybe_unused]] uint32_t fieldIdx) const { return true; }

    virtual int32_t decode_line([[maybe_unused]] uint8_t *dst, 
                                [[maybe_unused]] uint32_t len, 
                                [[maybe_unused]] uint8_t split_ch = '\n', 
                                [[maybe_unused]] bool need2hold = false) { return 0; }
    virtual int32_t decode_line([[maybe_unused]] uint8_t *dst, 
                                [[maybe_unused]] uint32_t len, 
                                [[maybe_unused]] uint8_t *rely = nullptr, 
                                [[maybe_unused]] uint8_t split_ch = '\n', 
                                [[maybe_unused]] bool need2hold = false)
    {
        return 0;
    }

    virtual void set_level(int32_t level)
    {
        io->set_level(level);
    }

protected:
    coder_io *io;
};

#endif
