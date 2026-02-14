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
