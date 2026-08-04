/*
 * coder.cpp - Source file for pbgz project
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

#include "coder.h"
#include "coder_fc.h"

coder_alloc_func alloc_proc = nullptr;
coder_realloc_func realloc_proc = nullptr ;
coder_logger_func logger_proc = nullptr;
coder_free_func free_proc = nullptr;

namespace coder_ns {
    void register_alloc_proc(coder_alloc_func proc) {
        alloc_proc = proc;
    }

    void register_realloc_proc(coder_realloc_func proc) {
        realloc_proc = proc;
    }

    void resister_logger_proc(coder_logger_func proc) {
        logger_proc = proc;
    }

    void register_free_func(coder_free_func proc) {
        free_proc = proc;
    }

    void initFcCoder() {
        int result = fcinit();
        check_exit(FC_OK == result, coder_ns::CODER_ERR_INNER, "fcinit failed (%d) : %d.", __LINE__, result);
    }
}

void* safe_realloc_init(uint32_t& size, uint8_t* ptr, size_t new_size, char ch) {
    if (realloc_proc == nullptr) {
        return ptr;
    }

    if (new_size > size) {
        uint8_t* temp_ptr = realloc_proc((size_t&)size, ptr, new_size);
        if (temp_ptr == nullptr) {
            return ptr;
        }

        ptr = temp_ptr;
        memset(ptr + size, ch, static_cast<size_t>(new_size - size));
        size = new_size;
    }

    return ptr;
}

void* safe_realloc(uint32_t& size, uint8_t* ptr, size_t new_size) {
    return safe_realloc_init(size, ptr, new_size, 0);
}


void* safe_alloc(size_t size) {
    if (alloc_proc == nullptr) {
        return nullptr;
    }

    return alloc_proc(size);
}

void* safe_alloc_init(size_t size, char ch) {
    if (alloc_proc == nullptr) {
        return nullptr;
    }

    void* temp_ptr =  alloc_proc(size);
    if (temp_ptr == nullptr) {
        return nullptr;
    }

    memset(temp_ptr, ch, size);
    return temp_ptr;
}

void* safe_calloc(size_t num, size_t size) {
    return safe_alloc_init(num * size, 0);
}

void safe_free(void** ptr) {
    if (free_proc != nullptr) {
        return free_proc(*ptr);
    }
    return;
}

void coder_logger(coder_ns::coder_log_level level, const char* log_format, ...) {
    char log_message[2048];
    va_list args;
    va_start(args, log_format);
    vsnprintf(log_message, sizeof(log_message), log_format, args);
    va_end(args);

    if (logger_proc != nullptr) {
        return logger_proc(level, log_message);
    }

    fprintf(stderr, "%s\n", log_message);
    return;
}

/*
 * 不再调用 exit_proc / exit()：错误改为异常上抛，由调用方（PbgzEngine 的单块
 * 处理边界）捕获后走正常的失败返回路径，保证每一个错误最终都有人处理。
 */
void coder_exit(int16_t exit_code, const char* exit_msg_format, ...) {
    char exit_message[2048];
    va_list args;
    va_start(args, exit_msg_format);
    vsnprintf(exit_message, sizeof(exit_message), exit_msg_format, args);
    va_end(args);

    coder_logger(coder_ns::ERROR, "%s", exit_message);
    throw coder_exception(exit_code, exit_message);
}

void check_exit(bool condition, int16_t exit_code, const char* exit_msg_format, ...) {
    if (condition) {
        return;
    }

    char exit_message[2048];
    va_list args;
    va_start(args, exit_msg_format);
    vsnprintf(exit_message, sizeof(exit_message), exit_msg_format, args);
    va_end(args);

    coder_logger(coder_ns::ERROR, "%s", exit_message);
    throw coder_exception(exit_code, exit_message);
}
