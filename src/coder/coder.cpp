#include "coder.h"

coder_exit_func exit_proc = nullptr;
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

    void register_exit_proc(coder_exit_func proc) {
        exit_proc = proc;
        return;
    }

    void resister_logger_proc(coder_logger_func proc) {
        logger_proc = proc;
    }

    void register_free_func(coder_free_func proc) {
        free_proc = proc;
    }
}

void* safe_realloc_init(int32_t& size, uint8_t* ptr, size_t new_size, char ch) {
    if (realloc_proc == nullptr) {
        return ptr;
    }

    if (new_size > size) {
        uint8_t* temp_ptr = alloc_proc(new_size);
        if (temp_ptr == nullptr) {
            return ptr;
        }

        ptr = temp_ptr;
        memset(ptr + size, ch, static_cast<size_t>(new_size - size));
        size = new_size;
    }

    return ptr;
}

void* safe_realloc(int32_t& size, uint8_t* ptr, size_t new_size) {
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

// The following functions will trigger process exit, not recommended for use within the coder library
void coder_exit(int16_t exit_code, const char* exit_msg_format, ...) {
    char exit_message[2048];
    va_list args;
    va_start(args, exit_msg_format);
    vsnprintf(exit_message, sizeof(exit_message), exit_msg_format, args);
    va_end(args);
    
    if (exit_proc != nullptr) {
        return exit_proc(exit_code, exit_message);
    }

    fprintf(stderr, "%s\n", exit_message);
    return exit(exit_code);
}

void check_exit(bool condition, int16_t exit_code, const char* exit_msg_format, ...) {
    if (!condition) {
        char exit_message[2048];
        va_list args;
        va_start(args, exit_msg_format);
        vsnprintf(exit_message, sizeof(exit_message), exit_msg_format, args);
        va_end(args);
        
        if (exit_proc != nullptr) {
            return exit_proc(exit_code, exit_message);
        }

        fprintf(stderr, "%s\n", exit_message);
        return exit(exit_code);
    }

    return;
}


