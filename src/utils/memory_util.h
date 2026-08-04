/*
 * memory_util.h - Header file for pbgz project
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
#include <initializer_list>
#include <any>
#include <log/logger.h>

namespace MemoryUtil {
    template<typename T>
    T* safeAlloc(size_t size) {
        do {
            T* ptr = nullptr;
            if (size > 0) {
                ptr = static_cast<T*>(calloc(size, sizeof(T)));
                if (ptr == nullptr) {
                    LOG_ERROR("Memory calloc failed");
                }
            } 
            return ptr;
        } while(0);
    }

    template <typename T>
    T* safeAlloInit(size_t size, char ch) {
        do {
            T* ptr = nullptr;
            if (size> 0) {
                ptr = static_cast<T*>(malloc(sizeof(T) * size));
                if (nullptr != ptr) {
                    memset(ptr, ch, size_t(size * sizeof(T)));
                } else {
                    LOG_ERROR("Memory calloc failed");
                }
            }
            return ptr;
        } while(0);
    }

    template<typename T>
    T* safeRealloc(size_t& size, T* ptr, size_t newSize) {
        do {
            T* tempPtr = nullptr;
            if (newSize > (size_t)size) {
                tempPtr = static_cast<T*>(realloc(ptr, newSize * sizeof(T)));
                if (tempPtr != nullptr) {
                    ptr = tempPtr;
                    memset(ptr + size, 0, (newSize - size) * sizeof(T));
                    size = newSize;
                } else {
                    LOG_ERROR("Memory realloc failed");
                }
            }

            return ptr;
        } while(0);
    }

    template <typename T>
    T* safeReallocInit(size_t& size, T* ptr, int newSize, char ch) {
        do {
            T* tempPtr = nullptr;
            if (newSize > size) {
                tempPtr = static_cast<T*>(realloc(ptr, newSize * sizeof(T)));
                if (tempPtr != nullptr) {
                    ptr = tempPtr;
                    memset(ptr + size, ch, (newSize - size) * sizeof(T));
                    size = newSize;
                } else {
                    LOG_ERROR("Memory realloc failed");
                }
            }

            return ptr;
        } while(0);
    }

    template <typename T>
    void safeFree(T*& ptr) {
        if (ptr != nullptr) {
            free(ptr);
            ptr = nullptr;
        }
    }

    template <typename T>
    T* safeNew(int size) {
        do {
            T* ptr = nullptr;
            if (size > 0) {
                ptr = new T[size];
                if (nullptr == ptr) {
                    LOG_ERROR("Memory not enough");
                }
            }
            return ptr;
        } while(0);
    }

    template <typename T>
    void safeDelete(T*& ptr) {
        do {
            if (ptr) {
                delete [] ptr;
                ptr = nullptr;
            }
        } while(0);
    }

    template <typename T,  typename... Args>
    T* safeNewClass(Args&&... args) {
        do {
            T* ptr = new T(std::forward<Args>(args)...);
            if (nullptr == ptr) {
                LOG_ERROR("Memory not enough");
            }
            return ptr;
        } while(0);
    }

    template <typename T>
    void safeDeleteClass(T*& ptr) {
        do {
            if (ptr) {
                delete ptr;
                ptr = nullptr;
            }
        } while(0);
    }
}
