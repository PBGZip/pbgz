#pragma once
#include <initializer_list>
#include <any>

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
    T* safeRealloc(int& size, T* ptr, size_t newSize) {
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
    T* safeReallocInit(int& size, T* ptr, int newSize, char ch) {
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
                ptr = std::shared_ptr<T>(new T[size], [](T* p) { delete[] p;});
                if (nullptr == ptr) {
                    LOG_ERROR("Memory not enough");
                }
            }
            return ptr;
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
}