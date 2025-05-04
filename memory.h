#ifndef PBGZ_MEMORY_H_
#define PBGZ_MEMORY_H_

#include <cstdlib>
#include <cstring>
#include <memory>
#include <cinttypes>  // For PRIu64 format specifier
#include "manager.h"

/**
 * @brief Safely allocate zero-initialized memory
 */
#define safe_alloc(n, type, p)                                                               \
    do {                                                                                     \
        if ((n) > 0) {                                                                       \
            (p) = static_cast<type*>(calloc((n), sizeof(type)));                             \
            if (!(p)) {                                                                      \
                fprintf(stderr, "Error: Insufficient memory: need %" PRIu64 " MB\n",         \
                       (static_cast<uint64_t>(sizeof(type)) * (n)) >> 20);                   \
                manage::instance().exit(ERR_MEM_NOENOUGH);                                   \
            }                                                                                \
        } else {                                                                             \
            (p) = nullptr;                                                                   \
        }                                                                                    \
    } while (0)

/**
 * @brief Safely allocate memory with custom initialization value
 */
#define safe_alloc_init(n, type, p, x)                                                       \
    do {                                                                                     \
        if ((n) > 0) {                                                                       \
            (p) = static_cast<type*>(malloc(static_cast<size_t>(n) * sizeof(type)));         \
            if (!(p)) {                                                                      \
                fprintf(stderr, "Error: Insufficient memory: need %" PRIu64 " MB\n",         \
                       (static_cast<uint64_t>(sizeof(type)) * (n)) >> 20);                   \
                manage::instance().exit(ERR_MEM_NOENOUGH);                                   \
            }                                                                                \
            memset((p), (x), static_cast<size_t>(n) * sizeof(type));                         \
        } else {                                                                             \
            (p) = nullptr;                                                                   \
        }                                                                                    \
    } while (0)

/**
 * @brief Safely reallocate memory with zero initialization for new portion
 */
#define safe_realloc(n, type, p, N)                                                          \
    do {                                                                                     \
        if ((N) > (n)) {                                                                     \
            type* temp = static_cast<type*>(realloc((p), static_cast<size_t>(N) * sizeof(type))); \
            if (!temp) {                                                                     \
                fprintf(stderr, "Error: Insufficient memory: need %" PRIu64 " MB\n",         \
                       (static_cast<uint64_t>(sizeof(type)) * (N)) >> 20);                   \
                manage::instance().exit(ERR_MEM_NOENOUGH);                                   \
            }                                                                                \
            (p) = temp;                                                                      \
            memset((p) + (n), 0, static_cast<size_t>((N) - (n)) * sizeof(type));             \
            (n) = (N);                                                                       \
        }                                                                                    \
    } while (0)

/**
 * @brief Safely reallocate memory with custom initialization for new portion
 */
#define safe_realloc_init(n, type, p, N, x)                                                  \
    do {                                                                                     \
        if ((N) > (n)) {                                                                     \
            type* temp = static_cast<type*>(realloc((p), static_cast<size_t>(N) * sizeof(type))); \
            if (!temp) {                                                                     \
                fprintf(stderr, "Error: Insufficient memory: need %" PRIu64 " MB\n",         \
                       (static_cast<uint64_t>(sizeof(type)) * (N)) >> 20);                   \
                manage::instance().exit(ERR_MEM_NOENOUGH);                                   \
            }                                                                                \
            (p) = temp;                                                                      \
            memset((p) + (n), (x), static_cast<size_t>((N) - (n)) * sizeof(type));           \
            (n) = (N);                                                                       \
        }                                                                                    \
    } while (0)

/**
 * @brief Safely create a new shared_ptr to array with zero initialization
 */
#define safe_new(n, type, p)                                                                 \
    do {                                                                                     \
        try {                                                                                \
            (p) = std::shared_ptr<type>(new type[(n)], [](type* ptr) { delete[] ptr; });     \
            if ((p)) {                                                                       \
                memset(static_cast<type*>((p).get()), 0, static_cast<size_t>(n) * sizeof(type)); \
            }                                                                                \
        } catch (const std::exception& e) {                                                  \
            fprintf(stderr, "Error: Memory allocation failed: %s (%" PRIu64 " MB)\n",        \
                   e.what(), (static_cast<uint64_t>(sizeof(type)) * (n)) >> 20);             \
            manage::instance().exit(ERR_MEM_NOENOUGH);                                       \
        } catch (...) {                                                                      \
            fprintf(stderr, "Error: Unknown allocation failure: %" PRIu64 " MB\n",           \
                   (static_cast<uint64_t>(sizeof(type)) * (n)) >> 20);                       \
            manage::instance().exit(ERR_MEM_NOENOUGH);                                       \
        }                                                                                    \
    } while (0)

/**
 * @brief Safely create a new class instance with constructor arguments
 * 
 * @param c Class type
 * @param init Constructor arguments
 * @param p Pointer that will hold the new class instance
 */
#define safe_new_class(c, init, p)                                                           \
    do {                                                                                     \
        try {                                                                                \
            (p) = new c(init);                                                               \
        } catch (const std::exception& e) {                                                  \
            fprintf(stderr, "Error: Class instantiation failed: %s\n", e.what());            \
            manage::instance().exit(ERR_MEM_NOENOUGH);                                       \
        } catch (...) {                                                                      \
            fprintf(stderr, "Error: Unknown class instantiation failure\n");                 \
            manage::instance().exit(ERR_MEM_NOENOUGH);                                       \
        }                                                                                    \
    } while (0)

#endif  // PBGZ_MEMORY_H_