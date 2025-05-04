#ifndef PBGZ_PLATFORM_COMPAT_H_
#define PBGZ_PLATFORM_COMPAT_H_

#include <stdio.h>
#include <sys/types.h>

// Platform detection
#if defined(__APPLE__) || defined(__MACH__)
#define PLATFORM_MACOS 1
#elif defined(__linux__) || defined(__linux) || defined(linux)
#define PLATFORM_LINUX 1
#elif defined(_WIN32) || defined(__WIN32__) || defined(WIN32)
#define PLATFORM_WINDOWS 1
#endif

// File I/O compatibility for large file support
#if defined(PLATFORM_MACOS)
// On macOS, fseeko/ftello already handle large files
#define fseeko64 fseeko
#define ftello64 ftello
#define fopen64 fopen
#elif defined(PLATFORM_WINDOWS)
// Windows specific handling would go here
#define fseeko64 _fseeki64
#define ftello64 _ftelli64
#define fopen64 fopen
#endif

#endif // PBGZ_PLATFORM_COMPAT_H_