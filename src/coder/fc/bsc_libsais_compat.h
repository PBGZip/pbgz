/*--

bsc_libsais_compat.h - bsc-style Burrows-Wheeler entry points for the FC coder.

The FC coder calls the bsc-style entry points transform_do() / untransform_do()
/ untransform(). These were historically provided by the vendored, bsc-merged
libsais.c (a stripped libsais 2.0 fork). PBGZ now vendors upstream libsais
2.8.7 (Apache-2.0) verbatim in libsais.c/.h and implements the three entry
points here, on top of the public libsais API, with byte-identical output.

SPDX-License-Identifier: MIT

--*/

#ifndef PBGZ_BSC_LIBSAIS_COMPAT_H
#define PBGZ_BSC_LIBSAIS_COMPAT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int32_t transform_do(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t r, int32_t * I);
int32_t untransform_do(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t r, const int32_t * I);
int32_t untransform(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t i);

#ifdef __cplusplus
}
#endif

#endif
