/*--

bsc_libsais_compat.cpp - see bsc_libsais_compat.h.

SPDX-License-Identifier: MIT

--*/

#include "bsc_libsais_compat.h"

#include "libsais.h"

int32_t transform_do(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t r, int32_t * I)
{
    if ((T == nullptr) || (U == nullptr) || (A == nullptr) || (n < 0) || (fs < 0) || (r < 2) || ((r & (r - 1)) != 0) || (I == nullptr))
    {
        return -1;
    }

    if (n <= 1)
    {
        if (n == 1) { U[0] = T[0]; }

        I[0] = n;
        return 0;
    }

    return libsais_bwt_aux(T, U, A, n, fs, nullptr, r, I);
}

int32_t untransform_do(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t r, const int32_t * I)
{
    if ((T == nullptr) || (U == nullptr) || (A == nullptr) || (n < 0) || ((r != n) && ((r < 2) || ((r & (r - 1)) != 0))) || (I == nullptr))
    {
        return -1;
    }

    if (n <= 1)
    {
        if (I[0] != n) { return -1; }
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    return libsais_unbwt_aux(T, U, A, n, nullptr, r, I);
}

int32_t untransform(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t i)
{
    if ((T == nullptr) || (U == nullptr) || (A == nullptr) || (n < 0) || (i <= 0) || (i > n))
    {
        return -1;
    }

    if (n <= 1)
    {
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    return libsais_unbwt(T, U, A, n, nullptr, i);
}
