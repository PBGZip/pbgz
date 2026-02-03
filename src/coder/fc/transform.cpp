/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Burrows Wheeler Transform                                 */
/*-----------------------------------------------------------*/

/*--

This file is a part of bsc and/or libbsc, a program and a library for
lossless, block-sorting data compression.

   Copyright (c) 2009-2021 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information and file AUTHORS
for full list of contributors.

See also the bsc and libbsc web site:
  http://libbsc.com/ for more information.

--*/
#include <cstdlib>

#include "fc_header.h"
#include "libsais.h"
#include "transform.h"

int fc_transform(unsigned char * T, int n, unsigned char * num_indexes, int * indexes)
{
    if (int * RESTRICT A = (int *)malloc(n * sizeof(int)))
    {
        int mod = n / 8;
        {
            mod |= mod >> 1;  mod |= mod >> 2;
            mod |= mod >> 4;  mod |= mod >> 8;
            mod |= mod >> 16; mod >>= 1;
        }

#ifdef FC_OPENMP
        check_exit(false, ERR_INTERNEL, "undefined...");
        int index = optbwt_aux_omp(T, T, A, n, 0, mod + 1, indexes, (features & FC_FEATURE_MULTITHREADING) > 0 ? 0 : 1);
#else
        int index = transform_do(T, T, A, n, 0, mod + 1, indexes);
#endif

        free(A);

        switch (index)
        {
            case -1 : return FC_BAD_ARGS;
            case -2 : return FC_LACK_OF_MEMORY;
        }

        num_indexes[0] = (unsigned char)((n - 1) / (mod + 1));
        index = indexes[0]; for (int t = 0; t < num_indexes[0]; ++t) indexes[t] = indexes[t + 1] - 1;

        return index;
    }
    return FC_LACK_OF_MEMORY;
}

int fc_untransform(unsigned char * T, int n, int index, unsigned char num_indexes, int * indexes)
{
    if ((T == NULL) || (n < 0) || (index <= 0) || (index > n))
    {
        return FC_BAD_ARGS;
    }
    if (n <= 1)
    {
        return FC_OK;
    }
    if (int * P = (int *)malloc((n + 1) * sizeof(int)))
    {
        int mod = n / 8;
        {
            mod |= mod >> 1;  mod |= mod >> 2;
            mod |= mod >> 4;  mod |= mod >> 8;
            mod |= mod >> 16; mod >>= 1;
        }

        if (num_indexes == (unsigned char)((n - 1) / (mod + 1)) && indexes != NULL)
        {
            int I[256]; I[0] = index; for (int t = 0; t < num_indexes; ++t) { I[t + 1] = indexes[t] + 1; }

#ifdef FC_OPENMP
            index = untransform_do_omp(T, T, P, n, mod + 1, I, (features & FC_FEATURE_MULTITHREADING) > 0 ? 0 : 1);
#else
            index = untransform_do(T, T, P, n, mod + 1, I);
#endif
        }
        else
        {
#ifdef FC_OPENMP
            index = untransform_omp(T, T, P, n, index, (features & FC_FEATURE_MULTITHREADING) > 0 ? 0 : 1);
#else
            index = untransform(T, T, P, n, index);
#endif
        }

        free(P);

        switch (index)
        {
            case -1 : return FC_BAD_ARGS;
            case -2 : return FC_LACK_OF_MEMORY;
        }       

        return FC_OK;
    };
    return FC_LACK_OF_MEMORY;
}

