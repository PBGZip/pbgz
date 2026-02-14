/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Statistical data compression model for QLFC               */
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

#include <stdlib.h>
#include <memory.h>

#include "fc_model.h"
#include "fc_header.h"
#include <shared_mutex>
#include <mutex>

fcModel1 gfcModel1;
fcModel2 gfcModel2;

static const int16_t tid_state_model0[256] =
{
    293,   434,  2047,   186,    78,   119,   791,  2028,  2058,  1023,  2049,    86,    38,  2048,  2048,  4018,
   2048,   128,    72,  2048,    18,  2048,  2048,    19,     6,   311,  2715,  2048,    28,  2048,  2048,   198,
   2048,  2048,  2048,    26,  2048,  2048,    22,  2048,  2048,    27,    33,   738,   177,  2048,    45,   807,
   2048,   510,    13,  1914,  2048,    20,   845,  2048,   358,  2048,  2048,  2048,    34,  2048,   617,  4058,
   2048,     9,  2048,  2039,  2048,  2048,  2048,  2048,    36,  2048,   764,   156,    12,  2048,  2048,  2048,
   2048,  2048,  3961,  2048,  2048,  2048,    97,  2048,  2048,   102,  2048,    29,    21,  2048,  2048,  2048,
   2048,  2048,   658,   151,  2048,  2048,  2048,  1353,  2048,  2048,  2048,   950,    14,    21,   561,    32,
     16,  2048,  2048,  2048,    24,  2048,  2048,  1443,    23,  2048,    64,  1042,   207,  2048,    44,  2048,
    105,  2048,  2048,  2048,  2048,  2048,  1209,  2048,  2048,    23,  2048,  2048,   997,  1282,  2048,  2048,
    137,  2048,  2048,  2048,  2048,  2048,  3248,  2048,  2048,  2048,  2048,  2048,   266,    17,  2048,  1583,
   2048,   123,  2048,    41,  2048,    18,     1,    18,    22,  2048,   113,  2048,  1089,  1157,  2048,  2048,
    237,  2048,  2048,  2048,    82,  2329,  2048,    33,  2048,  2048,  2048,    14,  2048,  2048,    17,  2048,
     57,  2048,    16,    47,  2048,  2048,    35,  2048,  2048,  2048,  2048,   108,   227,  1991,  2048,    92,
    116,    14,  2048,  2048,  2048,  2048,  1503,  2048,  2048,  2048,    51,    21,    30,   165,   922,  2048,
   4057,  2048,  2048,  2048,  2048,    26,  2048,    22,  2048,  2048,   970,  2048,  2048,  2048,  1789,  2048,
   2048,  2048,  2048,  2048,    11,  2048,  1054,  2048,  2048,   143,  1660,  2048,  2048,    15,  2048,   897,
};

static const int16_t run_state_model0[256] = 
{
   1945,  2048,  3480,  2723,  2048,  2065,  2060,  1845,  2585,  2068,  4008,  2048,  2048,  1750,   959,  3593,
   2048,  3087,   400,  4071,  1583,  1221,  2048,  2048,  2274,  2048,  2048,  2076,  2747,  2048,  2226,  2048,
   2110,  2285,  3946,  2943,  3258,  2373,  2592,  2023,  2048,  2046,  2301,  2794,  2080,  2056,  2048,  2043,
   2894,  3052,  3928,  1642,  1558,  2048,  3031,  2640,  2756,  2048,  1681,  2049,  2604,  2624,    35,  1077,
   3305,  2048,  2265,  2048,  2048,  2690,  2033,  2048,  2347,  2308,  2048,  2165,  2048,  3040,  2048,  2087,
   1885,  2048,  2044,  3181,  2048,  2328,  2048,  3074,  3988,  1836,  2048,  2048,  2048,  2048,  3835,  2048,
   2296,  2058,  2053,  2048,  1704,  2554,  2105,  2506,   525,  3115,  2048,  2802,  2157,  2173,  3388,  2088,
   2048,  2048,  2195,  2280,  2256,  2920,  2039,  2097,  2083,  2062,  2048,  1535,  2184,   162,  2048,  2048,
   2048,  2048,  2236,  2048,  2048,  2779,  2567,  1953,  2048,  2914,  2048,  2838,  2048,  2048,  1982,  2048,
   2656,  2130,  2143,  2213,   640,  2523,  3145,  1914,  2048,  4051,  2499,  2115,  2057,   415,  2864,  2242,
   3717,  2254,  2048,  4033,  2048,  3214,  2029,  2048,  1897,  3280,  2151,  2048,   661,  2048,  2048,  2968,
   2048,  2250,  2048,  2059,  2048,     2,  2048,  2048,  2544,  1932,  3272,  2806,  2048,  2147,  2259,  2467,
    843,  2048,  2813,  1967,  2283,  2426,  2048,  2289,  2048,  1873,  2048,  2048,  2048,  3350,  1143,  2048,
   2048,  2048,  2048,  2048,  2112,  2058,  2048,  2048,  2057,  2048,  2048,  2100,  2048,  2672,  2048,  2048,
   2048,  2048,  2005,  2060,  2048,  2048,  2048,  1358,  2048,  1414,  1815,  2069,  2048,  2048,  2048,  2094,
   2071,  2394,  2048,  3006,  2287,  2048,  2048,  2063,  2048,  1215,  2048,  2048,  2139,  2048,  2048,  1717,
};

void fcmemset(void * dst, int size, short v)
{
    for (int i = 0; i < size / 2; ++i) ((short *)dst)[i] = v;
}

static std::shared_mutex model_mutex;

int fcinit_model()
{
    std::unique_lock<std::shared_mutex> lock(model_mutex);

    for (int mixer = 0; mixer < CHAR_SIZE; ++mixer)
    {
        gfcModel1.tid_mixer[mixer].Init();
        gfcModel1.tid_mixerout_of_range[mixer].Init();
        gfcModel1.run_mixer[mixer].Init();
    }
    for (int bit = 0; bit < 8; ++bit)
    {
        gfcModel1.tid_mixerbits_value[bit].Init();
        for (int context = 0; context < 8; ++context)
            gfcModel1.tid_mixerbit_width[context][bit].Init();
    }
    for (int bit = 0; bit < 32; ++bit)
    {
        gfcModel1.run_mixerbits_value[bit].Init();
        for (int context = 0; context < 32; ++context)
            gfcModel1.run_mixerbit_width[context][bit].Init();
    }

    fcmemset(&gfcModel1.tid_t, sizeof(gfcModel1.tid_t), 2048);
    fcmemset(&gfcModel1.run_t, sizeof(gfcModel1.run_t), 2048);

    for (int s = 0; s < CHAR_SIZE; ++s)
    {
        gfcModel1.tid_t.state_model[s] = (short)tid_state_model0[s];
        gfcModel1.run_t.state_model[s] = (short)run_state_model0[s];
    }

    fcmemset(&gfcModel2.tid_t, sizeof(gfcModel2.tid_t), 4096);
    fcmemset(&gfcModel2.run_t, sizeof(gfcModel2.run_t), 1024);

    return FC_OK;
}

void fcinit_model(fcModel1 * model)
{
    std::shared_lock<std::shared_mutex> lock(model_mutex);
    memcpy(model, &gfcModel1, sizeof(fcModel1));
}

void fcinit_model(fcModel2 * model)
{
    std::shared_lock<std::shared_mutex> lock(model_mutex);
    memcpy(model, &gfcModel2, sizeof(fcModel2));
}
