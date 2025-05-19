/*-----------------------------------------------------------*/
/* Block Sorting, Lossless Data Compression Library.         */
/* Range coder                                               */
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

#ifndef _FC_CODER_RANGECODER_H
#define _FC_CODER_RANGECODER_H

#include <memory.h>
#include "fc_header.h"

class RangeCoder
{

private:

    union rc_data
    {
        struct u
        {
            unsigned int low32;
            unsigned int carry;
        } u;
        unsigned long long low;
    } rc_data;

    unsigned int rc_data_code;
    unsigned int rc_data_ffnum;
    unsigned int rc_data_cache;
    unsigned int rc_data_range;

    const unsigned short * RESTRICT rc_data_input;
          unsigned short * RESTRICT rc_data_output;
          unsigned short * RESTRICT rc_data_outputEOB;
          unsigned short * RESTRICT rc_data_outputStart;

    INLINE void OutputShort(unsigned short s)
    {
#if defined(FC_UNALIGN)
        *rc_data_output++ = s;
#else
        memcpy(rc_data_output++, &s, sizeof(unsigned short));
#endif
    };

    INLINE unsigned short InputShort()
    {
#if defined(FC_UNALIGN)
        return *rc_data_input++;
#else
        unsigned short ret;
        memcpy(&ret, rc_data_input++, sizeof(unsigned short));
        return ret;
#endif
    };

    NOINLINE unsigned int lshift()
    {
        if (rc_data.u.low32 < 0xffff0000U || rc_data.u.carry)
        {
            OutputShort(rc_data_cache + rc_data.u.carry);
            if (rc_data_ffnum)
            {
                unsigned short s = rc_data.u.carry - 1;
                do { OutputShort(s); } while (--rc_data_ffnum);
            }
            rc_data_cache = rc_data.u.low32 >> 16; rc_data.u.carry = 0;
        } else rc_data_ffnum++;
        rc_data.u.low32 <<= 16;

        return rc_data_range << 16;
    }

public:

    INLINE bool is_end()
    {
        return rc_data_output >= rc_data_outputEOB;
    }

    INLINE void InitEncoder(unsigned char * output, int outputSize)
    {
        rc_data_outputStart = (unsigned short *)output;
        rc_data_output      = (unsigned short *)output;
        rc_data_outputEOB   = (unsigned short *)(output + outputSize - 16);
        rc_data.low         = 0;
        rc_data_ffnum       = 0;
        rc_data_cache       = 0;
        rc_data_range       = 0xffffffff;
    };

    INLINE int FinishEncoder()
    {
        if (rc_data_range < 0x10000)
        {
            lshift();
        }

        lshift(); lshift(); lshift();
        return (int)(rc_data_output - rc_data_outputStart) * sizeof(rc_data_output[0]);
    }

    template <int P = 12> INLINE void EncodeBit0(int probability)
    {
        if (rc_data_range < 0x10000)
        {
            rc_data_range = lshift();
        }

        rc_data_range = (rc_data_range >> P) * probability;
    }

    template <int P = 12> INLINE void EncodeBit1(int probability)
    {
        if (rc_data_range < 0x10000)
        {
            rc_data_range = lshift();
        }

        unsigned int range = (rc_data_range >> P) * probability;
        rc_data.low += range; rc_data_range -= range;
    }

    template <int P = 12> INLINE void EncodeBit(unsigned int bit, int probability)
    {
        if (rc_data_range < 0x10000)
        {
            rc_data_range = lshift();
        }

        unsigned int range = (rc_data_range >> P) * probability;

        rc_data.low   = rc_data.low + ((~bit + 1u) & range);
        rc_data_range = range   + ((~bit + 1u) & (rc_data_range - range - range));
    }

    INLINE void EncodeBit(unsigned int bit)
    {
        if (bit) EncodeBit1(2048); else EncodeBit0(2048);
    };

    INLINE void EncodeByte(unsigned int byte)
    {
        for (int bit = 7; bit >= 0; --bit)
        {
            EncodeBit(byte & (1 << bit));
        }
    };

    INLINE void EncodeWord(unsigned int word)
    {
        for (int bit = 31; bit >= 0; --bit)
        {
            EncodeBit(word & (1 << bit));
        }
    };

    INLINE void InitDecoder(const unsigned char * input)
    {
        rc_data_input = (unsigned short *)input;
        rc_data_code  = 0;
        rc_data_range = 0xffffffff;
        rc_data_code  = (rc_data_code << 16) | InputShort();
        rc_data_code  = (rc_data_code << 16) | InputShort();
        rc_data_code  = (rc_data_code << 16) | InputShort();
    };

    template <int P = 12> INLINE int PeakBit(int probability)
    {
        if (rc_data_range < 0x10000)
        {
            rc_data_range <<= 16; rc_data_code = (rc_data_code << 16) | InputShort();
        }

        return rc_data_code >= (rc_data_range >> P) * probability;
    }

    template <int P = 12> INLINE int DecodeBit(int probability)
    {
        if (rc_data_range < 0x10000)
        {
            rc_data_range <<= 16; rc_data_code = (rc_data_code << 16) | InputShort();
        }

        unsigned int range = (rc_data_range >> P) * probability;
        int bit = rc_data_code >= range;

        rc_data_range = bit ? rc_data_range - range : range;
        rc_data_code  = bit ? rc_data_code  - range : rc_data_code;

        return bit;
    }

    template <int P = 12> INLINE void DecodeBit0(int probability)
    {
        rc_data_range = (rc_data_range >> P) * probability;
    }

    template <int P = 12> INLINE void DecodeBit1(int probability)
    {
        unsigned int range = (rc_data_range >> P) * probability;
        rc_data_code -= range; rc_data_range -= range;
    }

    INLINE unsigned int DecodeBit()
    {
        return DecodeBit(2048);
    }

    INLINE unsigned int DecodeByte()
    {
        unsigned int byte = 0;
        for (int bit = 7; bit >= 0; --bit)
        {
            byte += byte + DecodeBit();
        }
        return byte;
    }

    INLINE unsigned int DecodeWord()
    {
        unsigned int word = 0;
        for (int bit = 31; bit >= 0; --bit)
        {
            word += word + DecodeBit();
        }
        return word;
    }
};

#endif

/*-----------------------------------------------------------*/
/* End                                          rangecoder.h */
/*-----------------------------------------------------------*/
