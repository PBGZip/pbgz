/*
 * clr.h - Header file for pbgz project
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

/*
 * Note it is up to the calling code to ensure that no overruns on input and
 * output buffers occur.
 *
 * Call the input() and output() functions to set and query the current
 * buffer locations.
 */
#ifndef _CLR_H_
#define _CLR_H_

#include <stdlib.h>
#include <stdint.h>

#define  DO(n)       for (int _=0; _<n; _++)
#define  TOP       (1<<24)

#define CodecName "CLR"

typedef unsigned char uc;

class RangeCoder {
    uint64_t low;
    uint range, code;

public:

    uc *in_buf;
    uc *out_buf;

    void input(char *in) { out_buf = in_buf = (uc *) in; }

    void output(char *out) { in_buf = out_buf = (uc *) out; }

    char *input(void) { return (char *) in_buf; }

    char *output(void) { return (char *) out_buf; }

    int size_out(void) { return out_buf - in_buf; }

    int size_in(void) { return in_buf - out_buf; }

    void StartEncode(void) {
        low = 0;
        range = (uint) -1;
    }

    void StartDecode(void) {
        low = 0;
        range = (uint) -1;
        DO(8) code = (code << 8) | *in_buf++;
    }

    void FinishEncode(void) {
        DO(8) (*out_buf++ = low >> 56), low <<= 8;
    }

    void FinishDecode(void) {}

    void Encode(uint cumFreq, uint freq, uint totFreq) {
        low += cumFreq * (range /= totFreq);
        range *= freq;

        if (cumFreq + freq > totFreq)
            abort();

        while (range < TOP) {
            // range = 0x00ffffff..
            // low/high may be matching
            //       eg 88332211/88342211 (range 00010000)
            // or differing
            //       eg 88ff2211/89002211 (range 00010000)
            //
            // If the latter, we need to reduce range down
            // such that high=88ffffff.
            // Eg. top-1      == 00ffffff
            //     low|top-1  == 88ffffff
            //     ...-low    == 0000ddee
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint(low) | (TOP - 1)) - uint(low));
            *out_buf++ = low >> 56, range <<= 8, low <<= 8;
        }
    }

    uint GetFreq(uint totFreq) {
        return code / (range /= totFreq);
    }

    void Decode(uint cumFreq, uint freq) {
        uint temp = cumFreq * range;
        low += temp;
        code -= temp;
        range *= freq;

        while (range < TOP) {
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint(low) | (TOP - 1)) - uint(low));
            code = (code << 8) | *in_buf++, range <<= 8, low <<= 8;
        }
    }
};

#endif
