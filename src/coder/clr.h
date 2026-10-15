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

/* Call input()/output() to set and query the buffer locations; the caller must
   ensure no buffer overruns. */
#ifndef _CLR_H_
#define _CLR_H_

#include <stdlib.h>
#include <stdint.h>

#define  DO(n)       for (int _=0; _<n; _++)
#define  TOP       (1<<24)

#define CodecName "CLR"

typedef unsigned char uc;

/* Boundary/error handling: overflow sets a sticky err and stops advancing;
   errors are checked once in FinishEncode/FinishDecode. */
class RangeCoder {
    uint64_t low;
    uint range, code;

public:

    uc *in_buf;
    uc *out_buf;
    uc *in_end;
    uc *out_end;
    int err;

    RangeCoder() : low(0), range(0), code(0),
                   in_buf(nullptr), out_buf(nullptr),
                   in_end(nullptr), out_end(nullptr), err(0) {}

    /* Endpoints are required parameters: without knowing where it ends, the buffer cannot be set up. */
    void input(char *in, char *end) { out_buf = in_buf = (uc *) in; in_end = (uc *) end; }

    void output(char *out, char *end) { in_buf = out_buf = (uc *) out; out_end = (uc *) end; }

    char *input(void) { return (char *) in_buf; }

    char *output(void) { return (char *) out_buf; }

    int size_out(void) { return out_buf - in_buf; }

    int size_in(void) { return in_buf - out_buf; }

    void StartEncode(void) {
        low = 0;
        range = (uint) -1;
        err = 0;
    }

    void StartDecode(void) {
        low = 0;
        range = (uint) -1;
        err = 0;
        if (in_buf + 8 > in_end) {
            err = -1;
            in_buf = in_end;   /* Push to the end so all subsequent decoding stops immediately at the boundary check */
            return;
        }
        DO(8) code = (code << 8) | *in_buf++;
    }

    int FinishEncode(void) {
        DO(8) {
            if (out_buf >= out_end) { err = -1; break; }
            (*out_buf++ = low >> 56), low <<= 8;
        }
        return err;
    }

    int FinishDecode(void) { return err; }

    void Encode(uint cumFreq, uint freq, uint totFreq) {
        if (err) return;

        /* A contradictory model means corrupt state; set the error instead of aborting. */
        if (!totFreq || cumFreq + freq > totFreq) { err = -1; return; }

        low += cumFreq * (range /= totFreq);
        range *= freq;

        while (range < TOP) {
            /* Bring high down to 88ffffff when low and low+range differ in the top byte. */
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint(low) | (TOP - 1)) - uint(low));
            if (out_buf >= out_end) { err = -1; return; }
            *out_buf++ = low >> 56, range <<= 8, low <<= 8;
        }
    }

    uint GetFreq(uint totFreq) {
        /* A corrupt code stream can make totFreq 0 or collapse range; without this guard it would be a SIGFPE. */
        if (err || !totFreq || range < totFreq) { err = -1; return 0; }
        return code / (range /= totFreq);
    }

    void Decode(uint cumFreq, uint freq) {
        if (err) return;

        uint temp = cumFreq * range;
        low += temp;
        code -= temp;
        range *= freq;

        while (range < TOP) {
            if (uc((low ^ (low + range)) >> 56))
                range = ((uint(low) | (TOP - 1)) - uint(low));
            if (in_buf >= in_end) { err = -1; return; }
            code = (code << 8) | *in_buf++, range <<= 8, low <<= 8;
        }
    }
};

#endif
