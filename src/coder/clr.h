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

/*
 * Boundary and error handling follow htscodecs/c_range_coder.h: this class
 * shares its origin (Shelwien's original), but that one was hardened by
 * Bonfield while this one has always been the unhardened version—no end markers,
 * no error flag, reads and writes go straight past the boundary.
 *
 * Three disciplines: the endpoints live on the coder itself, not as a verbal
 * convention for the caller; overflow sets a sticky err and immediately stops
 * advancing; errors are picked up at exactly one place, FinishEncode/
 * FinishDecode. Because err, once latched, stops all further progress, checking
 * once at the end is equivalent to checking every symbol, so the hot path needs
 * no checks at all.
 */
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

        /*
         * This used to be abort(): killing the process deep inside a library
         * function, with no logging and bypassing the engine's error
         * aggregation. Changed to setting an error—a contradictory model means
         * the data or state is already corrupt, and it is left to the surrounding
         * code to decide.
         */
        if (!totFreq || cumFreq + freq > totFreq) { err = -1; return; }

        low += cumFreq * (range /= totFreq);
        range *= freq;

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
