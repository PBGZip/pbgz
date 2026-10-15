/*
 * seq_stream_util.h - byte-level helpers the SEQ columns of both formats share
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
 * Byte-level helpers shared by the SEQ columns of both formats: positions and lengths as base-128
 * varints, and the reference match stream either whole or split into zero runs plus values. They
 * know nothing about records, references or CIGARs. The varint helpers are also used by the SAM
 * field layouts and the QNAME column.
 */
#ifndef _SEQ_STREAM_UTIL_H_
#define _SEQ_STREAM_UTIL_H_

#include <cstdint>
#include <memory>
#include <vector>

/* LEB128 varint and zigzag helpers for the position, length and numeric sub-streams; zigzag keeps
   small negative deltas small by moving the sign into the low bit. */
static inline uint32_t writeVarint(uint8_t* p, uint32_t v) {
    uint32_t n = 0;
    while (v >= 0x80) {
        p[n++] = (uint8_t)(v | 0x80);
        v >>= 7;
    }
    p[n++] = (uint8_t)v;
    return n;
}

static inline uint32_t zigzag32(int32_t v) {
    return (uint32_t)((v << 1) ^ (v >> 31));
}

static inline int32_t unzigzag32(uint32_t v) {
    return (int32_t)((v >> 1) ^ (uint32_t)(-(int32_t)(v & 1)));
}

/* Appends value as a base-128 varint, the encoding the position and length lists use. */
static inline void appendVarint(std::vector<uint8_t>& out, uint32_t value)
{
    while (value >= 0x80) {
        out.push_back((uint8_t)(value | 0x80));
        value >>= 7;
    }
    out.push_back((uint8_t)value);
}

/* Same, for the signed values whose magnitude can exceed 32 bits once zigzagged. */
static inline void appendVarint64(std::vector<uint8_t>& out, uint64_t value)
{
    while (value >= 0x80) {
        out.push_back((uint8_t)(value | 0x80));
        value >>= 7;
    }
    out.push_back((uint8_t)value);
}

/* Bytes appendVarint64 would write for this value, without writing them. */
static inline uint32_t varintLen64(uint64_t value)
{
    uint32_t n = 1;
    while (value >= 0x80) {
        value >>= 7;
        ++n;
    }
    return n;
}

/* Zigzag: the sign of a signed delta in the low bit, so small negatives stay small. */
static inline uint64_t zigzag64(int64_t value)
{
    return ((uint64_t)value << 1) ^ (uint64_t)(value >> 63);
}

/* The inverse of zigzag64. */
static inline int64_t unzigzag64(uint64_t value)
{
    return (int64_t)((value >> 1) ^ (uint64_t)(-(int64_t)(value & 1)));
}

/* Read one varint from buf[off], advancing off; false on end or longer than five bytes. */
static inline bool readVarint(const uint8_t* buf, uint32_t len, uint32_t& off, uint32_t& value)
{
    uint32_t result = 0;
    int shift = 0;
    for (;;) {
        if (off >= len || shift > 28) {
            return false;
        }
        const uint8_t byte = buf[off++];
        result |= (uint32_t)(byte & 0x7F) << shift;
        shift += 7;
        if ((byte & 0x80) == 0) {
            value = result;
            return true;
        }
    }
}

/* readVarint for values needing more than 32 bits; false on end or longer than ten bytes. */
static inline bool readVarint64(const uint8_t* buf, uint32_t len, uint32_t& off, uint64_t& value)
{
    uint64_t result = 0;
    int shift = 0;
    for (;;) {
        if (off >= len || shift > 63) {
            return false;
        }
        const uint8_t byte = buf[off++];
        result |= (uint64_t)(byte & 0x7F) << shift;
        shift += 7;
        if ((byte & 0x80) == 0) {
            value = result;
            return true;
        }
    }
}

/* The two halves of a match stream split into runs and values: the varint run lengths of the zero
   runs, and the surviving non-zero values one byte each (the "m" and "mval" sub-streams). */
struct SeqRleSplit {
    bool useRle = false;
    std::unique_ptr<uint8_t[]> run;   /* varint run lengths of the zero runs */
    std::unique_ptr<uint8_t[]> val;   /* the surviving non-zero values, one byte each */
    uint32_t runLength = 0;
    uint32_t valLength = 0;
    uint32_t nonZero = 0;             /* how many bytes were not zero, i.e. what the threshold saw */
};

/* Whether the match stream is split, and the split itself when it is. The verdict is per block:
   RLE is taken when at least 98% of the bytes are zero. The decoder dispatches on whether the "m"
   sub-stream's meta carries the "rle" member. */
static inline SeqRleSplit splitSeqMatchStream(const uint8_t* match, uint32_t matchLen)
{
    SeqRleSplit split;
    uint32_t nNonZero = 0;
    for (uint32_t i = 0; i < matchLen; i++) {
        if (match[i] != 0) {
            nNonZero++;
        }
    }
    const uint32_t zeroCount = matchLen - nNonZero;
    split.nonZero = nNonZero;
    split.useRle = (matchLen > 0) &&
                   ((uint64_t)zeroCount * 100ull >= (uint64_t)matchLen * 98ull);
    if (!split.useRle) {
        return split;
    }

    split.run = std::make_unique<uint8_t[]>((nNonZero + 1) * 5 + 16);
    split.val = std::make_unique<uint8_t[]>(nNonZero + 16);
    uint32_t rp = 0;
    uint32_t vp = 0;
    uint32_t run = 0;
    for (uint32_t i = 0; i < matchLen; i++) {
        if (match[i] == 0) {
            run++;
        } else {
            uint32_t v = run;
            while (v >= 0x80) { split.run[rp++] = (uint8_t)(v | 0x80); v >>= 7; }
            split.run[rp++] = (uint8_t)v;
            split.val[vp++] = match[i];
            run = 0;
        }
    }
    uint32_t tail = run;                /* trailing run of zeros */
    while (tail >= 0x80) { split.run[rp++] = (uint8_t)(tail | 0x80); tail >>= 7; }
    split.run[rp++] = (uint8_t)tail;

    split.runLength = rp;               /* length of the run-length segment */
    split.valLength = vp;               /* length of the value segment */
    return split;
}

#endif
