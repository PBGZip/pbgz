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
 * The SEQ column of both formats is a byte stream that is written in the same handful of shapes:
 * positions and lengths as base-128 varints, and the reference match stream either whole or split
 * into zero runs and the values that follow them. Those shapes know nothing about records,
 * references or CIGARs, so they live here rather than inside either format's module - the SAM side
 * grew them first and the FASTQ side uses the very same ones, which is what this file exists for.
 *
 * The varint helpers are the more general half: they are what the SAM field layouts and the QNAME
 * column write and read their numbers with too, so they carry plain names despite having grown up
 * next to the TLEN payload that first needed them.
 */
#ifndef _SEQ_STREAM_UTIL_H_
#define _SEQ_STREAM_UTIL_H_

#include <cstdint>
#include <memory>
#include <vector>

/*
 * LEB128 varint and zigzag helpers, used by the position, length and numeric sub-streams (the
 * names date from the TLEN stream, which was the first to replace a fixed-width layout with them).
 * Both replace layouts whose high bytes were almost always 0, and zigzag keeps small negative
 * differences small by moving the sign into the low bit.
 */
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

/* How many bytes appendVarint64 would write for this value, without writing them: the size a
   column would take in a varint layout, for a caller deciding between two of them. */
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

/* Reads one varint from buf[off], advancing off; false when it runs off the end or is longer than
   five bytes. */
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

/* Same as readVarint, for the values that need more than 32 bits (the reference positions of the
   FASTQ SEQ column are deltas of 64-bit coordinates); false when it runs off the end or is longer
   than ten bytes. */
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

/*
 * The two halves of a match stream split into runs and values.
 *
 * The match stream is a sparse 0..3 byte stream in which a block whose reads all take their bases
 * from the reference is ~99% zeros. When the zero share is high enough to pay for it, the stream is
 * replaced by two independent sub-streams: this one (the varint run lengths of the zero runs) and
 * the surviving non-zero values, one byte each. The decision is splitSeqMatchStream's, and its two
 * halves are the "m" and "mval" sub-streams.
 */
struct SeqRleSplit {
    bool useRle = false;
    std::unique_ptr<uint8_t[]> run;   /* varint run lengths of the zero runs */
    std::unique_ptr<uint8_t[]> val;   /* the surviving non-zero values, one byte each */
    uint32_t runLength = 0;
    uint32_t valLength = 0;
    uint32_t nonZero = 0;             /* how many bytes were not zero, i.e. what the threshold saw */
};

/*
 * Whether the match stream is split, and the split itself when it is.
 *
 * The verdict is per block and depends only on the payload: RLE is taken when at least 98% of the
 * bytes are zero. The crossover is sharp and was measured on SAM blocks built with a known fraction
 * of reads taken from the reference: at 99.0% zeros the split wins by 2.3%, at 95.6% it already
 * loses by 12%, and the two kinds of block production actually hands over sit three orders of
 * magnitude away from the line (all reads from the reference: 99.2% zeros; none of them: 30-36%).
 * Below the threshold, splitting is worse than useless: each surviving value then costs one byte of
 * "mval" plus about one of run length, so the coders see 14%-39% more input, and the two halves are
 * no longer adjacent, which costs most where the data repeats across reads.
 *
 * A block whose match stream is not split is written byte-for-byte as it was before the split
 * existed, so nothing here needs a format version: the decoder dispatches on whether the "m"
 * sub-stream's meta carries the "rle" member.
 */
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
