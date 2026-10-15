/*
 * bam_columns.h - column-oriented BAM record storage
 *
 * Keeps BAM fields in structured columns instead of converting records to SAM text:
 * fixed-width scalars stay scalars, bases stay 4-bit packed, and variable-length
 * fields are concatenated with a per-record length table. Each column is encoded by
 * the same coder preprocessing picked for that SAM field.
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

#pragma once

#include <stdint.h>
#include <string.h>
#include <vector>

struct BamColumns {
    uint32_t nRecords = 0;
    /* False until a worker has parsed the raw records into the columns below.
       In fast mode each coder thread parses its own block. */
    bool columnsBuilt = false;
    /* Byte offset where the [u32 size][record] stream starts; non-zero only when a
       SAM header is prepended to the first data block (!splitHeader). */
    uint32_t recordStreamOffset = 0;

    /* Fixed-width columns, one entry per record. */
    std::vector<uint16_t> flag;
    std::vector<int32_t>  refId;
    std::vector<int32_t>  pos;        /* 0-based, as stored in BAM */
    std::vector<uint8_t>  mapq;
    std::vector<int32_t>  nextRefId;
    std::vector<int32_t>  nextPos;
    std::vector<int32_t>  tlen;

    /* Variable-length columns: payload + per-record length. */
    std::vector<uint8_t>  qname;      /* read names, no NUL terminator */
    std::vector<uint32_t> qnameLen;

    std::vector<uint32_t> cigar;      /* BAM cigar ops, (op<<4)|len */
    std::vector<uint32_t> cigarLen;   /* ops per record */

    std::vector<uint8_t>  seq;        /* 4-bit packed bases, (lSeq+1)/2 bytes */
    std::vector<uint32_t> seqLen;     /* packed bytes per record */
    std::vector<int32_t>  lSeq;       /* base count per record */

    std::vector<uint8_t>  qual;       /* raw quality bytes, lSeq per record */

    std::vector<uint8_t>  tags;       /* raw aux bytes */
    std::vector<uint32_t> tagsLen;

    void clear()
    {
        nRecords = 0;
        columnsBuilt = false;
        flag.clear();
        refId.clear();
        pos.clear();
        mapq.clear();
        nextRefId.clear();
        nextPos.clear();
        tlen.clear();
        qname.clear();
        qnameLen.clear();
        cigar.clear();
        cigarLen.clear();
        seq.clear();
        seqLen.clear();
        lSeq.clear();
        qual.clear();
        tags.clear();
        tagsLen.clear();
    }
};

/* Record/record-stream parsers used by the worker that builds one block's columns
   (BamCodecActuator in fast mode). */
namespace bamcol {

inline const char kBaseMap[16] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V',
                                  'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

inline uint16_t rU16(const uint8_t*& p) {
    uint16_t v;
    memcpy(&v, p, 2);
    p += 2;
    return v;
}
inline int32_t rI32(const uint8_t*& p) {
    int32_t v;
    memcpy(&v, p, 4);
    p += 4;
    return v;
}

/* Parse one BAM alignment record (data..data+size) into cols. Bases are expanded to
   characters; quality is raw phred with 0xFF rendered as '*'. */
inline int32_t parseRecord(const uint8_t* data, int32_t size, BamColumns& cols)
{
    const uint8_t* p = data;
    const uint8_t* end = data + size;

    const int32_t refID = rI32(p);
    const int32_t pos = rI32(p);
    const uint8_t lReadName = *p++;
    const uint8_t mapq = *p++;
    (void)rU16(p);              /* bin */
    const uint16_t nCigarOp = rU16(p);
    const uint16_t flag = rU16(p);
    const int32_t lSeq = rI32(p);
    const int32_t nextRefID = rI32(p);
    const int32_t nextPos = rI32(p);
    const int32_t tlen = rI32(p);

    if (lReadName == 0 || p > end) {
        return -1;
    }

    cols.flag.push_back(flag);
    cols.refId.push_back(refID);
    cols.pos.push_back(pos);
    cols.mapq.push_back(mapq);
    cols.nextRefId.push_back(nextRefID);
    cols.nextPos.push_back(nextPos);
    cols.tlen.push_back(tlen);

    const uint32_t nameLen = (uint32_t)lReadName - 1u;
    cols.qname.insert(cols.qname.end(), p, p + nameLen);
    cols.qnameLen.push_back(nameLen);
    p += lReadName;

    const uint8_t* cigarData = p;
    p += (size_t)nCigarOp * 4;
    cols.cigarLen.push_back((uint32_t)nCigarOp);
    for (uint16_t i = 0; i < nCigarOp; ++i) {
        const uint8_t* op = cigarData + (size_t)i * 4;   /* BAM is little-endian */
        cols.cigar.push_back((uint32_t)op[0] | ((uint32_t)op[1] << 8) |
                             ((uint32_t)op[2] << 16) | ((uint32_t)op[3] << 24));
    }

    const size_t packedLen = (size_t)((lSeq > 0 ? lSeq : 0) + 1) / 2;
    const size_t seqBase = cols.seq.size();
    cols.seq.resize(seqBase + (size_t)(lSeq > 0 ? lSeq : 0));
    for (int32_t b = 0; b < lSeq; ++b) {
        const uint8_t byte = p[(size_t)(b >> 1)];
        const uint8_t nib = (b & 1) ? (byte & 0x0f) : (uint8_t)((byte >> 4) & 0x0f);
        cols.seq[seqBase + (size_t)b] = (uint8_t)kBaseMap[nib];
    }
    cols.seqLen.push_back((uint32_t)(lSeq > 0 ? lSeq : 0));
    cols.lSeq.push_back(lSeq);
    p += packedLen;

    const size_t qualLen = (size_t)(lSeq > 0 ? lSeq : 0);
    const size_t qualBase = cols.qual.size();
    cols.qual.resize(qualBase + qualLen);
    for (size_t i = 0; i < qualLen; ++i) {
        cols.qual[qualBase + i] = (p[i] == 0xFF) ? (uint8_t)'*' : (uint8_t)(p[i] + 33);
    }
    p += qualLen;

    if (p > end) {
        return -1;
    }
    cols.tags.insert(cols.tags.end(), p, end);
    cols.tagsLen.push_back((uint32_t)(end - p));
    return 0;
}

/* Parse a block buffer of [u32 recordSize][record] entries. On success fills cols
   and marks it built; returns 0, or -1 on a corrupt stream. */
inline int32_t parseRecordStream(const uint8_t* buf, size_t len, BamColumns& cols)
{
    cols.clear();
    size_t off = 0;
    while (off + 4 <= len) {
        int32_t blockSize;
        memcpy(&blockSize, buf + off, 4);
        if (blockSize <= 0 || off + 4 + (size_t)blockSize > len) {
            return -1;
        }
        if (parseRecord(buf + off + 4, blockSize, cols) != 0) {
            return -1;
        }
        off += 4 + (size_t)blockSize;
    }
    if (off != len) {
        return -1;
    }
    cols.nRecords = (uint32_t)cols.lSeq.size();
    cols.columnsBuilt = true;
    return 0;
}

}  /* namespace bamcol */
