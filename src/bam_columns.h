/*
 * bam_columns.h - column-oriented BAM record storage
 *
 * The classic path turns every BAM record into a SAM text line (BamBlockReader
 * ::parseBamRecord) and hands the text to SamCodecActuator. That inflates the
 * data by ~2.2x (2.95 GB of BAM became 6.39 GB of SAM text on the measured
 * long-read file) and forces every numeric column through its decimal text
 * representation.
 *
 * BamColumns keeps the fields the BAM record already stores: fixed-width
 * scalars stay scalars (FLAG uint16, POS int32, MAPQ uint8 ...), bases stay
 * 4-bit packed, and variable-length fields are concatenated with a per-record
 * length table. The coders are unchanged - each column is still encoded by the
 * coder that preprocessing picked for that SAM field - only the bytes handed to
 * them are structured instead of textual.
 */

#pragma once

#include <stdint.h>
#include <string.h>
#include <vector>

struct BamColumns {
    uint32_t nRecords = 0;
    /* False until a worker has parsed the raw records in the block buffer into
       the payload columns below. In fast mode the reader thread only streams
       the raw BAM record bytes ([u32 size][record]...), and each coder thread
       parses its own block here - so record->column work runs in parallel
       across blocks instead of serially on the reader. */
    bool columnsBuilt = false;
    /* Byte offset in the block buffer where the [u32 size][record] stream
       starts. Zero for the usual split-header layout; non-zero only when a SAM
       header text is prepended to the first data block (!splitHeader). */
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

/* Record/record-stream parsers used by the worker that builds the columns of
   one block (BamCodecActuator in fast mode). They were moved out of
   BamBlockReader so the reader thread never parses: it only streams raw
   records, and parsing is fanned out over the coder threads, one block each. */
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

/* Parse one BAM alignment record (data..data+size) into cols. Columns keep
   their on-disk values: bases are expanded to characters, quality is raw
   phred with 0xFF rendered as '*' (same conversion as the SAM text path). */
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

/* Parse a whole block whose buffer is a sequence of [u32 recordSize][record]
   entries (exactly what BamBlockReader appends for structured BAM blocks).
   On success fills cols and marks it built; returns 0, -1 on a corrupt
   stream. */
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
