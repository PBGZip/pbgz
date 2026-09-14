/*
 * bam_record.h - SAM text line -> BAM binary record
 *
 * The single place that turns one SAM alignment line into the byte stream a
 * BAM file stores for that record. It sits outside BamWriter, the only
 * writer-thread component allowed to touch the output, so the conversion can run
 * in the decompression workers (in parallel across blocks) instead of every
 * block travelling as SAM text to the writer thread to be re-parsed there, which
 * would serialise text->binary conversion on a single thread.
 *
 * The builder is deliberately free of BGZF framing and of any writer state:
 * it appends a complete record (the leading u32 block_size included) to a
 * caller-owned buffer. Two callers use it:
 *   - BamCodecActuator::decompress builds the records inside the decompression
 *     workers, in parallel across blocks, and marks the block as already-BAM;
 *   - BamWriter keeps using it when it receives SAM text (the historical
 *     path, and the unit tests), so both paths produce identical bytes.
 *
 * BamRecordScratch carries the per-record temporaries so a caller that walks a
 * whole block reuses the same allocations instead of churning one set of
 * vectors per record.
 */

#pragma once

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <map>
#include <string>
#include <vector>

#include "log/logger.h"

namespace bamrec {

/* ---- Little-endian byte stream output ---- */
inline void putU8(std::vector<uint8_t>& out, uint8_t v) {
    out.push_back(v);
}
inline void putU16(std::vector<uint8_t>& out, uint16_t v) {
    out.push_back((uint8_t)(v & 0xFF));
    out.push_back((uint8_t)((v >> 8) & 0xFF));
}
inline void putI32(std::vector<uint8_t>& out, int32_t v) {
    uint32_t u = (uint32_t)v;
    out.push_back((uint8_t)(u & 0xFF));
    out.push_back((uint8_t)((u >> 8) & 0xFF));
    out.push_back((uint8_t)((u >> 16) & 0xFF));
    out.push_back((uint8_t)((u >> 24) & 0xFF));
}
inline void putU32(std::vector<uint8_t>& out, uint32_t v) {
    putI32(out, (int32_t)v);
}

/* ---- BAM CIGAR opcodes (indices match BAM_CIGAR_OPS "MIDNSHP=XB") ---- */
inline int cigarOpCode(char c) {
    switch (c) {
    case 'M': return 0;
    case 'I': return 1;
    case 'D': return 2;
    case 'N': return 3;
    case 'S': return 4;
    case 'H': return 5;
    case 'P': return 6;
    case '=': return 7;
    case 'X': return 8;
    default:  return -1;
    }
}

/* CIGAR operations that consume reference length: M/D/N/=/X */
inline bool cigarConsumesRef(int code) {
    return code == 0 || code == 2 || code == 3 || code == 7 || code == 8;
}

/* Parse a CIGAR string (returns true with no ops for "*" or empty); refSpan is the reference length consumed */
inline bool parseCigar(const char* s, size_t n, std::vector<uint32_t>& ops, int64_t& refSpan) {
    ops.clear();
    refSpan = 0;
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        while (j < n && s[j] >= '0' && s[j] <= '9') {
            ++j;
        }
        if (j == i) {
            return false;   // Missing length
        }
        uint32_t len = 0;
        for (size_t k = i; k < j; ++k) {
            len = len * 10 + (uint32_t)(s[k] - '0');
        }
        if (j >= n) {
            return false;   // Missing operator
        }
        const int code = cigarOpCode(s[j]);
        if (code < 0) {
            return false;
        }
        ops.push_back((len << 4) | (uint32_t)code);
        if (cigarConsumesRef(code)) {
            refSpan += len;
        }
        i = j + 1;
    }
    return true;
}

/* BAM 4-bit base encoding: char -> index (inverse of the read-side BAM_BASE_MAP) */
inline int baseNibble(char c) {
    switch (c) {
    case '=': return 0;
    case 'A': return 1;
    case 'C': return 2;
    case 'M': return 3;
    case 'G': return 4;
    case 'R': return 5;
    case 'S': return 6;
    case 'V': return 7;
    case 'T': return 8;
    case 'W': return 9;
    case 'Y': return 10;
    case 'H': return 11;
    case 'K': return 12;
    case 'D': return 13;
    case 'B': return 14;
    case 'N': return 15;
    default:  return 15;
    }
}

/* Classic reg2bin: beg is the 0-based start, end is the 0-based exclusive end */
inline uint16_t samReg2Bin(int64_t beg, int64_t end) {
    const int64_t e = end - 1;
    if ((beg >> 14) == (e >> 14)) {
        return (uint16_t)(((1 << 15) - 1) / 7 + (beg >> 14));
    }
    if ((beg >> 17) == (e >> 17)) {
        return (uint16_t)(((1 << 12) - 1) / 7 + (beg >> 17));
    }
    if ((beg >> 20) == (e >> 20)) {
        return (uint16_t)(((1 << 9) - 1) / 7 + (beg >> 20));
    }
    if ((beg >> 23) == (e >> 23)) {
        return (uint16_t)(((1 << 6) - 1) / 7 + (beg >> 23));
    }
    if ((beg >> 26) == (e >> 26)) {
        return (uint16_t)(((1 << 3) - 1) / 7 + (beg >> 26));
    }
    return 0;
}

/*
 * Convert one SAM optional field "TAG:TYPE:VALUE" into BAM aux bytes (including TAG and
 * type); returns -1 on failure (the caller skips that field). Named differently from the
 * read-side appendBamAux (BAM->SAM).
 */
inline int appendSamAuxToBam(const char* opt, size_t n, std::vector<uint8_t>& out) {
    size_t p = 0;
    while (p < n && opt[p] != ':') {
        ++p;
    }
    if (p != 2) {
        return -1;   // TAG must be 2 characters
    }
    if (p + 3 > n || opt[p + 2] != ':') {
        return -1;   // Expected "TAG:TYPE:VALUE"
    }
    const char type = opt[p + 1];
    const char* value = opt + p + 3;
    const size_t vlen = n - (p + 3);

    putU8(out, (uint8_t)opt[0]);
    putU8(out, (uint8_t)opt[1]);
    switch (type) {
    case 'A': {
        if (vlen != 1) {
            return -1;
        }
        putU8(out, (uint8_t)'A');
        putU8(out, (uint8_t)value[0]);
        break;
    }
    case 'i': {
        long long v = 0;
        try {
            v = std::stoll(std::string(value, vlen));
        } catch (...) {
            return -1;
        }
        if (v >= INT8_MIN && v <= INT8_MAX) {
            putU8(out, (uint8_t)'c');
            putU8(out, (uint8_t)(int8_t)v);
        } else if (v >= 0 && v <= UINT8_MAX) {
            putU8(out, (uint8_t)'C');
            putU8(out, (uint8_t)v);
        } else if (v >= INT16_MIN && v <= INT16_MAX) {
            putU8(out, (uint8_t)'s');
            putU16(out, (uint16_t)(int16_t)v);
        } else if (v >= 0 && v <= UINT16_MAX) {
            putU8(out, (uint8_t)'S');
            putU16(out, (uint16_t)v);
        } else if (v >= INT32_MIN && v <= INT32_MAX) {
            putU8(out, (uint8_t)'i');
            putI32(out, (int32_t)v);
        } else {
            putU8(out, (uint8_t)'I');
            putU32(out, (uint32_t)v);
        }
        break;
    }
    case 'f': {
        double d = strtod(std::string(value, vlen).c_str(), nullptr);
        float f = (float)d;
        putU8(out, (uint8_t)'f');
        uint8_t raw[4];
        memcpy(raw, &f, 4);
        out.insert(out.end(), raw, raw + 4);
        break;
    }
    case 'Z':
    case 'H': {
        putU8(out, (uint8_t)type);
        out.insert(out.end(), value, value + vlen);
        putU8(out, 0);
        break;
    }
    case 'B': {
        /* VALUE has the form "SUBTYPE,V1,V2,..." */
        if (vlen < 3 || value[1] != ',') {
            return -1;
        }
        const char sub = value[0];
        std::vector<double> vals;
        {
            size_t i = 2;
            while (i <= vlen) {
                size_t j = i;
                while (j < vlen && value[j] != ',') {
                    ++j;
                }
                try {
                    vals.push_back(std::stod(std::string(value + i, j - i)));
                } catch (...) {
                    return -1;
                }
                i = j + 1;
            }
        }
        putU8(out, (uint8_t)'B');
        putU8(out, (uint8_t)sub);
        const size_t countPos = out.size();
        putI32(out, 0);   // Count placeholder, backfilled below
        for (double dv : vals) {
            switch (sub) {
            case 'c': putU8(out, (uint8_t)(int8_t)dv); break;
            case 'C': putU8(out, (uint8_t)dv); break;
            case 's': putU16(out, (uint16_t)(int16_t)dv); break;
            case 'S': putU16(out, (uint16_t)dv); break;
            case 'i': putI32(out, (int32_t)dv); break;
            case 'I': putU32(out, (uint32_t)dv); break;
            case 'f': {
                float f = (float)dv;
                uint8_t raw[4];
                memcpy(raw, &f, 4);
                out.insert(out.end(), raw, raw + 4);
                break;
            }
            default:
                return -1;
            }
        }
        const int32_t cnt = (int32_t)vals.size();
        out[countPos] = (uint8_t)(cnt & 0xFF);
        out[countPos + 1] = (uint8_t)((cnt >> 8) & 0xFF);
        out[countPos + 2] = (uint8_t)((cnt >> 16) & 0xFF);
        out[countPos + 3] = (uint8_t)((cnt >> 24) & 0xFF);
        break;
    }
    default:
        return -1;
    }
    return 0;
}

/* Per-record temporaries reused across a whole block. */
struct BamRecordScratch {
    std::vector<const char*> fields;
    std::vector<size_t> flens;
    std::vector<uint32_t> cigarOps;
    std::vector<uint8_t> packedSeq;
    std::vector<uint8_t> quals;
    std::vector<uint8_t> rec;   /* the built record, block_size prefix included */

    void clear() {
        fields.clear();
        flens.clear();
        cigarOps.clear();
        packedSeq.clear();
        quals.clear();
        rec.clear();
    }
};

/*
 * Build one BAM record from a SAM alignment line. refIndex maps RNAME/RNEXT
 * names to reference ids (built from the @SQ header list).
 *
 * Returns 1 when a record was produced in scratch.rec; 0 when the line was
 * skipped (empty, too few fields, empty QNAME - a warning is logged); -1 on an
 * internal error. The record length in scratch.rec is variable, so callers
 * must read scratch.rec rather than assume a size.
 */
inline int buildBamRecordFromSamLine(const uint8_t* line, size_t len,
                                     const std::map<std::string, int32_t>& refIndex,
                                     BamRecordScratch& scratch) {
    /* Strip trailing newline/carriage return */
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        --len;
    }
    if (len == 0) {
        return 0;
    }

    /* Split fields on \t */
    std::vector<const char*>& fields = scratch.fields;
    std::vector<size_t>& flens = scratch.flens;
    fields.clear();
    flens.clear();
    {
        size_t start = 0;
        for (size_t i = 0; i <= len; ++i) {
            if (i == len || line[i] == '\t') {
                fields.push_back((const char*)line + start);
                flens.push_back(i - start);
                start = i + 1;
            }
        }
    }
    if (fields.size() < 11) {
        LOG_WARNING("BamRecord: SAM record with too few fields, skipped.");
        return 0;   // Too few fields; skip the invalid record
    }

    const char* qname = fields[0];
    const size_t qlen = flens[0];
    const char* rname = fields[2];
    const size_t rlen = flens[2];
    const char* rnext = fields[6];
    const size_t rnextLen = flens[6];
    const char* cigarStr = fields[5];
    const size_t clen = flens[5];
    const char* seqStr = fields[9];
    const size_t slen = flens[9];
    const char* qualStr = fields[10];
    const size_t qlen2 = flens[10];

    if (qlen == 0) {
        return 0;   // QNAME must not be empty; skip the invalid record
    }

    const uint16_t flag = (uint16_t)strtoul(fields[1], nullptr, 10);
    const uint8_t mapq = (uint8_t)strtoul(fields[4], nullptr, 10);
    const int64_t pos1 = strtoll(fields[3], nullptr, 10);      // 1-based
    const int64_t pnext1 = strtoll(fields[7], nullptr, 10);    // 1-based
    const int32_t btlen = (int32_t)strtol(fields[8], nullptr, 10);

    /* refID: RNAME '*' / empty -> -1, otherwise look up the reference sequence list */
    int32_t refId = -1;
    if (rlen > 0 && !(rlen == 1 && rname[0] == '*')) {
        const std::string nm(rname, rlen);
        std::map<std::string, int32_t>::const_iterator it = refIndex.find(nm);
        if (it != refIndex.end()) {
            refId = it->second;
        } else {
            LOG_WARNING("BamRecord: RNAME '%s' not in @SQ, treat as unmapped.", nm.c_str());
        }
    }
    int32_t bpos = (refId >= 0 && pos1 > 0) ? (int32_t)(pos1 - 1) : -1;

    int32_t nextRefId = -1;
    if (rnextLen == 1 && rnext[0] == '=') {
        nextRefId = refId;   // '=' means the same reference as RNAME
    } else if (rnextLen > 0 && !(rnextLen == 1 && rnext[0] == '*')) {
        const std::string nm(rnext, rnextLen);
        std::map<std::string, int32_t>::const_iterator it = refIndex.find(nm);
        if (it != refIndex.end()) {
            nextRefId = it->second;
        }
    }
    int32_t bnext = (nextRefId >= 0 && pnext1 > 0) ? (int32_t)(pnext1 - 1) : -1;

    /* CIGAR */
    std::vector<uint32_t>& cigarOps = scratch.cigarOps;
    int64_t refSpan = 0;
    if (!(clen == 1 && cigarStr[0] == '*')) {
        if (!parseCigar(cigarStr, clen, cigarOps, refSpan)) {
            LOG_WARNING("BamRecord: invalid CIGAR, treat as no CIGAR.");
            cigarOps.clear();
            refSpan = 0;
        }
    } else {
        cigarOps.clear();
        refSpan = 0;
    }
    const uint16_t nCigar = (uint16_t)cigarOps.size();

    /* bin */
    uint16_t bin = 0;
    if (refId >= 0 && bpos >= 0) {
        const int64_t end = (int64_t)bpos + (refSpan > 0 ? refSpan : 1);
        bin = samReg2Bin(bpos, end);
    }

    /* SEQ */
    int32_t lSeq = 0;
    std::vector<uint8_t>& packedSeq = scratch.packedSeq;
    packedSeq.clear();
    if (slen > 0 && !(slen == 1 && seqStr[0] == '*')) {
        lSeq = (int32_t)slen;
        packedSeq.resize((size_t)((slen + 1) / 2), 0);
        for (size_t i = 0; i < slen; ++i) {
            const int nib = baseNibble(seqStr[i]);
            if (i & 1) {
                packedSeq[i >> 1] |= (uint8_t)(nib & 0xF);
            } else {
                packedSeq[i >> 1] = (uint8_t)((nib & 0xF) << 4);
            }
        }
    }

    /* QUAL: all 0xFF when missing (*); otherwise per-byte ascii-33 */
    std::vector<uint8_t>& quals = scratch.quals;
    quals.clear();
    const bool qualMissing = (qlen2 == 0) || (qlen2 == 1 && qualStr[0] == '*');
    if (!qualMissing) {
        size_t qn = (qlen2 < (size_t)lSeq) ? qlen2 : (size_t)lSeq;
        quals.resize(qn);
        for (size_t i = 0; i < qn; ++i) {
            int phred = (int)(uint8_t)qualStr[i] - 33;
            if (phred < 0) {
                phred = 0;
            }
            if (phred > 93) {
                phred = 93;
            }
            quals[i] = (uint8_t)phred;
        }
    }

    /* read_name (BAM requires l_read_name <= 255, including the trailing \0) */
    size_t qn = (qlen < 254) ? qlen : 254;
    const uint8_t lReadName = (uint8_t)(qn + 1);

    /* Assemble the record */
    std::vector<uint8_t>& rec = scratch.rec;
    rec.clear();
    rec.reserve(32 + qn + cigarOps.size() * 4 + packedSeq.size() + quals.size() + 16);
    putI32(rec, 0);   // block_size placeholder
    putI32(rec, refId);
    putI32(rec, bpos);
    putU8(rec, lReadName);
    putU8(rec, mapq);
    putU16(rec, bin);
    putU16(rec, nCigar);
    putU16(rec, flag);
    putI32(rec, lSeq);
    putI32(rec, nextRefId);
    putI32(rec, bnext);
    putI32(rec, btlen);
    rec.insert(rec.end(), qname, qname + qn);
    rec.push_back(0);
    for (size_t i = 0; i < cigarOps.size(); ++i) {
        putU32(rec, cigarOps[i]);
    }
    rec.insert(rec.end(), packedSeq.begin(), packedSeq.end());
    /* QUAL: append per byte when present; pad with 0xFF when missing (*) or shorter than SEQ, keeping the record length consistent with the BAM layout */
    if (!quals.empty()) {
        rec.insert(rec.end(), quals.begin(), quals.end());
    }
    for (int32_t i = (int32_t)quals.size(); i < lSeq; ++i) {
        rec.push_back(0xFF);
    }
    for (size_t i = 11; i < fields.size(); ++i) {
        if (0 != appendSamAuxToBam(fields[i], flens[i], rec)) {
            LOG_WARNING("BamRecord: skip invalid SAM option field.");
        }
    }

    const int32_t blockSize = (int32_t)(rec.size() - 4);
    memcpy(rec.data(), &blockSize, 4);
    return 1;
}

}  // namespace bamrec
