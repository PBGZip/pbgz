/*
 * sam_field_layout.h - the byte layouts a SAM field payload can be in, and how a stream says which.
 *
 * The archive records a layout per sub-stream in its meta, and a decoder has to read every layout a
 * released revision wrote, not only the one this build writes. Those readers live here - pure
 * functions over a payload buffer - instead of inline in the record loop, so that the compatibility
 * reading is in one place and each layout can be pinned by a test without an archive (see
 * test/testcase/sam_field_layout_testcase.cpp).
 *
 * What this build writes:
 *   - TLEN exceptions: (delta line index, zigzag32 TLEN) varint pairs, marked "exc" = 1;
 *   - SEQ exception positions: forward deltas as varints ("nposd"), or runs ("nposr");
 *   - SEQ exception positions, absolute: "npos" together with a "ch".
 * What older revisions wrote and is still read:
 *   - TLEN exceptions as fixed int32 pairs, with no marker at all (up to 2026-08-19);
 *   - SEQ exception positions as one absolute list under "npos" with no "ch", holding 'N' and 'n'
 *     alike, counted by the field-level "ncount", and writing 'N' back at every position.
 */

#pragma once

#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include <json/json.h>

/*
 * LEB128 varint and zigzag helpers, used by the position, length and numeric sub-streams (the
 * names date from the TLEN stream, which was the first to replace a fixed-width layout with them).
 * Both replace layouts whose high bytes were almost always 0, and zigzag keeps small negative
 * differences small by moving the sign into the low bit.
 */
static inline uint32_t tlenPutVarint(uint8_t* p, uint32_t v) {
    uint32_t n = 0;
    while (v >= 0x80) {
        p[n++] = (uint8_t)(v | 0x80);
        v >>= 7;
    }
    p[n++] = (uint8_t)v;
    return n;
}

static inline uint32_t tlenZigzag32(int32_t v) {
    return (uint32_t)((v << 1) ^ (v >> 31));
}

static inline int32_t tlenUnzigzag32(uint32_t v) {
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

/* Zigzag: the sign of a signed delta in the low bit, so small negatives stay small. */
static inline uint64_t zigzag64(int64_t value)
{
    return ((uint64_t)value << 1) ^ (uint64_t)(value >> 63);
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

/* The marker the TLEN field meta carries, naming the layout its exception stream is in. */
enum { TLEN_EXC_LAYOUT_FIXED = 0, TLEN_EXC_LAYOUT_VARINT = 1 };

/*
 * The two layouts a TLEN exception stream arrives in. This build writes (delta line index, zigzag32
 * TLEN) as varint pairs and records that in the field meta under the short key "exc" - the marker
 * sits in the meta of every block that has exceptions, so its name is three characters and its
 * value a small integer rather than a string. The layout before it wrote (line, TLEN) as fixed
 * int32 pairs - eight bytes per entry, the line an absolute content index - and wrote no marker at
 * all, so an archive without "exc" holds that one.
 *
 * The marker is authoritative where it is present: the layout it names is the one decoded, and a
 * stream that does not decode under it is corruption rather than the other layout. Where it is
 * absent the fixed layout is assumed, as an archive predating the marker must - unless the stream's
 * size rules that out, since eight bytes per entry and the count the field declares are the one
 * thing the fixed layout must satisfy.
 */
static inline bool tlenDecodeVarints(const uint8_t* buf, uint32_t srclen, uint32_t declared,
                                     std::map<uint32_t, int32_t>& out)
{
    uint32_t off = 0, line = 0, entries = 0;
    while (off < srclen) {
        uint32_t delta = 0, zz = 0;
        if (!readVarint(buf, srclen, off, delta) || !readVarint(buf, srclen, off, zz)) {
            return false;
        }
        line += delta;
        out[line] = tlenUnzigzag32(zz);
        entries++;
    }
    return declared == 0 || entries == declared;
}

static inline bool tlenDecodeFixedPairs(const uint8_t* buf, uint32_t srclen, uint32_t declared,
                                        std::map<uint32_t, int32_t>& out)
{
    if (declared == 0 || srclen != declared * 2 * (uint32_t)sizeof(uint32_t)) {
        return false;
    }
    for (uint32_t i = 0; i < declared; ++i) {
        uint32_t line = 0, value = 0;
        memcpy(&line, buf + (size_t)i * 2 * sizeof(uint32_t), sizeof(uint32_t));
        memcpy(&value, buf + ((size_t)i * 2 + 1) * sizeof(uint32_t), sizeof(uint32_t));
        out[line] = (int32_t)value;
    }
    return true;
}

/*
 * The SEQ exception sub-streams: one per character present, each a strictly increasing list of
 * block offsets, and each named for the form its positions are written in. The names are the single
 * source of truth for the encoder and the reader - they are the meta, so both sides have to agree
 * on them exactly.
 */
static const char kSeqExcDeltaName[] = "nposd";
static const char kSeqExcAbsName[] = "npos";
static const char kSeqExcRunName[] = "nposr";

enum class SeqExcForm { None, Absolute, Delta, Run };

/* None when the name is not one of the forms; the decoder then knows the block's exception streams
   have ended, since they are written consecutively. */
static inline SeqExcForm seqExcFormOf(const std::string& sname)
{
    if (sname == kSeqExcDeltaName) {
        return SeqExcForm::Delta;
    }
    if (sname == kSeqExcAbsName) {
        return SeqExcForm::Absolute;
    }
    if (sname == kSeqExcRunName) {
        return SeqExcForm::Run;
    }
    return SeqExcForm::None;
}

/* What a SEQ exception sub-stream's meta says: the form, the character to write at every position it
   carries, how many positions, and for the run form how many runs. */
struct SeqExcLayout {
    SeqExcForm form;
    uint8_t ch;
    uint32_t count;
    uint32_t runs;
};

static inline SeqExcLayout seqExcLayoutOf(const Json::Value& fieldMeta, const Json::Value& streamMeta)
{
    SeqExcLayout out;
    out.form = seqExcFormOf(streamMeta["sname"].asString());
    /* The first layout's "npos" predates "ch": it held 'N' and 'n' alike and wrote 'N' back at
       every position, which is what a reader without the key has to do. */
    const Json::Value& chJson = streamMeta["ch"];
    out.ch = chJson.isUInt() ? (uint8_t)chJson.asUInt() : (uint8_t)'N';
    /* A stream's own count; "ncount" is the legacy form of it for that "npos". */
    const Json::Value& countJson = streamMeta["count"];
    out.count = countJson.isUInt() ? countJson.asUInt() : fieldMeta["ncount"].asUInt();
    out.runs = streamMeta["runs"].asUInt();
    return out;
}

/*
 * Expands a decoded payload into the block offsets it names. Absolute holds one 4-byte offset per
 * position; delta holds one varint per position, each the step from the previous one; run holds
 * both of its streams as varints, all the gaps first and then all the lengths (see the encoder), and
 * the runs it describes must add up to the declared count.
 */
static inline bool seqExcDecodePositions(const SeqExcLayout& layout, const uint8_t* buf,
                                         uint32_t srclen, std::vector<uint32_t>& out)
{
    out.clear();
    switch (layout.form) {
        case SeqExcForm::Absolute: {
            if (srclen != layout.count * (uint32_t)sizeof(uint32_t)) {
                return false;
            }
            out.resize(layout.count);
            if (layout.count > 0) {
                memcpy(out.data(), buf, srclen);
            }
            return true;
        }
        case SeqExcForm::Delta: {
            out.resize(layout.count);
            uint32_t off = 0, pos = 0, step = 0;
            for (uint32_t i = 0; i < layout.count; ++i) {
                if (!readVarint(buf, srclen, off, step)) {
                    return false;
                }
                pos += step;
                out[i] = pos;
            }
            return true;
        }
        case SeqExcForm::Run: {
            std::vector<uint32_t> gaps(layout.runs, 0);
            std::vector<uint32_t> lens(layout.runs, 0);
            uint32_t off = 0, value = 0;
            for (uint32_t i = 0; i < layout.runs; ++i) {
                if (!readVarint(buf, srclen, off, value)) {
                    return false;
                }
                gaps[i] = value;
            }
            for (uint32_t i = 0; i < layout.runs; ++i) {
                if (!readVarint(buf, srclen, off, value)) {
                    return false;
                }
                lens[i] = value;
            }
            out.reserve(layout.count);
            uint32_t end = 0;
            for (uint32_t i = 0; i < layout.runs; ++i) {
                end += gaps[i];
                for (uint32_t k = 0; k < lens[i]; ++k) {
                    out.push_back(end);
                    end++;
                }
            }
            return out.size() == layout.count;
        }
        case SeqExcForm::None:
        default:
            return false;
    }
}
