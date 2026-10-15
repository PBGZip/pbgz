/*
 * sam_field_layout.h - the byte layouts a SAM field payload can be in, and how a stream says
 * which. A decoder must read every layout a released revision wrote, so these pure readers over a
 * payload buffer live here rather than inline in the record loop (see
 * test/testcase/sam_field_layout_testcase.cpp).
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

#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include <json/json.h>

/* Varint/zigzag helpers live in seq_stream_util.h. */
#include "seq_stream_util.h"

/* The marker the TLEN field meta carries, naming the layout its exception stream is in. */
enum { TLEN_EXC_LAYOUT_FIXED = 0, TLEN_EXC_LAYOUT_VARINT = 1 };

/*
 * Decode a TLEN exception stream written as (delta line index, zigzag32 TLEN) varint pairs, the
 * layout this build writes and marks "exc" in the field meta. Returns false if the entries do not
 * match the declared count.
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
        out[line] = unzigzag32(zz);
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

/* SEQ exception sub-stream names: each names the form its positions are written in and is the
   meta key both encoder and reader agree on. */
static const char kSeqExcDeltaName[] = "nposd";
static const char kSeqExcAbsName[] = "npos";
static const char kSeqExcRunName[] = "nposr";

enum class SeqExcForm { None, Absolute, Delta, Run };

/* None when the name is not one of the forms (exception streams have ended). */
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

/* What a SEQ exception sub-stream's meta says: form, character, position count, run count. */
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
    /* The legacy "npos" had no "ch"; default to 'N'. */
    const Json::Value& chJson = streamMeta["ch"];
    out.ch = chJson.isUInt() ? (uint8_t)chJson.asUInt() : (uint8_t)'N';
    /* A stream's own count; "ncount" is the legacy form of it for that "npos". */
    const Json::Value& countJson = streamMeta["count"];
    out.count = countJson.isUInt() ? countJson.asUInt() : fieldMeta["ncount"].asUInt();
    out.runs = streamMeta["runs"].asUInt();
    return out;
}

/* Expand a decoded payload into the block offsets it names. Absolute: one 4-byte offset per
   position. Delta: one varint step per position. Run: all gaps as varints then all lengths; the
   runs must add up to the declared count. */
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
