/*
 * coder_fcv2.h - quality-value context-mixing coder for the SAM QUAL column.
 *
 * Builds several probability models of differing granularity per quality
 * symbol, mixes their predictions in the log-odds domain and feeds a binary
 * range coder. The interface is per record (strand and read length are
 * record-level). The implementation lives in the .cpp so this header can avoid
 * including fc/rangecoder.h (name clash with clr.h's RangeCoder).
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
#include <memory>
#include <vector>

#include "coder_io.h"

class fcv2_impl;

/*
 * Mixing models used by fcv2: m0..m4, always all five. The count is fixed and
 * is not carried in the stream header.
 */
const int FCV2_MODEL_COUNT = 5;

/* Huffman tree cap of the value alphabet; the coder's internal TREE_CAP. */
const int FCV2_TREE_CAP = 64;

/*
 * fcv2 stream header version, written as the stream's first symbol. Bump when
 * the parameter set, their order or the model layout changes.
 */
const uint8_t FCV2_STREAM_VERSION = 3;

/*
 * Parameter tiers of the fcv2 context model; chosen per data by QualSelector
 * and read back from the stream header.
 *   cycleMax / cycleBucket  position context upper bound and bin count
 *   prevShift               right shift quantizing the predecessor value in m3
 *   useDedup                enable adjacent duplicate read dedup
 */
struct Fcv2Cfg {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  prevShift   = 1;
    bool useDedup    = false;

    bool operator==(const Fcv2Cfg& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               prevShift == o.prevShift && useDedup == o.useDedup;
    }
    bool operator!=(const Fcv2Cfg& o) const { return !(*this == o); }
};

/*
 * The QUAL context presets worth trial-compressing against each other; the
 * selector trials them and keeps the smallest. minVolumeTier is the
 * levelVolumeTier() from which trialling the preset is worth the time. The
 * model count is fixed at FCV2_MODEL_COUNT for every preset.
 */
struct Fcv2Preset {
    const char* name;
    int         minVolumeTier;
    Fcv2Cfg     cfg;
};

const int FCV2_PRESET_COUNT = 3;

/* The three measured presets: ultra, ultra-80 and fine. */
const Fcv2Preset FCV2_PRESETS[FCV2_PRESET_COUNT] = {
    { "ultra",    0, {   96,  4, 0, false } },
    { "ultra-80", 0, {   80,  4, 0, false } },
    { "fine",     3, {   96, 24, 0, false } },
};

/* The name of a configuration, or nullptr when it is not one of the presets. */
inline const char* fcv2PresetName(const Fcv2Cfg& cfg)
{
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        const Fcv2Cfg& p = FCV2_PRESETS[i].cfg;
        if (cfg.cycleMax == p.cycleMax && cfg.cycleBucket == p.cycleBucket &&
            cfg.prevShift == p.prevShift && cfg.useDedup == p.useDedup) {
            return FCV2_PRESETS[i].name;
        }
    }
    return nullptr;
}

/*
 * The alphabet is data-derived: only the quality values that occur in the block
 * are admitted, in ascending byte order. A byte outside it is skipped by
 * encode_record, not rejected. The tree cap (FCV2_TREE_CAP) is a real ceiling:
 * a block exceeding it leaves alphaSize at 0 and the coder unusable.
 */
class coder_fcv2 {
public:
    /*
     * freqTable passes the occurrence counts of the quality values, indexed by
     * symbol value (raw byte minus '!'), and is used to build the Huffman tree.
     * It improves speed but not the compression ratio.
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable);

    /*
     * Whether the stream's alphabet contains this quality byte; a value outside
     * it cannot be coded (encoding one fails with coder_io::IO_UNCODABLE).
     */
    bool coversByte(uint8_t b) const;

    /* Version taking context parameter tiers; cfg is written to the stream header. */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg);

    /*
     * Creates a coder from a previously exported model snapshot. modelLoaded
     * reports whether the snapshot validated; on failure the model built from
     * freqTable stays in place.
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);

    /* Version taking context parameter tiers + a prior snapshot. */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);
    ~coder_fcv2();

    /*
     * Exports the currently learned model. Returns false if there is no valid
     * alphabet to export or the output buffer could not be allocated.
     */
    bool export_model(std::vector<uint8_t>& out) const;

    /*
     * This coder needs the length and strand of every record; only the QUAL
     * column of an aligned SAM can provide these. The applicability check lives
     * in CoderFactory::coderSupports.
     */

    /*
     * Encodes one record of quality values. rev comes from FLAG's 0x10 bit and
     * recovers the true sequencing-cycle index. seq is the base sequence used
     * as the condition for the base context model; nullptr falls into the
     * "unknown" bin.
     */
    void encode_record(const uint8_t* qual, uint32_t len, bool rev,
                       const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* End of encoding; returns the number of bytes written to io. */
    int32_t encode_flush();

    /*
     * Decodes one record of quality values into dst. The strand is read from
     * the stream; seq is the base sequence used as the condition for the base
     * context model and must agree with the encoding side.
     */
    int32_t decode_record(uint8_t* dst, uint32_t len,
                          const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* Called once before decoding; reads the stream header. */
    int32_t begin_decode();

private:
    std::unique_ptr<fcv2_impl> impl;

    coder_fcv2(const coder_fcv2&) = delete;
    coder_fcv2& operator=(const coder_fcv2&) = delete;
};
