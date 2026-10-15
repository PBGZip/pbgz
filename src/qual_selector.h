/*
 * qual_selector.h - encoder evaluation for the quality-value column. The column needs its own
 * path: its candidates (coder_qual, fcv2) need record-level input, and coder_qual.h cannot share
 * a translation unit with coder_fc.h.
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
#include <vector>
#include <string>

#include "preprocess_info.h"

/* One record from the sample; both evaluated candidates need this information. */
struct QualSampleRecord {
    std::string qual;   /* quality value */
    std::string seq;    /* corresponding base sequence, used by coder_qual as context */
    bool        rev;    /* strand direction, taken from bit 0x10 of FLAG, used by fcv2 to recover the cycle index */
};

/* Mean record length at which a file counts as long-read; fcv2 is dropped at or above it. */
extern const uint32_t QUAL_LONG_READ_LEN;

class QualSelector {
public:
    /* Trial-compress the candidates on the sampled records and return the smallest. freqByByte
       counts quality-value occurrences by raw byte value. Returns SKIPPED for a too-small sample
       and FAILED when no candidate runs; both leave the caller's default encoder in place. */
    static FieldCodecSelection select(const std::vector<QualSampleRecord>& records,
                                      const std::vector<uint32_t>& freqByByte,
                                      uint8_t compressLevel = 5);

    /* Same, with the candidate set given by the caller (how the FASTQ QUAL column is evaluated). */
    static FieldCodecSelection select(const std::vector<QualSampleRecord>& records,
                                      const std::vector<uint32_t>& freqByByte,
                                      uint8_t compressLevel,
                                      const std::vector<CoderType>& candidates);
};
