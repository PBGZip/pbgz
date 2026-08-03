/*
 * codec_selector.h - Codec pre-selection for the file preprocessing module
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
#include <string>
#include <vector>

#include "preprocess_info.h"
#include "qual_selector.h"

class RoughIOBlock;

/*
 * CodecSelector
 *
 * Codec pre-selection step of the file preprocessing module. It scans a sample
 * of the first data block, splits it into per-field byte streams (SAM columns
 * or FASTQ lines), trial-compresses each stream with every candidate coder and
 * records the best one into a PreprocessInfo.
 *
 * The selector is self-contained: it parses the block itself and does not rely
 * on any actuator state, which keeps it independently unit-testable.
 *
 * Failure policy: any field that is too small to give a reliable comparison, or
 * whose trial compression fails, is marked SKIPPED / FAILED. Actuators then fall
 * back to their built-in default coder, so a preprocessing failure never breaks
 * compression.
 */
class CodecSelector {
public:
    /*
     * Analyze a sample of the given block and fill info.
     * Returns 0 on success (including the "some fields skipped" case),
     * -1 only on a fatal error such as a null block.
     */
    static int32_t analyze(RoughIOBlock* block, PreprocessInfo& info);

    /*
     * Trial-compress one byte stream with every candidate coder and return the
     * selection result. Exposed for unit testing.
     */
    static FieldCodecSelection selectCoder(const uint8_t* data, uint32_t len);

private:
    static int32_t analyzeSam(RoughIOBlock* block, PreprocessInfo& info);
    static int32_t analyzeFastq(RoughIOBlock* block, PreprocessInfo& info);

    /* Extract per-field concatenated samples from a block. */
    /*
     * 单独采集质量值列的样本。
     *
     * 通用采样把整列拼成连续字节，记录边界就没了；而质量值的两个候选都需要记录级
     * 信息——fcv2 要每条记录的长度来还原测序循环序号，coder_qual 要对应的碱基序列
     * 作上下文，链方向则取自 FLAG 的 0x10 位。所以这里按记录收集，不做拼接。
     */
    static void extractQualSamples(RoughIOBlock* block,
                                   std::vector<QualSampleRecord>& records,
                                   std::vector<uint32_t>& freqByByte,
                                   uint32_t sampleBudget);

    static uint32_t extractSamFieldSamples(RoughIOBlock* block,
                                       std::vector<std::string>& fieldBufs,
                                       uint32_t sampleBudget);
    static uint32_t extractFastqFieldSamples(RoughIOBlock* block,
                                         std::vector<std::string>& fieldBufs,
                                         uint32_t sampleBudget);

    /* Choose the smallest coder_bwt_cm block level whose size fits the sample. */
    static int pickBwtLevel(uint32_t sampleLen);
};
