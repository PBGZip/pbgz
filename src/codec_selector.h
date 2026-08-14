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
 * 一个按行切开的字段样本。coder_affix_match 是逐行编码器，前后缀匹配的收益来自
 * 相邻行之间，试压时必须按行喂（含行尾 tab），整列连成一个字节流会丢掉行边界，
 * 测出来的数字和真实压缩完全不是一回事。
 */
struct LineSample {
    const uint8_t* data;
    uint32_t len;
};

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
    /* 先验训练最多使用 45 MB QUAL。实测到此处收益约为 1.84 个百分点，继续加到 90 MB
       只多约 0.10 个百分点；多出的时间却发生在切换到并行压缩之前的串行阶段，会直接
       推迟后续工作线程启动，所以必须在收益曲线变平前截断。 */
    static constexpr uint32_t QUAL_PRIOR_TRAIN_MAX = 45u * 1024u * 1024u; // TEMP

    /*
     * Analyze a sample of the given block and fill info.
     * Returns 0 on success (including the "some fields skipped" case),
     * -1 only on a fatal error such as a null block.
     *
     * inputTotalBytes 是整个输入的字节数，0 表示不可知（管道输入）。它只用于
     * 判断按全文件数据量摊销才划算的决策——目前是 QUAL 先验要不要训练和写出。
     * 做成入参而不是塞进 PreprocessInfo：后者是决策的产物，不该同时兼作输入。
     */
    static int32_t analyze(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info);

    /*
     * Trial-compress one byte stream with every candidate coder and return the
     * selection result. Exposed for unit testing.
     *
     * trialAffix 控制是否把 coder_affix_match 纳入候选。affix 是列式前后缀匹配编码器，
     * 只对部分 SAM 常规字段（FLAG/POS/MAPQ/CIGAR/PNEXT/TLEN）有收益，其他字段（如
     * SEQ、QNAME）参与比较只会浪费时间，甚至可能被误选。
     *
     * affix 必须按行喂数据（见 LineSample），lines 给出逐行样本；为空时退化为
     * 整段一次 encode_line，测出来的结果不可信，此时 affix 不参与比较。
     */
    static FieldCodecSelection selectCoder(const uint8_t* data, uint32_t len,
                                           bool trialAffix = false,
                                           const std::vector<LineSample>* lines = nullptr);

    /*
     * QUAL 先验的跨块训练样本累积器。由 CompressEngine 在读线程持有：前 N 个块的
     * QUAL 按记录累积进来（每块扫描到块尾，但收集到 QUAL_PRIOR_TRAIN_MAX 即止），
     * 读完目标块数或累积满后由 trainQualPriorModel 一次性训练并导出快照。
     *
     * 为什么跨块：首块通常只有 ~9 MB QUAL，而先验收益在 45 MB 处才基本收敛（见
     * QUAL_PRIOR_TRAIN_MAX）。多块累积让小块也能逼近这个训练量。
     */
    struct QualPriorAccum {
        std::vector<QualSampleRecord> records;
        std::vector<uint32_t> freqByByte;
        uint64_t collectedBytes = 0;

        bool full() const { return collectedBytes >= QUAL_PRIOR_TRAIN_MAX; }
    };

    /* 把一块的 QUAL 记录追加进累积器（上限 QUAL_PRIOR_TRAIN_MAX，超出部分丢弃）。 */
    static void accumulateQualPrior(RoughIOBlock* block, QualPriorAccum& acc);

    /* 在累积样本上训练 fcv2 并导出模型快照；trainedBytes 输出实际训练字节数。
     * params 是 QualSelector 选定的上下文参数档位（仅 fcv2 胜出时有意义），必须与
     * 实际压缩用的档位一致，否则先验快照加载会因布局不符而整体回退。
     * 失败或样本为空返回空向量。 */
    static std::vector<uint8_t> trainQualPriorModel(const QualPriorAccum& acc,
                                                    const QualFcv2Params& params,
                                                    uint64_t* trainedBytes);

private:
    static int32_t analyzeSam(RoughIOBlock* block, uint64_t inputTotalBytes, PreprocessInfo& info);
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
                                   uint32_t sampleBudget,
                                   uint64_t qualBudget = UINT64_MAX);

    static uint32_t extractSamFieldSamples(RoughIOBlock* block,
                                       std::vector<std::string>& fieldBufs,
                                       std::vector<std::vector<LineSample>>& fieldLines,
                                       uint32_t sampleBudget);
    static uint32_t extractFastqFieldSamples(RoughIOBlock* block,
                                         std::vector<std::string>& fieldBufs,
                                         uint32_t sampleBudget);

    /* Choose the smallest coder_bwt_cm block level whose size fits the sample. */
    static int pickBwtLevel(uint32_t sampleLen);
};
