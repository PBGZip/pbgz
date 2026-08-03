/*
 * preprocess_info.h - Data model for the file preprocessing module
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
#include <atomic>
#include <string>
#include <utility>
#include <vector>

#include "io_block.h"

/*
 * Candidate coders that the pre-selector may trial-compress with.
 * The magic strings match what each coder writes into coder_io::meta,
 * so the decompression side can dispatch on them unchanged.
 */
enum class CoderType : uint8_t {
    BWT_CM = 0,   /* coder_bwt_cm: BWT + context model (current default)   */
    FC,           /* coder_fc:     LZP + BWT + MTF + context range coder   */
    SIMPLE_RC,    /* coder_simple_rc: 有损，已从编解码路径移除，仅保留枚举占位 */
    FCV2,         /* coder_fcv2:   质量值上下文混合编码器，仅适用于 SAM 的 QUAL 列 */
    QUAL,         /* coder_qual:   quality-specific context model coder    */
    COUNT
};

/* Map a CoderType to the magic string stored in block meta. */
static inline const char* coderTypeToMagic(CoderType type)
{
    switch (type) {
    case CoderType::BWT_CM:    return "coder_bwt_cm";
    case CoderType::FC:        return "coder_fc";
    case CoderType::SIMPLE_RC: return "coder_simple_rc";
    case CoderType::FCV2:      return "coder_fcv2";
    case CoderType::QUAL:      return "coder_qual";
    default:                   return "unknown";
    }
}

/* Selection status for a single field. */
enum class FieldStatus : uint8_t {
    SELECTED = 0,  /* a coder was chosen successfully                  */
    FAILED,        /* trial compression failed; use the default coder  */
    SKIPPED        /* field not analyzed (empty / too small / n-a)     */
};

/* Codec selection result for one field. */
struct FieldCodecSelection {
    FieldStatus status;
    CoderType   selectedCoder;
    uint32_t    sampleLen;
    uint32_t    bestCompLen;

    uint32_t    trialBwtCmLen = 0;
    uint32_t    trialFcLen = 0;

    /*
     * 试压耗时，微秒。
     *
     * 只看压缩率来选编码器是不够的：两个编码器压缩率相差不到一个百分点、吞吐却差
     * 好几倍的情况很常见，这时候慢的那个未必划算。把耗时一并记下来，选择策略以后
     * 就能在压缩率和速度之间做取舍，而不是无条件取最小。
     *
     * 注意样本只有几百 KB 到一两 MB，单次计时的相对误差不小（同一二进制重复测量的
     * 波动在百分之几的量级），所以这个数字适合用来看数量级差异，不适合拿来比较
     * 两个相差百分之几的编码器。
     */
    uint32_t    trialBwtCmUs = 0;
    uint32_t    trialFcUs = 0;

    /*
     * 实际用于定案的样本字节数，以及为此跑了几轮。
     *
     * 评估不再是"把整个采样压一遍"，而是从小样本起步逐轮加倍，一旦领先者拉开足够
     * 差距就立刻收手。多数字段在头一两轮就能分出胜负，剩下的预算留给真正难分的字段。
     * decidedLen 小于 sampleLen 就说明这个字段提前定案了。
     */
    uint32_t    decidedLen = 0;
    uint32_t    rounds = 0;

    FieldCodecSelection()
        : status(FieldStatus::SKIPPED),
          selectedCoder(CoderType::BWT_CM),
          sampleLen(0),
          bestCompLen(0) {}

    /* Compression ratio of the selected coder on the sample (1.0 = no gain). */
    double ratio() const
    {
        return (sampleLen == 0) ? 1.0 : (double)bestCompLen / (double)sampleLen;
    }
};

/* SAM mandatory field indices (0-based, matching the SAM spec columns). */
enum SamField : uint32_t {
    SAM_QNAME = 0,
    SAM_FLAG,
    SAM_RNAME,
    SAM_POS,
    SAM_MAPQ,
    SAM_CIGAR,
    SAM_RNEXT,
    SAM_PNEXT,
    SAM_TLEN,
    SAM_SEQ,
    SAM_QUAL,
    SAM_FIELD_COUNT   /* 11 mandatory fields */
};

/* FASTQ field indices (one record = ID / SEQ / comment / QUAL lines). */
enum FastqField : uint32_t {
    FQ_ID = 0,
    FQ_SEQ,
    FQ_QUAL,
    FQ_COMMENT,
    FQ_FIELD_COUNT    /* 4 fields */
};

/*
 * Result of file preprocessing.
 *
 * Preprocessing runs exactly once, independent of the compression pipeline's
 * block-level concurrency.  The first worker to arrive for a SAM/FASTQ block
 * performs the trial-compression; other workers see "in progress" and simply
 * use their default coder without blocking, so the pipeline never stalls.
 *
 * The publish pattern is lock-free: writers set `state` to DONE with
 * memory_order_release; readers acquire-load `state` and only consult `fields`
 * after seeing DONE, which gives them a consistent snapshot without locks.
 *
 * The struct is intentionally extensible: future preprocessing results (ID
 * split-symbol statistics, quality-score range, read-length histogram, ...)
 * can be added as new members without changing the selection contract.
 */
enum class PreprocessState : uint8_t {
    IDLE = 0,
    RUNNING,
    DONE
};

struct PreprocessInfo {
    BlockType fileType;
    std::atomic<PreprocessState> state;

    /* 各字段样本长度之和，不含 tab 分隔符、换行符和头部行。 */
    uint32_t  sampleBytes;

    /* 为做分析实际扫过的原始数据字节数，也就是从数据块头部截下来的那一段的大小。 */
    uint32_t  scannedBytes;

    std::vector<FieldCodecSelection> fields;

    /*
     * fcv2 由首块完整 QUAL 训练出的模型快照及其训练量。
     *
     * 先验跟随预处理结果保存，而不是交给某个编码器实例：PreprocessInfo 是流水线和
     * 编码器决策之间唯一的交接点，编码器层因此无需知道数据块或文件偏移，继续保持
     * coder/ 不反向依赖 src/ 的分层。空快照明确表示本次没有可用先验，调用方可无条件
     * 退回既有的冷启动行为。
     */
    std::vector<uint8_t> qualPriorSnapshot;
    uint64_t qualPriorTrainingBytes;

    PreprocessInfo() : fileType(TYPE_UNKNOW), state(PreprocessState::IDLE), sampleBytes(0), scannedBytes(0),
                       qualPriorTrainingBytes(0) {}

    void reset(BlockType type)
    {
        fileType = type;
        state.store(PreprocessState::IDLE, std::memory_order_relaxed);
        sampleBytes = 0;
        scannedBytes = 0;
        fields.clear();
        qualPriorSnapshot.clear();
        qualPriorTrainingBytes = 0;
    }

    bool isDone() const
    {
        return state.load(std::memory_order_acquire) == PreprocessState::DONE;
    }

    /*
     * Atomically claim the right to run preprocessing.
     * Returns true if this caller won the claim (should run preprocessing),
     * false if another thread is already running or has finished.
     */
    bool tryClaim()
    {
        PreprocessState expected = PreprocessState::IDLE;
        return state.compare_exchange_strong(
            expected, PreprocessState::RUNNING,
            std::memory_order_acq_rel, std::memory_order_acquire);
    }

    void markDone()
    {
        state.store(PreprocessState::DONE, std::memory_order_release);
    }

    const FieldCodecSelection* getField(uint32_t idx) const
    {
        return (idx < fields.size()) ? &fields[idx] : nullptr;
    }

    /*
     * Return the coder to use for a field.  If preprocessing is DONE and
     * selected a coder, return it; otherwise return fallback (the actuator's
     * default).  Safe to call from any worker without blocking.
     */
    CoderType coderFor(uint32_t idx, CoderType fallback) const
    {
        if (!isDone()) {
            return fallback;
        }
        const FieldCodecSelection* sel = getField(idx);
        if (sel != nullptr && sel->status == FieldStatus::SELECTED) {
            return sel->selectedCoder;
        }
        return fallback;
    }

    const std::vector<uint8_t>& qualPrior() const
    {
        return qualPriorSnapshot;
    }

    uint64_t qualPriorTrainedBytes() const
    {
        return qualPriorTrainingBytes;
    }

    void setQualPrior(std::vector<uint8_t>&& snapshot, uint64_t trainedBytes)
    {
        qualPriorSnapshot = std::move(snapshot);
        qualPriorTrainingBytes = qualPriorSnapshot.empty() ? 0 : trainedBytes;
    }
};
