/*
 * field_coder_config.h - SAM 字段编码器配置表
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

#include <algorithm>
#include <cstdint>
#include <vector>

#include "preprocess_info.h"

/*
 * SAM 每字段的编码器配置表：候选编码器（预处理阶段逐字段试压比较）与默认编码器
 * （未选择 / 选择失败 / 未接通选择时执行器的兜底）。修改某个字段的候选或默认编码器
 * 只动这一张表，预处理试压范围与执行器默认值都从这里读，无需再改两处代码。
 *
 * 特例：
 *  - QUAL(10)：走 QualSelector 专用路径（coder_qual / fcv2 / bwt_cm），不参与
 *    通用试压，这里只登记默认值。
 *  - POS(3) / PNEXT(7) / TLEN(8)：固定走差分 / 推算压缩（见 compressPosFieldDelta /
 *    compressPNextFieldDelta / compressTLen），无通用候选，默认值是底层 bwt_cm。
 */

struct FieldCoderConfig {
    std::vector<CoderType> candidates; /* 可选编码器，按序试压比较 */
    CoderType fallback;                /* 默认编码器 */
};

/* 参与编码器选择的 SAM 字段数：11 个必选字段 + 第 12 列 OPTION。 */
static const uint32_t SAM_FIELD_COUNT_SELECT = SAM_FIELD_COUNT + 1;

inline const FieldCoderConfig kSamFieldCoderConfig[SAM_FIELD_COUNT_SELECT] = {
    /* QNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* FLAG  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* RNAME */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* POS   */ {{}, CoderType::BWT_CM},
    /* MAPQ  */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* CIGAR */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
    /* RNEXT */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::BWT_CM},
    /* PNEXT */ {{}, CoderType::BWT_CM},
    /* TLEN  */ {{}, CoderType::BWT_CM},
    /* SEQ   */ {{CoderType::BWT_CM, CoderType::FC}, CoderType::FC},
    /* QUAL  */ {{}, CoderType::QUAL},
    /* OPTION */ {{CoderType::BWT_CM, CoderType::FC, CoderType::AFFIX_MATCH}, CoderType::BWT_CM},
};

inline const FieldCoderConfig* samFieldCoderConfig(uint32_t fieldIdx)
{
    if (fieldIdx >= SAM_FIELD_COUNT_SELECT) {
        return nullptr;
    }
    return &kSamFieldCoderConfig[fieldIdx];
}

/* 某字段是否把指定编码器列为候选。 */
inline bool samFieldCandidate(uint32_t fieldIdx, CoderType type)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    if (cfg == nullptr) {
        return false;
    }
    return std::find(cfg->candidates.begin(), cfg->candidates.end(), type) != cfg->candidates.end();
}

/* 默认编码器；未登记时返回调用方给的兜底。 */
inline CoderType samFieldDefaultCoder(uint32_t fieldIdx, CoderType fallback)
{
    const FieldCoderConfig* cfg = samFieldCoderConfig(fieldIdx);
    return (cfg != nullptr) ? cfg->fallback : fallback;
}
