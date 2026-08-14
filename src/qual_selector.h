/*
 * qual_selector.h - 质量值列的编码器评估
 *
 * 质量值列不能走通用的评估路径，原因有三个。
 *
 * 一是候选不同。通用字段的候选是 coder_bwt_cm 和 coder_fc 这类字节流压缩器，而质量值
 * 实际使用的是 coder_qual 和 fcv2。此前预处理一直在拿 bwt_cm 和 fc 评估质量值列，
 * 而 compressQuality 用的是另外两个，评估结果实际上没有被使用。
 *
 * 二是输入形态不同。通用评估把整列拼成一段连续字节，记录边界就丢了；而这两个候选都
 * 需要记录级信息——fcv2 需要每条记录的长度来还原测序循环序号，coder_qual 需要对应的
 * 碱基序列作为上下文。
 *
 * 三是编译约束。coder_qual.h 会间接引入 clr.h，其中的 RangeCoder 与 coder_fc.h 引入的
 * fc/rangecoder.h 里的同名类冲突，两者不能出现在同一个编译单元，所以质量值的评估必须
 * 单独成一个文件。
 */

#pragma once

#include <stdint.h>
#include <vector>
#include <string>

#include "preprocess_info.h"

/* 采样得到的一条记录，评估两个候选都需要这些信息。 */
struct QualSampleRecord {
    std::string qual;   /* 质量值 */
    std::string seq;    /* 对应的碱基序列，coder_qual 用作上下文 */
    bool        rev;    /* 链方向，取自 FLAG 的 0x10 位，fcv2 用来还原循环序号 */
};

class QualSelector {
public:
    /*
     * 在采样记录上试压 coder_qual、fcv2（多个上下文参数档位）与 bwt_cm，返回压缩后
     * 更小的那个。freqByByte 是质量值的出现次数，下标为原始字节值。三个候选都需要它：
     * coder_qual 用它建自己的频率表，fcv2 用它建哈夫曼树。
     *
     * 评估走多轮收敛：从 64KB 起步逐轮翻倍样本量，领先者拉开 3% 以上差距即定案
     * （见 codec_selector.cpp 的 SETTLE_MARGIN）。选中的 fcv2 参数档位（若有）随
     * FieldCodecSelection.fcv2Params 返回，压缩端与先验训练据此保持一致；解码端从
     * fcv2 码流头部读回同一组参数。
     *
     * 样本太少时返回 SKIPPED，由调用方沿用默认编码器。评估失败（所有候选都跑不通）
     * 时返回 FAILED，同样退回默认，绝不会因为评估出问题而让压缩无法进行。
     */
    static FieldCodecSelection select(const std::vector<QualSampleRecord>& records,
                                      const std::vector<uint32_t>& freqByByte);
};
