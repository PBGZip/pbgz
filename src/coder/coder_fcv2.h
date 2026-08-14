/*
 * coder_fcv2.h - 质量值上下文混合编码器
 *
 * 针对 SAM 质量值列设计的专用编码器。它不是通用的字节流压缩器，依赖两项只有比对后
 * 的 SAM 才能提供的信息：每条记录的长度，以及该记录的链方向（FLAG 的 0x10 位）。
 * 因此它只能用于 SAM 的 QUAL 列，其他字段上不可用，见 supports()。
 *
 * 基本做法是给每个质量值符号建立多个不同粒度的概率模型，把它们的预测在对数几率域
 * 加权合并，再交给二值区间编码器。模型粒度从粗到细，粗的保证任何时候都有足够样本，
 * 细的在数据充分时提供额外精度，权重由梯度下降自行调整，不需要人工设定回退阈值。
 *
 * 关于接口粒度：对外按"一条记录"暴露，而不是按字节或按块。原因有两个。
 * 一是链方向和读长本来就是记录级的属性，按记录传最自然。二是实现全部藏在 .cpp 里
 * （原因见下），跨编译单元调用无法内联，按记录调用意味着这点开销被整条记录的字节数
 * 摊薄，而真正的热路径——逐比特的预测与更新——留在 .cpp 内部，照常内联。
 *
 * 关于为什么实现必须放在 .cpp：本编码器要用 coder/fc/rangecoder.h 里的二值区间编码器，
 * 而该文件定义的 RangeCoder 与 coder/clr.h 里的同名类是两个不同的东西（一个是二值的，
 * 一个是多符号频率的）。使用本编码器的 sam_actuator.cpp 同时需要 coder_qual.h，
 * 后者会引入 clr.h。两者出现在同一个编译单元就会重复定义，所以本头文件不能包含
 * rangecoder.h，实现只能放到 .cpp 中。
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

#include "coder_io.h"

class fcv2_impl;

/*
 * fcv2 上下文模型的参数档位。默认值即历史固定常量；QualSelector 按数据特征（质量值
 * 字母表大小、样本量、平均读长）选定一组参数参与试压，选中的档位随 PreprocessInfo
 * 传给编码端与先验训练，解码端从码流头部读回同一组参数（见 encode_record/begin_decode）。
 *
 * 含义：
 *   cycleMax   / cycleBucket  记录内位置（测序循环序号）的取值上限与分档数
 *   deltaMax   / deltaBucket  read 内质量跳变次数的取值上限与分档数（策略 1）
 *   prevShift  m3 里前驱质量值的量化右移位数，越大越粗（字母表大时上下文更稀疏）
 *   useDelta   是否启用 m5 跳变次数上下文；置 false 时 deltaBucket 归一为 1
 *   useDedup   是否启用相邻重复 read 去重（策略 3，fqzcomp 的 do_dedup）：
 *              每条记录先与上一条比对，完全相同就只写 1 bit，整条质量串跳过。
 *   useQa      是否启用 read 平均质量分 bin 上下文（策略 4，fqzcomp 的 do_qa）：
 *              每条记录先算平均质量、量化成 4 档写进码流，档位作为 m6 上下文。
 */
struct Fcv2Cfg {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  deltaMax    = 32;
    int  deltaBucket = 8;
    int  prevShift   = 1;
    bool useDelta    = true;
    bool useDedup    = false;
    bool useQa       = false;

    bool operator==(const Fcv2Cfg& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               deltaMax == o.deltaMax && deltaBucket == o.deltaBucket &&
               prevShift == o.prevShift && useDelta == o.useDelta &&
               useDedup == o.useDedup && useQa == o.useQa;
    }
    bool operator!=(const Fcv2Cfg& o) const { return !(*this == o); }
};

class coder_fcv2 {
public:
    /*
     * freqTable 传质量值的出现次数，下标是符号值（原始字节减去 '!'），值是计数。
     * 用途是建哈夫曼树：质量值分布高度倾斜，让高频符号的编码路径短一些，可以显著
     * 减少每个符号需要的二值编码次数，进而提升速度。注意这不影响压缩率——算术编码
     * 的代价只取决于预测概率，把一个符号拆成几步条件判断，代价总和是一样的。
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable);

    /* 带上下文参数档位的版本；cfg 会原样写进码流头部供解码端读回。 */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg);

    /*
     * 从先前导出的模型快照创建编码器。modelLoaded 非空时会写入实际加载结果：只有快照
     * 的版本、编译期尺寸参数、字母表和全部模型数组都通过校验时才为 true。任何损坏或不
     * 兼容的快照都会保留 freqTable 所构造的固定初值模型，绝不留下半恢复的状态；这让
     * 调用方可以把快照当作可选的性能优化，而不把错误处理扩散到压缩主流程。
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);

    /* 带上下文参数档位 + 先验快照的版本。 */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);
    ~coder_fcv2();

    /*
     * 导出当前学习到的模型。快照同时携带字母表和量化频率，而不是只保存计数器：计数器
     * 的下标是哈夫曼内部节点号，节点号只能由完全相同的字母表和频率重建。返回 false
     * 表示当前编码器没有可导出的有效字母表，或输出缓冲区无法分配。
     */
    bool export_model(std::vector<uint8_t>& out) const;

    /*
     * 本编码器需要每条记录的长度和链方向，这两样只有比对后的 SAM 的 QUAL 列能提供。
     *
     * 适用范围的判断放在 CoderFactory::coderSupports 里，而不是做成本类的成员。
     * 原因是那个判断要比较 BlockType 和 SamField，它们定义在 src 层；coder 层的
     * 编译目标只包含 coder/ 目录，反向依赖上层会破坏现有的分层。
     */

    /*
     * 编码一条记录的质量值。
     *
     * rev 取自该记录 FLAG 的 0x10 位。按 SAM 规范，该位置位时 SEQ 和 QUAL 在文件里
     * 存的是相对参考正链的形式，也就是相对测序仪原始读出顺序已经反转，测序仪的第一个
     * 循环对应存储中的最后一个字节。编码器据此还原真实的循环序号——质量值随测序循环
     * 系统性下降，这个规律只有在循环序号还原正确时才能用上。实测还原与不还原相差
     * 0.37 个百分点。
     *
     * 未比对的数据（FASTQ 转来的 uBAM）该位恒为 0，此时存储顺序就是原始顺序，
     * 传 false 即可。
     *
     * seq 是对应位置的碱基序列（ACGTN），seqLen 是其长度，用作第五个上下文模型的
     * 条件（当前碱基 + 循环序号）。反向链时按与循环序号一致的映射取碱基，即取存储序里
     * 位置 len-1-i 的碱基，因为那才是产生本质量值的那个循环。seq 传 nullptr 或
     * seqLen 不足以覆盖 mapped 位置时，该位置碱基上下文归入"未知"档，行为与无碱基
     * 版本一致。
     */
    void encode_record(const uint8_t* qual, uint32_t len, bool rev,
                       const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* 编码收尾，返回写入 io 的字节数。 */
    int32_t encode_flush();

    /*
     * 解码一条记录的质量值到 dst。链方向由码流自带、解码时自行读回，只需要长度，
     * 这样解压侧不必为了解 QUAL 去跟踪 FLAG 字段。
     *
     * seq 是已经解出的对应记录碱基序列（长度 seqLen），作为碱基上下文模型的条件；
     * 传 nullptr 时归入"未知"档，必须与编码侧一致（编码侧也应传空，否则两端上下文
     * 不一致）。
     */
    int32_t decode_record(uint8_t* dst, uint32_t len,
                          const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* 解码前调用一次，读取码流头部的字母表等信息。 */
    int32_t begin_decode();

private:
    std::unique_ptr<fcv2_impl> impl;

    coder_fcv2(const coder_fcv2&) = delete;
    coder_fcv2& operator=(const coder_fcv2&) = delete;
};
