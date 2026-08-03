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

class coder_fcv2 {
public:
    /*
     * freqTable 传质量值的出现次数，下标是符号值（原始字节减去 '!'），值是计数。
     * 用途是建哈夫曼树：质量值分布高度倾斜，让高频符号的编码路径短一些，可以显著
     * 减少每个符号需要的二值编码次数，进而提升速度。注意这不影响压缩率——算术编码
     * 的代价只取决于预测概率，把一个符号拆成几步条件判断，代价总和是一样的。
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable);

    /*
     * 从先前导出的模型快照创建编码器。modelLoaded 非空时会写入实际加载结果：只有快照
     * 的版本、编译期尺寸参数、字母表和全部模型数组都通过校验时才为 true。任何损坏或不
     * 兼容的快照都会保留 freqTable 所构造的固定初值模型，绝不留下半恢复的状态；这让
     * 调用方可以把快照当作可选的性能优化，而不把错误处理扩散到压缩主流程。
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable,
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
     */
    void encode_record(const uint8_t* qual, uint32_t len, bool rev);

    /* 编码收尾，返回写入 io 的字节数。 */
    int32_t encode_flush();

    /*
     * 解码一条记录的质量值到 dst。只需要长度，不需要再传链方向——它已经随码流写出，
     * 解码时自行读回。这样解压侧不必为了解 QUAL 去跟踪 FLAG 字段。
     */
    int32_t decode_record(uint8_t* dst, uint32_t len);

    /* 解码前调用一次，读取码流头部的字母表等信息。 */
    int32_t begin_decode();

private:
    std::unique_ptr<fcv2_impl> impl;

    coder_fcv2(const coder_fcv2&) = delete;
    coder_fcv2& operator=(const coder_fcv2&) = delete;
};
