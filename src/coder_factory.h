/*
 * coder_factory.h - 按类型或 magic 创建编码器/解码器
 *
 * 在此之前，各个执行器都是直接 new 出具体的编码器类，编码器类型在代码里写死。
 * 预处理阶段（CodecSelector）虽然会试压出每个字段的最优编码器并存进 PreprocessInfo，
 * 但执行器从不去读它，选择结果实际上没有生效。
 *
 * 这个工厂把"类型 -> 实例"的映射集中到一处，让执行器可以按预处理的选择结果动态
 * 创建编码器，同时保证任何异常情况下都有可用的兜底。
 *
 * 注意本头文件刻意不包含任何具体编码器的头文件，实现全部放在 coder_factory.cpp。
 * 原因是 coder_fc.h 会间接引入 clr.h，而 coder_qual.h 引入的 qual_model.h 里另有
 * 一个同名但接口不同的 RangeCoder，两者出现在同一个编译单元就会冲突。把实现挪进
 * .cpp 之后，包含本头文件的调用方只会看到声明，不会被这层依赖波及。
 */

#pragma once

#include <memory>
#include <string>

#include "preprocess_info.h"
#include "coder/coder.h"

class CoderFactory {
public:
    /*
     * 压缩侧：按预处理选出的类型创建编码器。
     *
     * 这里必须永远返回一个可用的编码器，不能返回空指针。压缩是单向过程，一旦某个
     * 字段因为拿不到编码器而被跳过，产出的文件就是残缺的，问题会推迟到解压时才暴露。
     * 所以无法识别的类型（包括将来新增、但本版本还不认识的类型）一律退回 BWT_CM，
     * 它是所有字段都能处理的通用编码器，压缩率未必最优但一定正确。
     *
     * CoderType::QUAL 也走这条兜底路径：coder_qual 不继承 coder 基类，构造时还需要
     * 额外的频率表，没法由本工厂统一创建，质量字段有自己单独的压缩函数。通用字段
     * 不应该被选成 QUAL（试压候选里根本没有它），万一出现就说明选择结果异常，
     * 退回 BWT_CM 是安全的。
     */
    static std::shared_ptr<coder> makeEncoder(CoderType type, coder_io* io);

    /*
     * 解压侧：按码流里记录的 magic 创建解码器。
     *
     * 这里和压缩侧相反，必须返回空指针而不是兜底。解压时 magic 是写在文件里的事实，
     * 表明这段数据当初是用哪个编码器压的。如果不认识这个 magic，说明文件是更新版本
     * 的 pbgz 写的，或者数据已经损坏。此时随便挑一个编码器去解，只会解出一堆垃圾，
     * 而且很可能不报错——这比直接失败危险得多。
     *
     * 调用方拿到空指针后必须报错中止，不能静默跳过。
     */
    static std::shared_ptr<coder> makeDecoder(const std::string& magic, coder_io* io);

    /*
     * 该类型能否由本工厂创建。
     *
     * 预处理阶段用它来过滤试压候选：如果某个类型工厂造不出来，就算试压赢了也没法
     * 在压缩时被真正使用，不如一开始就别参与比较。
     */
    static bool canMake(CoderType type);

    /*
     * 某个编码器能否用于指定文件类型的指定字段。
     *
     * 绝大多数编码器是通用字节流压缩器，放哪个字段上都能跑。少数有额外前置依赖：
     * fcv2 需要每条记录的长度和链方向，只有比对后的 SAM 的 QUAL 列能提供，用在别处
     * 拿不到这些信息。预处理试压前先问一句，不适用就不参与比较，免得选出一个实际
     * 用不了的编码器。
     *
     * 这个判断放在这里而不是各编码器自己身上，是因为它要比较 BlockType 和 SamField，
     * 而 coder 层的编译目标不包含 src 目录，不能反向依赖。
     */
    static bool coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx);
};
