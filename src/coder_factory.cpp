/*
 * coder_factory.cpp - 编码器工厂的实现
 *
 * 实现放在 .cpp 而不是头文件里，是为了不让具体编码器的头文件泄漏给调用方：
 * coder_fc.h 会间接引入 clr.h，其中定义了一个 RangeCoder；而 coder_qual.h 引入的
 * qual_model.h 里另有一个同名但接口不同的 RangeCoder。两者同时出现在一个编译单元
 * 就会冲突。执行器普遍需要 coder_qual.h，所以工厂头文件必须保持干净。
 *
 * 关于 coder_simple_rc：它已从编解码两侧全部移除。
 *
 * 实测这个编码器是有损的——单独测它的往返时，10 个数据块全部校验失败，游程长度信息
 * 丢了将近四成，它那个看起来很漂亮的压缩率正是丢数据换来的。
 *
 * 在接通预处理选择之前，执行器把编码器类型写死，它从来没有被真正调用过，所以问题
 * 一直没暴露，也不存在用它压出来的历史文件。接通之后风险就变实了：试压只比较压缩后
 * 的大小、并不验证能否原样解回来，只要它在某个字段上压得最小就会被选中，产出一个
 * 解不出原始数据的文件。既然没有历史包袱，编码解码两侧一并去掉最干净。
 */

#include "coder_factory.h"

#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"

std::shared_ptr<coder> CoderFactory::makeEncoder(CoderType type, coder_io* io)
{
    switch (type) {
    case CoderType::FC:
        return std::make_shared<coder_fc>(io);

    /*
     * QUAL 走兜底：coder_qual 不继承 coder 基类，构造时还需要额外的频率表，没法由
     * 本工厂统一创建，质量字段有自己单独的压缩函数。通用字段不应该被选成 QUAL
     * （试压候选里根本没有它），万一出现就说明选择结果异常，退回 BWT_CM 是安全的。
     */
    case CoderType::QUAL:
    case CoderType::SIMPLE_RC:
    case CoderType::BWT_CM:
    default:
        return std::make_shared<coder_bwt_cm>(io);
    }
}

std::shared_ptr<coder> CoderFactory::makeDecoder(const std::string& magic, coder_io* io)
{
    if (magic == "coder_bwt_cm") {
        return std::make_shared<coder_bwt_cm>(io);
    }
    if (magic == "coder_fc") {
        return std::make_shared<coder_fc>(io);
    }
    return nullptr;
}

bool CoderFactory::coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx)
{
    if (type == CoderType::FCV2) {
        return fileType == (uint32_t)SAM && fieldIdx == (uint32_t)SAM_QUAL;
    }
    return true;
}

bool CoderFactory::canMake(CoderType type)
{
    return type == CoderType::BWT_CM || type == CoderType::FC;
}
