/*
 * coder_bwt_cm 逐记录往返验证（用于定位 QUAL 走 bwt_cm 时解压失败的根因）。
 *
 * 要回答的问题：把 QUAL 按记录逐条 encode_line(ptr, len) 喂进 bwt_cm、
 * 再按同样的长度逐条 decode_line(out, len) 取回，能不能逐字节还原。
 *
 * 之所以必须先单测编码器本身，是因为端到端失败（解出 64 MB vs 原 19 MB）
 * 有两种可能：编码器的定长流接口本身不成立，或者 sam_actuator 的接线错了。
 * 这两者的修法完全不同，靠读代码分不出来，必须用一个不含任何 pbgz 上下文的
 * 最小用例把编码器单独钉死。
 *
 * 关键点：decode_line 传 split_ch = UINT8_MAX 时走的是纯定长流分支，
 * 不看分隔符，只按调用方给的 out_len 取字节。encode 侧无分隔符喂入
 * 与它是配套的，所以理论上应当成立——本测试就是验证这个"理论上"。
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "coder/coder.h"
#include "coder/coder_bwt_cm.h"

namespace {

void registerCoderCallbacks()
{
    coder_ns::register_alloc_proc([](size_t size) -> uint8_t* {
        return (uint8_t*)malloc(size);
    });
    coder_ns::register_realloc_proc(
        [](size_t& size, uint8_t* ptr, size_t newSize) -> uint8_t* {
            uint8_t* p = (uint8_t*)realloc(ptr, newSize);
            if (p != NULL) size = newSize;
            return p;
        });
    coder_ns::register_free_func([](void*& ptr) { free(ptr); ptr = NULL; });
    coder_ns::resister_logger_proc([](int, const char* msg) {
        fprintf(stderr, "coder 日志: %s\n", msg != NULL ? msg : "");
    });
    coder_ns::initFcCoder();
}

/*
 * 造一批长度不等的伪质量值记录。长度故意不统一，因为真实 SAM 里 QUAL 长度
 * 随 SEQ 变化，定长记录会掩盖掉"跨块边界时记录被劈开"这一类问题。
 */
struct Records {
    std::vector<uint8_t> flat;
    std::vector<uint32_t> lens;
};

Records makeRecords(uint32_t count, uint32_t baseLen)
{
    Records r;
    uint32_t seed = 12345;
    for (uint32_t i = 0; i < count; ++i) {
        uint32_t len = baseLen + (i % 37);
        r.lens.push_back(len);
        for (uint32_t j = 0; j < len; ++j) {
            seed = seed * 1103515245u + 12345u;
            /* 质量值集中在少数几个符号上，贴近真实分布，否则压不动看不出问题 */
            static const char alphabet[] = "!,:FI";
            r.flat.push_back((uint8_t)alphabet[(seed >> 16) % 5]);
        }
    }
    return r;
}

/* 返回 0 表示通过 */
int runCase(const char* name, uint32_t count, uint32_t baseLen, int level)
{
    Records rec = makeRecords(count, baseLen);
    const uint32_t srcLen = (uint32_t)rec.flat.size();

    /* 压缩缓冲给足，避免因为溢出把问题误判成编码器缺陷 */
    std::vector<uint8_t> comp(srcLen * 2 + (1u << 20), 0);

    uint32_t dstLen = 0;
    {
        coder_io io(comp.data(), (uint32_t)comp.size());
        io.meta["level"] = (Json::Value::Int)level;
        coder_bwt_cm enc(&io);
        uint32_t off = 0;
        for (uint32_t i = 0; i < rec.lens.size(); ++i) {
            enc.encode_line(rec.flat.data() + off, rec.lens[i]);
            off += rec.lens[i];
        }
        enc.encode_flush();
        dstLen = io.data_len;
    }

    std::vector<uint8_t> out(srcLen + 1024, 0);
    uint32_t outOff = 0;
    int bad = 0;
    {
        coder_io io(comp.data(), dstLen);
        coder_bwt_cm dec(&io);
        for (uint32_t i = 0; i < rec.lens.size(); ++i) {
            int32_t got = dec.decode_line(out.data() + outOff, rec.lens[i]);
            if (got != (int32_t)rec.lens[i]) {
                fprintf(stderr,
                        "  [%s] 第 %u 条记录长度不符: 期望 %u, 实际 %d\n",
                        name, i, rec.lens[i], got);
                bad = 1;
                break;
            }
            outOff += rec.lens[i];
        }
    }

    if (!bad && outOff != srcLen) {
        fprintf(stderr, "  [%s] 总长不符: 期望 %u, 实际 %u\n", name, srcLen, outOff);
        bad = 1;
    }
    if (!bad && memcmp(out.data(), rec.flat.data(), srcLen) != 0) {
        for (uint32_t i = 0; i < srcLen; ++i) {
            if (out[i] != rec.flat[i]) {
                fprintf(stderr, "  [%s] 内容首个不一致位置 %u: 期望 0x%02X, 实际 0x%02X\n",
                        name, i, rec.flat[i], out[i]);
                break;
            }
        }
        bad = 1;
    }

    printf("%-28s level=%d  %u 条 / %u 字节 -> %u 字节 (%.2f%%)  %s\n",
           name, level, count, srcLen, dstLen,
           srcLen == 0 ? 0.0 : (double)dstLen * 100.0 / (double)srcLen,
           bad ? "FAIL" : "PASS");
    return bad;
}

} // namespace

int main()
{
    registerCoderCallbacks();

    int fail = 0;
    /* 小数据：全部落在 encode_flush 的单块路径上 */
    fail |= runCase("单块-小数据", 2000, 100, 4);
    /* 大数据：必须跨越 bsize，触发 encode_line 里的整块提交 + 记录被劈开 */
    fail |= runCase("跨块-level1(1MB块)", 40000, 150, 1);
    /* 多次跨块 */
    fail |= runCase("多次跨块-level1", 200000, 150, 1);
    /* 产品默认档位 */
    fail |= runCase("默认档-level4", 50000, 150, 4);

    printf("\n%s\n", fail ? "结论: 有用例失败" : "结论: 全部通过");
    return fail;
}
