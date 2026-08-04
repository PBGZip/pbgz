/*
 * QUAL 列编码器的横向实测（仅用于评估，不参与产品构建）。
 *
 * 要回答的问题：质量值这一列到底该用哪个编码器。
 *
 * 现状是 QualSelector 只在 coder_qual 和 coder_fcv2 之间二选一，coder_bwt_cm
 * 压根不在候选集里——通用字段走的 CodecSelector 才会试 bwt_cm。所以 fcv2 目前
 * 只是"在两个候选里赢了"，不是"在全部候选里赢了"。这个缺口必须用数据补上：
 * 如果 bcm 在 QUAL 上更好，那么先验和分片这一整套就是对着错的编码器在做。
 *
 * 三个编码器都按记录逐条喂，与 sam_actuator 里的实际调用方式一致：
 * bwt_cm 和 qual 用 encode_line / encode_qual_gen2，fcv2 用 encode_record。
 * 之前一次性把几 MB 塞进 encode_line 是错的用法，那不是它的接口约定。
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <chrono>
#include <string>
#include <vector>

#include "coder/coder.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_qual.h"

namespace {

/*
 * coder 层的分配、释放、退出、日志四个回调必须先注册，否则用到 safe_alloc 的编码器
 * （coder_bwt_cm 初始化时要一次性拿约 96 MB：16 MB 缓冲 + 16 MB BWT + 64 MB 索引数组）
 * 会拿到空指针，紧接着 memcpy 就是段错误——进程无声无息地消失。
 * 注意 safe_alloc 走的是 alloc_proc，与 realloc_proc 是两个独立的回调，
 * 只注册后者不够。产品路径上这套注册在 PbgzEngine::init 里做，独立测试程序
 * 必须自己补上，否则测的不是编码器而是缺失的初始化。
 *
 * coder_fcv2 不依赖 safe_alloc，所以之前几个只测 fcv2 的工具没暴露这个问题。
 */
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
    coder_ns::resister_logger_proc([](int, const char*) { });
    coder_ns::initFcCoder();
}

struct Record {
    std::string qual;
    std::string seq;
    bool rev;
};

bool loadRecords(const char* path, size_t takeBytes,
                 std::vector<Record>& out, std::vector<uint32_t>& freq)
{
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) return false;

    freq.assign(256, 0);
    size_t taken = 0;
    char* line = NULL;
    size_t cap = 0;

    while (taken < takeBytes && getline(&line, &cap, fp) > 0) {
        if (line[0] == '@') continue;

        int col = 0;
        char* cur = line;
        char* flagCol = NULL;
        char* seqCol = NULL; size_t seqLen = 0;
        char* qualCol = NULL; size_t qualLen = 0;
        for (char* p = line; ; p++) {
            if (*p == '\t' || *p == '\n' || *p == '\0') {
                if (col == 1) flagCol = cur;
                if (col == 9) { seqCol = cur; seqLen = (size_t)(p - cur); }
                if (col == 10) { qualCol = cur; qualLen = (size_t)(p - cur); }
                col++;
                cur = p + 1;
                if (*p == '\n' || *p == '\0') break;
            }
        }
        if (qualCol == NULL || qualLen == 0) continue;

        Record r;
        r.qual.assign(qualCol, qualLen);
        r.seq.assign(seqCol != NULL ? seqCol : "", seqCol != NULL ? seqLen : 0);
        r.rev = (flagCol != NULL) && ((atoi(flagCol) & 0x10) != 0);
        for (size_t i = 0; i < qualLen; i++) freq[(uint8_t)r.qual[i]]++;
        out.push_back(r);
        taken += qualLen;
    }
    free(line);
    fclose(fp);
    return !out.empty();
}

/*
 * coder_qual 要的是按出现次数降序排列的 <符号, 计数> 表，符号是原始字节减 '!'。
 * 计数字段是 uint16_t，样本里的真实计数远超它的量程，所以按名次折算成一个单调递减
 * 的权重即可——这个表只用来定字母表顺序，不参与概率计算。sam_actuator 里构造
 * fcv2 的频率表时也是同样的处理。
 */
std::vector<std::pair<uint16_t, uint16_t> > buildQualFreqTable(const std::vector<uint32_t>& freq)
{
    std::vector<std::pair<uint32_t, uint32_t> > tmp;
    for (int b = 0; b < 256; b++) {
        if (freq[b] > 0) {
            tmp.push_back(std::make_pair((uint32_t)(b - '!'), freq[b]));
        }
    }
    for (size_t i = 0; i < tmp.size(); i++) {
        for (size_t j = i + 1; j < tmp.size(); j++) {
            if (tmp[j].second > tmp[i].second) std::swap(tmp[i], tmp[j]);
        }
    }
    std::vector<std::pair<uint16_t, uint16_t> > table;
    for (size_t i = 0; i < tmp.size(); i++) {
        table.push_back(std::make_pair((uint16_t)tmp[i].first,
                                       (uint16_t)(tmp.size() - i)));
    }
    return table;
}

void report(const char* name, size_t raw, size_t packed, double sec)
{
    if (packed == 0) {
        printf("  %-14s 失败\n", name);
        return;
    }
    printf("  %-14s %10zu -> %9zu  (%.4f%%)  %.2fs  %.1f MB/s\n",
           name, raw, packed, 100.0 * (double)packed / (double)raw, sec,
           (double)raw / (1024.0 * 1024.0) / sec);
}

} /* namespace */

int main(int argc, char** argv)
{
    if (argc < 3) {
        fprintf(stderr, "用法: %s <sam文件> <QUAL字节数>\n", argv[0]);
        return 1;
    }
    const char* samPath = argv[1];
    const size_t takeBytes = (size_t)atoll(argv[2]);

    registerCoderCallbacks();

    std::vector<Record> recs;
    std::vector<uint32_t> freq;
    if (!loadRecords(samPath, takeBytes, recs, freq)) {
        fprintf(stderr, "读取失败\n");
        return 1;
    }

    size_t raw = 0;
    for (size_t i = 0; i < recs.size(); i++) raw += recs[i].qual.size();
    printf("记录数 %zu，QUAL 原始 %zu 字节\n\n", recs.size(), raw);
    fflush(stdout);

    const size_t cap = raw * 2 + (1 << 20);

    /* coder_bwt_cm：按记录逐条 encode_line，与 sam_actuator 里的用法一致。 */
    {
        printf("  coder_bwt_cm 运行中...\n"); fflush(stdout);
        std::vector<uint8_t> out(cap, 0);
        coder_io io(out.data(), (int32_t)out.size());
        coder_bwt_cm coder(&io);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < recs.size(); i++) {
            coder.encode_line((uint8_t*)recs[i].qual.data(), (uint32_t)recs[i].qual.size());
        }
        coder.encode_flush();
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        report("coder_bwt_cm", raw, io.data_len > 0 ? (size_t)io.data_len : 0, sec);
    }

    {
        printf("  coder_qual 运行中...\n"); fflush(stdout);
        std::vector<uint8_t> out(cap, 0);
        coder_io io(out.data(), (int32_t)out.size());
        std::vector<std::pair<uint16_t, uint16_t> > table = buildQualFreqTable(freq);
        coder_qual coder(&io, true, table);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < recs.size(); i++) {
            coder.encode_qual_gen2((uint8_t*)recs[i].seq.data(),
                                   (uint8_t*)recs[i].qual.data(),
                                   (uint32_t)recs[i].qual.size());
        }
        coder.encode_flush();
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        report("coder_qual", raw, io.data_len > 0 ? (size_t)io.data_len : 0, sec);
    }

    {
        printf("  coder_fcv2 运行中...\n"); fflush(stdout);
        std::vector<uint8_t> out(cap, 0);
        coder_io io(out.data(), (int32_t)out.size());
        coder_fcv2 coder(&io, freq);
        const auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < recs.size(); i++) {
            coder.encode_record((const uint8_t*)recs[i].qual.data(),
                                (uint32_t)recs[i].qual.size(), recs[i].rev);
        }
        const int32_t n = coder.encode_flush();
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        report("coder_fcv2", raw, n > 0 ? (size_t)n : 0, sec);
    }

    return 0;
}
