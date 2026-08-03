/*
 * 先验收益实测工具（仅用于评估，不参与产品构建）。
 *
 * 要回答的问题只有一个：pbgz 树内这份 fcv2，先验到底值多少个百分点？
 *
 * benchmark 目录里那个 1.59pp 是另一套独立实现测出来的，模型结构、参数、字母表处理
 * 都未必一致，不能直接搬来当作接入决策的依据。在把先验块、相对偏移、共享缓存这一整套
 * 管道铺开之前，必须先在真实数据上把这个数字量出来；量不到就没有必要做后面的工作。
 *
 * 三种模式对照：
 *   A. 无先验、整段一个分片        —— 当前线上行为的等价物
 *   B. 无先验、按记录数切独立分片  —— 只加随机访问，不加先验，代价有多大
 *   C. 有先验、按记录数切独立分片  —— 目标形态
 *
 * 关键点：训练用文件靠前的一段，测量用**训练区之后**的另一段。否则等于拿训练数据
 * 考自己，先验的收益会被严重高估。
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "coder/coder_fcv2.h"

namespace {

struct QualRecord {
    std::vector<uint8_t> qual;
    bool rev;
};

/*
 * 从 SAM 里顺序读取 QUAL 记录。skipBytes 指定先跳过多少字节的 QUAL 再开始收集，
 * 用来把训练区和测量区分开。takeBytes 是收集的上限。
 */
bool loadRecords(const char* path, size_t skipBytes, size_t takeBytes,
                 std::vector<QualRecord>& out, std::vector<uint32_t>& freq)
{
    FILE* fp = fopen(path, "rb");
    if (fp == NULL) return false;

    freq.assign(256, 0);
    size_t skipped = 0;
    size_t taken = 0;
    char* line = NULL;
    size_t cap = 0;
    ssize_t len = 0;

    while (taken < takeBytes && (len = getline(&line, &cap, fp)) > 0) {
        if (line[0] == '@') continue;

        int col = 0;
        char* cur = line;
        char* flagCol = NULL;
        char* qualCol = NULL;
        size_t qualLen = 0;
        for (char* p = line; ; p++) {
            if (*p == '\t' || *p == '\n' || *p == '\0') {
                if (col == 1) flagCol = cur;
                if (col == 10) { qualCol = cur; qualLen = (size_t)(p - cur); }
                col++;
                cur = p + 1;
                if (*p == '\n' || *p == '\0') break;
            }
        }
        if (qualCol == NULL || qualLen == 0) continue;

        if (skipped < skipBytes) { skipped += qualLen; continue; }

        QualRecord rec;
        rec.qual.assign(qualCol, qualCol + qualLen);
        rec.rev = (flagCol != NULL) && ((atoi(flagCol) & 0x10) != 0);
        for (size_t i = 0; i < qualLen; i++) freq[rec.qual[i]]++;
        out.push_back(rec);
        taken += qualLen;
    }
    free(line);
    fclose(fp);
    return !out.empty();
}

size_t totalQualBytes(const std::vector<QualRecord>& recs)
{
    size_t n = 0;
    for (size_t i = 0; i < recs.size(); i++) n += recs[i].qual.size();
    return n;
}

/*
 * 按 sliceRecords 条一个分片压缩整段数据，每个分片都新建一个编码器。
 *
 * 新建编码器就意味着分片之间没有任何状态传递，这正是随机访问要求的独立性。
 * sliceRecords 为 0 表示不切分片，整段用一个编码器压完。
 * blob 非空时每个分片都从同一份先验初始化。
 *
 * 返回所有分片压缩结果的字节总和。分片各自带自己的流头，这部分开销也计入，
 * 因为真实格式里同样躲不掉。
 */
size_t compressAll(const std::vector<QualRecord>& recs,
                   const std::vector<uint32_t>& freq,
                   size_t sliceRecords,
                   const std::vector<uint8_t>* blob,
                   size_t* sliceCountOut,
                   size_t* priorFailures)
{
    const size_t step = (sliceRecords == 0) ? recs.size() : sliceRecords;
    size_t total = 0;
    size_t slices = 0;

    for (size_t begin = 0; begin < recs.size(); begin += step) {
        const size_t end = (begin + step < recs.size()) ? begin + step : recs.size();

        size_t rawLen = 0;
        for (size_t i = begin; i < end; i++) rawLen += recs[i].qual.size();

        std::vector<uint8_t> buf(rawLen * 2 + (1 << 16), 0);
        coder_io io(buf.data(), (int32_t)buf.size());

        int32_t produced = 0;
        if (blob != NULL) {
            bool loaded = false;
            coder_fcv2 coder(&io, freq, *blob, &loaded);
            if (!loaded && priorFailures != NULL) (*priorFailures)++;
            for (size_t i = begin; i < end; i++) {
                coder.encode_record(recs[i].qual.data(), (uint32_t)recs[i].qual.size(), recs[i].rev);
            }
            produced = coder.encode_flush();
        } else {
            coder_fcv2 coder(&io, freq);
            for (size_t i = begin; i < end; i++) {
                coder.encode_record(recs[i].qual.data(), (uint32_t)recs[i].qual.size(), recs[i].rev);
            }
            produced = coder.encode_flush();
        }

        if (produced <= 0) {
            fprintf(stderr, "分片 %zu 压缩失败\n", slices);
            return 0;
        }
        total += (size_t)produced;
        slices++;
    }

    if (sliceCountOut != NULL) *sliceCountOut = slices;
    return total;
}

void report(const char* label, size_t raw, size_t comp, size_t slices)
{
    printf("  %-34s %10zu -> %9zu  (%.4f%%)  分片数 %zu\n",
           label, raw, comp, 100.0 * (double)comp / (double)raw, slices);
}

} /* namespace */

int main(int argc, char** argv)
{
    if (argc < 5) {
        fprintf(stderr,
                "用法: %s <sam文件> <训练QUAL字节数> <测量QUAL字节数> <分片记录数>\n"
                "例:   %s /tmp/con_sorted.sam 9000000 45000000 10000\n",
                argv[0], argv[0]);
        return 1;
    }
    const char* samPath = argv[1];
    const size_t trainBytes = (size_t)atoll(argv[2]);
    const size_t evalBytes = (size_t)atoll(argv[3]);
    const size_t sliceRecords = (size_t)atoll(argv[4]);

    /* 第一步：读训练区，训练并导出先验。 */
    std::vector<QualRecord> trainRecs;
    std::vector<uint32_t> trainFreq;
    if (!loadRecords(samPath, 0, trainBytes, trainRecs, trainFreq)) {
        fprintf(stderr, "训练区读取失败\n");
        return 1;
    }
    const size_t trainRaw = totalQualBytes(trainRecs);

    std::vector<uint8_t> blob;
    {
        std::vector<uint8_t> buf(trainRaw * 2 + (1 << 20), 0);
        coder_io io(buf.data(), (int32_t)buf.size());
        coder_fcv2 coder(&io, trainFreq);
        for (size_t i = 0; i < trainRecs.size(); i++) {
            coder.encode_record(trainRecs[i].qual.data(),
                                (uint32_t)trainRecs[i].qual.size(), trainRecs[i].rev);
        }
        coder.encode_flush();
        if (!coder.export_model(blob)) {
            fprintf(stderr, "导出先验失败\n");
            return 1;
        }
    }

    /* 第二步：读测量区，必须跳过训练区，避免拿训练数据考自己。 */
    std::vector<QualRecord> evalRecs;
    std::vector<uint32_t> evalFreq;
    if (!loadRecords(samPath, trainRaw, evalBytes, evalRecs, evalFreq)) {
        fprintf(stderr, "测量区读取失败\n");
        return 1;
    }
    const size_t evalRaw = totalQualBytes(evalRecs);

    printf("训练区: %zu 条 / %zu 字节\n", trainRecs.size(), trainRaw);
    printf("测量区: %zu 条 / %zu 字节（已跳过训练区）\n", evalRecs.size(), evalRaw);
    printf("先验快照原始大小: %zu 字节\n", blob.size());
    printf("分片记录数: %zu\n\n", sliceRecords);

    size_t slicesA = 0, slicesB = 0, slicesC = 0;
    size_t priorFail = 0;

    const size_t compA = compressAll(evalRecs, evalFreq, 0, NULL, &slicesA, NULL);
    const size_t compB = compressAll(evalRecs, evalFreq, sliceRecords, NULL, &slicesB, NULL);
    const size_t compC = compressAll(evalRecs, evalFreq, sliceRecords, &blob, &slicesC, &priorFail);

    if (compA == 0 || compB == 0 || compC == 0) return 1;

    printf("测量结果:\n");
    report("A 无先验 / 整段一片", evalRaw, compA, slicesA);
    report("B 无先验 / 独立分片", evalRaw, compB, slicesB);
    report("C 有先验 / 独立分片", evalRaw, compC, slicesC);

    const double ratioA = 100.0 * (double)compA / (double)evalRaw;
    const double ratioB = 100.0 * (double)compB / (double)evalRaw;
    const double ratioC = 100.0 * (double)compC / (double)evalRaw;

    printf("\n结论:\n");
    printf("  分片的代价 (B-A): %+.4f pp\n", ratioB - ratioA);
    printf("  先验的收益 (B-C): %+.4f pp\n", ratioB - ratioC);
    printf("  相对当前行为 (C-A): %+.4f pp\n", ratioC - ratioA);
    if (priorFail > 0) {
        printf("  警告: 有 %zu 个分片未能加载先验\n", priorFail);
    }
    return 0;
}
