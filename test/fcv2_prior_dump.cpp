/*
 * 先验快照落盘工具（仅用于评估，不参与产品构建）。
 *
 * 目的：用真实 SAM 的质量值训练 fcv2，把训练后的模型快照导出到文件，
 * 以便实测各种通用压缩器对这份快照的压缩效果。
 *
 * 之所以必须用真实数据而不是合成数据：快照的可压缩性几乎完全取决于
 * "有多少上下文从未被访问过、仍保持初值"，而这个分布只有真实质量值才能反映。
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>
#include <vector>

#include "coder/coder_fcv2.h"

int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "用法: %s <sam文件> <训练字节数> <输出快照文件>\n", argv[0]);
        return 1;
    }
    const char* samPath = argv[1];
    const size_t trainBytes = (size_t)atoll(argv[2]);
    const char* outPath = argv[3];

    FILE* fp = fopen(samPath, "rb");
    if (fp == NULL) {
        fprintf(stderr, "打不开 %s\n", samPath);
        return 1;
    }

    /* 逐行读 SAM，取第 11 列（QUAL）和第 2 列（FLAG，用于链方向）。 */
    std::vector<std::vector<uint8_t> > records;
    std::vector<bool> revs;
    std::vector<uint32_t> freq(256, 0);
    size_t collected = 0;

    char* line = NULL;
    size_t cap = 0;
    ssize_t len = 0;
    while (collected < trainBytes && (len = getline(&line, &cap, fp)) > 0) {
        if (line[0] == '@') continue;
        /* 按制表符切出各列 */
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

        const int flag = (flagCol != NULL) ? atoi(flagCol) : 0;
        std::vector<uint8_t> rec(qualCol, qualCol + qualLen);
        for (size_t i = 0; i < qualLen; i++) freq[rec[i]]++;
        records.push_back(rec);
        revs.push_back((flag & 0x10) != 0);
        collected += qualLen;
    }
    free(line);
    fclose(fp);

    if (records.empty()) {
        fprintf(stderr, "没读到任何 QUAL 记录\n");
        return 1;
    }

    /* 训练：正常走一遍编码，编码结果丢弃，只要训练后的模型状态。 */
    std::vector<uint8_t> scratch(collected * 2 + (1 << 20), 0);
    coder_io io(scratch.data(), (int32_t)scratch.size());
    coder_fcv2 coder(&io, freq);
    for (size_t i = 0; i < records.size(); i++) {
        coder.encode_record(records[i].data(), (uint32_t)records[i].size(), revs[i]);
    }
    const int32_t streamLen = coder.encode_flush();

    std::vector<uint8_t> blob;
    if (!coder.export_model(blob)) {
        fprintf(stderr, "导出模型快照失败\n");
        return 1;
    }

    FILE* out = fopen(outPath, "wb");
    if (out == NULL) {
        fprintf(stderr, "打不开输出 %s\n", outPath);
        return 1;
    }
    fwrite(blob.data(), 1, blob.size(), out);
    fclose(out);

    /* 统计字母表大小与仍保持初值的计数器占比，用于解释快照的可压缩性。 */
    int alpha = 0;
    for (int i = 0; i < 256; i++) if (freq[i] > 0) alpha++;

    printf("训练记录数: %zu\n", records.size());
    printf("训练 QUAL 字节: %zu\n", collected);
    printf("字母表大小: %d\n", alpha);
    printf("训练后码流长度: %d (%.2f%%)\n", streamLen,
           100.0 * (double)streamLen / (double)collected);
    printf("快照原始大小: %zu\n", blob.size());
    return 0;
}
