/*
 * fcv2 往返测试：用真实 SAM 数据验证 coder_fcv2 编解码一致。
 * 独立于 pbgz 主流程，确认编码器本身正确后再接入。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>

#include "coder/coder_fcv2.h"

struct Rec { std::string qual; bool rev; };

int main(int argc, char** argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <in.sam> [maxrec]\n", argv[0]); return 1; }
    long maxRec = (argc >= 3) ? atol(argv[2]) : 200000;

    FILE* f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "open failed\n"); return 1; }

    std::vector<Rec> recs;
    std::vector<uint32_t> freq(256, 0);
    char* line = (char*)malloc(1 << 20);
    size_t totalQual = 0;

    while (fgets(line, 1 << 20, f) && (long)recs.size() < maxRec) {
        if (line[0] == '@') continue;
        char* fields[12]; int nf = 0; fields[nf++] = line;
        for (char* p = line; *p && nf < 12; p++) {
            if (*p == '\t') { *p = '\0'; fields[nf++] = p + 1; }
        }
        if (nf < 11) continue;
        long flag = strtol(fields[1], NULL, 10);
        char* q = fields[10];
        int qn = 0;
        while (q[qn] && q[qn] != '\n' && q[qn] != '\r' && q[qn] != '\t') qn++;
        if (qn <= 0 || (qn == 1 && q[0] == '*')) continue;
        Rec r; r.qual.assign(q, qn); r.rev = (flag & 16) != 0;
        for (int i = 0; i < qn; i++) freq[(uint8_t)q[i]]++;
        totalQual += qn;
        recs.push_back(r);
    }
    fclose(f); free(line);

    int alpha = 0; for (int i = 0; i < 256; i++) if (freq[i]) alpha++;
    printf("记录数 %zu  质量值字节 %zu  字母表 %d\n", recs.size(), totalQual, alpha);
    if (recs.empty()) return 1;

    size_t cap = totalQual + (totalQual >> 1) + (1 << 20);
    std::vector<uint8_t> out(cap, 0);

    coder_io encIo(out.data(), (int32_t)cap);
    coder_fcv2 enc(&encIo, freq);
    for (size_t i = 0; i < recs.size(); i++) {
        enc.encode_record((const uint8_t*)recs[i].qual.data(),
                          (uint32_t)recs[i].qual.size(), recs[i].rev);
    }
    int32_t compLen = enc.encode_flush();
    printf("压缩后 %d 字节, 比率 %.4f%%\n", compLen, 100.0 * compLen / totalQual);
    if (compLen <= 0) { printf("往返: FAIL (压缩输出为空)\n"); return 1; }

    coder_io decIo(out.data(), compLen);
    coder_fcv2 dec(&decIo, freq);
    if (dec.begin_decode() != 0) { printf("往返: FAIL (begin_decode)\n"); return 1; }

    std::vector<uint8_t> buf(1 << 16);
    size_t bad = 0;
    for (size_t i = 0; i < recs.size(); i++) {
        uint32_t n = (uint32_t)recs[i].qual.size();
        dec.decode_record(buf.data(), n);
        if (memcmp(buf.data(), recs[i].qual.data(), n) != 0) {
            if (bad < 3) printf("  记录 %zu 不一致\n", i);
            bad++;
        }
    }
    printf("往返: %s (%zu/%zu 条不一致)\n", bad ? "FAIL" : "PASS", bad, recs.size());
    return bad ? 1 : 0;
}
