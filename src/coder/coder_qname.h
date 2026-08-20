/*
 * coder_qname.h - QNAME 专用编码器
 *
 * 背景：HG00106 这类由多个 run 拼接、coordinate-sorted 的 SAM，QNAME 形如
 *   ERR015528.21801860
 * 固定 9 字符前缀只取 5 种值（熵约 2.3 bit），而同片段 R1/R2 两条记录的 QNAME
 * 完全相同、坐标排序后往往相邻（行距 <64 占 99%+）。coder_affix_match 只做
 * 相邻行共享前缀匹配，对随机数字 id 基本无效；本编码器针对这两点：
 *
 *   1) 跨行去重：维护"最近行哈希表"，QNAME 曾出现过且行距在 RING_SIZE 以内时，
 *      只写 1 bit 标记 + 行距（行距熵实测 ~5.1 bit），不再重复写整串；
 *   2) 按位置建模：未命中的行按"行长 + 每位置一个自适应模型"逐字符写，前缀段
 *      收敛到极小的字符集，数字段收敛到接近均匀数字的熵。
 *
 * 哈希表设计（编码端专用，解码端只维护环形缓冲、不建表）：
 *   - 直接映射：槽 = hash & MASK，无线性探测、无删除，天然有界、不会死循环。
 *   - 槽存"最近出现行号"；命中需二次校验（行在窗口内 + 内容一致）防哈希碰撞误判。
 *   - 碰撞代价仅是丢一次 copy 机会（同槽被窗口内的其它行覆盖），dist 小（本文件
 *     中位 19）时命中率损失可忽略。
 *
 * 输入约定：整条 QNAME 含末尾 '\t'（与 compressIdFieldInAll 一致），in_len<=255。
 */
#ifndef _CODER_QNAME_H_
#define _CODER_QNAME_H_

#include <stdint.h>
#include <cstring>

#include "simple_model.h"
#include "coder_io.h"
#include "coder.h"

class coder_qname : public coder
{
public:
    static const uint32_t DM_BITS = 18;
    static const uint32_t DM_SIZE = 1u << DM_BITS;
    static const uint32_t RING_SIZE = 4096;
    /* 支持的最长 QNAME(含 '\t')：SAM 规范 QNAME<=254，+tab=255。 */
    static const uint32_t MAX_POS = 256;

    coder_qname(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qname");
        this->flushed = false;
        this->line_idx = 0;
        this->last_len = 0;

        void* ptr = safe_alloc(MAX_POS * sizeof(SIMPLE_MODEL<256>));
        model_pos = static_cast<SIMPLE_MODEL<256>*>(ptr);
        check_exit(model_pos != nullptr, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
            "coder_qname: model_pos alloc failed");
        for (uint32_t n = 0; n < MAX_POS; n++) {
            model_pos[n].reset();
        }
        for (uint32_t n = 0; n < 256; n++) {
            model_len[n].reset();
        }
        model_copy.reset();
        model_dist.reset();
        memset(dm, 0xff, sizeof(dm)); /* line=0xffffffff 表示空 */
        memset(rlen, 0, sizeof(rlen));
    }

    virtual ~coder_qname()
    {
        if (model_pos) {
            safe_free((void**)&model_pos);
        }
        if (!flushed && io->m == coder_io::MENC) {
            encode_flush();
        }
    }

    /* 编码一行 QNAME（含末尾 '\t'）。 */
    void encode_line(const uint8_t *in, const uint32_t in_len, [[maybe_unused]] bool need2hold = false)
    {
        check_exit(in_len <= MAX_POS, coder_ns::CODER_ERR_BAD_ARGS,
            "coder_qname: QNAME too long (%u)", in_len);

        if (io->m != coder_io::MENC) {
            rc.output((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartEncode();
            io->m = coder_io::MENC;
        }

        uint32_t h = hash(in, in_len);
        const uint32_t ridx = line_idx % RING_SIZE;
        const uint32_t slot = h & (DM_SIZE - 1);

        /* 探测最近出现：直接映射槽命中 + 行在窗口内 + 内容一致 -> COPY。 */
        int32_t dist = -1;
        if (dm[slot].hash == h && dm[slot].line != UINT32_MAX) {
            uint32_t lastLine = dm[slot].line;
            if (lastLine < line_idx && line_idx - lastLine <= RING_SIZE - 1 &&
                rline[lastLine % RING_SIZE] == lastLine &&
                rlen[lastLine % RING_SIZE] == in_len &&
                memcmp(recent[lastLine % RING_SIZE], in, in_len) == 0) {
                dist = (int32_t)(line_idx - lastLine);
            }
        }

        if (dist >= 0) {
            model_copy.encodeSymbol(&rc, 1);
            model_dist.encodeSymbol(&rc, (uint16_t)dist);
        } else {
            model_copy.encodeSymbol(&rc, 0);
            model_len[last_len].encodeSymbolOrder(&rc, in_len);
            for (uint32_t p = 0; p < in_len; p++) {
                model_pos[p].encodeSymbol(&rc, in[p]);
            }
        }

        last_len = (int32_t)in_len;

        /* 写入环形缓冲与直接映射表。 */
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)in_len;
        memcpy(recent[ridx], in, in_len);

        dm[slot].hash = h;
        dm[slot].line = line_idx;

        line_idx++;
    }

    /* 解码一行 QNAME（含末尾 '\t'），返回写入 out 的字节数。 */
    int32_t decode_line(uint8_t *out, uint32_t out_len,
                        [[maybe_unused]] uint8_t split_ch = UINT8_MAX,
                        [[maybe_unused]] bool need2hold = false)
    {
        if (io->m != coder_io::MDEC) {
            rc.input((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartDecode();
            io->m = coder_io::MDEC;
        }

        uint16_t copy = model_copy.decodeSymbol(&rc);
        uint32_t len;
        if (copy) {
            uint32_t dist = model_dist.decodeSymbol(&rc);
            /* dist 可为 line_idx（复制第 0 行），但必须 >= 1 且 <= line_idx。 */
            if (dist == 0 || dist > line_idx) {
                return coder_ns::CODER_ERR_BAD_ARGS;
            }
            uint32_t ridx = (line_idx - dist) % RING_SIZE;
            len = rlen[ridx];
            if (len > out_len) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            memcpy(out, recent[ridx], len);
        } else {
            len = model_len[last_len].decodeSymbolOrder(&rc);
            if (len > out_len || len > MAX_POS) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            for (uint32_t p = 0; p < len; p++) {
                out[p] = (uint8_t)model_pos[p].decodeSymbol(&rc);
            }
        }
        last_len = (int32_t)len;

        /* 与编码端一致地更新环形缓冲（无哈希表：解码不需要探测）。 */
        uint32_t ridx = line_idx % RING_SIZE;
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)len;
        memcpy(recent[ridx], out, len);

        line_idx++;
        return (int32_t)len;
    }

    /* fake */
    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) { return 0; }

    void encode_flush()
    {
        if (flushed) {
            return;
        }
        if (io->m == coder_io::MENC) {
            if (rc.FinishEncode() < 0) {
                io->set_err(coder_io::IO_BUF_FULL);
            }
            io->data_len += rc.size_out();
        }
        flushed = true;
    }

private:
    static uint32_t hash(const uint8_t *s, uint32_t n)
    {
        uint32_t h = 2166136261u;
        for (uint32_t i = 0; i < n; i++) {
            h ^= s[i];
            h *= 16777619u;
        }
        return h;
    }

    RangeCoder rc;
    SIMPLE_MODEL<2> model_copy;
    SIMPLE_MODEL<4096> model_dist;
    SIMPLE_MODEL<256> model_len[256];
    SIMPLE_MODEL<256> *model_pos;

    struct HashSlot {
        uint32_t hash;
        uint32_t line;
    } dm[DM_SIZE];
    uint8_t rlen[RING_SIZE];
    uint32_t rline[RING_SIZE];
    uint8_t recent[RING_SIZE][MAX_POS];

    uint32_t line_idx;
    int32_t last_len;
    bool flushed;
};

#endif
