#ifndef _ACTUATOR_FASTQ_H_
#define _ACTUATOR_FASTQ_H_

#include "block.h"
#include "reference.h"
#include "coder/coder.h"
#include <json/json.h>
#include "coder/coder_qual.h"

/* fastq传动装置 */

// #define ANALYZ_STREAMS /* 分析流的压缩率 */

class actuator_fastq
{
private:
    enum comtype
    {
        CT_JUST_PLUS,  /* comment只有一个加号 */
        CT_SAME_AS_ID, /* comment与id内容相同 */
        CT_OTHER,      /* comment非上面两类 */
        CT_UNKNOW
    };

    typedef struct mapping
    {
        void set(uint8_t *_ps, uint32_t _ps_len, uint8_t *_psp, uint32_t _psp_len, uint32_t _offset)
        {
            ps[0] = _ps;
            ps_len[0] = _ps_len;
            ps[1] = _psp;
            ps_len[1] = _psp_len;
            offset = _offset;
        }
        inline void inc_offset() { ++offset; }
        uint8_t *get_squash(bool pair) const {return ((pair) ? (ps[1] - offset) : (ps[0] + offset));}

        /*  squash buffer , 0表示正向，1表示互补链 */
        uint8_t *ps[2];      /* 指向squash buffer */
        uint32_t ps_len[2];  /* 对应squash buffer的数据长度 */
        uint32_t ps_l_unalign_len[2]; /* 左边没有对齐的碱基长度 */
        uint8_t ps_l_unalign[2][3]; /* 左边没有对齐的碱基，这里存的是squash之后的数据 */
        uint32_t ps_r_unalign_len[2]; /* 右边没有对齐的碱基长度 */
        uint8_t ps_r_unalign[2][3]; /* 右边没有对齐的碱基，这里存的是squash之后的数据 */
        uint32_t offset;  /* 当前偏移 */
    } mapping_t;

public:
    actuator_fastq(const block_rough_ptr bptr_in, block_rough_ptr bptr_out, const reference *r);
    virtual ~actuator_fastq();

    /* 详细分析fastq，格式不符合时返回false */
    bool analyze_fastq();

    /* 压缩fastq */
    bool compress();

    /* 解压fastq */
    bool decompress();

private:
    /* 压缩id行 */
    bool compress_id();

    /* 解压id行初始化 */
    bool initialize_decode_id(std::vector<coder *> &id_decoders);

    /* 解压base行初始化 */
    bool initialize_decode_base(coder *&base_decoder);

    /* 解压comment行初始化 */
    bool initialize_decode_comment(coder *&com_decoder);

    /* 解压quality行初始化 */
    bool initialize_decode_quality(coder_qual *&qual_decoder);

    /* 压缩base行 */
    bool compress_base();

    /* 压缩初始化 */
    void compress_initialize();

    /* 压缩comment */
    bool compress_comment();

    /* 压缩质量行 */
    bool compress_quality();

    /* 解压一行base */
    uint32_t decode_base_line(uint8_t *out, uint32_t out_len, coder *base_decoder);

    /* mapping */
    void (actuator_fastq::*mapping)(const uint8_t *, uint32_t, uint8_t *&, uint32_t &, uint64_t &, uint8_t &);

    /* 二代数据匹配reference */
    inline void mapping_gen2(const uint8_t *base, uint32_t base_len, uint8_t *&out, uint32_t &out_len, uint64_t &mpos, uint8_t &mdir);

    /* 三代数据匹配reference */
    inline void mapping_gen3(const uint8_t *base, uint32_t base_len, uint8_t *&out, uint32_t &out_len, uint64_t &mpos, uint8_t &mdir);

private:
    block_rough_ptr indata;
    block_rough_ptr outdata;

    /* 参考基因组 */
    reference *refgene;

    /* 定义一个fastq块 */

    /* id的分割符 */
    const std::string idsplit = "/:= _.,-#\r\n\t";
    /* id第一行解析后分割符信息 */
    uint8_t *idsplit_syms;
    uint32_t idsplit_symslen;
    uint32_t *idsplit_minlen;
    uint32_t *idsplit_maxlen;

    /* id按照分割符切割后所有片段的位置信息信息，其中位置为相对该行开头的偏移 */
    uint16_t *idpos;
    uint32_t idpos_offset;

    /* 是否是二代数据 */
    bool is_gen2;

    /* comment行 */
    comtype ctype;
    
    /* 当前块base最小长度和最大长度 */
    uint32_t baselen_min;
    uint32_t baselen_max;

    /* 当前压缩的fastq块中base中N的个数 */
    uint32_t base_n_cnt;

    /* 当前block对应编码器的meta信息 */
    Json::Value meta;

    /* 质量数符号频率表 */
    std::vector<std::pair<uint16_t, uint16_t>> qual_freq_table;

private:
    /* some buffer for base mapping or encode */
    uint8_t *somebuffer;

    /* base的互补链 */
    uint8_t *somebuffer_basepair;

    /* base的互补链squash, 4->1byte */
    uint8_t *somebuffer_basepair_squash[4];

    /* base squash, 4->1byte */
    uint8_t *somebuffer_base_squash[4];

    /* base mapping完后比对结果的buffer */
    uint8_t *somebuffer_base_mapped;
    uint32_t len_base_mapped;

    /* base mapping之后的位置 */
    uint64_t *somebuffer_base_mpos;

    /* base 是否mapping到互补链 */
    uint8_t *somebuffer_base_mpair;

    /* base每行去除N后的存储buffer */
    uint8_t *somebuffer_base_stripn;

    /* base中N的位置 */
    uint32_t *somebuffer_base_npos;

    /* base每行的长度 */
    uint16_t *somebuffer_base_len2;
    uint32_t *somebuffer_base_len3;

    /* reference strech之后的buffer, 解压时使用 */
    uint8_t *somebuffer_refe_stretch;

    /* debug */
#ifdef ANALYZ_STREAMS
    int64_t mapping_cnt = 0;
    int64_t actg_total = 0;
    int64_t actg_matched = 0;
#endif
};

#endif