#ifndef _REFERENCE_H_
#define _REFERENCE_H_

/*  参考基因检索 */

#include <stdint.h>
#include <json/json.h>
#include "manager.h"
#include "bar.h"

#if __x86_64__  || __i386__ || _M_X64 || _M_IX86
#define XCHGBACKOFF // 0.9s
#else
#define PTHREAD
#endif

#ifdef XCHG
#include "spinlock-xchg.h"
#elif defined(XCHGBACKOFF)
#include "spinlock/spinlock-xchg-backoff.h"
#elif defined(K42)
#include "spinlock/spinlock-k42.h"
#elif defined(MCS)
#include "spinlock/spinlock-mcs.h"
#elif defined(TICKET)
#include "spinlock/spinlock-ticket.h"
#elif defined(PTHREAD)
#include "spinlock/spinlock-pthread.h"
#elif defined(CMPXCHG)
#include "spinlock/spinlock-cmpxchg.h"
#elif defined(RTM)
#include "spinlock/spinlock-xchg.h"
#include "spinlock/rtm.h"
#elif defined(HLE)
#include "spinlock/spinlock-xchg-hle.h"
#else
#error "must define a spinlock implementation"
#endif

typedef std::vector<std::pair<std::pair<uint32_t *, uint32_t>, uint32_t>> htable;

class reference
{
private:
    typedef struct
    {
        uint32_t basegroup_pos; /* 碱基在referene中的偏移的basegroup_step个数 */
        uint32_t hash_bucket;   /* 碱基对应到hash的buket */
    } bg_hash;

public:
    reference(const std::string &fasta);
    virtual ~reference();

    /* 对外接口，创建索引表 */
    bool make_index();

    /* 查询位置 */
    const uint32_t* query_pos(const uint32_t &hash, uint32_t &len);

    /* 获取建索引时取的碱基的长度 */
    const uint32_t get_bglen() const;

    /* 获取建索引时取的碱基的步长 */
    const uint32_t get_bgstep() const;

    /* 得到reference squash buffer */
    const uint8_t *get_squash() const;

    /* 得到reference squash buffer的长度 */
    const int64_t get_squashlen() const;

    /* 获取fasta文件名 */
    const std::string get_fasta_name() const;

    /* 获取fasta文件内容长度 */
    const int64_t get_fasta_len() const;

    /* 获取fasta文件内容md5 */
    const std::string get_fasta_md5() const;

    /* 获取ni文件路径 */
    const std::string get_ni_path() const;

    /* 将没有matched上的reference清零 */
    void reference_squash_sanitize(int64_t start_pos, int64_t len);

    /* 打包参考基因组有效信息 */
    void pack();

    /* 更新参考基因组匹配信息，以一个squash字节为单位更新，这样可以并发更新 */
    void update_matched(uint64_t actg_pos, uint32_t match_len);

    /* 得到reference对应的ni文件名 */
    void ni_from_reference(std::string &ni);

    /* 初始化参考基因squash buffer, 通过.ni索引文件 */
    bool initialize_squash_1();

    /* 初始化参考基因squash buffer, 通过stream */
    uint8_t *initialize_squash_2(int64_t squash_len);

    /* 获取指定位置对应长度的actg碱基 */
    void get_stretch_actg(uint8_t *out, uint32_t out_len, uint64_t actg_pos);

    /* 获取指定位置对应长度的squash碱基，即2个bits 放到一个字符的末尾*/
    void get_stretch_2bits1char(uint8_t *out, uint32_t out_len, uint64_t actg_pos);

    /*  根据每个字节末尾的2个bits，转换成actg */
    void get_actg_from2bits(const uint8_t *src_2bits, uint32_t src_2blits_len, uint8_t *dst_actg);

private:
    /* 制作索引表：从参考基因组中拿碱基*/
    void make_step1_fetchBG(bg_hash *&hash);

    /* 制作索引表：计算hash table的size */
    void make_step2_calcHT(htable &hash_table);

    /* 制作索引表：初始化hash table */
    void make_step3_initHT(const htable &hash_table, uint32_t *&hash_bucket_curpos);

    /* 制作索引表：将碱基的位置信息存入hash table*/
    void make_step4_buildHT(const bg_hash *hash, uint32_t *&hash_bucket_curpos);

    /* 制作索引表：将同一个bucket中的位置排序，保证并发生成的索引表都一致 */
    void make_step5_sortHT();

    /* .ni文件是否有效 */
    bool ni_valid(const std::string &ni);

    /* 创建.ni索引文件 */
    bool ni_make(const std::string &ni);

    /* 得到reference对应的ni文件名 */
    void ni2reference(std::string &ni);

    /* dump hash table to file */
    void dump_hash_table();

private:
    /* 进度条*/
    guard_bar *mbar;
    int64_t mbar_current, mbar_total;

    /* fasta文件内容的md5 */
    std::string fasta_md5;

    /* ni 文件的路径 */
    std::string path_ni;

    /* fasta文件文件内容的长度 */
    int64_t fasta_len;

    /* 建立索引表时取碱基作为key的长度，必须为奇数 */
    const uint32_t basegroup_len = 31;

    /* 建立索引表时碱基的步长 */
    const uint32_t basegroup_step = 32;

    /* 当前的并发 */
    uint32_t parallel;

    /*  hash table信息 */
    /* 每个hash buckets元素个数 */
    uint32_t *hash_bucket_cnt;

    /*  记录参考基因组碱基位置信息的hash table buffer*/
    uint32_t *hash_table_buffer;

    /* 参考基因组文件 */
    std::string reference_file;

    /* 参考基因actg编码后的buffer */
    uint8_t *reference_squash;

    /* 参考基因actg编码后的buffer 对应的长度 */
    int64_t reference_squashlen;

    /* 参考基因组actg编码后,匹配上了base的buffer */
    uint8_t *reference_squash_matched;

    /* 对应的长度 */
    uint64_t reference_squash_matchedlen;

    /* 缓存文件 */
    Json::Value ref2ni_cache;

    /* actg hash信息 */
    const int32_t hash_buckets = (32 << 20); /* 32M */

    /* hash buckets对应的mask */
    const int32_t hash_mask = hash_buckets - 1;

    /* 1 byte squashed actg stretch to 4 actg */
    const uint32_t actg_stretch[256] = {
        0x41414141, 0x43414141, 0x54414141, 0x47414141, 0x41434141, 0x43434141, 0x54434141, 0x47434141, 0x41544141, 0x43544141, 0x54544141, 0x47544141, 0x41474141, 0x43474141, 0x54474141, 0x47474141,
        0x41414341, 0x43414341, 0x54414341, 0x47414341, 0x41434341, 0x43434341, 0x54434341, 0x47434341, 0x41544341, 0x43544341, 0x54544341, 0x47544341, 0x41474341, 0x43474341, 0x54474341, 0x47474341,
        0x41415441, 0x43415441, 0x54415441, 0x47415441, 0x41435441, 0x43435441, 0x54435441, 0x47435441, 0x41545441, 0x43545441, 0x54545441, 0x47545441, 0x41475441, 0x43475441, 0x54475441, 0x47475441,
        0x41414741, 0x43414741, 0x54414741, 0x47414741, 0x41434741, 0x43434741, 0x54434741, 0x47434741, 0x41544741, 0x43544741, 0x54544741, 0x47544741, 0x41474741, 0x43474741, 0x54474741, 0x47474741,
        0x41414143, 0x43414143, 0x54414143, 0x47414143, 0x41434143, 0x43434143, 0x54434143, 0x47434143, 0x41544143, 0x43544143, 0x54544143, 0x47544143, 0x41474143, 0x43474143, 0x54474143, 0x47474143,
        0x41414343, 0x43414343, 0x54414343, 0x47414343, 0x41434343, 0x43434343, 0x54434343, 0x47434343, 0x41544343, 0x43544343, 0x54544343, 0x47544343, 0x41474343, 0x43474343, 0x54474343, 0x47474343,
        0x41415443, 0x43415443, 0x54415443, 0x47415443, 0x41435443, 0x43435443, 0x54435443, 0x47435443, 0x41545443, 0x43545443, 0x54545443, 0x47545443, 0x41475443, 0x43475443, 0x54475443, 0x47475443,
        0x41414743, 0x43414743, 0x54414743, 0x47414743, 0x41434743, 0x43434743, 0x54434743, 0x47434743, 0x41544743, 0x43544743, 0x54544743, 0x47544743, 0x41474743, 0x43474743, 0x54474743, 0x47474743,
        0x41414154, 0x43414154, 0x54414154, 0x47414154, 0x41434154, 0x43434154, 0x54434154, 0x47434154, 0x41544154, 0x43544154, 0x54544154, 0x47544154, 0x41474154, 0x43474154, 0x54474154, 0x47474154,
        0x41414354, 0x43414354, 0x54414354, 0x47414354, 0x41434354, 0x43434354, 0x54434354, 0x47434354, 0x41544354, 0x43544354, 0x54544354, 0x47544354, 0x41474354, 0x43474354, 0x54474354, 0x47474354,
        0x41415454, 0x43415454, 0x54415454, 0x47415454, 0x41435454, 0x43435454, 0x54435454, 0x47435454, 0x41545454, 0x43545454, 0x54545454, 0x47545454, 0x41475454, 0x43475454, 0x54475454, 0x47475454,
        0x41414754, 0x43414754, 0x54414754, 0x47414754, 0x41434754, 0x43434754, 0x54434754, 0x47434754, 0x41544754, 0x43544754, 0x54544754, 0x47544754, 0x41474754, 0x43474754, 0x54474754, 0x47474754,
        0x41414147, 0x43414147, 0x54414147, 0x47414147, 0x41434147, 0x43434147, 0x54434147, 0x47434147, 0x41544147, 0x43544147, 0x54544147, 0x47544147, 0x41474147, 0x43474147, 0x54474147, 0x47474147,
        0x41414347, 0x43414347, 0x54414347, 0x47414347, 0x41434347, 0x43434347, 0x54434347, 0x47434347, 0x41544347, 0x43544347, 0x54544347, 0x47544347, 0x41474347, 0x43474347, 0x54474347, 0x47474347,
        0x41415447, 0x43415447, 0x54415447, 0x47415447, 0x41435447, 0x43435447, 0x54435447, 0x47435447, 0x41545447, 0x43545447, 0x54545447, 0x47545447, 0x41475447, 0x43475447, 0x54475447, 0x47475447,
        0x41414747, 0x43414747, 0x54414747, 0x47414747, 0x41434747, 0x43434747, 0x54434747, 0x47434747, 0x41544747, 0x43544747, 0x54544747, 0x47544747, 0x41474747, 0x43474747, 0x54474747, 0x47474747};

/* 1 byte squashed actg stretch to 4 2bits */
const uint32_t actg_stretch_2bits[256] = {
        0x00000000, 0x01000000, 0x02000000, 0x03000000, 0x00010000, 0x01010000, 0x02010000, 0x03010000, 0x00020000, 0x01020000, 0x02020000, 0x03020000, 0x00030000, 0x01030000, 0x02030000, 0x03030000,
        0x00000100, 0x01000100, 0x02000100, 0x03000100, 0x00010100, 0x01010100, 0x02010100, 0x03010100, 0x00020100, 0x01020100, 0x02020100, 0x03020100, 0x00030100, 0x01030100, 0x02030100, 0x03030100,
        0x00000200, 0x01000200, 0x02000200, 0x03000200, 0x00010200, 0x01010200, 0x02010200, 0x03010200, 0x00020200, 0x01020200, 0x02020200, 0x03020200, 0x00030200, 0x01030200, 0x02030200, 0x03030200,
        0x00000300, 0x01000300, 0x02000300, 0x03000300, 0x00010300, 0x01010300, 0x02010300, 0x03010300, 0x00020300, 0x01020300, 0x02020300, 0x03020300, 0x00030300, 0x01030300, 0x02030300, 0x03030300,
        0x00000001, 0x01000001, 0x02000001, 0x03000001, 0x00010001, 0x01010001, 0x02010001, 0x03010001, 0x00020001, 0x01020001, 0x02020001, 0x03020001, 0x00030001, 0x01030001, 0x02030001, 0x03030001,
        0x00000101, 0x01000101, 0x02000101, 0x03000101, 0x00010101, 0x01010101, 0x02010101, 0x03010101, 0x00020101, 0x01020101, 0x02020101, 0x03020101, 0x00030101, 0x01030101, 0x02030101, 0x03030101,
        0x00000201, 0x01000201, 0x02000201, 0x03000201, 0x00010201, 0x01010201, 0x02010201, 0x03010201, 0x00020201, 0x01020201, 0x02020201, 0x03020201, 0x00030201, 0x01030201, 0x02030201, 0x03030201,
        0x00000301, 0x01000301, 0x02000301, 0x03000301, 0x00010301, 0x01010301, 0x02010301, 0x03010301, 0x00020301, 0x01020301, 0x02020301, 0x03020301, 0x00030301, 0x01030301, 0x02030301, 0x03030301,
        0x00000002, 0x01000002, 0x02000002, 0x03000002, 0x00010002, 0x01010002, 0x02010002, 0x03010002, 0x00020002, 0x01020002, 0x02020002, 0x03020002, 0x00030002, 0x01030002, 0x02030002, 0x03030002,
        0x00000102, 0x01000102, 0x02000102, 0x03000102, 0x00010102, 0x01010102, 0x02010102, 0x03010102, 0x00020102, 0x01020102, 0x02020102, 0x03020102, 0x00030102, 0x01030102, 0x02030102, 0x03030102,
        0x00000202, 0x01000202, 0x02000202, 0x03000202, 0x00010202, 0x01010202, 0x02010202, 0x03010202, 0x00020202, 0x01020202, 0x02020202, 0x03020202, 0x00030202, 0x01030202, 0x02030202, 0x03030202,
        0x00000302, 0x01000302, 0x02000302, 0x03000302, 0x00010302, 0x01010302, 0x02010302, 0x03010302, 0x00020302, 0x01020302, 0x02020302, 0x03020302, 0x00030302, 0x01030302, 0x02030302, 0x03030302,
        0x00000003, 0x01000003, 0x02000003, 0x03000003, 0x00010003, 0x01010003, 0x02010003, 0x03010003, 0x00020003, 0x01020003, 0x02020003, 0x03020003, 0x00030003, 0x01030003, 0x02030003, 0x03030003,
        0x00000103, 0x01000103, 0x02000103, 0x03000103, 0x00010103, 0x01010103, 0x02010103, 0x03010103, 0x00020103, 0x01020103, 0x02020103, 0x03020103, 0x00030103, 0x01030103, 0x02030103, 0x03030103,
        0x00000203, 0x01000203, 0x02000203, 0x03000203, 0x00010203, 0x01010203, 0x02010203, 0x03010203, 0x00020203, 0x01020203, 0x02020203, 0x03020203, 0x00030203, 0x01030203, 0x02030203, 0x03030203,
        0x00000303, 0x01000303, 0x02000303, 0x03000303, 0x00010303, 0x01010303, 0x02010303, 0x03010303, 0x00020303, 0x01020303, 0x02020303, 0x03020303, 0x00030303, 0x01030303, 0x02030303, 0x03030303};
};


#endif