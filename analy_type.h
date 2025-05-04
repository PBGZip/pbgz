#ifndef _ANALY_TYPE_H_
#define _ANALY_TYPE_H_

/* 数据类型分析识别 */

#include <stdint.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <memory>

typedef enum
{
    TYPE_UNKNOW = 0,
    GZIP = (1 << 0),
    BINARY = (1 << 1),
    FASTQ_GEN2 = (1 << 2),
    FASTQ_GEN3 = (1 << 3),
    BINARY_GZIP = (BINARY | GZIP),
    FASTQ_GEN2_GZIP = (FASTQ_GEN2 | GZIP),
    FASTQ_GEN3_GZIP = (FASTQ_GEN3 | GZIP),
    BAM = (1 << 4),
    PBGZFILE = UINT32_MAX
} filetype, blocktype;

/* 数据块类型是否是fastq类型 */
bool block_type_isfastq(const blocktype &btype);

/* 获取type对应的名称 */
std::shared_ptr<uint8_t> get_typename(const blocktype &btype);

/*  id: @开头，baselen == qualitylen, comment: +或者与id相同
 *  这里'ACTGN'没做检查，放到block中做，提升性能. 分析后将回车的位置记录下来
 *  成功时返回最大base的长度
 */
int64_t is_fastq(const uint8_t *buffer, const int32_t &len, std::vector<uint32_t> &npos, bool eof=false);

bool is_gzip(const uint8_t *buffer, int32_t size);

filetype judge_file_type(const std::string &file_name);

#endif