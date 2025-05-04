#include "analy_type.h"
#include "manager.h"
#include "io.h"

/* 获取type对应的名称 */
std::shared_ptr<uint8_t> get_typename(const blocktype &btype)
{
    uint32_t len = 1;
    std::shared_ptr<uint8_t> ptr_name;
    switch (btype)
    {
    case FASTQ_GEN2:
        len += strlen("fastq(gen2)");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "fastq(gen2)", len);
        break;
    case FASTQ_GEN2_GZIP:
        len += strlen("fastq.gz(gen2)");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "fastq.gz(gen2)", len);
        break;
    case FASTQ_GEN3:
        len += strlen("fastq(gen3)");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "fastq(gen3)", len);
        break;
    case FASTQ_GEN3_GZIP:
        len += strlen("fastq.gz(gen3)");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "fastq.gz(gen3)", len);
        break;
    case GZIP:
        len += strlen("gz");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "gz", len);
        break;
    case BINARY:
        len += strlen("binary");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "binary", len);
        break;
    case BINARY_GZIP:
        len += strlen("binary.gz");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "binary.gz", len);
        break;
    case TYPE_UNKNOW:
        len += strlen("unknow");
        safe_new(len, uint8_t, ptr_name);
        memcpy(ptr_name.get(), "unknow", len);
        break;
    default:
        check_exit(false, ERR_INTERNEL, "can not get typename");
        break;
    }
    return ptr_name;
}

/*  id: @开头，baselen == qualitylen, comment: +或者与id相同
 *  这里'ACTGN'没做检查，放到block中做，提升性能. 分析后将回车的位置记录下来
 *  成功时返回最大base的长度
 */
int64_t is_fastq(const uint8_t *buffer, const int32_t &len, std::vector<uint32_t> &npos, bool eof)
{
    int64_t i, j, line_cnt = 0;
    int64_t line_start = -1; /* 标识每行开始对应文件的偏移位置 */
    int64_t max_len_base = 0;
    int64_t id_pos = -1, id_len = -1;    /* 记录一条id信息 */
    int64_t base_pos = -1, base_len = 0; /* 记录一条base信息 */
    const std::string actg = "ACTGNactgn";
    bool flag = true;

    for(i = 0; i < len; i++) {
        if(*(buffer+ i) == '\n') {
            npos.push_back(i);
            const uint8_t *pline = buffer + line_start; /* 指向解析到的行数据 */
            switch ((++line_cnt) & 0x3)
            {
            case 1: /* judge line id */
                if (*pline != '@')
                    flag = false;
                id_pos = line_start;
                id_len = i - line_start;
                break;

            case 2: /* judge line base */
                base_pos = line_start;
                base_len = i - line_start; // 不包含'\n'
                max_len_base = std::max(max_len_base, base_len);
                if (actg.find_first_of(*pline) == std::string::npos)
                    return false;
                break;

            case 3: /* judge line comment */
                // if (*pline != '+') { /*  不为+时，判断是否与id行内容相同  */
                //     if (id_len != i - line_start)
                //         flag = false;
                //     else if (memcmp((uint8_t *)pline, buffer + id_pos, id_len))
                //         flag = false;
                // }
                break;

            case 0: /* judge line quality */
                if (base_len != i - line_start)
                    flag = false;
                break;

            default:
                break;
            }

            if (!flag)
                return 0;

            line_start = -1;
        } else {
            if(line_start == -1)
                line_start = i; 
        }
    }

    if (eof && len > 0 && npos.back() + 1 != len) /*  考虑文件末尾没有\n的情况*/
        npos.push_back(len - 1);
    
    if (npos.size() < 4) /* 不够4行不认为是一个fastq */
        max_len_base = 0;

    return max_len_base;
}

bool is_gzip(const uint8_t *buffer, int32_t size)
{
    return (buffer && size > 2 && (*(buffer) == 0x1F ) && (*(buffer + 1) == 0x8B));
}

filetype judge_file_type(const std::string &file_name) 
{
    filetype ftype = TYPE_UNKNOW;
    int len_judge[2] = {(1 << 20), (GENE3_MAX_BASE << 2)};
    for (int i = 0; i  < sizeof(len_judge)/sizeof(int); i++) {
        io read_io(file_name, mr, false);
        block_rough block(len_judge[i]);
        read_io.read_one_block(block);
        if (block.btype != TYPE_UNKNOW) {
            ftype = static_cast<filetype>(block.btype);
            break;
        }
    }
    return ftype;
}

/* 数据块类型是否是fastq类型 */
bool block_type_isfastq(const blocktype &btype)
{
    return (btype == FASTQ_GEN2 || btype == FASTQ_GEN2_GZIP ||
            btype == FASTQ_GEN3 || btype == FASTQ_GEN3_GZIP);
}