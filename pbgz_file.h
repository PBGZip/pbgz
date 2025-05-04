#ifndef _PBGZ_FILE_H_
#define _PBGZ_FILE_H_

#include "manager.h"
#include "io.h"

/* 读写pbgz格式文件 */

typedef enum pbgz_stream_type
{
    /* 文件meta流，如整个文件有多少条流，多少个块等，通常为json格式 */
    NST_META_FILE,

    /* 块的meta流，针对单个块，该流后面通过接对应的block */
    NST_META_BLOCK,

    /* fastq类型的数据块流 */
    NST_FASTQ_DATA_BLOCK,

    /* 未归类的数据块流 */
    NST_OTHER_DATA_BLOCK,

    /* reference的meta流，针对单个块，该流后面通过接对应的block */
    NST_META_REFE,

    /* reference数据流，譬如fastq的块压缩对应的流 */
    NST_DATA_REFE,

    /* 只用一个字符表示 */
    NST_UNKNOW = UINT8_MAX
} nstype;

typedef struct pbgz_stream_header
{
public:
    /* 分配一条pbgz stream格式的头 */
    pbgz_stream_header()
    {
        /* SID + 4个bytes的mainid标识(hash值)  + 8个字节的subid  + 1个byte的流类型(支持256种)  + 4个bytes的流长度(支持4G) */
        std::string sid = STREAM_ID;
        streamid_len = sid.length();
        buffer_len = streamid_len + mainid_len + subid_len + type_len + data_len;
        safe_alloc(buffer_len, uint8_t, buffer);
        memcpy(buffer, sid.c_str(), sid.length());
        this->alloced = true;
    }

    /* 将外部buffer转换为pbgz stream格式头解析 */
    pbgz_stream_header(uint8_t *header)
    {
        /* SID + 4个bytes的mainid标识(hash值)  + 8个字节的subid  + 1个byte的流类型(支持256种)  + 4个bytes的流长度(支持4G) */
        std::string sid = STREAM_ID;
        streamid_len = sid.length();
        buffer_len = streamid_len + mainid_len + subid_len + type_len + data_len;
        this->buffer = header;
        memcpy(buffer, sid.c_str(), sid.length());
        this->alloced = false;
    }

    virtual ~pbgz_stream_header()
    {
        if (alloced && buffer)
            free(buffer);
    }

    void set(uint32_t main_id, int64_t sub_id, uint8_t type, uint32_t data_len)
    {
        uint8_t *pmainid, *psubid, *ptype, *pdatalen;
        pmainid = buffer + streamid_len;
        *((uint32_t *)pmainid) = main_id;
        psubid = (uint8_t *)(pmainid) + mainid_len;
        *((int64_t *)psubid) = sub_id;
        ptype = (uint8_t *)(psubid) + subid_len;
        *((uint8_t *)ptype) = type;
        pdatalen = (uint8_t *)(ptype) + type_len;
        *((uint32_t *)pdatalen) = data_len;
    }

    /* 检查头标识是否ok */
    bool check_mark() const{
        if (get_bufferlen() < strlen(STREAM_ID))
            return false;
        return (memcmp(get_buffer(), STREAM_ID, strlen(STREAM_ID)) == 0);
    }

    uint8_t *get_buffer() const
    {
        return this->buffer;
    }

    const uint32_t get_bufferlen() const
    {
        return this->buffer_len;
    }

    const uint32_t get_mainid() const
    {
        return *((uint32_t *)(buffer + streamid_len));
    }

    const int64_t get_subid() const
    {
        return *((int64_t *)(buffer + streamid_len + mainid_len));
    }

    const uint8_t get_type() const
    {
        return *((uint8_t *)(buffer + streamid_len + mainid_len + subid_len));
    }

    const uint32_t get_datalen() const
    {
        return *((uint32_t *)(buffer + streamid_len + mainid_len + subid_len + type_len));
    }

    // /* 流的标识，如果后面流格式有修改，可以把改为其他字样  */
    // const std::string id = STREAM_ID;
    // /* 流的主id，一个pbgz文件所有流main_id相同，由文件名和创建时间组成生成hash值，尽量随机种子保证唯一性 */
    // uint32_t main_id;
    // /* 流的子id，在当前流中的偏移，从1开始 */
    // int64_t sub_id;
    // /* 流的类型，对应pbgz_stream_type */
    // uint8_t type;
    // /* 流对应的长度 */
    // uint32_t data_len;

private:
    uint8_t *buffer;
    uint32_t buffer_len;
    uint32_t streamid_len;
    bool alloced;

    const uint32_t mainid_len = sizeof(uint32_t);
    const uint32_t subid_len = sizeof(int64_t);
    const uint32_t type_len = sizeof(uint8_t);
    const uint32_t data_len = sizeof(uint32_t);
} nshead;

class pbgz_file
{
public:
    pbgz_file(const std::string file_name, bool readonly);
    virtual ~pbgz_file();

    /* 从pbgz文件中读一对流：data meta + data，如果读出错了内部会处理，返回读到的流的长度 */
    virtual int64_t read_pair_stream(block_rough &block);

    /* 写一条流，如果写出错内部会处理 */
    virtual void write_one_stream(const uint8_t *buffer, uint32_t buffer_len, nstype stream_type, bool eof = false);

    /* 刷磁盘 */
    virtual void flush();

    /* 获取文件当前的offset */
    virtual const int64_t get_offset() const;

    /* 获取当前的sub id */
    virtual const int64_t get_subid() const;

    /* 设置主备两条meta流的文件偏移地址 */
    virtual void set_meta_offset(int64_t offset_1, int64_t offset_2);

    /* 解析全局的文件meta流 */
    virtual void get_file_meta(Json::Value &file_meta);

    /* 得到所有流的offset */
    virtual const Json::Value get_stream_offset();

private:
    /* 流的标识，如果后面流格式有修改，可以把改为其他字样  */
    const std::string id_stream = STREAM_ID;

    /* 文件标识 */
    const std::string id_file = FILE_ID;

    std::string file;
    bool readonly;

    /* 写pbgz文件时对应的当前的main id */
    uint32_t mid2write;

    /* 写pbgz文件时对应的当前sub id */
    int64_t sid2write;

    /* 读写文件的io*/
    io *fileio;

    /* 当前文件指针的offset */
    int64_t file_offset;

    /* 第一条file meta流的offset */
    int64_t file_meta_offset;

    /* 读pbgz文件时当前读到的流的计数 */
    int64_t streams_readed;

    /* 流开始的位置 */
    int64_t stream_start;

    /* 流的header buffer */
    nshead *head;

    /* 所有流的offset */
    Json::Value meta_stream_offset;
};

#endif