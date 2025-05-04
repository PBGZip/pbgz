#ifndef _CODER_IO_H_
#define _CODER_IO_H_

struct coder_io
{
    enum mode
    {
        MENC,
        MDEC,
        MUNSET
    };

    coder_io(const uint8_t *buff, int32_t buff_len)
    {
        data = (uint8_t *)buff;
        data_capacity = buff_len;
        data_len = 0;
        meta.clear();
        m = MUNSET;
    }

    /* 追加coder标识 */
    void appen_magic(const std::string magic)
    {
        meta["magic"] = magic;
    }

    std::string get_magic() const
    {
        return meta["magic"].asString();
    }

    /* 设置level */
    void set_level(int32_t level)
    {
        meta["level"] = level;
    }

    const int get_level() const
    {
        return meta["level"].asInt();
    }

    /* 写一个字符 */
    void putc(uint8_t c)
    {
        *(data + data_len++) = c;
    }

    /* 读一个字符 */
    uint8_t getc()
    {
        return (data_len == data_capacity) ? '\0' : (*(data + data_len++));
    }

    /* io模式 */
    mode m;

    uint8_t *data;
    /* data的总长度 */
    int32_t data_capacity;
    /* 当前已经处理的长度 */
    int32_t data_len;
    /* 编码器参数传入和编码器输出的原信息由meta交互 */
    Json::Value meta;
};

#endif