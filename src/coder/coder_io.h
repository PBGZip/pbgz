#ifndef _CODER_IO_H_
#define _CODER_IO_H_

#include <stdint.h>
#include <json/json.h>

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

    /* Append coder identifier */
    void appen_magic(const std::string magic)
    {
        meta["magic"] = magic;
    }

    std::string get_magic() const
    {
        return meta["magic"].asString();
    }

    /* Set level */
    void set_level(int32_t level)
    {
        meta["level"] = level;
    }

    int get_level() const
    {
        return meta["level"].asInt();
    }

    /* Write one character */
    void putc(uint8_t c)
    {
        *(data + data_len++) = c;
    }

    /* Read one character */
    uint8_t getc()
    {
        return (data_len == data_capacity) ? '\0' : (*(data + data_len++));
    }

    /* IO mode */
    mode m;

    uint8_t *data;
    /* Total length of data */
    int32_t data_capacity;
    /* Currently processed length */
    int32_t data_len;
    /* Encoder parameter input and output metadata interaction through meta */
    Json::Value meta;
};

#endif