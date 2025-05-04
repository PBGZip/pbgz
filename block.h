#ifndef _BLOCK_H_
#define _BLOCK_H_

#include "memory.h"
#include "vendor.h"
#include <json/json.h>
#include "md5sum.h"

struct block_rough /*  粗糙的块*/
{
    block_rough(uint32_t len = BLOCK_SIZE) : buffer_size(len)
    {
        safe_alloc(len, uint8_t, this->buffer);
        reset();
    }
    virtual ~block_rough()
    {
        if (buffer)
            free(buffer);
    }

    void reset()
    {
        block_id = -1;
        current_len = 0;
        meta_encoded_len = 0;
        max_line_len = 0;
        btype = TYPE_UNKNOW;
        npos.clear();
    }

    const uint8_t* get_buffer() const {
        return buffer;
    }

    uint8_t * get_curr() {
        return buffer + current_len;
    }

    uint32_t get_remain() {
        return buffer_size - current_len;
    }

    uint8_t *buffer;
    uint32_t buffer_size;       /*  buffer的size */
    uint32_t current_len;       /* 当前buffer的实际长度，该长度为block数据压缩后长度 + block meta压缩后长度之和 */
    blocktype btype;            /*  当前块对应块类型 */
    std::vector<uint32_t> npos; /*  buffer中换行符的位置，从0开始 */
    int64_t block_id;
    uint32_t max_line_len;
    uint32_t meta_encoded_len; /* meta压缩后的数据会放在buffer一起，用该值标识meta压缩后流的长度 */
};

typedef block_rough *block_rough_ptr;

static inline void calc_md5(std::string &md5_get, const uint8_t *data, uint32_t len)
{
    MD5_CONTEXT md5;
    md5_init(&md5);
    md5_write(&md5, (uint8_t *)data, len);
    md5_final(&md5);
    md5_get = md5.hexstr();
}

#endif
