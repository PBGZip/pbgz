#ifndef _CODER_JSON_H_
#define _CODER_JSON_H_

#include <json/json.h>
#include <zstd.h>
#include "../err_code.h"

/* json编码器，采用zstd流式压缩算法 */
class coder_json
{
public:
    /* level [1, 7] */
    coder_json(int32_t level=6);
    virtual ~coder_json();

    /* 压缩json数据，返回压缩后的长度，如果为负数说明out空间不够，如返回-28,说明out_len长度还需再加28
     * 调该函数则假设out_len肯定是足够的，否则一旦out_len不足则会报错退出
     */
    virtual int64_t encoder(const Json::Value &in, uint8_t *out, const int64_t out_len);

    /* 压缩json数据，在内部申请空间，压缩后数据存在out中 */
    virtual void encoder(const Json::Value &in, std::string &out);

    /* 解压json格式数据 */
    virtual void decoder(const uint8_t *in, const int64_t in_len, Json::Value &out);

private:
    /* 对应的压缩级别 */
    int32_t cLevel;
};

#endif