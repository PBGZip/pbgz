#ifndef _CODER_JSON_H_
#define _CODER_JSON_H_

#include <json/json.h>
#include <zstd.h>

#include "coder.h"

/* JSON encoder using zstd streaming compression algorithm */
class coder_json
{
public:
    /* level [1, 7] */
    coder_json(int32_t level=6);
    virtual ~coder_json();

    /* Compress JSON data, return compressed length, negative value indicates insufficient out space, e.g., return -28 means out_len needs 28 more bytes
     * This function assumes out_len is sufficient, otherwise it will exit with error if out_len is insufficient
     */
    virtual int64_t encoder(const Json::Value &in, uint8_t *out, const int64_t out_len);

    /* Compress JSON data, allocate space internally, compressed data stored in out */
    virtual void encoder(const Json::Value &in, std::string &out);

    /* Decompress JSON format data */
    virtual void decoder(const uint8_t *in, const int64_t in_len, Json::Value &out);

private:
    /* Corresponding compression level */
    int32_t cLevel;
};

#endif