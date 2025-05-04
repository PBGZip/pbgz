#include "actuator_everything.h"
#include "manager.h"
#include "coder/coder_fc.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_json.h"
#include "pbgz_file.h"

actuator_everything::actuator_everything(const block_rough_ptr bptr_in, block_rough_ptr bptr_out)
{
    indata = bptr_in;
    outdata = bptr_out;
    this->outdata->block_id = this->indata->block_id;
    this->outdata->btype = this->indata->btype;
}

actuator_everything::~actuator_everything()
{
}

bool actuator_everything::compress()
{
    std::string md5;
    uint32_t src_len;
    int64_t len;
    Json::Value meta, sub_meta;
    info conf = manage::instance().get_zipinfo();

    coder_io block_io(outdata->get_curr(), outdata->get_remain());
    // block_io.set_level(conf.get_clevel());
    src_len = indata->current_len;

    // coder_ppmd block_cm(&block_io);
    // block_cm.encode(indata->buffer, indata->current_len);
    // block_cm.encode_flush();

    if (src_len <= FC_MIN_LEN || src_len >= FC_MAX_LEN) // 用bcm压缩
    {
        coder_bwt_cm coder(&block_io);
        coder.encode_line(indata->buffer, indata->current_len);
        coder.encode_flush();
    } else { // 用fc压缩
        coder_fc coder(&block_io);
        coder.encode_line(indata->buffer, indata->current_len);
        coder.encode_flush();
    }
    
    outdata->current_len += block_io.data_len;
    // fprintf(stderr, "\ndone in line %d, outdata->current_len %d....\n", __LINE__, outdata->current_len);
    sub_meta["srclen"] = (Json::Value::Int)src_len;
    sub_meta["dstlen"] = (Json::Value::Int)(block_io.data_len);
    sub_meta["coder"] = block_io.meta;

    /* 计算数据块的md5 */
    calc_md5(md5, indata->get_buffer(), indata->current_len);
    meta["stream"] = sub_meta;
    meta["md5"] = md5;

    /*  压缩block meta */
    coder_json cmeta;
    len = cmeta.encoder(meta, outdata->get_curr(), outdata->get_remain());
    check(len > 0, false, "block %ld meta info compress failed, get len: %ld", outdata->block_id, len);
    // std::cout << meta.dump() << std::endl;
    /*  设置meta压缩后的长度 */
    outdata->meta_encoded_len = len;
    outdata->current_len += len;

    if (outdata->current_len >= indata->current_len)
        check((indata->current_len * 2) <= outdata->buffer_size, ERR_INTERNEL,
              "encoder abnormal: block id %ld, type %s, encoded len %ld > origin len %ld\n",
              indata->block_id, get_typename(indata->btype).get(), outdata->current_len, indata->current_len);
    return true;
}

bool actuator_everything::decompress()
{
    int32_t len;
    coder_json cmeta;
    std::string md5;
    Json::Value meta;
    coder *decoder = nullptr;
    int32_t decode_len_expect;
    pbgz_stream_header header_meta(indata->buffer);

    /* 首先解析出当前block的meta信息，获取对应行的流信息和编码器等信息 */
    cmeta.decoder(indata->buffer + header_meta.get_bufferlen(), header_meta.get_datalen(), meta);

    indata->current_len = header_meta.get_bufferlen() + header_meta.get_datalen() + header_meta.get_bufferlen();

    if (meta["stream"]["coder"]["magic"] == "coder_bwt_cm")
    {
        decoder = new coder_bwt_cm(new coder_io(indata->get_curr(), indata->get_remain()));
        // decoder->set_level(meta["stream"]["coder"]["level"].asInt());

        coder_io block_io(indata->get_curr(), indata->get_remain());
        // block_io.set_level(1);
        // block_io.set_meta(meta);
        coder_bwt_cm coder(&block_io);
        decode_len_expect = meta["stream"]["srclen"].asInt();
        len = coder.decode_line((uint8_t *)(outdata->get_buffer()), decode_len_expect, UINT8_MAX, false);
    }

    else if (meta["stream"]["coder"]["magic"] == "coder_fc")
    {
        coder_io block_io(indata->get_curr(), indata->get_remain());
        block_io.meta = meta["stream"];
        block_io.meta["tot_dstlen"] = block_io.meta["dstlen"]; // 兼容字段，可以统一
        coder_fc coder(&block_io);
        decode_len_expect = meta["stream"]["srclen"].asInt();
        len = coder.decode_line((uint8_t *)(outdata->get_buffer()), decode_len_expect, UINT8_MAX, false);
    }

    // else if (meta["stream"]["coder"]["magic"] == "coder_ppmd")
    // {
    //     coder_ppmd  *ppmd = new coder_ppmd(new coder_io(indata->get_curr(), indata->get_remain()));
    //     // decoder->set_level(meta["stream"]["coder"]["level"].int_value());

    //     decode_len_expect = meta["stream"]["srclen"].int_value();
    //     len = ppmd->decode((uint8_t *)(outdata->get_buffer()), decode_len_expect);
    //     check_exit(len == decode_len_expect, ERR_INTERNEL, "check failed, expect %d, current decode len %d", decode_len_expect, len);
    // }

    else
        check_exit(false, ERR_INTERNEL, "undefine...");
    
    outdata->current_len += len;

    /* 检查md5 */
    calc_md5(md5, outdata->get_buffer(), outdata->current_len);
    check_exit(md5 == meta["md5"].asString(), ERR_INTERNEL, "md5 check failed");

    if (decoder)
        delete decoder;
    return true;
}