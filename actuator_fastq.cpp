#include <algorithm>
#include <memory> // Add this for std::shared_ptr
#include "actuator_fastq.h"
#include "manager.h"
#include "pbgz_file.h"
#include "actg.h"
#include "coder/simple_model.h"
#include "coder/clr.h"
#include "coder/coder_fc.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_json.h"
#include "coder/coder.h"
#include "coder/coder_ppmd.h"
#include "coder/coder_base.h"
#include "city.h"

actuator_fastq::actuator_fastq(const block_rough_ptr bptr_in, block_rough_ptr bptr_out, const reference *r)
{
    uint32_t n;
    this->indata = bptr_in;
    this->outdata = bptr_out;
    this->meta.clear();
    this->base_n_cnt = 0;
    this->baselen_min = UINT32_MAX;
    this->baselen_max = 0;
    this->idpos = nullptr;
    this->idpos_offset = 0;
    this->idsplit_syms = nullptr;
    this->idsplit_symslen = 0;
    this->idsplit_minlen = nullptr;
    this->idsplit_maxlen = nullptr;
    this->refgene = (reference *)r;
    this->ctype = CT_UNKNOW;
    this->meta.clear();
    this->outdata->block_id = this->indata->block_id;
    this->outdata->btype = this->indata->btype;
    this->meta["block_id"] = (Json::Value::Int)(bptr_out->block_id);
    this->somebuffer = nullptr;
}
actuator_fastq::~actuator_fastq()
{
    uint32_t n;
    if (idpos)
        free(idpos);
    if (idsplit_syms)
        free(idsplit_syms);
    if (idsplit_minlen)
        free(idsplit_minlen);
    if (idsplit_maxlen)
        free(idsplit_maxlen);
    if (somebuffer)
        free(somebuffer);
}

/* 详细分析fastq，格式不符合时返回false */
bool actuator_fastq::analyze_fastq()
{
    uint8_t ch, splitch, *p, *pline, *ptmp;
    uint32_t i, line, n, m, s = 0, e, l;
    uint32_t last, curr_len, curr_idx, last_tmp;
    std::vector<uint16_t> idpos_first_line;
    std::pair<uint8_t, uint32_t> qual_freq[256];

    line = indata->npos.size();
    p = indata->buffer;
    safe_alloc(indata->npos[0], uint8_t, idsplit_syms);
    for (n = 0; n < 256; n++) {
        qual_freq[n].first = n;
        qual_freq[n].second = 0;
    }

    for (i = 0; i < line; i++) {
        e = indata->npos[i];
        switch (i & 0x3)
        {
        case 0:                             /* id 行 , 先做域切割，域大于等于7时分开压缩 */
            if (UINT32_MAX != idpos_offset) /* 目前为止，id没有检查失败，继续分割 */
            {
                if (!idsplit_symslen) { /* 第一行id，还未解析出分隔符信息 */
                    for (n = s + 1; n <= e; n++) { /* 分割到换行符，需要考虑整行没有其他符号行的情形 */
                        ch = *(p + n);
                        if (idsplit.find_first_of(ch) != std::string::npos)
                        {
                            *(idsplit_syms + idsplit_symslen++) = ch;
                            idpos_first_line.push_back(n - s); /* 记录相对当前行开始的偏移 */
                        }
                    }
                    /* 将第一行id行分析信息拷贝至 idpos */
                    n = idsplit_symslen * ((line >> 2) + 1);
                    safe_alloc(n, uint16_t, idpos);
                    safe_alloc_init(idsplit_symslen * 2, uint32_t, idsplit_minlen, UINT32_MAX);
                    safe_alloc_init(idsplit_symslen * 2, uint32_t, idsplit_maxlen, 0);
                    last_tmp = 0;
                    for (n = 0; n < idpos_first_line.size(); n++) {
                        curr_len = idpos_first_line[n] - last_tmp - 1;
                        curr_idx = idpos_offset & (idsplit_symslen - 1);
                        if (curr_len < idsplit_minlen[curr_idx])
                            idsplit_minlen[curr_idx] = curr_len;
                        if (curr_len > idsplit_maxlen[curr_idx])
                             idsplit_maxlen[curr_idx] = curr_len;
                        *(idpos + idpos_offset++) = idpos_first_line[n];
                        last_tmp = idpos_first_line[n];
                    }

                } else { /* 非第一行，已经拿到分隔符模板，其他行都以模板为准去解析 */
                    
                    last = s + 1; /* 当前行上一个符对应的位置 */
                    last_tmp = 0;
                    for (m = 0; m < idsplit_symslen; m++) {
                        splitch = *(idsplit_syms + m);
                        for (n = last; n <= e; n++) {
                            if (*(p + n) == splitch) {
                                curr_len = (n - s) - last_tmp - 1;
                                curr_idx = m;
                                if (curr_len < idsplit_minlen[curr_idx])
                                    idsplit_minlen[curr_idx] = curr_len;
                                if (curr_len > idsplit_maxlen[curr_idx])
                                    idsplit_maxlen[curr_idx] = curr_len;
                                *(idpos + idpos_offset++) = n - s;
                                last_tmp = n - s;
                                last = n + 1;
                                break;
                            }
                        }
                        if (n > e) { /* 没有找到对应的分割符 */
                            idpos_offset = UINT32_MAX;
                            break;
                        }
                    }
                    /* 检查当前行是否找到了与模板对应对应长度的分割符 */
                    if (idpos_offset != (((i >> 2) + 1) *idsplit_symslen)) {
                        idpos_offset = UINT32_MAX;
                        break;
                    }
                }
            } 
            break;

        case 1: /* base 行分析 */
            l = e - s;
            if (l < baselen_min)
                baselen_min = l;
            if (l > baselen_max)
                baselen_max = l;
            for (n = s; n < e; n++)
            {
                switch (*(p + n))
                {
                case 'a':
                    break;
                case 'c':
                    break;
                case 't':
                    break;
                case 'g':
                    break;
                case 'n':
                    base_n_cnt++;
                    break;
                case 'A':
                    break;
                case 'C':
                    break;
                case 'T':
                    break;
                case 'G':
                    break;
                case 'N':
                    base_n_cnt++;
                    break;
                default:
                    return false;
                    break;
                }
            }
            break;

        case 2: /* comment 行分析 */
            pline = indata->buffer + indata->npos[i - 1] + 1;
            if (ctype != CT_OTHER) {
                if (ctype == CT_UNKNOW) {
                    if (*pline == '+' && indata->npos[i] - indata->npos[i - 1] == 2)
                        ctype = CT_JUST_PLUS;
                    else if((indata->npos[0] == indata->npos[i] - indata->npos[i - 1]) &&
                        !memcpy((uint8_t *)pline, indata->buffer,  indata->npos[0] - 1))
                        ctype = CT_SAME_AS_ID;
                    else
                        ctype == CT_OTHER;
                } else if (ctype == CT_JUST_PLUS) {
                    if (*pline != '+' || indata->npos[i] - indata->npos[i - 1] != 2)
                        ctype = CT_OTHER;
                } else if (ctype == CT_SAME_AS_ID) {
                    if((indata->npos[0] != indata->npos[i] - indata->npos[i - 1]) ||
                        memcpy((uint8_t *)pline, indata->buffer,  indata->npos[0] - 1))
                        ctype = CT_OTHER;
                }
            }
            break;

        case 3: /* quality分析 */
            for (n = s; n < e; n++)
                qual_freq[*(p + n)].second++;
            break;

        default:
            break;
        }
        s = e + 1;
    }

    /* 分析质量数 */
    /* 提升性能：排序后在模型中概率大的符号放在前面，减少减少匹配次数 */
    std::sort(qual_freq, qual_freq + 256, [](const std::pair<uint8_t, uint32_t> &e1, const std::pair<uint8_t, uint32_t> &e2) -> bool {
        return e1.second > e2.second;
    });

    for (i = 0; i < 256; i++) {
        if (qual_freq[i].second == 0)
            break;
        qual_freq_table.push_back(std::make_pair((qual_freq[i].first - '!'), 1));
    }
    return true;
}

/* 压缩id行 */
bool actuator_fastq::compress_id()
{
    Json::Value meta_id, meta_streams;
    Json::Value curr_stream;
    uint8_t *data, *pidlast;
    uint32_t i, m, n, line, curr_len, level;
    uint32_t curr_lineoffset; /* 指向当前line的偏移位置 */
    uint16_t *pcurr_idpos;    /* 指向当前idpos的位置 */
    uint32_t src_len, dst_len, total_srclen, total_dstlen;
    uint8_t idsplit_symbol[idsplit_symslen];

    level = manage::instance().get_zipinfo().get_clevel();
    line = indata->npos.size();
    total_srclen = total_dstlen = 0;

    if (UINT32_MAX != idpos_offset) {
        for (m = 0; m < idsplit_symslen; m++) { /* 检查是否有域内容为空 */
            idsplit_symbol[m] = *(indata->buffer + *(idpos + m));
            if (idsplit_minlen[m] == 0) {
                idpos_offset = UINT32_MAX;
                break;
            }
        }
    }
    if (UINT32_MAX == idpos_offset) { /* id检查无效，整块压缩 */
        src_len = 0;
        dst_len = outdata->current_len;
        std::shared_ptr<coder_io> id_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
        std::shared_ptr<coder_bwt_cm> id_cm = std::make_shared<coder_bwt_cm>(id_io.get());
        for (n = 0, i = 0; i < line;) {
            src_len += indata->npos[i] - n + 1;
            id_cm->encode_line(indata->buffer + n, indata->npos[i] - n + 1);
            i += 4;
            n = indata->npos[i - 1] + 1;
        }
        id_cm->encode_flush();
        outdata->current_len += id_io->data_len;
        dst_len = outdata->current_len - dst_len;
        total_srclen += src_len; total_dstlen += dst_len;

        curr_stream.clear();
        curr_stream["srclen"] = src_len;
        curr_stream["dstlen"] = dst_len;
        curr_stream["coder"] = id_io->meta;
        meta_streams.append(curr_stream);

        meta_id["tot_srclen"] = total_srclen;
        meta_id["tot_dstlen"] = total_dstlen;
        meta_id["streams"] = meta_streams;
        meta["id"] = meta_id;
        return (outdata->current_len <= outdata->buffer_size);
    }

    meta_id["split_sym"] = std::string((char *)idsplit_symbol, idsplit_symslen);

    line = indata->npos.size();
    for (i = 0; i < idsplit_symslen; i++) { /* 总共被分成了这么多条流 */
        src_len = 0;
        dst_len = outdata->current_len;

        /* id coder 1 */
        /* 判断是否是数字类型，这里只判断了首行 */
        {
            if (idsplit_maxlen[i] != idsplit_minlen[i]) { /* 不相等且第一行是数字 */
                pcurr_idpos = idpos + i;
                curr_lineoffset = 0;
                if (i == 0) {
                    data = indata->buffer + curr_lineoffset + 1; /* 1: 偏移开头的'@' */
                    curr_len = *pcurr_idpos - 1;
                } else {
                    data = indata->buffer + curr_lineoffset + *(pcurr_idpos - 1) + 1; /* 偏移当前分隔符 */
                    curr_len = *pcurr_idpos - *(pcurr_idpos - 1) - 1;
                }
                for (n = 0; n < curr_len; n++) {
                    if (*(data +n) < '0' || *(data +n) > '9')
                        break;
                }

                if (n >= curr_len) { /* 是数字类型 */
                    std::shared_ptr<coder_io>  id_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
                    std::shared_ptr<coder_bwt_cm> id_cm = std::make_shared<coder_bwt_cm>(id_io.get());

                    curr_lineoffset = 0;
                    pcurr_idpos = idpos + i;
                    for (n = 0; n < line;) { /* 取第n行对应i域的id */
                        if (i == 0) {
                            data = indata->buffer + curr_lineoffset;
                            curr_len = *pcurr_idpos + 1;
                        } else {
                            data = indata->buffer + curr_lineoffset + *(pcurr_idpos - 1) + 1;
                            curr_len = *pcurr_idpos - *(pcurr_idpos - 1);
                        }
                        id_cm->encode_line(data, curr_len);

                        src_len += curr_len;
                        pcurr_idpos += idsplit_symslen;
                        n += 4;
                        curr_lineoffset = indata->npos[n - 1] + 1; /* 1: 表示要偏移一个回车 */
                    }
                    id_cm->encode_flush();
                    outdata->current_len += id_io->data_len;
                    dst_len = outdata->current_len - dst_len;
                    total_srclen += src_len; total_dstlen += dst_len;

                    curr_stream.clear();
                    curr_stream["srclen"] = src_len;
                    curr_stream["dstlen"] = dst_len;
                    curr_stream["coder"] = id_io->meta;
                    meta_streams.append(curr_stream);
                    continue;
                }
            }
        }

        /* id coder 2 */
        std::shared_ptr<coder_io> id_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
        std::shared_ptr<coder_affix_match> id_am = std::make_shared<coder_affix_match>(id_io.get());
        
        curr_lineoffset = 0;
        pcurr_idpos = idpos + i;
        for (n = 0; n < line;) { /* 取第n行对应i域的id */
            if (i == 0) {
                data = indata->buffer + curr_lineoffset;
                curr_len = *pcurr_idpos + 1;
            } else {
                data = indata->buffer + curr_lineoffset + *(pcurr_idpos - 1) + 1;
                curr_len = *pcurr_idpos - *(pcurr_idpos - 1);
            }

            id_am->encode_line(data, curr_len);
            src_len += curr_len;
            pcurr_idpos += idsplit_symslen;
            n += 4;
            curr_lineoffset = indata->npos[n - 1] + 1; /* 1: 表示要偏移一个回车 */
        }

        id_am->encode_flush();
        outdata->current_len += id_io->data_len;
        dst_len = outdata->current_len - dst_len;
        total_srclen += src_len; total_dstlen += dst_len;

        curr_stream.clear();
        curr_stream["srclen"] = src_len;
        curr_stream["dstlen"] = dst_len;
        curr_stream["coder"] = id_io->meta;
        meta_streams.append(curr_stream);       
    }

    meta_id["tot_srclen"] = total_srclen;
    meta_id["tot_dstlen"] = total_dstlen;
    meta_id["streams"] = meta_streams;
    meta["id"] = meta_id;

    return (outdata->current_len <= outdata->buffer_size);
}
/* 解压id行初始化 */
bool actuator_fastq::initialize_decode_id(std::vector<coder *> &id_decoders)
{
    int32_t n;
    std::string coder_name;
    Json::Value id_streams;
    pbgz_stream_header header_meta(indata->buffer);

    /* 指向id待解压的数据位置 */
    indata->current_len = header_meta.get_bufferlen() + header_meta.get_datalen() + header_meta.get_bufferlen();
    outdata->current_len = 0;
    id_streams = meta["id"]["streams"];
    check(id_streams.isArray(), false, "prase id streams failed");

    for (n = 0; n < id_streams.size(); n++)
    {
        coder_name = id_streams[n]["coder"]["magic"].asString();
        if (coder_name == "coder_affix_match")
        {
            id_decoders.push_back(new coder_affix_match(new coder_io(indata->get_curr(), indata->get_remain())));
            id_decoders.back()->set_level(id_streams[n]["coder"]["level"].asInt());
        }
        else if (coder_name == "coder_bwt_cm")
        {
            id_decoders.push_back(new coder_bwt_cm(new coder_io(indata->get_curr(), indata->get_remain())));
            id_decoders.back()->set_level(id_streams[n]["coder"]["level"].asInt());
        }
        else
            check_exit(false, ERR_INTERNEL, "undefined coder: %s", coder_name.c_str());

        indata->current_len += id_streams[n]["dstlen"].asUInt();
    }
    return true;
}

/* 解压base行初始化 */
bool actuator_fastq::initialize_decode_base(coder *&base_decoder)
{
    uint8_t *p, *ps;
    uint32_t src_len, dst_len, id;
    Json::Value base_meta, meta_streams;
    uint32_t n, len_max_base, offset = 0;

    base_meta = meta["base"];

    if (refgene)
    {
        meta_streams = base_meta["streams"];
        len_max_base = base_meta["lenmax"].asUInt();
        n = len_max_base;  /* sub stream 1, strip n  */
        n += len_max_base; /* for somebuffer_refe_stretch */
        for (id = 1; id < meta_streams.size(); id++)
            n += meta_streams[id]["srclen"].asUInt();
        safe_alloc(n, uint8_t, somebuffer);
        ps = somebuffer;
        somebuffer_base_stripn = ps;
        ps += len_max_base;
        somebuffer_refe_stretch = ps;
        ps += len_max_base;

        /* check sub streams 1 */
        id = 0;
        check_exit(meta_streams[id]["sname"].asString() == "m",
                   ERR_INTERNEL, "check sub stream failed: %s", meta_streams[id]["sname"].asString().c_str());

        if (meta_streams[id]["coder"]["magic"].asString() == "coder_bwt_cm")
            base_decoder = new coder_bwt_cm(new coder_io(indata->get_curr(), indata->get_remain()));
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
        offset += meta_streams[id]["dstlen"].asInt();

        /* check sub streams 2 */
        id = 1;
        somebuffer_base_mpos = nullptr;
        check_exit(meta_streams[id]["sname"].asString() == "mpos",
                   ERR_INTERNEL, "check sub stream failed: %s", meta_streams[id]["sname"].asString().c_str());
        if (meta_streams[id]["coder"]["magic"].asString() == "coder_bwt_cm")
        {
            p = indata->get_curr() + offset;
            src_len = meta_streams[id]["srclen"].asUInt();
            dst_len = meta_streams[id]["dstlen"].asInt();
            somebuffer_base_mpos = (uint64_t *)ps;
            ps += src_len;
            coder_io pos_io(p, dst_len);
            coder_bwt_cm pos_cm(&pos_io);
            pos_cm.decode_line((uint8_t *)somebuffer_base_mpos, src_len, UINT8_MAX, false);
        }
        else if (meta_streams[id]["coder"]["magic"].asString() == "coder_fc")
        {
            p = indata->get_curr() + offset;
            src_len = meta_streams[id]["srclen"].asUInt();
            dst_len = meta_streams[id]["dstlen"].asInt();
            somebuffer_base_mpos = (uint64_t *)ps;
            ps += src_len;
            coder_io pos_io(p, dst_len);
            pos_io.meta = meta_streams[id];
            pos_io.meta["tot_dstlen"] = meta_streams[id]["dstlen"]; // 兼容字段，可以统一
            coder_fc pos_cm(&pos_io);
            pos_cm.decode_line((uint8_t *)somebuffer_base_mpos, src_len, UINT8_MAX, false);
        }
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
        offset += meta_streams[id]["dstlen"].asInt();

        /* check sub stream 3 */
        id = 2;
        somebuffer_base_mpair = nullptr;
        check_exit(meta_streams[id]["sname"].asString() == "mpair",
                   ERR_INTERNEL, "check sub stream failed: %s", meta_streams[id]["sname"].asString().c_str());
        if (meta_streams[id]["coder"]["magic"].asString() == "coder_bwt_cm")
        {
            p = indata->get_curr() + offset;
            src_len = meta_streams[id]["srclen"].asUInt();
            dst_len = meta_streams[id]["dstlen"].asInt();
            somebuffer_base_mpair = ps;
            ps += src_len;
            coder_io pair_io(p, dst_len);
            coder_bwt_cm pair_cm(&pair_io);
            pair_cm.decode_line(somebuffer_base_mpair, src_len, UINT8_MAX, false);
        }
        else if (meta_streams[id]["coder"]["magic"].asString() == "coder_fc")
        {
            p = indata->get_curr() + offset;
            src_len = meta_streams[id]["srclen"].asUInt();
            dst_len = meta_streams[id]["dstlen"].asInt();
            somebuffer_base_mpair = ps;
            ps += src_len;
            coder_io pair_io(p, dst_len);
            pair_io.meta = meta_streams[id];
            pair_io.meta["tot_dstlen"] = meta_streams[id]["dstlen"]; // 兼容字段，可以统一
            coder_fc pair_cm(&pair_io);
            pair_cm.decode_line(somebuffer_base_mpair, src_len, UINT8_MAX, false);
        }
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
        offset += meta_streams[id]["dstlen"].asInt();

        /* check sub stream 4 */
        somebuffer_base_npos = nullptr;
        if (base_meta["ncount"].asUInt())
        {
            id = 3;
            check_exit(meta_streams[id]["sname"].asString() == "npos",
                       ERR_INTERNEL, "check sub stream failed: %s", meta_streams[id]["sname"].asString().c_str());
            if (meta_streams[id]["coder"]["magic"].asString() == "coder_bwt_cm")
            {
                p = indata->get_curr() + offset;
                src_len = meta_streams[id]["srclen"].asUInt();
                dst_len = meta_streams[id]["dstlen"].asInt();
                somebuffer_base_npos = (uint32_t *)ps;
                ps += src_len;
                coder_io npos_io(p, dst_len);
                coder_bwt_cm npos_cm(&npos_io);
                npos_cm.decode_line((uint8_t *)somebuffer_base_npos, src_len, UINT8_MAX, false);
            }
            else
                check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
            offset += meta_streams[id]["dstlen"].asInt();
        }

        /* check sub stream 5 */
        somebuffer_base_len2 = nullptr;
        if (base_meta["lenmin"].asUInt() != base_meta["lenmax"].asUInt())
        {
            id = 4;
            check_exit(meta_streams[id]["sname"].asString() == "baselen",
                       ERR_INTERNEL, "check sub stream failed: %s", meta_streams[id]["sname"].asString().c_str());
            if (meta_streams[id]["coder"]["magic"].asString() == "coder_bwt_cm")
            {
                p = indata->get_curr() + offset;
                src_len = meta_streams[id]["srclen"].asUInt();
                dst_len = meta_streams[id]["dstlen"].asInt();
                somebuffer_base_len2 = (uint16_t *)ps;
                ps += src_len;
                coder_io len_io(p, dst_len);
                coder_bwt_cm len_cm(&len_io);
                len_cm.decode_line((uint8_t *)somebuffer_base_len2, src_len, UINT8_MAX, false);
            }
            else
                check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
            offset += meta_streams[id]["dstlen"].asInt();
        }
    }
    else
    {
        if (base_meta["coder"]["magic"].asString() == "coder_fc")
        {
            src_len = base_meta["tot_srclen"].asUInt();
            dst_len = base_meta["tot_dstlen"].asUInt();

            coder_io base_io(indata->get_curr(), dst_len);
            base_io.meta = base_meta;
            coder_fc base_cm(&base_io);
            check_exit(outdata->get_remain() > src_len,
                       ERR_INTERNEL, "lack of buffer: %d should > %d", outdata->get_remain(), src_len);
            base_cm.decode_line((uint8_t *)(outdata->get_buffer()) + outdata->buffer_size - src_len, src_len, UINT8_MAX, false);
        }
        else if (base_meta["coder"]["magic"].asString() == "coder_bwt_cm")
            base_decoder = new coder_bwt_cm(new coder_io(indata->get_curr(), indata->get_remain()));
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
    }
    indata->current_len += base_meta["tot_dstlen"].asUInt();
    return true;
}

/* 解压comment行初始化 */
bool actuator_fastq::initialize_decode_comment(coder *&com_decoder)
{
    Json::Value meta_comment, meta_streams;
    meta_comment = meta["comment"];

    com_decoder = nullptr;
    if (meta_comment["type"].asString() == "plus")
        ctype = CT_JUST_PLUS;
    else if (meta_comment["type"].asString() == "same")
        ctype = CT_SAME_AS_ID;
    else if (meta_comment["type"].asString() == "other")
    {
        ctype = CT_OTHER;
        meta_streams = meta_comment["streams"];
        check_exit(meta_streams.size() == 1, ERR_INTERNEL,
                   "check comment decoder count failed: %u", meta_streams.size());
        check_exit(meta_streams["sname"].asString() == "comment",
                   ERR_INTERNEL, "check sub stream failed: %s", meta_streams["sname"].asString().c_str());

        if (meta_streams["coder"]["magic"].asString() == "coder_bwt_cm")
        {
            com_decoder = new coder_bwt_cm(new coder_io(indata->get_curr(), indata->get_remain()));
            indata->current_len += meta_streams["dstlen"].asInt();
        }
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
    }
    else
        check_exit(false, ERR_INTERNEL,
                   "check comment type failed: %s", meta_comment["type"].asString().c_str());
    return true;
}

/* 解压quality行初始化 */
bool actuator_fastq::initialize_decode_quality(coder_qual *&qual_decoder)
{
    uint32_t len_src, len_get, flen, n;
    Json::Value meta_qual, meta_streams;
    meta_qual = meta["quality"];
    meta_streams = meta_qual["streams"];

    check_exit(meta_streams.size() == 2, ERR_INTERNEL,
               "check quality decoder count failed: %u", meta_streams.size());
    len_src = meta_streams[1]["srclen"].asUInt();
    flen = len_src / sizeof(uint16_t);
    uint16_t qual_freq_arr[flen];

    /* 先解析 qual_freq_table */
    if (meta_streams[1]["coder"]["magic"].asString() == "coder_bwt_cm")
    {
        coder_io qual_io(indata->get_curr() + meta_streams[0]["dstlen"].asInt(),
                         indata->get_remain() - meta_streams[0]["dstlen"].asInt());
        coder_bwt_cm *qual_freq_cm = new coder_bwt_cm(&qual_io);

        len_get = qual_freq_cm->decode_line((uint8_t *)qual_freq_arr, len_src, UINT8_MAX, false);
        check_exit(len_get == len_src, ERR_INTERNEL,
                   "decode quality failed: expect lenght %u, actual %u", len_src, len_get);
        for (n = 0; n < flen; n += 2)
            qual_freq_table.push_back(std::make_pair(qual_freq_arr[n], qual_freq_arr[n + 1]));
        delete qual_freq_cm;

        if (meta_streams[0]["coder"]["magic"].asString() == "coder_qual")
        {
            qual_decoder = new coder_qual(new coder_io(
                                              indata->get_curr(), indata->get_remain() - meta_streams[1]["dstlen"].asInt()),
                                          true, qual_freq_table);
        }
        else
            check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");

        indata->current_len += meta_qual["tot_dstlen"].asUInt();
    }
    else
        check_exit(false, ERR_INTERNEL, "check sub stream failed: coder name unmatch");
    return true;
}

/* 二代数据匹配reference */
inline void actuator_fastq::mapping_gen2(const uint8_t *base, uint32_t base_len, uint8_t *&out, uint32_t &out_len, uint64_t &mpos, uint8_t &mdir)
{
    mapping_t mt[4];
    uint32_t s, e, n, m, o, l, len_squash[2], loffset;
    uint32_t align4[2], align4_curr;
    uint8_t *psquash, *psquash_refe, *p, ch;
    const uint32_t bg_len = refgene->get_bglen();
    const uint32_t bg_mid = bg_len >> 1;
    const uint8_t *pseq[2] = {base, somebuffer_basepair};
    uint32_t hash32, pos_cnts, *pos_vals, total;
    const uint8_t bg_is_unalign4 = !!(bg_len & 0x3);
    const uint32_t bg_align4_len = (bg_len >> 2) << 2;
    const uint32_t len_bgs = (bg_len >> 2) + bg_is_unalign4;

    uint32_t base_squash_align4 = (base_len >> 2) + !!(base_len & 0x3);
    uint32_t match_pair, match_pair_origin;
    uint32_t match_pos, unmatchs = base_len;
    uint8_t *prefe_squash = (uint8_t *)(refgene->get_squash());
    int64_t refe_squashlen = refgene->get_squashlen();
    uint32_t best_pos_inrefe = UINT32_MAX, best_is_pair = 0; /* 记录当前最好匹配时对应的reference中的位置和正负链方向 */
    uint32_t best_unmatchs = UINT32_MAX, best_align4;
    uint32_t best_pos = UINT32_MAX; /* 转换为reference原始碱基对应的位置 */
    uint64_t xsquash, xsquash_match, xsquash_macth_refe;
    // const uint64_t xsquash_tab[2] = {0xFCFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFF3F};
    const uint64_t xsquash_tab[2] = {0xFCFFFFFFFFFFFFFF, 0xFCFFFFFFFFFFFFFF};
    const uint8_t bg_unalign4_len = 1; // (refgene->get_bgstep() - refgene->get_bglen()) & 0x3;

    /* case 1: base长度不大于reference索引对应的碱基长度 */
    if (base_len <= bg_len)
    {
        actg_encode(base, out, base_len);
        out_len = base_len;
        mpos = 0;
        mdir = 2; /*  解压时先判断mdir，如果为2说明没有匹配上 */
        return;
    }

    /* case 2: base长度大于reference索引对应的碱基长度 */
    e = base_len - bg_len;
    actg_pair(somebuffer_basepair, base, base_len);

    /* 计算align4  squash buffer and pair squash buffer */
    for (n = 0; n < 4; n++)
    {
        len_squash[0] = (base_len - n) >> 2;
        total = len_squash[0] << 2;
        mt[n].ps_l_unalign_len[0] = n;
        for (m = 0; m < n; m++)
            mt[n].ps_l_unalign[0][m] = ((*(pseq[0] + m)) >> 1) & 0x3; /* 左边没有4对齐的碱基squash值 */
        mt[n].ps_r_unalign_len[0] = base_len - n - total;
        for (m = 0; m < mt[n].ps_r_unalign_len[0]; m++)
            mt[n].ps_r_unalign[0][m] = ((*(pseq[0] + base_len - mt[n].ps_r_unalign_len[0] + m)) >> 1) & 0x3;
        actg_squash(pseq[0] + n, total, somebuffer_base_squash[n]);

        len_squash[1] = (base_len - n + 1) >> 2; /*  在最右边补一个字符使之与32对齐，如果match到mt[0]的pair，且mt[0] offset为0时需要处理最后一个字符 */
        total = len_squash[1] << 2;
        mt[n].ps_l_unalign_len[1] = base_len + 1 - n - total;
        for (m = 0; m < mt[n].ps_l_unalign_len[1]; m++)
            mt[n].ps_l_unalign[1][m] = ((*(pseq[1] + m)) >> 1) & 0x3;
        mt[n].ps_r_unalign_len[1] = (n == 0) ? 0 : (n - 1); /* 右边补一个到key对齐，所以需要减1 */
        for (m = 0; m < mt[n].ps_r_unalign_len[1]; m++)
            mt[n].ps_r_unalign[1][m] = ((*(pseq[1] + base_len - mt[n].ps_r_unalign_len[1] + m) >> 1) & 0x3);
        actg_squash(pseq[1] + mt[n].ps_l_unalign_len[1], total, somebuffer_basepair_squash[n]);

        mt[n].set(somebuffer_base_squash[n], len_squash[0], /* 建立base squash与对应pair base squash的对应关系 */
                  somebuffer_basepair_squash[n] + len_squash[1] - len_bgs, len_squash[1], 0);

        /* do mapping */
        align4_curr = n & 0x3;
        match_pair_origin = (*(pseq[0] + n + bg_mid) < *(pseq[1] + base_len - bg_len - n + bg_mid));

        psquash = mt[align4_curr].get_squash(match_pair_origin);
        xsquash = *((uint64_t *)(psquash));
        xsquash &= xsquash_tab[match_pair_origin];
        hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
        pos_vals = (uint32_t *)(refgene->query_pos(hash32, pos_cnts));

        for (o = 0; o < pos_cnts; o++)
        {
            match_pos = *pos_vals++;
            match_pair = ((match_pos & 0x80000000) >> 31) ^ match_pair_origin;
            match_pos = (match_pos & 0x7FFFFFFF) << 3; /* to squash reference pos */

            /* check left and right boundary simply */
            if (match_pos + base_squash_align4 >= refe_squashlen || match_pos < base_squash_align4)
                continue;

            loffset = (match_pair) ? (mt[align4_curr].ps_len[1] - mt[align4_curr].offset - len_bgs) : (mt[align4_curr].offset);

            /* caculate unmatch count */
            psquash = mt[align4_curr].get_squash(match_pair) - loffset;
            psquash_refe = prefe_squash + match_pos - loffset;

            xsquash_match = ((*((uint64_t *)(mt[align4_curr].get_squash(match_pair)))) & xsquash_tab[match_pair]);
            xsquash_macth_refe = ((*((uint64_t *)(prefe_squash + match_pos))) & xsquash_tab[match_pair]);
            if (xsquash_match != xsquash_macth_refe) /* key is not same, skip */
                continue;

            /* align 4 */
            unmatchs = actg_squash_diffcnt(psquash, psquash_refe, mt[align4_curr].ps_len[match_pair]);
#ifdef ANALYZ_STREAMS
            mapping_cnt++;
#endif
            if (unmatchs >= best_unmatchs)
                continue;

            /*  因为这种情况在最右边补一个字符使之与32对齐：如果match到mt[0]的pair，且mt[0] offset为0时需要处理最后一个字符 */
            unmatchs -= (match_pair && (mt[align4_curr].offset == 0)) ? ((xsquash_match & 0x80000000000000) != (xsquash_macth_refe & 0x80000000000000)) : 0;

            /* left unalign */
            l = mt[align4_curr].ps_l_unalign_len[match_pair];
            for (ch = *(psquash_refe - 1), m = 0; m < mt[align4_curr].ps_l_unalign_len[match_pair]; m++)
            {
                unmatchs += ((ch >> (m << 1)) & 0x3) != (mt[align4_curr].ps_l_unalign[match_pair][l - 1]);
                l--;
            }
            /* right unalign */
            l = mt[align4_curr].ps_len[match_pair];
            for (ch = *(psquash_refe + l), m = 0; m < mt[align4_curr].ps_r_unalign_len[match_pair]; m++)
                unmatchs += ((ch >> (6 - (m << 1)) & 0x3) != (mt[align4_curr].ps_r_unalign[match_pair][m]));
            if (unmatchs < best_unmatchs)
            {
                best_pos = (match_pos << 2) - (loffset << 2) - mt[align4_curr].ps_l_unalign_len[match_pair];
                best_pos_inrefe = match_pos - loffset;
                best_is_pair = match_pair;
                best_align4 = align4_curr;
                best_unmatchs = unmatchs;
            }
            if (best_unmatchs <= MAPPED_THRESHOLD_GEN2)
                break;
        }

        if (best_unmatchs <= MAPPED_THRESHOLD_GEN2)
            break;
        mt[align4_curr].inc_offset();
    }

    if (best_unmatchs > MAPPED_THRESHOLD_GEN2)
    { /* 继续mapping */
        for (n = 4; n <= e; n++)
        {
            /* do mapping */
            align4_curr = n & 0x3;
            match_pair_origin = (*(pseq[0] + n + bg_mid) < *(pseq[1] + base_len - bg_len - n + bg_mid));

            psquash = mt[align4_curr].get_squash(match_pair_origin);
            xsquash = *((uint64_t *)(psquash));
            xsquash &= xsquash_tab[match_pair_origin];
            hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
            pos_vals = (uint32_t *)(refgene->query_pos(hash32, pos_cnts));

            for (o = 0; o < pos_cnts; o++)
            {
                match_pos = *pos_vals++;
                match_pair = ((match_pos & 0x80000000) >> 31) ^ match_pair_origin;
                match_pos = (match_pos & 0x7FFFFFFF) << 3; /* to squash reference pos */

                /* check left and right boundary simply */
                if (match_pos + base_squash_align4 >= refe_squashlen || match_pos < base_squash_align4)
                    continue;

                loffset = (match_pair) ? (mt[align4_curr].ps_len[1] - mt[align4_curr].offset - len_bgs) : (mt[align4_curr].offset);

                /* caculate unmatch count */
                psquash = mt[align4_curr].get_squash(match_pair) - loffset;
                psquash_refe = prefe_squash + match_pos - loffset;

                xsquash_match = ((*((uint64_t *)(mt[align4_curr].get_squash(match_pair)))) & xsquash_tab[match_pair]);
                xsquash_macth_refe = ((*((uint64_t *)(prefe_squash + match_pos))) & xsquash_tab[match_pair]);
                if (xsquash_match != xsquash_macth_refe) /* key is not same, skip */
                    continue;

                /* align 4 */
                unmatchs = actg_squash_diffcnt(psquash, psquash_refe, mt[align4_curr].ps_len[match_pair]);
#ifdef ANALYZ_STREAMS
                mapping_cnt++;
#endif
                if (unmatchs >= best_unmatchs)
                    continue;

                /*  因为这种情况在最右边补一个字符使之与32对齐：如果match到mt[0]的pair，且mt[0] offset为0时需要处理最后一个字符 */
                unmatchs -= (match_pair && (mt[align4_curr].offset == 0)) ? ((xsquash_match & 0x80000000000000) != (xsquash_macth_refe & 0x80000000000000)) : 0;

                /* left unalign */
                l = mt[align4_curr].ps_l_unalign_len[match_pair];
                for (ch = *(psquash_refe - 1), m = 0; m < mt[align4_curr].ps_l_unalign_len[match_pair]; m++)
                {
                    unmatchs += ((ch >> (m << 1)) & 0x3) != (mt[align4_curr].ps_l_unalign[match_pair][l - 1]);
                    l--;
                }
                /* right unalign */
                l = mt[align4_curr].ps_len[match_pair];
                for (ch = *(psquash_refe + l), m = 0; m < mt[align4_curr].ps_r_unalign_len[match_pair]; m++)
                    unmatchs += ((ch >> (6 - (m << 1)) & 0x3) != (mt[align4_curr].ps_r_unalign[match_pair][m]));
                if (unmatchs < best_unmatchs)
                {
                    best_pos = (match_pos << 2) - (loffset << 2) - mt[align4_curr].ps_l_unalign_len[match_pair];
                    best_pos_inrefe = match_pos - loffset;
                    best_is_pair = match_pair;
                    best_align4 = align4_curr;
                    best_unmatchs = unmatchs;
                }
                if (best_unmatchs <= MAPPED_THRESHOLD_GEN2)
                    break;
            }

            if (best_unmatchs <= MAPPED_THRESHOLD_GEN2)
                break;
            mt[align4_curr].inc_offset();
        }
    }

    /* calc the result of base or base pair mapping with the match pos reference */
    out_len = 0;
    if (best_unmatchs != UINT32_MAX) /* get pos in reference table */
    {
        psquash = (best_is_pair) ? (mt[best_align4].ps[1] - (mt[best_align4].ps_len[1] - len_bgs)) : (mt[best_align4].ps[0]);
        psquash_refe = prefe_squash + best_pos_inrefe;

        /* left unalign */
        n = o = l = mt[best_align4].ps_l_unalign_len[best_is_pair];
        for (ch = *(psquash_refe - 1), m = 0; m < mt[best_align4].ps_l_unalign_len[best_is_pair]; m++)
        {
            /* xor  method*/
            out[--o] = (((ch >> (m << 1)) & 0x3) ^ (mt[best_align4].ps_l_unalign[best_is_pair][l - 1]));
            l--;
        }
        out_len += n; /* 注意字节顺序 */
        /* middle align */
        // out_len += actg_stretch_mapping(psquash, psquash_refe, mt[best_align4].ps_len[best_is_pair], out + out_len);
        /* xor  method*/
        out_len += actg_stretch_mapping_xor(psquash, psquash_refe, mt[best_align4].ps_len[best_is_pair], out + out_len);

        out_len -= (best_is_pair && best_align4 == 0);
        /* right unalign */
        l = mt[best_align4].ps_len[best_is_pair];
        for (ch = *(psquash_refe + l), m = 0; m < mt[best_align4].ps_r_unalign_len[best_is_pair]; m++)
        {
            // out[out_len++] = ((ch >> (6 - (m << 1)) & 0x3) == (mt[best_align4].ps_r_unalign[best_is_pair][m])) ? 0 : mt[best_align4].ps_r_unalign[best_is_pair][m];
            /* xor method */
            out[out_len++] = ((ch >> (6 - (m << 1)) & 0x3) ^ (mt[best_align4].ps_r_unalign[best_is_pair][m]));
        }
    }
    else
    { /* not match valid pos */
        actg_encode(base, out, base_len);
        out_len = base_len;

        best_pos = 0;
        best_is_pair = 2; ///*  解压时先判断mdir，如果为2说明没有匹配上 */
    }

    mpos = best_pos;
    mdir = best_is_pair;
}

/* 三代数据匹配reference */
inline void actuator_fastq::mapping_gen3(const uint8_t *base, uint32_t base_len, uint8_t *&out, uint32_t &out_len, uint64_t &mpos, uint8_t &mdir)
{
    check_exit(false, ERR_INTERNEL, "undefined");
} 

/* 压缩base行 */
bool actuator_fastq::compress_base()
{
    uint8_t *p, *pb, mdir, ch, *pbuff;
    uint32_t i, line, n, s = 0, e, out_len, offset = 0;
    int64_t mpos, noffset = 0, currpos = 0, base_src_len, l;
    uint32_t line4 = (indata->npos.size() >> 2);
    bool enc_baselen = (baselen_min != baselen_max);
    uint32_t src_len, dst_len, total_srclen, total_dstlen;
    Json::Value meta_base, meta_streams, meta_subs;

    line = indata->npos.size();
    p = indata->buffer;
    total_srclen = total_dstlen = 0;

    check_exit(is_gen2, ERR_INTERNEL, "not support gen3 fastq compress");
#ifdef ANALYZ_STREAMS
    timer cost_ms(true);
#endif

    if (refgene) { /* 使用reference压缩 */

        src_len = dst_len = 0;
        std::shared_ptr<coder_io> match_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
        std::shared_ptr<coder_bwt_cm> match_cm = std::make_shared<coder_bwt_cm>(match_io.get());

        for (i = 1; i < line; i+= 4)
        {
            e = indata->npos[i];
            s = indata->npos[i - 1] + 1;

            pbuff = somebuffer_base_stripn;
            for (n = s; n < e; n++)
            {
                ch = *(p + n);
                if (ch == 'n' || ch == 'N') {
                    *(somebuffer_base_npos + noffset) = currpos + n - s;
                    noffset++;
                } else {
                    *pbuff = ch;
                    pbuff++;
                }
            }

            (this->*actuator_fastq::mapping)(somebuffer_base_stripn, pbuff - somebuffer_base_stripn,
                                                somebuffer_base_mapped, out_len, somebuffer_base_mpos[offset], somebuffer_base_mpair[offset]);
            src_len += pbuff - somebuffer_base_stripn;
            match_cm->encode_line(somebuffer_base_mapped, out_len);

#ifdef ANALYZ_STREAMS
            for (int32_t _ = 0; _ < out_len; _++)
                actg_matched += !(somebuffer_base_mapped[_]);
            actg_total += out_len;
#endif

            this->refgene->update_matched(somebuffer_base_mpos[offset], out_len);

            if (enc_baselen) 
                somebuffer_base_len2[offset] = e - s - baselen_min;
            offset++;
            currpos += e - s;
        }
        /* 第一条子流：reads与reference的match流 */
        match_cm->encode_flush();
        meta_subs.clear();
        meta_subs["srclen"] =  (Json::Value::UInt)src_len;
        meta_subs["dstlen"] = (Json::Value::Int)(match_io->data_len);
        meta_subs["coder"] = match_io->meta;
        meta_subs["sname"] = "m";
        meta_streams.append(meta_subs);
        outdata->current_len += match_io->data_len;
        total_srclen += src_len;
        total_dstlen += match_io->data_len;
#ifdef ANALYZ_STREAMS
        fprintf(stderr, "\n\tsub streams: ---[match xor]--- src size %llu, dst size %llu, cost ms %llu\n\n",
               src_len, match_io->data_len, cost_ms.elapsed());
        cost_ms.reset();
#endif

        /* 第二条子流：reads与reference的match的位置流 */
        std::shared_ptr<coder_io> pos_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
        src_len = (line4 << 3);
        if (src_len > FC_MIN_LEN && src_len < FC_MAX_LEN)
        {
            std::shared_ptr<coder_fc> pos_cm = std::make_shared<coder_fc>(pos_io.get());
            pos_cm->encode_line((uint8_t *)somebuffer_base_mpos, src_len);
            pos_cm->encode_flush();
        } else {
            std::shared_ptr<coder_bwt_cm> pos_cm = std::make_shared<coder_bwt_cm>(pos_io.get());
            pos_cm->encode_line((uint8_t *)somebuffer_base_mpos, src_len);
            pos_cm->encode_flush();
        }
        meta_subs.clear();
        meta_subs["srclen"] =  (Json::Value::UInt)src_len;
        meta_subs["dstlen"] =  (Json::Value::Int)pos_io->data_len;
        meta_subs["coder"] = pos_io->meta;
        meta_subs["sname"] = "mpos";
        meta_streams.append(meta_subs);
        outdata->current_len += pos_io->data_len;
        total_srclen += src_len;
        total_dstlen += pos_io->data_len;
#ifdef ANALYZ_STREAMS
        printf("\n\tsub streams: ---[match pos]--- src size %llu, dst size %llu, cost ms %llu\n\n",
               src_len, pos_io->data_len, cost_ms.elapsed());
        cost_ms.reset();
#endif

        /* 第三条子流：reads与reference的match的pair标识流 */
        auto pair_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
        src_len = line4;
        if (src_len > FC_MIN_LEN && src_len < FC_MAX_LEN)
        {
            std::shared_ptr<coder_fc> sub_coder = std::make_shared<coder_fc>(pair_io.get());
            sub_coder->encode_line((uint8_t *)somebuffer_base_mpair, src_len);
            sub_coder->encode_flush();
        } else {
            std::shared_ptr<coder_bwt_cm> sub_coder = std::make_shared<coder_bwt_cm>(pair_io.get());
            sub_coder->encode_line((uint8_t *)somebuffer_base_mpair, src_len);
            sub_coder->encode_flush();
        }
        meta_subs.clear();
        meta_subs["srclen"] = (Json::Value::UInt)src_len;
        meta_subs["dstlen"] = (Json::Value::Int)pair_io->data_len;
        meta_subs["coder"] = pair_io->meta;
        meta_subs["sname"] = "mpair";
        meta_streams.append(meta_subs);
        outdata->current_len += pair_io->data_len;
        total_srclen += src_len;
        total_dstlen += pair_io->data_len;

        /* 第四条子流：reads中所有N的位置 */
        if (base_n_cnt > 0) {
            std::shared_ptr<coder_io> npos_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
            src_len = (base_n_cnt << 2);
            if (false)
            {
                std::shared_ptr<coder_fc> sub_coder = std::make_shared<coder_fc>(npos_io.get());
                sub_coder->encode_line((uint8_t *)somebuffer_base_npos, src_len);
                sub_coder->encode_flush();
            } else {
                std::shared_ptr<coder_bwt_cm> sub_coder = std::make_shared<coder_bwt_cm>(npos_io.get());
                sub_coder->encode_line((uint8_t *)somebuffer_base_npos, src_len);
                sub_coder->encode_flush();
            }
            meta_subs.clear();
            meta_subs["srclen"] = (Json::Value::UInt)src_len;
            meta_subs["dstlen"] = (Json::Value::Int)npos_io->data_len;
            meta_subs["coder"] = npos_io->meta;
            meta_subs["sname"] = "npos";
            meta_streams.append(meta_subs);
            outdata->current_len += npos_io->data_len;
            total_srclen += src_len;
            total_dstlen += npos_io->data_len;
        }
        meta_base["ncount"] = (Json::Value::UInt)base_n_cnt;

        /* 第五条流：每行base的长度 */
        if (enc_baselen) {
            std::shared_ptr<coder_io> len_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
            src_len = (line4 << 1);
            if (false)
            {
                std::shared_ptr<coder_fc> sub_coder = std::make_shared<coder_fc>(len_io.get());
                sub_coder->encode_line((uint8_t *)somebuffer_base_len2, src_len);
                sub_coder->encode_flush();
            } else {
                std::shared_ptr<coder_bwt_cm> sub_coder = std::make_shared<coder_bwt_cm>(len_io.get());
                sub_coder->encode_line((uint8_t *)somebuffer_base_len2, src_len);
                sub_coder->encode_flush();
            }
            meta_subs.clear();
            meta_subs["srclen"] = (Json::Value::UInt)src_len;
            meta_subs["dstlen"] = (Json::Value::Int)len_io->data_len;
            meta_subs["coder"] = len_io->meta;
            meta_subs["sname"] = "baselen";
            meta_streams.append(meta_subs);
            outdata->current_len += len_io->data_len;
            total_srclen += src_len;
            total_dstlen += len_io->data_len;
        }
        meta_base["lenmin"] = (Json::Value::UInt)baselen_min;
        meta_base["lenmax"] = (Json::Value::UInt)baselen_max;

        meta_base["tot_srclen"] = (Json::Value::UInt)total_srclen;
        meta_base["tot_dstlen"] = (Json::Value::UInt)total_dstlen;
        meta_base["streams"] = meta_streams;
        meta["base"] = meta_base;

    } else { /* 不使用reference压缩 */

        int add_enter = !!enc_baselen;
        base_src_len = 0;
        for (i = 1; i < line; i+= 4)
        {
            e = indata->npos[i];
            s = indata->npos[i - 1] + 1;
            l = e - s + add_enter;
            base_src_len += l;
        }

        if (base_src_len <= FC_MIN_LEN || base_src_len >= FC_MAX_LEN) // 用bcm压缩
        {
            check_exit(false, ERR_INTERNEL,
                "base total len is invalid: %d, should be (%d, %d)", base_src_len, FC_MIN_LEN, FC_MAX_LEN); 
            std::shared_ptr<coder_io> base_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
            std::shared_ptr<coder_bwt_cm> base_cm = std::make_shared<coder_bwt_cm>(base_io.get());
            for (i = 1; i < line; i+= 4)
            {
                e = indata->npos[i];
                s = indata->npos[i - 1] + 1;
                l = e - s + add_enter;
                src_len += l;
                base_cm->encode_line(p + s, l);
            }
            base_cm->encode_flush();
            dst_len = base_io->data_len;
            meta_subs = base_io->meta;
        }
        else // fc效果基本都比bcm效果好， fc整块压缩
        {
            // 为了节约内存，这里将base临时存放在outdata的末端
            check_exit(outdata->get_remain() > base_src_len, 
                ERR_INTERNEL, "current need this error, but not");
            
            uint8_t *pbase_tmp;
            src_len = 0;
            pbase_tmp = (uint8_t*)(outdata->get_buffer()) + outdata->buffer_size - base_src_len;
            for (i = 1; i < line; i+= 4)
            {
                e = indata->npos[i];
                s = indata->npos[i - 1] + 1;
                l = e - s + add_enter;
                memcpy(pbase_tmp + src_len, p + s, l);
                src_len += l;
            }

            std::shared_ptr<coder_io> base_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain() - base_src_len);
            std::shared_ptr<coder_fc> base_cm = std::make_shared<coder_fc>(base_io.get());
            base_cm->encode_line(pbase_tmp, base_src_len);
            base_cm->encode_flush();
            dst_len = base_io->data_len;
            meta_subs = base_io->meta;

            check_exit(base_io->data_len <= (outdata->get_remain() - base_src_len),
                   ERR_INTERNEL, "lack of buffer: %d should <= %d", base_io->data_len, outdata->get_remain() - base_src_len);
        }

        total_srclen += base_src_len;
        total_dstlen += dst_len;
        outdata->current_len += dst_len;

        meta_base.clear();
        meta_base["lenmin"] = (Json::Value::UInt)baselen_min;
        meta_base["lenmax"] = (Json::Value::UInt)baselen_max;
        
        meta_base["coder"] = meta_subs;
        meta_base["tot_srclen"] = (Json::Value::UInt)total_srclen;
        meta_base["tot_dstlen"] = (Json::Value::UInt)total_dstlen;
        meta["base"] = meta_base;
    }
    return true;
}

/* 压缩comment */
bool actuator_fastq::compress_comment()
{
    Json::Value meta_comment, meta_subs, meta_streams;
    uint32_t i, line, s, e, src_len = 0, l;
    uint8_t *p = indata->buffer;
    std::shared_ptr<coder_io> c_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
    std::shared_ptr<coder_bwt_cm> c_cm = std::make_shared<coder_bwt_cm>(c_io.get());

    switch (ctype)
    {
    case CT_JUST_PLUS:
        meta_comment["type"] = "plus";
        break;
    case CT_SAME_AS_ID:
        meta_comment["type"] = "same";
        break;
    case CT_OTHER:
        meta_comment["type"] = "other";
        line = indata->npos.size();
        for (i = 2; i < line; i += 4)
        {
            e = indata->npos[i];
            s = indata->npos[i - 1] + 1;
            l = e - s;
            c_cm->encode_line(p + s, l + 1); /* encode '\n'*/

            src_len += l;
        }
        c_cm->encode_flush();
        meta_subs.clear();
        meta_subs["srclen"] = (Json::Value::UInt)src_len;
        meta_subs["dstlen"] = (Json::Value::Int)c_io->data_len;
        meta_subs["coder"] = c_io->meta;
        meta_subs["sname"] = "comment";
        meta_streams.append(meta_subs);
        outdata->current_len += c_io->data_len;
        meta_comment["streams"] = meta_streams;
        break;
    default:
        break;
    }
    meta["comment"] = meta_comment;
    return true;
}

/* 压缩质量行 */
bool actuator_fastq::compress_quality()
{
    Json::Value meta_qual, meta_subs, meta_streams;
    uint32_t i, line, s, e, src_len = 0, l, flen = (qual_freq_table.size()) << 1;
    uint32_t total_srclen, total_dstlen;
    uint16_t qual_freq_arr[flen];
    uint8_t *p = indata->buffer;
    std::shared_ptr<coder_io> qual_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
    std::shared_ptr<coder_qual> qual_cm = std::make_shared<coder_qual>(qual_io.get(), true, qual_freq_table);

    total_srclen = total_dstlen = 0;
    line = indata->npos.size();
    for (i = 3; i < line; i += 4)
    {
        e = indata->npos[i];
        s = indata->npos[i - 1] + 1;
        l = e - s;
        qual_cm->encode_qual_gen2(p + indata->npos[i - 3] + 1, p + s, l);
        src_len += l;
    }
    qual_cm->encode_flush();

    /* 第一条流 */
    meta_subs.clear();
    meta_subs["srclen"] = (Json::Value::UInt)src_len;
    meta_subs["dstlen"] = (Json::Value::Int)qual_io->data_len;
    meta_subs["coder"] = qual_io->meta;
    meta_subs["sname"] = "qual";
    meta_streams.append(meta_subs);
    outdata->current_len += qual_io->data_len;
    total_srclen += src_len;
    total_dstlen += qual_io->data_len;

    /* 第二条流, 质量数符号频率表 */
    std::shared_ptr<coder_io> qual_freq_io = std::make_shared<coder_io>(outdata->get_curr(), outdata->get_remain());
    std::shared_ptr<coder_bwt_cm> qual_freq_cm = std::make_shared<coder_bwt_cm>(qual_freq_io.get());
    for (i = 0; i < qual_freq_table.size(); i++) {
        l = (i << 1);
        qual_freq_arr[l] = qual_freq_table[i].first;
        qual_freq_arr[l + 1] = qual_freq_table[i].second;
    }
    src_len = flen * sizeof(uint16_t);
    qual_freq_cm->encode_line((uint8_t *)qual_freq_arr, src_len);
    qual_freq_cm->encode_flush();
    meta_subs.clear();
    meta_subs["srclen"] = (Json::Value::UInt)src_len;
    meta_subs["dstlen"] = (Json::Value::Int)qual_freq_io->data_len;
    meta_subs["coder"] = qual_freq_io->meta;
    meta_subs["sname"] = "qual_freq";
    meta_streams.append(meta_subs);
    outdata->current_len += qual_freq_io->data_len;
    total_srclen += src_len;
    total_dstlen += qual_freq_io->data_len;

    meta_qual["tot_srclen"] = (Json::Value::UInt)total_srclen;
    meta_qual["tot_dstlen"] = (Json::Value::UInt)total_dstlen;
    meta_qual["streams"] = meta_streams;
    meta["quality"] = meta_qual;
    return true;
}



/* 压缩初始化 */
void actuator_fastq::compress_initialize()
{
    uint8_t *p;
    uint32_t n, line4 = (indata->npos.size() >> 2), baselen_len;
    const uint32_t lmax = indata->max_line_len + 4; /* 4个预留给碱基key不以4对齐的情况 */
    const uint32_t lsquash = (lmax >> 2) + !!(lmax & 0x3); /* squash 长度 */

    this->is_gen2 = (indata->btype == FASTQ_GEN2) || (indata->btype == FASTQ_GEN2_GZIP);
    len_base_mapped = (is_gen2) ? lmax : (lmax << 1);
    baselen_len = ((baselen_min == baselen_max) ? 0 : (is_gen2 ? (line4 << 1) : (line4 << 2)));

    if (refgene) {
      /* base pair + 4 base pair squash + 4 base squash + base mapped + base N pos in block +
        * base delete N + mapped pos + mapped pair + baselen each line
        */
        n = lmax + (lsquash << 3) + len_base_mapped + (base_n_cnt << 2);
        n += lmax + (line4 << 3) + line4 + baselen_len;
        safe_alloc(n, uint8_t, somebuffer);
        // fprintf(stderr, "\n--- [some buffer] ---len %ld\n", n);
        p = somebuffer;
        somebuffer_basepair = p;
        p += lmax;
        for (n = 0; n < 4; n++) {
            somebuffer_basepair_squash[n] = p;
            p += lsquash;
            somebuffer_base_squash[n] = p;
            p += lsquash;
        }

        somebuffer_base_mapped = p;
        p += len_base_mapped;

        somebuffer_base_npos = (uint32_t *)p;
        p += (base_n_cnt << 2);

        somebuffer_base_stripn = p;
        p += lmax;

        somebuffer_base_mpos = (uint64_t *)p;
        p += (line4 << 3);

        somebuffer_base_mpair = p;
        p += line4;

        if (baselen_len) {
            if (is_gen2)
                somebuffer_base_len2 = (uint16_t *)p;
            else
                somebuffer_base_len3 = (uint32_t *)p;
            p += baselen_len;
        }

        mapping = (is_gen2) ? (&actuator_fastq::mapping_gen2) : (&actuator_fastq::mapping_gen3);
    }
    
}

/* 压缩fastq */
bool actuator_fastq::compress()
{
    bool res;
    int64_t len;
    std::string md5;

    compress_initialize();

    check_exit(is_gen2, ERR_INTERNEL, "currently, not support gen3 compress.");

    /* 压缩 id */
    res = compress_id();
    check(res, false, "id compress failed");
    meta["id_lines"] = (Json::Value::Int)(indata->npos.size() >> 2);

    /* 压缩base */
    res = compress_base();
    check(res, false, "base compress failed");

    /* 压缩comment */
    res = compress_comment();
    check(res, false, "comment compress failed");

    /* 压缩quality*/
    res = compress_quality();
    check(res, false, "quality compress failed");

    /* 计算数据块的md5 */
    calc_md5(md5, indata->get_buffer(), indata->current_len);
    meta["md5"] = md5;

    /*  压缩block meta */
    coder_json cmeta;
    len = cmeta.encoder(meta, outdata->get_curr(), outdata->get_remain());
    check(len > 0, false, "block %ld meta info compress failed, get len: %ld", outdata->block_id, len);
    // std::cout << meta << std::endl;
    /*  设置meta压缩后的长度 */
    outdata->meta_encoded_len = len;
    outdata->current_len += len;
    // fprintf(stderr, "\n\n------------------[block meta]---------------- encoded len %lu\n",  len);

    if (outdata->current_len >= indata->current_len)
        check((indata->current_len * 2) <= outdata->buffer_size, ERR_INTERNEL,
              "encoder abnormal: block id %ld, type %s, encoded len %ld > origin len %ld\n",
              indata->block_id, get_typename(indata->btype).get(), outdata->current_len, indata->current_len);

#ifdef ANALYZ_STREAMS
    fprintf(stderr, "\n\n\tmapping_cnt [ %lld ] actg_matched [ %llu ] actg_total [ %llu ]\n",
            mapping_cnt, actg_matched, actg_total);
#endif
    return true;
}

/* 解压fastq */
bool actuator_fastq::decompress()
{
    uint8_t *pout, *pdata, *pid, *pseq;
    uint32_t base_maxlen, base_minlen, base_ncount, total_baselen = 0;
    int64_t enclen_meta, enclen_data, len_id;
    int32_t lines, n, l, len, len_actual, o, len_stripn, out_base_len, out_quality_len;
    int32_t base_lines = 0, ncnt_curr_line, npos_offset = 0;
    coder_json cmeta;
    std::string id_splitsym, md5;
    pbgz_stream_header header_meta(indata->buffer);
    std::vector<coder *> id_decoders;
    coder *base_decoder, *com_decoder;
    coder_qual *qual_decoder;
    const uint8_t actg4[4] = {'A', 'C', 'T', 'G'};

    /* 首先解析出当前block的meta信息，获取对应行的流信息和编码器等信息 */
    cmeta.decoder(indata->buffer + header_meta.get_bufferlen(), header_meta.get_datalen(), meta);

    base_maxlen = meta["base"]["lenmax"].asUInt();
    base_minlen = meta["base"]["lenmin"].asUInt();
    base_ncount = meta["base"]["ncount"].asUInt();
    outdata->block_id = meta["block_id"].asInt64();
    id_splitsym = meta["id"]["split_sym"].asString();
    lines = meta["id_lines"].asInt();

    initialize_decode_id(id_decoders);
    initialize_decode_base(base_decoder);
    initialize_decode_comment(com_decoder);
    initialize_decode_quality(qual_decoder);

    uint8_t *pbase_out = nullptr;
    const uint8_t *pend = outdata->buffer + outdata->buffer_size;
    if (meta["base"]["coder"]["magic"].asString() == "coder_fc") // 指向整个base块
        pbase_out = (uint8_t*)(outdata->get_buffer()) + outdata->buffer_size - meta["base"]["tot_srclen"].asUInt();

    for (l = 0; l < lines; l++)
    {
        /* decode id */
        pid = outdata->get_curr();
        for (n = 0; n < id_decoders.size(); n++) {
            len = id_decoders[n]->decode_line(outdata->get_curr(), outdata->get_remain(), id_splitsym[n]);
            outdata->current_len += len;
        }
        len_id = outdata->get_curr() - pid;

        /* decode base */
        pseq = outdata->get_curr();
        if (refgene) {
            /* decode base mapping stream with strip N */
            len_actual = (somebuffer_base_len2) ? (somebuffer_base_len2[base_lines] + base_minlen) : base_maxlen;
            pout = outdata->get_curr();

            if (base_ncount && npos_offset < base_ncount) { /* 该block有N，且未处理完所有N */
                /* 计算当前base行N的个数 */
                n = total_baselen+ len_actual;
                for (o = npos_offset; o < base_ncount; o++) {
                    if (somebuffer_base_npos[o] + 1 > n)
                        break;
                }
                ncnt_curr_line = o - npos_offset; 

                /* 解压当前base行的mapping流 */
                len_stripn = len_actual - ncnt_curr_line;
                len = base_decoder->decode_line(somebuffer_base_stripn, len_stripn, UINT8_MAX, false);
                check_exit(len_stripn == len, ERR_INTERNEL,
                           "base decode failed in block %llu, expect len %u, actual %u", meta["block_id"].asInt64(), len_stripn, len);

                if (2 == somebuffer_base_mpair[base_lines]) {
                    for (o = 0; o < len_stripn; o++)
                        somebuffer_base_stripn[o] = actg4[somebuffer_base_stripn[o]];
                    pdata = somebuffer_base_stripn;
                } else {
                    /* 得到对应位置的reference */
                    refgene->get_stretch_2bits1char(somebuffer_refe_stretch, len_stripn, somebuffer_base_mpos[base_lines]);

                    /* 还原mapping的base squash之后的流 */
                    actg_xor(somebuffer_base_stripn, somebuffer_refe_stretch, somebuffer_base_stripn, len_stripn);

                    refgene->get_actg_from2bits(somebuffer_base_stripn, len_stripn, somebuffer_refe_stretch);
                    if (somebuffer_base_mpair[base_lines]) {
                        actg_pair(somebuffer_base_stripn, somebuffer_refe_stretch, len_stripn);
                        pdata = somebuffer_base_stripn;
                    } else
                        pdata = somebuffer_refe_stretch;
                }

                /* 将还原数据复制到输出buffer */
                for (o = 0, n = 0; n < len_actual; n++) {
                    if (npos_offset < base_ncount && (total_baselen + n) == somebuffer_base_npos[npos_offset]) { /* 当前位置为N */
                        *pout++ = 'N' ;
                        npos_offset++;
                    } else /* 当前位置不为N */
                        *pout++ = (pdata[o++]) ;
                }

            } else { /* 块中没有N */
                /* 解压当前base行的mapping流 */
                len_stripn = len_actual;
                len = base_decoder->decode_line(somebuffer_base_stripn, len_stripn, UINT8_MAX, false);
                check_exit(len_stripn == len, ERR_INTERNEL,
                           "base decode failed in block %llu, expect len %u, actual %u", meta["block_id"].asInt64(), len_stripn, len);

                 if (2 == somebuffer_base_mpair[base_lines]) {
                    for (o = 0; o < len_stripn; o++)
                        pout[o] = actg4[somebuffer_base_stripn[o]];
                } else {
                    /* 得到对应位置的reference */
                    refgene->get_stretch_2bits1char(somebuffer_refe_stretch, len_stripn, somebuffer_base_mpos[base_lines]);

                    /* 还原mapping的base squash之后的流 */
                    actg_xor(somebuffer_base_stripn, somebuffer_refe_stretch, somebuffer_base_stripn, len_stripn);

                    if (somebuffer_base_mpair[base_lines]) {
                        refgene->get_actg_from2bits(somebuffer_base_stripn, len_stripn, somebuffer_refe_stretch);
                        actg_pair(pout, somebuffer_refe_stretch, len_stripn);
                    } else
                        refgene->get_actg_from2bits(somebuffer_base_stripn, len_stripn, pout);
                }
            }
            total_baselen += len_actual;
            outdata->current_len += len_actual;

            *(outdata->get_curr()) = '\n'; outdata->current_len++;
        } else {
            if (base_maxlen == base_minlen)
            { 
                len_actual = base_maxlen;
                if (meta["base"]["coder"]["magic"].asString() == "coder_fc")
                {
                    /* 检查是否还有pbase_out数据 ，并且outdata buffer是否与pbase_out覆盖了 */
                    check_exit(pbase_out + len_actual <= pend && outdata->buffer < pbase_out,
                               ERR_INTERNEL, "decode failed !");
                    memcpy(pseq, pbase_out, len_actual);
                    pbase_out += len_actual;
                }
                else if (meta["base"]["coder"]["magic"].asString() == "coder_bwt_cm")
                    base_decoder->decode_line(pseq, len_actual, UINT8_MAX, false);
                else
                    check_exit(false, ERR_INTERNEL, "undefined decoder type");

                outdata->current_len += len_actual;
                *(outdata->get_curr()) = '\n'; outdata->current_len++; // 相等时没有压缩回车

            } else {
                if (meta["base"]["coder"]["magic"].asString() == "coder_fc")
                {
                    uint8_t *p, *p1 = pseq;
                    for (p = pbase_out; p < pend; p++)
                    {
                        *p1++ = *p;
                        if (*p == '\n')break;
                    }
                    len_actual = p - pbase_out + 1;
                    /* 检查是否还有pbase_out数据 ，并且outdata buffer是否与pbase_out覆盖了 */
                    check_exit(pbase_out + len_actual <= pend && outdata->buffer < pbase_out,
                               ERR_INTERNEL, "decode failed !");
                    pbase_out += len_actual;
                }
                else if (meta["base"]["coder"]["magic"].asString() == "coder_bwt_cm")
                    len_actual = base_decoder->decode_line(pseq, base_maxlen, '\n', false);
                else
                    check_exit(false, ERR_INTERNEL, "undefined decoder type");

                outdata->current_len += len_actual;
                len_actual -= 1;
            }
        }

        /* decode comment */
        pout = outdata->get_curr();
        switch (ctype)
        {
        case CT_JUST_PLUS:
            *pout++ = '+';
            *pout = '\n';
            outdata->current_len += 2;
            break;
        case CT_SAME_AS_ID:
            memcpy(pout, pid, len_id); /*  have included '\n' */
            outdata->current_len += len_id;
            break;
        case CT_OTHER:
            n = com_decoder->decode_line(pout, 0, '\n', false);
            outdata->current_len += n;
            break;
        default:
            break;
        }

        /* decode quality */
        qual_decoder->decode_qual_gen2(pseq, outdata->get_curr(), len_actual);
        outdata->current_len += len_actual;
        *(outdata->get_curr()) = '\n'; outdata->current_len++;

        base_lines++;
    }

    if (refgene)
        out_base_len = meta["base"]["streams"][0]["srclen"].asUInt() + base_ncount;
    else
        out_base_len = meta["base"]["tot_srclen"].asUInt() - ((base_maxlen == base_minlen) ? (0) : (lines));
    out_quality_len = meta["quality"]["streams"][0]["srclen"].asUInt();
    if (out_base_len == out_quality_len + 1)
        outdata->current_len -= 1; // strip n
    else
        check_exit(out_base_len == out_quality_len,
                   ERR_INTERNEL, "check falied: base out len %u != quality out len %u", out_base_len, out_quality_len);

    /* 检查md5 */
    calc_md5(md5, outdata->get_buffer(), outdata->current_len);
    check_exit(md5 == meta["md5"].asString(), ERR_INTERNEL, "md5 check failed in block %ld", outdata->block_id);
    for(auto &c: id_decoders)
        delete c;
    delete base_decoder;
    if (com_decoder)
        delete com_decoder;
    delete qual_decoder;

    return true;
}