#ifndef _CODER_AFFIX_MATCH_H_
#define _CODER_AFFIX_MATCH_H_

#include "simple_model.h"
#include "../manager.h"
#include "coder_io.h"
#include "coder.h"

/* 基于前后缀匹配的字符串编码器 */
class coder_affix_match : public coder
{
public:
    coder_affix_match(coder_io *io)
    {
        int32_t n;
        safe_alloc(8192, SIMPLE_MODEL<128>, model_middle);
        for (n = 0; n < 8192; n++)
            model_middle[n].reset();
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_affix_match");
        this->plast = nullptr;

        last_len = last_prelen = last_suflen = 0;
        this->flushed = false;

        level = (io->meta["level"].isInt()) ? (io->meta["level"].asInt()) : 2;
        io->meta["level"] = (Json::Value::Int)level;
        check_exit(level <= 2 && level >= 1, ERR_INTERNEL, "coder level should in [1, 2], current is %d", level);
    }

    virtual ~coder_affix_match()
    {
        if (model_middle)
            free(model_middle);
        if (this->plast)
            free(this->plast);
        if (!flushed)
            if (io->m == coder_io::MENC)
                encode_flush();
    }

    /* 对外压缩接口, 如果调用encode后in会被释放，那么need2hold需要设置为true */
    void encode_line(const uint8_t *in, const int32_t in_len, bool need2hold = false)
    {
         int32_t pre_len, suf_len;
        int32_t i, j, k, last_ch;
        int32_t len2, match;

        if (io->m != coder_io::MENC)
        {
            rc.output((char *)(this->io->data));
            rc.StartEncode();
            io->m = coder_io::MENC;
            this->last_capacity = in_len;
            safe_alloc_init(this->last_capacity, uint8_t, last, ' ');
            this->plast = last;
        }

        /* 获取前缀长度 */
        for (i = 0; i < in_len && i < last_len; i++)
        {
            if (in[i] != last[i])
                break;
        }
        pre_len = i;

        /* 获取后缀长度 */
        for (i = in_len - 1, j = last_len - 1; i >= 0 && j >= 0; i--, j--)
        {
            if (in[i] != last[j])
                break;
        }
        suf_len = in_len - 1 - i;
        if (in_len - suf_len - pre_len < 0)
            suf_len = in_len - pre_len;

        if (encode_higher())
        { /* 更高压缩率压缩 */
            model_prefix[last_prelen].encodeSymbolOrder(&rc, pre_len);
            model_suffix[last_suflen].encodeSymbolOrder(&rc, suf_len);
            model_len[last_len].encodeSymbolOrder(&rc, in_len);
        }
        else
        {
            model_prefix[last_prelen].encodeSymbol(&rc, pre_len);
            model_suffix[last_suflen].encodeSymbol(&rc, suf_len);
            model_len[last_len].encodeSymbol(&rc, in_len);
        }

        last_prelen = pre_len;
        last_suflen = suf_len;
        len2 = in_len - suf_len;
        match = !!pre_len;
        for (i = j = pre_len, k = 0; i < len2; i++, j++, k++)
        {
            last_ch = (((last[j] - 32) << 1) + match + (k << 6)) & 0x1FFF;
            if (encode_higher()) /* 更高压缩率压缩 */
                model_middle[last_ch].encodeSymbolOrder(&rc, in[i] & 0x7f);
            else
                model_middle[last_ch].encodeSymbol(&rc, in[i] & 0x7f);

            if (in[i] == ' ' && last[j] != ' ') j++;
            if (in[i] != ' ' && in[j] == ' ') j--;
            if (in[i] == ':' && in[j] != ':') j++;
            if (in[i] != ':' && in[j] == ':') j--;
            if (in[i] == ':' || in[i] == ' ') k = (k + 3) >> 2 << 2;

            match = (in[i] == last[j]);
        }

        if (need2hold)
        {
            safe_realloc(this->last_capacity, uint8_t, this->last, in_len);
            memcpy(this->last, in, in_len);
            last_len = in_len;
        }
        else
        {
            last = (uint8_t *)in;
            last_len = in_len;
        }
    }

    /* 对外解压接口，如果调用decode后in会被释放，那么need2hold需要设置为true */
    int32_t decode_line(uint8_t *out, int32_t out_len, uint8_t split_ch = UINT8_MAX, bool need2hold = false)
    {
        uint8_t c;
        int32_t pre_len, suf_len, len;
        int32_t i, j, k;
        int32_t last_ch;
        int32_t len2, match;

        if (io->m != coder_io::MDEC)
        {
            rc.input((char *)(this->io->data));
            rc.StartDecode();
            io->m = coder_io::MDEC;
            this->last_capacity = 1024;
            safe_alloc_init(this->last_capacity, uint8_t, last, ' ');
            this->plast = last;
        }

        if (encode_higher()) {
            pre_len = model_prefix[last_prelen].decodeSymbolOrder(&rc);
            suf_len = model_suffix[last_suflen].decodeSymbolOrder(&rc);
            len = model_len[last_len].decodeSymbolOrder(&rc);
        } else {
            pre_len = model_prefix[last_prelen].decodeSymbol(&rc);
            suf_len = model_suffix[last_suflen].decodeSymbol(&rc);
            len = model_len[last_len].decodeSymbol(&rc);
        }

        last_prelen = pre_len;
        last_suflen = suf_len;

        for (i = 0; i < pre_len; i++)
            out[i] = last[i];

        len2 = len - suf_len, match = pre_len ? 1 : 0;
        for (i = j = pre_len, k = 0; i < len2; i++, j++, k++) {
            last_ch = (((last[j] - 32) << 1) + match + (k << 6)) & 0x1FFF;
            if (encode_higher())
                c = model_middle[last_ch].decodeSymbolOrder(&rc);
            else
                c = model_middle[last_ch].decodeSymbol(&rc);

            out[i] = c;

            if (c == ' ' && last[j] != ' ') j++;
            if (c != ' ' && last[j] == ' ') j--;
            if (c == ':' && last[j] != ':') j++;
            if (c != ':' && last[j] == ':') j--;
            if (out[i] == ':' || out[i] == ' ') k = (k + 3) >> 2 << 2;

            match = c == last[j];
        }

        for (j = last_len - suf_len; i < len; i++, j++)
            out[i] = last[j];

        if (need2hold)
        {
            safe_realloc(this->last_capacity, uint8_t, this->last, len);
            memcpy(this->last, out, len);
            last_len = len;
        }
        else
        {
            last = (uint8_t *)out;
            last_len = len;
        }

        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t *out, int32_t out_len, uint8_t *rely = nullptr, uint8_t split_ch = UINT8_MAX, bool need2hold = false)
    {
        return 0;
    }

    /* 对外接口，获取当前已经消耗数据的长度 */
    int32_t decode_inlen()
    {
        return rc.size_in();
    }

    /* 编码完未处理的数据 */
    void encode_flush()
    {
        if (flushed)
            return;
        if (io->m == coder_io::MENC)
        {
            rc.FinishEncode();
            io->data_len += rc.size_out();
        }
        flushed = true;
    }

private:
    /* 是否启用更高压缩 */
    const bool encode_higher() const
    {
        return this->level > 1;
    }

private:
    RangeCoder rc;
    uint8_t *last, *plast; /* 上一个压缩的字符串 */
    int32_t last_capacity;
    int32_t last_len;    /* 上一个压缩的的长度 */
    int32_t last_prelen; /*上一次匹配时相同的前缀长度 */
    int32_t last_suflen; /*上一次匹配时相同的后缀长度 */

    SIMPLE_MODEL<256> model_prefix[256];
    SIMPLE_MODEL<256> model_suffix[256];
    SIMPLE_MODEL<256> model_len[256];
    SIMPLE_MODEL<128> *model_middle;
    int32_t level;
    bool flushed;
};

#endif