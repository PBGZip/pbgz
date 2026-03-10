/*
 * coder_affix_match.h - Header file for pbgz project
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _CODER_AFFIX_MATCH_H_
#define _CODER_AFFIX_MATCH_H_

#include <stdint.h>
#include <cinttypes> 

#include "simple_model.h"
#include "coder_io.h"
#include "coder.h"


/* String encoder based on prefix and suffix matching */
class coder_affix_match : public coder
{
public:
    coder_affix_match(coder_io *io)
    {
        int32_t n;
        // safe_alloc(8192, SIMPLE_MODEL<128>, model_middle);
        void* ptr = safe_alloc(8192 * sizeof(SIMPLE_MODEL<128>));
        model_middle = static_cast<SIMPLE_MODEL<128>*>(ptr);
        check_exit((model_middle != nullptr),  coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Error: Insufficient memory: need %" PRIu64 " MB", 8192 * sizeof(SIMPLE_MODEL<128>) >> 20);

        for (n = 0; n < 8192; n++)
            model_middle[n].reset();
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_affix_match");
        this->plast = nullptr;
        this->last = nullptr;

        last_len = last_prelen = last_suflen = 0;
        this->flushed = false;

        level = (io->meta["level"].isInt()) ? (io->meta["level"].asInt()) : 2;
        io->meta["level"] = (Json::Value::Int)level;
        check_exit(level <= 2 && level >= 1,  coder_ns::CODER_ERR_BAD_ARGS, "coder level should in [1, 2], current is %d", level);
    }

    virtual ~coder_affix_match()
    {
        if (model_middle)
            safe_free((void**)&model_middle);
        if (this->plast)
            safe_free((void**)&(this->plast));
        if (!flushed)
            if (io->m == coder_io::MENC)
                encode_flush();
    }

    /* External compression interface, set need2hold to true if input will be freed after encode call */
    void encode_line(const uint8_t *in, const uint32_t in_len, bool need2hold = false)
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
            this->last = static_cast<uint8_t*>(safe_alloc_init(this->last_capacity, ' '));
            this->last_len = 0;
            this->plast = last;
        }
        /* Get prefix length */
        for (i = 0; i < (int32_t)in_len && i < last_len; i++)
        {
            if (in[i] != last[i])
                break;
        }
        pre_len = i;

        /* Get suffix length */
        for (i = in_len - 1, j = last_len - 1; i >= 0 && j >= 0; i--, j--)
        {
            if (in[i] != last[j])
                break;
        }
        suf_len = in_len - 1 - i;
        if ((int32_t)in_len < suf_len + pre_len)
            suf_len = in_len - pre_len;

        if (encode_higher())
        { /* Higher compression rate */
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
            uint8_t ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            last_ch = (((ch - 32) << 1) + match + (k << 6)) & 0x1FFF;
            if (encode_higher()) /* Higher compression rate */
                model_middle[last_ch].encodeSymbolOrder(&rc, in[i] & 0x7f);
            else
                model_middle[last_ch].encodeSymbol(&rc, in[i] & 0x7f);

            ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            if (in[i] == ' ' && ch != ' ') j++;

            ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            if (in[i] != ' ' && ch == ' ') j--;

            ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            if (in[i] == ':' && ch != ':') j++;

            ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            if (in[i] != ':' && ch == ':') j--;

            if (in[i] == ':' || in[i] == ' ') k = (k + 3) >> 2 << 2;

            ch = (j < (int32_t)last_len && j >= 0) ? last[j] : ' ';
            match = (in[i] == ch);
        }

        if (need2hold)
        {
            /* Allocate new memory to avoid in-place reallocation issues */
            if (in_len > this->last_capacity) 
            {
                uint8_t* new_last = static_cast<uint8_t*>(safe_alloc(in_len));
                if (new_last) {
                    memcpy(new_last, in, in_len);
                    if (last != in && this->last) 
                    {
                        safe_free((void**)&this->last);
                    }
                    this->last = new_last;
                    this->last_len = in_len;
                    this->last_capacity = in_len;
                    this->plast = last;
                }
            } 
            else 
            {
                memcpy(this->last, in, in_len);
                this->last_len = in_len;
            }
        }
        else
        {
            last = (uint8_t *)in;
            last_len = in_len;
        }
    }

    /* External decompression interface, set need2hold to true if input will be freed after decode call */
    int32_t decode_line(uint8_t *out, uint32_t __attribute__((unused)) out_len,
                        uint8_t __attribute__((unused)) split_ch = UINT8_MAX, 
                        bool __attribute__((unused)) need2hold = false)
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
            last = static_cast<uint8_t*>(safe_alloc_init(this->last_capacity, ' '));
            //check_exit(last, coder_ns::CODER_ERR_MEM_ALLOC_FAIL, "Error: Insufficient memory: need %" PRIu64 " MB\n",  static_cast<uint64_t>(this->last_capacity) >> 20);
            if (!last) 
            {
                coder_logger(coder_ns::ERROR, "Error: Insufficient memory: need %" PRIu64 " MB\n",  
                    static_cast<uint64_t>(this->last_capacity) >> 20);
                return coder_ns::CODER_ERR_MEM_ALLOC_FAIL;
            }
            this->last_len = 0;
            this->plast = last;
        }

        if (encode_higher()) 
        {
            pre_len = model_prefix[last_prelen].decodeSymbolOrder(&rc);
            suf_len = model_suffix[last_suflen].decodeSymbolOrder(&rc);
            len = model_len[last_len].decodeSymbolOrder(&rc);
        } 
        else 
        {
            pre_len = model_prefix[last_prelen].decodeSymbol(&rc);
            suf_len = model_suffix[last_suflen].decodeSymbol(&rc);
            len = model_len[last_len].decodeSymbol(&rc);
        }

        last_prelen = pre_len;
        last_suflen = suf_len;

        for (i = 0; i < pre_len; i++)
            out[i] = last[i];

        len2 = len - suf_len, match = pre_len ? 1 : 0;
        for (i = j = pre_len, k = 0; i < len2; i++, j++, k++) 
        {
            uint8_t ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            last_ch = (((ch - 32) << 1) + match + (k << 6)) & 0x1FFF;
            if (encode_higher())
                c = model_middle[last_ch].decodeSymbolOrder(&rc);
            else
                c = model_middle[last_ch].decodeSymbol(&rc);

            out[i] = c;

            ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            if (c == ' ' && ch != ' ') j++;
            
            ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            if (c != ' ' && ch == ' ') j--;
            
            ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            if (c == ':' && ch != ':') j++;
            
            ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            if (c != ':' && ch == ':') j--;
            
            if (out[i] == ':' || out[i] == ' ') k = (k + 3) >> 2 << 2;
            
            ch = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
            match = c == ch;
        }

        for (j = last_len - suf_len; i < len; i++, j++) 
        {
            out[i] = (j < (int32_t)this->last_len && j >= 0) ? last[j] : ' ';
        }

        if (need2hold)
        {
            if (len > (int32_t)this->last_capacity) {
                uint8_t* new_last = static_cast<uint8_t*>(safe_alloc(len));
                if (new_last) 
                {
                    memcpy(new_last, out, len);
                    if (this->last) 
                    {
                        safe_free((void**)&this->last);
                    }
                    this->last = new_last;
                    this->last_len = len;
                    this->last_capacity = len;
                    this->plast = last;
                }
            } 
            else 
            {
                memcpy(this->last, out, len);
                last_len = len;
            }
        }
        else
        {
            last = (uint8_t *)out;
            last_len = len;
        }

        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool)
    {
        return 0;
    }

    /* External interface to get the length of currently consumed data */
    int32_t decode_inlen()
    {
        return rc.size_in();
    }

    /* Flush unprocessed encoded data */
    void encode_flush()
    {
        if (flushed) {
            return;
        }
           
        if (io->m == coder_io::MENC)
        {
            rc.FinishEncode();
            io->data_len += rc.size_out();
        }

        flushed = true;
    }

private:
    /* Whether to enable higher compression */
    bool encode_higher() const
    {
        return this->level > 1;
    }

private:
    RangeCoder rc;
    uint8_t *last, *plast; /* Previously compressed string */
    uint32_t last_capacity;
    int32_t last_len;    /* Length of previous compression */
    int32_t last_prelen; /* Length of same prefix in last match */
    int32_t last_suflen; /* Length of same suffix in last match */

    SIMPLE_MODEL<256> model_prefix[256];
    SIMPLE_MODEL<256> model_suffix[256];
    SIMPLE_MODEL<256> model_len[256];
    SIMPLE_MODEL<128> *model_middle;
    int32_t level;
    bool flushed;
};

#endif
