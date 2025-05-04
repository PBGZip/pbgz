#ifndef _CODER_BASE_H_
#define _CODER_BASE_H_

#include "coder.h"

template <int32_t bits_valid_cnt, int32_t model_order, int32_t rate>
class coder_base
{
    template <int32_t rate_bit_counter>
    class bit_counter
    {
    public:
        bit_counter(void) : p(1 << 15) {}
        inline void update_bit0(void)
        {
            p -= p >> rate_bit_counter;
        }
        inline void update_bit1(void)
        {
            p += (p ^ 0xFFFF) >> rate_bit_counter;
        }
        uint16_t p;
    };

    class bit_coder
    {
    public:
        bit_coder(coder_io *io)
        {
            low = 0;
            high = UINT32_MAX;
            code = 0;
            this->io = io;
            flushed = false;
        }

        virtual ~bit_coder()
        {
            if (!flushed && io->m == coder_io::MENC)
                encode_flush();
        }

        void encode_flush()
        {
            int32_t i;
            if (flushed)
                return;
            for (int32_t i = 0; i < 4; ++i)
            {
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
            }
            flushed = true;
        }

        inline void encode_bit(int32_t bit, uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG); /* 这里是把高概率符号放前面了，低概率放后面，P表示高概率符号的概率 */

            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* 归一化 */
            while ((low ^ high) < (1 << 24))
            {
                *(io->data + io->data_len++) = (low >> 24);
                low <<= 8;
                high = (high << 8) + 255;
            }
        }

        inline int32_t decode_bit(uint32_t p)
        {
            const uint32_t mid = low + ((uint64_t(high - low) * p) >> P_LOG);

            const int32_t bit = (code <= mid);
            if (bit)
                high = mid;
            else
                low = mid + 1;

            /* 归一化 */
            while ((low ^ high) < (1 << 24))
            {
                low <<= 8;
                high = (high << 8) + 255;
                code = (code << 8) + *(io->data + io->data_len++);
            }
            return bit;
        }

        void decode_init()
        {
            int32_t i;
            for (i = 0; i < 4; ++i)
                code = (code << 8) + *(io->data + io->data_len++);
        }

    private:
        uint32_t low;
        uint32_t high;
        uint32_t code;
        bool flushed;
        coder_io *io;
        static const int16_t P_LOG = 16;
    };

    template <int32_t rate_char_coder>
    class char_coder
    {
    public:
        char_coder(void)
        {
        }

        inline void encode_8bits(int32_t c, bit_coder *coder)
        {
            int32_t ctx = 1;
            while (this->context_cnt > ctx)
            {
                const int32_t bit = c & (this->context_cnt >> 1);
                const int32_t p = this->counter[ctx].p;
                c += c;
                if (bit)
                {
                    coder->encode_bit(1, p);
                    this->counter[ctx].update_bit1();
                    ctx += ctx + 1;
                }
                else
                {
                    coder->encode_bit(0, p);
                    this->counter[ctx].update_bit0();
                    ctx += ctx;
                }
            }
        }

        inline int32_t decode_8bits(bit_coder *coder)
        {
            int32_t ctx = 1;
            while (this->context_cnt > ctx)
            {
                const int32_t p = this->counter[ctx].p;
                const int32_t bit = coder->decoder_bit(p);

                if (bit)
                {
                    this->counter[ctx].update_bit1();
                    ctx += ctx + 1;
                }
                else
                {
                    this->counter[ctx].update_bit0();
                    ctx += ctx;
                }
            }
            return ctx & context_mask;
        }

    private:
        static const int32_t context_cnt = (1 << bits_valid_cnt);
        static const int32_t context_mask = (context_cnt - 1);

        bit_counter<rate_char_coder> counter[context_cnt];
    };

    template <int32_t rate_model, int32_t order_model>
    class model
    {
    public:
        model() : hash(0)
        {
            this->coders = new char_coder<rate_model>[this->hash_mask + 1];
        }

        ~model()
        {
            delete[] this->coders;
        }

        void encode_line(const uint8_t *in, int32_t in_len, bit_coder *coder)
        {
            int32_t n, symbol;
            this->hash = 0;
            for (n = 0; n < in_len; n++)
            {
                symbol = in[n];
                this->coders[this->get_hash()].encode_8bits(symbol, coder);
                this->update_hash(symbol);
            }
        }

        int32_t decode_line(uint8_t *out, int32_t out_len, bit_coder *coder)
        {
            int32_t n, symbol;
            this->hash = 0;
            for (n = 0; n < out_len; n++)
            {
                symbol = this->coders[this->get_hash()].decode_8bits(coder);
                this->update_hash(symbol);
                out[n] = (uint8_t)symbol;
            }
            return out_len;
        }

    private:
        void update_hash(int32_t symbol)
        {
            this->hash <<= this->valid_bits_cnt;
            this->hash |= (uint32_t)symbol;
            this->hash &= this->hash_mask;
        }

        int64_t get_hash()
        {
            return this->hash;
        }

    private:
        static const int32_t valid_bits_cnt = bits_valid_cnt;
        static const int32_t order = order_model;
        static const uint64_t hash_mask = ((1ul << (order * valid_bits_cnt)) - 1ul);

    private:
        uint64_t hash;
        char_coder<rate_model> *coders;
    };

public:
    coder_base(coder_io *io)
    {
        this->io = io;
        this->io->appen_magic("coder_base");

        coder = new bit_coder(io);
        check_exit(coder, ERR_MEM_NOENOUGH, "coder initialize failed: no enough memory");
        check_exit(io->m != coder_io::MUNSET, ERR_INTERNEL,"unset code mode");
    }
    ~coder_base()
    {
        if (coder)
            delete coder;
    }

    void encode_line(const uint8_t *in, int32_t in_len)
    {
        return m.encode_line(in, in_len, coder);
    }

    int32_t decode_line(uint8_t *out, int32_t out_len)
    {
        return m.decode_line(out, out_len);
    }

    void encode_flush(){
        io->meta["bits_valid"] = (Json::Value::Int)bits_valid_cnt;
        io->meta["order"] = (Json::Value::Int)model_order;
        io->meta["rate"] = (Json::Value::Int)rate;

        this->coder->encode_flush();
    }

private:
    coder_io *io;
    bit_coder *coder;
    model<model_order, rate> m;
};
#endif