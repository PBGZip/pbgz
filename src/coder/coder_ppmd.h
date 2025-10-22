#ifndef _CODER_PPMD_H_
#define _CODER_PPMD_H_

#include "coder_io.h"
#include "coder.h"
#include "ppmd/Ppmd8.h"

class coder_ppmd : public coder
{
public:
    coder_ppmd(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_ppmd");

        flushed = false;
    }

    ~coder_ppmd()
    {
        if (io->m != coder_io::MENC && !flushed)
            encode_flush();
        Ppmd8_Free(&ppmd, &ialloc);
    }

    int64_t encode(const uint8_t *src, int64_t src_len)
    {
        int64_t n;
        if (io->m != coder_io::MENC)
        {
            opt_mem = 8;
            opt_order = 6;
            opt_restore = 0;
            info = (opt_order - 1) | ((opt_mem - 1) << 4) | (('I' - 'A') << 12);
            fnlen = 1;

            cw = {Write, io};
            ppmd.Stream.Out = reinterpret_cast<IByteOut *>(&cw);
            Ppmd8_Construct(&ppmd);
            Ppmd8_Alloc(&ppmd, opt_mem << 20, &ialloc);
            Ppmd8_RangeEnc_Init(&ppmd);
            Ppmd8_Init(&ppmd, this->opt_order, 0);

            io->m = coder_io::MENC;
        }
        for (n = 0; n < src_len; n++)
            Ppmd8_EncodeSymbol(&ppmd, src[n]);
        return src_len;
    }

    int64_t decode(uint8_t *out, int64_t out_len)
    {
        int32_t ch;
        int64_t len = 0;
        if (io->m != coder_io::MDEC)
        {
            info = io->meta["info"].asUInt();
            fnlen = io->meta["fnlen"].asUInt();

            opt_restore = fnlen >> 14;
            opt_order = (info & 0xf) + 1;
            opt_mem = ((info >> 4) & 0xff) + 1;

            struct CharReader cr = {Read, io};
            ppmd.Stream.In = (IByteIn *)&cr;
            Ppmd8_Construct(&ppmd);
            Ppmd8_Alloc(&ppmd, opt_mem << 20, &ialloc);
            Ppmd8_RangeDec_Init(&ppmd);
            Ppmd8_Init(&ppmd, opt_order, opt_restore);

            io->m = coder_io::MDEC;
        }

        for (;;)
        {
            ch = Ppmd8_DecodeSymbol(&ppmd);
            if (ch < 0)
                break;
            *(out + len++) = ch;
            if (len == out_len)
                return len;
        }
        return len;
    }

    /* 编码完成后需要调用它将缓存中的数据压缩 */
    void encode_flush()
    {
        if (io->m != coder_io::MENC || flushed)
            return;
        Ppmd8_EncodeSymbol(&ppmd, -1); /* EndMark */
        Ppmd8_RangeEnc_FlushData(&ppmd);
        io->meta["info"] = (Json::Value::UInt)info;
        io->meta["fnlen"] = (Json::Value::UInt)fnlen;
        flushed = true;
    }

private:
    struct CharWriter
    {
        /* Inherits from IByteOut */
        void (*Write)(void *p, Byte b);
        coder_io *out;
    };

    static void Write(void *p, Byte b)
    {
        struct CharWriter *cw = (struct CharWriter *)p;
        cw->out->putc(b);
    }

    struct CharReader
    {
        /* Inherits from IByteIn */
        Byte (*Read)(void *p);
        coder_io *in;
    };


    static Byte Read(void *p)
    {
        struct CharReader *cr = (struct CharReader *)p;
        
        int c = cr->in->getc();
        return (c == '\0') ? 0 : c;
    }

    static void *pmalloc(void *p, size_t size)
    {
        (void)p;
        return malloc(size);
    }

    static void pfree(void *p, void *addr)
    {
        (void)p;
        free(addr);
    }

    ISzAlloc ialloc = {pmalloc, pfree};

private:
    struct CharWriter cw;

    CPpmd8 ppmd;
    uint16_t info;
    uint16_t fnlen;

    int32_t opt_mem;
    int32_t opt_order;
    int32_t opt_restore;

    bool flushed;
};

#endif