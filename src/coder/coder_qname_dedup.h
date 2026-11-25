/*
 * coder_qname_dedup.h - the repeat stream of the QNAME dedup layout.
 *
 * One payload per block: `lines` flag bytes (0 = first occurrence, 1 = repeat) followed
 * by three planes holding the repeat distances (24-bit little-endian, plane by plane).
 * Flags and each distance plane get their own adaptive model, which is what makes a
 * repeat cost about ten bits instead of a whole name. The shape travels in the stream
 * meta ("lines"/"reps"), so decoder and encoder agree without extra framing.
 *
 * Copyright (C) 2025 PBGZip
 */
#ifndef _CODER_QNAME_DEDUP_H_
#define _CODER_QNAME_DEDUP_H_

#include <stdint.h>

#include "simple_model.h"
#include "coder_io.h"
#include "coder.h"

class coder_qname_dedup : public coder
{
public:
    coder_qname_dedup(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qname_dedup");
        model_flag.reset();
        for (uint32_t p = 0; p < 3; ++p) {
            model_plane[p].reset();
        }
        flushed = false;
    }

    virtual ~coder_qname_dedup()
    {
        if (!flushed && io->m == coder_io::MENC) {
            encode_flush();
        }
    }

    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        uint32_t lines = 0;
        uint32_t reps = 0;
        if (!shape(lines, reps) || (uint64_t)lines + 3ull * reps != in_len) {
            io->set_err(coder_io::IO_READ_EMPTY);
            return;
        }
        if (io->m != coder_io::MENC) {
            rc.output((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartEncode();
            io->m = coder_io::MENC;
        }
        for (uint32_t i = 0; i < lines; ++i) {
            model_flag.encodeSymbol(&rc, in[i] ? 1 : 0);
        }
        for (uint32_t p = 0; p < 3; ++p) {
            const uint8_t *plane = in + lines + (size_t)p * reps;
            for (uint32_t k = 0; k < reps; ++k) {
                model_plane[p].encodeSymbol(&rc, plane[k]);
            }
        }
    }

    int32_t decode_line(uint8_t *out, uint32_t out_len,
                        [[maybe_unused]] uint8_t split_ch = UINT8_MAX,
                        [[maybe_unused]] bool need2hold = false) override
    {
        uint32_t lines = 0;
        uint32_t reps = 0;
        if (!shape(lines, reps) || (uint64_t)lines + 3ull * reps != out_len) {
            return coder_ns::CODER_ERR_BAD_ARGS;
        }
        if (io->m != coder_io::MDEC) {
            rc.input((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartDecode();
            io->m = coder_io::MDEC;
        }
        for (uint32_t i = 0; i < lines; ++i) {
            out[i] = (uint8_t)model_flag.decodeSymbol(&rc);
        }
        for (uint32_t p = 0; p < 3; ++p) {
            uint8_t *plane = out + lines + (size_t)p * reps;
            for (uint32_t k = 0; k < reps; ++k) {
                plane[k] = (uint8_t)model_plane[p].decodeSymbol(&rc);
            }
        }
        return (int32_t)out_len;
    }

    /* fake */
    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) { return 0; }

    void encode_flush() override
    {
        if (flushed) {
            return;
        }
        if (io->m == coder_io::MENC) {
            if (rc.FinishEncode() < 0) {
                io->set_err(coder_io::IO_BUF_FULL);
            }
            io->data_len += rc.size_out();
        }
        flushed = true;
    }

private:
    bool shape(uint32_t &lines, uint32_t &reps) const
    {
        if (!io->meta.isMember("lines") || !io->meta.isMember("reps")) {
            return false;
        }
        lines = io->meta["lines"].asUInt();
        reps = io->meta["reps"].asUInt();
        return true;
    }

    RangeCoder rc;
    SIMPLE_MODEL<2> model_flag;
    SIMPLE_MODEL<256> model_plane[3];
    bool flushed;
};

#endif
