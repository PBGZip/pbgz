/*
 * coder_golomb.h - Golomb-Rice byte-stream coder
 *
 * A byte-stream coder with the same external interface as coder_bwt_cm and
 * coder_affix_match (encode_line / encode_flush / decode_line / decode_inlen),
 * and with the same block-buffering behavior: the compression level selects
 * one "compress at once" block size, and input smaller than one block is
 * cached until a whole block can be committed (see encode_line).
 *
 * Golomb-Rice with per-4096-byte-segment adaptive k. Symbols are bytes
 * (0..255). Each byte x is coded with parameter k: q = x >> k as a unary code
 * (q '1' bits + a '0' separator), then the k low bits of x verbatim (MSB
 * first). Per segment the encoder scans for the k in [0, 8] that minimizes
 * the total code length and writes it into the segment header, so the decoder
 * needs no extra statistics. If even the best k does not beat 8 bits/byte,
 * the whole block is stored raw.
 *
 * Golomb codes are optimal for geometric distributions (small values
 * frequent, large values rare), which is exactly the shape of the LEB128
 * varint POS-delta streams produced by SamCodecActuator::compressPosFieldDelta;
 * the segmented k adapts to slowly changing means within a block.
 *
 * Bitstream block format (each block is independent):
 *   [u32 LE]  source length N, 0 = end of stream
 *   [u8]      flags: 0x00 = raw, 0x01 = Golomb
 *   [payload] N source bytes (raw); for Golomb, per 4096-byte segment a
 *             [u8 k] header followed by the segment bitstream padded to a
 *             byte boundary at block end.
 *
 * Error convention is the same as the other coders: decode_line returns < 0 on
 * a corrupt/truncated stream, 0 on clean EOF and > 0 on the decoded length.
 */

#ifndef _CODER_GOLOMB_H_
#define _CODER_GOLOMB_H_

#include <stdint.h>
#include <algorithm>

#include "coder_io.h"
#include "coder.h"

namespace {
constexpr int32_t kGolombSegSize = 4096; /* Golomb k-segment length */
}

class coder_golomb : public coder
{
public:
    coder_golomb(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_golomb");

        coder_buff = nullptr;
        buff_capacity = 0;
        curr_incache = 0;
        curr_out_offset = 0;
        curr_blk_len = 0;
        bsize = 0;
        bitbuf = 0;
        nbits = 0;
        flushed = false;
    }

    virtual ~coder_golomb()
    {
        if (!flushed && io->m == coder_io::MENC)
            encode_flush();
        if (coder_buff)
            safe_free((void**)&coder_buff);
    }

    /* External compression interface. */
    void encode_line(const uint8_t *in, const uint32_t in_len,
                     [[maybe_unused]] bool need2hold = false) override
    {
        if (io->m != coder_io::MENC)
        { /* Initialize when not initialized */
            const int32_t tab[10] =
                {
                    0,
                    1 << 20,    // -1 - 1 MB
                    1 << 22,    // -2 - 4 MB
                    1 << 23,    // -3 - 8 MB
                    0x00FFFFFF, // -4 - ~16 MB (Default)
                    1 << 25,    // -5 - 32 MB
                    1 << 26,    // -6 - 64 MB
                    1 << 27,    // -7 - 128 MB
                    1 << 28,    // -8 - 256 MB
                    0x7FFFFFFF, // -9 - ~2 GB
                };
            const int32_t BLOCK_SIZE = (268435456);

            io->m = coder_io::MENC;
            int32_t level = (io->meta["level"].isInt()) ? (io->meta["level"].asInt()) : 4;
            io->meta["level"] = (Json::Value::Int)level;
            /* level 0 would leave bsize = 0 and break the buffering; engine
               levels are always 1..9, clamp defensively. */
            if (level < 1) level = 1;
            if (level > 9) level = 9;

            bsize = std::min(tab[level], BLOCK_SIZE); /* Block size */
            coder_buff = static_cast<uint8_t*>(safe_alloc(bsize));
            check_exit(coder_buff, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
                       "coder_golomb: insufficient memory for %d-byte block buffer", bsize);
            buff_capacity = bsize;
        }

        uint32_t off = 0;
        while (off < in_len)
        {
            /* Fill up to one block at a time; a single encode_line call may
               carry more than one block's worth of data (e.g. a whole SAM
               block assembled in memory), so commit whole blocks as they are
               filled and keep only the tail cached. */
            const int32_t len2add = bsize - curr_incache;
            const uint32_t take = std::min((uint32_t)(in_len - off), (uint32_t)len2add);
            memcpy(coder_buff + curr_incache, in + off, take);
            curr_incache += (int32_t)take;
            off += take;
            if (curr_incache == bsize)
            {
                encode_one_block(bsize);
                curr_incache = 0;
            }
        }
    }

    /* Call this after encoding to compress cached data */
    void encode_flush() override
    {
        if (io->m != coder_io::MENC || flushed)
            return;
        /* First compress all cached data */
        if (curr_incache > 0)
            encode_one_block(curr_incache);
        /* Then set EOF marker */
        put32(0);
        flushed = true;
    }

    /* External decompression interface, returns actual decompressed length,
       exits when encountering split_ch during decompression. */
    int32_t decode_line(uint8_t *out, uint32_t out_len, uint8_t split_ch = UINT8_MAX,
                        bool need2hold [[maybe_unused]] = false) override
    {
        uint8_t ch;
        int32_t len = 0, n;

        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;

        if (io->m != coder_io::MDEC)
        {
            io->m = coder_io::MDEC;
            int32_t ret = decode_one_block();
            if (ret < 0)
                return ret;
            if (ret == 0)
                return 0;
        }

        if (split_ch != UINT8_MAX)
        {
            for (;;)
            {
                if (curr_out_offset == curr_blk_len)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < curr_blk_len; n++)
                {
                    ch = coder_buff[n];
                    *(out + len++) = ch;
                    if (ch == split_ch)
                    {
                        curr_out_offset = ++n;
                        return len;
                    }
                    if (len >= (int32_t)out_len)
                    {
                        /* Filled up without finding the separator: the output
                           buffer is too small. */
                        curr_out_offset = ++n;
                        io->set_err(coder_io::IO_BUF_FULL);
                        return coder_ns::CODER_ERR_BUF_SMALL;
                    }
                }
                curr_out_offset = n;
            }
        }
        else
        {
            for (;;)
            {
                if (curr_out_offset == curr_blk_len)
                { /* Current block decompression completed, need to decompress again */
                    int32_t ret = decode_one_block();
                    if (ret < 0)
                        return ret;
                    if (ret == 0)
                        goto stream_end;
                }
                for (n = curr_out_offset; n < curr_blk_len; n++)
                {
                    *(out + len++) = coder_buff[n];
                    if (len >= (int32_t)out_len)
                    {
                        /* Fixed-length mode: filling out_len is a normal end; the rest is left for the next call. */
                        curr_out_offset = ++n;
                        return len;
                    }
                }
                curr_out_offset = n;
            }
        }

stream_end:
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        return len;
    }

    /* fake */
    int32_t decode_line(uint8_t*, uint32_t, uint8_t*, uint8_t, bool) override
    {
        return 0;
    }

    /* External interface to get the length of currently consumed data */
    int32_t decode_inlen()
    {
        return this->io->data_len;
    }

private:
    /* Best Golomb parameter and its bit length for a segment. */
    struct SegPick {
        int32_t k;
        uint64_t bits;
    };

    static SegPick pickK(const uint8_t *buf, int32_t len)
    {
        SegPick best;
        best.k = 0;
        best.bits = UINT64_MAX;
        for (int32_t kk = 0; kk <= 8; kk++)
        {
            uint64_t bits = 0;
            for (int32_t i = 0; i < len; i++)
            {
                bits += (buf[i] >> kk) + 1 + kk;
                if (bits >= best.bits)
                    break; /* prune */
            }
            if (bits < best.bits)
            {
                best.bits = bits;
                best.k = kk;
            }
        }
        return best;
    }

    /* Write one block of n source bytes plus its header. */
    void encode_one_block(int32_t n)
    {
        put32((uint32_t)n);

        /* Golomb: sum the segment-wise optimal bit lengths to decide raw. */
        uint64_t totalBits = 0;
        for (int32_t off = 0; off < n; off += kGolombSegSize)
            totalBits += pickK(coder_buff + off, std::min(kGolombSegSize, n - off)).bits;

        if (totalBits >= (uint64_t)n * 8)
        {
            io->putc(0x00); /* flags: raw */
            for (int32_t i = 0; i < n; i++)
                io->putc(coder_buff[i]);
            return;
        }

        io->putc(0x01); /* flags: golomb */
        bitbuf = 0;
        nbits = 0;
        for (int32_t off = 0; off < n; off += kGolombSegSize)
        {
            const int32_t segLen = std::min(kGolombSegSize, n - off);
            const int32_t k = pickK(coder_buff + off, segLen).k;
            io->putc((uint8_t)k);
            for (int32_t i = off; i < off + segLen; i++)
            {
                const int32_t x = coder_buff[i];
                const int32_t q = x >> k;
                const int32_t r = x & ((1 << k) - 1);
                for (int32_t j = 0; j < q; j++)
                    put_bit(1);
                put_bit(0);
                for (int32_t j = k - 1; j >= 0; j--)
                    put_bit((r >> j) & 1);
            }
            while (nbits > 0) /* pad every segment to a byte boundary so the
                                 next segment's k header stays byte-aligned */
                put_bit(0);
        }
    }

    /* Decode one block from the code stream into coder_buff; returns the
       source length, 0 on the end marker, < 0 on a corrupt stream. */
    int32_t decode_one_block()
    {
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        int32_t n = (int32_t)get32();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END; /* stream truncated while reading length */
        if (n == 0)
            return 0;
        if (n < 0 || n > (1 << 30)) /* sanity bound for a corrupt header */
        {
            coder_logger(coder_ns::ERROR, "coder_golomb: block size %d is invalid", n);
            return coder_ns::CODER_ERR_INNER;
        }

        ensure_buffer((uint32_t)n);

        const int32_t flags = io->getc();
        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        if (flags == 0x00)
        {
            for (int32_t i = 0; i < n; i++)
                coder_buff[i] = io->getc();
        }
        else if (flags == 0x01)
        {
            bitbuf = 0;
            nbits = 0;
            int32_t off = 0;
            while (off < n)
            {
                const int32_t k = io->getc();
                if (io->err != coder_io::IO_OK)
                    return coder_ns::CODER_ERR_STREAM_END;
                if (k < 0 || k > 8)
                {
                    coder_logger(coder_ns::ERROR, "coder_golomb: k %d is invalid", k);
                    return coder_ns::CODER_ERR_INNER;
                }
                const int32_t segLen = std::min(kGolombSegSize, n - off);
                for (int32_t i = 0; i < segLen; i++)
                {
                    int32_t q = 0;
                    while (get_bit() == 1)
                    {
                        if (++q > 256) /* a byte needs at most 256 unary bits (k = 0, x = 255) */
                        {
                            coder_logger(coder_ns::ERROR, "coder_golomb: corrupt unary code");
                            return coder_ns::CODER_ERR_INNER;
                        }
                    }
                    int32_t r = 0;
                    for (int32_t j = 0; j < k; j++)
                        r = (r << 1) | get_bit();
                    coder_buff[off + i] = (uint8_t)((q << k) | r);
                }
                bitbuf = 0; /* skip the segment padding, back to byte alignment */
                nbits = 0;
                off += segLen;
            }
        }
        else
        {
            coder_logger(coder_ns::ERROR, "coder_golomb: invalid block flags %d", flags);
            return coder_ns::CODER_ERR_INNER;
        }

        if (io->err != coder_io::IO_OK)
            return coder_ns::CODER_ERR_STREAM_END;
        curr_blk_len = n;
        curr_out_offset = 0;
        return n;
    }

    /* Grow the internal buffer (capacity only ever increases). */
    void ensure_buffer(uint32_t need)
    {
        if (need <= (uint32_t)buff_capacity)
            return;
        uint8_t *nb = static_cast<uint8_t*>(safe_alloc(need));
        if (!nb)
        {
            coder_logger(coder_ns::ERROR, "coder_golomb: insufficient memory for %u-byte block buffer", need);
            io->set_err(coder_io::IO_BUF_FULL);
            return;
        }
        if (coder_buff)
        {
            memcpy(nb, coder_buff, buff_capacity);
            safe_free((void**)&coder_buff);
        }
        coder_buff = nb;
        buff_capacity = (int32_t)need;
    }

    inline void put32(uint32_t x)
    {
        io->putc((uint8_t)(x & 0xff));
        io->putc((uint8_t)((x >> 8) & 0xff));
        io->putc((uint8_t)((x >> 16) & 0xff));
        io->putc((uint8_t)((x >> 24) & 0xff));
    }

    inline uint32_t get32()
    {
        uint32_t x = 0;
        x |= (uint32_t)io->getc();
        x |= (uint32_t)io->getc() << 8;
        x |= (uint32_t)io->getc() << 16;
        x |= (uint32_t)io->getc() << 24;
        return x;
    }

    inline void put_bit(int32_t bit)
    {
        bitbuf = (bitbuf << 1) | (bit & 1);
        if (++nbits == 8)
        {
            io->putc((uint8_t)bitbuf);
            bitbuf = 0;
            nbits = 0;
        }
    }

    inline int32_t get_bit()
    {
        if (nbits == 0)
        {
            bitbuf = io->getc();
            nbits = 8;
        }
        return (int32_t)((bitbuf >> (--nbits)) & 1);
    }

private:
    uint8_t *coder_buff;  /* Block cache for encoding / decoded block for decoding */
    int32_t buff_capacity; /* Capacity of coder_buff */
    int32_t curr_incache; /* Length of currently cached encode data */
    int32_t curr_out_offset; /* Offset position during decompression */
    int32_t curr_blk_len;    /* Source length of the currently decoded block */
    int32_t bsize;           /* Encode block size selected by level */
    uint32_t bitbuf;         /* Bit accumulator for the Golomb bitstream */
    int32_t nbits;           /* Number of valid bits in bitbuf */
    bool flushed;
};
#endif
