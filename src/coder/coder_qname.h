/*
 * coder_qname.h - Codec dedicated to QNAME
 *
 * Background: in coordinate-sorted SAM files such as HG00106 assembled from
 * multiple runs, QNAMEs look like
 *   ERR015528.21801860
 * The fixed 9-character prefix takes only 5 distinct values (entropy ~2.3 bit),
 * while the QNAMEs of the R1/R2 records of the same fragment are identical and,
 * after coordinate sorting, usually adjacent (distance <64 for 99%+ of them).
 * coder_affix_match only matches shared prefixes between adjacent lines and is
 * essentially ineffective against random numeric ids; this codec targets those
 * two points:
 *
 *   1) Cross-line deduplication: it maintains a "recent-line hash table". When a
 *      QNAME has appeared before and the line distance is within RING_SIZE, it
 *      writes only a 1-bit marker plus the distance (distance entropy empirically
 *      ~5.1 bit) instead of repeating the whole string;
 *   2) Position-based modeling: for lines that miss, it writes each character
 *      using "line length + one adaptive model per position"; the prefix segment
 *      converges to a very small character set, and the numeric segment
 *      converges to entropy close to that of uniform digits.
 *
 * Hash table design (encoder-only; the decoder maintains only a ring buffer and
 * does not build the table):
 *   - Direct mapping: slot = hash & MASK, no linear probing, no deletion;
 *     naturally bounded and cannot loop forever.
 *   - Each slot stores the "most recent line number"; a hit requires a second
 *     verification (line within the window + identical content) to guard against
 *     hash-collision misjudgments.
 *   - The cost of a collision is only a lost copy opportunity (the slot is
 *     overwritten by another line within the window); when dist is small (median
 *     19 in this file), the hit-rate loss is negligible.
 *
 * Input convention: the whole QNAME includes the trailing '\t' (consistent with
 * compressIdFieldInAll), and in_len<=255.
 */
#ifndef _CODER_QNAME_H_
#define _CODER_QNAME_H_

#include <stdint.h>
#include <cstring>

#include "simple_model.h"
#include "coder_io.h"
#include "coder.h"

class coder_qname : public coder
{
public:
    static const uint32_t DM_BITS = 18;
    static const uint32_t DM_SIZE = 1u << DM_BITS;
    static const uint32_t RING_SIZE = 4096;
    /* Longest supported QNAME (including '\t'): SAM spec QNAME<=254, +tab=255. */
    static const uint32_t MAX_POS = 256;

    coder_qname(coder_io *io)
    {
        this->io = io;
        this->io->m = coder_io::MUNSET;
        this->io->appen_magic("coder_qname");
        this->flushed = false;
        this->line_idx = 0;
        this->last_len = 0;

        void* ptr = safe_alloc(MAX_POS * sizeof(SIMPLE_MODEL<256>));
        model_pos = static_cast<SIMPLE_MODEL<256>*>(ptr);
        check_exit(model_pos != nullptr, coder_ns::CODER_ERR_MEM_ALLOC_FAIL,
            "coder_qname: model_pos alloc failed");
        for (uint32_t n = 0; n < MAX_POS; n++) {
            model_pos[n].reset();
        }
        for (uint32_t n = 0; n < 256; n++) {
            model_len[n].reset();
        }
        model_copy.reset();
        model_dist.reset();
        memset(dm, 0xff, sizeof(dm)); /* line=0xffffffff means empty */
        memset(rlen, 0, sizeof(rlen));
    }

    virtual ~coder_qname()
    {
        if (model_pos) {
            safe_free((void**)&model_pos);
        }
        if (!flushed && io->m == coder_io::MENC) {
            encode_flush();
        }
    }

    /* Encode one line of QNAME (including the trailing '\t'). */
    void encode_line(const uint8_t *in, const uint32_t in_len, [[maybe_unused]] bool need2hold = false)
    {
        check_exit(in_len <= MAX_POS, coder_ns::CODER_ERR_BAD_ARGS,
            "coder_qname: QNAME too long (%u)", in_len);

        if (io->m != coder_io::MENC) {
            rc.output((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartEncode();
            io->m = coder_io::MENC;
        }

        uint32_t h = hash(in, in_len);
        const uint32_t ridx = line_idx % RING_SIZE;
        const uint32_t slot = h & (DM_SIZE - 1);

        /* Probe for the most recent occurrence: direct-mapped slot hit + line within window + identical content -> COPY. */
        int32_t dist = -1;
        if (dm[slot].hash == h && dm[slot].line != UINT32_MAX) {
            uint32_t lastLine = dm[slot].line;
            if (lastLine < line_idx && line_idx - lastLine <= RING_SIZE - 1 &&
                rline[lastLine % RING_SIZE] == lastLine &&
                rlen[lastLine % RING_SIZE] == in_len &&
                memcmp(recent[lastLine % RING_SIZE], in, in_len) == 0) {
                dist = (int32_t)(line_idx - lastLine);
            }
        }

        if (dist >= 0) {
            model_copy.encodeSymbol(&rc, 1);
            model_dist.encodeSymbol(&rc, (uint16_t)dist);
        } else {
            model_copy.encodeSymbol(&rc, 0);
            model_len[last_len].encodeSymbolOrder(&rc, in_len);
            for (uint32_t p = 0; p < in_len; p++) {
                model_pos[p].encodeSymbol(&rc, in[p]);
            }
        }

        last_len = (int32_t)in_len;

        /* Update the ring buffer and the direct-mapping table. */
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)in_len;
        memcpy(recent[ridx], in, in_len);

        dm[slot].hash = h;
        dm[slot].line = line_idx;

        line_idx++;
    }

    /* Decode one line of QNAME (including the trailing '\t'), returning the number of bytes written to out. */
    int32_t decode_line(uint8_t *out, uint32_t out_len,
                        [[maybe_unused]] uint8_t split_ch = UINT8_MAX,
                        [[maybe_unused]] bool need2hold = false)
    {
        if (io->m != coder_io::MDEC) {
            rc.input((char *)io->data, (char *)io->data + io->data_capacity);
            rc.StartDecode();
            io->m = coder_io::MDEC;
        }

        uint16_t copy = model_copy.decodeSymbol(&rc);
        uint32_t len;
        if (copy) {
            uint32_t dist = model_dist.decodeSymbol(&rc);
            /* dist can be line_idx (copying line 0), but must be >= 1 and <= line_idx. */
            if (dist == 0 || dist > line_idx) {
                return coder_ns::CODER_ERR_BAD_ARGS;
            }
            uint32_t ridx = (line_idx - dist) % RING_SIZE;
            len = rlen[ridx];
            if (len > out_len) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            memcpy(out, recent[ridx], len);
        } else {
            len = model_len[last_len].decodeSymbolOrder(&rc);
            if (len > out_len || len > MAX_POS) {
                return coder_ns::CODER_ERR_BUF_SMALL;
            }
            for (uint32_t p = 0; p < len; p++) {
                out[p] = (uint8_t)model_pos[p].decodeSymbol(&rc);
            }
        }
        last_len = (int32_t)len;

        /* Update the ring buffer consistently with the encoder (no hash table: decoding needs no probing). */
        uint32_t ridx = line_idx % RING_SIZE;
        rline[ridx] = line_idx;
        rlen[ridx] = (uint8_t)len;
        memcpy(recent[ridx], out, len);

        line_idx++;
        return (int32_t)len;
    }

    /* fake */
    int32_t decode_line(uint8_t *, uint32_t, uint8_t *, uint8_t, bool) { return 0; }

    void encode_flush()
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
    static uint32_t hash(const uint8_t *s, uint32_t n)
    {
        uint32_t h = 2166136261u;
        for (uint32_t i = 0; i < n; i++) {
            h ^= s[i];
            h *= 16777619u;
        }
        return h;
    }

    RangeCoder rc;
    SIMPLE_MODEL<2> model_copy;
    SIMPLE_MODEL<4096> model_dist;
    SIMPLE_MODEL<256> model_len[256];
    SIMPLE_MODEL<256> *model_pos;

    struct HashSlot {
        uint32_t hash;
        uint32_t line;
    } dm[DM_SIZE];
    uint8_t rlen[RING_SIZE];
    uint32_t rline[RING_SIZE];
    uint8_t recent[RING_SIZE][MAX_POS];

    uint32_t line_idx;
    int32_t last_len;
    bool flushed;
};

#endif
