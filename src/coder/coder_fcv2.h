/*
 * coder_fcv2.h - quality-value context-mixing coder
 *
 * A dedicated coder designed for the SAM quality-value column. It is not a
 * generic byte-stream compressor; it depends on two pieces of information only
 * an aligned SAM can provide: the length of each record, and that record's
 * strand (the 0x10 bit of FLAG). It can therefore only be used on the QUAL
 * column of SAM and not on other fields, see supports().
 *
 * The basic approach is to build, for each quality-value symbol, several
 * probability models of differing granularity, combine their predictions with
 * weights in the log-odds domain, and hand the result to a binary range coder.
 * The model granularities run from coarse to fine: the coarse ones guarantee
 * enough samples at all times, the fine ones provide extra precision when the
 * data is plentiful, and the weights are adjusted by gradient descent, so no
 * manual fallback threshold is needed.
 *
 * On interface granularity: the interface is exposed per "record" rather than
 * per byte or per block. There are two reasons. First, strand and read length
 * are record-level attributes anyway, so passing them per record is most
 * natural. Second, the implementation is entirely hidden in the .cpp (see
 * below for why), so calls across translation units cannot be inlined; calling
 * per record means this overhead is amortized over the record's byte count,
 * while the true hot path - the bit-by-bit prediction and update - stays inside
 * the .cpp and is inlined as usual.
 *
 * On why the implementation must live in the .cpp: this coder needs the binary
 * range coder from coder/fc/rangecoder.h, and the RangeCoder defined there and
 * the same-named class in coder/clr.h are two different things (one is binary,
 * the other multi-symbol with frequencies). sam_actuator.cpp, which uses this
 * coder, also needs coder_qual.h, and the latter pulls in clr.h. If both
 * appeared in the same translation unit they would be redefined, so this header
 * cannot include rangecoder.h and the implementation can only go in the .cpp.
 */

#pragma once

#include <stdint.h>
#include <memory>
#include <vector>

#include "coder_io.h"

class fcv2_impl;

/*
 * Parameter tiers of the fcv2 context model. The default values are the
 * historical fixed constants; QualSelector chooses a parameter set for trial
 * compression based on data characteristics (quality-value alphabet size,
 * sample count, average read length), and the chosen tiers are passed with
 * PreprocessInfo to the coder and the prior training; the decoder reads back
 * the same set from the stream header (see encode_record/begin_decode).
 *
 * Meaning:
 *   cycleMax   / cycleBucket  upper bound and bin count for the in-record
 *                             position (sequencing-cycle index)
 *   deltaMax   / deltaBucket  upper bound and bin count for the in-read
 *                             quality-transition count (strategy 1)
 *   prevShift  right shift used to quantize the predecessor quality value in
 *              m3; larger is coarser (context is sparser with a large alphabet)
 *   useDelta   whether the m5 transition-count context is enabled; when false,
 *              deltaBucket normalizes to 1
 *   useDedup   whether adjacent duplicate read dedup is enabled (strategy 3,
 *              fqzcomp's do_dedup): each record is first compared with the
 *              previous one; if identical, only 1 bit is written and the whole
 *              quality string is skipped.
 *   useQa      whether the read average-quality bin context is enabled
 *              (strategy 4, fqzcomp's do_qa): each record first computes its
 *              average quality, quantizes it to 4 tiers written into the
 *              stream, and the tier serves as the m6 context.
 */
struct Fcv2Cfg {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  deltaMax    = 32;
    int  deltaBucket = 8;
    int  prevShift   = 1;
    bool useDelta    = true;
    bool useDedup    = false;
    bool useQa       = false;

    bool operator==(const Fcv2Cfg& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               deltaMax == o.deltaMax && deltaBucket == o.deltaBucket &&
               prevShift == o.prevShift && useDelta == o.useDelta &&
               useDedup == o.useDedup && useQa == o.useQa;
    }
    bool operator!=(const Fcv2Cfg& o) const { return !(*this == o); }
};

class coder_fcv2 {
public:
    /*
     * freqTable passes the occurrence counts of the quality values; the index
     * is the symbol value (raw byte minus '!') and the value is the count. It
     * is used to build the Huffman tree: quality-value distributions are
     * highly skewed, so giving high-frequency symbols shorter coding paths
     * significantly reduces the number of binary codings per symbol, which
     * improves speed. Note this does not affect the compression ratio - the
     * cost of arithmetic coding depends only on the predicted probability, and
     * splitting a symbol into several conditional decisions leaves the total
     * cost the same.
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable);

    /* Version taking context parameter tiers; cfg is written verbatim into the
       stream header for the decoder to read back. */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg);

    /*
     * Creates a coder from a previously exported model snapshot. When
     * modelLoaded is non-null, the actual load result is written into it: it
     * is only true when the snapshot's version, compile-time size parameters,
     * alphabet, and all model arrays pass validation. Any corrupted or
     * incompatible snapshot leaves the fixed initial model built from
     * freqTable in place, never leaving a half-restored state; this lets the
     * caller treat the snapshot as an optional performance optimization
     * without spreading error handling into the main compression flow.
     */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);

    /* Version taking context parameter tiers + a prior snapshot. */
    coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg,
               const std::vector<uint8_t>& modelBlob, bool* modelLoaded);
    ~coder_fcv2();

    /*
     * Exports the currently learned model. The snapshot carries the alphabet
     * and quantized frequencies along with it, rather than saving only the
     * counters: the counters are indexed by Huffman internal node number, and
     * node numbers can only be rebuilt from exactly the same alphabet and
     * frequencies. Returning false means the current coder has no valid
     * alphabet to export, or the output buffer could not be allocated.
     */
    bool export_model(std::vector<uint8_t>& out) const;

    /*
     * This coder needs the length and strand of every record; only the QUAL
     * column of an aligned SAM can provide these.
     *
     * The applicability check lives in CoderFactory::coderSupports rather than
     * as a member of this class. The reason is that the check must compare
     * BlockType and SamField, which are defined at the src layer; the coder
     * layer's build target contains only the coder/ directory, and a reverse
     * dependency on the upper layer would break the existing layering.
     */

    /*
     * Encodes one record of quality values.
     *
     * rev comes from the 0x10 bit of the record's FLAG. Per the SAM spec, when
     * that bit is set, SEQ and QUAL are stored in the file relative to the
     * forward strand of the reference, i.e. reversed relative to the
     * sequencer's original readout order, so the sequencer's first cycle
     * corresponds to the last byte in storage. The coder uses this to recover
     * the true cycle index - quality systematically declines with the
     * sequencing cycle, and this regularity can only be exploited if the cycle
     * index is recovered correctly. Empirically, recovering it vs. not differs
     * by 0.37 percentage points.
     *
     * For unaligned data (uBAM converted from FASTQ) this bit is always 0; the
     * storage order is then the original order, so pass false.
     *
     * seq is the base sequence (ACGTN) at the corresponding positions, seqLen
     * its length, used as the condition for the fifth context model (current
     * base + cycle index). On the reverse strand the base is taken with the
     * same mapping as the cycle index, i.e. the base at storage position
     * len-1-i, because that is the cycle that produced this quality value. When
     * seq is nullptr or seqLen is insufficient to cover a mapped position,
     * that position's base context falls into the "unknown" bin, behaving the
     * same as the version without a base sequence.
     */
    void encode_record(const uint8_t* qual, uint32_t len, bool rev,
                       const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* End of encoding; returns the number of bytes written to io. */
    int32_t encode_flush();

    /*
     * Decodes one record of quality values into dst. The strand is carried in
     * the stream and read back during decoding, so only the length is needed;
     * the decompression side need not track the FLAG field to decode QUAL.
     *
     * seq is the base sequence of the corresponding record already decoded
     * (length seqLen), used as the condition for the base context model; when
     * nullptr it falls into the "unknown" bin, and it must agree with the
     * encoding side (the encoder should also pass null, otherwise the two
     * sides' contexts disagree).
     */
    int32_t decode_record(uint8_t* dst, uint32_t len,
                          const uint8_t* seq = nullptr, uint32_t seqLen = 0);

    /* Called once before decoding; reads the alphabet and other info from the
       stream header. */
    int32_t begin_decode();

private:
    std::unique_ptr<fcv2_impl> impl;

    coder_fcv2(const coder_fcv2&) = delete;
    coder_fcv2& operator=(const coder_fcv2&) = delete;
};
