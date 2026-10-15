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
 * Mixing models available to fcv2: m0..m4, always all five. The model count is
 * fixed - there is no per-file or per-read-length variation - so neither
 * Fcv2Cfg nor the stream header carries a count: encoder and decoder agree on
 * the mix by construction.
 *
 * The two models beyond m4 are intentionally gone: m6 (read average-quality
 * bin, strategy 4) duplicates m1 while its switch is off - its context index
 * degenerates to m1's - and m5 (in-read quality-transition count, strategy 1)
 * overlaps what m2/m3 already carry. Each model costs a mix-and-update per bit,
 * which is the hot loop, so a smaller mix is faster; it can also compress
 * better, because fewer weights are easier to fit. The configuration parameters
 * that shaped m5/m6 (deltaMax, deltaBucket, useDelta, useQa) are removed along
 * with them, and the model count is pinned at 5 (m0..m4), which measured
 * against the old 6-model mix is another ~15%/~12% faster and compresses about
 * the same or a little better on short reads.
 */
const int FCV2_MODEL_COUNT = 5;

/*
 * Huffman tree cap of the value alphabet. A block with more distinct quality
 * values than this cannot be coded by fcv2 at all (the constructor leaves go
 * unused), so a caller can tell from the sample's alphabet size whether the coder
 * is applicable. The coder itself uses it as its internal TREE_CAP.
 */
const int FCV2_TREE_CAP = 64;

/*
 * Version of the fcv2 stream header (see encode_record / begin_decode), written as the stream's
 * first symbol.
 *
 * The header carries the context tiers, which fix the model layout (the mix is
 * always all FCV2_MODEL_COUNT models, so its count is implied), so a reader that
 * assumes a different layout misparses the rest of the stream instead of failing
 * at the header. Bump this whenever the parameter set, their order or the model
 * state layout changes. Compatibility with old layouts is deliberately not
 * maintained: the decoder reads back whatever this build wrote, and anything
 * else is answered with CODER_ERR_UNSUPPORTED_VERSION rather than a guess.
 *
 * History:
 *   none   alpha, then cycleMax, cycleBucket, deltaMax, deltaBucket, prevShift, useDelta,
 *          useDedup, useQa, with the mix always using every model (MODEL_COUNT = 7 then).
 *   1      version byte, then alpha, then the same eight tiers, then modelCount.
 *   2      version byte, then alpha, then cycleMax, cycleBucket, prevShift, useDedup,
 *          then modelCount. The m5 transition-count tiers (deltaMax, deltaBucket, useDelta) and
 *          the m6 read-average-quality switch (useQa) are gone with their models.
 *   3      current: version byte, then alpha, then cycleMax, cycleBucket, prevShift, useDedup.
 *          The model count no longer travels in the header because it is fixed at
 *          FCV2_MODEL_COUNT.
 */
const uint8_t FCV2_STREAM_VERSION = 3;

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
 *   prevShift  right shift used to quantize the predecessor quality value in
 *              m3; larger is coarser (context is sparser with a large alphabet)
 *   useDedup   whether adjacent duplicate read dedup is enabled (strategy 3,
 *              fqzcomp's do_dedup): each record is first compared with the
 *              previous one; if identical, only 1 bit is written and the whole
 *              quality string is skipped.
 *
 * The mixing-model count is deliberately not a tier: it is fixed at
 * FCV2_MODEL_COUNT (m0..m4), same on both sides, so it appears nowhere here.
 */
struct Fcv2Cfg {
    int  cycleMax    = 96;
    int  cycleBucket = 16;
    int  prevShift   = 1;
    bool useDedup    = false;

    bool operator==(const Fcv2Cfg& o) const
    {
        return cycleMax == o.cycleMax && cycleBucket == o.cycleBucket &&
               prevShift == o.prevShift && useDedup == o.useDedup;
    }
    bool operator!=(const Fcv2Cfg& o) const { return !(*this == o); }
};

/*
 * The QUAL context presets, as a table rather than as code in the selector.
 *
 * They are the configurations worth trial-compressing against each other; which
 * of them wins is a property of the data, not something a rule predicts, so the
 * selector trials them and keeps the smallest. This table is what makes that set
 * explicit, shared with the measurement tools, and printable under its own name.
 *
 *   cycleMax     reach of the positional context; cycles at or beyond it collapse
 *                into one slot, and it also sets the bucket width (cycleMax /
 *                cycleBucket cycles per bin)
 *   cycleBucket  position bins for m3 and for the mixing-weight rows
 *   prevShift    quantization of m3's predecessor quality values (0 = none)
 *
 * minVolumeTier is the levelVolumeTier() value from which trial-compressing this
 * preset starts being worth the time; see levelVolumeTier in qual_selector.cpp.
 * The fine preset needs enough per-block QUAL for its fine buckets to stay
 * populated, and below that it only spends trial time.
 *
 * The mixing-model count is not part of a preset: it is fixed at FCV2_MODEL_COUNT
 * on every one of them.
 */
struct Fcv2Preset {
    const char* name;
    int         minVolumeTier;
    Fcv2Cfg     cfg;
};

const int FCV2_PRESET_COUNT = 3;

/*
 * Naming is by what the preset changes, not by a ranking: they are alternatives.
 *
 * These are the three measured to be worth trialling; the other two former rows
 * were dropped because they never win on the data that reaches the trial:
 *
 *   - "coarse" (96/8/ps2) is beaten by ultra at every measured block volume and
 *     read length, so trialling it only spent time.
 *   - "default" (96/16/ps1, the Fcv2Cfg fallback) trails the three below at
 *     every measured volume (its best showing is still a couple of points
 *     worse), so it too never wins.
 *
 * The remaining rows span the ranked set on short reads, which is where fcv2 is
 * used at all (QualSelector drops the whole coder on long reads and on a quality
 * alphabet wider than its tree cap):
 *
 *   - ultra (4/ps0) leads at every measured small to medium block volume (0.4%-
 *     1.2% over the next tier).
 *   - ultra-80 (80/4/ps0) is the same bins over a shorter cycle span; it is
 *     second or third across volumes and costs almost nothing to keep.
 *   - fine (24/ps0) needs a per-block QUAL volume large enough for its fine
 *     buckets to stay populated (it leads from a few MB of QUAL per block up), so
 *     it is tried only at -l 8/9 (minVolumeTier 3, the 100k reads/block tier).
 *
 * The former "short-trained" row is gone: its entire point was an eight-times
 * finer m5 transition-count context (deltaMax 64, deltaBucket 16) and its own
 * model count of 6. With m5 and the delta tiers removed, and the model count
 * fixed at FCV2_MODEL_COUNT, it collapses to a duplicate of ultra-80, so
 * trialling it would only spend time.
 */
const Fcv2Preset FCV2_PRESETS[FCV2_PRESET_COUNT] = {
    { "ultra",    0, {   96,  4, 0, false } },
    { "ultra-80", 0, {   80,  4, 0, false } },
    { "fine",     3, {   96, 24, 0, false } },
};

/* The name of a configuration, or nullptr when it is not one of the presets. */
inline const char* fcv2PresetName(const Fcv2Cfg& cfg)
{
    for (int i = 0; i < FCV2_PRESET_COUNT; ++i) {
        const Fcv2Cfg& p = FCV2_PRESETS[i].cfg;
        if (cfg.cycleMax == p.cycleMax && cfg.cycleBucket == p.cycleBucket &&
            cfg.prevShift == p.prevShift && cfg.useDedup == p.useDedup) {
            return FCV2_PRESETS[i].name;
        }
    }
    return nullptr;
}

/*
 * The alphabet is data-derived: only the quality values that occur in the block
 * are admitted, kept in ascending byte order. That keeps the Huffman tree
 * shallow, but it also means the coder can only code the values it was built
 * from - a byte outside the alphabet is skipped by encode_record, not rejected -
 * so a coder is only ever used on the data it was built for.
 *
 * The tree cap (FCV2_TREE_CAP) is a real ceiling in practice.
 * Measured quality-byte ranges fall into two conventions, which matters when
 * reasoning about how close a file is to that ceiling:
 *
 *   Phred+33   bytes 33..93    (measured: 34..80 and 34..93, i.e. phred 1..60)
 *   Phred+64   bytes 64..105   (measured: 66..105, i.e. phred 2..41 - legacy
 *                               Illumina/Solexa encoding, still present in data
 *                               that predates the +33 switch)
 *
 * A file in one convention stays under the cap, but a file spanning both would
 * need 73 symbols, which does not fit: the constructor then leaves alphaSize at 0
 * and the coder is unusable, which the upper layer is expected to notice (see
 * CoderFactory::coderSupports).
 */
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

    /*
     * Whether the stream's alphabet contains this quality byte. The alphabet is fixed once the
     * coder is built - from the frequency table, or from the prior's own alphabet when a prior was
     * loaded - and a value outside it cannot be coded (encoding one fails the stream with
     * coder_io::IO_UNCODABLE). The caller checks the block's values first and routes the column to
     * another codec when the prior does not cover them.
     */
    bool coversByte(uint8_t b) const;

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
