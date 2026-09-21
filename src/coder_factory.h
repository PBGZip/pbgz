/*
 * coder_factory.h - coder registry and encoder/decoder factory
 *
 * This header is the single source of truth for "which coders exist and what
 * are they good at". Every coder is described by one CoderDescriptor row:
 * its bitstream magic, the compression profiles (-m) it may be selected in,
 * whether it supports line-by-line accumulation, how it takes part in trial
 * compression, and how the engine level maps onto it.
 *
 * Adding a coder therefore means implementing a coder subclass and adding one
 * row to the descriptor table in coder_factory.cpp plus, where relevant, its
 * name to a field's candidate list in field_coder_config.h. No actuator change
 * is required: generic fields obtain their encoder through makeEncoder() driven
 * by the preprocessing result, and the preprocessing trial walks the configured
 * candidate list through the same table.
 *
 * Note that this header deliberately includes no concrete encoder headers; all
 * implementations live in coder_factory.cpp. The reason is that coder_fc.h
 * brings in clr.h indirectly, while qual_model.h (brought in by coder_qual.h)
 * has another RangeCoder of the same name but a different interface; the two
 * conflict when they appear in the same translation unit. After moving the
 * implementation into the .cpp, callers that include this header only see
 * declarations and are not affected by that dependency.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "pbgz_types.h"
#include "preprocess_info.h"
#include "coder/coder.h"

/*
 * The quality column's record-level encoder.
 *
 * QUAL is the one field whose coders are not interchangeable byte-stream compressors:
 * coder_qual codes each record against its SEQ, coder_fcv2 codes it against its sequencing
 * cycle and strand, and neither inherits coder. This interface is what lets
 * makeQualEncoder() hand back something the caller can feed record by record without
 * knowing which of them a trial picked.
 */
class qual_record_encoder {
public:
    virtual ~qual_record_encoder() = default;

    /* Append one record: its QUAL bytes, the SEQ it belongs to (may be null) and the strand
       direction. Only a coder that restores the sequencing cycle reads the last two. */
    virtual void encode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq, uint32_t seqLen,
                               bool rev) = 0;

    /* Commit the stream. */
    virtual void flush() = 0;
};

/*
 * What a trial's verdict hands over, beyond the type it picked: the material the coder's
 * constructor wants. Every member may be absent - preprocessing did not run, or the sample
 * was too small to be trusted - in which case the coder's own defaults stand.
 */
struct QualCoderArgs {
    /* The quality alphabet in frequency-descending order, as preprocessing collected it. */
    const std::vector<std::pair<uint16_t, uint16_t>>* freqTable = nullptr;
    /* The context tier the trial settled on for fcv2; null = the coder's own defaults. */
    const QualFcv2Params* fcv2Params = nullptr;
    /* The engine's trained prior, and where the loaded flag is reported back. */
    const std::vector<uint8_t>* priorBlob = nullptr;
    bool* priorLoaded = nullptr;
};

/*
 * The decoding twin of qual_record_encoder: one call fetches a record, whichever coder the
 * stream's meta says wrote it. What the actuator has to supply is the record's target, its
 * length, and - for the coders that use it - the SEQ that was already decoded before the
 * quality column (the same material the encoding side fed encode_record).
 */
class qual_record_decoder {
public:
    virtual ~qual_record_decoder() = default;

    /* Fetch one record into `qual`, with the record's decoded SEQ as context (may be null).
       Returns 0 on success, negative on a corrupt or truncated stream. */
    virtual int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq,
                                  uint32_t seqLen) = 0;
};

/*
 * What a field's stream may hand its decoder beyond the stream itself: the level the field
 * was written with, and the file-level prior a coder_arith stream needs to stay on the
 * encoder's model. Both are optional; a stream that carries neither gets the defaults.
 */
struct FieldDecoderArgs {
    /* The level the stream was written with; -1 = its meta recorded none, so the coder's own
       is left alone (which is not the same as level 0: a stream that says 0 had set_level(0)
       applied to its encoder, and this side must do the same). */
    int level = -1;
    /* coder_arith only: the position prior, in the form set_prior expects. */
    const std::vector<uint8_t>* posPrior = nullptr;
};

/* What the stream's meta and the engine hand the decoder's constructor. */
struct QualDecoderArgs {
    /* The alphabet, as the stream's own frequency-table sub-stream carried it. */
    const std::vector<std::pair<uint16_t, uint16_t>>* freqTable = nullptr;
    /* Whether the stream was written from a trained prior (its meta carries the address). */
    bool priorRequired = false;
    int64_t priorAddress = -1;
    /* The prior snapshot itself; a load failure is fatal on this side (see the factory). */
    const std::vector<uint8_t>* priorBlob = nullptr;
};

/*
 * Compression profile bitmask: a coder declares which -m modes may select it.
 * fast holds the speed-oriented coders, archive the ratio-oriented ones; a coder
 * may belong to both.
 */
enum CoderProfile : uint32_t {
    CODER_PROFILE_NONE    = 0u,
    CODER_PROFILE_FAST    = 1u << 0,
    CODER_PROFILE_ARCHIVE = 1u << 1,
};

/* How the engine compression level maps onto a coder. */
enum class CoderLevelPolicy : uint8_t {
    NONE = 0,      /* the coder does not consume the level */
    SET_LEVEL,     /* forward the level through coder_io::set_level */
};

/*
 * One coder's metadata. The create/supports pointers live in the .cpp so that
 * concrete coder headers do not leak to this header (see the file comment).
 *
 *   profiles        bitmask of CoderProfile; the -m gate
 *   trialLineBased  must be fed line by line during trial compression to be
 *                   measured faithfully (its gain comes from adjacent lines;
 *                   only coder_affix_match sets this)
 *   trialCandidate  may take part in trial-compression selection
 *   trialUsesLevel  trial-compression should pick the coder's own block level
 *   minLevel        lowest engine level at which it may be trialled (0/1 = always)
 *   trialPriority   order inside a trial; smaller is tried first and wins ties
 *   levelPolicy     how applyLevel() treats the coder at real encoding time
 *   create/supports the implementation row (null when this build cannot build it)
 *
 * The last two are the traits a caller needs to know about a *decoder* without holding one:
 * it has the magic the stream carries, which is exactly what the descriptor is keyed by
 * (see descriptorByMagic). They exist so that no caller has to compare magic strings to find
 * out how a coder wants to be fed - the answer belongs to the coder, not to the caller.
 *
 *   holdsCallerBuffer  this coder's decoder points at the caller's buffer and keeps pointing
 *                      there until the next call unless it is told to hold
 *                      (coder_affix_match: that reference *is* its cross-line context, so a
 *                      caller that reads a record before decoding the next one must set
 *                      need2hold). Every other coder keeps its state inside itself.
 *   wholeBlockOnly     this coder's decoder can only decode its whole stream in one call
 *                      (coder_fc), so the caller stages the block and slices it per record
 *                      itself. Every other coder reads record by record.
 *
 * Both default to false, so a row that does not name them is an ordinary line coder.
 */
struct CoderDescriptor {
    CoderType        type;
    uint32_t         profiles;
    bool             trialLineBased;
    bool             trialCandidate;
    bool             trialUsesLevel;
    uint8_t          minLevel;
    uint8_t          trialPriority;
    CoderLevelPolicy levelPolicy;
    std::shared_ptr<coder> (*create)(coder_io* io);
    bool             (*supports)(uint32_t fileType, uint32_t fieldIdx);
    bool             holdsCallerBuffer;
    bool             wholeBlockOnly;
};

class CoderFactory {
public:
    /* Metadata row for a coder type; nullptr when the type has no implementation
     * (e.g. the removed coder_simple_rc). */
    static const CoderDescriptor* descriptor(CoderType type);

    /* Metadata row for a bitstream magic; nullptr when unknown. */
    static const CoderDescriptor* descriptorByMagic(const std::string& magic);

    /*
     * The decoder traits of a stream the caller has not built a decoder for yet, looked up
     * by the magic that stream carries (see CoderDescriptor). An unknown magic answers
     * false, the same way the caller's own decoder lookup would fail.
     */
    static bool decoderHoldsCallerBuffer(const std::string& magic);
    static bool decoderIsWholeBlock(const std::string& magic);

    /* Whether this run's profile (-m mode) may select the coder. */
    static bool eligibleProfile(CoderType type, uint8_t mode);

    /* Profile and level gate: the coder may be trialled at this engine level. */
    static bool eligible(CoderType type, uint8_t mode, uint8_t compressLevel);

    /*
     * Compression side: create an encoder by the type picked during
     * preprocessing.
     *
     * This must always return a usable encoder, never a null pointer.
     * Compression is a one-way process: if a field is skipped because no encoder
     * could be obtained, the produced file is incomplete and the problem only
     * surfaces at decompression time. So unrecognized types (including ones
     * added in the future but unknown to this version) all fall back to BWT_CM,
     * a general-purpose encoder that handles every field; its compression ratio
     * may not be optimal but it is always correct.
     *
     * CoderType::QUAL also goes down this fallback path: coder_qual does not
     * inherit from the coder base class and also needs an extra frequency table
     * at construction, so this factory cannot create it uniformly; the quality
     * field has its own dedicated compression function. Generic fields should
     * never be picked as QUAL (it is not even in the trial candidates), and if
     * it ever happens the selection result is anomalous, so falling back to
     * BWT_CM is safe.
     */
    static std::shared_ptr<coder> makeEncoder(CoderType type, coder_io* io);

    /*
     * Decompression side: create a decoder by the magic recorded in the
     * bitstream.
     *
     * Unlike the compression side, this must return a null pointer rather than a
     * fallback. At decompression time the magic is a fact written in the file,
     * stating which encoder originally compressed this data. If the magic is
     * unrecognized, the file was written by a newer version of pbgz, or the
     * data is corrupt. Picking some encoder at random and decompressing with it
     * would only produce garbage, and most likely without any error — far more
     * dangerous than failing outright.
     *
     * On receiving a null pointer the caller must error out and abort; it must
     * not silently skip.
     */
    /*
     * The encoder for the coder a trial selected, with the parameters that trial settled on
     * (see qual_record_encoder / QualCoderArgs). The compression level is applied to the
     * coder that goes through the registry, after it is constructed - the coders that read
     * it do so on first use.
     *
     * Implemented in qual_coder_factory.cpp, not here: coder_qual.h and coder_fcv2.h cannot
     * be included in this translation unit (see the file comment above), so the two
     * record-level coders are built there and everything else is delegated to makeEncoder.
     */
    static std::shared_ptr<qual_record_encoder> makeQualEncoder(CoderType picked, coder_io* io,
                                                               const QualCoderArgs& args,
                                                               uint8_t compressLevel);

    /*
     * The decoder for the QUAL stream whose meta carries `magic` - the value the encoder
     * wrote there, which is what the stream itself says it is (see qual_record_decoder).
     * Returns null when the magic is unknown or the stream cannot be opened; the reason is
     * logged there, where the coder's own error codes are understood.
     *
     * Implemented in qual_coder_factory.cpp for the same reason as makeQualEncoder.
     */
    static std::shared_ptr<qual_record_decoder> makeQualDecoder(const std::string& magic,
                                                               coder_io* io,
                                                               const QualDecoderArgs& args);

    /*
     * The decoder for a field stream whose meta carries `magic`, with the arguments that
     * stream needs (see FieldDecoderArgs). This is what a caller that only knows the meta
     * asks for - it no longer names coder classes, which is the point: a coder's decoder is
     * chosen by the magic the encoder wrote, and adding one never touches an actuator again.
     *
     * Returns null for a magic this build has no decoder for; the caller decides whether
     * that is fatal, because only it knows which field it is assembling.
     */
    static std::shared_ptr<coder> makeFieldDecoder(const std::string& magic, coder_io* io,
                                                   const FieldDecoderArgs& args);

    static std::shared_ptr<coder> makeDecoder(const std::string& magic, coder_io* io);

    /*
     * Before creating an encoder, convert the engine's compression level into
     * the level accepted by that encoder type and write it into coder_io's
     * meta. The encoder reads this value at construction / first encoding and
     * writes it into the block meta (the decoding side replays it from
     * meta["coder"]["level"]).
     *
     * Which coders consume the level is declared by their descriptor
     * (levelPolicy); the rest ignore it, so setting it is a no-op for them.
     */
    static void applyLevel(coder_io* io, CoderType type, uint8_t compressLevel);

    /*
     * Whether this factory can create the given type (i.e. it has an
     * implementation row with a create function). Preprocessing uses it to
     * filter trial candidates: an uncreatable type could never be used even if
     * it won a trial.
     */
    static bool canMake(CoderType type);

    /*
     * Whether an encoder can be used for a given field of a given file type.
     *
     * Most encoders are generic byte-stream compressors that work on any field.
     * A few have extra prerequisites: fcv2 needs each record's length and
     * strand direction, which only the QUAL column of an aligned SAM provides;
     * elsewhere that information is unavailable. Ask this question before
     * trial-compressing so an unusable encoder is not picked.
     *
     * This check lives here rather than on each encoder itself because it
     * compares BlockType and SamField, and the coder layer's compilation target
     * does not include the src directory, so it cannot depend upward.
     */
    static bool coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx);
};
