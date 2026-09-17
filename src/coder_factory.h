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

#include "pbgz_types.h"
#include "preprocess_info.h"
#include "coder/coder.h"

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
};

class CoderFactory {
public:
    /* Metadata row for a coder type; nullptr when the type has no implementation
     * (e.g. the removed coder_simple_rc). */
    static const CoderDescriptor* descriptor(CoderType type);

    /* Metadata row for a bitstream magic; nullptr when unknown. */
    static const CoderDescriptor* descriptorByMagic(const std::string& magic);

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
