/*
 * coder_factory.h - create encoders/decoders by type or magic
 *
 * Previously each actuator directly new'd a concrete encoder class, with the
 * encoder type hard-coded in the code. The preprocessing phase (CodecSelector)
 * trial-compressed to find the best encoder per field and stored it in
 * PreprocessInfo, but actuators never read it, so the selection result never
 * took effect.
 *
 * This factory centralizes the "type -> instance" mapping so actuators can
 * create encoders dynamically according to the preprocessing result, while
 * guaranteeing a usable fallback under any abnormal condition.
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

#include "preprocess_info.h"
#include "coder/coder.h"

class CoderFactory {
public:
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
     * Currently only coder_bwt_cm truly consumes the level (internal BWT block
     * size, 0-9), to which compressLevel 1-9 maps directly; the parameter only
     * takes effect on the encoding side, since on the decoding side the block
     * size is read from the bitstream, so changing it does not affect
     * compatibility. coder_affix_match also has a level (1-2), but its decoding
     * side calls set_level only after construction (ineffective, always the
     * default 2), and once the encoding side sets it to 1 it would disagree
     * with the decoding side and corrupt data, so it stays at the default 2 and
     * does not follow compressLevel. The remaining encoders (fc/fcv2/qual) do
     * not consume the level, so setting it is a no-op.
     */
    static void applyLevel(coder_io* io, CoderType type, uint8_t compressLevel);

    /*
     * Whether this factory can create the given type.
     *
     * Preprocessing uses it to filter trial candidates: if the factory cannot
     * build a type, even a trial win cannot be actually used at compression
     * time, so it is better not to include it in the comparison at all.
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
