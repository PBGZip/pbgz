/*
 * coder_factory.h - coder registry and encoder/decoder factory
 *
 * Single source of truth for which coders exist and their traits. Each coder is
 * described by one CoderDescriptor row (magic, profiles, trial behavior, level
 * policy). Adding a coder means implementing it and adding a row in
 * coder_factory.cpp, plus its name in field_coder_config.h where relevant.
 *
 * Deliberately includes no concrete encoder headers; all implementations live in
 * coder_factory.cpp to avoid a RangeCoder name clash between coder_fc.h and
 * coder_qual.h in the same translation unit.
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "preprocess_info.h"
#include "coder/coder.h"

/*
 * The quality column's record-level encoder. QUAL is the one field whose coders are not
 * interchangeable byte-stream compressors (coder_qual codes against SEQ, coder_fcv2 against
 * sequencing cycle and strand, neither inheriting coder). Lets makeQualEncoder() return
 * something feedable record by record without knowing which coder a trial picked.
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

    /*
     * Whether this coder can code the given quality byte. The stream's fixed alphabet (from the
     * frequency table or loaded prior) determines it, so a value outside it cannot be coded
     * (see coder_io::IO_UNCODABLE). The caller asks before feeding a block.
     */
    virtual bool coversByte(uint8_t b) const { (void)b; return true; }
};

/* Material a trial's verdict hands the coder's constructor beyond its type. Every member may
   be absent (preprocessing did not run or the sample was too small), so the coder's defaults
   stand. */
struct QualCoderArgs {
    /* The quality alphabet in frequency-descending order, as preprocessing collected it. */
    const std::vector<std::pair<uint16_t, uint16_t>>* freqTable = nullptr;
    /* The context tier the trial settled on for fcv2; null = the coder's own defaults. */
    const QualFcv2Params* fcv2Params = nullptr;
    /* The engine's trained prior, and where the loaded flag is reported back. */
    const std::vector<uint8_t>* priorBlob = nullptr;
    bool* priorLoaded = nullptr;
};

/* Decoding twin of qual_record_encoder: one call fetches a record, whichever coder the
   stream's meta says wrote it. Callers supply the target, its length, and the already-decoded
   SEQ for coders that use it. */
class qual_record_decoder {
public:
    virtual ~qual_record_decoder() = default;

    /* Fetch one record into `qual`, with the record's decoded SEQ as context (may be null).
       Returns 0 on success, negative on a corrupt or truncated stream. */
    virtual int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq,
                                  uint32_t seqLen) = 0;
};

/* What a field's stream hands its decoder beyond the stream itself: the level it was written
   with, and the file-level prior a coder_arith stream needs to stay on the encoder's model.
   Both optional; a stream carrying neither gets defaults. */
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

/* Compression profile bitmask: which -m modes may select a coder. A coder may belong
   to both. */
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
 * One coder's metadata. The create/supports pointers live in the .cpp so concrete coder
 * headers do not leak here (see the file comment).
 *
 *   profiles       bitmask of CoderProfile; the -m gate
 *   trialLineBased must be fed line by line during trial compression to be measured faithfully
 *                  (only coder_affix_match sets this)
 *   trialCandidate may take part in trial-compression selection
 *   trialUsesLevel trial compression should pick the coder's own block level
 *   minLevel       lowest engine level at which it may be trialled (0/1 = always)
 *   trialPriority  order inside a trial; smaller is tried first and wins ties
 *   levelPolicy    how applyLevel() treats the coder at real encoding time
 *   create/supports the implementation row (null when this build cannot build it)
 *   holdsCallerBuffer / wholeBlockOnly  decoder traits looked up by magic (see
 *                  descriptorByMagic), both false for an ordinary line coder
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

    /* Decoder traits of a stream, looked up by its magic (see CoderDescriptor). Unknown
       magic answers false. */
    static bool decoderHoldsCallerBuffer(const std::string& magic);
    static bool decoderIsWholeBlock(const std::string& magic);

    /* Whether this run's profile (-m mode) may select the coder. */
    static bool eligibleProfile(CoderType type, uint8_t mode);

    /* Profile and level gate: the coder may be trialled at this engine level. */
    static bool eligible(CoderType type, uint8_t mode, uint8_t compressLevel);

    /*
     * Compression side: create an encoder by the preprocessing-picked type. Always returns a
     * usable encoder, never null; unrecognized types fall back to BWT_CM (general-purpose,
     * always correct). QUAL also lands here (it cannot be created uniformly by this factory).
     */
    static std::shared_ptr<coder> makeEncoder(CoderType type, coder_io* io);

    /* Resolve a field's encoder for the active -m mode: fast always uses the field's fast
       default (coder_rans), archive keeps the caller's historical coder. */
    static CoderType fieldCoderForMode(uint32_t fieldIdx, uint8_t mode, CoderType archiveCoder);

    /* Build an encoder through fieldCoderForMode(), applying the level and replacing a coder
       that cannot be fed line by line with the mode's line fallback. This is the single entry
       point for the SAM field paths that used to instantiate coder_bwt_cm/coder_arith directly. */
    static std::shared_ptr<coder> makeProfileEncoder(uint32_t fieldIdx, uint8_t mode, coder_io* io,
                                                     CoderType archiveCoder, bool lineMode,
                                                     uint8_t compressLevel);

    /*
     * Decompression side: create a decoder by the magic recorded in the bitstream. Returns null
     * rather than a fallback on an unknown magic (the file is newer or corrupt); the caller must
     * error out, not silently skip.
     */
    /*
     * Encoder for the coder a trial selected, with the parameters that trial settled on (see
     * qual_record_encoder / QualCoderArgs). The compression level is applied after construction.
     * Implemented in qual_coder_factory.cpp (see the file comment).
     */
    static std::shared_ptr<qual_record_encoder> makeQualEncoder(CoderType picked, coder_io* io,
                                                               const QualCoderArgs& args,
                                                               uint8_t compressLevel);

    /* Decoder for the QUAL stream whose meta carries `magic`. Returns null when the magic is
       unknown or the stream cannot be opened. Implemented in qual_coder_factory.cpp. */
    static std::shared_ptr<qual_record_decoder> makeQualDecoder(const std::string& magic,
                                                               coder_io* io,
                                                               const QualDecoderArgs& args);

    /* Decoder for a field stream whose meta carries `magic`, chosen by magic rather than by
       coder class. Returns null for an unknown magic; the caller decides whether that is
       fatal. */
    static std::shared_ptr<coder> makeFieldDecoder(const std::string& magic, coder_io* io,
                                                   const FieldDecoderArgs& args);

    static std::shared_ptr<coder> makeDecoder(const std::string& magic, coder_io* io);

    /* Convert the engine level into the level the encoder type accepts and write it into
       coder_io's meta. Coders that do not consume the level (per levelPolicy) ignore it. */
    static void applyLevel(coder_io* io, CoderType type, uint8_t compressLevel);

    /* Whether this factory has an implementation row with a create function for the type;
       preprocessing filters trial candidates with it. */
    static bool canMake(CoderType type);

    /* Whether an encoder can be used for a given field of a given file type; most are generic,
       a few (e.g. fcv2) have prerequisites. Lives here because it compares BlockType and
       SamField, which the coder layer cannot depend on. */
    static bool coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx);
};
