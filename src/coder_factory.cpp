/*
 * coder_factory.cpp - coder registry and encoder/decoder factory
 *
 * The implementation lives in the .cpp rather than the header so that concrete
 * encoder headers do not leak to callers: coder_fc.h brings in clr.h
 * indirectly, which defines a RangeCoder; and qual_model.h (brought in by
 * coder_qual.h) has another RangeCoder of the same name but a different
 * interface. The two conflict if they appear in the same translation unit.
 * Actuators commonly need coder_qual.h, so the factory header must stay clean.
 *
 * The kCoderDescriptors table below is the registry: one row per coder, and the
 * only place that knows how to build it and which profile it belongs to. New
 * coders are added here (plus their name in field_coder_config.h); no actuator
 * needs to change.
 *
 * About coder_simple_rc: it has been removed entirely from both the encoding
 * and decoding sides.
 *
 * Measurements showed this encoder is lossy — in a standalone round-trip test,
 * all 10 data blocks failed verification and nearly 40% of run-length
 * information was lost; its attractive compression ratio was bought with lost
 * data.
 *
 * Before preprocessing selection was wired in, actuators hard-coded the encoder
 * type and it was never actually invoked, so the problem never surfaced and no
 * historical files compressed with it exist. After wiring it in, the risk
 * became real: trial compression only compares compressed size and does not
 * verify lossless round-trip, so as long as it compressed some field the
 * smallest it would be picked, producing a file that cannot be decompressed
 * back to the original data. Since there is no historical baggage, removing it
 * from both sides is the cleanest.
 */

#include "coder_factory.h"

#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_arith.h"
#include "coder/coder_rans.h"
#include "coder/coder_qcm.h"
#include "coder/coder_qname.h"
#include "field_coder_config.h"

namespace {

/*---------------------------------------------------------------------------
 * Coder implementations
 *-------------------------------------------------------------------------*/

std::shared_ptr<coder> createBwtCm(coder_io* io)     { return std::make_shared<coder_bwt_cm>(io); }
std::shared_ptr<coder> createFc(coder_io* io)        { return std::make_shared<coder_fc>(io); }
std::shared_ptr<coder> createAffix(coder_io* io)     { return std::make_shared<coder_affix_match>(io); }
std::shared_ptr<coder> createArith(coder_io* io)     { return std::make_shared<coder_arith>(io); }
std::shared_ptr<coder> createRans(coder_io* io)      { return std::make_shared<coder_rans>(io); }
std::shared_ptr<coder> createQcm(coder_io* io)       { return std::make_shared<coder_qcm>(io); }
std::shared_ptr<coder> createQname(coder_io* io)     { return std::make_shared<coder_qname>(io); }

/*---------------------------------------------------------------------------
 * Applicability checks
 *-------------------------------------------------------------------------*/

/* Generic byte-stream coders work on every field. */
bool supportsAny(uint32_t /*fileType*/, uint32_t /*fieldIdx*/) { return true; }

/*
 * coder_fcv2 needs each record's read length and strand direction, which only
 * the QUAL column of an aligned SAM/BAM provides.
 */
bool supportsFcv2(uint32_t fileType, uint32_t fieldIdx)
{
    return (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) &&
           fieldIdx == (uint32_t)SAM_QUAL;
}

/*
 * Whether coder_affix_match is a candidate is decided uniformly by the config
 * table, avoiding maintaining the trial-compression scope in two places.
 */
bool supportsAffix(uint32_t fileType, uint32_t fieldIdx)
{
    return (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) &&
           samFieldCandidate(fieldIdx, CoderType::AFFIX_MATCH);
}

/*---------------------------------------------------------------------------
 * The registry
 *-------------------------------------------------------------------------*/

/*
 * Trial priority preserves the historical order inside a trial: coder_fc is the
 * baseline, coder_affix_match and coder_bwt_cm only take the lead on a strictly
 * smaller result. Keeping it here rather than in field_coder_config.h lets a
 * candidate list be written in any order.
 */
const CoderDescriptor kCoderDescriptors[] = {
    /* type,                  profiles,               trialLineBased, trialCandidate, trialUsesLevel, minLevel, priority, levelPolicy,                 create,       supports */
    {CoderType::BWT_CM,       CODER_PROFILE_ARCHIVE,  false,          true,           true,           8,        30,       CoderLevelPolicy::SET_LEVEL, createBwtCm,  supportsAny,  false, false},
    {CoderType::FC,           CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        10,       CoderLevelPolicy::NONE,      createFc,     supportsAny,  false, true},
    {CoderType::AFFIX_MATCH,  CODER_PROFILE_ARCHIVE,  true,           true,           false,          1,        20,       CoderLevelPolicy::NONE,      createAffix,  supportsAffix, true, false},
    {CoderType::ARITH,        CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        40,       CoderLevelPolicy::SET_LEVEL, createArith,  supportsAny,  false, false},
    {CoderType::RANS,         CODER_PROFILE_FAST,     false,          true,           false,          1,        50,       CoderLevelPolicy::NONE,      createRans,   supportsAny,  false, false},

    /*
     * coder_lzma lived here. It mirrored CRAM 3.1's LZMA method, the one CRAM picks for its
     * base (SEQ) series on long-read data, where an LZ77 coder with a large window follows
     * the inter-read repetition a per-block BWT has to rediscover. It won on size - 1.7519
     * bits/base against coder_bwt_cm's 1.7590 on ERR11436629 SEQ - but at 0.82 MB/s against
     * coder_bwt_cm's ~14 MB/s, and a sweep over its own settings found no way out (hc4 buys
     * 6.4x for +7.4% of size, mode=FAST +5.3% for +18%, nice_len and depth each trade a
     * fraction of a percent for nothing). At -l9 that made it 6.3x the run time of -l8 for
     * 0.056% of the file, so it was dropped outright - coder and all, not just as a
     * candidate: the format it wrote is not something this build has to read back.
     */

    /*
     * coder_qcm is an order-2 context model written for the quality column (see
     * coder_qcm.h for the model and where the numbers come from). It is an
     * ordinary byte-stream coder - no position context, no strand - so unlike
     * the record-level quality coders it goes through the generic interface and
     * is applicable to any field; only the QUAL candidate list offers it.
     *
     * Measured on ERR11436629 (ONT, one 28.09 MB block of quality values):
     * 13,597,325 bytes against coder_bwt_cm's 13,656,746, i.e. 0.43% smaller,
     * and it is also the faster of the two (15.4 MB/s encoding, 12.5 MB/s
     * decoding, against 9.4 and 4.5). The reason coder_bwt_cm can be beaten
     * here at all is in the header note: on this column the BWT performs exactly
     * like an ideal order-1 model, so the space that is left is in order 2.
     *
     * It sits last in the trial order (priority 60) so a tie keeps the
     * incumbent, and the trial decides: a loss costs only trial time.
     */
    {CoderType::QCM,          CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        60,       CoderLevelPolicy::NONE,      createQcm,    supportsAny,  false, false},

    /*
     * coder_qname is the QNAME column's own coder: it codes a whole name per line, which is
     * why that field has a dedicated path (see compressIdFieldQname). It is never a trial
     * candidate - the field's layout decides when it is used - so its row exists for the
     * same reason the QUAL and FCV2 rows do: a stream carrying its magic can be built by
     * magic alone, with no actuator naming the class (see makeFieldDecoder).
     */
    {CoderType::QNAME,        CODER_PROFILE_ARCHIVE,  false,          false,          false,          1,        0,        CoderLevelPolicy::NONE,      createQname,  supportsAny,  false, false},

    /* QUAL (coder_qual) and FCV2 (coder_fcv2) do not inherit coder and need
     * record-level input, so they have no create function here: the QUAL column
     * has its own dedicated compression/evaluation path. Their rows still exist
     * so that profile/candidate bookkeeping can tell they are archive coders. */
    {CoderType::FCV2,         CODER_PROFILE_ARCHIVE,  false,          false,          false,          1,        0,        CoderLevelPolicy::NONE,      nullptr,      supportsFcv2, false, false},
    {CoderType::QUAL,         CODER_PROFILE_ARCHIVE,  false,          false,          false,          1,        0,        CoderLevelPolicy::NONE,      nullptr,      supportsAny,  false, false},

    /* coder_simple_rc was removed (lossy); kept only as an enum placeholder. */
    {CoderType::SIMPLE_RC,    CODER_PROFILE_NONE,     false,          false,          false,          1,        0,        CoderLevelPolicy::NONE,      nullptr,      supportsAny,  false, false},
};

const CoderDescriptor* findDescriptor(CoderType type)
{
    for (const CoderDescriptor& d : kCoderDescriptors) {
        if (d.type == type) {
            return &d;
        }
    }
    return nullptr;
}

} /* namespace */

const CoderDescriptor* CoderFactory::descriptor(CoderType type)
{
    return findDescriptor(type);
}

const CoderDescriptor* CoderFactory::descriptorByMagic(const std::string& magic)
{
    /* Reuse coderTypeToMagic so the decode dispatch key is defined in exactly
     * one place (preprocess_info.h). */
    for (const CoderDescriptor& d : kCoderDescriptors) {
        if (magic == coderTypeToMagic(d.type)) {
            return &d;
        }
    }
    return nullptr;
}

bool CoderFactory::eligibleProfile(CoderType type, uint8_t mode)
{
    const CoderDescriptor* d = findDescriptor(type);
    if (d == nullptr) {
        return false;
    }
    const uint32_t bit = (mode == PBGZ_MODE_FAST) ? CODER_PROFILE_FAST : CODER_PROFILE_ARCHIVE;
    return (d->profiles & bit) != 0;
}

bool CoderFactory::eligible(CoderType type, uint8_t mode, uint8_t compressLevel)
{
    const CoderDescriptor* d = findDescriptor(type);
    if (d == nullptr || !eligibleProfile(type, mode)) {
        return false;
    }
    return compressLevel >= d->minLevel;
}

std::shared_ptr<coder> CoderFactory::makeEncoder(CoderType type, coder_io* io)
{
    const CoderDescriptor* d = findDescriptor(type);
    if (d != nullptr && d->create != nullptr) {
        return d->create(io);
    }
    /* Unknown/uncreatable type: see the header for why BWT_CM is the safe
     * fallback on the compression side (QUAL and FCV2 land here). */
    return std::make_shared<coder_bwt_cm>(io);
}

std::shared_ptr<coder> CoderFactory::makeDecoder(const std::string& magic, coder_io* io)
{
    const CoderDescriptor* d = descriptorByMagic(magic);
    if (d == nullptr || d->create == nullptr) {
        return nullptr;
    }
    return d->create(io);
}

bool CoderFactory::decoderHoldsCallerBuffer(const std::string& magic)
{
    const CoderDescriptor* d = descriptorByMagic(magic);
    return (d != nullptr) && d->holdsCallerBuffer;
}

bool CoderFactory::decoderIsWholeBlock(const std::string& magic)
{
    const CoderDescriptor* d = descriptorByMagic(magic);
    return (d != nullptr) && d->wholeBlockOnly;
}

std::shared_ptr<coder> CoderFactory::makeFieldDecoder(const std::string& magic, coder_io* io,
                                                      const FieldDecoderArgs& args)
{
    std::shared_ptr<coder> decoder = makeDecoder(magic, io);
    if (decoder == nullptr) {
        return nullptr;
    }
    if (args.level >= 0) {
        decoder->set_level((uint8_t)args.level);
    }
    /*
     * A coder_arith decoder must use the same file-level prior the encoder did, or the
     * arithmetic decode diverges on the first symbol. The magic was checked by makeDecoder
     * above, so the cast is safe.
     */
    if (magic == "coder_arith" && args.posPrior != nullptr && !args.posPrior->empty()) {
        static_cast<coder_arith*>(decoder.get())
            ->set_prior(args.posPrior->data(), (uint32_t)args.posPrior->size());
    }
    return decoder;
}

void CoderFactory::applyLevel(coder_io* io, CoderType type, uint8_t compressLevel)
{
    if (io == nullptr || compressLevel < 1 || compressLevel > 9) {
        return;
    }
    const CoderDescriptor* d = findDescriptor(type);
    if (d != nullptr && d->levelPolicy == CoderLevelPolicy::SET_LEVEL) {
        io->set_level(compressLevel);
    }
}

bool CoderFactory::coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx)
{
    const CoderDescriptor* d = findDescriptor(type);
    if (d == nullptr || d->supports == nullptr) {
        return true;
    }
    return d->supports(fileType, fieldIdx);
}

bool CoderFactory::canMake(CoderType type)
{
    const CoderDescriptor* d = findDescriptor(type);
    return (d != nullptr) && (d->create != nullptr);
}
