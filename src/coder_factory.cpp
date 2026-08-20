/*
 * coder_factory.cpp - implementation of the encoder factory
 *
 * The implementation lives in the .cpp rather than the header so that concrete
 * encoder headers do not leak to callers: coder_fc.h brings in clr.h
 * indirectly, which defines a RangeCoder; and qual_model.h (brought in by
 * coder_qual.h) has another RangeCoder of the same name but a different
 * interface. The two conflict if they appear in the same translation unit.
 * Actuators commonly need coder_qual.h, so the factory header must stay clean.
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
#include "field_coder_config.h"

std::shared_ptr<coder> CoderFactory::makeEncoder(CoderType type, coder_io* io)
{
    switch (type) {
    case CoderType::FC:
        return std::make_shared<coder_fc>(io);

    case CoderType::AFFIX_MATCH:
        return std::make_shared<coder_affix_match>(io);

    /*
     * QUAL goes through the fallback: coder_qual does not inherit from the
     * coder base class and also needs an extra frequency table at construction,
     * so this factory cannot create it uniformly; the quality field has its own
     * dedicated compression function. Generic fields should never be picked as
     * QUAL (it is not even in the trial candidates), and if it ever happens the
     * selection result is anomalous, so falling back to BWT_CM is safe.
     */
    case CoderType::QUAL:
    case CoderType::SIMPLE_RC:
    case CoderType::BWT_CM:
    default:
        return std::make_shared<coder_bwt_cm>(io);
    }
}

std::shared_ptr<coder> CoderFactory::makeDecoder(const std::string& magic, coder_io* io)
{
    if (magic == "coder_bwt_cm") {
        return std::make_shared<coder_bwt_cm>(io);
    }
    if (magic == "coder_fc") {
        return std::make_shared<coder_fc>(io);
    }
    if (magic == "coder_affix_match") {
        return std::make_shared<coder_affix_match>(io);
    }
    return nullptr;
}

void CoderFactory::applyLevel(coder_io* io, CoderType type, uint8_t compressLevel)
{
    if (io == nullptr || compressLevel < 1 || compressLevel > 9) {
        return;
    }
    switch (type) {
    case CoderType::BWT_CM:
        io->set_level(compressLevel);
        break;
    case CoderType::FC:
    case CoderType::AFFIX_MATCH:
    case CoderType::QUAL:
    case CoderType::FCV2:
    case CoderType::SIMPLE_RC:
    default:
        break;
    }
}

bool CoderFactory::coderSupports(CoderType type, uint32_t fileType, uint32_t fieldIdx)
{
    if (type == CoderType::FCV2) {
        return (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) && fieldIdx == (uint32_t)SAM_QUAL;
    }
    if (type == CoderType::AFFIX_MATCH) {
        /* Whether it is a candidate is decided uniformly by the config table,
         * avoiding maintaining the trial-compression scope in two places. */
        return (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) && samFieldCandidate(fieldIdx, CoderType::AFFIX_MATCH);
    }
    return true;
}

bool CoderFactory::canMake(CoderType type)
{
    return type == CoderType::BWT_CM || type == CoderType::FC ||
           type == CoderType::AFFIX_MATCH;
}
