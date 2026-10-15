/*
 * coder_factory.cpp - coder registry and encoder/decoder factory
 *
 * Implementation lives in the .cpp so concrete encoder headers do not leak to callers:
 * coder_fc.h and coder_qual.h each define a RangeCoder of the same name, which conflict
 * in one translation unit. The kCoderDescriptors table below is the registry: one row
 * per coder. New coders are added here plus their name in field_coder_config.h.
 *
 * coder_simple_rc was removed from both sides (it was lossy), and coder_lzma was
 * removed as too slow.
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

#include "coder_factory.h"

#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_arith.h"
#include "coder/coder_arith_k.h"
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
std::shared_ptr<coder> createArithK(coder_io* io)    {
    /* htscodecs rANS order-1 on the A/C/G/T bases packed four per byte. */
    return std::make_shared<coder_arith_k>(io, /*order=*/1, /*mapBases=*/true,
                                           /*useRans=*/true, /*ransPack=*/true);
}
std::shared_ptr<coder> createRans(coder_io* io)      { return std::make_shared<coder_rans>(io); }
std::shared_ptr<coder> createQcm(coder_io* io)       { return std::make_shared<coder_qcm>(io); }
std::shared_ptr<coder> createQname(coder_io* io)     { return std::make_shared<coder_qname>(io); }

/*---------------------------------------------------------------------------
 * Applicability checks
 *-------------------------------------------------------------------------*/

/* Generic byte-stream coders work on every field. */
bool supportsAny(uint32_t /*fileType*/, uint32_t /*fieldIdx*/) { return true; }

/* A FASTQ block type carries one of the two generation bits, possibly with GZIP set. */
bool isFastqType(uint32_t fileType)
{
    return (fileType & ((uint32_t)FASTQ_GEN2 | (uint32_t)FASTQ_GEN3)) != 0;
}

/* coder_fcv2 needs each record's read length and strand direction, available for SAM/BAM
   QUAL and FASTQ QUAL (FASTQ hands a forward strand; the cycle index comes from the record
   length). */
bool supportsFcv2(uint32_t fileType, uint32_t fieldIdx)
{
    if (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) {
        return fieldIdx == (uint32_t)SAM_QUAL;
    }
    return isFastqType(fileType) && fieldIdx == (uint32_t)FQ_QUAL;
}

/* Whether coder_affix_match is a candidate is decided by the config table. */
bool supportsAffix(uint32_t fileType, uint32_t fieldIdx)
{
    if (fileType == (uint32_t)SAM || fileType == (uint32_t)BAM) {
        return samFieldCandidate(fieldIdx, CoderType::AFFIX_MATCH);
    }
    return isFastqType(fileType) && fastqFieldCandidate(fieldIdx, CoderType::AFFIX_MATCH,
                                                        PBGZ_MODE_ARCHIVE);
}

/*---------------------------------------------------------------------------
 * The registry
 *-------------------------------------------------------------------------*/

/* Trial priority: coder_fc is the baseline; coder_affix_match and coder_bwt_cm only take
   the lead on a strictly smaller result. Kept here so candidate lists can be any order. */
const CoderDescriptor kCoderDescriptors[] = {
    /* type,                  profiles,               trialLineBased, trialCandidate, trialUsesLevel, minLevel, priority, levelPolicy,                 create,       supports */
    {CoderType::BWT_CM,       CODER_PROFILE_ARCHIVE,  false,          true,           true,           8,        30,       CoderLevelPolicy::SET_LEVEL, createBwtCm,  supportsAny,  false, false},
    {CoderType::FC,           CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        10,       CoderLevelPolicy::NONE,      createFc,     supportsAny,  false, true},
    {CoderType::AFFIX_MATCH,  CODER_PROFILE_ARCHIVE,  true,           true,           false,          1,        20,       CoderLevelPolicy::NONE,      createAffix,  supportsAffix, true, false},
    {CoderType::ARITH,        CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        40,       CoderLevelPolicy::SET_LEVEL, createArith,  supportsAny,  false, false},
    {CoderType::ARITH_K,      CODER_PROFILE_ARCHIVE,  false,          false,          false,          1,        41,       CoderLevelPolicy::NONE,      createArithK, supportsAny,  false, true},
    {CoderType::RANS,         CODER_PROFILE_FAST,     false,          true,           false,          1,        50,       CoderLevelPolicy::NONE,      createRans,   supportsAny,  false, false},

    /* coder_lzma was removed as too slow (see the header comment). */

    /* coder_qcm is an order-2 context model for the quality column. It is an ordinary
       byte-stream coder applicable to any field; only the QUAL candidate list offers it.
       Last in trial order (priority 60) so a tie keeps the incumbent. */
    {CoderType::QCM,          CODER_PROFILE_ARCHIVE,  false,          true,           false,          1,        60,       CoderLevelPolicy::NONE,      createQcm,    supportsAny,  false, false},

    /* coder_qname codes a whole name per line (dedicated QNAME path, see
       compressIdFieldQname). Never a trial candidate; its row lets a stream be rebuilt by
       magic alone (see makeFieldDecoder). */
    {CoderType::QNAME,        CODER_PROFILE_ARCHIVE,  false,          false,          false,          1,        0,        CoderLevelPolicy::NONE,      createQname,  supportsAny,  false, false},

    /* QUAL (coder_qual) and FCV2 (coder_fcv2) do not inherit coder and need record-level
       input, so they have no create function here; their rows only record that they are
       archive coders. */
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
    /* Reuse coderTypeToMagic so the dispatch key is defined in one place. */
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
    /* Unknown/uncreatable type: BWT_CM is the safe fallback (QUAL/FCV2 land here). */
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
    /* A coder_arith decoder must use the encoder's file-level prior, or decode diverges.
       The magic was checked by makeDecoder; the cast is safe. */
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
