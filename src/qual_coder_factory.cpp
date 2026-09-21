/*
 * qual_coder_factory.cpp - building the quality column's encoder from a trial verdict
 *
 * The quality column is the one field whose coders are not interchangeable byte-stream
 * compressors: coder_qual codes each record against its SEQ, coder_fcv2 codes it against
 * its sequencing cycle and strand, and neither inherits coder. Everything that decides
 * which of them is used - the trial's verdict, the quality alphabet, the fcv2 context tier,
 * the trained prior - used to be assembled inside SamCodecActuator::compressQuality, which
 * left the actuator knowing the constructor signature of every coder and carrying one
 * encoder member per possibility.
 *
 * That knowledge lives here now. The actuator passes the verdict and the material (see
 * QualCoderArgs in coder_factory.h) and gets back a qual_record_encoder, whose single
 * encode_record() feeds a record whichever coder won.
 *
 * This is a translation unit of its own because coder_qual.h pulls in qual_model.h, which
 * has a RangeCoder of the same name as coder_fc.h's but a different interface; the registry
 * header deliberately includes neither (see the note in coder_factory.h), so this is the one
 * place that must include the record-level pair.
 */

#include "coder_factory.h"

#include "log/logger.h"

#include "coder/coder_fcv2.h"
#include "coder/coder_qual.h"

namespace {

/*
 * coder_qual: the original quality coder, which codes each record against its SEQ. It was
 * the only path before selection was wired in, and it stays the fallback when preprocessing
 * produced no verdict (see coderFor / samFieldDefaultCoder at the call site).
 */
class qual_seq_encoder : public qual_record_encoder {
public:
    explicit qual_seq_encoder(std::shared_ptr<coder_qual> c) : coder(std::move(c)) {}

    void encode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq, uint32_t, bool) override
    {
        coder->encode_qual_gen2(seq, qual, qualLen);
    }

    void flush() override { coder->encode_flush(); }

private:
    std::shared_ptr<coder_qual> coder;
};

/* coder_fcv2: codes each record against its sequencing cycle and strand direction. */
class qual_fcv2_encoder : public qual_record_encoder {
public:
    explicit qual_fcv2_encoder(std::shared_ptr<coder_fcv2> c) : coder(std::move(c)) {}

    void encode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq, uint32_t seqLen,
                       bool rev) override
    {
        coder->encode_record(qual, qualLen, rev, seq, seqLen);
    }

    void flush() override { coder->encode_flush(); }

private:
    std::shared_ptr<coder_fcv2> coder;
};

/*
 * A byte-stream coder (coder_bwt_cm, coder_qcm, ...): the record is a line like any other,
 * and the SEQ context and strand direction are simply not used.
 */
class qual_stream_encoder : public qual_record_encoder {
public:
    explicit qual_stream_encoder(std::shared_ptr<coder> c) : stream(std::move(c)) {}

    void encode_record(uint8_t* qual, uint32_t qualLen, uint8_t*, uint32_t, bool) override
    {
        stream->encode_line(qual, qualLen);
    }

    void flush() override { stream->encode_flush(); }

private:
    std::shared_ptr<coder> stream;
};

}  // namespace

std::shared_ptr<qual_record_encoder> CoderFactory::makeQualEncoder(CoderType picked, coder_io* io,
                                                                  const QualCoderArgs& args,
                                                                  uint8_t compressLevel)
{
    if (picked == CoderType::FCV2) {
        /* The tier the trial settled on; without a verdict the coder's own defaults stand,
           which is the cold-start behavior. */
        Fcv2Cfg cfg;
        if (args.fcv2Params != nullptr) {
            cfg.cycleMax = args.fcv2Params->cycleMax;
            cfg.cycleBucket = args.fcv2Params->cycleBucket;
            cfg.deltaMax = args.fcv2Params->deltaMax;
            cfg.deltaBucket = args.fcv2Params->deltaBucket;
            cfg.prevShift = args.fcv2Params->prevShift;
            cfg.useDelta = args.fcv2Params->useDelta;
            cfg.useDedup = args.fcv2Params->useDedup;
            cfg.useQa = args.fcv2Params->useQa;
            cfg.modelCount = args.fcv2Params->modelCount;
        }
        /*
         * coder_fcv2 keeps no counts of its own: the table only shapes its Huffman tree, so
         * the rank the preprocessing carried is turned into a monotonically decreasing weight
         * here.
         */
        std::vector<uint32_t> freqByByte(256, 0);
        if (args.freqTable != nullptr) {
            for (uint32_t i = 0; i < args.freqTable->size(); ++i) {
                const uint32_t b = (uint32_t)(*args.freqTable)[i].first + (uint32_t)'!';
                if (b < 256) {
                    freqByByte[b] = (uint32_t)(args.freqTable->size() - i);
                }
            }
        }
        /*
         * Start from the trained prior when the engine has one: the model need not relearn
         * from fixed initial values, and the gain is larger for smaller blocks. A load failure
         * is benign - the fixed-initial model stands and only the ratio suffers - but the
         * address is then reported as none (the caller drops it), so the decoding side is
         * never promised a snapshot it cannot obtain.
         */
        std::shared_ptr<coder_fcv2> coder;
        bool loaded = false;
        if (args.priorBlob != nullptr && !args.priorBlob->empty()) {
            coder = std::make_shared<coder_fcv2>(io, freqByByte, cfg, *args.priorBlob, &loaded);
        } else {
            coder = std::make_shared<coder_fcv2>(io, freqByByte, cfg);
        }
        if (args.priorLoaded != nullptr) {
            *args.priorLoaded = loaded;
        }
        return std::make_shared<qual_fcv2_encoder>(coder);
    }

    if (picked == CoderType::QUAL) {
        static const std::vector<std::pair<uint16_t, uint16_t>> kNoTable;
        const std::vector<std::pair<uint16_t, uint16_t>>& table =
            (args.freqTable != nullptr) ? *args.freqTable : kNoTable;
        return std::make_shared<qual_seq_encoder>(
            std::make_shared<coder_qual>(io, true, table));
    }

    /* Everything else takes part through the registry: build it there, apply the level the
       engine asked for, and present it as an ordinary line stream. */
    std::shared_ptr<coder> streamCoder = makeEncoder(picked, io);
    applyLevel(io, picked, compressLevel);
    return std::make_shared<qual_stream_encoder>(streamCoder);
}

/*
 * ---------------------------------------------------------------------------
 * Decoding
 *
 * The same division as above, in the other direction: what reads the stream's meta turns
 * the magic it finds into a decoder, and the actuator only has to say where the record goes.
 * ---------------------------------------------------------------------------
 */

namespace {

/* coder_qual: the record is decoded against its SEQ, which the caller has already decoded
   ahead of the quality column. */
class qual_seq_decoder : public qual_record_decoder {
public:
    explicit qual_seq_decoder(std::shared_ptr<coder_qual> c) : coder(std::move(c)) {}

    int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq, uint32_t) override
    {
        coder->decode_qual_gen2(seq, qual, qualLen);
        return 0;
    }

private:
    std::shared_ptr<coder_qual> coder;
};

/* coder_fcv2: the strand direction travels with the stream, so the only context to hand
   over here is this record's decoded bases. */
class qual_fcv2_decoder : public qual_record_decoder {
public:
    explicit qual_fcv2_decoder(std::shared_ptr<coder_fcv2> c) : coder(std::move(c)) {}

    int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t* seq, uint32_t seqLen) override
    {
        return coder->decode_record(qual, qualLen, seq, seqLen);
    }

private:
    std::shared_ptr<coder_fcv2> coder;
};

/* A byte-stream coder (coder_bwt_cm, coder_qcm, ...): the record is a line like any other,
   and the length is the one the caller asks for - the same length the encoding side fed it,
   with no delimiter in between. */
class qual_stream_decoder : public qual_record_decoder {
public:
    explicit qual_stream_decoder(std::shared_ptr<coder> c) : stream(std::move(c)) {}

    int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t*, uint32_t) override
    {
        /*
         * UINT8_MAX is the "no split symbol" value these coders have always been given on
         * this path - the record boundaries come from the caller's length, there is no
         * delimiter in the stream. It is spelled out rather than left to a default: the base
         * class declares two decode_line overloads whose defaults differ from the derived
         * ones, so an argument-less call is both ambiguous and would pick '\n'.
         */
        return (stream->decode_line(qual, qualLen, UINT8_MAX, false) < 0) ? -1 : 0;
    }

private:
    std::shared_ptr<coder> stream;
};

}  // namespace

std::shared_ptr<qual_record_decoder> CoderFactory::makeQualDecoder(const std::string& magic,
                                                                  coder_io* io,
                                                                  const QualDecoderArgs& args)
{
    if (magic == "coder_qual") {
        static const std::vector<std::pair<uint16_t, uint16_t>> kNoTable;
        const std::vector<std::pair<uint16_t, uint16_t>>& table =
            (args.freqTable != nullptr) ? *args.freqTable : kNoTable;
        return std::make_shared<qual_seq_decoder>(std::make_shared<coder_qual>(io, true, table));
    }

    if (magic == "coder_fcv2") {
        /*
         * fcv2's alphabet and per-symbol frequencies travel in its own stream header;
         * begin_decode reads them back to rebuild the Huffman tree, so an empty frequency
         * table at construction is enough.
         */
        const std::vector<uint32_t> emptyFreq(256, 0);
        std::shared_ptr<coder_fcv2> coder;
        if (args.priorRequired) {
            /*
             * If the encoder started from a prior, this side must start from the same
             * snapshot or the model diverges immediately. Unlike the benign fallback on the
             * compression side, a failure here is fatal: silently falling back to fixed
             * initial values would decode a stream of seemingly valid wrong data, far more
             * dangerous than failing outright. The address comes from the block meta rather
             * than sequential inference, so random access works as well.
             */
            if (args.priorBlob == nullptr || args.priorBlob->empty()) {
                LOG_ERROR("Qual prior at offset %lld required but unavailable",
                          (long long)args.priorAddress);
                return nullptr;
            }
            bool priorLoaded = false;
            coder = std::make_shared<coder_fcv2>(io, emptyFreq, *args.priorBlob, &priorLoaded);
            if (!priorLoaded) {
                LOG_ERROR("Qual prior at offset %lld failed to load", (long long)args.priorAddress);
                return nullptr;
            }
        } else {
            coder = std::make_shared<coder_fcv2>(io, emptyFreq);
        }
        const int32_t ret = coder->begin_decode();
        if (ret == coder_ns::CODER_ERR_UNSUPPORTED_VERSION) {
            /* Neither the current header layout nor the one from before the version byte
               fits this stream, so its layout is unknown - naming that is more useful than
               reporting it as corruption (see FCV2_STREAM_VERSION). */
            LOG_ERROR("QUAL stream fits neither the current fcv2 header layout (v%u) nor "
                      "the layout from before it: the archive was written by an "
                      "incompatible version and cannot be decoded",
                      (unsigned)FCV2_STREAM_VERSION);
            return nullptr;
        }
        if (ret != 0) {
            LOG_ERROR("fcv2 begin_decode failed");
            return nullptr;
        }
        return std::make_shared<qual_fcv2_decoder>(coder);
    }

    /* Everything else is a registered coder and is built by magic: coder_bwt_cm, coder_qcm
       and any later candidate the quality column's list gains. */
    std::shared_ptr<coder> streamCoder = makeDecoder(magic, io);
    if (streamCoder == nullptr) {
        return nullptr;
    }
    return std::make_shared<qual_stream_decoder>(streamCoder);
}
