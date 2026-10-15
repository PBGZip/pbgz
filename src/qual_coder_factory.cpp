/*
 * qual_coder_factory.cpp - builds the quality column's encoder/decoder from a trial verdict.
 * Its own translation unit because coder_qual.h cannot share one with coder_fc.h.
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

#include "log/logger.h"

#include "coder/coder_fcv2.h"
#include "coder/coder_qual.h"

namespace {

/* coder_qual: codes each record against its SEQ; the fallback with no verdict. */
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

    bool coversByte(uint8_t b) const override { return coder->coversByte(b); }

private:
    std::shared_ptr<coder_fcv2> coder;
};

/* A byte-stream coder (coder_bwt_cm, coder_qcm, ...): the record is a line. */
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
        /* The tier the trial settled on; defaults stand without a verdict. */
        Fcv2Cfg cfg;
        if (args.fcv2Params != nullptr) {
            cfg.cycleMax = args.fcv2Params->cycleMax;
            cfg.cycleBucket = args.fcv2Params->cycleBucket;
            cfg.prevShift = args.fcv2Params->prevShift;
            cfg.useDedup = args.fcv2Params->useDedup;
        }
        /* Turn the carried frequency rank into a monotonically decreasing weight. */
        std::vector<uint32_t> freqByByte(256, 0);
        if (args.freqTable != nullptr) {
            for (uint32_t i = 0; i < args.freqTable->size(); ++i) {
                const uint32_t b = (uint32_t)(*args.freqTable)[i].first + (uint32_t)'!';
                if (b < 256) {
                    freqByByte[b] = (uint32_t)(args.freqTable->size() - i);
                }
            }
        }
        /* Start from the trained prior when the engine has one; a load failure is benign but
           reported as none. */
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

    /* Everything else comes through the registry as an ordinary line stream. */
    std::shared_ptr<coder> streamCoder = makeEncoder(picked, io);
    applyLevel(io, picked, compressLevel);
    return std::make_shared<qual_stream_encoder>(streamCoder);
}

/*
 * ---------------------------------------------------------------------------
 * Decoding: turn the stream's magic into a decoder.
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

/* A byte-stream coder (coder_bwt_cm, coder_qcm, ...): the caller's length delimits the record. */
class qual_stream_decoder : public qual_record_decoder {
public:
    explicit qual_stream_decoder(std::shared_ptr<coder> c) : stream(std::move(c)) {}

    int32_t decode_record(uint8_t* qual, uint32_t qualLen, uint8_t*, uint32_t) override
    {
        /* UINT8_MAX is the "no split symbol" value these coders have always been given here. */
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
        /* fcv2's alphabet and frequencies travel in its stream header; an empty table suffices. */
        const std::vector<uint32_t> emptyFreq(256, 0);
        std::shared_ptr<coder_fcv2> coder;
        if (args.priorRequired) {
            /* A prior-encoded stream must start from the same snapshot; failure here is fatal. */
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
            /* Neither the current nor the pre-version byte header layout fits this stream. */
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

    /* Everything else is a registered coder, built by magic. */
    std::shared_ptr<coder> streamCoder = makeDecoder(magic, io);
    if (streamCoder == nullptr) {
        return nullptr;
    }
    return std::make_shared<qual_stream_decoder>(streamCoder);
}
