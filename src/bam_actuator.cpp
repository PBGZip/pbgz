/*
 * bam_actuator.cpp - column-wise BAM compression without SAM text
 *
 * See bam_actuator.h. Every column is encoded by the coder that preprocessing
 * picked for the matching SAM field, so the encoder of no field changes; what
 * changes is that the bytes are the values BAM already stores instead of their
 * textual rendering.
 *
 * Decompression is not implemented yet - this actuator exists to measure what
 * the structured representation does to the compression ratio.
 */

#include <map>
#include "bam_actuator.h"

#include <string.h>

#include "bam_record.h"

#include <string>
#include <vector>
#include <algorithm>

#include "coder/coder_bwt_cm.h"
#include "coder/coder_fcv2.h"
#include "coder/coder_json.h"
#include "coder/coder_qual.h"
#include "coder_factory.h"
#include "field_coder_config.h"
#include "sam_field_rules.h"
#include "sam_info.h"
#include "log/logger.h"
#include "preprocess_info.h"
#include "profile_stats.h"
#include "utils/md5_util.h"

namespace {

const char kSeqBaseTable[] = "=ACMGRSVTWYHKDBN";

/* Append v as n little-endian bytes. */
void appendLe(std::vector<uint8_t>& out, uint32_t v, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        out.push_back((uint8_t)((v >> (8 * i)) & 0xff));
    }
}

}  // namespace

namespace {

static const char kBamRefActg4[4] = {'A', 'C', 'T', 'G'};

static inline uint32_t cigLen(uint32_t v) { return v >> 4; }
static inline uint32_t cigOp(uint32_t v) { return v & 0xf; }

/* Shared by encoder and decoder so the per-record reference decision is
   symmetric. cigar/ops may be empty for unaligned records. */
static bool seqRefEligible(uint32_t flag, int32_t refId, int32_t pos,
                           const uint32_t* ops, size_t nOps, int32_t lSeq,
                           int64_t& chrStart, int64_t& refBase,
                           uint32_t& refConsumed)
{
    chrStart = -1;
    refBase = -1;
    refConsumed = 0;
    if (refId < 0 || (flag & 0x4u) || pos < 0 || nOps == 0)
        return false;
    const int64_t start = SamInfo::getInstance().getPositionByIndex((uint32_t)refId);
    if (start < 0)
        return false;
    int64_t readLen = 0;
    int64_t rc = 0;
    for (size_t i = 0; i < nOps; ++i) {
        const uint32_t op = cigOp(ops[i]);
        const uint32_t n = cigLen(ops[i]);
        switch (op) {
        case 0: case 7: case 8: readLen += n; rc += n; break; /* M = X */
        case 1: case 4:         readLen += n; break;         /* I S */
        case 2: case 3:         rc += n; break;              /* D N */
        default: break;
        }
    }
    if (readLen != (int64_t)lSeq)
        return false;
    chrStart = start;
    refBase = start + pos;
    if (refBase < 0)
        return false;
    refConsumed = (uint32_t)rc;
    return true;
}

}  // namespace

/* Reference bases consumed by a record's CIGAR (M/=/X/D/N). */
static uint32_t bamCigRefSpan(const uint32_t* ops, size_t n)
{
    uint32_t s = 0;
    for (size_t i = 0; i < n; ++i) {
        const uint32_t op = cigOp(ops[i]);
        const uint32_t l = cigLen(ops[i]);
        switch (op) {
        case 0: case 7: case 8: case 2: case 3:
            s += l;
            break;
        default:
            break;
        }
    }
    return s;
}

/* What a BAM record has to satisfy before TLEN can be inferred, and then the shared rule itself
   (samTemplateLen in sam_field_rules.h), which the SAM actuator applies to the same column.
   mateSpan = 0 when the mate is not in this block. */
static int32_t bamTlenCompute(int64_t flag, int64_t refId, int64_t pos,
                              int64_t nextRefId, int64_t nextPos,
                              uint32_t selfSpan, uint32_t mateSpan, bool minusOne)
{
    if (!(flag & 0x1) || (flag & 0x4) || (flag & 0x8))
        return 0;
    if (refId < 0 || pos < 0 || nextPos < 0)
        return 0;
    if (!(nextRefId == refId || nextRefId >= 0))
        return 0;   /* '=' or a real reference */
    return samTemplateLen(pos, nextPos, selfSpan, mateSpan, minusOne);
}

static inline int varintLenU(uint64_t v)
{
    int n = 1;
    while (v >= 0x80) { v >>= 7; ++n; }
    return n;
}

BamCodecActuator::BamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine,
                                   Reference* pRef)
    : CodecActuator(inPtr, outPtr, engine), cols(nullptr), pRefeGene(pRef)
{
}

void BamCodecActuator::putUvarint(std::vector<uint8_t>& out, uint64_t v)
{
    while (v >= 0x80) {
        out.push_back((uint8_t)(v | 0x80));
        v >>= 7;
    }
    out.push_back((uint8_t)v);
}

void BamCodecActuator::putSvarint(std::vector<uint8_t>& out, int64_t v)
{
    const uint64_t zz = (v < 0) ? (~((uint64_t)v) << 1) | 1u : ((uint64_t)v << 1);
    putUvarint(out, zz);
}

uint32_t BamCodecActuator::pickedFor(uint32_t fieldIdx, CoderType fallback) const
{
    const PreprocessInfo* preInfo = (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
    if (preInfo != nullptr) {
        return (uint32_t)preInfo->coderFor(fieldIdx, fallback);
    }
    return (uint32_t)fallback;
}

int32_t BamCodecActuator::encodeColumn(uint32_t fieldIdx, CoderType fallback, const uint8_t* data,
                                       size_t size, Json::Value& streamMeta, const char* name,
                                       const char* mode)
{
    if (size == 0) {
        return 0;
    }
    std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), name);
    std::shared_ptr<coder> coder = makeFieldEncoder(fieldIdx, fallback, io.get(), false);
    coder->encode_line(data, (uint32_t)size);
    coder->encode_flush();
    if (io->err != coder_io::IO_OK) {
        LOG_ERROR("BAM column %s: encode overflow", name);
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + io->data_len);

    Json::Value m;
    m["field"] = fieldIdx;
    m["name"] = name;
    m["mode"] = mode;
    m["srclen"] = (Json::Value::UInt)size;
    m["dstlen"] = (Json::Value::UInt)io->data_len;
    m["coder"] = io->meta;
    streamMeta.append(m);
    return (int32_t)io->data_len;
}


int32_t BamCodecActuator::encodeQnameColumn(const uint8_t* qnameStream, size_t qnameStreamSize,
                                            uint32_t nRecords, Json::Value& streamMeta)
{
    (void)nRecords;
    /* QNAME is concatenated with '\n' separators and coded as one whole block. */
    return encodeRawBwtColumn(SAM_QNAME, qnameStream, qnameStreamSize, streamMeta, "QNAME");
}

int32_t BamCodecActuator::encodeRawBwtColumn(uint32_t fieldIdx, const uint8_t* data, size_t size,
                                             Json::Value& streamMeta, const char* name)
{
    std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), name);
    std::shared_ptr<coder> coder = makeFieldEncoder(fieldIdx, CoderType::BWT_CM, io.get(), false);
    coder->encode_line(data, (uint32_t)size);
    coder->encode_flush();
    if (io->err != coder_io::IO_OK) {
        LOG_ERROR("BAM column %s: encode overflow", name);
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + io->data_len);

    Json::Value m;
    m["field"] = fieldIdx;
    m["name"] = name;
    m["mode"] = "block";
    m["srclen"] = (Json::Value::UInt)size;
    m["dstlen"] = (Json::Value::UInt)io->data_len;
    m["coder"] = io->meta;
    streamMeta.append(m);
    return (int32_t)io->data_len;
}

int32_t BamCodecActuator::encodePosColumn(Json::Value& streamMeta)
{
    /* POS: zigzag-varint delta against the previous record, baseline reset on a
       chromosome switch - the same transform the SAM path applies. */
    std::vector<uint8_t> buf;
    buf.reserve(cols->pos.size() * 2);
    int64_t prev = 0;
    int32_t prevRef = -1;
    for (size_t i = 0; i < cols->pos.size(); ++i) {
        if (cols->refId[i] != prevRef) {
            prev = 0;
            prevRef = cols->refId[i];
        }
        const int64_t d = (int64_t)cols->pos[i] - prev;
        putSvarint(buf, d);
        prev = cols->pos[i];
    }
    return encodeColumn(SAM_POS, CoderType::BWT_CM, buf.data(), buf.size(), streamMeta, "POS", "varint-delta");
}

int32_t BamCodecActuator::encodePnextColumn(Json::Value& streamMeta)
{
    std::vector<uint8_t> buf;
    buf.reserve(cols->nextPos.size() * 2);
    int64_t prev = 0;
    for (size_t i = 0; i < cols->nextPos.size(); ++i) {
        const int64_t d = (int64_t)cols->nextPos[i] - prev;
        putSvarint(buf, d);
        prev = cols->nextPos[i];
    }
    return encodeColumn(SAM_PNEXT, CoderType::BWT_CM, buf.data(), buf.size(), streamMeta, "PNEXT", "varint-delta");
}

int32_t BamCodecActuator::encodeTlenColumn(Json::Value& streamMeta)
{
    /* Decoder-side TLEN rebuild (same rule as the archive path): TLEN is only
       stored as exceptions. Paired sorted reads whose mate is in the same
       block are reconstructed from POS/PNEXT/CIGAR under the per-block
       convention (minusOne / plain); everything else is an exception. If the
       exceptions are not smaller than plain varints, fall back to plain. */
    std::vector<uint32_t> span(cols->nRecords, 0);
    {
        size_t cOff = 0;
        for (uint32_t r = 0; r < cols->nRecords; ++r) {
            span[r] = bamCigRefSpan(cols->cigar.data() + cOff, cols->cigarLen[r]);
            cOff += cols->cigarLen[r];
        }
    }
    std::map<std::pair<int64_t, int64_t>, uint32_t> mateMap;
    for (uint32_t r = 0; r < cols->nRecords; ++r) {
        const int64_t p = cols->pos[r];
        const int64_t np = cols->nextPos[r];
        if (p >= 0 && np >= 0)
            mateMap[{p, np}] = r;
    }
    auto mateOf = [&](uint32_t r) -> uint32_t {
        auto it = mateMap.find({cols->nextPos[r], cols->pos[r]});
        return (it != mateMap.end()) ? it->second : cols->nRecords;
    };
    auto computeFor = [&](uint32_t r, bool minusOne) -> int32_t {
        const uint32_t m = mateOf(r);
        const uint32_t ms = (m < cols->nRecords) ? span[m] : 0;
        return bamTlenCompute(cols->flag[r], cols->refId[r], cols->pos[r],
                              cols->nextRefId[r], cols->nextPos[r],
                              span[r], ms, minusOne);
    };
    uint64_t convM1 = 0, convP = 0;
    for (uint32_t r = 0; r < cols->nRecords; ++r) {
        const int32_t v = cols->tlen[r];
        if (computeFor(r, true) == v) convM1++;
        if (computeFor(r, false) == v) convP++;
    }
    const bool minusOne = (convM1 >= convP);

    std::vector<uint32_t> excLines;
    std::vector<int32_t> excVals;
    for (uint32_t r = 0; r < cols->nRecords; ++r) {
        const int32_t v = cols->tlen[r];
        if (computeFor(r, minusOne) != v) {
            excLines.push_back(r);
            excVals.push_back(v);
        }
    }
    /* payload: [u8 conv][ (varint lineDelta)(svarint value) ]* */
    std::vector<uint8_t> opt;
    size_t optBytes = 1;
    uint32_t prevLine = 0;
    for (size_t i = 0; i < excLines.size(); ++i) {
        optBytes += (size_t)varintLenU(excLines[i] - prevLine);
        prevLine = excLines[i];
        uint64_t zz = (excVals[i] < 0) ? (~((uint64_t)excVals[i]) << 1) | 1u
                                       : ((uint64_t)excVals[i] << 1);
        optBytes += (size_t)varintLenU(zz);
    }
    size_t plainBytes = 0;
    for (uint32_t r = 0; r < cols->nRecords; ++r) {
        const int64_t v = cols->tlen[r];
        const uint64_t zz = (v < 0) ? (~((uint64_t)v) << 1) | 1u : ((uint64_t)v << 1);
        plainBytes += (size_t)varintLenU(zz);
    }
    if (optBytes < plainBytes) {
        opt.push_back((uint8_t)(minusOne ? 1 : 0));
        prevLine = 0;
        for (size_t i = 0; i < excLines.size(); ++i) {
            putUvarint(opt, excLines[i] - prevLine);
            prevLine = excLines[i];
            putSvarint(opt, excVals[i]);
        }
        return encodeColumn(SAM_TLEN, CoderType::BWT_CM, opt.data(), opt.size(),
                            streamMeta, "TLEN_R", "rebuild");
    }
    std::vector<uint8_t> buf;
    buf.reserve(cols->tlen.size() * 2);
    for (size_t i = 0; i < cols->tlen.size(); ++i)
        putSvarint(buf, cols->tlen[i]);
    return encodeColumn(SAM_TLEN, CoderType::BWT_CM, buf.data(), buf.size(),
                        streamMeta, "TLEN", "varint");
}

int32_t BamCodecActuator::encodeCigarColumn(Json::Value& streamMeta)
{
    /* One varint per BAM cigar op, plus the per-record op count as its own
       column (the payload alone would not be splittable). */
    std::vector<uint8_t> buf;
    buf.reserve(cols->cigar.size() * 2);
    for (size_t i = 0; i < cols->cigar.size(); ++i) {
        putUvarint(buf, cols->cigar[i]);
    }

    std::vector<uint8_t> lens;
    lens.reserve(cols->cigarLen.size() * 2);
    for (size_t i = 0; i < cols->cigarLen.size(); ++i) {
        putUvarint(lens, cols->cigarLen[i]);
    }
    if (encodeColumn(SAM_FIELD_COUNT, CoderType::BWT_CM, lens.data(), lens.size(), streamMeta, "CIGAR_N", "varint") < 0) {
        return -1;
    }
    return encodeColumn(SAM_CIGAR, CoderType::BWT_CM, buf.data(), buf.size(), streamMeta, "CIGAR", "varint");
}

int32_t BamCodecActuator::encodeSeqColumnRef(Json::Value& streamMeta)
{
    const bool noRef = (pRefeGene == nullptr);
    const int64_t refMaxBase = noRef ? 0 : ((int64_t)pRefeGene->getSquashLength() << 2);

    /* pass 1: does any record qualify for the reference path? */
    bool anyRef = false;
    {
        size_t cOff = 0;
        for (uint32_t r = 0; r < cols->nRecords && !anyRef; ++r) {
            const uint32_t nOps = cols->cigarLen[r];
            const uint32_t len = (uint32_t)cols->lSeq[r];
            if (len == 0 || noRef) { cOff += nOps; continue; }
            int64_t cs, rb;
            uint32_t rc;
            if (seqRefEligible(cols->flag[r], cols->refId[r], cols->pos[r],
                               cols->cigar.data() + cOff, nOps, (int32_t)len,
                               cs, rb, rc)) {
                if (rb + (int64_t)rc + 8 <= refMaxBase)
                    anyRef = true;
            }
            cOff += nOps;
        }
    }
    if (!anyRef) {
        return encodeRawBwtColumn(SAM_SEQ, cols->seq.data(), cols->seq.size(),
                                  streamMeta, "SEQ");
    }

    std::vector<uint8_t> diff;
    diff.reserve(cols->seq.size());
    std::vector<uint8_t> refBuf;
    size_t sOff = 0;
    size_t cOff = 0;
    for (uint32_t r = 0; r < cols->nRecords; ++r) {
        const uint32_t len = (uint32_t)cols->lSeq[r];
        const uint32_t nOps = cols->cigarLen[r];
        const uint32_t* ops = cols->cigar.data() + cOff;
        int64_t chrStart = -1, refBase = -1;
        uint32_t refConsumed = 0;
        bool canRef = false;
        if (len > 0 && !noRef) {
            canRef = seqRefEligible(cols->flag[r], cols->refId[r], cols->pos[r],
                                    ops, nOps, (int32_t)len,
                                    chrStart, refBase, refConsumed);
            if (canRef && (refBase + (int64_t)refConsumed + 8 > refMaxBase))
                canRef = false;
            if (canRef) {
                /*
                 * Mark the reference span this record matched against. Without
                 * it the packed reference is sanitised to zero on the way out
                 * (sanitizeRefSquash keeps only matched bytes), so decompression
                 * would rebuild SEQ from an all-zero reference - i.e. all 'A'.
                 * Same contract as the SAM path (SamCodecActuator calls
                 * updateMatchedGene with the CIGAR reference-consumed span).
                 */
                pRefeGene->updateMatchedGene((uint64_t)refBase, refConsumed);
            }
        }
        uint32_t rd = 0;
        int64_t rl = canRef ? refBase : 0;
        for (uint32_t oi = 0; oi < nOps; ++oi) {
            const uint32_t op = cigOp(ops[oi]);
            uint32_t n = cigLen(ops[oi]);
            switch (op) {
            case 0: case 7: case 8: { /* M = X */
                if (rd + n > len) n = len - rd;
                if (canRef) {
                    if (refBuf.size() < (size_t)n + 4) refBuf.resize(n + 4);
                    pRefeGene->getStretch2Bits1Char(refBuf.data(), n, (uint64_t)rl);
                    for (uint32_t k = 0; k < n; ++k) {
                        const uint8_t ch = cols->seq[sOff + rd + k];
                        const uint8_t r2 = refBuf[k] & 3;
                        diff.push_back((ch == (uint8_t)kBamRefActg4[r2]) ? 0 : ch);
                    }
                    rl += n;
                } else {
                    for (uint32_t k = 0; k < n; ++k)
                        diff.push_back(cols->seq[sOff + rd + k]);
                }
                rd += n;
                break;
            }
            case 1: case 4: { /* I S */
                if (rd + n > len) n = len - rd;
                for (uint32_t k = 0; k < n; ++k)
                    diff.push_back(cols->seq[sOff + rd + k]);
                rd += n;
                break;
            }
            case 2: case 3: /* D N */
                if (canRef) rl += n;
                break;
            default:
                break;
            }
        }
        while (rd < len)
            diff.push_back(cols->seq[sOff + rd++]);
        sOff += len;
        cOff += nOps;
    }
    return encodeColumn(SAM_SEQ, CoderType::BWT_CM, diff.data(), diff.size(),
                        streamMeta, "SEQR", "ref");
}

int32_t BamCodecActuator::encodeQualColumn(Json::Value& streamMeta)
{
    /* fast mode codes QUAL as one whole column in payload order. Missing
       records were expanded to '*' during parsing, so lSeq sums to qual.size();
       the decoder splits it back per record and re-detects missing quals. */
    return encodeColumn(SAM_QUAL, CoderType::BWT_CM, cols->qual.data(), cols->qual.size(),
                        streamMeta, "QUAL", "block");
}

int32_t BamCodecActuator::compress()
{
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter for BAM compression");
        return -1;
    }

    cols = inBlockPtr->getBamColumns().get();
    if (cols == nullptr || cols->nRecords == 0) {
        LOG_ERROR("BamActuator: block %ld carries no structured columns", (long)inBlockPtr->getBlockId());
        return -1;
    }

    /*
     * Fast mode: the reader thread only streamed the raw BAM records
     * ([u32 size][record]...) into the block buffer; the columns are parsed
     * here, on the worker that compresses this block. Every worker parses its
     * own block concurrently, so record->column work does not serialise on the
     * reader (where it was the second-largest cost after inflate).
     */
    if (!cols->columnsBuilt) {
        PBGZ_PROF_SCOPE(pbgzprof::CODER_PARSE);
        const uint8_t* stream = inBlockPtr->getBuffer() + cols->recordStreamOffset;
        const size_t streamLen = (size_t)inBlockPtr->getDataLen() - cols->recordStreamOffset;
        if (0 != bamcol::parseRecordStream(stream, streamLen, *cols)) {
            LOG_ERROR("BamActuator: block %ld failed to parse its raw record stream",
                      (long)inBlockPtr->getBlockId());
            return -1;
        }
    }

    /* The output buffer is sized from the raw block; structured columns are
       smaller than the SAM text they replace, so this is a safe upper bound. */
    const size_t need = (size_t)inBlockPtr->getDataLen() + 4096;
    if (outBlockPtr->getBufferSize() < need) {
        if (0 != outBlockPtr->ensureCapacity(need)) {
            LOG_ERROR("BamActuator: preallocate output buffer failed");
            return -1;
        }
    }

    Json::Value streamMeta;
    uint32_t totalSrc = 0;

    /* QNAME: NUL-free by construction, so a '\n' separator makes the stream
       self-delimiting without a separate length column. */
    std::vector<uint8_t> qnameStream;
    qnameStream.reserve(cols->qname.size() + cols->nRecords);
    size_t qoff = 0;
    for (uint32_t i = 0; i < cols->nRecords; ++i) {
        qnameStream.insert(qnameStream.end(), cols->qname.begin() + qoff,
                           cols->qname.begin() + qoff + cols->qnameLen[i]);
        qnameStream.push_back('\n');
        qoff += cols->qnameLen[i];
    }
    totalSrc += (uint32_t)qnameStream.size();
    if (encodeQnameColumn(qnameStream.data(), qnameStream.size(), cols->nRecords, streamMeta) < 0) {
        return -1;
    }

    /* FLAG: uint16 stream. */
    std::vector<uint8_t> flagStream;
    flagStream.reserve(cols->flag.size() * 2);
    for (size_t i = 0; i < cols->flag.size(); ++i) {
        appendLe(flagStream, cols->flag[i], 2);
    }
    totalSrc += (uint32_t)flagStream.size();
    if (encodeColumn(SAM_FLAG, CoderType::BWT_CM, flagStream.data(), flagStream.size(), streamMeta,
                     "FLAG", "u16") < 0) {
        return -1;
    }

    /* RNAME / RNEXT: reference ids, delta-varint (records are sorted). */
    std::vector<uint8_t> refStream, nextRefStream;
    int32_t prevRef = 0, prevNextRef = 0;
    for (size_t i = 0; i < cols->refId.size(); ++i) {
        putSvarint(refStream, (int64_t)cols->refId[i] - prevRef);
        prevRef = cols->refId[i];
        putSvarint(nextRefStream, (int64_t)cols->nextRefId[i] - prevNextRef);
        prevNextRef = cols->nextRefId[i];
    }
    totalSrc += (uint32_t)(refStream.size() + nextRefStream.size());
    if (encodeColumn(SAM_RNAME, CoderType::BWT_CM, refStream.data(), refStream.size(), streamMeta,
                     "RNAME", "varint-delta") < 0) {
        return -1;
    }
    if (encodeColumn(SAM_RNEXT, CoderType::BWT_CM, nextRefStream.data(), nextRefStream.size(), streamMeta,
                     "RNEXT", "varint-delta") < 0) {
        return -1;
    }

    if (encodePosColumn(streamMeta) < 0) {
        return -1;
    }
    totalSrc += (uint32_t)(cols->pos.size() * 4);

    /* MAPQ: uint8 stream. */
    totalSrc += (uint32_t)cols->mapq.size();
    if (encodeColumn(SAM_MAPQ, CoderType::BWT_CM, cols->mapq.data(), cols->mapq.size(), streamMeta,
                     "MAPQ", "u8") < 0) {
        return -1;
    }

    /* Read length is needed to split SEQ/QUAL back into records. */
    std::vector<uint8_t> lseqStream;
    for (size_t i = 0; i < cols->lSeq.size(); ++i) {
        putUvarint(lseqStream, (uint64_t)cols->lSeq[i]);
    }
    if (encodeColumn(SAM_FIELD_COUNT, CoderType::BWT_CM, lseqStream.data(), lseqStream.size(), streamMeta,
                     "L_SEQ", "varint") < 0) {
        return -1;
    }

    if (encodeCigarColumn(streamMeta) < 0) {
        return -1;
    }
    if (encodePnextColumn(streamMeta) < 0) {
        return -1;
    }
    if (encodeTlenColumn(streamMeta) < 0) {
        return -1;
    }

    /* SEQ: encode against the reference when available (per-base diff stream),
       otherwise the whole text column. */
    totalSrc += (uint32_t)cols->seq.size();
    if (encodeSeqColumnRef(streamMeta) < 0) {
        return -1;
    }

    if (encodeQualColumn(streamMeta) < 0) {
        return -1;
    }
    totalSrc += (uint32_t)cols->qual.size();

    /* Aux tags: payload + per-record length. */
    std::vector<uint8_t> tagLenStream;
    for (size_t i = 0; i < cols->tagsLen.size(); ++i) {
        putUvarint(tagLenStream, cols->tagsLen[i]);
    }
    if (encodeColumn(SAM_FIELD_COUNT, CoderType::BWT_CM, tagLenStream.data(), tagLenStream.size(), streamMeta,
                     "TAGS_N", "varint") < 0) {
        return -1;
    }
    totalSrc += (uint32_t)cols->tags.size();
    if (encodeColumn(SAM_FIELD_COUNT, CoderType::BWT_CM, cols->tags.data(), cols->tags.size(), streamMeta,
                     "TAGS", "bytes") < 0) {
        return -1;
    }

    std::string md5;
    if (inBlockPtr->hasMd5()) {
        md5 = inBlockPtr->getMd5();
    } else {
        calcMd5sum(md5, inBlockPtr->getBuffer(), (uint32_t)inBlockPtr->getDataLen());
    }

    Json::Value bamMeta;
    bamMeta["records"] = cols->nRecords;
    bamMeta["streams"] = streamMeta;
    meta["bam"] = bamMeta;
    meta["md5"] = md5;

    coder_json metaCoder;
    const int32_t metaLen = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
    if (metaLen <= 0) {
        LOG_ERROR("BamActuator: encode meta failed");
        return -1;
    }
    outBlockPtr->setMetaLen(metaLen);
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());

    return 0;
}


namespace {
/* read one LEB128 varint */
inline bool readVarint(const uint8_t*& p, const uint8_t* end, uint64_t& out)
{
    out = 0;
    unsigned shift = 0;
    while (p < end && shift < 64) {
        const uint8_t b = *p++;
        out |= (uint64_t)(b & 0x7f) << shift;
        if ((b & 0x80) == 0) {
            return true;
        }
        shift += 7;
    }
    return false;
}

/* zig-zag decode of a signed varint */
inline int64_t readSvarint(const uint8_t*& p, const uint8_t* end, bool& ok)
{
    uint64_t zz = 0;
    ok = readVarint(p, end, zz);
    return (int64_t)((zz >> 1) ^ -(int64_t)(zz & 1));
}

/* one BAM cigar op (len<<4 | code) -> "100M" */
inline void cigarOpText(uint32_t op, std::string& out)
{
    static const char* kOps = "MIDNSHP=XB";
    const uint32_t len = op >> 4;
    const unsigned code = op & 0xf;
    out += std::to_string(len);
    out += (code < 11) ? kOps[code] : 'M';
}

/* Render raw BAM aux bytes (2-char tag, type, payload...) into SAM tag text
   appended to 'line', mirroring BamBlockReader::appendBamAux. */
void renderBamTags(const uint8_t* p, const uint8_t* end, std::string& line)
{
    while (p + 3 <= end) {
        char tag[3] = {(char)p[0], (char)p[1], 0};
        const char type = (char)p[2];
        p += 3;
        line += "\t";
        line += tag;
        line += ":";
        line += 'i';            /* SAM integer subtypes are always printed as 'i' */
        line += ":";
        switch (type) {
        case 'A':
            line[line.size() - 2] = 'A';
            line += (char)*p;
            p += 1;
            break;
        case 'c':
            line += std::to_string((int)(int8_t)*p);
            p += 1;
            break;
        case 'C':
            line += std::to_string((int)*p);
            p += 1;
            break;
        case 's': { int16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
        case 'S': { uint16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
        case 'i': { int32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
        case 'I': { uint32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
        case 'f': { float v; memcpy(&v, p, 4); char bf[64]; line[line.size() - 2] = 'f'; snprintf(bf, sizeof(bf), "%g", (double)v); line += bf; p += 4; } break;
        case 'Z':
        case 'H': {
            const size_t n = strnlen((const char*)p, end - p);
            line[line.size() - 2] = type;
            line.append((const char*)p, n);
            p += n + 1;
            break;
        }
        case 'B': {
            const char sub = (char)*p++;
            int32_t count = 0;
            if (p + 4 <= end) {
                memcpy(&count, p, 4);
            }
            p += 4;
            line[line.size() - 2] = 'B';
            line += sub;
            for (int32_t i = 0; i < count && p + 1 <= end; ++i) {
                line += ',';
                switch (sub) {
                case 'c': line += std::to_string((int)(int8_t)*p); p += 1; break;
                case 'C': line += std::to_string((int)*p); p += 1; break;
                case 's': { int16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
                case 'S': { uint16_t v; memcpy(&v, p, 2); line += std::to_string(v); p += 2; } break;
                case 'i': { int32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
                case 'I': { uint32_t v; memcpy(&v, p, 4); line += std::to_string(v); p += 4; } break;
                case 'f': { float v; memcpy(&v, p, 4); char bf[64]; snprintf(bf, sizeof(bf), "%g", (double)v); line += bf; p += 4; } break;
                default: return;
                }
            }
            break;
        }
        default:
            return;
        }
    }
}
}  // namespace

int32_t BamCodecActuator::decompress()
{
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter for BAM decompression");
        return -1;
    }

    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), (int32_t)inBlockPtr->getMetaLen(), meta);
    if (!meta.isMember("bam")) {
        LOG_ERROR("BamCodecActuator::decompress: not a structured BAM block");
        return -1;
    }
    const Json::Value& bamMeta = meta["bam"];
    const Json::Value& streams = bamMeta["streams"];
    const uint32_t nRec = bamMeta["records"].asUInt();
    if (nRec == 0) {
        return 0;
    }

    /* chromosome name table: header decoding populated SamInfo from @SQ lines in
       the same order BAM stores the references. */
    std::vector<std::string> chrNames;
    {
        const std::vector<ChromosomeInfo>& all =
            SamInfo::getInstance().getAllChromosomeInfo();
        chrNames.reserve(all.size());
        for (size_t i = 0; i < all.size(); ++i) {
            chrNames.push_back(all[i].name);
        }
    }

    /*
     * `-b` output: turn the records into BAM binary right here, on the
     * decompression worker that owns this block, instead of emitting SAM text
     * for the single writer thread to re-parse. That conversion used to
     * dominate `-b` (see BamWriter). The name->id map mirrors the @SQ table
     * BamWriter builds from the header block, so both paths resolve RNAME and
     * RNEXT identically.
     */
    const bool bamOut = (pbgzEngine != nullptr) && pbgzEngine->getParameter().isDecToBam;
    std::map<std::string, int32_t> refIndex;
    if (bamOut) {
        for (size_t i = 0; i < chrNames.size(); ++i) {
            refIndex[chrNames[i]] = (int32_t)i;
        }
    }

    const uint8_t* dataBase = inBlockPtr->getBuffer();

    /* per-record containers */
    std::vector<std::string> qname(nRec);
    std::vector<uint16_t> flag(nRec);
    std::vector<int32_t> refId(nRec);
    std::vector<uint32_t> pos(nRec);
    std::vector<uint8_t> mapq(nRec);
    std::vector<int32_t> nextRefId(nRec);
    std::vector<int32_t> nextPos(nRec);
    std::vector<int32_t> tlen(nRec);
    std::vector<uint32_t> lSeq(nRec);
    std::vector<std::vector<uint32_t>> cigar(nRec);
    std::vector<uint32_t> tagsLen(nRec);

    /* TLEN rebuild state (TLEN_R column) */
    bool tlenHasConv = false;
    bool tlenConvMinusOne = false;
    std::vector<uint32_t> tlenExcLine;
    std::vector<int32_t> tlenExcVal;

    std::vector<uint8_t> seqAll;   /* concatenated base text */
    std::vector<uint32_t> seqOff(nRec + 1, 0);
    std::vector<uint8_t> qualAll;  /* concatenated quality text */
    std::vector<uint8_t> tagAll;   /* concatenated raw aux bytes */

    size_t off = 0;
    for (Json::ArrayIndex si = 0; si < streams.size(); ++si) {
        const Json::Value& sm = streams[si];
        const std::string name = sm["name"].asString();
        const size_t dst = sm["dstlen"].asUInt64();
        const size_t src = sm["srclen"].asUInt64();
        const std::string magic = sm["coder"]["magic"].asString();
        if (dst == 0) {
            continue;
        }
        std::shared_ptr<coder_io> io = std::make_shared<coder_io>(
            (uint8_t*)dataBase + off, (int32_t)dst);
        off += dst;


        std::shared_ptr<coder> dec = CoderFactory::makeDecoder(magic, io.get());
        if (dec == nullptr) {
            LOG_ERROR("BamCodecActuator::decompress: no decoder for magic %s", magic.c_str());
            return -1;
        }

        /* Whole-block columns: encoded with a single encode_line, decoded in one shot. */
        std::vector<uint8_t> raw(src);
        const int32_t got = dec->decode_line(raw.data(), (uint32_t)raw.size(), UINT8_MAX, false);
        if (got < 0 || (size_t)got != src) {
            LOG_ERROR("BamCodecActuator::decompress: column %s decode failed (got %d of %llu)",
                      name.c_str(), got, (unsigned long long)src);
            return -1;
        }
        const uint8_t* p = raw.data();
        const uint8_t* end = raw.data() + raw.size();

        if (name == "QNAME") {
            size_t line = 0;
            size_t seg = 0;
            for (size_t i = 0; i < raw.size() && line < nRec; ++i) {
                if (raw[i] == '\n') {
                    qname[line].assign((const char*)raw.data() + seg, i - seg);
                    seg = i + 1;
                    ++line;
                }
            }
            if (line < nRec) {
                qname[line].assign((const char*)raw.data() + seg, raw.size() - seg);
            }
        } else if (name == "FLAG") {
            for (uint32_t r = 0; r < nRec; ++r) {
                flag[r] = (uint16_t)(raw[r * 2] | (raw[r * 2 + 1] << 8));
            }
        } else if (name == "RNAME") {
            int64_t prev = 0;
            for (uint32_t r = 0; r < nRec; ++r) {
                bool ok = false;
                prev += readSvarint(p, end, ok);
                if (!ok) { return -1; }
                refId[r] = (int32_t)prev;
            }
        } else if (name == "RNEXT") {
            int64_t prev = 0;
            for (uint32_t r = 0; r < nRec; ++r) {
                bool ok = false;
                prev += readSvarint(p, end, ok);
                if (!ok) { return -1; }
                nextRefId[r] = (int32_t)prev;
            }
        } else if (name == "POS") {
            int64_t prev = 0;
            int32_t prevRef = -1;
            for (uint32_t r = 0; r < nRec; ++r) {
                if (refId[r] != prevRef) {
                    prev = 0;
                    prevRef = refId[r];
                }
                bool ok = false;
                prev += readSvarint(p, end, ok);
                if (!ok) { return -1; }
                pos[r] = (uint32_t)prev;
            }
        } else if (name == "MAPQ") {
            for (uint32_t r = 0; r < nRec; ++r) {
                mapq[r] = raw[r];
            }
        } else if (name == "L_SEQ") {
            for (uint32_t r = 0; r < nRec; ++r) {
                uint64_t v = 0;
                if (!readVarint(p, end, v)) { return -1; }
                lSeq[r] = (uint32_t)v;
            }
        } else if (name == "CIGAR_N") {
            for (uint32_t r = 0; r < nRec; ++r) {
                uint64_t v = 0;
                if (!readVarint(p, end, v)) { return -1; }
                cigar[r].resize((size_t)v);
            }
        } else if (name == "CIGAR") {
            uint64_t v = 0;
            size_t used = 0;
            for (uint32_t r = 0; r < nRec; ++r) {
                for (size_t k = 0; k < cigar[r].size(); ++k) {
                    if (!readVarint(p, end, v)) { return -1; }
                    cigar[r][k] = (uint32_t)v;
                    used++;
                }
            }
            (void)used;
        } else if (name == "PNEXT") {
            int64_t prev = 0;
            for (uint32_t r = 0; r < nRec; ++r) {
                bool ok = false;
                prev += readSvarint(p, end, ok);
                if (!ok) { return -1; }
                nextPos[r] = (int32_t)prev;
            }
        } else if (name == "TLEN") {
            for (uint32_t r = 0; r < nRec; ++r) {
                bool ok = false;
                const int64_t v = readSvarint(p, end, ok);
                if (!ok) { return -1; }
                tlen[r] = (int32_t)v;
            }
        } else if (name == "TLEN_R") {
            if (raw.empty()) { return -1; }
            tlenConvMinusOne = (raw[0] != 0);
            tlenHasConv = true;
            size_t xi = 1;
            uint64_t d = 0;
            uint32_t prev = 0;
            while (xi < raw.size()) {
                bool ok = false;
                d = 0;
                /* varint (line delta) */
                uint32_t shift = 0;
                while (xi < raw.size()) {
                    const uint8_t c = raw[xi++];
                    d |= (uint64_t)(c & 0x7f) << shift;
                    if (!(c & 0x80)) break;
                    shift += 7;
                }
                prev += (uint32_t)d;
                ok = false;
                int64_t val = 0;
                shift = 0;
                while (xi < raw.size()) {
                    const uint8_t c = raw[xi++];
                    val |= (int64_t)(c & 0x7f) << shift;
                    if (!(c & 0x80)) { ok = true; break; }
                    shift += 7;
                }
                if (!ok) { return -1; }
                const uint64_t zz = (uint64_t)val;
                const int32_t v = (int32_t)((zz & 1) ? (~(zz >> 1)) : (zz >> 1));
                tlenExcLine.push_back(prev);
                tlenExcVal.push_back(v);
            }
        } else if (name == "QUAL") {
            qualAll.assign(raw.begin(), raw.end());
        } else if (name == "SEQ") {
            seqAll.assign(raw.begin(), raw.end());
            for (uint32_t r = 0; r < nRec; ++r) {
                seqOff[r + 1] = seqOff[r] + lSeq[r];
            }
        } else if (name == "SEQR") {
            if (pRefeGene == nullptr) {
                LOG_ERROR("BamCodecActuator::decompress: SEQR needs the reference genome");
                return -1;
            }
            const int64_t refMaxBase = (int64_t)pRefeGene->getSquashLength() << 2;
            seqAll.clear();
            seqAll.reserve(raw.size());
            std::vector<uint8_t> refBuf;
            size_t dfi = 0;
            size_t cOff = 0;
            for (uint32_t r = 0; r < nRec; ++r) {
                const uint32_t len = (uint32_t)lSeq[r];
                seqOff[r + 1] = seqOff[r] + len;
                if (len == 0)
                    continue;
                const uint32_t nOps = (uint32_t)cigar[r].size();
                int64_t chrStart = -1, refBase = -1;
                uint32_t refConsumed = 0;
                bool canRef = seqRefEligible(flag[r], refId[r], pos[r],
                                             cigar[r].data(), nOps, (int32_t)len,
                                             chrStart, refBase, refConsumed);
                if (canRef && (refBase + (int64_t)refConsumed + 8 > refMaxBase))
                    canRef = false;
                uint32_t rd = 0;
                int64_t rl = canRef ? refBase : 0;
                for (uint32_t oi = 0; oi < nOps; ++oi) {
                    const uint32_t op = cigOp(cigar[r][oi]);
                    uint32_t n = cigLen(cigar[r][oi]);
                    switch (op) {
                    case 0: case 7: case 8: {
                        if (rd + n > len) n = len - rd;
                        if (canRef) {
                            if (refBuf.size() < (size_t)n + 4) refBuf.resize(n + 4);
                            pRefeGene->getStretch2Bits1Char(refBuf.data(), n, (uint64_t)rl);
                            for (uint32_t k = 0; k < n; ++k) {
                                const uint8_t d = raw[dfi++];
                                seqAll.push_back(d ? d : (uint8_t)kBamRefActg4[refBuf[k] & 3]);
                            }
                            rl += n;
                        } else {
                            for (uint32_t k = 0; k < n; ++k)
                                seqAll.push_back(raw[dfi++]);
                        }
                        rd += n;
                        break;
                    }
                    case 1: case 4: {
                        if (rd + n > len) n = len - rd;
                        for (uint32_t k = 0; k < n; ++k)
                            seqAll.push_back(raw[dfi++]);
                        rd += n;
                        break;
                    }
                    case 2: case 3:
                        if (canRef) rl += n;
                        break;
                    default:
                        break;
                    }
                }
                while (rd < len) {
                    /* The encoder's tail loop advances the read offset; the
                       decoder must too, or it spins on the same base until the
                       diff stream runs out. */
                    if (dfi >= raw.size()) {
                        LOG_ERROR("BamCodecActuator::decompress: SEQR diff stream exhausted");
                        return -1;
                    }
                    seqAll.push_back(raw[dfi++]);
                    ++rd;
                }
                cOff += nOps;
            }
            if (dfi != raw.size()) {
                LOG_ERROR("BamCodecActuator::decompress: SEQR stream length mismatch");
                return -1;
            }
        } else if (name == "TAGS_N") {
            for (uint32_t r = 0; r < nRec; ++r) {
                uint64_t v = 0;
                if (!readVarint(p, end, v)) { return -1; }
                tagsLen[r] = (uint32_t)v;
            }
        } else if (name == "TAGS") {
            tagAll.assign(raw.begin(), raw.end());
        } else {
            LOG_ERROR("BamCodecActuator::decompress: unknown column %s", name.c_str());
            return -1;
        }
    }

    /* Decoder-side TLEN rebuild (when the TLEN_R column is present). */
    if (!tlenExcLine.empty() || tlenHasConv) {
        std::vector<uint32_t> span(nRec, 0);
        for (uint32_t r = 0; r < nRec; ++r) {
            span[r] = bamCigRefSpan(cigar[r].data(), cigar[r].size());
        }
        std::map<std::pair<int64_t, int64_t>, uint32_t> mateMap;
        for (uint32_t r = 0; r < nRec; ++r) {
            const int64_t p = (int32_t)pos[r];
            const int64_t np = nextPos[r];
            if (p >= 0 && np >= 0)
                mateMap[{p, np}] = r;
        }
        for (uint32_t r = 0; r < nRec; ++r) {
            auto it = mateMap.find({nextPos[r], (int32_t)pos[r]});
            const uint32_t m = (it != mateMap.end()) ? it->second : nRec;
            const uint32_t ms = (m < nRec) ? span[m] : 0;
            tlen[r] = bamTlenCompute(flag[r], refId[r], (int64_t)(int32_t)pos[r],
                                     nextRefId[r], nextPos[r],
                                     span[r], ms, tlenConvMinusOne);
        }
        for (size_t i = 0; i < tlenExcLine.size(); ++i)
            tlen[tlenExcLine[i]] = tlenExcVal[i];
    }

    /* Render SAM text lines, byte-identical to the textual BAM->SAM path.
       With bamOut the line is only an intermediate: bamScratch turns it into
       the BAM bytes appended to the block, and the writer thread just frames
       them into BGZF (BamWriter::writeBlock). */
    bamrec::BamRecordScratch bamScratch;
    size_t tagOff = 0;
    size_t qualReadOff = 0;
    for (uint32_t r = 0; r < nRec; ++r) {
        std::string line;
        line.reserve(64 + 2 * lSeq[r]);
        line += qname[r];
        line += "\t";
        line += std::to_string(flag[r]);
        line += "\t";
        if (refId[r] < 0 || (size_t)refId[r] >= chrNames.size()) {
            line += "*";
        } else {
            line += chrNames[refId[r]];
        }
        line += "\t";
        line += (pos[r] == 0xFFFFFFFFu) ? "0" : std::to_string(pos[r] + 1);
        line += "\t";
        line += std::to_string(mapq[r]);
        line += "\t";
        if (cigar[r].empty()) {
            line += "*";
        } else {
            for (size_t k = 0; k < cigar[r].size(); ++k) {
                cigarOpText(cigar[r][k], line);
            }
        }
        line += "\t";
        if (nextRefId[r] < 0 || (size_t)nextRefId[r] >= chrNames.size()) {
            line += "*";
        } else if (nextRefId[r] == refId[r]) {
            line += "=";
        } else {
            line += chrNames[nextRefId[r]];
        }
        line += "\t";
        line += (nextPos[r] < 0) ? "0" : std::to_string((uint32_t)nextPos[r] + 1);
        line += "\t";
        line += std::to_string(tlen[r]);
        line += "\t";
        if (lSeq[r] > 0) {
            line.append((const char*)seqAll.data() + seqOff[r], lSeq[r]);
        } else {
            line += "*";
        }
        line += "\t";
        bool qualMissing = false;
        if (lSeq[r] > 0) {
            const uint8_t* qseg = qualAll.data() + qualReadOff;
            bool allStar = true;
            for (uint32_t qi = 0; qi < lSeq[r]; ++qi) {
                if (qseg[qi] != '*') {
                    allStar = false;
                    break;
                }
            }
            qualMissing = allStar;
            qualReadOff += lSeq[r];
        }
        if (lSeq[r] == 0 || qualMissing) {
            line += "*";
        } else {
            line.append((const char*)qualAll.data() + qualReadOff - lSeq[r], lSeq[r]);
        }
        if (tagsLen[r] > 0) {
            renderBamTags(tagAll.data() + tagOff, tagAll.data() + tagOff + tagsLen[r], line);
            tagOff += tagsLen[r];
        }
        line += "\n";

        const uint8_t* payload = (const uint8_t*)line.data();
        size_t payloadLen = line.size();
        if (bamOut) {
            const int32_t built = bamrec::buildBamRecordFromSamLine(
                (const uint8_t*)line.data(), line.size(), refIndex, bamScratch);
            if (built < 0) {
                LOG_ERROR("BamCodecActuator::decompress: BAM record build failed");
                return -1;
            }
            if (built == 0) {
                continue;   /* invalid record: the text path drops it in BamWriter too */
            }
            payload = bamScratch.rec.data();
            payloadLen = bamScratch.rec.size();
        }

        if ((size_t)outBlockPtr->getDataLen() + payloadLen > outBlockPtr->getBufferSize()) {
            if (0 != outBlockPtr->ensureCapacity((size_t)outBlockPtr->getDataLen() + payloadLen)) {
                LOG_ERROR("BamCodecActuator::decompress: output buffer grow failed");
                return -1;
            }
        }
        memcpy(outBlockPtr->getBuffer() + outBlockPtr->getDataLen(), payload, payloadLen);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + (int64_t)payloadLen);
    }
    if (bamOut) {
        outBlockPtr->setBamPrebuilt(true);
    }
    return 0;
}
