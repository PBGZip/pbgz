/*
 * bam_actuator.h - column-wise BAM compression without SAM text
 *
 * BamCodecActuator consumes the BamColumns produced by BamBlockReader
 * (see bam_columns.h) instead of SAM text lines. Each column is encoded by the
 * very same coder that preprocessing selected for the corresponding SAM field
 * (makeFieldEncoder), so no field changes its encoder - only the representation
 * of the bytes fed to it:
 *
 *   - scalars stay scalars (FLAG u16, POS i32, MAPQ u8, RNEXT/PNEXT/TLEN i32),
 *     delta/varint transformed before entropy coding;
 *   - bases stay 4-bit packed (half the bytes of the textual SEQ column);
 *   - QUAL keeps the raw bytes and is still fed record by record, which is what
 *     fcv2 needs (per-record length + strand).
 */

#pragma once

#include <json/json.h>

#include "bam_columns.h"
#include "codec_actuator.h"
#include "reference.h"

class BamCodecActuator : public CodecActuator {
public:
    BamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr,
                     Reference* pRef = nullptr);
    virtual ~BamCodecActuator() override {}

    int32_t compress();
    int32_t decompress();

private:
    /*
     * Encode one already-serialised column with the coder preprocessing picked
     * for fieldIdx, and append its meta to streamMeta. Used for every column
     * whose payload can be handed to the coder as a single stream.
     */
    int32_t encodeColumn(uint32_t fieldIdx, CoderType fallback, const uint8_t* data, size_t size,
                         Json::Value& streamMeta, const char* name, const char* mode);

    /* QNAME is a column-wise prefix/suffix-matching coder: it must be fed line
       by line or it degrades to an ordinary context model, which loses the
       shared-prefix compression of the textual path. */
    int32_t encodeQnameColumn(const uint8_t* qnameStream, size_t qnameStreamSize,
                              uint32_t nRecords, Json::Value& streamMeta);

    /* Whole-column coding with a plain bwt_cm block, used for the columns whose
       preselected coder (coder_fc / affix) has a fragile whole-block decoder on
       the structured side. Deterministic and symmetric to decode. */
    int32_t encodeRawBwtColumn(uint32_t fieldIdx, const uint8_t* data, size_t size,
                               Json::Value& streamMeta, const char* name);

    /* LEB128 helpers for the numeric columns. */
    static void putUvarint(std::vector<uint8_t>& out, uint64_t v);
    static void putSvarint(std::vector<uint8_t>& out, int64_t v);

    /* Columns that need per-record handling or a transform. */
    int32_t encodePosColumn(Json::Value& streamMeta);
    int32_t encodePnextColumn(Json::Value& streamMeta);
    int32_t encodeTlenColumn(Json::Value& streamMeta);
    int32_t encodeCigarColumn(Json::Value& streamMeta);
    int32_t encodeQualColumn(Json::Value& streamMeta);

    /* SEQ against a reference: per-base diff stream (0 = matches reference,
       otherwise the read base character). Falls back to plain text SEQ when no
       reference is usable. */
    int32_t encodeSeqColumnRef(Json::Value& streamMeta);

    uint32_t pickedFor(uint32_t fieldIdx, CoderType fallback) const;

    BamColumns* cols;
    Reference* pRefeGene;
};
