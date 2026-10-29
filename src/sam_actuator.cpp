/*
 * sam_actuator.cpp - Source file for pbgz project
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

#include <memory>
#include <sstream>
#include "sam_actuator.h"
#include "profile_stats.h"
#include "sam_field_layout.h"
#include "sam_field_rules.h"
#include "bam_actuator.h"
#include "block_wrapper.h"
#include "coder/coder_io.h"
#include "coder/coder_fc.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_arith.h"
#include "utils/checksum.h"
#include "city.h"
#include "coder/coder_json.h"

#include <cstdlib>
#include <cstring>
#include "log/logger.h"
#include "sam_info.h"
#include "utils/path_util.h"
#include "pbgz_index.h"
#include "compress_engine.h"
#include "pbgz_stat.h"
#include "config_manager.h"
#include "field_coder_config.h"

namespace {
    uint16_t mapFieldIdxToStatUnitId(uint16_t fieldIdx) {
        switch (fieldIdx) {
            case 0: return StatObjectId::SAM_QNAME;
            case 1: return StatObjectId::SAM_FLAG;
            case 2: return StatObjectId::SAM_RNAME;
            case 3: return StatObjectId::SAM_POS;
            case 4: return StatObjectId::SAM_MAPQ;
            case 5: return StatObjectId::SAM_CIGAR;
            case 6: return StatObjectId::SAM_RNEXT;
            case 7: return StatObjectId::SAM_PNEXT;
            case 8: return StatObjectId::SAM_TLEN;
            case 9: return StatObjectId::SAM_SEQ;
            case 10: return StatObjectId::SAM_QUAL;
            case 11: return StatObjectId::SAM_OPTION;
            default: return 0;
        }
    }

    void recordFieldStats(PbgzEngine* engine, uint32_t fieldIdx, uint32_t fieldSrcLen, uint32_t fieldDstLen) {
        if (!engine) return;

        auto compressEngine = dynamic_cast<CompressEngine*>(engine);
        if (!compressEngine || !compressEngine->getStats()) return;

        auto samStat = dynamic_cast<SamStat*>(compressEngine->getStats());
        if (!samStat) return;

        uint16_t statObjectId = mapFieldIdxToStatUnitId(fieldIdx);
        if (statObjectId != 0) {
            samStat->addMetricValue(StatUnitIds::COMPRESSION_RATIO, statObjectId, StatMetricIds::ORIGINAL_SIZE, fieldSrcLen);
            samStat->addMetricValue(StatUnitIds::COMPRESSION_RATIO, statObjectId, StatMetricIds::COMPRESSED_SIZE, fieldDstLen);
        }
    }

    /* Basename used to compare an alignment's reference with the -r file. Both sides may be
       gzip-compressed and one may be a path in an @PG command line. Only the actual -r file is
       inspected for the gzip header's original name; a path parsed from the header is not opened. */
    std::string referenceBaseName(const std::string& path, bool inspectGzipHeader = false) {
        if (path.empty()) {
            return path;
        }
        std::string name = (inspectGzipHeader && PathUtil::isGzFile(path))
                               ? PathUtil::getFileNameFromGz(path)
                               : PathUtil::getFileName(path);
        if (name.size() > 3 && name.compare(name.size() - 3, 3, ".gz") == 0) {
            name.resize(name.size() - 3);
        }
        return name;
    }

    /* bwa mem options that consume the following token when the value is not attached
       (e.g. "-t 8", "-B 9", "-R @RG..."). Flags such as -M/-V/-Y/-q take no value. */
    bool bwaMemOptionTakesValue(char opt) {
        switch (opt) {
            case 't': case 'k': case 'w': case 'd': case 'r': case 'y': case 'c':
            case 'D': case 'W': case 'm': case 'A': case 'B': case 'O': case 'E':
            case 'L': case 'U': case 'x': case 'R': case 'H': case 'o': case 'j':
            case 'K': case 'v': case 'T': case 'h': case 'I':
                return true;
            default:
                return false;
        }
    }

    /* First positional argument of a "bwa mem ..." command line: the index/reference prefix.
       Option values are skipped so a value such as the "8" in "-t 8" is not mistaken for it. */
    std::string extractBwaMemReference(const std::string& command) {
        std::stringstream ss(command);
        std::string item;
        std::vector<std::string> tokens;
        while (std::getline(ss, item, ' ')) {
            if (!item.empty()) {
                tokens.push_back(item);
            }
        }
        if (tokens.size() < 3 || tokens[0] != "bwa" || tokens[1] != "mem") {
            return std::string();
        }
        for (size_t i = 2; i < tokens.size(); ++i) {
            const std::string& token = tokens[i];
            if (token.empty()) {
                continue;
            }
            if (token[0] != '-') {
                return token;
            }
            if (token == "--") {
                continue;
            }
            if (token.size() > 2 && token[1] == '-') {
                /* Long option: only the "--name=value" form is self-contained here. */
                continue;
            }
            if (token.size() <= 1) {
                continue;
            }
            if (bwaMemOptionTakesValue(token[1]) && token.size() == 2) {
                ++i;   /* skip the option's value */
            }
        }
        return std::string();
    }
}

SamCodecActuator::SamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine, Reference* pReferene): CodecActuator(inPtr, outPtr, engine) {
    pRefeGene = pReferene;
    headEndLine = 0;
    idAnalysis.posLength = 0;
    headerSrcLen = 0;
    headerDstLen = 0;
    readOffset = 0;
    baseNCount = 0;
    qualDecoder = nullptr;
    baseLengthBuffer = nullptr;
    baseSquashBuffer = nullptr;
    baseDiffSquashBuffer = nullptr;
    refeStrecchBuffer = nullptr;
    samLine = 0;
    refPosChrIndex = 65535;
    refPosBegin = 0;
    refPosEnd = 0;
    optionCacheEmpty = true;
    optionRecText.clear();
    optionRecOff.clear();
    optionRecPresent.clear();
    optionRecIds.clear();
    optionTagNames.clear();
    optionTagTypes.clear();
    optionTagCols.clear();
}

SamCodecActuator::~SamCodecActuator() {
    MemoryUtil::safeFree(baseLengthBuffer);
    MemoryUtil::safeFree(baseSquashBuffer);
    MemoryUtil::safeFree(baseDiffSquashBuffer);
    MemoryUtil::safeFree(refeStrecchBuffer);
    MemoryUtil::safeFree(matchBlockBuffer);

    // Release idDecoders
    idDecoders.clear();
    clearIdNumericState();

    // Release fieldDecoders
    fieldDecoders.clear();

    // Release the QUAL decoder
    qualDecoder.reset();

    ioVector.clear();
}

int32_t SamCodecActuator::preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLen) {
    return analyzeQnameFirstLine(pBuffer, bufLen, idAnalysis);
}

int32_t SamCodecActuator::preAnalysisIdLine(uint8_t* pBuffer, uint32_t bufferLen) {
    return analyzeQnameLine(pBuffer, bufferLen, idAnalysis);
}

/* The QNAME column's lines as the layout encoders want them (see sam_qname_column.h). */
QnameColumnInput SamCodecActuator::qnameColumnInput() const
{
    QnameColumnInput in;
    in.buffer = inBlockPtr->getBuffer();
    in.npos = &inBlockPtr->getNpos();
    in.headEndLine = headEndLine;
    in.fieldTabs = &contentPos;
    return in;
}

int32_t SamCodecActuator::compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    /* 0 = the whole block (the encoders' line bound is the trial's). */
    return encodeQnameSplit(qnameColumnInput(), idAnalysis, outBlockPtr, &ioErrSink, 0, fieldMeta,
                            fieldSrcLen);
}

int32_t SamCodecActuator::compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    return encodeQnameByQname(qnameColumnInput(), outBlockPtr, &ioErrSink, 0, fieldMeta,
                              fieldSrcLen);
}


/* Analyse one SAM header line: capture @SQ chromosomes and check the reference/@PG names. Returns 0, or -1 on error. */
int32_t SamCodecActuator::parseHeaderLine(const std::string& line)
{
    if (line.length() < 3) {
        return -1;
    }
    headEndLine++;
    if (line.substr(0, 3) == "@HD" || line.substr(0, 3) == "@RG" || line.substr(0, 3) == "@CO") {
        return 0;
    } else if (line.substr(0, 3) == "@SQ") {
        // Parse chromosome information from @SQ line
        if (SamUtil::parseChromosomeInfo(line) != 0) {
            LOG_ERROR("Failed to parse chromosome info from line: %s", line.c_str());
            return -1;
        }

        // Extract reference file path from @SQ line
        // Format: @SQ SN:ref_name LN:length UR:file_path
        static bool isCheckUR = false;
        if (!isCheckUR) {
            size_t pos = line.find("UR:");
            if (pos != std::string::npos) {
                pos += 3;
                size_t endPos = line.find("\t", pos);
                if (endPos == std::string::npos) {
                    endPos = line.length();
                }
                std::string refFilePath = line.substr(pos, endPos - pos);
                std::string refFileName = referenceBaseName(refFilePath);

                LOG_INFO("Reference fasta name: %s.", refFilePath.c_str());
                if (pRefeGene != nullptr) {
                    std::string inputFastq = referenceBaseName(pRefeGene->getFastaFileName(), true);
                    if (refFileName != inputFastq) {
                        fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fasta is %s. \n", refFileName.c_str(), inputFastq.c_str());
                        isCheckUR = true;
                    }
                }
            }
        }
        return 0;
    } else if (line.substr(0, 3) == "@PG") {
        static bool isCheckPG = false;
        if (!isCheckPG) {
            size_t pos = line.find("CL:");
            if (pos != std::string::npos) {
                pos += 3;
                size_t endPos = line.find("\t", pos);
                if (endPos == std::string::npos) {
                    endPos = line.length();
                }

                std::string command = line.substr(pos, endPos - pos);
                std::string refGeneName = extractBwaMemReference(command);
                if (!refGeneName.empty()) {
                    LOG_INFO("Reference gene name extracted from @PG CL: %s", refGeneName.c_str());
                    if (pRefeGene != nullptr) {
                        std::string inputFastq = referenceBaseName(pRefeGene->getFastaFileName(), true);
                        const std::string samFasta = referenceBaseName(refGeneName);
                        if (samFasta != inputFastq) {
                            fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fasta is %s \n", samFasta.c_str(), inputFastq.c_str());
                            isCheckPG = true;
                        }
                    }
                }
            }
        }
    } else {
        LOG_ERROR("Unexpected header %s", line.c_str());
        return -1;
    }
    return 0;
}

/* Analyse one record line: ID to the ID analysis, RNAME check, CIGAR/SEQ/QUAL length checks, base count and quality histogram. */
int32_t SamCodecActuator::scanDataLine(const uint8_t* line, uint32_t lineLen, uint32_t idx,
                                       std::pair<uint8_t, uint32_t>* qualityFrequnce)
{
    /* The line's tab positions (see LineTabs). */
    LineTabs linePos;
    uint32_t baseLen = 0;
    bool lineCigarMatchFlag = false;
    uint32_t baseFieldLen = 0;
    /* Walk the line through a raw pointer and hoisted length; this runs once per line of every block. */
    const char* const src = (const char*)line;
    const uint32_t srcLen = lineLen;
    for (uint32_t i = 0; i < srcLen; ++i) {
        if (src[i] == '\t' || src[i] == '\n') {
            // First tab before is ID column, need to split and analyze ID column
            if (linePos.empty()) {
                if (contentPos.empty()) {
                    // Pass ID column content (from line start to first tab position)
                    preAnalysisIdFirstLine((uint8_t*)src, i + 1);
                } else {
                    // Pass ID column content (from line start to first tab position)
                    preAnalysisIdLine((uint8_t*)src, i + 1);
                }
            } else if (linePos.size() == 2) {  // RNAME is the 3th field
                uint32_t chrLen = i - (uint32_t)linePos[1] - 1;
                const char* chrBegin = src + (size_t)linePos[1] + 1;
                const bool isStar = (chrLen == 1 && *chrBegin == '*');
                const bool isEq = (chrLen == 1 && *chrBegin == '=');
                if (!isStar && !isEq) {
                    /* Mapped read: build the (rare) name only to look it up. */
                    std::string chrName(chrBegin, chrLen);
                    if (SamInfo::getInstance().getChrNameIndex(chrName) == 65535) {
                        return -1;
                    }
                }
            } else if (linePos.size() == 5) {  // CIGAR is the 6th field
                uint32_t cigarLen = i - (uint32_t)linePos[4] - 1;
                if (cigarLen > 1) {
                    lineCigarMatchFlag = true;
                    uint8_t* cigarBegin  = (uint8_t*)src + (size_t)linePos[4] + 1;
                    baseFieldLen = cigarSeqConsumed(cigarBegin, cigarLen);
                } else {
                    lineCigarMatchFlag = false;
                }
            } else if (linePos.size() == 9) {   // Base is the 10th field
                baseLen = i - (uint32_t)linePos[8] - 1;
                if (baseLen > maxBaseLength) {
                    maxBaseLength = baseLen;
                }
                if (!lineCigarMatchFlag) {
                    if (baseLen < minBaseLength) {
                        minBaseLength =  baseLen;
                    }
                } else {
                    if (baseFieldLen != baseLen) {
                        return -1;
                    }
                }
            } else if (linePos.size() == 10) {  // Quality value is the 11th field
                uint32_t qualityLen = i - (uint32_t)linePos[9] - 1;
                if (baseLen != qualityLen) {
                    /* A QUAL of a single '*' denotes missing quality; its length need not equal the SEQ length */
                    const bool missingQual = (qualityLen == 1 &&
                                              src[(size_t)linePos[9] + 1] == '*');
                    if (!missingQual) {
                        LOG_ERROR("Not a valid sam data, baselen = %u, qualityLen = %u ", baseLen, qualityLen);
                        return -1;
                    }
                }
            }
            linePos.push_back(i);
            // Optional fields compressed as whole block
            if (linePos.size() == 11) {
                break;
            }
        } else {
            if (linePos.size() == 9) {
                char ch = src[i];
                if (ch == 'N' || ch == 'n') {
                    baseNCount++;
                }
            } else if (linePos.size() == 10) {
                /* Index by uint8_t: quality bytes are histogrammed in a 256-entry
                   table and char signedness differs between targets. */
                const uint8_t ch = line[i];
                qualityFrequnce[ch].second++;
            }
        }
    }
    // Each line must have at least 10 tabs, less than 10 means not a SAM file
    if (linePos.size() < 10) {
        LOG_ERROR("Sam field line invalid, size = %d", linePos.size());
        return -1;
    }
    if (linePos.size() > maxFieldSize) {
        maxFieldSize = linePos.size();
    }
    lineFiledCount.push_back(std::make_pair(idx, linePos.size()));
    contentPos.push_back(linePos);
    return 0;
}

int32_t SamCodecActuator::preAnalysis() {
    if (inBlockPtr == nullptr) {
        LOG_ERROR("Input block pointer is null");
        return -1;
    }

    /* Archive-mode reader offload: render a raw-record BAM block into SAM text before the usual passes. */
    if (inBlockPtr->isRawBam()) {
        if (0 != bamRenderRawBlockToText(inBlockPtr)) {
            LOG_ERROR("SAM block(%ld): failed to render raw BAM records to text.",
                      (long)inBlockPtr->getBlockId());
            return -1;
        }
    }

    /* A non-empty block with no line position is not SAM/BAM text; fail so the caller can store it verbatim as binary. */
    if (inBlockPtr->getDataLen() > 0 && inBlockPtr->getNpos().empty()) {
        LOG_ERROR("SAM block(%ld) carries %ld bytes but no line positions, not SAM/BAM text.",
                  (long)inBlockPtr->getBlockId(), (long)inBlockPtr->getDataLen());
        return -1;
    }

    /* SEQ's file-level decisions (match coder, N form) come from the codec pre-selection; this pass only fills column stats. */
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Initialize quality frequency analysis similar to FastqActuator
    std::pair<uint8_t, uint32_t> qualityFrequnce[256];
    for (int32_t i = 0; i < 256; ++i) {
        qualityFrequnce[i].first = i;
        qualityFrequnce[i].second = 0;
    }

    /* One reusable line buffer for the whole block (assign reuses the previous line's storage). */
    std::string line;
    // Parse SAM file header
    for (uint32_t idx = 0; idx < lineNum; ++idx) {
        uint32_t begin = (idx == 0 ? idx : npos[idx - 1]  + 1);
        uint32_t end = (idx < npos.size() ? npos[idx] : strlen((char*)buffer + begin));
        if (begin >= end) {
            return -1;
        }
        if (buffer[begin] == '@') {
            /* Header lines are few (split into their own block), so a copy here is fine. */
            line.assign((char*)buffer + begin, end - begin);
            if (parseHeaderLine(line) != 0) {
                return -1;
            }
        } else {
            /* Scan the record straight out of the block buffer; no per-line copy. */
            if (scanDataLine(buffer + begin, end - begin, idx, qualityFrequnce) != 0) {
                return -1;
            }
        }
    }

    // All match scenario
    if (minBaseLength == UINT32_MAX) {
        minBaseLength = maxBaseLength;
    }

    LOG_DEBUG("minBaseLength = %u, maxBaseLength = %u", minBaseLength, maxBaseLength);

    inBlockPtr->setMaxLineLen(maxBaseLength);

    std::sort(qualityFrequnce, qualityFrequnce + 256,
        [](const std::pair<uint8_t, uint32_t> &a, const std::pair<uint8_t, uint32_t> &b) { return  a.second > b.second;});
    for (int i = 0; i < 256; i++) {
        if (qualityFrequnce[i].second == 0) {
            continue;
        }
        qualFreqTable.push_back(std::make_pair(qualityFrequnce[i].first - '!', 1));
    }

    // Compression strategy judgment logic - refer to FastqActuator implementation
    if (idAnalysis.posLength != UINT32_MAX) {
        for (uint32_t idx = 0; idx < idAnalysis.symbols.size(); ++idx) {
            if (idAnalysis.minLen[idx] == 0) {
                idAnalysis.posLength = UINT32_MAX;
                break;
            }
        }
    }

    return 0;
}

int32_t SamCodecActuator::compress() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter, inBlockPtr or outBlockPtr is nullptr for SAM compression");
        return -1;
    }

    /* Same raw-record render as preAnalysis, for callers that compress without it. */
    if (inBlockPtr->isRawBam()) {
        if (0 != bamRenderRawBlockToText(inBlockPtr)) {
            LOG_ERROR("SAM block(%ld): failed to render raw BAM records to text.",
                      (long)inBlockPtr->getBlockId());
            return -1;
        }
    }

    /* Grow the output block to the input block's real size: -l is only an allocation hint, never a cap. */
    if (pbgzEngine != nullptr) {
        size_t bs = pbgzEngine->getFileBlockSize();
        if (bs == 0) {
            bs = ConfigManager::getInstance().getBlockSizeByCompressLevel(pbgzEngine->getParameter().compressLevel);
        }
        size_t outCapacity = bs * 2;
        const size_t inDataLen = (size_t)inBlockPtr->getDataLen();
        if (inDataLen > outCapacity) {
            outCapacity = inDataLen;
        }
        if (outBlockPtr->getBufferSize() < outCapacity) {
            if (0 != outBlockPtr->ensureCapacity(outCapacity)) {
                LOG_ERROR("Preallocate output buffer failed, need=%zu bytes", outCapacity);
                return -1;
            }
        }
    }

    if (headEndLine > 0) {
        PBGZ_PROF_SCOPE(pbgzprof::FIELD_HDR);
        if (0 != compressSamHeader()) {
            LOG_ERROR("Compress SAM header failed.");
            return -1;
        }
    }

    /* Same guard as preAnalysis: a block with bytes but no line position must not be written header-only. */
    if (inBlockPtr->getDataLen() > 0 && inBlockPtr->getNpos().empty()) {
        LOG_ERROR("SAM block(%ld) carries %ld bytes but no line positions, cannot code as text.",
                  (long)inBlockPtr->getBlockId(), (long)inBlockPtr->getDataLen());
        return -1;
    }

    if (inBlockPtr->getNpos().size() <= (size_t)headEndLine) {
        LOG_DEBUG("SAM head only");
        // For header-only case, we still need to write the metadata
        // Calculate checksum of original data (header only)
        std::string blockHash;
        if (inBlockPtr->hasMd5()) {
            /* Already hashed by the reader thread (CompressEngine::preDispatchBlock). */
            blockHash = inBlockPtr->getMd5();
        } else {
            blockHash = checksumHex(ChecksumType::FarmHash128,
                                    inBlockPtr->getBuffer(), (size_t)inBlockPtr->getDataLen());
        }
        setChecksumMeta(meta, ChecksumType::FarmHash128, blockHash);

        // Compress and write metadata
        coder_json metaCoder;
        int32_t metaLen = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
        if (metaLen <= 0) {
            LOG_ERROR("Failed to encode meta information for SAM compression (header only)");
            return -1;
        }
        outBlockPtr->setMetaLen(metaLen);

        // Set block information
        outBlockPtr->setBlockId(inBlockPtr->getBlockId());
        outBlockPtr->setBlockType(inBlockPtr->getBlockType());

        return 0;
    }

    // Check if preAnalysis was successful
    if (contentPos.empty()) {
        LOG_ERROR("preAnalysis() must be called before compress()");
        return -1;
    }

    if (0 != compressSamByFields()) {
        LOG_ERROR("Compress Sam Fields failed.");
        return -1;
    }

    if (buildSamIndex() != 0) {
        PBGZ_PROF_SCOPE(pbgzprof::FIELD_INDEX);
        LOG_ERROR("Build Sam index failed.");
    }

    return 0;
}

int32_t SamCodecActuator::compressSamHeader() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter for SAM header compression");
        return -1;
    }

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Clear previous chromosome information but keep chromosome ID counter for multi-block compatibility
    // Note: Do not clear chromosome information as multiple blocks may need to process chromosome info
    headerSrcLen = 0;
    // Create SAM file header compressor
    std::shared_ptr<coder_io> headerIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SAM header");
    std::shared_ptr<coder> headerCoder =
        CoderFactory::makeProfileEncoder(SAM_QNAME, compressMode(), headerIo.get(),
                                         CoderType::BWT_CM, true, engineCompressLevel());

    // Process SAM file header line by line
    for (uint32_t lineId = 0; lineId < headEndLine; ++lineId) {
        uint32_t begin = (lineId == 0 ? 0 : npos[lineId - 1] + 1);
        uint32_t end = npos[lineId];
        if (begin >= end) {
            continue;
        }

        uint32_t lineLength = end - begin + 1;
        // Compress current line (including newline character)
        headerCoder->encode_line(buffer + begin, lineLength);
        headerSrcLen += lineLength;
    }

    // Complete compression
    headerCoder->encode_flush();
    if (headerIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode SAM header overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + headerIo->data_len);
    headerDstLen = headerIo->data_len;

    // Set SAM file header metadata
    Json::Value headerMeta;
    headerMeta["srclen"] = headerSrcLen;
    headerMeta["dstlen"] = headerDstLen;
    headerMeta["lines"] = headEndLine;
    headerMeta["coder"] = headerIo->meta;
    meta["header"] = headerMeta;

    LOG_INFO("SAM header compression completed: %u lines, %u bytes -> %u bytes, compress ratio = %.2f%%",
             headEndLine, headerSrcLen, headerDstLen, (double)(headerDstLen* 100)/(double)headerSrcLen);

    return 0;
}

/* The QNAME column (field 0), written with the layout the codec pre-selection picked (affix or coder_qname). Negative on failure. */
int32_t SamCodecActuator::encodeQnameColumn(uint32_t& fieldSrcLen, Json::Value& fieldMeta)
{
    // ID field: compress based on analysis result
    if (idAnalysis.posLength == UINT32_MAX) {
        LOG_DEBUG("Id will compress in all.");
        return compressIdFieldInAll(fieldSrcLen, fieldMeta);
    }

    /* Fast mode writes QNAME as one column through coder_rans; the affix/qname layouts are
       archive-profile coders and would be rejected at encoder construction anyway. */
    if (compressMode() == PBGZ_MODE_FAST) {
        LOG_DEBUG("Id will compress in all (fast mode).");
        return compressIdFieldInAll(fieldSrcLen, fieldMeta);
    }

    /* Which layout writes this column is the file-level verdict (PreprocessInfo::qnameUseAffix); affix is the fallback. */
    const PreprocessInfo* preInfo = preprocessInfoMut();
    const int32_t verdict =
        (preInfo != nullptr) ? preInfo->qnameUseAffix.load(std::memory_order_relaxed) : -1;
    const bool useAffix = (verdict != 0);
    LOG_DEBUG("QNAME: file verdict %d -> %s", verdict, useAffix ? "affix" : "qname");
    return useAffix ? compressIdFieldSplit(fieldSrcLen, fieldMeta)
                    : compressIdFieldQname(fieldSrcLen, fieldMeta);
}

/* Close the block: write the SAM meta and per-field stream metas, hash the block, encode the meta. */
int32_t SamCodecActuator::finalizeSamBlockMeta(Json::Value& samMeta,
                                               const Json::Value& streamMeta, uint32_t lineNumber,
                                               uint32_t fieldCount, uint32_t totalSrcLen,
                                               uint32_t totalDstLen)
{
    // Set SAM metadata
    samMeta["lines"] = lineNumber;
    samMeta["fieldcount"] = fieldCount;
    samMeta["totalsrclen"] = totalSrcLen;
    samMeta["totaldstlen"] = totalDstLen;
    /* Exact decoded-text size of this block, used to size the decompressor's output buffer. */
    samMeta["textlen"] = (Json::Value::UInt64)inBlockPtr->getDataLen();
    samMeta["streams"] = streamMeta;
    meta["sam"] = samMeta;

    PBGZ_PROF_SCOPE(pbgzprof::FIELD_META);

    // Calculate checksum of original data
    std::string blockHash;
    if (inBlockPtr->hasMd5()) {
        /* Already hashed by the reader thread (CompressEngine::preDispatchBlock). */
        blockHash = inBlockPtr->getMd5();
    } else {
        blockHash = checksumHex(ChecksumType::FarmHash128,
                                inBlockPtr->getBuffer(), (size_t)inBlockPtr->getDataLen());
    }
    setChecksumMeta(meta, ChecksumType::FarmHash128, blockHash);

    // Compress metadata
    coder_json metaCoder;
    int32_t metaLen = metaCoder.encoder(meta, outBlockPtr->getMetaBuffer(), outBlockPtr->getRemain());
    if (metaLen <= 0) {
        LOG_ERROR("Failed to encode meta information for SAM compression");
        return -1;
    }
    outBlockPtr->setMetaLen(metaLen);

    // Set block information
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());
    return 0;
}

int32_t SamCodecActuator::compressSamByFields() {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size() - headEndLine;

    // Initialize metadata
    Json::Value samMeta;
    Json::Value streamMeta;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;

    // Compress each field (column) of SAM format separately
    // SAM format has at least 11 required fields, we'll compress each column
    uint32_t fieldCount = maxFieldSize + 1; // +1 for the last field after last tab

    // Compress each field separately
    for (uint32_t fieldIdx = 0; fieldIdx < fieldCount; ++fieldIdx) {
        uint32_t fieldSrcLen = 0;
        Json::Value fieldMeta;
        uint32_t fieldDstLen = 0;
        const std::chrono::steady_clock::time_point profFieldT0 =
            pbgzprof::enabled() ? std::chrono::steady_clock::now()
                                : std::chrono::steady_clock::time_point();
        switch (fieldIdx) {
            case 0: // QNAME as FQ:ID
                fieldDstLen = encodeQnameColumn(fieldSrcLen, fieldMeta);
                break;
            case 1: // FLAG
                if (pickedCoderFor(fieldIdx, CoderType::BWT_CM) == CoderType::AFFIX_MATCH) {
                    fieldDstLen = compressRegularField(fieldIdx, fieldSrcLen, fieldMeta);
                } else {
                    fieldDstLen = compressNumber<uint16_t>(fieldIdx, fieldSrcLen, fieldMeta);
                }
                break;
            case 2: // RNAME
                fieldDstLen = compressChrName(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 3: // POS
                /* POS as unsigned varint (LEB128) deltas against the previous line, baseline reset at each chromosome switch. */
                fieldDstLen = compressPosFieldDelta(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 4: // MAPQ
                if (pickedCoderFor(fieldIdx, CoderType::BWT_CM) == CoderType::AFFIX_MATCH) {
                    fieldDstLen = compressRegularField(fieldIdx, fieldSrcLen, fieldMeta);
                } else {
                    fieldDstLen = compressNumber<uint8_t>(fieldIdx, fieldSrcLen, fieldMeta);
                }
                break;
            case 5: // CIGAR
                fieldDstLen = compressCigar(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 6: // RNEXT
                fieldDstLen = compressChrName(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 7: // PNEXT
                /* PNEXT as the delta against POS; always the delta path, excluded from coder selection. */
                fieldDstLen = compressPNextFieldDelta(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 8: // TLEN
                /* TLEN reconstructed from POS/PNEXT/CIGAR on decompression, with exceptions only for lines that cannot be rebuilt. */
                fieldDstLen = compressTLen(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 9: // SEQ
                /* A reference the @SQ list cannot address is dropped for SEQ; the verdict is taken after the header is parsed. */
                if (!referenceUsableForSeq()) {
                    /* No usable reference: code against a block-built consensus where possible, else store the bases (see compressBaseSelfRef). */
                    fieldDstLen = compressBaseSelfRef(fieldIdx, fieldSrcLen, fieldMeta);
                } else {
                    LOG_DEBUG("Base will compress with reference");
                    fieldDstLen = compressBaseWithRef(fieldIdx, fieldSrcLen, fieldMeta);
                }
                break;
            case 10: // QUAL
                fieldDstLen = compressQuality(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 11: // Optional fields (all tags)
            {
                /* No trial: the per-tag layout is taken when it applies (MD:Z/NM:i regenerable,
                   where it wins by ~30%); otherwise the stripped layout is used, and the plain
                   whole column only if that fails too. */
                const int64_t savedDataLen = outBlockPtr->getDataLen();
                const int32_t splitDst = compressOptionField(fieldSrcLen, fieldMeta);
                if (splitDst >= 0) {
                    fieldDstLen = (uint32_t)splitDst;
                } else {
                    outBlockPtr->setDataLen(savedDataLen);
                    uint32_t stripDst = 0;
                    if (compressOptionStrip(fieldIdx, fieldSrcLen, fieldMeta, stripDst)) {
                        fieldDstLen = stripDst;
                    } else {
                        outBlockPtr->setDataLen(savedDataLen);
                        fieldDstLen = compressRegularField(fieldIdx, fieldSrcLen, fieldMeta);
                    }
                }
                break;
            }
        }

        // Record statistics for this field
        if (pbgzprof::enabled()) {
            const uint64_t profUs = (uint64_t)std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - profFieldT0).count();
            pbgzprof::add(pbgzprof::FIELD_BASE + (int)fieldIdx, profUs);
        }
        
        recordFieldStats(pbgzEngine, fieldIdx, fieldSrcLen, fieldDstLen);

        streamMeta.append(fieldMeta);
        totalSrcLen += fieldSrcLen;
        totalDstLen += fieldDstLen;
    }

    if (finalizeSamBlockMeta(samMeta, streamMeta, lineNum, fieldCount, totalSrcLen, totalDstLen) != 0) {
        return -1;
    }

    // Set block information
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());

    return 0;
}

CoderType SamCodecActuator::pickedCoderFor(uint32_t fieldIdx, CoderType fallback) const
{
    const PreprocessInfo* preInfo =
        (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
    if (preInfo != nullptr) {
        return preInfo->coderFor(fieldIdx, fallback);
    }
    return fallback;
}

/* Where one field of one content line sits: line start, tab positions, 1-based field index, and the two lengths. */
struct SamFieldRange {
    uint8_t* start = nullptr;
    /* The delimiter that ends the field: the next tab, or the end of the line for the last. */
    uint8_t* tab = nullptr;

    uint32_t lengthWithoutTab() const { return (uint32_t)(tab - start); }
    uint32_t lengthWithTab() const { return (uint32_t)(tab - start) + 1; }
};

static SamFieldRange samFieldRange(uint8_t* line, const LineTabs& tabs,
                                   uint32_t fieldIdx, uint32_t lineEnd)
{
    SamFieldRange range;
    const uint32_t prevTabPos = (uint32_t)tabs[fieldIdx - 1];
    const uint32_t currTabPos = (fieldIdx < tabs.size()) ? (uint32_t)tabs[fieldIdx] : lineEnd;
    range.start = line + prevTabPos + 1;
    range.tab = line + currTabPos;
    return range;
}
PreprocessInfo* SamCodecActuator::preprocessInfoMut() const
{
    return (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfoMut() : nullptr;
}

/* The coder for the SEQ match stream, decided once per file by the codec pre-selection; coder_bwt_cm is the fallback. */
CoderType SamCodecActuator::seqMatchCoderFor()
{
    const CoderType modeFallback =
        CoderFactory::fieldCoderForMode(SAM_SEQ, compressMode(), CoderType::BWT_CM);
    PreprocessInfo* preInfo = preprocessInfoMut();
    if (preInfo != nullptr) {
        const int32_t pinned = preInfo->seqMatchCoder.load(std::memory_order_relaxed);
        if (pinned >= 0) {
            const CoderType pinnedType = (CoderType)pinned;
            return CoderFactory::eligibleProfile(pinnedType, compressMode())
                       ? pinnedType : modeFallback;
        }
        /* No reference payload was trialled (a file with no reference): the pre-selection's verdict
           for the SEQ column itself - coder_fc by default, coder_fc/coder_bwt_cm at -l9 - applies. */
        const FieldCodecSelection* sel = preInfo->getField(SAM_SEQ);
        if (sel != nullptr && sel->status == FieldStatus::SELECTED &&
            CoderFactory::eligibleProfile(sel->selectedCoder, compressMode())) {
            return sel->selectedCoder;
        }
    }

    /* No verdict covers this file; the mode's fast vector or coder_bwt_cm. */
    return modeFallback;
}

/* Encode PNEXT's exception stream: (contentIdx, delta) pairs of records that cannot be rebuilt from their mate. */
uint32_t SamCodecActuator::writePnextExceptions(const std::vector<std::pair<uint32_t, int64_t>>& exc,
                                                  Json::Value& fieldMeta, Json::Value& metaStreams,
                                                  uint32_t& totalDstLen)
{
    Json::Value metaSubs;
    if (!exc.empty()) {
        std::shared_ptr<coder_io> excIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "PNEXT exceptions");
        std::shared_ptr<coder> excCoder =
            CoderFactory::makeProfileEncoder(SAM_PNEXT, compressMode(), excIo.get(),
                                             CoderType::BWT_CM, false, /* compressLevel */ 0);

        /* contentIdx as forward deltas, delta as zigzag varints; a constant delta moves into the meta and the column is dropped. */
        std::vector<uint8_t> excStream;
        excStream.reserve(exc.size() * 2 + 16);
        uint32_t prevOrdinal = 0;
        bool sameDelta = true;
        const int64_t firstDelta = exc[0].second;
        for (const auto& e : exc) {
            appendVarint(excStream, e.first - prevOrdinal);
            prevOrdinal = e.first;
            if (e.second != firstDelta) {
                sameDelta = false;
            }
        }
        if (sameDelta) {
            fieldMeta["exc_delta"] = (Json::Value::UInt64)zigzag64(firstDelta);
        } else {
            for (const auto& e : exc) {
                appendVarint64(excStream, zigzag64(e.second));
            }
        }
        excCoder->encode_line(excStream.data(), (uint32_t)excStream.size());
        excCoder->encode_flush();
        if (excIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode PNEXT exceptions overflow: output buffer too small");
            return -1;
        }
        metaSubs["srclen"] = (Json::Value::UInt)excStream.size();
        metaSubs["dstlen"] = excIo->data_len;
        metaSubs["coder"] = excIo->meta;
        metaSubs["sname"] = "pnextexc";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + excIo->data_len);
        totalDstLen += excIo->data_len;
        fieldMeta["exc_enc"] = "varint2";
    } else {
        fieldMeta["exc_enc"] = "none";
    }
    return totalDstLen;
}

int32_t SamCodecActuator::compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    fieldSrcLen = 0;
    uint32_t totalSrcLen = 0;

    /* PNEXT is exception-only: paired records with a mutually-mapped mate are rebuilt from the mate's POS; the rest are stored. */

    // Pass 1: collect (lineIdx -> qname, pos, pnext, flag) for every data line.
    /* Record table indexed by data-line number; QNAME kept as a slice of the block buffer. */
    struct RecInfo {
        uint32_t qoff = 0;           // QNAME slice in the block buffer
        uint32_t qlen = 0;
        int64_t pos = 0;
        int64_t pnext = 0;
        uint16_t flag = 0;
        uint32_t group = UINT32_MAX; // index into qnameGroups, when this line groups by name
        bool has = false;            // this slot was walked (a '@' line is skipped outright)
        bool valid = false;          // paired (0x1) and mate mapped (not 0x8)
        bool hasPnext = false;       // pnext != 0/ *
    };
    /* One group per distinct QNAME: valid-record count and the two lines that could pair up, through a hash -> index chain. */
    struct QnameGroup {
        uint32_t count = 0;          // valid records with this name
        uint32_t line[2] = {0, 0};   // the first two of them
        uint32_t qoff = 0;           // the name's slice, to confirm a hash hit
        uint32_t qlen = 0;
        uint32_t next = UINT32_MAX;  // next group with the same hash
    };
    const uint32_t slotCount = lineNum - headEndLine;
    std::vector<RecInfo> records(slotCount);
    std::vector<QnameGroup> qnameGroups;
    std::unordered_map<uint64_t, uint32_t> qnameHead;

    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        uint32_t contentIdx = lineIdx - headEndLine;

        // QNAME is field 0 (line start .. first tab).
        uint32_t qnameLen = contentPos[contentIdx].empty()
            ? 0 : contentPos[contentIdx][0];
        RecInfo& ri = records[contentIdx];
        ri.has = true;
        ri.qoff = lineStart;
        ri.qlen = qnameLen;

        // FLAG (field 1)
        auto flagIt = mappedFlag.find(lineIdx);
        if (flagIt != mappedFlag.end()) {
            ri.flag = flagIt->second;
        }
        // POS (field 3) is already decoded by compressPosFieldDelta into mappedPos.
        auto posIt = mappedPos.find(lineIdx);
        if (posIt != mappedPos.end()) {
            ri.pos = posIt->second;
        }
        // PNEXT (field 7): fieldLength includes trailing tab.
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();
        totalSrcLen += fieldLength;

        ri.valid = (ri.flag & 0x1) && !(ri.flag & 0x8);
        if (fieldLength > 1) {
            std::string pnStr((char*)fieldStart, fieldLength - 1);
            if (pnStr != "0" && pnStr != "*") {
                ri.pnext = (int64_t)std::stoll(pnStr);
                ri.hasPnext = true;
            }
        }
        nextMappedPos[lineIdx] = ri.hasPnext ? ri.pnext : 0;

        // Group by QNAME using the same criterion as the decoder (valid =
        // paired && mate mapped), so exception decisions match exactly.
        if (ri.valid) {
            const uint64_t h = CityHash64((const char*)line, qnameLen);
            uint32_t g = UINT32_MAX;
            auto headIt = qnameHead.find(h);
            if (headIt != qnameHead.end()) {
                for (uint32_t cand = headIt->second; cand != UINT32_MAX; cand = qnameGroups[cand].next) {
                    if (qnameGroups[cand].qlen == qnameLen &&
                        memcmp(buffer + qnameGroups[cand].qoff, line, qnameLen) == 0) {
                        g = cand;
                        break;
                    }
                }
            }
            if (g == UINT32_MAX) {
                g = (uint32_t)qnameGroups.size();
                QnameGroup ng;
                ng.qoff = lineStart;
                ng.qlen = qnameLen;
                qnameGroups.push_back(ng);
                if (headIt != qnameHead.end()) {
                    qnameGroups[g].next = headIt->second;
                    headIt->second = g;
                } else {
                    qnameHead.emplace(h, g);
                }
            }
            ri.group = g;
            QnameGroup& gr = qnameGroups[g];
            if (gr.count < 2) {
                gr.line[gr.count] = lineIdx;
            }
            gr.count++;
        }
    }

    // Pass 2: decide rebuildable vs exception.
    std::vector<std::pair<uint32_t, int64_t>> pnextExceptions; // (contentIdx, delta = pnext - pos)
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        RecInfo& ri = records[lineIdx - headEndLine];
        if (!ri.has) {
            continue;
        }

        bool rebuildable = false;
        if (ri.valid && ri.hasPnext && ri.group != UINT32_MAX) {
            const QnameGroup& gr = qnameGroups[ri.group];
            // For the decoder to uniquely locate the mate from the QNAME group,
            // this QNAME must contain exactly two mutually-mapped records.
            if (gr.count == 2) {
                const uint32_t mateLine = (gr.line[0] == lineIdx) ? gr.line[1] : gr.line[0];
                const RecInfo& mr = records[mateLine - headEndLine];
                if (mr.valid && mr.hasPnext &&
                    mr.pnext == ri.pos && ri.pnext == mr.pos) {
                    rebuildable = true;
                }
            }
        }
        if (!rebuildable) {
            /* A record with no real PNEXT (field "0"/"*") and not a mapped pair derives PNEXT 0; no exception is stored. */
            if (ri.valid || ri.hasPnext) {
                uint32_t contentIdx = lineIdx - headEndLine;
                /* Decoder reconstructs pnext = delta + pos; a record with no real PNEXT needs delta = -pos to reproduce 0. */
                int64_t delta = ri.hasPnext ? (ri.pnext - ri.pos) : (0 - (int64_t)ri.pos);
                pnextExceptions.emplace_back(contentIdx, delta);
            }
        }
    }

    // Encode the exception stream: (contentIdx, delta) pairs, like TLEN.
    Json::Value metaStreams;
    uint32_t totalDstLen = 0;

    /* (see writePnextExceptions) */
    totalDstLen = writePnextExceptions(pnextExceptions, fieldMeta, metaStreams, totalDstLen);

    // fieldSrcLen is the original PNEXT byte length (with trailing tab)
    // summed over all lines; the exception-stream size is only in the sub-stream meta.
    fieldSrcLen = totalSrcLen;
    fieldMeta["srclen"] = totalSrcLen;
    fieldMeta["dstlen"] = totalDstLen;
    fieldMeta["streams"] = metaStreams;
    fieldMeta["field"] = fieldIdx;
    fieldMeta["mode"] = "pnext_qname_rebuild";
    fieldMeta["exceptions"] = (Json::UInt64)pnextExceptions.size();

    LOG_INFO("SAM field(%d) (PNEXT) qname-rebuild compression completed: %u src -> %u dst, %u exceptions, ratio = %.2f%%",
        fieldIdx, totalSrcLen, totalDstLen, (uint32_t)pnextExceptions.size(),
        totalSrcLen ? (double)(totalDstLen * 100)/(double)totalSrcLen : 0.0);

    return (int32_t)totalDstLen;
}

int32_t SamCodecActuator::compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "POS delta");
    const CoderType archiveCoder =
        (pickedCoderFor(fieldIdx, CoderType::BWT_CM) == CoderType::ARITH)
            ? CoderType::ARITH : CoderType::BWT_CM;
    std::shared_ptr<coder> fieldCoder =
        CoderFactory::makeProfileEncoder(fieldIdx, compressMode(), fieldIo.get(), archiveCoder,
                                         true, /* compressLevel */ 0);
    /* The arithmetic backend takes the file-level delta prior to skip the model cold start; bwt_cm ignores it. */
    if (coder_arith* arithCoder = dynamic_cast<coder_arith*>(fieldCoder.get())) {
        const PreprocessInfo* preInfo =
            (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
        if (preInfo != nullptr && !preInfo->posPrior().empty()) {
            arithCoder->set_prior(preInfo->posPrior().data(),
                                  (uint32_t)preInfo->posPrior().size());
        }
    }

    fieldSrcLen = 0;
    uint32_t deltaLength = 0;
    int64_t prevPos = 0;

    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }

        uint32_t contentIdx = lineIdx - headEndLine;

        /* Detect chromosome switches via the per-line mappedChr index; the decoder replays the same comparison. */
        if (contentIdx > 0) {
            const auto curChrIt = mappedChr.find(lineIdx);
            const auto prevChrIt = mappedChr.find(lineIdx - 1);
            const uint16_t curChr = (curChrIt != mappedChr.end()) ? curChrIt->second : 0xFFFF;
            const uint16_t prevChr = (prevChrIt != mappedChr.end()) ? prevChrIt->second : 0xFFFF;
            if (curChr != prevChr) {
                prevPos = 0;
            }
        }

        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();

        int64_t delta = 0;
        if (fieldLength > 1) {
            std::string posStr = std::string((char*)fieldStart, fieldLength - 1);
            int64_t posValue = (int64_t)std::stoll(posStr);
            mappedPos[lineIdx] = posValue;
            delta = posValue - prevPos;
            fieldSrcLen += fieldLength;
            prevPos = posValue;
        } else {
            /* Empty/invalid values are treated as a delta of 0. */
            mappedPos[lineIdx] = 0;
            fieldSrcLen += 2;
        }
        /* Delta as unsigned LEB128 varint; after the baseline resets, sorted deltas are non-negative. */
        uint64_t u = (uint64_t)delta;
        uint8_t vbuf[10];
        uint32_t vlen = 0;
        do {
            uint8_t b = (uint8_t)(u & 0x7f);
            u >>= 7;
            if (u) b |= 0x80;
            vbuf[vlen++] = b;
        } while (u);
        fieldCoder->encode_line(vbuf, vlen);
        deltaLength += vlen;
    }

    fieldCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode POS delta overflow: output buffer too small");
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    fieldMeta["srclen"] = deltaLength;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;
    fieldMeta["mode"] = "pos_delta";

    LOG_INFO("SAM field(%d) (POS) delta compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

void SamCodecActuator::clearIdNumericState() {
    for (size_t i = 0; i < idNumericBufs.size(); ++i) {
        MemoryUtil::safeFree(idNumericBufs[i]);
    }
    idNumericBufs.clear();
    idNumericLens.clear();
    idNumericPos.clear();
    idNumericAcc.clear();
    idIntCoders.clear();
    idIntModes.clear();
    idIntActive.clear();
    idConstTexts.clear();
    idDictEntries.clear();
}


int32_t SamCodecActuator::compressIdFieldInAll(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create encoder for ID field whole compression
    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QNAME");
    std::shared_ptr<coder> fieldCoder = makeFieldEncoder(SAM_QNAME, samFieldDefaultCoder(SAM_QNAME, CoderType::BWT_CM), fieldIo.get(), true);

    fieldSrcLen = 0;

    // Process each line and compress ID field as whole
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;

        // Skip header lines (starting with @)
        if (buffer[lineStart] == '@') {
            continue;
        }

        // Extract ID field (from line start to first tab)
        uint8_t* idStart = buffer + lineStart;
        uint32_t idLength = contentPos[lineIdx - headEndLine][0] + 1;

        // Encode the ID field data
        if (idLength > 0) {
            fieldCoder->encode_line(idStart, idLength);
            fieldSrcLen += idLength;
        }
    }

    // Flush the encoder for ID field
    fieldCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode id field overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    Json::Value tempIdMeta;
    tempIdMeta["srclen"] = fieldSrcLen;
    tempIdMeta["dstlen"] = fieldIo->data_len;
    tempIdMeta["coder"] = fieldIo->meta;

    Json::Value streamMeta;
    streamMeta.append(tempIdMeta);
    // Set field metadata
    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["splitsym"] = "\t";
    fieldMeta["field"] = 0;

    LOG_INFO("SAM ID compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

int32_t SamCodecActuator::compressChrName(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create encoder for regular field compression
    std::shared_ptr<coder_io> chrIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "RNAME");
    std::shared_ptr<coder> chrCoder = makeFieldEncoder(fieldIdx, samFieldDefaultCoder(fieldIdx, CoderType::BWT_CM), chrIo.get(), true);

    fieldSrcLen = 0;
    uint32_t srcLen = 0;
    // Process each line and extract the current field
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;

        // Skip header lines (starting with @)
        if (*line == '@') {
            continue;
        }

        // Middle fields: between tabs
        uint32_t contentId = lineIdx - headEndLine;
        const SamFieldRange field = samFieldRange(line, contentPos[contentId], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithoutTab();

        std::string str = std::string((char*)fieldStart, fieldLength);
        uint16_t chrIndex = 0xFFFF;
        if (str == "*") {
            chrIndex = 0xFFFF;
        } else if (str == "=") {
            chrIndex = 0xFFFE;
        } else {
            chrIndex = SamInfo::getInstance().getChrNameIndex(str);
        }

        if (fieldIdx == 2) {
            mappedChr[lineIdx] = chrIndex;
        } else if (fieldIdx == 6) {
            nextMappedChr[lineIdx] = chrIndex;
        }

        if (fieldIdx == 6) {
            /* RNEXT stored as a one-byte code: 0x00 "=", 0x01 "*", 0x02..0xFD contig+2, 0xFF escape + raw uint16. */
            uint8_t code = 0;
            if (str == "=") {
                code = 0x00;
            } else if (str == "*") {
                code = 0x01;
            } else if (chrIndex >= 252) {
                code = 0xFF;
                chrCoder->encode_line(&code, 1);
                const uint8_t le[2] = {(uint8_t)(chrIndex & 0xff), (uint8_t)(chrIndex >> 8)};
                chrCoder->encode_line(le, 2);
                srcLen += 3;
                fieldSrcLen += str.length();
                continue;
            } else {
                code = (uint8_t)(chrIndex + 2);
            }
            chrCoder->encode_line(&code, 1);
            srcLen += 1;
        } else {
            // RNAME: encode the chromosome index (fixed 2 bytes).
            chrCoder->encode_line(reinterpret_cast<const uint8_t*>(&chrIndex), sizeof(chrIndex));
            srcLen += sizeof(chrIndex);
        }
        fieldSrcLen += str.length();
    }

    // Flush the encoder
    chrCoder->encode_flush();
    if (chrIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode chr name overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + chrIo->data_len);

    // Set field metadata
    fieldMeta["srclen"] = srcLen;
    fieldMeta["dstlen"] = chrIo->data_len;
    fieldMeta["coder"] = chrIo->meta;
    fieldMeta["field"] = fieldIdx;
    if (fieldIdx == 6) {
        /* RNEXT uses the one-byte code stream (see compressChrName). */
        fieldMeta["rn_enc"] = "byte";
    }

    LOG_INFO("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldIdx, fieldSrcLen, chrIo->data_len, (double)(chrIo->data_len * 100)/(double)fieldSrcLen);

    return chrIo->data_len;
 }

int32_t SamCodecActuator::compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create encoder for regular field compression
    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SAM field");
    std::shared_ptr<coder> fieldCoder = makeFieldEncoder(fieldIdx, samFieldDefaultCoder(fieldIdx, CoderType::BWT_CM), fieldIo.get(), true);

    fieldSrcLen = 0;

    /* Placeholder for a line without this field; outside the loop so coder_affix_match's
       context pointer stays valid. */
    uint8_t ch = '\n';

    // Process each line and extract the current field
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;
        // Skip header lines (starting with @)
        if (*line == '@') {
            continue;
        }

        uint32_t contentIdx = lineIdx - headEndLine;
        // Middle fields: between tabs
        if (fieldIdx > contentPos[contentIdx].size()) {
            fieldCoder->encode_line(&ch, 1);
            fieldSrcLen += 1;
            continue;
        }
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();

        // Encode the field data
        if (fieldLength > 0) {
            fieldCoder->encode_line(fieldStart, fieldLength);
            fieldSrcLen += fieldLength;

            /* Textual FLAG/POS/PNEXT must still fill the tracking maps here, keeping both sides consistent. */
            if (fieldLength > 1 && (fieldIdx == 1 || fieldIdx == 3 || fieldIdx == 7)) {
                std::string str = std::string((char*)fieldStart, fieldLength - 1);
                if (fieldIdx == 1) {
                    mappedFlag[lineIdx] = (uint16_t)std::stoll(str);
                } else if (fieldIdx == 3) {
                    mappedPos[lineIdx] = (int64_t)std::stoll(str);
                } else if (fieldIdx == 7) {
                    nextMappedPos[lineIdx] = (int64_t)std::stoll(str);
                }
            }
        }
    }

    // Flush the encoder for this field
    fieldCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode field overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    // Set field metadata
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;
    fieldMeta["mode"] = "string";

    LOG_INFO("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

/* The OPTION field (column 12 on) is compressed by CRAM-style tag columnization. */
/* Parse the OPTION column: tag dictionary (first-seen order), each line's tag sequence, each tag's value column.
   Takes the block views explicitly so the file-level layout trial can share it. */
static int32_t parseOptionTags(const uint8_t* buffer, const std::vector<size_t>& npos,
                               uint32_t headEndLine, const std::vector<LineTabs>& contentPos,
                               OptionParseState& st)
{
    uint32_t lineNum = npos.size();
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;
        const uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        uint32_t contentIdx = lineIdx - headEndLine;
        /* A record may carry fewer than the 12 columns (no OPTION). contentPos stores every
           '\t' and the trailing '\n', so position 10 is the 11th separator only when an OPTION
           column follows; otherwise there is nothing to parse. */
        const LineTabs& tabs = contentPos[contentIdx];
        if (tabs.size() < 11 || (uint64_t)tabs[10] >= (uint64_t)lineEnd) {
            st.recIds.push_back({});
            continue;
        }
        uint32_t fieldStartPos = (uint32_t)tabs[10] + 1;
        uint32_t fieldLen = lineEnd - fieldStartPos;
        if (fieldLen == 0) {
            st.recIds.push_back({});
            continue;
        }
        st.srcLen += fieldLen;

        std::vector<uint8_t> ids;
        ids.reserve(fieldLen / 5 + 1);
        const uint8_t* p = line + fieldStartPos;
        uint32_t segStart = 0;
        for (uint32_t i = 0; i <= fieldLen; ++i) {
            if (i != fieldLen && p[i] != '\t') {
                continue;
            }

            const uint32_t segLen = i - segStart;
            if (segLen > 0 && p[segStart] != '\n') {
                const uint8_t* ts = p + segStart;
                const uint8_t* segEnd = ts + segLen;
                const uint8_t* colon1 = nullptr;
                const uint8_t* colon2 = nullptr;
                if (segLen > 1) {
                    colon1 = static_cast<const uint8_t*>(
                        memchr(ts + 1, ':', (size_t)(segLen - 1)));
                }
                if (colon1 != nullptr) {
                    colon2 = static_cast<const uint8_t*>(
                        memchr(colon1 + 1, ':', (size_t)(segEnd - (colon1 + 1))));
                }
                if (colon2 == nullptr) {
                    /* Not TAG:TYPE:VALUE: the per-tag layout cannot represent it losslessly. */
                    return -1;
                }

                const size_t nameLen = (size_t)(colon1 - ts);
                const size_t typeLen = (size_t)(colon2 - colon1 - 1);
                const char* namePtr = reinterpret_cast<const char*>(ts);
                const char* typePtr = reinterpret_cast<const char*>(colon1 + 1);

                /* The tag set is small; a direct comparison avoids a temporary
                   std::string and a hash-map lookup for every tag of every line. */
                int id = -1;
                for (size_t t = 0; t < st.tagDict.size(); ++t) {
                    const std::string& known = st.tagDict[t].first;
                    if (known.size() == nameLen &&
                        memcmp(known.data(), namePtr, nameLen) == 0) {
                        id = (int)t;
                        break;
                    }
                }
                if (id < 0) {
                    id = (int)st.tagDict.size();
                    st.tagDict.emplace_back(std::string(namePtr, nameLen),
                                            std::string(typePtr, typeLen));
                    st.tagVals.emplace_back();
                } else {
                    const std::string& knownType = st.tagDict[id].second;
                    if (knownType.size() != typeLen ||
                        memcmp(knownType.data(), typePtr, typeLen) != 0) {
                        /* Same tag name with two types cannot share one column. */
                        return -1;
                    }
                }
                ids.push_back((uint8_t)id);
                st.tagVals[id].emplace_back(reinterpret_cast<const char*>(colon2 + 1),
                                            (size_t)(segEnd - (colon2 + 1)));
            }
            segStart = i + 1;
        }
        st.recIds.push_back(std::move(ids));
    }
    return 0;
}

int32_t SamCodecActuator::collectOptionTags(OptionParseState& st)
{
    return parseOptionTags(inBlockPtr->getBuffer(), inBlockPtr->getNpos(),
                           (uint32_t)headEndLine, contentPos, st);
}

/* Whether the block has a tag segment starting with `seg` (e.g. "MD:Z:"). A cheap byte scan:
   the OPTION column always follows a tab (after QUAL) and every tag starts after a tab, and a
   single-value field cannot contain a tab, so this only matches real OPTION segments. Used to
   gate the per-tag layout before paying for the per-record tag parse.

   Only a prefix is scanned: a block's tag set is uniform, and a scan of the whole block just to
   find that a tag is absent is what made no-MD files (e.g. con_sorted.sam) slower. */
static bool blockHasTagSegment(const uint8_t* data, size_t len, const char* seg)
{
    static const size_t kScanBytes = 1u << 20;   /* 1 MB: many records even with a large header */
    const size_t n = strlen(seg);
    if (data == nullptr || len < n) {
        return false;
    }
    const size_t end = (len < kScanBytes) ? len : kScanBytes;
    for (size_t i = 0; i + n <= end; ++i) {
        if (data[i] == '\t' && memcmp(data + i + 1, seg, n) == 0) {
            return true;
        }
    }
    return false;
}

/* SAM/BAM integer aux-tag types (c/C/s/S/i/I); every other type is a string/blob. */
static bool samTagIsInteger(const std::string& type)
{
    if (type.size() != 1) {
        return false;
    }
    switch (type[0]) {
        case 'c': case 'C': case 's': case 'S': case 'i': case 'I':
            return true;
        default:
            return false;
    }
}

/* CRAM-style OPTION serialization: one column per tag, each encoded on its own (integer tags try
   text/varint/delta-varint; string tags are line-oriented text; MD/NM are regenerated). */

/* Recompute MD/NM from the reference + CIGAR + SEQ (as CRAM does) and mark, per value,
   whether the stored one equals the regenerated one.  Takes the block views explicitly so
   the file-level layout trial can share it. */
static void deriveOptionTags(uint8_t* buffer, const std::vector<size_t>& npos,
                             uint32_t headEndLine, const std::vector<LineTabs>& contentPos,
                             Reference* ref, const OptionParseState& opt, int mdId, int nmId,
                             std::vector<std::vector<uint8_t>>& derivedFlags,
                             std::vector<std::string>& derivedExcs,
                             std::vector<bool>& derivedReady)
{
        std::vector<uint8_t> refScratch;
        std::vector<CigarOp> ops;
        size_t mdCur = 0, nmCur = 0;
        derivedFlags[mdId].assign(opt.tagVals[mdId].size(), 1);
        derivedFlags[nmId].assign(opt.tagVals[nmId].size(), 1);
        for (size_t r = 0; r < opt.recIds.size(); ++r) {
            bool recHasMd = false, recHasNm = false;
            for (uint8_t id : opt.recIds[r]) {
                if ((int)id == mdId) recHasMd = true;
                if ((int)id == nmId) recHasNm = true;
            }
            if (!recHasMd && !recHasNm) {
                continue;
            }

            const uint32_t lineIdx = headEndLine + (uint32_t)r;
            const uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
            const uint32_t lineEnd = (uint32_t)npos[lineIdx] - lineStart;
            uint8_t* line = buffer + lineStart;
            const LineTabs& tabs = contentPos[r];
            const SamFieldRange fl = samFieldRange(line, tabs, 1, lineEnd);
            const SamFieldRange rn = samFieldRange(line, tabs, 2, lineEnd);
            const SamFieldRange po = samFieldRange(line, tabs, 3, lineEnd);
            const SamFieldRange ci = samFieldRange(line, tabs, 5, lineEnd);
            const SamFieldRange sq = samFieldRange(line, tabs, 9, lineEnd);

            uint32_t flagVal = 0;
            for (uint32_t i = 0; i < fl.lengthWithoutTab(); ++i) {
                const char c = (char)fl.start[i];
                if (c >= '0' && c <= '9') flagVal = flagVal * 10 + (uint32_t)(c - '0');
            }
            uint64_t posVal = 0;
            for (uint32_t i = 0; i < po.lengthWithoutTab(); ++i) {
                const char c = (char)po.start[i];
                if (c >= '0' && c <= '9') posVal = posVal * 10 + (uint64_t)(c - '0');
            }
            const std::string chrName((const char*)rn.start, rn.lengthWithoutTab());
            const uint16_t chrId = SamInfo::getInstance().getChrNameIndex(chrName);

            std::string derivedMd;
            int64_t derivedNm = 0;
            bool derivedOk = false;
            if (chrId != 0xFFFF && chrId != 0xFFFE && posVal != 0 &&
                ci.lengthWithoutTab() > 0 && sq.lengthWithoutTab() > 0 &&
                !(sq.lengthWithoutTab() == 1 && sq.start[0] == '*')) {
                const int64_t chrStart = SamInfo::getInstance().getPositionByIndex(chrId);
                if (chrStart >= 0) {
                    ops.clear();
                    parseCigarOps(ci.start, ci.lengthWithoutTab(), ops);
                    derivedOk = deriveMdNm(ref, chrStart, posVal, (uint16_t)flagVal, ops,
                                           sq.start, sq.lengthWithoutTab(), refScratch, derivedMd,
                                           derivedNm);
                }
            }
            if (recHasMd) {
                const size_t idx = mdCur++;
                if (derivedOk && derivedMd == opt.tagVals[mdId][idx]) {
                    derivedFlags[mdId][idx] = 0;
                } else {
                    derivedExcs[mdId].append(opt.tagVals[mdId][idx]);
                    derivedExcs[mdId].push_back('\n');
                }
            }
            if (recHasNm) {
                const size_t idx = nmCur++;
                if (derivedOk && std::to_string((long long)derivedNm) == opt.tagVals[nmId][idx]) {
                    derivedFlags[nmId][idx] = 0;
                } else {
                    derivedExcs[nmId].append(opt.tagVals[nmId][idx]);
                    derivedExcs[nmId].push_back('\n');
                }
            }
        }
    derivedReady[mdId] = true;
    derivedReady[nmId] = true;
}

/* Build one tag's layout candidates: text, the integer forms, and the regenerated form. */
static void buildOptionCandidates(const std::vector<std::string>& vals, bool isInt,
                                  bool derivedReady,
                                  const std::vector<uint8_t>& derivedFlags,
                                  const std::string& derivedExc,
                                  std::vector<std::pair<const char*, std::vector<uint8_t>>>& out)
{
    out.clear();
    std::vector<uint8_t> textCol;
    for (size_t k = 0; k < vals.size(); ++k) {
        textCol.insert(textCol.end(), vals[k].begin(), vals[k].end());
        textCol.push_back('\n');
    }
    out.emplace_back("text", std::move(textCol));
    if (isInt) {
        std::vector<uint8_t> rawCol, deltaCol;
        rawCol.reserve(vals.size() * 2 + 16);
        deltaCol.reserve(vals.size() * 2 + 16);
        int64_t prev = 0;
        for (size_t k = 0; k < vals.size(); ++k) {
            const int64_t v = strtoll(vals[k].c_str(), nullptr, 10);
            appendVarint64(rawCol, (uint64_t)v);
            appendVarint64(deltaCol, zigzag64(v - prev));
            prev = v;
        }
        out.emplace_back("varint", std::move(rawCol));
        out.emplace_back("dvarint", std::move(deltaCol));
    }
    if (derivedReady) {
        /* flags (one per value) followed by the exception values, line by line. */
        std::vector<uint8_t> derivedCol;
        derivedCol.reserve(derivedFlags.size() + derivedExc.size());
        derivedCol.insert(derivedCol.end(), derivedFlags.begin(), derivedFlags.end());
        derivedCol.insert(derivedCol.end(), derivedExc.begin(), derivedExc.end());
        out.emplace_back("derived", std::move(derivedCol));
    }
}

/* Compressed size of one layout, encoded exactly as a block would encode it.  Returns 0 when
   the scratch buffer could not hold it (the tag is then left unplanned). */
static uint32_t trialLayoutSize(uint8_t mode, const std::vector<uint8_t>& blob,
                                uint8_t* scratch, uint32_t scratchLen)
{
    if (scratch == nullptr || scratchLen == 0) {
        return 0;
    }
    coder_io io(scratch, (int32_t)scratchLen, nullptr, "OPTION layout trial");
    std::shared_ptr<coder> c = CoderFactory::makeProfileEncoder(SAM_FIELD_COUNT, mode, &io,
                                                                CoderType::BWT_CM, false, 0);
    c->encode_line(blob.data(), (uint32_t)blob.size());
    c->encode_flush();
    return (io.err == coder_io::IO_OK) ? (uint32_t)io.data_len : 0;
}

/* File head line count and per-line tab positions for a block, the views parseOptionTags and
   deriveOptionTags take.  Returns false when a line is not a SAM record. */
static bool buildOptionLineViews(const uint8_t* buffer, const std::vector<size_t>& npos,
                                 std::vector<LineTabs>& contentPos, uint32_t& headEndLine)
{
    contentPos.clear();
    contentPos.reserve(npos.size());
    headEndLine = 0;
    bool header = true;
    for (size_t lineIdx = 0; lineIdx < npos.size(); ++lineIdx) {
        const uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
        const uint32_t lineEnd = (uint32_t)npos[lineIdx] - lineStart;
        if (lineEnd == 0) {
            return false;   /* empty line */
        }
        if (header && buffer[lineStart] == '@') {
            ++headEndLine;
            continue;
        }
        header = false;
        LineTabs tabs;
        for (uint32_t i = 0; i < lineEnd; ++i) {
            const uint8_t c = buffer[lineStart + i];
            if (c == '\t' || c == '\n') {
                tabs.push_back(i);
                if (tabs.size() == LineTabs::CAP) {
                    break;
                }
            }
        }
        if (tabs.size() < 10) {
            return false;   /* not a SAM record */
        }
        contentPos.push_back(tabs);
    }
    return !contentPos.empty();
}

/* Build the file-level OPTION layout plan: trial every tag's layouts once and record the
   winner.  Runs once, on a reader-side thread. */
int32_t buildSamOptionPlan(uint8_t* buffer, uint32_t dataLen, const std::vector<size_t>& npos,
                           Reference* ref, uint8_t mode,
                           std::unordered_map<std::string, std::string>& out)
{
    out.clear();
    if (ref == nullptr || buffer == nullptr || dataLen == 0 || npos.empty()) {
        return 0;
    }
    if (!blockHasTagSegment(buffer, (size_t)dataLen, "MD:Z:") ||
        !blockHasTagSegment(buffer, (size_t)dataLen, "NM:i:")) {
        return 0;
    }

    std::vector<LineTabs> contentPos;
    uint32_t headEndLine = 0;
    if (!buildOptionLineViews(buffer, npos, contentPos, headEndLine)) {
        return 0;
    }
    OptionParseState opt;
    if (parseOptionTags(buffer, npos, headEndLine, contentPos, opt) != 0) {
        return 0;
    }

    const uint32_t nTag = (uint32_t)opt.tagDict.size();
    int mdId = -1, nmId = -1;
    for (uint32_t t = 0; t < nTag; ++t) {
        if (opt.tagDict[t].first == "MD" && opt.tagDict[t].second == "Z") mdId = (int)t;
        if (opt.tagDict[t].first == "NM" &&
            (opt.tagDict[t].second == "i" || opt.tagDict[t].second == "I")) nmId = (int)t;
    }
    if (mdId < 0 || nmId < 0) {
        return 0;
    }

    std::vector<std::vector<uint8_t>> derivedFlags(nTag);
    std::vector<std::string> derivedExcs(nTag);
    std::vector<bool> derivedReady(nTag, false);
    deriveOptionTags(buffer, npos, headEndLine, contentPos, ref, opt, mdId, nmId,
                     derivedFlags, derivedExcs, derivedReady);

    /* Same gate as the block path: most MD/NM values must be regenerable, or the per-tag
       layout is not used at all and there is nothing to plan. */
    size_t total = 0, derived = 0;
    for (uint32_t t = 0; t < nTag; ++t) {
        if (!derivedReady[t]) {
            continue;
        }
        total += derivedFlags[t].size();
        for (uint8_t f : derivedFlags[t]) {
            if (f == 0) {
                ++derived;
            }
        }
    }
    if (total == 0 || derived * 2 < total) {
        return 0;
    }

    /* Generous scratch: a layout whose encode does not fit leaves the tag unplanned. */
    const size_t scratchLen = (size_t)dataLen * 2 + (1u << 20);
    std::vector<uint8_t> scratch(scratchLen);
    int32_t planned = 0;
    for (uint32_t t = 0; t < nTag; ++t) {
        std::vector<std::pair<const char*, std::vector<uint8_t>>> cands;
        buildOptionCandidates(opt.tagVals[t], samTagIsInteger(opt.tagDict[t].second),
                              derivedReady[t], derivedFlags[t], derivedExcs[t], cands);
        std::vector<uint32_t> sizes(cands.size(), 0);
        uint32_t best = UINT32_MAX;
        bool failed = false;
        for (size_t c = 0; c < cands.size(); ++c) {
            sizes[c] = trialLayoutSize(mode, cands[c].second, scratch.data(),
                                       (uint32_t)scratchLen);
            if (sizes[c] == 0) {
                failed = true;
                break;      /* no plan for this tag */
            }
            if (sizes[c] < best) {
                best = sizes[c];
            }
        }
        if (failed) {
            continue;
        }
        /* The winner plus any exact tie: a tie must still be decided per block. */
        std::string entry;
        for (size_t c = 0; c < cands.size(); ++c) {
            if (sizes[c] == best) {
                if (!entry.empty()) {
                    entry.push_back(',');
                }
                entry += cands[c].first;
            }
        }
        std::string key = opt.tagDict[t].first;
        key.push_back('\x1f');
        key += opt.tagDict[t].second;
        out[std::move(key)] = std::move(entry);
        ++planned;
    }
    return planned;
}

int32_t SamCodecActuator::compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {

    /* Gate the per-tag layout on a usable reference and the two regenerable tags *before* the
       per-record parse. That parse allocates a string per tag per record; a file with many tags
       but no MD/NM (e.g. con_sorted.sam) would otherwise pay it and then fall back to the whole
       column anyway. */
    if (pRefeGene == nullptr || !referenceUsableForSeq()) {
        return -1;
    }
    {
        const uint8_t* data = inBlockPtr->getBuffer();
        const size_t len = (size_t)inBlockPtr->getDataLen();
        if (!blockHasTagSegment(data, len, "MD:Z:") ||
            !blockHasTagSegment(data, len, "NM:i:")) {
            return -1;
        }
    }

    /* Parse OPTION line by line: see collectOptionTags. */
    OptionParseState opt;
    if (collectOptionTags(opt) != 0) {
        return -1;
    }
    fieldSrcLen = opt.srcLen;

    const uint32_t nTag = (uint32_t)opt.tagDict.size();

    /* The block scan above is only a hint; confirm the parsed dictionary really carries MD:Z and
       NM:i, since the layout's win comes from regenerating them. */
    int mdId = -1, nmId = -1;
    for (uint32_t t = 0; t < nTag; ++t) {
        if (opt.tagDict[t].first == "MD" && opt.tagDict[t].second == "Z") mdId = (int)t;
        if (opt.tagDict[t].first == "NM" &&
            (opt.tagDict[t].second == "i" || opt.tagDict[t].second == "I")) nmId = (int)t;
    }
    if (mdId < 0 || nmId < 0) {
        return -1;
    }

    /* File-level layout plan (see buildSamOptionPlan): built before this block was dispatched,
       so waiting here costs nothing in normal runs.  Without a plan every tag is trialled. */
    std::unordered_map<std::string, std::string> layoutPlan;
    if (pbgzEngine != nullptr) {
        SamOptionPlan& plan = pbgzEngine->getSamOptionPlan();
        std::unique_lock<std::mutex> lock(plan.mutex);
        if (plan.started) {
            plan.cv.wait(lock, [&plan] { return plan.ready; });
            layoutPlan = plan.layout;
        }
    }

    Json::Value streams(Json::arrayValue);
    uint32_t totalDst = 0;

    Json::Value tags(Json::arrayValue);
    for (uint32_t t = 0; t < nTag; ++t) {
        Json::Value e(Json::arrayValue);
        e.append(opt.tagDict[t].first);
        e.append(opt.tagDict[t].second);
        tags.append(e);
    }

    /* MD:Z and NM:i are recomputed from the reference + CIGAR + SEQ (as CRAM does), so only the
       records whose stored value differs from the recomputed one have to carry it. Each derived
       column is a byte per value (0 = derive, 1 = stored) followed by the exception values. */
    std::vector<std::vector<uint8_t>> derivedFlags(nTag);
    std::vector<std::string> derivedExcs(nTag);
    std::vector<bool> derivedReady(nTag, false);
    deriveOptionTags(inBlockPtr->getBuffer(), inBlockPtr->getNpos(), (uint32_t)headEndLine,
                     contentPos, pRefeGene, opt, mdId, nmId, derivedFlags, derivedExcs,
                     derivedReady);

    /* If most MD/NM values could not be regenerated (an aligner whose MD/NM differ from calmd),
       the per-tag layout loses its advantage; fall back to the plain column now, before writing
       anything. */
    {
        size_t total = 0, derived = 0;
        for (uint32_t t = 0; t < nTag; ++t) {
            if (!derivedReady[t]) continue;
            total += derivedFlags[t].size();
            for (uint8_t f : derivedFlags[t]) {
                if (f == 0) ++derived;
            }
        }
        if (total > 0 && derived * 2 < total) {
            return -1;
        }
    }

    /* Each column is encoded line by line with delimiters and compressed with bwt_cm's line-wise mode. */

    /* id column: one line per record, ids comma-separated, line terminated by '\n'. */
    {
        std::vector<uint8_t> idStream;
        for (size_t r = 0; r < opt.recIds.size(); ++r) {
            for (size_t k = 0; k < opt.recIds[r].size(); ++k) {
                if (k) idStream.push_back(',');
                uint32_t v = opt.recIds[r][k];
                char buf[4];
                int n = snprintf(buf, sizeof(buf), "%u", v);
                for (int b = 0; b < n; ++b) idStream.push_back((uint8_t)buf[b]);
            }
            idStream.push_back('\n');
        }
        std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "OPTION ids");
        std::shared_ptr<coder> c =
            CoderFactory::makeProfileEncoder(SAM_FIELD_COUNT, compressMode(), io.get(),
                                             CoderType::BWT_CM, true, /* compressLevel */ 0);
        size_t pos = 0;
        for (size_t r = 0; r < opt.recIds.size(); ++r) {
            /* Feed line by line to preserve bwt_cm's line-mode semantics. */
            size_t start = pos;
            while (pos < idStream.size() && idStream[pos] != '\n') pos++;
            pos++;
            if (pos > start) c->encode_line(idStream.data() + start, (uint32_t)(pos - start));
        }
        c->encode_flush();
        if (io->err != coder_io::IO_OK) return -1;
        Json::Value s;
        s["sname"] = "ids";
        s["srclen"] = (Json::Value::UInt)idStream.size();
        s["dstlen"] = (Json::Value::UInt)io->data_len;
        s["coder"] = io->meta;
        streams.append(s);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + io->data_len);
        totalDst += io->data_len;
    }

    /* Each tag gets its own column. Integer tags (AS/NM/XS/...) are tried as decimal text,
       raw varint and delta+zigzag varint, and the smallest encoding wins; string tags (MD/XA/
       SA/...) stay line-oriented. This is the CRAM layout - per-tag series with a per-tag
       method choice - instead of one mixed value stream. */
    auto encodeColumn = [&](const std::vector<uint8_t>& blob, std::vector<uint8_t>& bytes,
                            Json::Value& coderMeta, uint32_t& dst) -> bool {
        const int64_t saved = outBlockPtr->getDataLen();
        std::shared_ptr<coder_io> io =
            makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "OPTION tag");
        std::shared_ptr<coder> c =
            CoderFactory::makeProfileEncoder(SAM_FIELD_COUNT, compressMode(), io.get(),
                                             CoderType::BWT_CM, false, /* compressLevel */ 0);
        /* A byte-stream coder, so feeding the whole blob is correct. */
        c->encode_line(blob.data(), (uint32_t)blob.size());
        c->encode_flush();
        const bool ok = (io->err == coder_io::IO_OK);
        dst = ok ? (uint32_t)io->data_len : 0;
        if (ok) {
            const uint8_t* base = outBlockPtr->getBuffer() + saved;
            bytes.assign(base, base + dst);
            coderMeta = io->meta;
        }
        outBlockPtr->setDataLen(saved);
        return ok;
    };

    /* Each tag's layouts are compared by compressing them in full.  The coder emits a
       whole BWT block only on flush, so a partial encode cannot give a smaller size and
       any estimate here could change the chosen layout (and the archive size). */
    for (uint32_t t = 0; t < nTag; ++t) {
        const std::vector<std::string>& vals = opt.tagVals[t];

        std::vector<std::pair<const char*, std::vector<uint8_t>>> candidates;
        buildOptionCandidates(vals, samTagIsInteger(opt.tagDict[t].second), derivedReady[t],
                              derivedFlags[t], derivedExcs[t], candidates);

        const char* bestEnc = "text";
        uint32_t bestDst = UINT32_MAX, bestSrc = 0;
        Json::Value bestCoder;
        std::vector<uint8_t> bestBytes;

        /* Layouts the file-level plan still wants compared here: the winner plus any exact
           tie.  A tag without a plan entry (or a block where none of the planned layouts can
           be built) is trialled in full. */
        std::vector<const char*> planned;
        if (!layoutPlan.empty()) {
            std::string key = opt.tagDict[t].first;
            key.push_back('\x1f');
            key += opt.tagDict[t].second;
            auto it = layoutPlan.find(key);
            if (it != layoutPlan.end()) {
                size_t pos = 0;
                while (pos <= it->second.size()) {
                    const size_t comma = it->second.find(',', pos);
                    std::string name = it->second.substr(
                        pos, (comma == std::string::npos ? it->second.size() : comma) - pos);
                    if (!name.empty()) {
                        for (auto& cand : candidates) {
                            if (name == cand.first) {
                                planned.push_back(cand.first);
                                break;
                            }
                        }
                    }
                    if (comma == std::string::npos) {
                        break;
                    }
                    pos = comma + 1;
                }
            }
        }
        for (auto& cand : candidates) {
            if (!planned.empty()) {
                bool trialIt = false;
                for (const char* name : planned) {
                    if (strcmp(name, cand.first) == 0) {
                        trialIt = true;
                        break;
                    }
                }
                if (!trialIt) {
                    continue;
                }
            }
            std::vector<uint8_t> bytes;
            Json::Value coder;
            uint32_t dst = 0;
            if (!encodeColumn(cand.second, bytes, coder, dst)) {
                continue;
            }
            if (dst < bestDst) {
                bestDst = dst;
                bestSrc = (uint32_t)cand.second.size();
                bestEnc = cand.first;
                bestCoder = coder;
                bestBytes = std::move(bytes);
            }
        }
        if (bestDst == UINT32_MAX) {
            return -1;
        }

        const int64_t saved = outBlockPtr->getDataLen();
        if (bestDst > 0) {
            memcpy(outBlockPtr->getBuffer() + saved, bestBytes.data(), bestDst);
        }
        outBlockPtr->setDataLen(saved + bestDst);

        Json::Value s;
        s["sname"] = "tag";
        s["tag"] = opt.tagDict[t].first;
        s["type"] = opt.tagDict[t].second;
        s["enc"] = bestEnc;
        s["srclen"] = (Json::Value::UInt)bestSrc;
        s["dstlen"] = (Json::Value::UInt)bestDst;
        s["coder"] = bestCoder;
        if (std::string(bestEnc) == "derived") {
            s["n"] = (Json::Value::UInt)vals.size();   /* number of flag bytes at the blob's head */
        }
        streams.append(s);
        totalDst += bestDst;
        LOG_DEBUG("OPTION tag %s:%s enc=%s raw=%u dst=%u",
                  opt.tagDict[t].first.c_str(), opt.tagDict[t].second.c_str(), bestEnc,
                  (unsigned)bestSrc, (unsigned)bestDst);
    }

    fieldMeta["mode"] = "tag_split";
    fieldMeta["tags"] = tags;
    fieldMeta["streams"] = streams;
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = totalDst;
    fieldMeta["field"] = 11;

    LOG_INFO("SAM OPTION tag-split compression completed: %u bytes -> %u bytes, %u tags, ratio = %.2f%%",
        fieldSrcLen, totalDst, nTag, (double)(totalDst * 100) / (double)(fieldSrcLen ? fieldSrcLen : 1));
    return (int32_t)totalDst;
}

/* The OPTION column with tag name/type left out of the stream: a per-block table of tag
   sequences in the meta, plus one index and one line of values per record. Returns false
   (writes nothing) when the data does not suit it. */
bool SamCodecActuator::compressOptionStrip(uint32_t fieldIdx, uint32_t& fieldSrcLen,
                                           Json::Value& fieldMeta, uint32_t& dstLen)
{
    /* A block's tag set is small and constant; the plain layout is for the rest. */
    const size_t kMaxSequences = 4096;

    struct SeqEntry {
        std::string name;      /* empty when verbatim */
        std::string type;
        bool verbatim = false;
    };
    /* One line's segments as offsets into the block buffer, so a record costs no allocation. */
    struct SegView {
        uint32_t nameOff = 0, nameLen = 0;
        uint32_t typeOff = 0, typeLen = 0;
        uint32_t valueOff = 0, valueLen = 0;
        uint32_t segOff = 0, segLen = 0;   /* the whole segment, for a verbatim entry */
        bool verbatim = false;
    };

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    const uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::vector<std::vector<SeqEntry>> sequences;
    std::vector<uint8_t> indexText;
    std::vector<uint8_t> valueText;
    std::vector<SegView> segs;
    uint32_t srcLen = 0;

    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        const uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        const uint32_t lineEnd = npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        const uint32_t contentIdx = lineIdx - headEndLine;

        /* The same field view compressRegularField feeds: the record's fields 12.. as one column, or a bare newline. */
        const uint8_t* text = nullptr;
        uint32_t len = 0;
        bool hasField = true;
        if ((uint32_t)fieldIdx > contentPos[contentIdx].size()) {
            static const uint8_t kNewline = '\n';
            text = &kNewline;
            len = 1;
            hasField = false;
        } else {
            const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
            text = field.start;
            len = field.lengthWithTab();
            if (len == 0) {
                hasField = false;
                static const uint8_t kNewline = '\n';
                text = &kNewline;
                len = 1;
            }
        }
        srcLen += len;

        /* Blocks split at line ends, so the column ends with the record's newline (left implicit here). */
        if (len == 0 || text[len - 1] != '\n') {
            return false;
        }
        --len;

        /* Split the line into tab-separated segments; a present but empty last column is one empty verbatim segment. */
        segs.clear();
        if (len == 0 && hasField) {
            SegView v;
            v.verbatim = true;
            segs.push_back(v);
        } else if (len > 0) {
            uint32_t segStart = 0;
            for (uint32_t i = 0; i <= len; ++i) {
                if (i != len && text[i] != '\t') {
                    continue;
                }
                SegView v;
                v.segOff = segStart;
                v.segLen = i - segStart;
                const uint8_t* s = text + segStart;
                uint32_t c1 = UINT32_MAX, c2 = UINT32_MAX;
                for (uint32_t j = 1; j < v.segLen; ++j) {
                    if (s[j] == ':' && c1 == UINT32_MAX) {
                        c1 = j;
                    } else if (s[j] == ':' && c2 == UINT32_MAX) {
                        c2 = j;
                        break;
                    }
                }
                if (c2 == UINT32_MAX) {
                    v.verbatim = true;
                } else {
                    v.nameOff = segStart;
                    v.nameLen = c1;
                    v.typeOff = segStart + c1 + 1;
                    v.typeLen = c2 - c1 - 1;
                    v.valueOff = segStart + c2 + 1;
                    v.valueLen = v.segLen - c2 - 1;
                }
                segs.push_back(v);
                segStart = i + 1;
            }
        }

        /* Find this line's sequence by comparing segments. */
        int32_t seqIdx = -1;
        for (size_t s = 0; s < sequences.size() && seqIdx < 0; ++s) {
            const std::vector<SeqEntry>& seq = sequences[s];
            if (seq.size() != segs.size()) {
                continue;
            }
            bool same = true;
            for (size_t k = 0; k < segs.size(); ++k) {
                const SeqEntry& e = seq[k];
                if (e.verbatim != segs[k].verbatim) {
                    same = false;
                    break;
                }
                if (e.verbatim) {
                    continue;   /* the text is carried by the value stream per record */
                }
                if (e.name.size() != segs[k].nameLen ||
                    memcmp(e.name.data(), text + segs[k].nameOff, segs[k].nameLen) != 0 ||
                    e.type.size() != segs[k].typeLen ||
                    memcmp(e.type.data(), text + segs[k].typeOff, segs[k].typeLen) != 0) {
                    same = false;
                    break;
                }
            }
            if (same) {
                seqIdx = (int32_t)s;
            }
        }
        if (seqIdx < 0) {
            if (sequences.size() >= kMaxSequences) {
                return false;
            }
            std::vector<SeqEntry> seq;
            seq.reserve(segs.size());
            for (const SegView& v : segs) {
                SeqEntry e;
                e.verbatim = v.verbatim;
                if (!v.verbatim) {
                    e.name.assign((const char*)text + v.nameOff, v.nameLen);
                    e.type.assign((const char*)text + v.typeOff, v.typeLen);
                }
                seq.push_back(std::move(e));
            }
            sequences.push_back(std::move(seq));
            seqIdx = (int32_t)(sequences.size() - 1);
        }

        char idxBuf[16];
        const int idxLen = snprintf(idxBuf, sizeof(idxBuf), "%d\n", (int)seqIdx);
        indexText.insert(indexText.end(), idxBuf, idxBuf + idxLen);

        /* The values in sequence order, tab-separated, ending with the newline. */
        for (const SegView& v : segs) {
            if (v.verbatim) {
                valueText.insert(valueText.end(), text + v.segOff, text + v.segOff + v.segLen);
            } else {
                valueText.insert(valueText.end(), text + v.valueOff, text + v.valueOff + v.valueLen);
            }
            valueText.push_back('\t');
        }
        if (!segs.empty()) {
            valueText.back() = '\n';
        } else {
            valueText.push_back('\n');
        }
    }

    /* Commit: both sub-streams written before the meta; a failure rolls the block back. */
    const uint32_t savedDst = outBlockPtr->getDataLen();
    Json::Value streams(Json::arrayValue);
    uint32_t totalDst = 0;
    bool ok = true;
    {
        const char* const streamNames[2] = { "sidx", "vals" };
        const std::vector<uint8_t>* const texts[2] = { &indexText, &valueText };
        for (int i = 0; i < 2 && ok; ++i) {
            std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(),
                                                       outBlockPtr->getRemain(), streamNames[i]);
            std::shared_ptr<coder> c = makeFieldEncoder(fieldIdx,
                                                        samFieldDefaultCoder(fieldIdx, CoderType::BWT_CM),
                                                        io.get(), true);
            if (c == nullptr) {
                ok = false;
                break;
            }
            c->encode_line(texts[i]->data(), (uint32_t)texts[i]->size());
            c->encode_flush();
            if (io->err != coder_io::IO_OK) {
                ok = false;
                break;
            }
            Json::Value s;
            s["sname"] = streamNames[i];
            s["srclen"] = (Json::Value::UInt)texts[i]->size();
            s["dstlen"] = (Json::Value::UInt)io->data_len;
            s["coder"] = io->meta;
            streams.append(s);
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + io->data_len);
            totalDst += io->data_len;
        }
    }
    if (!ok) {
        outBlockPtr->setDataLen(savedDst);
        return false;
    }

    Json::Value seqJson(Json::arrayValue);
    for (const std::vector<SeqEntry>& seq : sequences) {
        Json::Value segJson(Json::arrayValue);
        for (const SeqEntry& e : seq) {
            Json::Value one(Json::arrayValue);
            one.append(e.verbatim ? 1 : 0);
            one.append(e.name);
            one.append(e.type);
            segJson.append(one);
        }
        seqJson.append(segJson);
    }

    fieldMeta["mode"] = "opt_strip";
    fieldMeta["seqs"] = seqJson;
    fieldMeta["streams"] = streams;
    fieldMeta["srclen"] = srcLen;
    fieldMeta["dstlen"] = totalDst;
    fieldMeta["field"] = fieldIdx;
    fieldSrcLen = srcLen;
    dstLen = totalDst;

    LOG_INFO("SAM OPTION tag-strip compression: %u bytes -> %u bytes (%zu sequences, %zu records), ratio = %.2f%%",
             srcLen, totalDst, sequences.size(), (size_t)(lineNum - headEndLine),
             (double)(totalDst * 100) / (double)(srcLen ? srcLen : 1));
    return true;
}

int32_t SamCodecActuator::compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create encoder for regular field compression
    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "CIGAR");
    std::shared_ptr<coder> fieldCoder = makeFieldEncoder(fieldIdx, samFieldDefaultCoder(fieldIdx, CoderType::BWT_CM), fieldIo.get(), true);

    fieldSrcLen = 0;

    uint32_t lineCount = lineNum - headEndLine;
    baseLengthBuffer =  MemoryUtil::safeAlloc<uint32_t>(lineCount);
    cigarOpList.resize(lineCount);

    // Process each line and extract the current field
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;
        // Skip header lines (starting with @)
        if (*line == '@') {
            continue;
        }

        uint32_t contentIdx = lineIdx - headEndLine;
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();

        // Parse CIGAR field, remove hard clipping sequence length
        if (fieldLength == 2 && *fieldStart == '*' ) {
            baseLengthBuffer[lineIdx - headEndLine] = 0;
            cigarReadLen[lineIdx] = 0;
            cigarOpList[contentIdx].clear();
        } else {
            uint32_t sequeceLength = cigarSeqConsumed(fieldStart, fieldLength);
            baseLengthBuffer[lineIdx - headEndLine] = sequeceLength;
            /* Reference span for TLEN reconstruction (field 8 is compressed after this field). */
            cigarReadLen[lineIdx] = cigarRefConsumed(fieldStart, fieldLength);
            /* Store the parsed operation list for the SEQ reference rebuild. */
            parseCigarOps(fieldStart, fieldLength, cigarOpList[contentIdx]);
        }

        // Encode the field data
        if (fieldLength > 0) {
            fieldCoder->encode_line(fieldStart, fieldLength);
            fieldSrcLen += fieldLength;
        }
    }

    // Flush the encoder for this field
    fieldCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode field overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    // Set field metadata
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_INFO("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

/* TLEN per the SAM spec: |TLEN| = rightmost mapped base - leftmost mapped base + 1; left fragment positive, right negative. */
bool SamCodecActuator::gatherTLEN(uint32_t lineIdx, TlenInput& in) {
    /* Not paired (FLAG bit 0) or either end unmapped (bits 2/3): TLEN is set to 0. */
    auto flagIt = mappedFlag.find(lineIdx);
    if (flagIt == mappedFlag.end() || !(flagIt->second & 0x1) ||
        (flagIt->second & 0x4) || (flagIt->second & 0x8)) {
        return false;
    }

    /* Reference sequence unavailable or mate on a different reference: TLEN is set to 0. */
    auto chrIt = mappedChr.find(lineIdx);
    if (chrIt == mappedChr.end() || chrIt->second == 0xFFFF) {
        return false;
    }
    auto nextChrIt = nextMappedChr.find(lineIdx);
    if (nextChrIt == nextMappedChr.end() || nextChrIt->second == 0xFFFF) {
        return false;
    }
    if (nextChrIt->second != 0xFFFE && nextChrIt->second != chrIt->second) {
        return false;
    }

    auto posIt = mappedPos.find(lineIdx);
    if (posIt == mappedPos.end()) {
        return false;
    }
    int64_t pos = posIt->second;

    auto pnextIt = nextMappedPos.find(lineIdx);
    if (pnextIt == nextMappedPos.end()) {
        return false;
    }
    int64_t pnext = pnextIt->second;

    /* This read's reference span comes from CIGAR. */
    auto readLenIt = cigarReadLen.find(lineIdx);
    uint32_t refSpan = (readLenIt != cigarReadLen.end()) ? readLenIt->second : 0;

    /* The mate's span is looked up in the full index via (pnext, pos); if not found, 0 is used (a missing span only hurts the compression ratio). */
    uint32_t mateRefSpan = 0;
    auto mateIt = tlenMateIndex.find(std::make_pair(pnext, pos));
    if (mateIt != tlenMateIndex.end()) {
        auto spanIt = cigarReadLen.find(mateIt->second);
        if (spanIt != cigarReadLen.end()) {
            mateRefSpan = spanIt->second;
        }
    }

    in.pos = pos;
    in.pnext = pnext;
    in.refSpan = refSpan;
    in.mateRefSpan = mateRefSpan;
    return true;
}

int32_t SamCodecActuator::computeTLEN(uint32_t lineIdx, bool minusOne) {
    TlenInput in;
    if (!gatherTLEN(lineIdx, in)) {
        return 0;
    }
    /* The rule and the two conventions it covers live in sam_field_rules.h, next to the BAM
       actuator that applies the same one to the same column. */
    return (int32_t)samTemplateLen(in.pos, in.pnext, in.refSpan, in.mateRefSpan, minusOne);
}

int32_t SamCodecActuator::compressTLen(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    Json::Value metaSubs;
    Json::Value metaStreams;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;

    /* Lines that cannot be reconstructed store (relative line number, actual value). */
    std::vector<std::pair<uint32_t, int32_t>> tlenExceptions;

    /* POS/PNEXT are already compressed in earlier fields; build the full mate index here. */
    tlenMateIndex.clear();
    for (const auto& entry : nextMappedPos) {
        auto posIt2 = mappedPos.find(entry.first);
        if (posIt2 == mappedPos.end()) {
            continue;
        }
        tlenMateIndex[std::make_pair(posIt2->second, entry.second)] = entry.first;
    }

    /* First pass: count matches for each TLEN convention (bwa minus 1 / minimap2 plain); the winner is stored in meta. */
    uint64_t convMinus1 = 0, convPlain = 0;
    /* Both conventions' values, kept for the second pass. */
    std::vector<int32_t> tlenMinus1(lineNum, 0);
    std::vector<int32_t> tlenPlain(lineNum, 0);
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        uint32_t contentIdx = lineIdx - headEndLine;
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();
        if (fieldLength > 1) {
            std::string tlenStr = std::string((char*)fieldStart, fieldLength - 1);
            int32_t currentTLEN = (int32_t)std::stoll(tlenStr);
            TlenInput in;
            const bool computable = gatherTLEN(lineIdx, in);
            const int32_t tlenM1 = computable
                ? (int32_t)samTemplateLen(in.pos, in.pnext, in.refSpan, in.mateRefSpan, true) : 0;
            const int32_t tlenPl = computable
                ? (int32_t)samTemplateLen(in.pos, in.pnext, in.refSpan, in.mateRefSpan, false) : 0;
            tlenMinus1[lineIdx] = tlenM1;
            tlenPlain[lineIdx] = tlenPl;
            if (tlenM1 == currentTLEN) convMinus1++;
            if (tlenPl == currentTLEN) convPlain++;
        }
    }
    const bool minusOne = (convMinus1 >= convPlain);
    fieldMeta["tlen_conv"] = minusOne ? 1 : 0;
    LOG_INFO("TLEN convention: minusOne=%d (matches %llu vs %llu)",
             (int)minusOne, (unsigned long long)convMinus1, (unsigned long long)convPlain);

    /* Second pass: collect the lines that cannot be reconstructed under the chosen convention. */
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }

        uint32_t contentIdx = lineIdx - headEndLine;
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();
        fieldSrcLen += fieldLength;

        if (fieldLength > 1) {
            std::string tlenStr = std::string((char*)fieldStart, fieldLength - 1);
            int32_t currentTLEN = (int32_t)std::stoll(tlenStr);
            /* The vote's winner, from the values the first pass already computed. */
            const int32_t computedTLEN = minusOne ? tlenMinus1[lineIdx] : tlenPlain[lineIdx];
            if (computedTLEN != currentTLEN) {
                tlenExceptions.push_back(std::make_pair(contentIdx, currentTLEN));
            }
        }
    }

    if (!tlenExceptions.empty()) {
        std::shared_ptr<coder_io> tlenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "TLEN exceptions");
        std::shared_ptr<coder> tlenCoder =
            CoderFactory::makeProfileEncoder(SAM_TLEN, compressMode(), tlenIo.get(),
                                             CoderType::BWT_CM, false, /* compressLevel */ 0);

        /* Worst case is 5 bytes per varint (32-bit value, 7 payload bits each). */
        const uint32_t cap = (uint32_t)(tlenExceptions.size() * 2 * 5);
        uint8_t* tlenExcBuffer = MemoryUtil::safeAlloc<uint8_t>(cap);
        if (tlenExcBuffer == nullptr) {
            return -1;
        }
        uint32_t tlenExcSrcLen = 0;
        uint32_t prevLine = 0;
        for (uint32_t i = 0; i < tlenExceptions.size(); ++i) {
            const uint32_t lineIdx = tlenExceptions[i].first;
            /* Indices are collected in increasing order -> the delta stays small. */
            tlenExcSrcLen += writeVarint(tlenExcBuffer + tlenExcSrcLen, lineIdx - prevLine);
            prevLine = lineIdx;
            tlenExcSrcLen += writeVarint(tlenExcBuffer + tlenExcSrcLen,
                                           zigzag32(tlenExceptions[i].second));
        }
        tlenCoder->encode_line((uint8_t*)tlenExcBuffer, tlenExcSrcLen);
        tlenCoder->encode_flush();
        MemoryUtil::safeFree(tlenExcBuffer);
        if (tlenIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode TLEN exceptions overflow: output buffer too small");
            return -1;
        }

        metaSubs["srclen"] = tlenExcSrcLen;
        metaSubs["dstlen"] = tlenIo->data_len;
        metaSubs["coder"] = tlenIo->meta;
        metaSubs["sname"] = "tlenexc";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + tlenIo->data_len);
        totalSrcLen += tlenExcSrcLen;
        totalDstLen += tlenIo->data_len;
        /* Which layout the stream is in (see tlenDecodeVarints); set here where a stream exists. */
        fieldMeta["exc"] = (Json::Value::UInt)TLEN_EXC_LAYOUT_VARINT;
    }

    fieldMeta["srclen"] = totalSrcLen;
    fieldMeta["dstlen"] = totalDstLen;
    fieldMeta["streams"] = metaStreams;
    fieldMeta["field"] = fieldIdx;
    fieldMeta["mode"] = "string";
    fieldMeta["optimized"] = true;
    fieldMeta["exceptions"] = (Json::UInt64)tlenExceptions.size();

    LOG_INFO("SAM TLEN field compression with optimization: %u bytes -> %u bytes, %u exceptions compressed, compress ratio = %.2f%%",
        fieldSrcLen, totalDstLen, (uint32_t)tlenExceptions.size(), (double)(totalDstLen * 100)/(double)fieldSrcLen);
    return totalDstLen;
}

int32_t SamCodecActuator::compressBaseWithoutRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    /* Same layout and same coder selection as the reference path, only with no reference to address:
       compressBaseWithRef sees a null file reference, so every record stores its own bases. */
    seqNoReference = true;
    const int32_t ret = compressBaseWithRef(fieldIdx, fieldSrcLen, fieldMeta);
    seqNoReference = false;
    return ret;
}

/* The SEQ match stream split into (run lengths, surviving values), or carried whole when the split would not pay. */
/* Where one field of one content line sits: line start, tab positions, 1-based index; range excludes surrounding tabs. */
/* Defined near the top of the file, before its first user. */

/* Record that `ch` occurred at block position `pos` (see writeSeqExceptionStreams for the emitter). */
static void addSeqException(std::vector<SeqExceptionClass>& classes,
                            std::vector<uint8_t>& present, uint8_t ch, uint32_t pos)
{
    SeqExceptionClass& exc = classes[ch];
    if (exc.count == 0) {
        present.push_back(ch);
        if (ch != (uint8_t)'N') {
            exc.varint.reserve(64);
        }
    }
    exc.add(pos);
}

/* Append one record's payload to the block's match stream, refusing to run past `capacity`. */
static bool appendSeqMatch(uint8_t* matchBuffer, uint32_t& matchLen, uint32_t capacity,
                           const uint8_t* bytes, uint32_t length)
{
    if (matchLen + length > capacity) {
        return false;
    }
    std::memcpy(matchBuffer + matchLen, bytes, length);
    matchLen += length;
    return true;
}

/* Base-length auxiliary stream: reads whose length is not implied by their CIGAR, as (ordinal, length) forward-delta varints. */
int32_t SamCodecActuator::writeSeqBaseLengthStream(Json::Value& metaSubs, Json::Value& metaStreams,
                                                  uint32_t& totalSrcLen, uint32_t& totalDstLen)
{
    if (minBaseLength != maxBaseLength) {
        std::shared_ptr<coder_io> lenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ length");
        /* The reads whose length is not implied by their CIGAR, as (ordinal, length) forward-delta varints. */
        std::vector<uint8_t> lenBuf;
        lenBuf.reserve((unmapedReadLength.size() << 1) + 16);
        uint32_t prevOrdinal = 0;
        for (const auto& entry : unmapedReadLength) {
            appendVarint(lenBuf, entry.first - prevOrdinal);
            prevOrdinal = entry.first;
        }
        for (const auto& entry : unmapedReadLength) {
            appendVarint(lenBuf, entry.second);
        }

        std::shared_ptr<coder> lenCoder =
            CoderFactory::makeProfileEncoder(SAM_SEQ, compressMode(), lenIo.get(),
                                             CoderType::BWT_CM, false, engineCompressLevel());
        lenCoder->encode_line(lenBuf.data(), (uint32_t)lenBuf.size());
        lenCoder->encode_flush();
        if (lenIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode base length stream overflow: output buffer too small");
            return -1;
        }

        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)lenBuf.size();
        metaSubs["dstlen"] = lenIo->data_len;
        metaSubs["coder"] = lenIo->meta;
        metaSubs["sname"] = "baselen";
        /* The delta layout and the two run lengths; without them the stream is absolute pairs. */
        metaSubs["count"] = (Json::Value::UInt)unmapedReadLength.size();
        metaSubs["delta"] = (Json::Value::UInt)1;
        LOG_DEBUG("SEQ baselen stream: entries=%u src=%u dst=%u", (uint32_t)unmapedReadLength.size(),
                  (uint32_t)lenBuf.size(), (uint32_t)lenIo->data_len);
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + lenIo->data_len);
        totalSrcLen += (uint32_t)lenBuf.size();
        totalDstLen += lenIo->data_len;
    }
    return 0;
}

/* Write the SEQ match stream: "m" (run lengths under RLE, else whole) and, under RLE, "mval" (surviving values). */
int32_t SamCodecActuator::writeSeqMatchStreams(const SeqRleSplit& rle, const uint8_t* matchBuffer,
                                               uint32_t matchLen, uint32_t srcLen, coder_io* matchIo,
                                               Json::Value& metaSubs, Json::Value& metaStreams,
                                               uint32_t& totalDstLen)
{
    {
        const uint32_t payLen = rle.useRle ? rle.runLength : matchLen;
        const uint8_t* payBuf = rle.useRle ? rle.run.get() : matchBuffer;

        /* An all-ACGT whole match stream is coded with coder_arith_k (rANS order-1), faster than coder_fc's LZP+BWT+CM. baseNCount is the pre-analysis's count of non-ACGT bases
           (N/n) in this block: when it is non-zero the stream is not all-ACGT, so coder_arith_k cannot pack it and would only rANS its raw bytes, and coder_fc is kept.
           At the extreme-ratio level (-l9) coder_arith_k is never taken: the pre-selection's verdict (coder_fc/coder_bwt_cm) is used as it stands. */
        const bool fastMode = (compressMode() == PBGZ_MODE_FAST);
        CoderType payType = seqMatchCoderFor();
        if (!fastMode && !rle.useRle && baseNCount == 0 && engineCompressLevel() < 9) {
            payType = CoderType::ARITH_K;
        }
        /* A pinned or sampled verdict can still name a segment too short for coder_fc (FC_MIN_LEN bytes or less). */
        if (!fastMode && payType == CoderType::FC &&
            (payLen <= FC_MIN_LEN || payLen >= (uint32_t)FC_MAX_LEN)) {
            payType = CoderType::BWT_CM;
        }
        CoderFactory::applyLevel(matchIo, payType, engineCompressLevel());
        std::shared_ptr<coder> payCoder = CoderFactory::makeEncoder(payType, matchIo);
        payCoder->encode_line(const_cast<uint8_t*>(payBuf), payLen);
        payCoder->encode_flush();
        LOG_DEBUG("SEQ match stream: matchLen=%u nonZero=%u (%.2f%%) rle=%d coder=%s src=%u dst=%u",
                  matchLen, rle.nonZero,
                  (double)rle.nonZero * 100.0 / (double)(matchLen ? matchLen : 1),
                  (int)rle.useRle, coderTypeToMagic(payType), payLen,
                  (uint32_t)matchIo->data_len);
    }
    if (matchIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode base match stream overflow: output buffer too small");
        return -1;
    }
    // First sub-stream: run lengths (RLE) / whole match stream (legacy layout)
    metaSubs.clear();
    metaSubs["srclen"] = rle.useRle ? (Json::Value::UInt)rle.runLength : (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] = matchIo->data_len;
    metaSubs["coder"] = matchIo->meta;
    metaSubs["sname"] = "m";
    if (rle.useRle) {
        /* orgrawlen: the expanded original length (base count); rlelen: this sub-stream's decode output length. */
        metaSubs["rle"] = (Json::Value::UInt)1;
        metaSubs["orgrawlen"] = (Json::Value::UInt)srcLen;
        metaSubs["rlelen"] = (Json::Value::UInt)rle.runLength;
    }
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + matchIo->data_len);
    totalDstLen += matchIo->data_len;

    /* Second sub-stream: the surviving non-zero values (RLE only). */
    if (rle.useRle && rle.valLength > 0) {
        std::shared_ptr<coder_io> valIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match val");
        /* Decoding dispatches on the recorded magic, so the mode picks the coder here too. */
        {
            std::shared_ptr<coder> valCoder =
                CoderFactory::makeProfileEncoder(SAM_SEQ, compressMode(), valIo.get(),
                                                 CoderType::BWT_CM, false, engineCompressLevel());
            valCoder->encode_line(rle.val.get(), rle.valLength);
            valCoder->encode_flush();
        }
        if (valIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ match value stream overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)rle.valLength;
        metaSubs["dstlen"] = valIo->data_len;
        metaSubs["coder"] = valIo->meta;
        metaSubs["sname"] = "mval";
        LOG_DEBUG("SEQ mval stream: src=%u dst=%u", rle.valLength, (uint32_t)valIo->data_len);
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + valIo->data_len);
        totalDstLen += valIo->data_len;
    }
    return 0;
}

bool SamCodecActuator::referenceUsableForSeq()
{
    if (pRefeGene == nullptr) {
        return false;
    }
    if (!refeUsableChecked) {
        refeUsableChecked = true;
        refeUsableForSeq = referenceLinesUpWithHeader(pRefeGene);
        if (!refeUsableForSeq) {
            warnReferenceNotUsableForSeq();
        }
    }
    return refeUsableForSeq;
}

/* The 2-bit code a SEQ character stands for, or -1 (order matches actgSquash / atcg4). */
static int atcgIndexOf(uint8_t ch)
{
    switch (ch) {
        case 'A': return 0;
        case 'C': return 1;
        case 'T': return 2;
        case 'G': return 3;
        default: return -1;
    }
}

bool SamCodecActuator::buildBlockConsensus(uint32_t fieldIdx, BlockConsensus& out)
{
    /* Bounds on building a consensus: window size cap (counters/squash memory) and a coverage floor. */
    static const int64_t kMaxWindowBases = 2 << 20;   /* 2 Mb: ~8 MB of counters, ~0.5 MB of squash */
    static const int64_t kMinCoverageShift = 3;       /* at least 1/8 of the window covered */

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    const uint32_t lineNum = (uint32_t)npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    /* Pass 1: the window, over every record that could be walked against a reference at all. */
    int64_t winStart = INT64_MAX;
    int64_t winEnd = 0;
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        const uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        const uint16_t chrId = mappedChr.find(lineIdx) == mappedChr.end() ? 0xFFFF : mappedChr[lineIdx];
        if (chrId == 0xFFFF || chrId == 0xFFFE) {
            continue;
        }
        const uint16_t flag = mappedFlag.find(lineIdx) == mappedFlag.end() ? 4 : mappedFlag[lineIdx];
        if (flag & 0x04) {
            continue;
        }
        const int64_t chrStart = SamInfo::getInstance().getPositionByIndex(chrId);
        if (chrStart == -1) {
            continue;
        }
        const auto posIt = mappedPos.find(lineIdx);
        if (posIt == mappedPos.end()) {
            continue;
        }
        const int64_t refStart = chrStart + (int64_t)posIt->second - 1;
        const auto crlIt = cigarReadLen.find(lineIdx);
        const uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
        if (refStart < 0 || refConsumed == 0) {
            continue;
        }
        if (refStart < winStart) {
            winStart = refStart;
        }
        if (refStart + (int64_t)refConsumed > winEnd) {
            winEnd = refStart + (int64_t)refConsumed;
        }
    }
    if (winStart == INT64_MAX) {
        return false;
    }
    const int64_t windowBases = winEnd - winStart;
    if (windowBases <= 0 || windowBases > kMaxWindowBases) {
        return false;
    }

    /* Pass 2: one saturating counter per base. */
    std::vector<uint8_t> count[4];
    for (int k = 0; k < 4; ++k) {
        count[k].assign((size_t)windowBases, 0);
    }
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        const uint32_t lineStart = (lineIdx == 0) ? 0 : (uint32_t)npos[lineIdx - 1] + 1;
        const uint32_t lineEnd = (uint32_t)npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        const uint32_t contentIdx = lineIdx - headEndLine;
        const uint16_t chrId = mappedChr.find(lineIdx) == mappedChr.end() ? 0xFFFF : mappedChr[lineIdx];
        const uint16_t flag = mappedFlag.find(lineIdx) == mappedFlag.end() ? 4 : mappedFlag[lineIdx];
        const auto posIt = mappedPos.find(lineIdx);
        if (chrId == 0xFFFF || chrId == 0xFFFE || (flag & 0x04) || posIt == mappedPos.end()) {
            continue;
        }
        const int64_t chrStart = SamInfo::getInstance().getPositionByIndex(chrId);
        if (chrStart == -1) {
            continue;
        }
        const int64_t refStart = chrStart + (int64_t)posIt->second - 1;
        const auto crlIt = cigarReadLen.find(lineIdx);
        if (crlIt == cigarReadLen.end() || crlIt->second == 0) {
            continue;
        }
        const SamFieldRange seq = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        const uint8_t* seqStart = seq.start;
        const uint32_t seqLength = seq.lengthWithoutTab();
        if (seqLength == 0 || contentIdx >= cigarOpList.size() || cigarOpList[contentIdx].empty()) {
            continue;
        }
        /* The same CIGAR walk as the coding side, so the vote lands on the bases the walk will compare. */
        const std::vector<CigarOp>& ops = cigarOpList[contentIdx];
        uint32_t readPos = 0;
        int64_t refPos = refStart - winStart;
        for (size_t oi = 0; oi < ops.size() && readPos < seqLength; ++oi) {
            const CigarOp& op = ops[oi];
            if (op.op == 'M' || op.op == '=' || op.op == 'X') {
                for (uint32_t i = 0; i < op.len && readPos < seqLength; ++i, ++readPos, ++refPos) {
                    const int k = atcgIndexOf(seqStart[readPos]);
                    if (k >= 0 && refPos >= 0 && refPos < windowBases) {
                        uint8_t& c = count[k][(size_t)refPos];
                        if (c < 255) {
                            ++c;
                        }
                    }
                }
            } else if (op.op == 'I' || op.op == 'S') {
                readPos += op.len;
            } else if (op.op == 'D' || op.op == 'N') {
                refPos += op.len;
            }
        }
    }

    /* Majority rule; an uncovered position stays at code 0. */
    const uint64_t squashLen = (uint64_t)(windowBases >> 2) + !!(windowBases & 0x3) + 1;
    std::vector<uint8_t> squash((size_t)squashLen, 0);
    int64_t covered = 0;
    for (int64_t i = 0; i < windowBases; ++i) {
        int best = 4;
        uint8_t bestCount = 0;
        for (int k = 0; k < 4; ++k) {
            if (count[k][(size_t)i] > bestCount) {
                bestCount = count[k][(size_t)i];
                best = k;
            }
        }
        if (best >= 4) {
            continue;
        }
        ++covered;
        const uint64_t bit = (uint64_t)i;
        squash[(size_t)(bit >> 2)] |= (uint8_t)(best << (6 - 2 * (bit & 0x3)));
    }
    if ((covered << kMinCoverageShift) < windowBases) {
        return false;
    }

    out.windowStart = winStart;
    out.windowBases = windowBases;
    out.squash.swap(squash);
    return true;
}

int32_t SamCodecActuator::compressBaseSelfRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta)
{
    BlockConsensus cons;
    const auto consT0 = pbgzprof::nowOrZero();
    if (!buildBlockConsensus(fieldIdx, cons)) {
        /* Nothing to build a reference from: the plain path is what is left. */
        return compressBaseWithoutRef(fieldIdx, fieldSrcLen, fieldMeta);
    }
    pbgzprof::addSince(pbgzprof::SEQ_MATCH_CONSENSUS, consT0);
    LOG_INFO("SEQ block %lld coded against a consensus built from its own reads: window %lld + %lld bases",
             (long long)inBlockPtr->getBlockId(), (long long)cons.windowStart,
             (long long)cons.windowBases);
    seqSelfRef = &cons;
    const int32_t ret = compressBaseWithRef(fieldIdx, fieldSrcLen, fieldMeta);
    seqSelfRef = nullptr;
    return ret;
}

int32_t SamCodecActuator::compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    /* The reference this block codes against: the file's reference, a block-built consensus
       (seqSelfRef), or none at all (compressBaseWithoutRef) - with none every record stores its own
       bases, and the "streams" layout still applies unchanged. */
    Reference* fileRef = seqNoReference ? nullptr : pRefeGene;

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    /* SEQ characters the 2-bit path cannot carry, one accumulator per character (see SeqExceptionClass); positions are delta-coded. */
    std::vector<SeqExceptionClass> excClasses(256);
    std::vector<uint8_t> excPresent;
    excClasses[(uint8_t)'N'].keepForms = true;
    excClasses[(uint8_t)'N'].abs.reserve((size_t)baseNCount + 16);
    excClasses[(uint8_t)'N'].varint.reserve((size_t)baseNCount + 16);

    // Initialize mapping buffers similar to FastqActuator
    const uint32_t baseMaxLength = inBlockPtr->getMaxLineLen() + 4;
    const uint32_t lsquash = (baseMaxLength >> 2) + !!(baseMaxLength & 0x3);

    uint32_t baseMappedLength = (baseMaxLength << 1);

    std::unique_ptr<uint8_t[]> basePairBuffer = std::make_unique<uint8_t[]>(baseMaxLength);
    baseSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(lsquash);
    std::unique_ptr<uint8_t[]> baseMappedBuffer = std::make_unique<uint8_t[]>(baseMappedLength);
    /* The N positions are collected into nposVarint below, which needs no upfront
       allocation of the (up to 4x larger) absolute-offset array. */
    /* Scratch for the reference 2-bit sequence in the CIGAR M/=/X rebuild; getStretch2Bits1Char can overshoot by 3 bytes. */
    std::unique_ptr<uint8_t[]> ref2bitBuf = std::make_unique<uint8_t[]>(baseMaxLength + 8);
    /* Empty op list for reads with a `*` CIGAR (no operations). */
    static const std::vector<CigarOp> emptyCigarOps;

    // Create metadata structure
    Json::Value metaSubs;
    Json::Value metaStreams;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;

    // Second pass: compress with reference
    std::shared_ptr<coder_io> matchIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match");
    /* The SEQ match stream (per-base 2-bit: M/=/X XOR the reference, I/S raw). Accumulated whole-block, then encoded once. */
    std::unique_ptr<uint8_t[]> matchBuffer = std::make_unique<uint8_t[]>(inBlockPtr->getDataLen() + 1);
    uint32_t matchLen = 0;
    int64_t srcLen = 0;
    uint32_t totalBaseLength = 0;

    const auto seqBuildT0 = pbgzprof::nowOrZero();
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;
        // Skip header lines (starting with @)
        if (*line == '@') {
            continue;
        }

        uint32_t contentIdx = lineIdx - headEndLine;
        const SamFieldRange seq = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* seqStart = seq.start;
        uint32_t seqLength = seq.lengthWithoutTab();

        if (seqLength == 0)  {
            continue;
        }

        if (baseLengthBuffer && baseLengthBuffer[contentIdx] == 0) {
            baseLengthBuffer[contentIdx] = seqLength;
            unmapedReadLength.push_back(std::make_pair(contentIdx, seqLength));
        }

        if (baseLengthBuffer && baseLengthBuffer[contentIdx] != seqLength) {
            LOG_WARNING("Warnning: sequece length(%u) not match cigar(%u)", seqLength, baseLengthBuffer[contentIdx]);
            baseLengthBuffer[contentIdx] = seqLength;
        }

        // Get mapping information from SAM fields
        uint16_t chrId = mappedChr.find(lineIdx) == mappedChr.end() ? 0xFFFF : mappedChr[lineIdx];
        uint64_t startPos = mappedPos.find(lineIdx) == mappedPos.end() ? 0 : mappedPos[lineIdx];
        // Extract FLAG field to determine strand
        uint16_t flag = mappedFlag.find(lineIdx) == mappedFlag.end() ? 4 : mappedFlag[lineIdx];

        /* CIGAR-segment reference rebuild: M/=/X are 2-bit XOR against the reference, I/S stored directly, D/N only advance. */
        static const std::vector<CigarOp> emptyCigarOps;
        /* The reference span this record's CIGAR consumes, from the parse the caller did. */
        const auto crlIt = cigarReadLen.find(lineIdx);
        const uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
        const std::vector<CigarOp>& ops =
            (contentIdx < cigarOpList.size()) ? cigarOpList[contentIdx] : emptyCigarOps;

        /* Which reference this record is walked against: the block consensus when in play, otherwise the file's. */
        const SeqRecordPayload payload =
            (seqSelfRef != nullptr)
                ? buildSeqSquashRecordPayload(chrId, flag, startPos, refConsumed, ops, seqStart,
                                              seqLength, seqSelfRef->squash.data(),
                                              seqSelfRef->windowStart, seqSelfRef->windowBases,
                                              baseMappedBuffer.get(), ref2bitBuf.get())
                : buildSeqRecordPayload(chrId, flag, startPos, refConsumed, ops, seqStart, seqLength,
                                        fileRef, baseMappedBuffer.get(), ref2bitBuf.get());
        const uint32_t outLen = payload.length;
        const bool rawBases = payload.rawBases;
        /* A consensus is not a file reference and carries no match statistics to update. */
        if (payload.usedReference && seqSelfRef == nullptr && fileRef != nullptr) {
            fileRef->updateMatchedGene((uint64_t)payload.refPos,
                (crlIt != cigarReadLen.end()) ? crlIt->second : seqLength);
        }
        /* Characters the payload cannot carry are recorded, one list per character (see SeqExceptionClass). */
        for (uint32_t n = 0; n < seqLength; n++) {
            const uint8_t ch = seqStart[n];
            const bool carried = (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') ||
                                 (rawBases && ch == 'N');
            if (!carried) {
                addSeqException(excClasses, excPresent, ch, totalBaseLength + n);
            }
        }
        totalBaseLength += seqLength;

        // Encode the mapped data
        if (outLen > 0) {
            /* Accumulate instead of encoding line-by-line: the coder is picked by
               the total length after the loop (coder_fc is block-only). */
            if (!appendSeqMatch(matchBuffer.get(), matchLen, inBlockPtr->getDataLen(),
                                payload.bytes, payload.length)) {
                LOG_ERROR("SEQ match buffer overflow in block %llu: %u + %u > %u",
                          inBlockPtr->getBlockId(), matchLen, seqLength, inBlockPtr->getDataLen());
                return -1;
            }
            srcLen += seqLength;
        }
    }
    /* RLE preprocessing: split the match stream into "m" (varint zero-run lengths) and "mval" (non-zero values), applied only when at least 98% zeros. */
    const auto seqSplitT0 = pbgzprof::nowOrZero();
    SeqRleSplit rle = splitSeqMatchStream(matchBuffer.get(), matchLen);
    pbgzprof::addSince(pbgzprof::SEQ_MATCH_BUILD, seqBuildT0);
    pbgzprof::addSince(pbgzprof::SEQ_MATCH_SPLIT, seqSplitT0);

    /* Sub-stream "m": run lengths under RLE, otherwise the whole match stream; the coder is the per-file verdict (seqMatchCoderFor). */
    /* Sub-stream "m": run-length segment under RLE, otherwise the whole match stream. */
    {
    const auto seqEncT0 = pbgzprof::nowOrZero();
    if (writeSeqMatchStreams(rle, matchBuffer.get(), matchLen, srcLen, matchIo.get(), metaSubs,
                             metaStreams, totalDstLen) != 0) {
        return -1;
    }
    pbgzprof::addSince(pbgzprof::SEQ_MATCH_ENCODE, seqEncT0);
    totalSrcLen += srcLen;
    }


    /* SEQ exception sub-streams: one per character present, in ascending character order (see writeSeqExceptionStreams). */
    if (writeSeqExceptionStreams(excClasses, excPresent, metaStreams, totalSrcLen, totalDstLen) != 0) {
        return -1;
    }

    if (writeSeqBaseLengthStream(metaSubs, metaStreams, totalSrcLen, totalDstLen) != 0) {
        return -1;
    }

    /* The block's own consensus, appended last: its squash plus the window it covers. */
    if (seqSelfRef != nullptr) {
        std::shared_ptr<coder_io> selfIo = makeCoderIo(outBlockPtr->getCurrent(),
                                                       outBlockPtr->getRemain(), "SEQ selfref");
        std::shared_ptr<coder> selfCoder =
            CoderFactory::makeProfileEncoder(SAM_SEQ, compressMode(), selfIo.get(),
                                             CoderType::BWT_CM, false, engineCompressLevel());
        selfCoder->encode_line(seqSelfRef->squash.data(), (uint32_t)seqSelfRef->squash.size());
        selfCoder->encode_flush();
        if (selfIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ self reference overflow: output buffer too small");
            return -1;
        }
        Json::Value selfMeta;
        selfMeta["srclen"] = (Json::Value::UInt)seqSelfRef->squash.size();
        selfMeta["dstlen"] = selfIo->data_len;
        selfMeta["coder"] = selfIo->meta;
        selfMeta["sname"] = "selfref";
        metaStreams.append(selfMeta);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + selfIo->data_len);
        totalDstLen += selfIo->data_len;
        fieldMeta["selfref"] = (Json::Value::UInt)1;
        fieldMeta["selfref_wstart"] = (Json::Value::Int64)seqSelfRef->windowStart;
        fieldMeta["selfref_wbases"] = (Json::Value::Int64)seqSelfRef->windowBases;
    }

    // Set metadata
    /* Count of the 'N' stream, for older archives. */
    fieldMeta["ncount"] = (Json::Value::UInt)excClasses[(uint8_t)'N'].count;
    fieldMeta["litbases"] = (Json::Value::UInt)1;
    /* No reference (file's or consensus): the decoder must not try to rebuild any record against one. */
    if (fileRef == nullptr && seqSelfRef == nullptr) {
        fieldMeta["noref"] = (Json::Value::UInt)1;
    }

    /* How to read a payload byte above the 2-bit range: as one of the record's own characters. */
    fieldMeta["minlen"] = minBaseLength;
    fieldMeta["maxlen"] = maxBaseLength;
    fieldMeta["totalsrclen"] = totalSrcLen;
    fieldMeta["totaldstlen"] = totalDstLen;
    fieldMeta["streams"] = metaStreams;
    fieldMeta["field"] = fieldIdx;
    /* fieldSrcLen is corrected to the match stream's source length only; meta's totalsrclen keeps the auxiliary streams. */
    fieldSrcLen = (uint32_t)srcLen;

    LOG_INFO("SAM base field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        totalSrcLen, totalDstLen, (double)(totalDstLen * 100)/(double)totalSrcLen);

    return totalDstLen;
}

int32_t SamCodecActuator::writeSeqExceptionStream(SeqExceptionClass& exc, uint8_t ch,
                                                   const char* sname, uint32_t srcLen,
                                                   const uint8_t* src, uint32_t runCount,
                                                   Json::Value& streamMeta, uint32_t& totalSrcLen,
                                                   uint32_t& totalDstLen)
{
    std::shared_ptr<coder_io> excIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ npos");
    std::shared_ptr<coder> subCoder =
        CoderFactory::makeProfileEncoder(SAM_SEQ, compressMode(), excIo.get(),
                                         CoderType::BWT_CM, false, engineCompressLevel());
    subCoder->encode_line(src, srcLen);
    subCoder->encode_flush();
    if (excIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode SEQ exception positions overflow: output buffer too small");
        return -1;
    }
    Json::Value metaSubs;
    metaSubs["srclen"] = srcLen;
    metaSubs["dstlen"] = excIo->data_len;
    metaSubs["coder"] = excIo->meta;
    metaSubs["sname"] = sname;
    metaSubs["ch"] = (Json::Value::UInt)ch;
    metaSubs["count"] = exc.count;
    LOG_DEBUG("SEQ exc ch=%u '%c' form=%s src=%u dst=%u count=%u runs=%u", (unsigned)ch,
              (char)ch, sname, srcLen, (uint32_t)excIo->data_len, exc.count, runCount);
    if (runCount > 0) {
        /* The run form: how many runs the two varint sections below hold. */
        metaSubs["runs"] = (Json::Value::UInt)runCount;
    }
    streamMeta.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + excIo->data_len);
    totalSrcLen += srcLen;
    totalDstLen += excIo->data_len;
    return 0;
}

int32_t SamCodecActuator::writeSeqExceptionStreams(std::vector<SeqExceptionClass>& excClasses,
                                                   std::vector<uint8_t>& excPresent,
                                                   Json::Value& streamMeta, uint32_t& totalSrcLen,
                                                   uint32_t& totalDstLen)
{
    std::sort(excPresent.begin(), excPresent.end());
    for (size_t ei = 0; ei < excPresent.size(); ++ei) {
        const uint8_t ch = excPresent[ei];
        SeqExceptionClass& exc = excClasses[ch];

        /* The 'N' class's form is a file-level verdict from the codec pre-selection; the absolute form is the fallback. */
        int32_t form = 1;
        if (ch == (uint8_t)'N') {
            PreprocessInfo* preInfo = preprocessInfoMut();
            const int32_t decided =
                (preInfo != nullptr) ? preInfo->nposForm.load(std::memory_order_relaxed) : -1;
            form = (decided >= 0) ? decided : 0;
        }

        std::vector<uint8_t> runBuf;
        uint32_t runCount = 0;
        const char* sname = kSeqExcDeltaName;
        const uint8_t* src = exc.varint.data();
        uint32_t srcLen = (uint32_t)exc.varint.size();
        if (form == 0) {
            sname = kSeqExcAbsName;
            src = (const uint8_t*)exc.abs.data();
            srcLen = (uint32_t)(exc.abs.size() * sizeof(uint32_t));
        } else if (form == 2) {
            exc.closeRun();
            runBuf.reserve((exc.runGaps.size() << 1) + 16);
            buildRunForm(exc, runBuf);
            runCount = (uint32_t)exc.runLens.size();
            sname = kSeqExcRunName;
            src = runBuf.data();
            srcLen = (uint32_t)runBuf.size();
        }
        if (writeSeqExceptionStream(exc, ch, sname, srcLen, src, runCount, streamMeta,
                                    totalSrcLen, totalDstLen) != 0) {
            return -1;
        }
    }
    return 0;
}

/* Encode the QUAL column record by record through the chosen coder; missing-quality records are expanded to seqLen '*' first. */
uint32_t SamCodecActuator::encodeQualRecords(qual_record_encoder* encoder, uint32_t lineNum,
                                                   uint32_t fieldIdx, bool needStrand)
{
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint8_t* buffer = inBlockPtr->getBuffer();
    uint32_t streamSrcLen = 0;
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;

        uint8_t* line = buffer + lineStart;

        // Skip header lines (starting with @)
        if (*line == '@') {
            continue;
        }

        // Extract QUAL field (field 11)
        uint32_t contentIdx = lineIdx - headEndLine;
        if (fieldIdx > contentPos[contentIdx].size()) {
            continue;
        }

        const SamFieldRange qual = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* qualStart = qual.start;
        uint32_t qualLength = qual.lengthWithoutTab();

        if (qualLength == 0) {
            continue;
        }

        // Get SEQ field for quality compression context (field 9)
        uint8_t* seqStart = nullptr;
        uint32_t seqLen = 0;
        if (fieldIdx >= 1 && contentPos[contentIdx].size() >= fieldIdx) {
            seqStart = line + contentPos[contentIdx][fieldIdx - 2] + 1;
            seqLen = contentPos[contentIdx][fieldIdx - 1] - contentPos[contentIdx][fieldIdx - 2] - 1;
        }

        /* Missing quality (single '*') is expanded to seqLen '*' so the record length matches the decoder's fetch. */
        if (qualLength == 1 && qualStart[0] == '*' && seqLen > 0) {
            qualMissingBuf.assign(seqLen, '*');
            qualStart = qualMissingBuf.data();
            qualLength = seqLen;
        }

        /* Strand direction from FLAG bit 0x10, parsed by hand (not NUL-terminated, hot path). */
        bool rev = false;
        if (needStrand && contentPos[contentIdx].size() > 1) {
            uint32_t flagBeg = contentPos[contentIdx][0] + 1;
            uint32_t flagEnd = contentPos[contentIdx][1];
            long flagVal = 0;
            for (uint32_t fp = flagBeg; fp < flagEnd; ++fp) {
                uint8_t ch = line[fp];
                if (ch < '0' || ch > '9') break;
                flagVal = flagVal * 10 + (ch - '0');
            }
            rev = (flagVal & 16) != 0;
        }
        encoder->encode_record(qualStart, qualLength, seqStart, seqLen, rev);
        streamSrcLen += qualLength;
    }
    return streamSrcLen;
}

/* Write the QUAL frequency table as an auxiliary stream: alphabet in frequency-descending order, two uint16 per entry, BWT_CM. */
int32_t SamCodecActuator::writeQualFreqStream(Json::Value& subMeta, Json::Value& streamMeta,
                                               uint32_t& totalSrcLength, uint32_t& totalDstLength)
{
    std::shared_ptr<coder_io> qualityFreqIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL freq table");
    std::shared_ptr<coder> qualityFreqCoder =
        CoderFactory::makeProfileEncoder(SAM_QUAL, compressMode(), qualityFreqIo.get(),
                                         CoderType::BWT_CM, false, engineCompressLevel());
    std::shared_ptr<uint16_t[]> qualityFreqArray(new uint16_t[qualFreqTable.size() << 1]);
    for (uint32_t i = 0; i < qualFreqTable.size(); ++i) {
        int idx = i << 1;
        qualityFreqArray[idx] = qualFreqTable[i].first;
        qualityFreqArray[idx + 1] = qualFreqTable[i].second;
    }

    uint32_t freqSrcLen = (qualFreqTable.size() << 1) * sizeof(uint16_t);
    qualityFreqCoder->encode_line((uint8_t*)qualityFreqArray.get(), freqSrcLen);
    qualityFreqCoder->encode_flush();
    if (qualityFreqIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode quality freq overflow: output buffer too small");
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityFreqIo->data_len);

    subMeta.clear();
    subMeta["srclen"] = freqSrcLen;
    subMeta["dstlen"] = qualityFreqIo->data_len;
    subMeta["coder"] = qualityFreqIo->meta;
    subMeta["streamname"] = "qualityfreq";
    streamMeta.append(subMeta);

    /* The frequency-table auxiliary stream still counts toward meta's totalsrclen (original accounting preserved). */
    totalSrcLength += freqSrcLen;
    totalDstLength += qualityFreqIo->data_len;
    return 0;
}

int32_t SamCodecActuator::compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::chrono::steady_clock::time_point qualT0, qualT1, qualT2;
    qualT0 = pbgzprof::nowOrZero();
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();

    // Create quality encoder similar to FastqActuator
    std::shared_ptr<coder_io> qualityIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL");

    /* QUAL may use coder_qual (with SEQ context) or coder_fcv2 (context mixing); the pre-selection's
       trial decides. In fast mode the trial picks the fast default (coder_rans). */
    int64_t qualPriorAddress = -1;
    CoderType pickedQualCoder = samFieldDefaultCoder(SAM_QUAL, compressMode(), CoderType::QUAL);
    const PreprocessInfo* qualPreInfo =
        (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
    if (qualPreInfo != nullptr) {
        const CoderType picked =
            qualPreInfo->coderFor(SAM_QUAL, pickedQualCoder);
        if (CoderFactory::eligibleProfile(picked, compressMode())) {
            pickedQualCoder = picked;
        }
    }
    /* Turning the verdict into an encoder is CoderFactory::makeQualEncoder's job; the prior's absolute address goes into the meta. */
    AuxPayloadPtr qualPriorBlob =
        (pbgzEngine != nullptr) ? pbgzEngine->getQualPrior(0) : AuxPayloadPtr();
    bool qualPriorLoaded = false;
    QualCoderArgs qualCoderArgs;
    qualCoderArgs.freqTable = qualFreqTable.empty() ? nullptr : &qualFreqTable;
    const FieldCodecSelection* qualSel =
        (qualPreInfo != nullptr) ? qualPreInfo->getField(SAM_QUAL) : nullptr;
    qualCoderArgs.fcv2Params =
        (qualSel != nullptr && qualSel->selectedCoder == CoderType::FCV2) ? &qualSel->fcv2Params
                                                                         : nullptr;
    qualCoderArgs.priorBlob = qualPriorBlob.get();
    qualCoderArgs.priorLoaded = &qualPriorLoaded;
    qualPriorAddress = (pbgzEngine != nullptr) ? pbgzEngine->getQualPriorAddress() : -1;

    qualT1 = pbgzprof::nowOrZero();
    std::shared_ptr<qual_record_encoder> qualEncoder =
        CoderFactory::makeQualEncoder(pickedQualCoder, qualityIo.get(), qualCoderArgs,
                                      engineCompressLevel());
    pbgzprof::addSince(pbgzprof::QUAL_PREP, qualT0);

    /* A block whose quality values left the prior's alphabet cannot use that prior; it is dropped for that block. */
    if (pickedQualCoder == CoderType::FCV2 && qualPriorBlob.get() != nullptr && qualEncoder != nullptr) {
        bool covered = true;
        for (size_t i = 0; i < qualFreqTable.size(); ++i) {
            if (!qualEncoder->coversByte((uint8_t)(qualFreqTable[i].first + '!'))) {
                covered = false;
                break;
            }
        }
        if (!covered) {
            LOG_INFO("QUAL prior does not cover block %lld's quality values, coding it without the prior.",
                     (long long)inBlockPtr->getBlockId());
            qualCoderArgs.priorBlob = nullptr;
            std::shared_ptr<coder_io> plainIo =
                makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL");
            qualEncoder = CoderFactory::makeQualEncoder(pickedQualCoder, plainIo.get(), qualCoderArgs,
                                                        engineCompressLevel());
            qualityIo = plainIo;
            qualPriorLoaded = false;
        }
    }
    if (!qualPriorLoaded) {
        qualPriorAddress = -1;
    }
    /* The measurement above covers the constructor's cost on every path. */
    pbgzprof::addSince(pbgzprof::QUAL_CTOR, qualT1);
    qualT2 = pbgzprof::nowOrZero();

    /* Only coder_fcv2 restores the sequencing cycle from the strand direction. */
    const bool needStrand = (pickedQualCoder == CoderType::FCV2);

    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;
    Json::Value streamMeta;

    // Encode quality data (see encodeQualRecords); the encoder is committed below.
    const uint32_t streamSrcLen = encodeQualRecords(qualEncoder.get(), lineNum, fieldIdx, needStrand);

    qualEncoder->flush();
    pbgzprof::addSince(pbgzprof::QUAL_CODEC, qualT2);
    if (qualityIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode quality overflow: output buffer too small");
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityIo->data_len);

    Json::Value subMeta;
    subMeta["srclen"] = streamSrcLen;
    subMeta["dstlen"] = qualityIo->data_len;
    subMeta["coder"] = qualityIo->meta;
    if (qualPriorAddress >= 0) {
        /* Absolute file offset of the prior block container header; the decoding side uses it to retrieve the same snapshot. */
        subMeta["prior"] = (Json::Value::Int64)qualPriorAddress;
    }
    streamMeta.append(subMeta);

    totalSrcLength += streamSrcLen;
    totalDstLength += qualityIo->data_len;

    // Encode quality frequency table
    /* The quality frequency table travels as an auxiliary stream: see writeQualFreqStream. */
    if (writeQualFreqStream(subMeta, streamMeta, totalSrcLength, totalDstLength) != 0) {
        return -1;
    }

    // Set field metadata
    fieldMeta["totalsrclen"] = totalSrcLength;
    fieldMeta["totaldstlen"] = totalDstLength;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["field"] = fieldIdx;

    /* fieldSrcLen is corrected to the quality text length only; meta's totalsrclen keeps the frequency-table stream. */
    fieldSrcLen = streamSrcLen;

    /* Verbose mode: print per block which QUAL coder was used and this block's ratio (single fprintf per line for thread safety). */
    if (pbgzEngine != nullptr && pbgzEngine->getParameter().verbose) {
        const double qualRatio = (totalSrcLength == 0)
            ? 0.0
            : (double)(totalDstLength * 100) / (double)totalSrcLength;
        fprintf(stderr,
                "[block %lu] QUAL -> %-12s  %u -> %u (%.2f%%)\n",
                (unsigned long)inBlockPtr->getBlockId(),
                coderTypeToMagic(pickedQualCoder),
                totalSrcLength,
                totalDstLength,
                qualRatio);
    }

    LOG_INFO("SAM quality field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        totalSrcLength, totalDstLength, (double)(totalDstLength * 100)/(double)totalSrcLength);

    return totalDstLength;
}

void SamCodecActuator::initMetaInfo() {
    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);
    if (meta.isMember("header")) {
        headEndLine = meta["header"]["lines"].asInt64();
    }
    if (meta.isMember("sam")) {
        samLine =  meta["sam"]["lines"].asUInt();
    }
    return;
}

int32_t SamCodecActuator::decompress() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter, inBlockPtr or outBlockPtr is nullptr for SAM decompression");
        return -1;
    }
    // Reset read offset before decompression
    readOffset = 0;
    optionCacheEmpty = true;
    optionRecText.clear();
    optionRecOff.clear();
    optionRecPresent.clear();
    // Parse meta information
    initMetaInfo();

    /* -m fast blocks carry structured "bam" meta; hand them to the column decoder. */
    if (meta.isMember("bam")) {
        BamCodecActuator bamActuator(inBlockPtr, outBlockPtr, pbgzEngine, pRefeGene);
        int32_t ret = bamActuator.decompress();
        return ret;
    }

    /* Pre-allocate block_size*2 from the file header (fallback getBlockSize()), the primary out-of-bounds defense. */
    size_t bs = pbgzEngine->getFileBlockSize();
    if (bs == 0) {
        bs = ConfigManager::getInstance().getBlockSizeByCompressLevel(pbgzEngine->getParameter().compressLevel);
    }
    /* Pre-allocate as the larger of block_size*2 and this block's actual data length. */
    size_t outCapacity = bs * 2;
    if ((size_t)inBlockPtr->getDataLen() > outCapacity) {
        outCapacity = (size_t)inBlockPtr->getDataLen();
    }
    /* The block's textlen is the exact output bound. */
    if (meta.isMember("sam") && meta["sam"].isMember("textlen")) {
        const uint64_t textLen = meta["sam"]["textlen"].asUInt64();
        if ((uint64_t)outCapacity < textLen) {
            outCapacity = (size_t)textLen;
        }
    }
    if (outBlockPtr->ensureCapacity(outCapacity) != 0) {
        LOG_ERROR("preallocate output buffer failed, need=%zu", outCapacity);
        return -1;
    }

    // Set block information
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());

    RoughIOBlock* targeBlock = outBlockPtr;

    const PbgzParameter& parameter = pbgzEngine->getParameter();
    if (!parameter.refeGenePos.empty()) {
        do {
            // Parse the p parameter, format like chr1:100-200, convert the part before : to chrID, the part before and after - for the reference gene start position
            size_t colonPos = parameter.refeGenePos.find(':');
            if (colonPos == std::string::npos) {   // Format with only chromosome name
                refPosChrIndex = SamInfo::getInstance().getChrNameIndex(parameter.refeGenePos);
                break;
            }

            refPosChrIndex = SamInfo::getInstance().getChrNameIndex(parameter.refeGenePos.substr(0, colonPos));
            size_t dashPos = parameter.refeGenePos.find('-');
            if (dashPos == std::string::npos) {
                break;
            }

            refPosBegin = std::stoi(parameter.refeGenePos.substr(colonPos + 1, dashPos - colonPos - 1));
            refPosEnd = std::stoi(parameter.refeGenePos.substr(dashPos + 1));
        } while(0);

        LOG_DEBUG("refPosChrIndex = %d, refPosBegin = %d, refPosEnd = %d", refPosChrIndex, refPosBegin, refPosEnd);
        if (refPosChrIndex != 65535) {
            targeBlock = MemoryUtil::safeNewClass<RoughIOBlock>(outBlockPtr->getBlockSize());
        }
    }

    if (meta.isMember("header")) {
        if (0 != decompressHeader(targeBlock)) {
            LOG_ERROR("Decompress header failed. block id = %d.", inBlockPtr->getBlockId());
            return -1;
        }
    } else {
        LOG_DEBUG("No header info for block: %d", inBlockPtr->getBlockId());
    }

    if (0 != decompressSamByFields(targeBlock)) {
        LOG_ERROR("Decompress fields failed. block id = %d", inBlockPtr->getBlockId());
        return -1;
    }

    if (refPosChrIndex == 65535) {
        // Verify with the algorithm recorded in the block meta (absent/md5 => legacy MD5)
        ChecksumType checksumType = ChecksumType::MD5;
        const std::string expected = getChecksumMeta(meta, checksumType);
        const std::string actual = checksumHex(
            checksumType, outBlockPtr->getBuffer(), (size_t)outBlockPtr->getDataLen());
        if (expected.empty() || actual != expected) {
            LOG_ERROR("checksum check failed for SAM data, blockid = %d, method = %s, expected: %s, got: %s",
                      outBlockPtr->getBlockId(), checksumTypeName(checksumType),
                      expected.c_str(), actual.c_str());
            return -1;
        }
    }
    if (targeBlock != outBlockPtr) {
        MemoryUtil::safeDeleteClass(targeBlock);
    }

    return 0;
}

/* Decode one line's fields in order, appending each; SEQ's landing spot and length travel back in `st` for QUAL. */

int32_t SamCodecActuator::decompressSamByFields(RoughIOBlock* outputBlock) {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr || outputBlock == nullptr) {
        LOG_ERROR("Invalid parameter, inBlockPtr or outputBlock is nullptr for SAM field-by-field decompression");
        return -1;
    }

    if (!meta.isMember("sam")) {
        LOG_INFO("No SAM info for field-by-field decompression");
        return 0;
    }

    // ensureCapacity was already called at the start of decompress(); not repeated here

    // Initialize decoders based on compression metadata
    if (0 != initDecoder(outputBlock)) {
        LOG_ERROR("Init decoder failed.");
        return -1;
    }

    Json::Value& samMeta = meta["sam"];
    Json::Value& streams = samMeta["streams"];
    uint32_t fieldCount = samMeta["fieldcount"].asUInt();
    uint8_t* pBaseEnd = outputBlock->getBuffer() + outputBlock->getBufferSize();

    baseSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    baseDiffSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    /* +8 slack: getStretch2Bits1Char may write outLen+3 bytes. */
    refeStrecchBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength + 8);
    uint32_t totalBaseLen = 0;
    /* The SEQ exception streams and their cursors live in seqExc (see SeqExcStream). */

    uint8_t* pBaseOut = nullptr;
    if (CoderFactory::decoderIsWholeBlock(streams[9]["coder"]["magic"].asString())) {
        /* A whole-block coder (coder_fc) decodes SEQ in one go, staged at the tail of the output block while the head appends lines. */
        pBaseOut = pBaseEnd - streams[9]["totalsrclen"].asUInt();
    }

    /* First pre-decode POS/CIGAR/PNEXT so TLEN reconstruction has the full mate index. */
    posDeltaPrev = 0;
    if (0 != preDecodeForTLEN()) {
        LOG_ERROR("Pre-decode for TLEN failed.");
        return -1;
    }

    /* Copy the pre-decoded field bytes to avoid decoding again; returns -1 when not cached. */

    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        // Decode each field for this line
        uint8_t* basePtr = nullptr;
        uint32_t actualBaseLen = 0;
    auto copyPreDecodedField = [&](uint32_t fieldIdx, uint32_t lineNo) -> int32_t {
        auto preIt = tlenPreDecodedFields.find(fieldIdx);
        if (preIt != tlenPreDecodedFields.end() && lineNo < preIt->second.size() && !preIt->second[lineNo].empty()) {
            const std::string& s = preIt->second[lineNo];
            memcpy(outputBlock->getCurrent(), s.data(), s.length());
            outputBlock->setDataLen(outputBlock->getDataLen() + s.length());
            return (int32_t)s.length();
        }
        return -1;
    };
        for (uint32_t fieldIdx = 0; fieldIdx < fieldCount; ++fieldIdx) {
            int32_t decoderLen = 0;
            if (fieldIdx == 0) {    /// ID
                decoderLen = decompressIdField(fieldIdx, streams[fieldIdx], outputBlock);
            } else if (fieldIdx == 1) {  /// FLAG
                if (streams[1]["mode"].asString() == "string") {
                    decoderLen = decompressRegularField(fieldIdx, lineNo, '\t', outputBlock);
                } else {
                    decoderLen = decompressNumber<uint16_t>(fieldIdx, lineNo, outputBlock);
                }
                /* Store FLAG value for PNEXT validation - retrieve from output */
                if (decoderLen > 1) {
                    uint8_t* flagOutput = outputBlock->getCurrent() - decoderLen;
                    std::string flagStr((char*)flagOutput, (size_t)decoderLen - 1);
                    try {
                        mappedFlag[lineNo] = (uint16_t)std::stoll(flagStr);
                    } catch (...) {
                        mappedFlag[lineNo] = 0;
                    }
                }
            } else if (fieldIdx == 2) {  /// RNAME
                decoderLen = decompressChrName(fieldIdx, lineNo, outputBlock);
            } else if (fieldIdx == 3) {  /// POS
                decoderLen = copyPreDecodedField(3, lineNo);
                if (decoderLen < 0) {
                    if (streams[3]["mode"].asString() == "pos_delta") {
                        decoderLen = decompressPosFieldDelta(fieldIdx, lineNo, '\t', outputBlock);
                    } else if (streams[3]["mode"].asString() == "string") {
                        decoderLen = decompressRegularField(fieldIdx, lineNo, '\t', outputBlock);
                    } else {
                        decoderLen = decompressNumber<uint32_t>(fieldIdx, lineNo, outputBlock);
                    }
                }
            } else if (fieldIdx == 4) {  /// MAPQ
                if (streams[4]["mode"].asString() == "string") {
                    decoderLen = decompressRegularField(fieldIdx, lineNo, '\t', outputBlock);
                } else {
                    decoderLen = decompressNumber<uint8_t>(fieldIdx, lineNo, outputBlock);
                }
            } else if (fieldIdx == 5) {  /// CIGAR
                decoderLen = copyPreDecodedField(5, lineNo);
                if (decoderLen < 0) {
                    decoderLen = decompressCigar(fieldIdx, '\t', lineNo, outputBlock);
                }
            } else if (fieldIdx == 6) {  /// RNEXT
                decoderLen = decompressChrName(fieldIdx, lineNo, outputBlock);
            } else if (fieldIdx == 7) {  /// PNEXT
                decoderLen = copyPreDecodedField(7, lineNo);
                if (decoderLen < 0) {
                    if (streams[7]["mode"].asString() == "pnext_delta") {
                        decoderLen = decompressPNextFieldDelta(fieldIdx, lineNo, '\t', outputBlock);
                    } else if (streams[7]["mode"].asString() == "string") {
                        decoderLen = decompressRegularField(fieldIdx, lineNo, '\t', outputBlock);
                    } else {
                        decoderLen = decompressNumber<uint32_t>(fieldIdx, lineNo, outputBlock);
                    }
                }
            } else if (fieldIdx == 8) {  /// TLEN
                decoderLen = decompressTLen(fieldIdx, lineNo, '\t', outputBlock, streams[8]);
            } else if (fieldIdx == 9) {  /// SEQ
                basePtr = outputBlock->getCurrent();
                decoderLen = decompressBase(fieldIdx, streams[fieldIdx], pBaseOut, lineNo, totalBaseLen, outputBlock);
                actualBaseLen = decoderLen;
            } else if (fieldIdx == 10 ) {  /// QUAL
                decoderLen = decompressQuality(basePtr, actualBaseLen, outputBlock);
                // No optional fields scenario, replace appended \t with \n
                if (fieldIdx + 1 == fieldCount) {
                    uint8_t* pEnd = outputBlock->getCurrent();
                    *(pEnd - 1) = '\n';
                }
            } else if (fieldIdx == 11) {   /// OPTION (all tags)
                /* This record's alignment, for regenerating derived MD/NM. */
                optCtxChrId = (mappedChr.find(lineNo) != mappedChr.end()) ? mappedChr[lineNo] : 0xFFFF;
                optCtxPos = (mappedPos.find(lineNo) != mappedPos.end()) ? (uint64_t)mappedPos[lineNo] : 0;
                optCtxFlag = (mappedFlag.find(lineNo) != mappedFlag.end()) ? mappedFlag[lineNo] : 4;
                optCtxSeq = basePtr;
                optCtxSeqLen = actualBaseLen;
                optCtxContentIdx = (int32_t)(lineNo - headEndLine);
                decoderLen = decompressOptionField(lineNo, '\n', outputBlock, streams[11]);
            } else {   /// Optional fields
                // Decode field until tab or end
                if (fieldIdx + 1 == fieldCount) {
                    decoderLen = decompressRegularField(fieldIdx, lineNo, '\n', outputBlock);
                    // Single newline means it's appended, need to change \t after quality to \n and remove appended \n
                    if (decoderLen == 1) {
                        uint8_t* pEnd = outputBlock->getCurrent();
                        *(pEnd - 2) = '\n';
                        outputBlock->setDataLen(outputBlock->getDataLen() - 1);
                    }
                } else {
                    decoderLen = decompressRegularField(fieldIdx, lineNo, '\t', outputBlock);
                }
            }
            if (decoderLen < 0) {
                LOG_ERROR("Decode field(%d) failed. lineNo = %d", fieldIdx, lineNo);
                return -1;
            }
        }

        if (refPosChrIndex != 65535) {
            if (mappedChr[lineNo] == refPosChrIndex) {
                if ((refPosBegin == 0 && refPosEnd == 0) || (mappedPos[lineNo]  >= refPosBegin && mappedPos[lineNo] <= refPosEnd)) {
                    memcpy(outBlockPtr->getCurrent(), outputBlock->getBuffer(), outputBlock->getDataLen());
                    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + outputBlock->getDataLen());
                }
            }
            outputBlock->reset();
        }
    }

    return 0;
}

int32_t SamCodecActuator::decompressHeader(RoughIOBlock* outputBlock) {
    if (inBlockPtr == nullptr || outputBlock == nullptr) {
        LOG_ERROR("Invalid parameter for SAM header decompression");
        return -1;
    }

    Json::Value& headerMeta = meta["header"];
    if (!headerMeta.isMember("srclen") || !headerMeta.isMember("dstlen") ||
        !headerMeta.isMember("lines") || !headerMeta.isMember("coder")) {
        LOG_ERROR("Invalid SAM header metadata for decompression");
        return -1;
    }

    uint32_t dstLen = headerMeta["dstlen"].asUInt();
    if (refPosChrIndex == 65535) {  // not set position paramter
        // Create SAM file header decompressor from the magic the writer recorded.
        std::shared_ptr<coder_io> headerIo = makeCoderIo(inBlockPtr->getBuffer(), dstLen, "SAM header");
        FieldDecoderArgs headerArgs;
        if (headerMeta["coder"].isMember("level")) {
            headerArgs.level = headerMeta["coder"]["level"].asInt();
        }
        std::shared_ptr<coder> headerDecoder =
            CoderFactory::makeFieldDecoder(headerMeta["coder"]["magic"].asString(),
                                           headerIo.get(), headerArgs);
        if (headerDecoder == nullptr) {
            LOG_ERROR("Unsupported SAM header coder: %s",
                      headerMeta["coder"]["magic"].asString().c_str());
            return -1;
        }

        // Decompress SAM file header data
        uint32_t lineCount = 0;
        uint32_t decoderTotalLen = 0;
        while (lineCount < headEndLine) {
            // Decompress one line of data
            int32_t decodedLen = headerDecoder->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(), '\n', false);
            if (decodedLen < 0) {
                LOG_ERROR("Decode SAM header failed: %d", decodedLen);
                return -1;
            }
            if (decodedLen == 0) {
                break; // No more data
            }
            std::string headStr = std::string((char*)outputBlock->getCurrent(), decodedLen);
            if (headStr.substr(0, 3) == "@SQ") {
                SamUtil::parseChromosomeInfo(headStr);
            }

            outputBlock->setDataLen(outputBlock->getDataLen() + decodedLen);
            lineCount++;
            decoderTotalLen += decodedLen;
        }
        LOG_DEBUG("SAM header decompression completed: %u lines, %u bytes - > %u bytes.", headEndLine, dstLen, decoderTotalLen);
    }
    readOffset += dstLen;
    return 0;
}

/* A "cnt" sub-stream's deltas expanded into the zigzag varint form the numeric layout carries. */
static int32_t expandCounterSubStream(const uint8_t* payload, uint32_t payloadLen,
                                      const Json::Value& streamMeta, uint8_t* out, uint32_t outCap)
{
    RangeCoder rc;
    rc.input((char*)payload, (char*)(payload + payloadLen));
    rc.StartDecode();
    if (rc.err != 0) {
        return -1;
    }
    const uint32_t numValues = streamMeta["numv"].asUInt();
    id_int::Counter counter;
    counter.reset(streamMeta["intw"].asUInt(), streamMeta["intk"].asUInt());
    uint32_t wp = 0;
    for (uint32_t k = 0; k < numValues; ++k) {
        const int64_t dv = counter.decode(rc);
        /* The numeric layout's zigzag form: magnitude in the high bits, sign in the low one. */
        const uint64_t mag = (dv >= 0) ? (uint64_t)dv : (uint64_t)(-dv);
        const uint32_t zz = (uint32_t)((dv >= 0) ? (mag << 1) : ((mag << 1) | 1u));
        if (wp + 5 > outCap) {
            return -1;
        }
        wp += writeVarint(out + wp, zz);
    }
    return (int32_t)wp;
}

/* Build the ID/QNAME column: one decoder per split sub-stream, plus the payload state each layout needs. */
int32_t SamCodecActuator::initIdFieldDecoders(Json::Value& idMeta)
{
    std::string idSplit = idMeta["splitsym"].asString();
    for (uint32_t i = 0; i < idSplit.length(); ++i) {
        idAnalysis.symbols.push_back(idSplit.c_str()[i]);
    }
    Json::Value& idStreamMeta = idMeta["streams"];
    if (idStreamMeta.size() != idAnalysis.symbols.size()) {
        LOG_ERROR("id streams not match id split, expect:%u, actual:%u", idAnalysis.symbols.size(), idStreamMeta.size());
        return -1;
    }

    // ID decoders
    idStreamOffsets.clear();
    idStreamDstLens.clear();
    idStreamCoders.clear();
    clearIdNumericState();
    for (uint32_t i = 0; i < idStreamMeta.size(); ++i) {
        std::string coderName = idStreamMeta[i]["coder"]["magic"].asString();
        uint32_t dstLength = idStreamMeta[i]["dstlen"].asUInt();
        idStreamOffsets.push_back(readOffset);
        idStreamDstLens.push_back(dstLength);
        idStreamCoders.push_back(coderName);
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "QNAME sub-stream");
        ioVector.push_back(io);
        const bool isNumeric = idStreamMeta[i].isMember("mode") &&
                               idStreamMeta[i]["mode"].asString() == "numeric";
        const bool isInt = idStreamMeta[i].isMember("mode") &&
                           idStreamMeta[i]["mode"].asString() == "intd";
        const bool isConst = idStreamMeta[i].isMember("mode") &&
                             idStreamMeta[i]["mode"].asString() == "const";
        const bool isCnt = idStreamMeta[i].isMember("mode") &&
                           idStreamMeta[i]["mode"].asString() == "cnt";
        const bool isIntu = idStreamMeta[i].isMember("mode") &&
                            idStreamMeta[i]["mode"].asString() == "intu";
        const bool isDict = idStreamMeta[i].isMember("mode") &&
                            idStreamMeta[i]["mode"].asString() == "dict";
        const bool isHex = idStreamMeta[i].isMember("mode") &&
                           idStreamMeta[i]["mode"].asString() == "hexd";
        if (isConst) {
            /* The payload is the segment's own text, stored once for the whole block. */
            idDecoders.push_back(nullptr);
        } else if (isDict) {
            /* The payload is a dictionary plus an index per line: no coder of its own. */
            idDecoders.push_back(nullptr);
        } else if (isHex) {
            /* The payload is uniform character codes: no coder of its own either. */
            idDecoders.push_back(nullptr);
        } else if (isCnt) {
            /* The payload is a coded delta sequence: neither text nor a varint buffer. */
            idDecoders.push_back(nullptr);
        } else if (isInt || isIntu) {
            /* The payload is the value-domain stream itself (no outer coder). */
            idDecoders.push_back(nullptr);
        } else {
            /* The sub-stream's decoder comes from its own meta (see CoderFactory::makeFieldDecoder).
               This covers the archive layouts and the fast profile's coder_rans alike. */
            FieldDecoderArgs idArgs;
            if (idStreamMeta[i]["coder"].isMember("level")) {
                idArgs.level = idStreamMeta[i]["coder"]["level"].asInt();
            }
            std::shared_ptr<coder> idDec =
                CoderFactory::makeFieldDecoder(coderName, io.get(), idArgs);
            if (idDec == nullptr) {
                LOG_ERROR("Unsupport coder name:%s", coderName.c_str());
                return -1;
            }
            if (coderName == "coder_qname") {
                idUsesQnameCoder = true;
            }
            idDecoders.push_back(idDec);
        }

        /* Numeric sub-streams have no per-line terminator, so the whole varint stream is decoded once and served per line. */
        uint8_t* nbuf = nullptr;
        uint32_t nlen = 0;
        RangeCoder intRc;
        id_int::Model intModel;
        uint8_t intActive = 0;
        std::vector<std::string> dictEntries;
        IdHexAlphabet hexAlpha;
        RangeCoder hexRc;
        if (isNumeric) {
            const uint32_t numSrcLen = idStreamMeta[i]["srclen"].asUInt();
            nbuf = MemoryUtil::safeAlloc<uint8_t>(numSrcLen + 1);
            if (nbuf == nullptr) {
                return -1;
            }
            int32_t dl = idDecoders[i]->decode_line(nbuf, numSrcLen, UINT8_MAX, false);
            if (dl < 0) {
                MemoryUtil::safeFree(nbuf);
                LOG_ERROR("Decode id numeric sub-stream(%u) failed: %d", i, dl);
                return -1;
            }
            nlen = (uint32_t)dl;
        } else if (isInt || isIntu) {
            nbuf = MemoryUtil::safeAlloc<uint8_t>(dstLength + 1);
            if (nbuf == nullptr) {
                return -1;
            }
            memcpy(nbuf, inBlockPtr->getBuffer() + readOffset, dstLength);
            nlen = dstLength;
            if (isIntu) {
                /* The values' exact range, with no model at all (see Model::resetUniform). */
                intModel.resetUniform(idStreamMeta[i]["intb"].asUInt(),
                                      idStreamMeta[i]["intn"].asUInt());
            } else {
                intModel.reset(idStreamMeta[i]["intw"].asUInt(), idStreamMeta[i]["intk"].asUInt());
            }
            intRc.input((char*)nbuf, (char*)nbuf + nlen);
            intRc.StartDecode();
            if (intRc.err != 0) {
                MemoryUtil::safeFree(nbuf);
                LOG_ERROR("Decode id value-domain sub-stream(%u) failed", i);
                return -1;
            }
            intActive = 1;
        } else if (isCnt) {
            /* A counter's deltas expanded into the same zigzag varint form as the numeric layout. */
            const uint32_t numSrcLen = idStreamMeta[i]["srclen"].asUInt();
            nbuf = MemoryUtil::safeAlloc<uint8_t>(numSrcLen + 8);
            if (nbuf == nullptr) {
                return -1;
            }
            const int32_t wp = expandCounterSubStream(inBlockPtr->getBuffer() + readOffset,
                                                      dstLength, idStreamMeta[i], nbuf,
                                                      numSrcLen + 8);
            if (wp < 0) {
                MemoryUtil::safeFree(nbuf);
                LOG_ERROR("Decode id counter sub-stream(%u) failed", i);
                return -1;
            }
            nlen = (uint32_t)wp;
        } else if (isDict) {
            /* Low-cardinality segment: distinct texts at the head of the payload, one index per line through id_int::Model. */
            const uint8_t* p = inBlockPtr->getBuffer() + readOffset;
            uint32_t pos = 0;
            uint32_t num = 0;
            if (!readVarint(p, dstLength, pos, num) || num < 2 || num > kMaxDictEntries) {
                LOG_ERROR("Decode id dictionary sub-stream(%u): bad entry count", i);
                return -1;
            }
            for (uint32_t e = 0; e < num; ++e) {
                uint32_t len = 0;
                if (!readVarint(p, dstLength, pos, len) || len > dstLength - pos) {
                    LOG_ERROR("Decode id dictionary sub-stream(%u): truncated entry", i);
                    return -1;
                }
                dictEntries.push_back(std::string((const char*)p + pos, len));
                pos += len;
            }
            /* The index stream is a counter layout; its deltas are expanded into one absolute index per line. */
            const uint32_t numValues = idStreamMeta[i]["numv"].asUInt();
            const uint8_t* ip = p + pos;
            const uint32_t ilen = dstLength - pos;
            RangeCoder idxRc;
            idxRc.input((char*)ip, (char*)(ip + ilen));
            idxRc.StartDecode();
            if (idxRc.err != 0) {
                LOG_ERROR("Decode id dictionary sub-stream(%u) failed", i);
                return -1;
            }
            const uint32_t ibufLen = numValues * 5 + 8;
            nbuf = MemoryUtil::safeAlloc<uint8_t>(ibufLen);
            if (nbuf == nullptr) {
                return -1;
            }
            id_int::Counter idxCounter;
            idxCounter.reset(idStreamMeta[i]["intw"].asUInt(), idStreamMeta[i]["intk"].asUInt());
            uint32_t wp = 0;
            int64_t acc = 0;
            for (uint32_t k = 0; k < numValues; ++k) {
                acc += idxCounter.decode(idxRc);
                if (acc < 0 || (uint64_t)acc >= dictEntries.size() || wp + 5 > ibufLen) {
                    MemoryUtil::safeFree(nbuf);
                    LOG_ERROR("Decode id dictionary sub-stream(%u): index out of range", i);
                    return -1;
                }
                wp += writeVarint(nbuf + wp, (uint32_t)acc);
            }
            nlen = wp;
        } else if (isHex) {
            /* Fixed-alphabet characters: one uniform code per character, plus a length per line when lengths vary. */
            hexAlpha.chars = idStreamMeta[i]["hexa"].asString();
            hexAlpha.minLen = idStreamMeta[i]["hexn"].asUInt();
            hexAlpha.maxLen = idStreamMeta[i]["hexm"].asUInt();
            if (hexAlpha.chars.size() < 2 || hexAlpha.maxLen < hexAlpha.minLen ||
                hexAlpha.maxLen > kMaxHexLen) {
                LOG_ERROR("Decode id fixed-alphabet sub-stream(%u): bad alphabet", i);
                return -1;
            }
            nbuf = MemoryUtil::safeAlloc<uint8_t>(dstLength + 1);
            if (nbuf == nullptr) {
                return -1;
            }
            memcpy(nbuf, inBlockPtr->getBuffer() + readOffset, dstLength);
            nlen = dstLength;
            hexRc.input((char*)nbuf, (char*)nbuf + nlen);
            hexRc.StartDecode();
            if (hexRc.err != 0) {
                MemoryUtil::safeFree(nbuf);
                LOG_ERROR("Decode id fixed-alphabet sub-stream(%u) failed", i);
                return -1;
            }
        }
        idNumericBufs.push_back(nbuf);
        idNumericLens.push_back(nlen);
        idNumericPos.push_back(0);
        idNumericAcc.push_back(0);
        idIntCoders.push_back(intRc);
        idIntModes.push_back(intModel);
        idIntActive.push_back(intActive);
        if (isConst) {
            /* No coder and no per-line data: the segment text is the payload, verbatim. */
            idConstTexts.push_back(
                std::string((const char*)inBlockPtr->getBuffer() + readOffset, dstLength));
        } else {
            idConstTexts.push_back(std::string());
        }
        idDictEntries.push_back(std::move(dictEntries));
        idHexAlphabets.push_back(std::move(hexAlpha));
        idHexCoders.push_back(hexRc);

        readOffset += dstLength;
    }
    return 0;
}

/* Build the SEQ column: reference layout (match, exception, base-length streams) or the single whole-column stream. */
int32_t SamCodecActuator::initSeqFieldDecoders(Json::Value& baseMeta, uint32_t idx,
                                              RoughIOBlock* outputBlock)
{
    maxBaseLength = baseMeta["maxlen"].asUInt();
    minBaseLength = baseMeta["minlen"].asUInt();
    LOG_DEBUG("maxBaseLen = %d, minBaseLen = %d", maxBaseLength, minBaseLength);
    /* A block-consensus column uses the same layout as the reference one but needs no external reference. */
    baseSelfRefActive = false;
    /* The multi-stream layout is identified by the meta itself (a "selfref" consensus, or a reference
       layout that may be a file reference or, when "noref", no reference at all). */
    const bool usesStreams = baseMeta.isMember("selfref") || baseMeta.isMember("streams");
    if (usesStreams && !baseMeta.isMember("selfref") &&
        !(baseMeta.isMember("noref") && baseMeta["noref"].asUInt() != 0) && pRefeGene == nullptr) {
        LOG_ERROR("SEQ block %lld was coded against a reference but none was supplied",
                  (long long)inBlockPtr->getBlockId());
        return -1;
    }
    if (!usesStreams) {
        /* One stream for the whole column; initSeqWholeBlockDecoder builds its decoder or decodes it outright. */
        const int32_t consumed = initSeqWholeBlockDecoder(baseMeta, idx, outputBlock);
        if (consumed < 0) {
            return -1;
        }
        readOffset += (uint32_t)consumed;
    } else {
        uint32_t id = 0;
        // Scenario using reference genome
        Json::Value& baseMetaStreams = baseMeta["streams"];
        if (baseMetaStreams[id]["sname"] != "m") {
            LOG_ERROR("check sub stream failed:%s", baseMetaStreams[id]["sname"].asString().c_str());
            return -1;
        }

        uint32_t dstLength = baseMetaStreams[id]["dstlen"].asUInt();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "SEQ match");
        /* The sub-stream's meta travels with the view: a whole-block coder (coder_fc) needs it to decode. */
        io->meta = baseMetaStreams[id];
        ioVector.push_back(io);
        const bool rleEnabled = baseMetaStreams[id].isMember("rle") &&
                                baseMetaStreams[id]["rle"].asUInt() != 0;
        if (rleEnabled) {
            if (predecodeSeqMatchStream(baseMetaStreams, id, dstLength, io.get(), idx,
                                        matchExtraOffset) != 0) {
                return -1;
            }
        } else {
            if (predecodeSeqWholeMatchStream(baseMetaStreams, id, io.get(), idx) != 0) {
                return -1;
            }
        }

        readOffset += dstLength + matchExtraOffset;   /* +mval sub-stream under RLE */
        matchExtraOffset = 0;
        /* SEQ exception streams (see SeqExceptionClass): one per character present, a strictly increasing offset list. */
        if (predecodeSeqExceptionStreams(baseMeta, baseMetaStreams, id) != 0) {
            return -1;
        }

        if (minBaseLength != maxBaseLength) {
            id++;
            if (predecodeSeqBaseLengthStream(baseMetaStreams, id) != 0) {
                return -1;
            }
        }

        /* The block's own consensus, when used: its squash is the last sub-stream and the meta gives its window. */
        if (baseMeta.isMember("selfref")) {
            id++;
            if (id >= (uint32_t)baseMetaStreams.size() ||
                predecodeSeqSelfRefStream(baseMetaStreams, id) != 0) {
                LOG_ERROR("SEQ self reference stream is missing or unreadable");
                return -1;
            }
            baseSelfRefWindowStart = baseMeta["selfref_wstart"].asInt64();
            baseSelfRefWindowBases = baseMeta["selfref_wbases"].asInt64();
            baseSelfRefActive = true;
        }
    }
    return 0;
}

/* Record the stream offset of a lazily-decoded field and report whether this was such a field. */
bool SamCodecActuator::recordDeferredFieldOffset(uint32_t idx, Json::Value& fieldMeta)
{
    if (idx == 11 && fieldMeta.isMember("mode") &&
        (fieldMeta["mode"].asString() == "tag_split" ||
         fieldMeta["mode"].asString() == "opt_strip")) {
        /* OPTION is the last field, so its streams need no advancing. */
        fieldIoStart[idx] = readOffset;
        return true;
    }
    if (idx == 7 && fieldMeta.isMember("mode") &&
        fieldMeta["mode"].asString() == "pnext_qname_rebuild") {
        fieldIoStart[idx] = readOffset;
        if (fieldMeta.isMember("streams") && fieldMeta["streams"].size() > 0) {
            const uint32_t excDstLen = fieldMeta["streams"][0]["dstlen"].asUInt();
            readOffset += excDstLen;
        }
        return true;
    }
    return false;
}

int32_t SamCodecActuator::initDecoder(RoughIOBlock* outputBlock) {
    if (outputBlock == nullptr) {
        return -1;
    }
    Json::Value& streamMeta = meta["sam"]["streams"];
    if (!streamMeta.isArray()) {
        return -1;
    }

    uint32_t lineNumber = meta["sam"]["lines"].asUInt();
    LOG_DEBUG("Line number = %d", lineNumber);
    baseLengthBuffer = MemoryUtil::safeAlloc<uint32_t>(lineNumber);

    for (uint32_t idx = 0; idx < streamMeta.size(); ++idx) {
        if (idx == 0) {
            if (initIdFieldDecoders(streamMeta[idx]) != 0) {
                return -1;
            }
        } else if (idx == 9) {
            if (initSeqFieldDecoders(streamMeta[idx], idx, outputBlock) != 0) {
                return -1;
            }
        } else if (idx == 10) {
            /* The whole column: its frequency table, then the decoder its magic names. */
            if (initQualFieldDecoders(streamMeta[idx]) != 0) {
                return -1;
            }
        } else if (idx == 8) {
            /* TLEN: the reconstruction layout keeps only an exception stream (see
               initTlenFieldDecoders). */
            if (initTlenFieldDecoders(streamMeta[idx], idx) != 0) {
                return -1;
            }
        } else if (recordDeferredFieldOffset(idx, streamMeta[idx])) {
            continue;   /* decoded lazily: only the stream's offset is kept */
        } else {
            if (initFieldDecoder(idx, streamMeta[idx]) != 0) {
                return -1;
            }
        }
    }
    return 0;
}

int32_t SamCodecActuator::decompressRegularField(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock) {
    /* A negative decode_line return is an error code (corrupted stream / insufficient buffer); pass it through, never treat it as a length. */
    int32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(), splitFlag, false);
    if (fieldLen < 0) {
        LOG_ERROR("Decode regular field(%u) failed: %d", fieldIdx, fieldLen);
        return -1;
    }
    outputBlock->setDataLen(outputBlock->getDataLen() + fieldLen);

    /* Textual FLAG/POS/PNEXT must repopulate the tracking maps here, consistent with decompressNumber. */
    if (fieldLen > 1 && (fieldIdx == 1 || fieldIdx == 3 || fieldIdx == 7)) {
        uint8_t* vs = outputBlock->getBuffer() + outputBlock->getDataLen() - fieldLen;
        std::string str((char*)vs, (size_t)fieldLen - 1);
        if (fieldIdx == 1) {
            mappedFlag[lineNo] = (uint16_t)std::stoull(str);
        } else if (fieldIdx == 3) {
            mappedPos[lineNo] = (int64_t)std::stoull(str);
        } else if (fieldIdx == 7) {
            nextMappedPos[lineNo] = (int64_t)std::stoull(str);
        }
    }
    return fieldLen;
}

/* OPTION field (disabled tag-columnization path): decode all columns once and cache them, then emit per line. */
int32_t SamCodecActuator::decompressOptionField(uint32_t lineNo, uint8_t splitFlag,
                                                RoughIOBlock* outputBlock,
                                                const Json::Value& fieldMeta) {
    /* The stripped layout: rebuild each record from its sequence index and values (see compressOptionStrip). */
    if (fieldMeta.isMember("mode") && fieldMeta["mode"].asString() == "opt_strip") {
        if (optionCacheEmpty) {
            if (decodeOptionStripColumn(fieldMeta) != 0) {
                return -1;
            }
        }
        /* lineNo counts data lines here, one per record. */
        const uint32_t recIdx = lineNo;
        if (recIdx + 1 < optionRecOff.size() && recIdx < optionRecPresent.size() &&
            optionRecPresent[recIdx] != 0) {
            const uint32_t from = optionRecOff[recIdx];
            const uint32_t len = optionRecOff[recIdx + 1] - from;
            memcpy(outputBlock->getCurrent(), optionRecText.data() + from, len);
            outputBlock->setDataLen(outputBlock->getDataLen() + len);
            *(outputBlock->getCurrent()) = splitFlag;
            outputBlock->setDataLen(outputBlock->getDataLen() + 1);
            return (int32_t)len + 1;
        }
        /* This line has no OPTION: turn the '\t' just appended after QUAL back into '\n'. */
        uint8_t* pEnd = outputBlock->getCurrent();
        if (pEnd > outputBlock->getBuffer()) {
            *(pEnd - 1) = '\n';
        }
        return 0;
    }

    if (!fieldMeta.isMember("tags") || fieldMeta["mode"].asString() != "tag_split") {
        /* affix form: the whole block's OPTION is one column; a line with no OPTION yielded a bare delimiter. */
        const int32_t fieldLen = decompressRegularField(11, lineNo, splitFlag, outputBlock);
        if (fieldLen == 1) {
            uint8_t* pEnd = outputBlock->getCurrent();
            if (pEnd - 2 >= outputBlock->getBuffer()) {
                *(pEnd - 2) = '\n';
                outputBlock->setDataLen(outputBlock->getDataLen() - 1);
            }
        }
        return fieldLen;
    }

    /* The per-tag columns are decoded once per block; each line is reassembled here so a derived
       MD/NM can be recomputed from this record's decoded alignment. */
    if (optionCacheEmpty) {
        if (0 != decodeOptionColumn(fieldMeta)) {
            return -1;
        }
    }
    const std::vector<uint8_t>* ids =
        (lineNo < optionRecIds.size()) ? &optionRecIds[lineNo] : nullptr;
    if (ids == nullptr || ids->empty()) {
        /* This line has no OPTION: turn the '\t' just appended after QUAL back into '\n'. */
        uint8_t* pEnd = outputBlock->getCurrent();
        if (pEnd > outputBlock->getBuffer()) {
            *(pEnd - 1) = splitFlag;
        }
        return 0;
    }

    std::string line;
    bool first = true;
    for (uint8_t tid : *ids) {
        if (tid >= optionTagCols.size()) {
            continue;
        }
        OptionTagCol& col = optionTagCols[tid];
        std::string value;
        if (col.enc == "derived") {
            if (col.cursor >= col.flags.size()) {
                continue;
            }
            const uint8_t flag = col.flags[col.cursor++];
            if (flag == 0) {
                if (!deriveOptionTagValue(tid, value)) {
                    LOG_ERROR("OPTION tag %s: cannot regenerate at decode (block %lld line %u)",
                              optionTagNames[tid].c_str(), (long long)inBlockPtr->getBlockId(),
                              lineNo);
                    return -1;
                }
            } else {
                if (col.excCursor >= col.excs.size()) {
                    LOG_ERROR("OPTION derived column %s ran out of exceptions",
                              optionTagNames[tid].c_str());
                    return -1;
                }
                value = col.excs[col.excCursor++];
            }
        } else {
            if (col.cursor >= col.vals.size()) {
                continue;
            }
            value = col.vals[col.cursor++];
        }
        if (!first) {
            line.push_back('\t');
        }
        first = false;
        line.append(optionTagNames[tid]);
        line.push_back(':');
        line.append(optionTagTypes[tid]);
        line.push_back(':');
        line.append(value);
    }
    if (!line.empty()) {
        memcpy(outputBlock->getCurrent(), line.data(), line.size());
        outputBlock->setDataLen(outputBlock->getDataLen() + line.size());
    }
    *(outputBlock->getCurrent()) = splitFlag;
    outputBlock->setDataLen(outputBlock->getDataLen() + 1);
    return (int32_t)line.size() + 1;
}

int32_t SamCodecActuator::decodeOptionColumn(const Json::Value& fieldMeta) {
    optionRecIds.clear();
    optionTagNames.clear();
    optionTagTypes.clear();
    optionTagCols.clear();
    optionCacheEmpty = false;

    /* Start of this field's stream: recorded by initDecoder; readOffset cannot be used (it has been advanced by line-by-line decoding). */
    uint32_t optBase = readOffset;
    auto optIt = fieldIoStart.find(11);
    if (optIt != fieldIoStart.end()) {
        optBase = optIt->second;
    }

    const Json::Value& streams = fieldMeta["streams"];
    if (!streams.isArray()) {
        return -1;
    }

    const Json::Value& tags = fieldMeta["tags"];
    const uint32_t nTag = (uint32_t)tags.size();
    optionTagNames.resize(nTag);
    optionTagTypes.resize(nTag);
    optionTagCols.assign(nTag, OptionTagCol());
    for (uint32_t t = 0; t < nTag; ++t) {
        optionTagNames[t] = tags[t][0].asString();
        optionTagTypes[t] = tags[t][1].asString();
    }
    if (nTag == 0) {
        /* A block with no tags at all: every record has no OPTION column. */
        return 0;
    }

    /* Decode the id-sequence column (line-wise text, comma-separated). */
    uint32_t idStreamDst = 0;
    const Json::Value* idsMeta = nullptr;
    for (uint32_t i = 0; i < streams.size(); ++i) {
        if (streams[i]["sname"].asString() == "ids") {
            idStreamDst = streams[i]["dstlen"].asUInt();
            idsMeta = &streams[i];
            break;
        }
    }
    {
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + optBase, idStreamDst, "OPTION ids");
        ioVector.push_back(io);
        FieldDecoderArgs args;
        if (idsMeta != nullptr && (*idsMeta)["coder"].isMember("level")) {
            args.level = (*idsMeta)["coder"]["level"].asInt();
        }
        const std::string magic = (idsMeta != nullptr)
                                      ? (*idsMeta)["coder"]["magic"].asString() : std::string();
        std::shared_ptr<coder> c = CoderFactory::makeFieldDecoder(magic, io.get(), args);
        if (c == nullptr) {
            LOG_ERROR("Unsupported OPTION ids coder: %s", magic.c_str());
            return -1;
        }
        uint8_t tmp[1 << 16];
        int32_t n;
        std::string line;
        while ((n = c->decode_line(tmp, sizeof(tmp), '\n', false)) > 0) {
            line.assign((char*)tmp, (size_t)n);
            if (!line.empty() && line.back() == '\n') line.pop_back();
            std::vector<uint8_t> ids;
            size_t p = 0;
            while (p < line.size()) {
                size_t q = line.find(',', p);
                if (q == std::string::npos) q = line.size();
                std::string tok = line.substr(p, q - p);
                if (!tok.empty()) ids.push_back((uint8_t)atoi(tok.c_str()));
                p = q + 1;
            }
            optionRecIds.push_back(std::move(ids));
        }
        optBase += idStreamDst;
    }

    /* Decode each tag-value column: line-oriented text, a varint blob, or the derived layout
       (one flag byte per value followed by the exception values). */
    for (uint32_t i = 0; i < streams.size(); ++i) {
        const Json::Value& s = streams[i];
        if (s["sname"].asString() != "tag") continue;
        uint32_t tagIdx = 0;
        for (uint32_t t = 0; t < nTag; ++t) {
            if (optionTagNames[t] == s["tag"].asString()) { tagIdx = t; break; }
        }
        uint32_t dstlen = s["dstlen"].asUInt();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + optBase, dstlen, "OPTION tag");
        ioVector.push_back(io);
        FieldDecoderArgs args;
        if (s["coder"].isMember("level")) {
            args.level = s["coder"]["level"].asInt();
        }
        const std::string magic = s["coder"]["magic"].asString();
        std::shared_ptr<coder> c = CoderFactory::makeFieldDecoder(magic, io.get(), args);
        if (c == nullptr) {
            LOG_ERROR("Unsupported OPTION tag coder: %s", magic.c_str());
            return -1;
        }
        const std::string enc = s["enc"].isString() ? s["enc"].asString() : std::string("text");
        OptionTagCol& col = optionTagCols[tagIdx];
        col.enc = enc;
        col.cursor = 0;
        col.excCursor = 0;
        if (enc == "text") {
            uint8_t tmp[1 << 16];
            int32_t n;
            while ((n = c->decode_line(tmp, sizeof(tmp), '\n', false)) > 0) {
                size_t tlen = (size_t)n;
                if (tlen > 0 && tmp[tlen-1] == '\n') tlen--;
                col.vals.emplace_back((char*)tmp, tlen);
            }
        } else {
            /* Whole-blob column (varint or derived). */
            const uint32_t srclen = s["srclen"].asUInt();
            std::vector<uint8_t> dec(srclen);
            uint32_t got = 0;
            while (got < srclen) {
                const int32_t n = c->decode_line(dec.data() + got, srclen - got, UINT8_MAX, false);
                if (n <= 0) break;
                got += (uint32_t)n;
            }
            if (got != srclen) {
                LOG_ERROR("OPTION column '%s' is short: %u of %u bytes",
                          s["tag"].asString().c_str(), got, srclen);
                return -1;
            }
            if (enc == "derived") {
                const uint32_t n = s["n"].asUInt();
                if (n > srclen) {
                    LOG_ERROR("OPTION derived column '%s' flag count %u exceeds %u bytes",
                              s["tag"].asString().c_str(), n, srclen);
                    return -1;
                }
                col.flags.assign(dec.begin(), dec.begin() + n);
                size_t off = n;
                while (off < srclen) {
                    size_t e = off;
                    while (e < srclen && dec[e] != '\n') ++e;
                    col.excs.emplace_back((char*)dec.data() + off, e - off);
                    off = e + 1;
                }
            } else {
                uint32_t off = 0;
                int64_t prev = 0;
                while (off < srclen) {
                    uint64_t zz = 0;
                    if (!readVarint64(dec.data(), srclen, off, zz)) {
                        LOG_ERROR("OPTION integer column '%s' has a malformed varint",
                                  s["tag"].asString().c_str());
                        return -1;
                    }
                    const int64_t v = (enc == "dvarint") ? (prev + unzigzag64(zz)) : (int64_t)zz;
                    prev = v;
                    char num[24];
                    const int nlen = snprintf(num, sizeof(num), "%lld", (long long)v);
                    col.vals.emplace_back(num, (size_t)nlen);
                }
            }
        }
        optBase += dstlen;
    }
    return 0;
}

/* Recompute a derived tag (MD:Z / NM:i) for the record currently being decoded. Returns false if
   the tag is not derivable or the alignment cannot be recomputed. */
bool SamCodecActuator::deriveOptionTagValue(uint32_t tagIdx, std::string& out)
{
    if (tagIdx >= optionTagNames.size()) {
        return false;
    }
    const std::string& name = optionTagNames[tagIdx];
    const std::string& type = optionTagTypes[tagIdx];
    const bool isMd = (name == "MD" && type == "Z");
    const bool isNm = (name == "NM" && (type == "i" || type == "I"));
    if ((!isMd && !isNm) || pRefeGene == nullptr || optCtxSeq == nullptr ||
        optCtxContentIdx < 0 || optCtxChrId == 0xFFFF || optCtxChrId == 0xFFFE) {
        return false;
    }
    const int64_t chrStart = SamInfo::getInstance().getPositionByIndex(optCtxChrId);
    if (chrStart < 0) {
        return false;
    }
    static const std::vector<CigarOp> kEmptyOps;
    const std::vector<CigarOp>& ops =
        ((size_t)optCtxContentIdx < cigarOpList.size()) ? cigarOpList[optCtxContentIdx] : kEmptyOps;
    std::string md;
    int64_t nm = 0;
    if (!deriveMdNm(pRefeGene, chrStart, optCtxPos, optCtxFlag, ops, optCtxSeq, optCtxSeqLen,
                    optMdScratch, md, nm)) {
        return false;
    }
    out = isMd ? md : std::to_string((long long)nm);
    return true;
}

/* Read the stripped OPTION layout back: one sequence index and one line of values per record, both line-wise. */
int32_t SamCodecActuator::decodeOptionStripColumn(const Json::Value& fieldMeta)
{
    optionRecText.clear();
    optionRecOff.clear();
    optionRecPresent.clear();
    optionCacheEmpty = false;

    const Json::Value& streams = fieldMeta["streams"];
    const Json::Value& seqs = fieldMeta["seqs"];
    if (!streams.isArray() || !seqs.isArray()) {
        LOG_ERROR("OPTION stripped layout is missing its sequence table or streams");
        return -1;
    }

    /* This field's stream start, recorded by recordDeferredFieldOffset. */
    uint32_t off = readOffset;
    auto startIt = fieldIoStart.find(11);
    if (startIt != fieldIoStart.end()) {
        off = startIt->second;
    }

    /* The two sub-streams sit one after the other, in the order the writer emitted them. */
    uint32_t sidxOff = 0, sidxLen = 0, valsOff = 0, valsLen = 0;
    const Json::Value* sidxMeta = nullptr;
    const Json::Value* valsMeta = nullptr;
    for (uint32_t i = 0; i < streams.size(); ++i) {
        const uint32_t dst = streams[i]["dstlen"].asUInt();
        const std::string name = streams[i]["sname"].asString();
        if (name == "sidx") {
            sidxOff = off;
            sidxLen = dst;
            sidxMeta = &streams[i];
        } else if (name == "vals") {
            valsOff = off;
            valsLen = dst;
            valsMeta = &streams[i];
        }
        off += dst;
    }
    if (sidxMeta == nullptr || valsMeta == nullptr) {
        LOG_ERROR("OPTION stripped layout is missing a sub-stream");
        return -1;
    }

    /* The level travels as in initFieldDecoder. */
    auto openSub = [&](const Json::Value& meta, uint32_t at, uint32_t len) -> std::shared_ptr<coder> {
        const std::string coderName = meta["coder"]["magic"].asString();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + at, len, "OPTION stream");
        ioVector.push_back(io);
        FieldDecoderArgs args;
        if (coderName == "coder_affix_match" ||
            (coderName == "coder_arith" && meta["coder"].isMember("level"))) {
            args.level = meta["coder"]["level"].asInt();
        }
        return CoderFactory::makeFieldDecoder(coderName, io.get(), args);
    };

    std::shared_ptr<coder> sidxCoder = openSub(*sidxMeta, sidxOff, sidxLen);
    std::shared_ptr<coder> valsCoder = openSub(*valsMeta, valsOff, valsLen);
    if (sidxCoder == nullptr || valsCoder == nullptr) {
        LOG_ERROR("Unsupported OPTION sub-stream coder");
        return -1;
    }

    /* The sequence table's texts, read out of the meta once per sequence (seqText[i][k] is segment k of sequence i). */
    struct SeqSegment {
        bool verbatim = false;
        std::string prefix;   /* "NAME:TYPE:", empty when the segment is verbatim */
    };
    std::vector<std::vector<SeqSegment>> seqText((size_t)seqs.size());
    std::vector<uint8_t> seqReady((size_t)seqs.size(), 0);

    /* optionRecOff is a prefix sum; its leading zero precedes the records. */
    optionRecOff.push_back(0);
    uint8_t tmp[1 << 16];
    std::vector<uint32_t> valueStarts;
    while (true) {
        const int32_t idxLen = sidxCoder->decode_line(tmp, sizeof(tmp), '\n', false);
        if (idxLen <= 0) {
            break;   /* end of the column */
        }
        uint32_t idxEnd = (uint32_t)idxLen;
        if (idxEnd > 0 && tmp[idxEnd - 1] == '\n') {
            --idxEnd;
        }
        uint32_t seqIdx = 0;
        for (uint32_t k = 0; k < idxEnd; ++k) {
            if (tmp[k] < '0' || tmp[k] > '9') {
                LOG_ERROR("OPTION sequence index is not a number");
                return -1;
            }
            seqIdx = seqIdx * 10 + (uint32_t)(tmp[k] - '0');
        }
        if (seqIdx >= (uint32_t)seqs.size()) {
            LOG_ERROR("OPTION sequence index %u is out of range", seqIdx);
            return -1;
        }

        const int32_t valsGot = valsCoder->decode_line(tmp, sizeof(tmp), '\n', false);
        if (valsGot <= 0) {
            LOG_ERROR("OPTION values line is missing");
            return -1;
        }
        uint32_t valsEnd = (uint32_t)valsGot;
        if (valsEnd > 0 && tmp[valsEnd - 1] == '\n') {
            --valsEnd;
        }

        /* The value count is fixed by the sequence, not the line's length. */
        const Json::Value& seq = seqs[seqIdx];
        if (!seqReady[seqIdx]) {
            const uint32_t segs = (uint32_t)seq.size();
            seqText[seqIdx].resize(segs);
            for (uint32_t k = 0; k < segs; ++k) {
                SeqSegment& e = seqText[seqIdx][k];
                e.verbatim = (seq[k][0].asInt() == 1);
                if (!e.verbatim) {
                    e.prefix = seq[k][1].asString() + ':' + seq[k][2].asString() + ':';
                }
            }
            seqReady[seqIdx] = 1;
        }
        valueStarts.clear();
        if (seq.size() > 0) {
            uint32_t start = 0;
            for (uint32_t k = 0; k <= valsEnd; ++k) {
                if (k != valsEnd && tmp[k] != '\t') {
                    continue;
                }
                valueStarts.push_back(start);
                start = k + 1;
            }
        }

        if (seq.size() != valueStarts.size()) {
            LOG_ERROR("OPTION sequence has %u segments but the record carries %zu values",
                      (unsigned)seq.size(), valueStarts.size());
            return -1;
        }
        for (uint32_t k = 0; k < (uint32_t)valueStarts.size(); ++k) {
            const uint32_t from = valueStarts[k];
            const uint32_t to = (k + 1 < (uint32_t)valueStarts.size())
                ? valueStarts[k + 1] - 1 : valsEnd;   /* minus the tab that separates them */
            if (k) {
                optionRecText.push_back('\t');
            }
            const SeqSegment& e = seqText[seqIdx][k];
            if (!e.verbatim) {
                optionRecText.insert(optionRecText.end(), e.prefix.begin(), e.prefix.end());
            }
            optionRecText.insert(optionRecText.end(), tmp + from, tmp + to);
        }
        optionRecOff.push_back((uint32_t)optionRecText.size());
        optionRecPresent.push_back(seq.size() > 0 ? 1 : 0);
    }
    return 0;
}

int32_t SamCodecActuator::decompressPNextFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock) {
    uint8_t deltaBuffer[32] = {0};
    int32_t deltaLen = fieldDecoders[fieldIdx]->decode_line(deltaBuffer, sizeof(deltaBuffer), splitFlag, true);
    if (deltaLen < 0) {
        LOG_ERROR("Decode PNEXT delta failed at line %u", lineNo);
        return -1;
    }

    /* Check if PNEXT is valid based on FLAG bits (same logic as compression) */
    auto flagIt = mappedFlag.find(lineNo);
    bool isPNextValid = false;
    if (flagIt != mappedFlag.end()) {
        uint16_t flag = flagIt->second;
        /* PNEXT is valid only when FLAG 0x1 is set and 0x8 is clear. */
        isPNextValid = ((flag & 0x1) != 0) && ((flag & 0x8) == 0);
    }

    int64_t pNext = 0;
    if (isPNextValid && (uint32_t)deltaLen > 1) {
        /* Valid PNEXT: decode as delta and reconstruct original value */
        int64_t pos = mappedPos.find(lineNo) == mappedPos.end() ? 0 : mappedPos[lineNo];
        std::string pNextDeltaStr = std::string((char*)deltaBuffer, (size_t)deltaLen - 1);
        try {
            pNext = (int64_t)std::stoll(pNextDeltaStr) + pos;
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to parse delta value '%s' for line %d: %s", pNextDeltaStr.c_str(), lineNo, e.what());
            return -1;
        }
    } else {
        /* Invalid PNEXT: decode as original value directly */
        if ((uint32_t)deltaLen > 1) {
            std::string pNextStr = std::string((char*)deltaBuffer, (size_t)deltaLen - 1);
            try {
                pNext = (int64_t)std::stoll(pNextStr);
            } catch (const std::exception& e) {
                LOG_ERROR("Failed to parse PNEXT value '%s' for line %d: %s", pNextStr.c_str(), lineNo, e.what());
                return -1;
            }
        }
    }
    
    nextMappedPos[lineNo] = pNext;
    char buff[32];
    int pNextLen = snprintf(buff, sizeof(buff), "%" PRId64 "\t", pNext);
    memcpy(outputBlock->getCurrent(), buff, pNextLen);
    outputBlock->setDataLen(outputBlock->getDataLen() + pNextLen);
    return pNextLen;
}

int32_t SamCodecActuator::decompressPosFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock) {
    /* Decode an unsigned LEB128 varint (mirrors compressPosFieldDelta); the accumulator resets on RNAME change. */
    if (lineNo > 0) {
        const auto curChrIt = mappedChr.find(lineNo);
        const auto prevChrIt = mappedChr.find(lineNo - 1);
        const uint16_t curChr = (curChrIt != mappedChr.end()) ? curChrIt->second : 0xFFFF;
        const uint16_t prevChr = (prevChrIt != mappedChr.end()) ? prevChrIt->second : 0xFFFF;
        if (curChr != prevChr) {
            posDeltaPrev = 0;
        }
    }
    uint64_t u = 0;
    int32_t shift = 0;
    int32_t deltaLen;
    do {
        uint8_t b = 0;
        deltaLen = fieldDecoders[fieldIdx]->decode_line(&b, 1, UINT8_MAX, false);
        if (deltaLen != 1) {
            LOG_ERROR("Decode POS delta failed at line %u", lineNo);
            return -1;
        }
        u |= (uint64_t)(b & 0x7f) << shift;
        if ((b & 0x80) == 0) break;
        shift += 7;
        if (shift >= 64) {
            LOG_ERROR("Decode POS delta overlong varint at line %u", lineNo);
            return -1;
        }
    } while (true);
    int64_t delta = (int64_t)u;

    posDeltaPrev += delta;
    int64_t pos = posDeltaPrev;
    mappedPos[lineNo] = pos;
    char buff[32];
    int posLen = snprintf(buff, sizeof(buff), "%" PRId64 "%c", pos, splitFlag);
    memcpy(outputBlock->getCurrent(), buff, posLen);
    outputBlock->setDataLen(outputBlock->getDataLen() + posLen);
    return posLen;
}

int32_t SamCodecActuator::decompressTLen(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock, const Json::Value& fieldMeta) {
    bool isOptimized = fieldMeta.isMember("optimized") && fieldMeta["optimized"].asBool();
    int32_t val;
    if (isOptimized) {
        auto cacheIt = tlenCache.find(lineNo);
        if (cacheIt != tlenCache.end()) {
            val = cacheIt->second;
        } else {
            bool minusOne = !fieldMeta.isMember("tlen_conv") || fieldMeta["tlen_conv"].asInt() != 0;
            val = computeTLEN(lineNo, minusOne);
        }
    } else {
        std::string mode = fieldMeta.isMember("mode") ? fieldMeta["mode"].asString() : "";
        if (mode == "string") {
            return decompressRegularField(fieldIdx, lineNo, splitFlag, outputBlock);
        }
        return decompressNumber<int32_t>(fieldIdx, lineNo, outputBlock);
    }

    char buffer[16];
    int32_t len = snprintf(buffer, sizeof(buffer), "%d", val);
    memcpy(outputBlock->getCurrent(), buffer, len);
    char* currentPos = (char*)outputBlock->getCurrent();
    currentPos[len] = splitFlag;
    outputBlock->setDataLen(outputBlock->getDataLen() + len + 1);
    return len + 1;
}

/* Pre-decode POS/CIGAR/PNEXT and cache per line, so TLEN has the full mate index. */
/* SEQ with no reference: one stream; build its decoder or decode it outright. */
int32_t SamCodecActuator::initSeqWholeBlockDecoder(const Json::Value& baseMeta, uint32_t idx,
                                                   RoughIOBlock* outputBlock)
{
    const std::string coderName = baseMeta["coder"]["magic"].asString();
    const uint32_t dstLength = baseMeta["totaldstlen"].asUInt();
    const uint32_t srcLength = baseMeta["totalsrclen"].asUInt();
    LOG_DEBUG("srclen = %d, dstlen = %d", srcLength, dstLength);

    /* Whether the stream is whole-block-only is the coder's own trait (CoderFactory::decoderIsWholeBlock). */
    if (CoderFactory::decoderIsWholeBlock(coderName)) {
        /* The whole-block SEQ is staged at the output block's tail, then moved to the head line by line. */
        coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstLength, &ioErrSink, "SEQ");
        baseIo.meta = baseMeta;
        baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
        coder_fc baseDecoder = coder_fc(&baseIo);
        if (baseDecoder.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcLength,
                                    srcLength, UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode SEQ by coder_fc failed, srclen = %u", srcLength);
            return -1;
        }
    } else {
        /* A record-by-record coder's stream is only prepared here; decompressBase pulls one record at a time. */
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "SEQ");
        io->meta = baseMeta;
        io->meta["dstlen"] = dstLength;

        FieldDecoderArgs args;
        if (baseMeta["coder"].isMember("level")) {
            args.level = baseMeta["coder"]["level"].asInt();
        }
        std::shared_ptr<coder> decoder = CoderFactory::makeFieldDecoder(coderName, io.get(), args);
        if (decoder == nullptr) {
            LOG_ERROR("Unsupported coder name:%s", coderName.c_str());
            return -1;
        }
        ioVector.push_back(io);
        fieldDecoders[idx] = decoder;
    }
    return (int32_t)dstLength;
}

/* The reference layout's match stream: under RLE, "m" (run lengths) and "mval" (values) expanded back to the sparse stream. */
int32_t SamCodecActuator::predecodeSeqMatchStream(const Json::Value& streams, uint32_t& streamId,
                                                  uint32_t firstDstLength, coder_io* matchIo,
                                                  uint32_t idx, uint32_t& extraOffset)
{
    const std::string matchCoderName = streams[streamId]["coder"]["magic"].asString();
    uint32_t runLen = streams[streamId]["srclen"].asUInt();
    uint32_t srcLen = streams[streamId]["orgrawlen"].asUInt();
    std::shared_ptr<coder> runDecoder =
        CoderFactory::makeFieldDecoder(matchCoderName, matchIo, FieldDecoderArgs());
    if (runDecoder == nullptr) {
        LOG_ERROR("check sub stream failed, coder name not match: %s", matchCoderName.c_str());
        return -1;
    }
    std::unique_ptr<uint8_t[]> runBuf = std::make_unique<uint8_t[]>(runLen + 1);
    if (runDecoder->decode_line(runBuf.get(), runLen, UINT8_MAX, false) < 0) {
        LOG_ERROR("Decode SEQ match run stream failed");
        return -1;
    }

    /* The value sub-stream, if any, sits after the run sub-stream; a fully matching block writes none. */
    uint32_t valDstLen = 0;
    uint32_t valLen = 0;
    std::unique_ptr<uint8_t[]> valBuf;
    if (streamId + 1 < streams.size() && streams[streamId + 1]["sname"].asString() == "mval") {
        streamId++;                       /* step past "m" onto the "mval" sub-stream */
        valDstLen = streams[streamId]["dstlen"].asUInt();
        valLen = streams[streamId]["srclen"].asUInt();
        std::shared_ptr<coder_io> valIo = makeCoderIo(
            inBlockPtr->getBuffer() + readOffset + firstDstLength, valDstLen, "SEQ match val");
        valIo->meta = streams[streamId];
        ioVector.push_back(valIo);
        /* Each sub-stream carries its own coder name; the two halves may differ. */
        const std::string valCoderName = streams[streamId]["coder"]["magic"].asString();
        std::shared_ptr<coder> valDecoder =
            CoderFactory::makeFieldDecoder(valCoderName, valIo.get(), FieldDecoderArgs());
        if (valDecoder == nullptr) {
            LOG_ERROR("check sub stream failed, coder name not match: %s", valCoderName.c_str());
            return -1;
        }
        valBuf = std::make_unique<uint8_t[]>(valLen + 1);
        if (valDecoder->decode_line(valBuf.get(), valLen, UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode SEQ match value stream failed");
            return -1;
        }
        extraOffset = valDstLen;   /* skipped past by the caller */
    } else {
        /* No surviving values in this fully-matching block. */
        extraOffset = 0;
    }

    /* Expand: each run is followed by one surviving value; calloc zero-fills the rest. */
    MemoryUtil::safeFree(matchBlockBuffer);
    matchBlockBuffer = MemoryUtil::safeAlloc<uint8_t>(srcLen + 1);
    if (matchBlockBuffer == nullptr) {
        LOG_ERROR("Alloc SEQ match block buffer failed");
        return -1;
    }
    uint32_t rp = 0, vp = 0, out = 0;
    for (uint32_t k = 0; k < valLen; k++) {
        uint32_t run = 0, shift = 0;
        while (rp < runLen) {
            uint8_t b = runBuf[rp++];
            run |= (uint32_t)(b & 0x7F) << shift;
            shift += 7;
            if ((b & 0x80) == 0) {
                break;
            }
        }
        out += run;
        if (out < srcLen) {
            matchBlockBuffer[out++] = valBuf[vp++];
        }
    }
    /* The trailing varint (trailing zeros) needs no action. */
    matchBlockLength = srcLen;
    matchBlockOffset = 0;
    matchBlockDecode = true;
    fieldDecoders[idx] = runDecoder;    /* kept for error reporting / reuse */
    return 0;
}

/* The reference layout's match stream when written whole (legacy, no "rle" member). */
int32_t SamCodecActuator::predecodeSeqWholeMatchStream(const Json::Value& streams, uint32_t streamId,
                                                       coder_io* matchIo, uint32_t idx)
{
    const std::string matchCoderName = streams[streamId]["coder"]["magic"].asString();
    /* The decoder comes from the stream's own magic; whether it is whole-block is the coder's trait. */
    fieldDecoders[idx] = CoderFactory::makeFieldDecoder(matchCoderName, matchIo, FieldDecoderArgs());
    if (fieldDecoders[idx] == nullptr) {
        LOG_ERROR("check sub stream failed, coder name not match: %s", matchCoderName.c_str());
        return -1;
    }
    matchBlockDecode = CoderFactory::decoderIsWholeBlock(matchCoderName);
    if (matchBlockDecode) {
        matchBlockLength = streams[streamId]["srclen"].asUInt();
        matchBlockOffset = 0;
        MemoryUtil::safeFree(matchBlockBuffer);
        matchBlockBuffer = MemoryUtil::safeAlloc<uint8_t>(matchBlockLength + 1);
        if (matchBlockBuffer == nullptr) {
            LOG_ERROR("Alloc SEQ match block buffer failed");
            return -1;
        }
        if (fieldDecoders[idx]->decode_line(matchBlockBuffer, matchBlockLength, UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode SEQ match stream (coder_fc) failed");
            return -1;
        }
    }
    return 0;
}

/* The baselen sub-stream: (ordinal, length) forward deltas, or absolute 4-byte pairs in older archives. */
int32_t SamCodecActuator::predecodeSeqBaseLengthStream(const Json::Value& streams, uint32_t streamId)
{
    if (streams[streamId]["sname"].asString() != "baselen") {
        LOG_ERROR("check sub stream failed. sname not match");
        return -1;
    }
    const uint32_t dstlen = streams[streamId]["dstlen"].asUInt();
    const uint32_t srclen = streams[streamId]["srclen"].asUInt();
    uint8_t* baseLenBuffer = MemoryUtil::safeAlloc<uint8_t>(srclen);

    {
        coder_io baseLenIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ length");
        FieldDecoderArgs args;
        if (streams[streamId]["coder"].isMember("level")) {
            args.level = streams[streamId]["coder"]["level"].asInt();
        }
        std::shared_ptr<coder> baseLenCoder =
            CoderFactory::makeFieldDecoder(streams[streamId]["coder"]["magic"].asString(),
                                           &baseLenIo, args);
        if (baseLenCoder == nullptr ||
            baseLenCoder->decode_line((uint8_t*)baseLenBuffer, srclen, UINT8_MAX, false) < 0) {
            MemoryUtil::safeFree(baseLenBuffer);
            LOG_ERROR("Decode base lengths failed (coder %s)",
                      streams[streamId]["coder"]["magic"].asString().c_str());
            return -1;
        }
    }

    if (streams[streamId]["delta"].isUInt()) {
        const uint32_t count = streams[streamId]["count"].asUInt();
        std::vector<uint32_t> ordinals(count, 0);
        uint32_t off = 0;
        uint32_t value = 0;
        uint32_t ordinal = 0;
        for (uint32_t i = 0; i < count; ++i) {
            if (!readVarint(baseLenBuffer, srclen, off, value)) {
                MemoryUtil::safeFree(baseLenBuffer);
                LOG_ERROR("Decode base lengths failed: malformed ordinal");
                return -1;
            }
            ordinal += value;
            ordinals[i] = ordinal;
        }
        for (uint32_t i = 0; i < count; ++i) {
            if (!readVarint(baseLenBuffer, srclen, off, value)) {
                MemoryUtil::safeFree(baseLenBuffer);
                LOG_ERROR("Decode base lengths failed: malformed length");
                return -1;
            }
            unmapedReadLength.push_back(std::make_pair(ordinals[i], value));
        }
    } else {
        const uint32_t* baseLenPtr = (const uint32_t*)baseLenBuffer;
        const uint32_t baseLenCount = srclen >> 2;
        for (uint32_t i = 0; i + 1 < baseLenCount; i += 2) {
            unmapedReadLength.push_back(std::make_pair(baseLenPtr[i], baseLenPtr[i + 1]));
        }
    }
    readOffset += dstlen;
    MemoryUtil::safeFree(baseLenBuffer);
    return 0;
}

int32_t SamCodecActuator::predecodeSeqSelfRefStream(const Json::Value& streams, uint32_t streamId)
{
    if (streams[streamId]["sname"].asString() != "selfref") {
        LOG_ERROR("check sub stream failed. sname not match: %s",
                  streams[streamId]["sname"].asString().c_str());
        return -1;
    }
    const uint32_t srclen = streams[streamId]["srclen"].asUInt();
    const uint32_t dstlen = streams[streamId]["dstlen"].asUInt();
    if (srclen == 0) {
        LOG_ERROR("SEQ self reference stream is empty");
        return -1;
    }

    std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstlen,
                                               "SEQ selfref");
    io->meta = streams[streamId];
    std::shared_ptr<coder> decoder =
        CoderFactory::makeFieldDecoder(streams[streamId]["coder"]["magic"].asString(), io.get(),
                                       FieldDecoderArgs());
    if (decoder == nullptr) {
        LOG_ERROR("check sub stream failed, coder name not match: %s",
                  streams[streamId]["coder"]["magic"].asString().c_str());
        return -1;
    }
    baseSelfRefSquash.assign(srclen, 0);
    if (decoder->decode_line(baseSelfRefSquash.data(), srclen, UINT8_MAX, false) <= 0) {
        LOG_ERROR("Decode SEQ self reference failed");
        return -1;
    }
    readOffset += dstlen;
    return 0;
}

/* SEQ exception streams: one per character present (stream meta "ch"); the older layout's lone "npos" list is also read. */
int32_t SamCodecActuator::predecodeSeqExceptionStreams(const Json::Value& baseMeta,
                                                       const Json::Value& streams,
                                                       uint32_t& streamId)
{
    seqExc.clear();
    for (;;) {
        const std::string name = (streams.size() > (size_t)streamId + 1)
                                     ? streams[streamId + 1]["sname"].asString()
                                     : std::string();
        /* Exception streams come next; stop at the first other name. */
        if (seqExcFormOf(name) == SeqExcForm::None) {
            break;   /* this block has no further exception stream */
        }
        streamId++;
        const SeqExcLayout layout = seqExcLayoutOf(baseMeta, streams[streamId]);
        const uint8_t ch = layout.ch;
        const uint32_t dstlen = streams[streamId]["dstlen"].asUInt();
        const uint32_t srclen = streams[streamId]["srclen"].asUInt();
        coder_io excIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ npos");
        FieldDecoderArgs args;
        if (streams[streamId]["coder"].isMember("level")) {
            args.level = streams[streamId]["coder"]["level"].asInt();
        }
        std::shared_ptr<coder> excCoder =
            CoderFactory::makeFieldDecoder(streams[streamId]["coder"]["magic"].asString(),
                                           &excIo, args);
        /* The payload carries positions alone; the character is in the stream meta's "ch". */
        std::vector<uint8_t> raw(srclen);
        std::vector<uint32_t> pos;
        if (excCoder == nullptr ||
            excCoder->decode_line(raw.data(), srclen, UINT8_MAX, false) < 0 ||
            !seqExcDecodePositions(layout, raw.data(), srclen, pos)) {
            LOG_ERROR("Decode SEQ exception positions failed: %u bytes, %u positions, form %d",
                      srclen, layout.count, (int)layout.form);
            return -1;
        }
        readOffset += dstlen;

        seqExc.push_back(SeqExcStream());
        seqExc.back().byte = ch;
        seqExc.back().pos = std::move(pos);
        LOG_DEBUG("SEQ exc read ch=%u '%c' positions=%u first=%u", (unsigned)ch, (char)ch,
                  (uint32_t)seqExc.back().pos.size(),
                  seqExc.back().pos.empty() ? 0u : seqExc.back().pos[0]);
    }
    return 0;
}

/* The QUAL column: decode its frequency table first, then build the value decoder from the stream's magic. */
int32_t SamCodecActuator::initQualFieldDecoders(const Json::Value& qualMeta)
{
    const Json::Value& qualStreamMeta = qualMeta["streams"];
    if (qualStreamMeta.size() != 2) {
        LOG_ERROR("quality streams check failded, size = %d", qualStreamMeta.size());
        return -1;
    }
    const uint32_t qualDstLength = qualStreamMeta[0]["dstlen"].asUInt();
    const uint32_t freqDstLength = qualStreamMeta[1]["dstlen"].asUInt();
    coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualDstLength, freqDstLength, &ioErrSink, "QUAL freq table");
    FieldDecoderArgs freqArgs;
    if (qualStreamMeta[1]["coder"].isMember("level")) {
        freqArgs.level = qualStreamMeta[1]["coder"]["level"].asInt();
    }
    std::shared_ptr<coder> qualFreqCoder =
        CoderFactory::makeFieldDecoder(qualStreamMeta[1]["coder"]["magic"].asString(),
                                       &qualFreqIo, freqArgs);
    const uint32_t qualFreqSrcLength = qualStreamMeta[1]["srclen"].asUInt();
    /* Counting with uint8_t wraps when the alphabet exceeds 127 symbols. */
    const uint32_t qualFreqArrLength = qualFreqSrcLength / sizeof(uint16_t);
    uint16_t* qualFreqArr = new uint16_t[qualFreqArrLength];
    const int32_t qualFreq = (qualFreqCoder == nullptr)
                                 ? -1
                                 : qualFreqCoder->decode_line((uint8_t*)qualFreqArr,
                                                             qualFreqSrcLength, UINT8_MAX, false);
    if (qualFreq < 0 || (uint32_t)qualFreq != qualFreqSrcLength) {
        LOG_ERROR("Decode quality frequncy failed");
        delete [] qualFreqArr;
        return -1;
    }
    for (uint32_t i = 0; i < qualFreqArrLength; i += 2) {
        qualFreqTable.push_back(std::make_pair(qualFreqArr[i], qualFreqArr[i + 1]));
    }
    delete [] qualFreqArr;

    std::shared_ptr<coder_io> qualIo = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualDstLength, "QUAL");
    ioVector.push_back(qualIo);
    /* The magic-to-decoder step is CoderFactory::makeQualDecoder's; this side supplies the alphabet and the prior. */
    QualDecoderArgs qualArgs;
    qualArgs.freqTable = qualFreqTable.empty() ? nullptr : &qualFreqTable;
    qualArgs.priorRequired = qualStreamMeta[0].isMember("prior");
    qualArgs.priorAddress = qualArgs.priorRequired ? qualStreamMeta[0]["prior"].asInt64() : -1;
    AuxPayloadPtr qualPriorBlob =
        (qualArgs.priorRequired && pbgzEngine != nullptr)
            ? pbgzEngine->getQualPrior(inBlockPtr->getPackageIndex())
            : AuxPayloadPtr();
    qualArgs.priorBlob = qualPriorBlob.get();
    qualDecoder = CoderFactory::makeQualDecoder(
        qualStreamMeta[0]["coder"]["magic"].asString(), qualIo.get(), qualArgs);
    if (qualDecoder == nullptr) {
        LOG_ERROR("Unsupport coder type: %s",
                  qualStreamMeta[0]["coder"]["magic"].asString().c_str());
        return -1;
    }
    readOffset += (qualDstLength + freqDstLength);
    return 0;
}

/* Any ordinary column: one stream, decoder built from its magic; its position is recorded for preDecodeForTLEN. */
int32_t SamCodecActuator::initFieldDecoder(uint32_t fieldIdx, const Json::Value& fieldMeta)
{
    const std::string coderName = fieldMeta["coder"]["magic"].asString();
    const uint32_t dstLen = fieldMeta["dstlen"].asUInt();
    fieldIoStart[fieldIdx] = readOffset;
    fieldIoDstLen[fieldIdx] = dstLen;
    if (fieldMeta["coder"].isMember("level")) {
        fieldIoLevel[fieldIdx] = fieldMeta["coder"]["level"].asInt();
    }

    std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "SAM field");
    ioVector.push_back(io);

    /* The level is handed over only for coders whose stream records one. */
    FieldDecoderArgs args;
    if (coderName == "coder_affix_match" ||
        (coderName == "coder_arith" && fieldMeta["coder"].isMember("level"))) {
        args.level = fieldMeta["coder"]["level"].asInt();
    }
    AuxPayloadPtr posPrior;
    if (coderName == "coder_arith") {
        posPrior = (pbgzEngine != nullptr) ? pbgzEngine->getPosPrior() : AuxPayloadPtr();
        args.posPrior = posPrior.get();
    }

    fieldDecoders[fieldIdx] = CoderFactory::makeFieldDecoder(coderName, io.get(), args);
    if (fieldDecoders[fieldIdx] == nullptr) {
        LOG_ERROR("Unsupport coder type: %s", coderName.c_str());
        return -1;
    }
    readOffset += dstLen;
    return 0;
}

/* The TLEN exception stream: records the pair walk cannot produce; its layout is read from the field meta. */
int32_t SamCodecActuator::decodeTlenExceptionStream(const Json::Value& stream,
                                                    const Json::Value& tlenMeta)
{
    const uint32_t srclen = stream["srclen"].asUInt();
    const uint32_t dstlen = stream["dstlen"].asUInt();
    uint8_t* excBuffer = MemoryUtil::safeAlloc<uint8_t>(srclen);
    if (excBuffer == nullptr) {
        return -1;
    }
    coder_io tlenIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "TLEN exceptions");
    FieldDecoderArgs args;
    if (stream["coder"].isMember("level")) {
        args.level = stream["coder"]["level"].asInt();
    }
    std::shared_ptr<coder> tlenDecoder =
        CoderFactory::makeFieldDecoder(stream["coder"]["magic"].asString(), &tlenIo, args);
    if (tlenDecoder == nullptr ||
        tlenDecoder->decode_line(excBuffer, srclen, UINT8_MAX, false) < 0) {
        MemoryUtil::safeFree(excBuffer);
        LOG_ERROR("Decode TLEN exceptions failed (coder %s)",
                  stream["coder"]["magic"].asString().c_str());
        return -1;
    }

    const uint32_t declared = tlenMeta["exceptions"].asUInt();
    const Json::Value& layout = tlenMeta["exc"];
    std::map<uint32_t, int32_t> decoded;
    bool decodedOk;
    if (layout.isUInt()) {
        decodedOk = (layout.asUInt() == (Json::Value::UInt)TLEN_EXC_LAYOUT_VARINT)
                        ? tlenDecodeVarints(excBuffer, srclen, declared, decoded)
                        : tlenDecodeFixedPairs(excBuffer, srclen, declared, decoded);
    } else {
        decodedOk = tlenDecodeFixedPairs(excBuffer, srclen, declared, decoded) ||
                    tlenDecodeVarints(excBuffer, srclen, declared, decoded);
    }
    if (!decodedOk) {
        MemoryUtil::safeFree(excBuffer);
        LOG_ERROR("Corrupt TLEN exception stream: neither the varint layout nor the "
                  "fixed-pair layout fits (%u bytes, %u exceptions)", srclen, declared);
        return -1;
    }
    tlenCache.insert(decoded.begin(), decoded.end());
    MemoryUtil::safeFree(excBuffer);
    readOffset += dstlen;
    return 0;
}

/* The TLEN column: the exception stream under reconstruction, else an ordinary column. */
int32_t SamCodecActuator::initTlenFieldDecoders(const Json::Value& tlenMeta, uint32_t fieldIdx)
{
    if (tlenMeta.isMember("optimized") && tlenMeta["optimized"].asBool()) {
        if (tlenMeta.isMember("streams") && tlenMeta["streams"].isArray()) {
            const Json::Value& tlenStreams = tlenMeta["streams"];
            for (Json::Value::const_iterator it = tlenStreams.begin(); it != tlenStreams.end(); ++it) {
                const Json::Value& stream = *it;
                if (stream["sname"].asString() != "tlenexc") {
                    continue;
                }
                if (decodeTlenExceptionStream(stream, tlenMeta) != 0) {
                    return -1;
                }
            }
        }
        return 0;
    }

    const std::string coderName = tlenMeta["coder"]["magic"].asString();
    const uint32_t dstLen = tlenMeta["dstlen"].asUInt();
    std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "TLEN");
    ioVector.push_back(io);
    FieldDecoderArgs args;
    if (tlenMeta["coder"].isMember("level")) {
        args.level = tlenMeta["coder"]["level"].asInt();
    }
    fieldDecoders[fieldIdx] = CoderFactory::makeFieldDecoder(coderName, io.get(), args);
    if (fieldDecoders[fieldIdx] == nullptr) {
        LOG_ERROR("Unsupported TLEN coder type: %s", coderName.c_str());
        return -1;
    }
    readOffset += dstLen;
    return 0;
}

/* The QNAME column's own decoders, rebuilt from the recorded offsets so the main loop's are not consumed. */
int32_t SamCodecActuator::predecodeIdStreams(IdPredecodeState& state)
{
    Json::Value& streams = meta["sam"]["streams"];
    uint8_t* buffer = inBlockPtr->getBuffer();

    state.streams.reserve(idStreamOffsets.size());
    state.decoders.reserve(idStreamOffsets.size());
    for (uint32_t si = 0; si < idStreamOffsets.size(); ++si) {
        std::shared_ptr<coder_io> io = makeCoderIo(buffer + idStreamOffsets[si], idStreamDstLens[si], "QNAME predecode");
        state.streams.push_back(io);
        int32_t lvl = -1;
        if (streams[0].isMember("streams") && streams[0]["streams"].isValidIndex(si) &&
            streams[0]["streams"][si]["coder"].isMember("level")) {
            lvl = streams[0]["streams"][si]["coder"]["level"].asInt();
        }
        FieldDecoderArgs idPreArgs;
        if (lvl >= 0) {
            idPreArgs.level = lvl;
        }
        state.decoders.push_back(
            CoderFactory::makeFieldDecoder(idStreamCoders[si], io.get(), idPreArgs));
    }

    state.numericBufs.assign(idStreamOffsets.size(), nullptr);
    state.numericLens.assign(idStreamOffsets.size(), 0);
    state.numericPos.assign(idStreamOffsets.size(), 0);
    state.numericAcc.assign(idStreamOffsets.size(), 0);
    for (uint32_t si = 0; si < idStreamOffsets.size(); ++si) {
        if (si >= state.decoders.size()) {
            continue;
        }
        if (!(streams[0].isMember("streams") && streams[0]["streams"].isValidIndex(si) &&
              streams[0]["streams"][si].isMember("mode"))) {
            continue;
        }
        const std::string preMode = streams[0]["streams"][si]["mode"].asString();
        const uint32_t srclen = streams[0]["streams"][si]["srclen"].asUInt();
        if (preMode == "cnt") {
            uint8_t* cb = MemoryUtil::safeAlloc<uint8_t>(srclen + 8);
            if (cb == nullptr) {
                return -1;
            }
            const int32_t cl = expandCounterSubStream(buffer + idStreamOffsets[si],
                                                      idStreamDstLens[si],
                                                      streams[0]["streams"][si], cb, srclen + 8);
            if (cl < 0) {
                MemoryUtil::safeFree(cb);
                LOG_ERROR("Predecode id counter sub-stream(%u) failed", si);
                return -1;
            }
            state.numericBufs[si] = cb;
            state.numericLens[si] = (uint32_t)cl;
            continue;
        }
        if (!state.decoders[si] || preMode != "numeric") {
            continue;
        }
        uint8_t* b = MemoryUtil::safeAlloc<uint8_t>(srclen + 8);
        if (b == nullptr) {
            return -1;
        }
        int32_t dl = state.decoders[si]->decode_line(b, srclen, UINT8_MAX, false);
        if (dl < 0) {
            MemoryUtil::safeFree(b);
            LOG_ERROR("Predecode id numeric sub-stream(%u) failed: %d", si, dl);
            return -1;
        }
        state.numericBufs[si] = b;
        state.numericLens[si] = (uint32_t)dl;
    }
    return 0;
}

/* Pre-decode FLAG (field 1): PNEXT's delta form depends on the flag bits (mirrors compressPNextFieldDelta). */
int32_t SamCodecActuator::predecodeFlagColumn()
{
    Json::Value& streams = meta["sam"]["streams"];
    auto flagStartIt = fieldIoStart.find(1);
    if (flagStartIt == fieldIoStart.end() || !streams.isValidIndex(1) ||
        !streams[1].isMember("coder")) {
        return 0;
    }

    const uint32_t off = flagStartIt->second;
    const uint32_t dstlen = fieldIoDstLen[1];
    const std::string coderName = streams[1]["coder"]["magic"].asString();
    const std::string mode = streams[1].isMember("mode") ? streams[1]["mode"].asString() : "";

    uint8_t tmpBuf[1024];
    std::shared_ptr<coder_io> tmpIo = makeCoderIo(inBlockPtr->getBuffer() + off, dstlen, "FLAG predecode");
    FieldDecoderArgs flagArgs;
    auto flagLvIt = fieldIoLevel.find(1);
    if (flagLvIt != fieldIoLevel.end()) {
        flagArgs.level = flagLvIt->second;
    }
    std::shared_ptr<coder> tmpDec = CoderFactory::makeFieldDecoder(coderName, tmpIo.get(), flagArgs);
    if (!tmpDec) {
        return 0;
    }

    const bool binary = (mode == "number" || mode == "");
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        const bool need2hold = CoderFactory::decoderHoldsCallerBuffer(coderName);
        int32_t len;
        if (binary) {
            len = tmpDec->decode_line(tmpBuf, sizeof(uint16_t), UINT8_MAX, need2hold);
            mappedFlag[lineNo] = (len >= (int32_t)sizeof(uint16_t)) ? *(uint16_t*)tmpBuf : 0;
        } else {
            len = tmpDec->decode_line(tmpBuf, (uint32_t)sizeof(tmpBuf), '\t', need2hold);
            if (len > 1) {
                const std::string flagStr((char*)tmpBuf, (size_t)len - 1);
                try {
                    mappedFlag[lineNo] = (uint16_t)std::stoll(flagStr);
                } catch (...) {
                    mappedFlag[lineNo] = 0;
                }
            } else {
                mappedFlag[lineNo] = 0;
            }
        }
    }
    return 0;
}

/* Pre-decode RNAME (field 2) so the POS delta chain can detect chromosome switches. */
int32_t SamCodecActuator::predecodeRnameColumn()
{
    Json::Value& streams = meta["sam"]["streams"];
    auto rnameStartIt = fieldIoStart.find(2);
    if (rnameStartIt == fieldIoStart.end() || !streams.isValidIndex(2) ||
        !streams[2].isMember("coder")) {
        return 0;
    }

    const uint32_t rnOff = rnameStartIt->second;
    const uint32_t rnDst = fieldIoDstLen[2];
    const std::string rnCoderName = streams[2]["coder"]["magic"].asString();

    std::shared_ptr<coder_io> rnIo = makeCoderIo(inBlockPtr->getBuffer() + rnOff, rnDst, "RNAME predecode");
    FieldDecoderArgs rnArgs;
    auto rnLvIt = fieldIoLevel.find(2);
    if (rnLvIt != fieldIoLevel.end()) {
        rnArgs.level = rnLvIt->second;
    }
    std::shared_ptr<coder> rnDec = CoderFactory::makeFieldDecoder(rnCoderName, rnIo.get(), rnArgs);
    if (!rnDec) {
        return 0;
    }

    const bool rnHold = CoderFactory::decoderHoldsCallerBuffer(rnCoderName);
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        uint16_t chrIndex = 0;
        const int32_t len = rnDec->decode_line((uint8_t*)&chrIndex, sizeof(chrIndex), UINT8_MAX, rnHold);
        mappedChr[lineNo] = (len >= (int32_t)sizeof(chrIndex)) ? chrIndex : 0xFFFF;
    }
    return 0;
}

/* Build the pre-decoder for a TLEN-reconstruction column from the stream's magic. Returns false when there is nothing to pre-decode. */
bool SamCodecActuator::buildTlenColumnPredecoder(uint32_t fieldIdx, TlenColumnPredecoder& out)
{
    Json::Value& streams = meta["sam"]["streams"];
    auto startIt = fieldIoStart.find(fieldIdx);
    if (startIt == fieldIoStart.end() || !streams.isValidIndex(fieldIdx) ||
        !streams[fieldIdx].isMember("coder")) {
        return false;
    }

    const uint32_t off = startIt->second;
    const uint32_t dstlen = fieldIoDstLen[fieldIdx];
    out.coderName = streams[fieldIdx]["coder"]["magic"].asString();
    out.mode = streams[fieldIdx].isMember("mode") ? streams[fieldIdx]["mode"].asString() : "";
    /* Column source length: bounds the widest line (sizes the decode buffer in rebuildTlenColumn). */
    out.srcLen = streams[fieldIdx].isMember("srclen") ? streams[fieldIdx]["srclen"].asUInt() : 0;

    out.stream = makeCoderIo(inBlockPtr->getBuffer() + off, dstlen, "TLEN predecode");
    FieldDecoderArgs args;
    auto lvIt = fieldIoLevel.find(fieldIdx);
    if (lvIt != fieldIoLevel.end()) {
        args.level = lvIt->second;
    }
    AuxPayloadPtr posPrior;
    if (out.coderName == "coder_arith") {
        posPrior = (pbgzEngine != nullptr) ? pbgzEngine->getPosPrior() : AuxPayloadPtr();
        args.posPrior = posPrior.get();
    }
    out.decoder = CoderFactory::makeFieldDecoder(out.coderName, out.stream.get(), args);
    return out.decoder != nullptr;
}

/* Rebuild one TLEN column from its pre-decoder, line by line, into tlenPreDecodedFields and the per-column maps. */
int32_t SamCodecActuator::rebuildTlenColumn(uint32_t fieldIdx, const TlenColumnPredecoder& pre)
{
    const std::string& coderName = pre.coderName;
    const std::string& mode = pre.mode;
    std::vector<std::string>& fieldCache = tlenPreDecodedFields[fieldIdx];
    fieldCache.clear();
    fieldCache.reserve(samLine);

    /* Per-line decode buffer, sized from the column's source length; left uninitialized (a
       vector would zero-fill it, and the buffer is overwritten by every decode_line). */
    const uint32_t decodeCap = std::max<uint32_t>(1024, pre.srcLen + 1);
    std::unique_ptr<uint8_t[]> tmpStore(new uint8_t[decodeCap]);
    uint8_t* tmpBuf = tmpStore.get();
    /* Delta chain of pos_delta mode: the absolute POS reconstructed from the previous line. */
    int64_t posPrev = 0;
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        /* coder_affix_match points at the caller's buffer unless need2hold is set; set it here or the context is overwritten. */
        const bool need2hold = CoderFactory::decoderHoldsCallerBuffer(coderName);
        /* CIGAR is always textual; POS/PNEXT are binary unless mode is textual. POS uses a zigzag varint read byte by byte. */
        if (fieldIdx == 3 && mode == "pos_delta") {
            /* Unsigned varint decode, mirroring compressPosFieldDelta. */
            if (lineNo > 0) {
                const auto curChrIt = mappedChr.find(lineNo);
                const auto prevChrIt = mappedChr.find(lineNo - 1);
                const uint16_t curChr = (curChrIt != mappedChr.end()) ? curChrIt->second : 0xFFFF;
                const uint16_t prevChr = (prevChrIt != mappedChr.end()) ? prevChrIt->second : 0xFFFF;
                if (curChr != prevChr) {
                    posPrev = 0;
                }
            }
            uint64_t u = 0;
            int32_t shift = 0;
            int32_t vlen;
            do {
                uint8_t b = 0;
                vlen = pre.decoder->decode_line(&b, 1, UINT8_MAX, need2hold);
                if (vlen != 1) {
                    LOG_ERROR("Decode POS delta failed at line %u", lineNo);
                    return -1;
                }
                u |= (uint64_t)(b & 0x7f) << shift;
                if ((b & 0x80) == 0) break;
                shift += 7;
                if (shift >= 64) {
                    LOG_ERROR("Decode POS delta overlong varint at line %u", lineNo);
                    return -1;
                }
            } while (true);
            int64_t delta = (int64_t)u;
            int64_t pos = posPrev + delta;
            posPrev = pos;
            mappedPos[lineNo] = pos;
            fieldCache.emplace_back(std::to_string(pos) + '\t');
            continue;
        }

        bool binary = (fieldIdx != 5) && (mode == "number" || mode == "");
        int32_t len;
        if (binary) {
            len = pre.decoder->decode_line(tmpBuf, 4, UINT8_MAX, need2hold);
        } else {
            len = pre.decoder->decode_line(tmpBuf, decodeCap, '\t', need2hold);
        }
        if (len <= 1) {
            fieldCache.emplace_back();
            continue;
        }

        switch (fieldIdx) {
            case 3:
                if (binary) {
                    mappedPos[lineNo] = (int64_t)(uint32_t)(*(uint32_t*)tmpBuf);
                    fieldCache.emplace_back(std::to_string(mappedPos[lineNo]) + '\t');
                } else {
                    mappedPos[lineNo] = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                    fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                }
                break;
            case 5:
                cigarReadLen[lineNo] = cigarRefConsumed(tmpBuf, len);
                baseLengthBuffer[lineNo] = cigarSeqConsumed(tmpBuf, len);
                if (cigarOpList.size() <= lineNo) cigarOpList.resize(lineNo + 1);
                parseCigarOps(tmpBuf, (uint32_t)len, cigarOpList[lineNo]);
                fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                break;
            case 7:
                if (binary) {
                    int64_t pnextBin = (int64_t)(uint32_t)(*(uint32_t*)tmpBuf);
                    nextMappedPos[lineNo] = pnextBin;
                    fieldCache.emplace_back(std::to_string(pnextBin) + '\t');
                } else if (mode == "pnext_delta") {
                    /* Check if PNEXT is valid based on FLAG bits (same logic as compression) */
                    auto flagIt = mappedFlag.find(lineNo);
                    bool isPNextValid = false;
                    if (flagIt != mappedFlag.end()) {
                        uint16_t flag = flagIt->second;
                        /* PNEXT is valid only when FLAG 0x1 is set and 0x8 is clear. */
                        isPNextValid = ((flag & 0x1) != 0) && ((flag & 0x8) == 0);
                    }

                    if (isPNextValid) {
                        /* Valid PNEXT: stored as delta against POS, reconstruct */
                        int64_t delta = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                        int64_t pos = mappedPos.count(lineNo) ? mappedPos[lineNo] : 0;
                        int64_t pnext = delta + pos;
                        nextMappedPos[lineNo] = pnext;
                        fieldCache.emplace_back(std::to_string(pnext) + '\t');
                    } else {
                        /* Invalid PNEXT: stored as original value directly */
                        if ((uint32_t)len > 1) {
                            nextMappedPos[lineNo] = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                            fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                        } else {
                            nextMappedPos[lineNo] = 0;
                            fieldCache.emplace_back("0\t");
                        }
                    }
                } else {
                    nextMappedPos[lineNo] = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                    fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                }
                break;
            default:
                fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                break;
        }
    }
    return 0;
}

/* Pre-decode QNAME (field 0) into decodedQnames, the pairing key for pnext_qname_rebuild. */
int32_t SamCodecActuator::rebuildQnameColumn()
{
    decodedQnames.clear();
    decodedQnames.resize(samLine);
    if (idStreamOffsets.empty()) {
        return 0;
    }

    IdPredecodeState preState;
    if (predecodeIdStreams(preState) != 0) {
        return -1;
    }
    /* The rebuild loop below speaks in terms of the state's vectors, as it always has. */
    auto& preIdDec = preState.decoders;
    auto& numBufs = preState.numericBufs;
    auto& numLens = preState.numericLens;
    auto& numPos = preState.numericPos;
    auto& numAcc = preState.numericAcc;

    char qnameBuf[512];
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        std::string qname;
        if (idUsesQnameCoder && !preIdDec.empty()) {
            // Single-stream coder_qname: one QNAME per line (with trailing '\t').
            bool hold = CoderFactory::decoderHoldsCallerBuffer(idStreamCoders[0]);
            int32_t len = preIdDec[0]->decode_line((uint8_t*)qnameBuf, sizeof(qnameBuf), UINT8_MAX, hold);
            if (len > 0) {
                int32_t strip = (qnameBuf[len - 1] == '\t') ? 1 : 0;
                qname.assign(qnameBuf, (size_t)len - strip);
            }
        } else {
            // Reconstruct QNAME from split sub-streams, appending each decoded segment (mirrors decompressIdField).
            for (uint32_t si = 0; si < preIdDec.size(); ++si) {
                if (!preIdDec[si] && numBufs[si] == nullptr) continue;
                if (numBufs[si] != nullptr) {
                    /* Numeric sub-stream: one varint per line, no separator byte. */
                    if (numPos[si] >= numLens[si]) {
                        break;
                    }
                    uint32_t zz = 0;
                    if (!readVarint(numBufs[si], numLens[si], numPos[si], zz)) {
                        LOG_ERROR("Corrupt id numeric sub-stream(%u) in predecode", si);
                        return -1;
                    }
                    const uint32_t delta = zz >> 1;
                    numAcc[si] = ((zz & 1u) == 0) ? (numAcc[si] + delta) : (numAcc[si] - delta);
                    char ntmp[24];
                    int nn = snprintf(ntmp, sizeof(ntmp), "%llu",
                                      (unsigned long long)numAcc[si]);
                    if (nn > 0) {
                        qname.append(ntmp, (size_t)nn);
                    }
                    continue;
                }
                if (!idConstTexts.empty() && !idConstTexts[si].empty()) {
                    /* Constant segment: the same text every line, stored once in the block. */
                    qname.append(idConstTexts[si]);
                    continue;
                }
                if (!idIntActive.empty() && idIntActive[si]) {
                    /* Value-domain segment: one value per line, no separator. */
                    char ntmp[24];
                    const uint64_t v = idIntModes[si].decode(idIntCoders[si]);
                    const int nn = snprintf(ntmp, sizeof(ntmp), "%llu", (unsigned long long)v);
                    if (nn > 0) {
                        qname.append(ntmp, (size_t)nn);
                    }
                    continue;
                }
                if (!preIdDec[si]) {
                    /* Every coder-backed layout is handled above; an unknown mode must not be dereferenced. */
                    LOG_ERROR("Predecode id sub-stream(%u) has no decoder", si);
                    return -1;
                }
                uint8_t sep = (si < idAnalysis.symbols.size()) ? idAnalysis.symbols[si] : UINT8_MAX;
                // coder_affix_match keeps cross-line context only when need2hold is set;
                // without it the next line would corrupt the prefix reference.
                bool hold = CoderFactory::decoderHoldsCallerBuffer(idStreamCoders[si]);
                int32_t len = preIdDec[si]->decode_line((uint8_t*)qnameBuf, sizeof(qnameBuf), sep, hold);
                if (len > 0) {
                    int32_t strip = (qnameBuf[len - 1] == '\t') ? 1 : 0;
                    qname.append(qnameBuf, (size_t)len - strip);
                }
            }
        }
        decodedQnames[lineNo] = qname;
    }
    for (uint32_t si = 0; si < numBufs.size(); ++si) {
        MemoryUtil::safeFree(numBufs[si]);
    }
    return 0;
}

/* PNEXT (field 7) in pnext_qname_rebuild mode: exceptions from streams[7]["streams"][0], the rest rebuilt by QNAME pairing. */
int32_t SamCodecActuator::rebuildPnextByQname()
{
    Json::Value& streams = meta["sam"]["streams"];
    std::unordered_map<uint32_t, int64_t> pnextCache;   /* line -> value, for the exception pairs */
    pnextCache.clear();
    // Decode the exception stream (pairs of contentIdx, delta).
    if (streams[7].isMember("streams") && streams[7]["streams"].size() > 0) {
        Json::Value& excMeta = streams[7]["streams"][0];
        uint32_t excOff = fieldIoStart[7];
        uint32_t excDstLen = excMeta["dstlen"].asUInt();
        std::string excCoderName = excMeta["coder"]["magic"].asString();
        std::shared_ptr<coder_io> excIo = makeCoderIo(inBlockPtr->getBuffer() + excOff, excDstLen, "PNEXT exc predecode");
        std::shared_ptr<coder> excDec =
            CoderFactory::makeFieldDecoder(excCoderName, excIo.get(), FieldDecoderArgs());
        if (excDec) {
            // Exception coder level: prefer the value stored in the sub-stream meta.
            if (excMeta["coder"].isMember("level")) {
                excDec->set_level(excMeta["coder"]["level"].asInt());
            } else {
                auto lvIt = fieldIoLevel.find(7);
                if (lvIt != fieldIoLevel.end()) excDec->set_level(lvIt->second);
            }
            const int32_t pairCount = streams[7]["exceptions"].asUInt();
            const std::string excEnc = streams[7].isMember("exc_enc")
                                           ? streams[7]["exc_enc"].asString()
                                           : std::string();
            /* One buffer size for both varint layouts: a u32 varint (<=5 B) plus a zigzag i64 varint (<=10 B). */
            std::vector<uint8_t> raw((size_t)pairCount * 16, 0);
            size_t rawOff = 0;
            auto nextVarint = [&raw, &rawOff](uint64_t& out) -> bool {
                out = 0;
                int shift = 0;
                while (rawOff < raw.size()) {
                    const uint8_t b = raw[rawOff++];
                    out |= (uint64_t)(b & 0x7f) << shift;
                    if ((b & 0x80) == 0) return true;
                    shift += 7;
                    if (shift >= 64) return false;
                }
                return false;
            };
            if (excEnc == "varint2") {
                /* Ordinals as forward deltas, then deltas as zigzag varints; a constant delta lives in the meta. */
                excDec->decode_line(raw.data(), (uint32_t)raw.size(), UINT8_MAX, false);
                const bool constDelta = streams[7].isMember("exc_delta");
                const uint64_t constZ = constDelta ? streams[7]["exc_delta"].asUInt64() : 0;
                std::vector<uint32_t> ordinals((size_t)pairCount, 0);
                uint32_t ordinal = 0;
                for (int32_t i = 0; i < pairCount; ++i) {
                    uint64_t step = 0;
                    if (!nextVarint(step)) {
                        LOG_ERROR("PNEXT exception varint stream truncated");
                        return -1;
                    }
                    ordinal += (uint32_t)step;
                    ordinals[i] = ordinal;
                }
                for (int32_t i = 0; i < pairCount; ++i) {
                    uint64_t z = constZ;
                    if (!constDelta && !nextVarint(z)) {
                        LOG_ERROR("PNEXT exception varint stream truncated");
                        return -1;
                    }
                    const int64_t delta = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); /* de-zigzag */
                    const uint32_t lineNo = ordinals[i];
                    const int64_t basePos = mappedPos.count(lineNo) ? mappedPos[lineNo] : 0;
                    pnextCache[lineNo] = delta + basePos;
                }
            } else if (excEnc == "varint") {
                /* Absolute contentIdx then a zigzag-LEB128 delta per pair (older layout). */
                excDec->decode_line(raw.data(), (uint32_t)raw.size(), UINT8_MAX, false);
                for (int32_t i = 0; i < pairCount; ++i) {
                    uint64_t ci = 0, z = 0;
                    if (!nextVarint(ci) || !nextVarint(z)) {
                        LOG_ERROR("PNEXT exception varint stream truncated");
                        return -1;
                    }
                    const int64_t delta = (int64_t)(z >> 1) ^ -(int64_t)(z & 1); /* de-zigzag */
                    const uint32_t lineNo = (uint32_t)ci;
                    const int64_t basePos = mappedPos.count(lineNo) ? mappedPos[lineNo] : 0;
                    pnextCache[lineNo] = delta + basePos;
                }
            } else {
                /* Legacy layout: fixed int32 x 2 per pair. */
                std::vector<int32_t> excBuf((size_t)pairCount * 2, 0);
                excDec->decode_line((uint8_t*)excBuf.data(),
                    (uint32_t)((size_t)pairCount * 2 * sizeof(int32_t)), UINT8_MAX, false);
                for (int32_t i = 0; i < pairCount; ++i) {
                    uint32_t contentIdx = (uint32_t)excBuf[2 * i];
                    int32_t delta = excBuf[2 * i + 1];
                    /* The exception stores the 0-based data-line index (lineIdx - headEndLine), which is also lineNo here. */
                    uint32_t lineNo = contentIdx;
                    int64_t pos = mappedPos.count(lineNo) ? mappedPos[lineNo] : 0;
                    pnextCache[lineNo] = delta + pos;
                }
            }
        }
    }

    // Build qname -> lines mapping from decodedQnames (only mapped, paired lines).
    std::unordered_map<std::string, std::vector<uint32_t>> qnameToLines;
    qnameToLines.reserve(samLine);
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        auto flagIt = mappedFlag.find(lineNo);
        bool paired = flagIt != mappedFlag.end() && (flagIt->second & 0x1) && !(flagIt->second & 0x8);
        if (paired && !decodedQnames[lineNo].empty()) {
            qnameToLines[decodedQnames[lineNo]].push_back(lineNo);
        }
    }

    // Rebuild PNEXT for every line.
    std::vector<std::string>& pnextCacheOut = tlenPreDecodedFields[7];
    pnextCacheOut.clear();
    pnextCacheOut.reserve(samLine);
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        int64_t pnext = 0;
        auto cacheIt = pnextCache.find(lineNo);
        if (cacheIt != pnextCache.end()) {
            pnext = cacheIt->second; // exception: stored value
        } else {
            // Rebuild from mate: the encoder guarantees a non-exception line
            // belongs to a QNAME group of exactly two mutually-mapped
            // records, so the other record in the group is the mate.
            auto flagIt = mappedFlag.find(lineNo);
            bool paired = flagIt != mappedFlag.end() && (flagIt->second & 0x1) && !(flagIt->second & 0x8);
            if (paired && !decodedQnames[lineNo].empty()) {
                const auto& mates = qnameToLines[decodedQnames[lineNo]];
                if (mates.size() == 2) {
                    for (uint32_t ml : mates) {
                        if (ml == lineNo) continue;
                        pnext = mappedPos.count(ml) ? mappedPos[ml] : 0;
                        break;
                    }
                }
            }
        }
        nextMappedPos[lineNo] = pnext;
        pnextCacheOut.emplace_back(std::to_string(pnext) + '\t');
    }
    return 0;
}

/* Full mate index: built at once, once all (pos, pnext) pairs are known. */
void SamCodecActuator::buildTlenMateIndex()
{
    tlenMateIndex.clear();
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        auto posIt = mappedPos.find(lineNo);
        auto pnextIt = nextMappedPos.find(lineNo);
        if (posIt != mappedPos.end() && pnextIt != nextMappedPos.end()) {
            tlenMateIndex[std::make_pair(posIt->second, pnextIt->second)] = lineNo;
        }
    }
}

int32_t SamCodecActuator::preDecodeForTLEN() {
    if (samLine == 0 || !meta.isMember("sam")) {
        return 0;
    }

    Json::Value& streams = meta["sam"]["streams"];

    /* QNAME (field 0) first: the PNEXT mode below pairs records by it. See rebuildQnameColumn. */
    if (rebuildQnameColumn() != 0) {
        return -1;
    }

    /* FLAG (field 1) first: the flag bits decide whether PNEXT holds a delta (see predecodeFlagColumn). */
    if (predecodeFlagColumn() != 0) {
        return -1;
    }

    /* Pre-decode RNAME (field 2) so the POS delta chain can detect chromosome switches; a separate decoder. */
    if (predecodeRnameColumn() != 0) {
        return -1;
    }

    /* Fields needed to reconstruct TLEN: POS(3), CIGAR(5), PNEXT(7). */
    static const uint32_t tlenFields[] = {3, 5, 7};

    for (uint32_t f : tlenFields) {
        TlenColumnPredecoder pre;
        if (!buildTlenColumnPredecoder(f, pre)) {
            continue;   /* nothing to pre-decode here (no stream, or an unreadable coder) */
        }
        if (rebuildTlenColumn(f, pre) != 0) {
            return -1;
        }
    }

    /* PNEXT (7) in pnext_qname_rebuild mode: only exceptions in the stream; the rest from QNAME groups. */
    std::string pnextMode = streams.isValidIndex(7) && streams[7].isMember("mode")
        ? streams[7]["mode"].asString() : "";
    if (pnextMode == "pnext_qname_rebuild" && rebuildPnextByQname() != 0) {
        return -1;
    }

    buildTlenMateIndex();
    return 0;
}


/* Reconstruct one split segment of this line's ID by whatever layout wrote it, and append it. Returns 0, or -1 on a corrupt stream. */
int32_t SamCodecActuator::reconstructIdSegment(uint32_t splitIdx, Json::Value& splitMeta,
                                                RoughIOBlock* outputBlock, uint32_t& idLength)
{
    const uint32_t splitDstLen = splitMeta["dstlen"].asUInt();
    std::string coderName = splitMeta["coder"]["magic"].asString();

    /* "cnt" is the numeric layout from here on: initDecoder already expanded its deltas. */
    const std::string idMode = splitMeta.isMember("mode") ? splitMeta["mode"].asString()
                                                         : std::string();
    const bool numericMode = (idMode == "numeric" || idMode == "cnt");
    if (numericMode) {
        /* Numeric layout: serve one value per line, render to decimal text, re-append the split symbol. */
        readOffset += splitDstLen;
        idLength += splitDstLen;

        if (splitIdx >= idNumericBufs.size() || idNumericBufs[splitIdx] == nullptr ||
            idNumericPos[splitIdx] >= idNumericLens[splitIdx]) {
            LOG_ERROR("id numeric segment(%u) exhausted", splitIdx);
            return -1;
        }
        uint8_t* numBuf = idNumericBufs[splitIdx];
        uint32_t& npos = idNumericPos[splitIdx];
        uint64_t& nacc = idNumericAcc[splitIdx];
        const uint32_t nlen = idNumericLens[splitIdx];

        uint32_t zz = 0;
        if (!readVarint(numBuf, nlen, npos, zz)) {
            LOG_ERROR("Corrupt id numeric segment(%u): truncated varint", splitIdx);
            return -1;
        }
        const uint32_t delta = zz >> 1;
        nacc = ((zz & 1u) == 0) ? (nacc + delta) : (nacc - delta);

        const uint8_t sep = (splitIdx < idAnalysis.symbols.size())
                            ? (uint8_t)idAnalysis.symbols[splitIdx] : (uint8_t)'\t';
        char tmp[24];
        int n = snprintf(tmp, sizeof(tmp), "%llu", (unsigned long long)nacc);
        if (n <= 0 || outputBlock->getRemain() < (uint32_t)n + 1) {
            LOG_ERROR("Reconstruct id numeric segment(%u) overflow", splitIdx);
            return -1;
        }
        memcpy(outputBlock->getCurrent(), tmp, (uint32_t)n);
        outputBlock->setDataLen(outputBlock->getDataLen() + (uint32_t)n);
        outputBlock->getCurrent()[0] = sep;
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return 0;
    }
    if (!idDictEntries.empty() && !idDictEntries[splitIdx].empty()) {
        /* Dictionary segment: one index per line; initDecoder turned it into one absolute index per line. */
        if (splitIdx >= idNumericBufs.size() || idNumericBufs[splitIdx] == nullptr ||
            idNumericPos[splitIdx] >= idNumericLens[splitIdx]) {
            LOG_ERROR("id dictionary segment(%u) exhausted", splitIdx);
            return -1;
        }
        uint32_t idx = 0;
        if (!readVarint(idNumericBufs[splitIdx], idNumericLens[splitIdx],
                        idNumericPos[splitIdx], idx) ||
            idx >= idDictEntries[splitIdx].size()) {
            LOG_ERROR("Reconstruct id dictionary segment(%u): bad index", splitIdx);
            return -1;
        }
        const std::string& text = idDictEntries[splitIdx][(size_t)idx];
        if (outputBlock->getRemain() < text.size()) {
            LOG_ERROR("Reconstruct id dictionary segment(%u) overflow", splitIdx);
            return -1;
        }
        memcpy(outputBlock->getCurrent(), text.data(), text.size());
        outputBlock->setDataLen(outputBlock->getDataLen() + (uint32_t)text.size());
        return 0;
    }
    if (!idHexAlphabets.empty() && !idHexAlphabets[splitIdx].chars.empty()) {
        /* Fixed-alphabet segment: the length when they vary, one uniform char code each, then the split symbol. */
        const std::string& alpha = idHexAlphabets[splitIdx].chars;
        RangeCoder& hrc = idHexCoders[splitIdx];
        uint32_t len = idHexAlphabets[splitIdx].minLen;
        if (idHexAlphabets[splitIdx].maxLen > len) {
            const uint32_t range = idHexAlphabets[splitIdx].maxLen - len + 1;
            const uint32_t v = hrc.GetFreq(range);
            hrc.Decode(v, 1);
            len += v;
        }
        const uint8_t sep = (splitIdx < idAnalysis.symbols.size())
                            ? (uint8_t)idAnalysis.symbols[splitIdx] : (uint8_t)'\t';
        if (hrc.err != 0 || len == 0 || len > 256) {
            LOG_ERROR("Reconstruct id fixed-alphabet segment(%u) failed", splitIdx);
            return -1;
        }
        if (outputBlock->getRemain() < len + 1) {
            LOG_ERROR("Reconstruct id fixed-alphabet segment(%u) overflow", splitIdx);
            return -1;
        }
        char tmp[257];
        for (uint32_t c = 0; c < len; ++c) {
            const uint32_t v = hrc.GetFreq((uint32_t)alpha.size());
            hrc.Decode(v, 1);
            if (hrc.err != 0 || v >= alpha.size()) {
                LOG_ERROR("Reconstruct id fixed-alphabet segment(%u) failed", splitIdx);
                return -1;
            }
            tmp[c] = alpha[v];
        }
        memcpy(outputBlock->getCurrent(), tmp, len);
        outputBlock->setDataLen(outputBlock->getDataLen() + len);
        outputBlock->getCurrent()[0] = sep;
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return 0;
    }
    if (!idConstTexts.empty() && !idConstTexts[splitIdx].empty()) {
        /* Constant segment: stored once; its text carries the trailing split symbol. */
        const std::string& text = idConstTexts[splitIdx];
        if (outputBlock->getRemain() < text.size()) {
            LOG_ERROR("Reconstruct id constant segment(%u) overflow", splitIdx);
            return -1;
        }
        memcpy(outputBlock->getCurrent(), text.data(), text.size());
        outputBlock->setDataLen(outputBlock->getDataLen() + (uint32_t)text.size());
        return 0;
    }
    if (!idIntActive.empty() && idIntActive[splitIdx]) {
        /* Value-domain segment: one value per line, then the split symbol. */
        const uint8_t sep = (splitIdx < idAnalysis.symbols.size())
                            ? (uint8_t)idAnalysis.symbols[splitIdx] : (uint8_t)'\t';
        const uint64_t v = idIntModes[splitIdx].decode(idIntCoders[splitIdx]);
        char tmp[24];
        int n = snprintf(tmp, sizeof(tmp), "%llu", (unsigned long long)v);
        if (n <= 0 || outputBlock->getRemain() < (uint32_t)n + 1) {
            LOG_ERROR("Reconstruct id value-domain segment(%u) overflow", splitIdx);
            return -1;
        }
        memcpy(outputBlock->getCurrent(), tmp, (uint32_t)n);
        outputBlock->setDataLen(outputBlock->getDataLen() + (uint32_t)n);
        outputBlock->getCurrent()[0] = sep;
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return 0;
    }

    // Decode segment (legacy textual layout)
    int32_t segmentLen = idDecoders[splitIdx]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(),
        (splitIdx < idAnalysis.symbols.size()) ? idAnalysis.symbols[splitIdx] : UINT8_MAX, false);
    if (segmentLen < 0) {
        LOG_ERROR("Decode id field segment(%u) failed: %d", splitIdx, segmentLen);
        return -1;
    }
    readOffset += splitDstLen;
    idLength += splitDstLen;
    outputBlock->setDataLen(outputBlock->getDataLen() + segmentLen);
    return 0;
    }

int32_t SamCodecActuator::decompressIdField(uint32_t fieldIdx, Json::Value& fieldMeta, RoughIOBlock* outputBlock) {
    if (fieldIdx != 0) {
        return -1;
    }
    // Handle ID field with split compression
    Json::Value& idStreams = fieldMeta["streams"];
    uint32_t idLength = 0;

    if (idUsesQnameCoder && !idDecoders.empty()) {
        /* Single-stream coder_qname: decode one whole QNAME per line (including the trailing '\t'). */
        int32_t segLen = idDecoders[0]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(),
            UINT8_MAX, false);
        if (segLen < 0) {
            LOG_ERROR("Decode QNAME (coder_qname) failed: %d", segLen);
            return -1;
        }
        readOffset += idStreams[0]["dstlen"].asUInt();
        outputBlock->setDataLen(outputBlock->getDataLen() + segLen);
        return (int32_t)idStreams[0]["dstlen"].asUInt();
    }

    // Reconstruct ID from split segments
    for (uint32_t splitIdx = 0; splitIdx < idStreams.size(); ++splitIdx) {
        if (reconstructIdSegment(splitIdx, idStreams[splitIdx], outputBlock, idLength) != 0) {
            return -1;
        }
    }
    return idLength;
}

int32_t SamCodecActuator::decompressChrName(uint32_t fieldIdx, uint32_t lineNo, RoughIOBlock* outputBlock) {
    uint16_t chrIndex = 0;
    bool byteCode = false;
    if (fieldIdx == 6) {
        /* RNEXT compact one-byte code stream (see compressChrName). */
        const Json::Value& st = meta["sam"]["streams"];
        if (st.isValidIndex(6) && st[6].isMember("rn_enc") &&
            st[6]["rn_enc"].asString() == "byte") {
            byteCode = true;
        }
    }
    if (byteCode) {
        uint8_t code = 0;
        if (fieldDecoders[fieldIdx]->decode_line(&code, 1, UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode RNEXT code failed, lineNo = %u", lineNo);
            return -1;
        }
        if (code == 0x00) {
            chrIndex = 0xFFFE; /* "=" */
        } else if (code == 0x01) {
            chrIndex = 0xFFFF; /* "*" */
        } else if (code == 0xFF) {
            uint8_t le[2] = {0, 0};
            if (fieldDecoders[fieldIdx]->decode_line(le, 2, UINT8_MAX, false) < 0) {
                LOG_ERROR("Decode RNEXT escape index failed, lineNo = %u", lineNo);
                return -1;
            }
            chrIndex = (uint16_t)(le[0] | (le[1] << 8));
        } else {
            chrIndex = (uint16_t)(code - 2);
        }
    } else {
        if (fieldDecoders[fieldIdx]->decode_line((uint8_t*)&chrIndex, sizeof(chrIndex), UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode chr name field(%u) failed, lineNo = %u", fieldIdx, lineNo);
            return -1;
        }
    }
    /* Record the chromosome index per line ("*"=0xFFFF, "="=0xFFFE), mirroring compressChrName. */
    if (fieldIdx == 2) {
        mappedChr[lineNo] = chrIndex;
    } else if (fieldIdx == 6) {
        nextMappedChr[lineNo] = chrIndex;
    }
    if (chrIndex == 0xFFFF) {
        *outputBlock->getCurrent() = '*';
        *(outputBlock->getCurrent() + 1) = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 2);
        return 2;
    } else if (chrIndex == 0xFFFE) {
        *outputBlock->getCurrent() = '=';
        *(outputBlock->getCurrent() + 1) = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 2);
        return 2;
    } else {
        std::string chrName = SamInfo::getInstance().getChromosomeInfo(chrIndex).name;
        memcpy(outputBlock->getCurrent(), chrName.c_str(), chrName.length());
        outputBlock->setDataLen(outputBlock->getDataLen() + chrName.length());
        *outputBlock->getCurrent() = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return  chrName.length() + 1;
    }
}

int32_t SamCodecActuator::readMatchLine(uint32_t fieldIdx, uint8_t* dst, uint32_t len, uint32_t lineNo)
{
    if (matchBlockDecode) {
        /* coder_fc: the whole match stream was decoded up front; slice this record out. */
        if (matchBlockOffset + len > matchBlockLength) {
            LOG_ERROR("SEQ match block exhausted in block %llu, line %u (need %u, left %u)",
                      inBlockPtr->getBlockId(), lineNo, len, matchBlockLength - matchBlockOffset);
            return -1;
        }
        std::memcpy(dst, matchBlockBuffer + matchBlockOffset, len);
        matchBlockOffset += len;
        return (int32_t)len;
    }
    return fieldDecoders[fieldIdx]->decode_line(dst, len, UINT8_MAX, false);
}

int32_t SamCodecActuator::decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& totalBaseLen, RoughIOBlock* outputBlock) {
    /* ensureCapacity must not be called here: the caller captured basePtr for QUAL, and a realloc would dangle it. */
    /* Three layouts: the whole-column stream ("coder"), the multi-stream reference/consensus layout
       ("streams"), and the same multi-stream layout with no reference at all ("streams" + "noref"). */
    const bool isSelfReference = baseSelfRefActive && fieldMeta.isMember("selfref");
    const bool isStreamsLayout = fieldMeta.isMember("streams");
    const bool isNoRefBlock = isStreamsLayout && !isSelfReference &&
                              fieldMeta.isMember("noref") && fieldMeta["noref"].asUInt() != 0;
    const bool isUserReference = isSelfReference || isStreamsLayout;
    uint32_t actualBaseLen = 0;
    if (!isUserReference) {
        /* Each row length is known, so records are taken by exact length (split_ch = UINT8_MAX); a whole-block coder is copied from staging. */
        const std::string seqCoder = fieldMeta["coder"]["magic"].asString();
        if (minBaseLength == maxBaseLength) {
            actualBaseLen = baseLengthBuffer[lineNo] == 0 ? maxBaseLength : baseLengthBuffer[lineNo];
            if (CoderFactory::decoderIsWholeBlock(seqCoder)) {
                memcpy(outputBlock->getCurrent(), pBaseOut, actualBaseLen);
                pBaseOut += actualBaseLen;
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
            } else {
                int32_t decLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), actualBaseLen, UINT8_MAX, false);
                if (decLen < 0 || (uint32_t)decLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %lld, line %d, expect len %d, actual len %d", (long long)inBlockPtr->getBlockId(), lineNo, actualBaseLen, decLen);
                    return -1;
                }
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
            }

            *(outputBlock->getCurrent()) = '\t';
            outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        } else {
            if (CoderFactory::decoderIsWholeBlock(seqCoder)) {
                /* A whole-block coder's stream is already decoded into staging; copy this record up to its tab. */
                uint8_t* pBaseTmp = outputBlock->getCurrent();
                uint8_t* ptr = pBaseOut;
                uint8_t* pBaseEnd = outputBlock->getBuffer() + outputBlock->getBufferSize();
                for (; ptr < pBaseEnd; ++ptr) {
                    *pBaseTmp++ = *ptr;
                    if (*ptr == '\t') {
                        break;
                    }
                }
                actualBaseLen = ptr - pBaseOut + 1;
                pBaseOut += actualBaseLen;
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
                actualBaseLen -= 1; // Remove \t length
            } else {
                /* The per-row limit must be maxBaseLength + 1: the row stores its terminating tab and split mode reports BUF_SMALL at out_len chars. */
                int32_t decLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), maxBaseLength + 1, '\t', false);
                if (decLen <= 0) {
                    LOG_ERROR("base decode failed in block %lld, line %d: %d", (long long)inBlockPtr->getBlockId(), lineNo, decLen);
                    return -1;
                }
                actualBaseLen = decLen;
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
                actualBaseLen -= 1; // Remove \t length
            }
        }
    } else {
        actualBaseLen =  baseLengthBuffer[lineNo];
        if (actualBaseLen == 0) {
            for (auto pairIter = unmapedReadLength.begin(); pairIter < unmapedReadLength.end(); ++pairIter) {
                if (pairIter->first == lineNo) {
                    actualBaseLen = pairIter->second;
                    break;
                }
            }
            if (actualBaseLen == 0) {
                actualBaseLen = maxBaseLength;
            }
        }
        Json::Value& baseStream = fieldMeta["streams"];
        /* A payload byte above 2-bit range is one of the record's own characters (record does not use the reference). */
        const bool litBases = fieldMeta.isMember("litbases");
        /* `src` is the buffer this record's payload was read into. */
        const auto writeDirectBases = [&](const uint8_t* src, int32_t len) {
            if (litBases) {
                memcpy(outputBlock->getCurrent(), src, (size_t)len);
            } else {
                for (int32_t o = 0; o < len; ++o) {
                    outputBlock->getCurrent()[o] = atcg4[src[o]];
                }
            }
        };
        /* The layout is identified by its sub-stream names; initDecoder has already expanded both RLE halves. */
        if (baseStream[0]["sname"].asString() == "m") {
            int32_t decoderLen = 0;
            uint16_t mapFlag = mappedFlag.find(lineNo) == mappedFlag.end() ? 4 : mappedFlag[lineNo];
            // Not matched, or the block was coded with no reference at all: the payload is the bases themselves
            if (isNoRefBlock || (mapFlag & 0x04)) {
                decoderLen = readMatchLine(fieldIdx, baseSquashBuffer, actualBaseLen, lineNo);
                if ((uint32_t)decoderLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                    return -1;
                }
                writeDirectBases(baseSquashBuffer, decoderLen);
            } else {
                /* CIGAR-segment reference rebuild, mirroring compression: M/=/X from the reference, I/S direct 2-bit, D/N advance the reference. */
                bool findMappedPos = false;
                int64_t refeMappedPos = 0;
                if (isSelfReference) {
                    /* The consensus covers only the block's window, so the position is rebased on it (see seqRecordUsesSquashWindow). */
                    const uint16_t chrIdx = mappedChr.find(lineNo) == mappedChr.end() ? 0xFFFF : mappedChr[lineNo];
                    const auto crlIt = cigarReadLen.find(lineNo);
                    const uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
                    const auto posIt = mappedPos.find(lineNo);
                    const uint64_t pos = (posIt != mappedPos.end()) ? (uint64_t)posIt->second : 0;
                    const uint16_t mapFlag = mappedFlag.find(lineNo) == mappedFlag.end() ? 4 : mappedFlag[lineNo];
                    findMappedPos = seqRecordUsesSquashWindow(chrIdx, mapFlag, pos, refConsumed,
                                                              baseSelfRefWindowStart,
                                                              baseSelfRefWindowBases, refeMappedPos);
                } else {
                do {
                    uint16_t chrIdx = mappedChr.find(lineNo) == mappedChr.end() ? 0xFFFF : mappedChr[lineNo];
                    if (chrIdx == 0xFFFF || chrIdx == 0xFFFE) {
                        break;
                    }
                    int64_t refeChrPos = SamInfo::getInstance().getPositionByIndex(chrIdx);
                    if (refeChrPos == -1) {
                        break;
                    }
                    refeMappedPos = refeChrPos + mappedPos[lineNo] - 1;
                    auto crlIt = cigarReadLen.find(lineNo);
                    uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
                    uint64_t needSquash = (uint64_t)(refConsumed >> 2) + !!(refConsumed & 0x3) + 1;
                    if (refeMappedPos < 0 || ((uint64_t)(refeMappedPos >> 2)) + needSquash > (uint64_t)pRefeGene->getSquashLength()) {
                        break;
                    }
                    findMappedPos = true;
                } while(0);
                }

                if (!findMappedPos) {
                    decoderLen = readMatchLine(fieldIdx, baseSquashBuffer, actualBaseLen, lineNo);
                    if ((uint32_t)decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    writeDirectBases(baseSquashBuffer, decoderLen);
                } else {
                    decoderLen = readMatchLine(fieldIdx, baseDiffSquashBuffer, actualBaseLen, lineNo);
                    if ((uint32_t)decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d,expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    uint8_t* out = outputBlock->getCurrent();
                    /* No op list (or empty): the record took the encoder's fallback and its bases are not reference-coded. */
                    const bool hasOps = (lineNo < cigarOpList.size()) && !cigarOpList[lineNo].empty();
                    if (!hasOps) {
                        writeDirectBases(baseDiffSquashBuffer, (int32_t)actualBaseLen);
                    } else {
                        /* The same shared CIGAR walk the writer runs (see seqWalkCigar). */
                        const std::vector<CigarOp>& ops = cigarOpList[lineNo];
                        const SeqWalkEnd end = seqWalkCigar(ops, actualBaseLen, refeMappedPos,
                            [&](uint32_t readPos, int64_t pos, uint32_t len) {
                                if (isSelfReference) {
                                    /* The consensus is read the same way as the file's reference. */
                                    seqSquashStretch(baseSelfRefSquash.data(), refeStrecchBuffer, len,
                                                     (uint64_t)pos);
                                    for (uint32_t i = 0; i < len; ++i) {
                                        out[readPos + i] =
                                            atcg4[(refeStrecchBuffer[i] ^ baseDiffSquashBuffer[readPos + i]) & 0x3];
                                    }
                                } else {
                                    pRefeGene->getStretch2Bits1Char(refeStrecchBuffer, len, pos);
                                    for (uint32_t i = 0; i < len; ++i) {
                                        baseSquashBuffer[i] = refeStrecchBuffer[i] ^ baseDiffSquashBuffer[readPos + i];
                                    }
                                    pRefeGene->getActgFrom2Bits(baseSquashBuffer, len, out + readPos);
                                }
                                return SeqWalkAction::Continue;
                            },
                            [&](uint32_t readPos, uint32_t len) {
                                for (uint32_t i = 0; i < len; ++i) {
                                    out[readPos + i] = atcg4[baseDiffSquashBuffer[readPos + i] & 0x3];
                                }
                                return SeqWalkAction::Continue;
                            },
                            [](uint32_t) { return SeqWalkAction::Continue; });
                        /* The ops did not describe this read: the writer fell back to actgEncode()
                           for the whole payload, so the bytes are 2-bit base codes, not ASCII. */
                        if (!end.complete) {
                            for (int32_t i = 0; i < (int32_t)actualBaseLen; ++i) {
                                out[i] = atcg4[baseDiffSquashBuffer[i] & 0x3];
                            }
                        }
                    }
                }
            }

            /* Write the exception characters back, one cursor per stream (see SeqExceptionClass). */
            for (uint32_t n = 0; n < actualBaseLen; ++n) {
                const uint32_t pos = totalBaseLen + n;
                for (size_t s = 0; s < seqExc.size(); ++s) {
                    SeqExcStream& stream = seqExc[s];
                    if (stream.off < stream.pos.size() && stream.pos[stream.off] == pos) {
                        *(outputBlock->getCurrent() + n) = (char)stream.byte;
                        stream.off++;
                        break;
                    }
                }
            }
            totalBaseLen += actualBaseLen;
            outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);

            *(outputBlock->getCurrent()) = '\t';
            outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        } else {
            LOG_ERROR("Not supported coder name:%s", fieldMeta["coder"]["magic"].asString().c_str());
            return -1;
        }
    }
    return actualBaseLen;
}

int32_t SamCodecActuator::decompressQuality(uint8_t* basePtr, uint32_t actualBaseLen, RoughIOBlock* outputBlock) {
    uint8_t* dst = outputBlock->getCurrent();
    /* One call, whichever coder wrote the column: the record goes to dst, with the record's decoded SEQ as context. */
    if (qualDecoder->decode_record(dst, actualBaseLen, basePtr, actualBaseLen) < 0) {
        LOG_ERROR("Decode quality failed, len = %u", actualBaseLen);
        return -1;
    }

    /* Missing quality (single '*') was expanded to seqLen '*' on compression; fold it back here. */
    bool missing = true;
    for (uint32_t i = 0; i < actualBaseLen; ++i) {
        if (dst[i] != '*') {
            missing = false;
            break;
        }
    }
    if (missing) {
        dst[0] = '*';
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
    } else {
        outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
    }
    *(outputBlock->getCurrent()) = '\t';
    outputBlock->setDataLen(outputBlock->getDataLen() + 1);
    return actualBaseLen;
}

int32_t SamCodecActuator::decompressCigar(uint32_t fieldIdx, uint8_t splitFlag, uint32_t lineIdx, RoughIOBlock* outputBlock) {
    int32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(), splitFlag, false);
    if (fieldLen < 0) {
        LOG_ERROR("Decode cigar field(%u) failed: %d", fieldIdx, fieldLen);
        return -1;
    }
    /* cigarOpList is indexed by the 0-based data line, always non-negative here. */
    const size_t contentIdx = static_cast<size_t>(lineIdx - headEndLine);
    if (fieldLen > 1) {
        uint32_t seqLength = cigarSeqConsumed(outputBlock->getCurrent(), fieldLen);
        baseLengthBuffer[lineIdx] = seqLength;
        cigarReadLen[lineIdx] = cigarRefConsumed(outputBlock->getCurrent(), fieldLen);
        if (cigarOpList.size() <= contentIdx) cigarOpList.resize(contentIdx + 1);
        parseCigarOps(outputBlock->getCurrent(), (uint32_t)fieldLen, cigarOpList[contentIdx]);
    } else {
        baseLengthBuffer[lineIdx] = 0;
        cigarReadLen[lineIdx] = 0;
        if (cigarOpList.size() <= contentIdx) cigarOpList.resize(contentIdx + 1);
        cigarOpList[contentIdx].clear();
    }

    outputBlock->setDataLen(outputBlock->getDataLen() + fieldLen);
    return fieldLen;
}

int32_t SamCodecActuator::buildSamIndex() {
    if (!pbgzEngine->getParameter().isMakeIndex) {
        return 0;
    }

    uint32_t totalLineNum = inBlockPtr->getNpos().size();
    struct SortKey {
        SortKey(uint16_t c = 0, int64_t p = 0) : chrIndex(c), mapPos(p) {}
        uint16_t chrIndex;
        int64_t mapPos;
        bool operator<(const SortKey& other) const {
            if (chrIndex != other.chrIndex) {
                return chrIndex < other.chrIndex;
            }
            return mapPos < other.mapPos;
        }
    };

    SortKey lastKey = {0, 0};
    std::map<uint16_t, std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>> chrBlockStats;
    const uint32_t MAX_SPLITS_PER_CHR = 1000;

    for (uint32_t lineNo = headEndLine; lineNo < totalLineNum; ++lineNo) {
        if (mappedFlag.find(lineNo) == mappedFlag.end()) {
            continue;
        }
        if (mappedChr.find(lineNo) == mappedChr.end()) {
            continue;
        }
        if (mappedPos.find(lineNo) == mappedPos.end()) {
            continue;
        }

        uint16_t flag = mappedFlag[lineNo];
        uint16_t chrIndex = mappedChr[lineNo];
        int64_t mapPos = mappedPos[lineNo];

        if ((flag & 0x04) == 0 && chrIndex != 0xFFFF && chrIndex != 0xFFFE) {
            SortKey currentKey = {chrIndex, mapPos};

            if (currentKey < lastKey) {
                LOG_ERROR("SAM block %d is not sorted by chrIndex and mapPos: line %d (chr=%u,pos=%ld) < line %d-1 (chr=%u,pos=%ld)",
                    inBlockPtr->getBlockId(), lineNo, chrIndex, mapPos, lineNo,
                    lastKey.chrIndex, lastKey.mapPos);
                static bool isPrint = true;
                static std::mutex mutex;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    if (isPrint) {
                        isPrint = false;
                        fprintf(stderr, "The SAM file is unsorted, which may cause the index file not to be created.\n");
                        pbgzEngine->parameter.isMakeIndex = false;
                    }
                }
                return -1;
            }

            lastKey = currentKey;

            auto& items = chrBlockStats[chrIndex];
            if (items.empty()) {
                items.emplace_back(mapPos, mapPos, 1);
            } else {
                uint32_t& lastPos = std::get<1>(items.back());
                uint32_t& count = std::get<2>(items.back());

                if (count >= MAX_SPLITS_PER_CHR && mapPos != lastPos) {
                    items.emplace_back(mapPos, mapPos, 1);
                } else {
                    lastPos = mapPos;
                    count++;
                }
            }
        }
    }

    for (const auto& chrPair : chrBlockStats) {
        uint16_t chrIndex = chrPair.first;
        const auto& items = chrPair.second;

        int64_t refPos = SamInfo::getInstance().getPositionByIndex(chrIndex);
        if (refPos == -1) {
            continue;
        }

        for (const auto& item : items) {
            uint32_t firstPos = std::get<0>(item);
            uint32_t count = std::get<2>(item);
            SamIndex::getInstance().addSamIndex(chrIndex, firstPos, count, inBlockPtr->getBlockId());
        }
    }

    return 0;
}