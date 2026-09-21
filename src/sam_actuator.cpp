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
#include "coder/coder_io.h"
#include "coder/coder_fc.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_arith.h"
#include "coder/coder_qual.h"
#include "utils/md5_util.h"
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
    optionRecLines.clear();
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

/*
 * Where the QNAME column's lines are, as the layout encoders want them (see sam_qname_column.h):
 * the block being compressed, the tab positions the parse recorded, and how many leading lines
 * are header.
 */
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
    /* 0 = the whole block: the line bound the module's encoders take is the trial's, and this is
       not the trial (see the declaration). */
    return encodeQnameSplit(qnameColumnInput(), idAnalysis, outBlockPtr, &ioErrSink, 0, fieldMeta,
                            fieldSrcLen);
}

int32_t SamCodecActuator::compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    return encodeQnameByQname(qnameColumnInput(), outBlockPtr, &ioErrSink, 0, fieldMeta,
                              fieldSrcLen);
}


/*
 * Analyse one SAM header line (it starts with '@'): capture the chromosome table from @SQ,
 * check the reference path it names against the one in use, and read the bwa command line out
 * of @PG for the same check.
 *
 * Answers 0 when the line has been dealt with, which is what the loop's `continue` used to
 * mean; an unknown header type, an unparsable @SQ line or a chromosome missing from the index
 * is an error.
 */
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
                // Parse filename from refFilePath, compatible with local and network paths
                std::string refFileName;
                size_t lastSlash = refFilePath.find_last_of("/\\");
                if (lastSlash != std::string::npos) {
                    refFileName = refFilePath.substr(lastSlash + 1);
                } else {
                    refFileName = refFilePath;
                }

                LOG_INFO("Reference fasta name: %s.", refFilePath.c_str());
                if (pRefeGene != nullptr) {
                    std::string inputFastq = PathUtil::getFileName(pRefeGene->getFastaFileName());
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
                // Split command by spaces and take the third as reference gene name
                std::stringstream ss(command);
                std::string item;
                std::vector<std::string> tokens;
                while (std::getline(ss, item, ' ')) {
                    if (!item.empty()) {
                        tokens.push_back(item);
                    }
                }
                if (tokens.size() >= 3) {
                    std::string refGeneName;
                    if (tokens.size() >= 2 && tokens[0] == "bwa" && tokens[1] == "mem") {
                        int nonOptionCount = 0;
                        for (size_t i = 2; i < tokens.size(); i++) {
                            if (tokens[i].substr(0, 1) != "-") {
                                nonOptionCount++;
                                if (nonOptionCount == 1) {
                                    refGeneName = tokens[i];
                                    break;
                                }
                            }
                        }
                    }

                    if (!refGeneName.empty()) {
                        LOG_INFO("Reference gene name extracted from @PG CL: %s", refGeneName.c_str());
                        if (pRefeGene != nullptr) {
                            std::string inputFastq = PathUtil::getFileName(pRefeGene->getFastaFileName());
                            if (PathUtil::isGzFile(pRefeGene->getFastaFileName())) {
                                inputFastq = PathUtil::getFileNameFromGz(pRefeGene->getFastaFileName());
                            }
                            if (refGeneName != inputFastq) {
                                fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fasta is %s \n", refGeneName.c_str(), inputFastq.c_str());
                                isCheckPG = true;
                            }
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

/*
 * Analyse one record line (it does not start with '@'): hand the ID column to the ID analysis,
 * check RNAME against the chromosome index, and measure the CIGAR-derived length against the
 * SEQ column, the SEQ length against QUAL, the base count and the quality histogram. The
 * field offset table it builds is kept in contentPos.
 *
 * `qualityFrequnce` is the caller's 256-entry histogram, indexed by the quality byte.
 */
int32_t SamCodecActuator::scanDataLine(const std::string& line, uint32_t idx,
                                      std::pair<uint8_t, uint32_t>* qualityFrequnce)
{
    std::vector<int64_t> linePos;
    uint32_t baseLen = 0;
    bool lineCigarMatchFlag = false;
    uint32_t baseFieldLen = 0;
    for (uint32_t i = 0; i < line.length(); ++i) {
        if (line.at(i) == '\t' || line.at(i) == '\n') {
            // First tab before is ID column, need to split and analyze ID column
            if (linePos.empty()) {
                if (contentPos.empty()) {
                    // Pass ID column content (from line start to first tab position)
                    preAnalysisIdFirstLine((uint8_t*)line.data(), i + 1);
                } else {
                    // Pass ID column content (from line start to first tab position)
                    preAnalysisIdLine((uint8_t*)line.data(), i + 1);
                }
            } else if (linePos.size() == 2) {  // RNAME is the 3th field
                uint32_t chrLen = i - linePos.at(1) - 1;
                std::string chrName = line.substr(linePos.at(1) + 1, chrLen);
                if (chrName != "*" && chrName != "=") {
                    if (SamInfo::getInstance().getChrNameIndex(chrName) == 65535) {
                        if (headEndLine >= 0) {
                            SamInfo::getInstance().clearChromosomeInfo();
                        }
                        return -1;
                    }
                }
            } else if (linePos.size() == 5) {  // CIGAR is the 6th field
                uint32_t cigarLen = i - linePos.at(4) - 1;
                if (cigarLen > 1) {
                    lineCigarMatchFlag = true;
                    uint8_t* cigarBegin  = (uint8_t*)line.c_str() + linePos.at(4) + 1;
                    baseFieldLen = parseCigar(cigarBegin, cigarLen);
                } else {
                    lineCigarMatchFlag = false;
                }
            } else if (linePos.size() == 9) {   // Base is the 10th field
                baseLen = i - linePos.at(8) - 1;
                if (baseLen > maxBaseLength) {
                    maxBaseLength = baseLen;
                }
                if (!lineCigarMatchFlag) {
                    if (baseLen < minBaseLength) {
                        minBaseLength =  baseLen;
                    }
                } else {
                    if (baseFieldLen != baseLen) {
                        if (headEndLine >= 0) {
                            SamInfo::getInstance().clearChromosomeInfo();
                        }
                        return -1;
                    }
                }
            } else if (linePos.size() == 10) {  // Quality value is the 11th field
                uint32_t qualityLen = i - linePos.at(9) - 1;
                if (baseLen != qualityLen) {
                    /* A QUAL of a single '*' denotes missing quality; its length need not equal the SEQ length */
                    const bool missingQual = (qualityLen == 1 &&
                                              line.at(linePos.at(9) + 1) == '*');
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
                char ch = line[i];
                if (ch == 'N' || ch == 'n') {
                    baseNCount++;
                }
            } else if (linePos.size() == 10) {
                char ch = line[i];
                if (ch >= 0) {
                    qualityFrequnce[(uint8_t)ch].second++;
                }
            }
        }
    }
    // Each line must have at least 10 tabs, less than 10 means not a SAM file
    if (linePos.size() < 10) {
        LOG_ERROR("Sam field line invalid, size = %d", linePos.size());
        if (headEndLine >= 0) {
            SamInfo::getInstance().clearChromosomeInfo();
        }
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

    /*
     * A SAM block's line positions are recorded by the format's own reader (SamBlockReader /
     * BamGzBlockReader): they split the input on newlines and keep every record boundary, so a
     * non-empty block with no line position at all is not SAM/BAM text. The pbgz-input path is
     * the case that matters here - it hands over the archive's already-coded payloads (see
     * PbgzBlockReader::readBlock) - but any such block must not be coded as text.
     *
     * Failing here is what makes the bytes survive: the caller sees it and swaps in a
     * BinaryCompressActuator (CompressEngine::actuatorPreProc), which stores the block verbatim.
     * Returning 0 instead would let compress() take its "header only" branch and write an empty
     * block, dropping the whole block's data without a word.
     */
    if (inBlockPtr->getDataLen() > 0 && inBlockPtr->getNpos().empty()) {
        LOG_ERROR("SAM block(%ld) carries %ld bytes but no line positions, not SAM/BAM text.",
                  (long)inBlockPtr->getBlockId(), (long)inBlockPtr->getDataLen());
        return -1;
    }

    /*
     * The SEQ column's file-level decisions (its match coder, and the form the N positions travel
     * in) are not taken here: they are the codec pre-selection's, which measures them on the first
     * block before any actuator exists (see CodecSelector::selectSeqReferenceCoder). This pass only
     * fills the column stats below, and the consumption side reads the pinned verdicts
     * (seqMatchCoderFor, and the N-form choice inside writeSeqExceptionStreams).
     */
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Initialize quality frequency analysis similar to FastqActuator
    std::pair<uint8_t, uint32_t> qualityFrequnce[256];
    for (int32_t i = 0; i < 256; ++i) {
        qualityFrequnce[i].first = i;
        qualityFrequnce[i].second = 0;
    }

    // Parse SAM file header
    for (uint32_t idx = 0; idx < lineNum; ++idx) {
        uint32_t begin = (idx == 0 ? idx : npos[idx - 1]  + 1);
        uint32_t end = (idx < npos.size() ? npos[idx] : strlen((char*)buffer + begin));
        if (begin >= end) {
            return -1;
        }
        std::string line((char*)buffer + begin, end - begin);
        if (buffer[begin] == '@') {
            if (parseHeaderLine(line) != 0) {
                return -1;
            }
        } else {
            if (scanDataLine(line, idx, qualityFrequnce) != 0) {
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

    /*
     * Grow the output block to the real size of this input block before any
     * stream is written. The engine pool hands out output blocks whose initial
     * buffer only follows -l; with long reads a SAM/BAM data block (bounded by
     * the read-count/base-count slicer in BlockReader) can be far larger than
     * that, and a fixed -l-sized coder_io would overflow and abort compression.
     * So -l is only an allocation hint here, never a cap on the encoded block.
     * The formula mirrors the decompression side (see SamCodecActuator:decompress).
     */
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

    /*
     * The same guard as preAnalysis(), for the callers that compress without running the
     * pre-analysis first: a block with bytes but no line position cannot be "header only" -
     * taking that branch would write an empty block and drop the data. The pipeline always
     * runs preAnalysis() first (CompressEngine::actuatorPreProc), where the block is turned
     * into a binary one instead, so this is the belt to that branch's braces: fail loudly.
     */
    if (inBlockPtr->getDataLen() > 0 && inBlockPtr->getNpos().empty()) {
        LOG_ERROR("SAM block(%ld) carries %ld bytes but no line positions, cannot code as text.",
                  (long)inBlockPtr->getBlockId(), (long)inBlockPtr->getDataLen());
        return -1;
    }

    if (inBlockPtr->getNpos().size() <= (size_t)headEndLine) {
        LOG_DEBUG("SAM head only");
        // For header-only case, we still need to write the metadata
        // Calculate MD5 of original data (header only)
        std::string md5;
        if (inBlockPtr->hasMd5()) {
            /* Already hashed by the reader thread (CompressEngine::preDispatchBlock). */
            md5 = inBlockPtr->getMd5();
        } else {
            calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
        }
        meta["md5"] = md5;

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
    std::shared_ptr<coder_bwt_cm> headerCoder = std::make_shared<coder_bwt_cm>(headerIo.get());
    CoderFactory::applyLevel(headerIo.get(), CoderType::BWT_CM, engineCompressLevel());

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

/*
 * The QNAME column (field 0), written in whichever of its two layouts the codec pre-selection
 * picked for this file: affix segmentation, or coder_qname's cross-line deduplication. Both are
 * encoders, so the choice is made by running them (see CodecSelector::selectQnameLayout), and it
 * is a property of the file's naming scheme rather than of a block - which is why it is taken
 * before any block is written and not here. This reads that verdict (PreprocessInfo::qnameUseAffix)
 * and encodes with it; a column with no split to work from is written whole instead.
 *
 * Returns the column's encoded size, negative on failure.
 */
int32_t SamCodecActuator::encodeQnameColumn(uint32_t& fieldSrcLen, Json::Value& fieldMeta)
{
    // ID field: compress based on analysis result
    if (idAnalysis.posLength == UINT32_MAX) {
        LOG_DEBUG("Id will compress in all.");
        return compressIdFieldInAll(fieldSrcLen, fieldMeta);
    }

    /*
     * Which of the two layouts writes this column is a file-level verdict, and it arrives from the
     * codec pre-selection, which ran both encoders on a sample of the first block before any block
     * was dispatched (see CodecSelector::selectQnameLayout). Nothing is measured here: writing a
     * block is not a place to be comparing encoders.
     *
     * A file the pre-selection could not judge - the ID analysis above failed on the deciding
     * block, or the actuator is driven without an engine - keeps the affix layout, which is the
     * segmentation this column has always been written with.
     */
    const PreprocessInfo* preInfo = preprocessInfoMut();
    const int32_t verdict =
        (preInfo != nullptr) ? preInfo->qnameUseAffix.load(std::memory_order_relaxed) : -1;
    const bool useAffix = (verdict != 0);
    LOG_DEBUG("QNAME: file verdict %d -> %s", verdict, useAffix ? "affix" : "qname");
    return useAffix ? compressIdFieldSplit(fieldSrcLen, fieldMeta)
                    : compressIdFieldQname(fieldSrcLen, fieldMeta);
}

/*
 * Close the block: write the SAM meta (line and field counts, the two length totals, the exact
 * decoded text size the decompressor sizes its buffer from, then the per-field stream metas), hash
 * the block and encode the whole meta into the meta buffer.
 *
 * Returns 0 on success; the caller has already set the block id and type.
 */
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
    /* Exact decoded-text size of this block, so the decompressor can size its
       output buffer from the real block length instead of guessing from -l or
       from the compressed payload length (which under-allocates for large
       long-read blocks and caused out-of-bounds writes on decode). */
    samMeta["textlen"] = (Json::Value::UInt64)inBlockPtr->getDataLen();
    samMeta["streams"] = streamMeta;
    meta["sam"] = samMeta;

    PBGZ_PROF_SCOPE(pbgzprof::FIELD_META);

    // Calculate MD5 of original data
    std::string md5;
    if (inBlockPtr->hasMd5()) {
        /* Already hashed by the reader thread (CompressEngine::preDispatchBlock). */
        md5 = inBlockPtr->getMd5();
    } else {
        calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
    }
    meta["md5"] = md5;

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
                /*
                 * POS is compressed as unsigned varint (LEB128) deltas against the
                 * previous line's POS; the baseline resets at each chromosome
                 * switch, so deltas stay small and non-negative. The underlying
                 * entropy coder is chosen by preprocessing on the delta-varint
                 * stream itself (bwt_cm vs order-0 arithmetic, see
                 * CodecSelector::selectPosDeltaCoder); bwt_cm is the fallback.
                 */
                if (pickedCoderFor(fieldIdx, CoderType::BWT_CM) == CoderType::ARITH) {
                    fieldDstLen = compressPosFieldDelta<coder_arith>(fieldIdx, fieldSrcLen, fieldMeta);
                } else {
                    fieldDstLen = compressPosFieldDelta<coder_bwt_cm>(fieldIdx, fieldSrcLen, fieldMeta);
                }
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
                /*
                 * PNEXT is compressed as the delta against POS. The deltas of
                 * consecutive lines are far smaller than the raw coordinates;
                 * empirically bwt_cm compresses the delta text to ~0.86
                 * bytes/line, far better than compressing raw PNEXT or affix.
                 * The trial-based selection is based on raw PNEXT text, which
                 * does not match the delta encoding actually used here, so this
                 * field always takes the delta path and is excluded from coder
                 * selection.
                 */
                fieldDstLen = compressPNextFieldDelta<coder_bwt_cm>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 8: // TLEN
                /*
                 * TLEN is not stored verbatim: it is reconstructed from
                 * POS/PNEXT/CIGAR on decompression, storing exceptions only for
                 * lines that cannot be reconstructed. CIGAR/FLAG/RNAME/RNEXT/
                 * POS/PNEXT are all compressed before this field, so the
                 * tracking maps needed for reconstruction are already populated.
                 */
                fieldDstLen = compressTLen<coder_bwt_cm>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 9: // SEQ
                if (pRefeGene == nullptr) {
                    LOG_DEBUG("Base will compress without reference");
                    fieldDstLen = compressBaseWithoutRef(fieldIdx, fieldSrcLen, fieldMeta);
                } else {
                    LOG_DEBUG("Base will compress with reference");
                    fieldDstLen = compressBaseWithRef(fieldIdx, fieldSrcLen, fieldMeta);
                }
                break;
            case 10: // QUAL
                fieldDstLen = compressQuality(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 11: // Optional fields (all tags)
                fieldDstLen = compressRegularField(fieldIdx, fieldSrcLen, fieldMeta);
                break;
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

/*
 * Where one field of one content line sits.
 *
 * `line` is the start of the line, `tabs` the tab positions of that line (contentPos for this
 * record) and `fieldIdx` the 1-based SAM column; `lineEnd` is the line's length, used for the
 * last field, which no tab follows.
 *
 * The two lengths are the point of this struct. The column walks in this file disagree on
 * whether a field includes its trailing tab - the bytes handed to a coder usually do, since
 * the tab is what the decoder splits on, while a length or a hash usually does not - and
 * hiding that behind one `length` is exactly how the two conventions drifted apart at the
 * dozen sites that used to compute this inline.
 */
struct SamFieldRange {
    uint8_t* start = nullptr;
    /* The delimiter that ends the field: the next tab, or the end of the line for the last. */
    uint8_t* tab = nullptr;

    uint32_t lengthWithoutTab() const { return (uint32_t)(tab - start); }
    uint32_t lengthWithTab() const { return (uint32_t)(tab - start) + 1; }
};

static SamFieldRange samFieldRange(uint8_t* line, const std::vector<int64_t>& tabs,
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

/*
 * Which coder codes the SEQ match stream: BWT_CM or FC, and the smallest encoded sample wins
 * (see CodecSelector::trialSeqMatchCoder).
 *
 * The question is a file-level one, asked once by the codec pre-selection before any block is
 * dispatched and only read here. There is no per-block trial: that existed for coder_lzma, whose
 * per-block variation was worth bytes, and it cost a sample encode of every block - measured at -l8
 * on ERR11436629 it was 43% of the run time for a byte-identical archive. With coder_lzma gone from
 * the candidates (it wins on size but costs ~17x the time, see its registry row) there is nothing
 * left that varies enough per block to pay for that.
 *
 * A caller driving the actuator without an engine has no verdict to read, and a run without a
 * reference never got one (the column is not written as the match stream then): both fall back to
 * coder_bwt_cm, which is what this stream was coded with before the trial existed.
 */
CoderType SamCodecActuator::seqMatchCoderFor()
{
    PreprocessInfo* preInfo = preprocessInfoMut();
    if (preInfo != nullptr) {
        const int32_t pinned = preInfo->seqMatchCoder.load(std::memory_order_relaxed);
        if (pinned >= 0) {
            return (CoderType)pinned;
        }
    }

    /*
     * No verdict covers this file: the pre-selection found no SEQ payload to measure, or a caller
     * is driving the actuator without an engine. coder_bwt_cm is the coder this stream used before
     * the trial existed.
     */
    return CoderType::BWT_CM;
}

/*
 * Encode PNEXT's exception stream: the (contentIdx, delta) pairs of the records that cannot be
 * rebuilt from their mate.
 *
 * The contentIdx column goes out as forward deltas and the delta column as zigzag varints, and
 * where every delta is the same value that value moves into the meta and the column is dropped
 * outright: a block whose records are all unpaired - or all unmapped, like ERR14949932 - has
 * every line here with the same delta, and there the absolute layout this replaced spent
 * 937,737 B of stream on 3,343,586 entries whose delta is uniformly 0.
 *
 * contentIdx is the block-internal data-line index (small, increasing); delta = pnext - pos is
 * signed (the mate offset, negative when the record had no real PNEXT), so it goes through
 * zigzag. Records in fieldMeta how the stream was written, and returns its encoded size.
 */
template<typename CoderType>
uint32_t SamCodecActuator::writePnextExceptions(const std::vector<std::pair<uint32_t, int64_t>>& exc,
                                                  Json::Value& fieldMeta, Json::Value& metaStreams,
                                                  uint32_t& totalDstLen)
{
    Json::Value metaSubs;
    if (!exc.empty()) {
        std::shared_ptr<coder_io> excIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "PNEXT exceptions");
        std::shared_ptr<CoderType> excCoder = std::make_shared<CoderType>(excIo.get());

        /*
         * Exception stream encoding: the contentIdx column as forward deltas and the delta
         * column as zigzag varints, and where every delta is the same value that value moves
         * into the meta and the column is dropped outright. A block whose records are all
         * unpaired - or all unmapped, like ERR14949932 - has every line in this stream with
         * the same delta, and there the old layout (absolute contentIdx then delta, 1-3 B + 1 B
         * per entry) spent 937,737 B of stream on 3,343,586 entries whose delta is uniformly 0.
         *
         * contentIdx is the block-internal data-line index (small, increasing); delta =
         * pnext - pos is signed (mate offset, negative when the record had no real PNEXT), so
         * it goes through zigzag.
         */
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

template<typename CoderType>
int32_t SamCodecActuator::compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    fieldSrcLen = 0;
    uint32_t totalSrcLen = 0;

    /*
     * PNEXT is stored as "exception-only": for a paired record whose mate is
     * present and mutually mapped (same QNAME), the PNEXT can be rebuilt on the
     * decoder from the mate's POS, so it need not be stored. Only records that
     * cannot be rebuilt (unpaired / mate unmapped / supplementary / mate
     * missing from the block / asymmetric coords) are stored as exceptions.
     *
     * The decoder rebuilds non-exception PNEXT by grouping records with the
     * same QNAME (a QNAME with exactly two mutually-mapped records is a mate
     * pair). This preserves the block line order, so POS delta coding and the
     * output order are unaffected.
     */

    // Pass 1: collect (lineIdx -> qname, pos, pnext, flag) for every data line.
    struct RecInfo {
        std::string qname;
        int64_t pos = 0;
        int64_t pnext = 0;
        uint16_t flag = 0;
        bool valid = false;   // paired (0x1) and mate mapped (not 0x8)
        bool hasPnext = false; // pnext != 0/ *
    };
    std::map<uint32_t, RecInfo> records;
    std::unordered_map<std::string, std::vector<uint32_t>> qnameToLines;

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
        RecInfo ri;
        ri.qname.assign((char*)line, qnameLen);

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

        records[lineIdx] = std::move(ri);
        // Group by QNAME using the same criterion as the decoder (valid =
        // paired && mate mapped), so exception decisions match exactly.
        if (records[lineIdx].valid) {
            qnameToLines[records[lineIdx].qname].push_back(lineIdx);
        }
    }

    // Pass 2: decide rebuildable vs exception.
    std::vector<std::pair<uint32_t, int64_t>> pnextExceptions; // (contentIdx, delta = pnext - pos)
    for (auto& kv : records) {
        uint32_t lineIdx = kv.first;
        RecInfo& ri = kv.second;

        bool rebuildable = false;
        if (ri.valid && ri.hasPnext) {
            const auto& mates = qnameToLines[ri.qname];
            // For the decoder to uniquely locate the mate from the QNAME group,
            // this QNAME must contain exactly two mutually-mapped records.
            if (mates.size() == 2) {
                for (uint32_t ml : mates) {
                    if (ml == lineIdx) continue;
                    const RecInfo& mr = records[ml];
                    if (mr.valid && mr.hasPnext &&
                        mr.pnext == ri.pos && ri.pnext == mr.pos) {
                        rebuildable = true;
                        break;
                    }
                }
            }
        }
        if (!rebuildable) {
            uint32_t contentIdx = lineIdx - headEndLine;
            /* Reconstruction on the decoder is pnext = delta + pos. When the record
               had no real PNEXT (field was "0" or "*"), hasPnext is false and the
               original value is 0, so delta must be -pos (not 0) for the decoder to
               reproduce 0. */
            int64_t delta = ri.hasPnext ? (ri.pnext - ri.pos) : (0 - (int64_t)ri.pos);
            pnextExceptions.emplace_back(contentIdx, delta);
        }
    }

    // Encode the exception stream: (contentIdx, delta) pairs, like TLEN.
    Json::Value metaStreams;
    uint32_t totalDstLen = 0;

    /* (see writePnextExceptions) */
    totalDstLen = writePnextExceptions<CoderType>(pnextExceptions, fieldMeta, metaStreams, totalDstLen);

    // fieldSrcLen (output) must be the original PNEXT field byte length (with
    // trailing tab) summed over all lines, so the -s statistics and
    // recordFieldStats print the true source length. totalSrcLen accumulates
    // exactly that (the exception-stream size is recorded only in the sub-stream
    // meta, not folded into the field source length).
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

template<typename CoderType>
int32_t SamCodecActuator::compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "POS delta");
    std::shared_ptr<CoderType> fieldCoder = std::make_shared<CoderType>(fieldIo.get());
    /* The arithmetic backend consumes the file-level delta prior to skip each
     * block's model cold start. bwt_cm (the other CoderType this template is
     * instantiated with) ignores it: dynamic_cast yields null. The prior was
     * produced by preprocessing only when arith won the POS trial, so it is
     * consistent with this instance being coder_arith. */
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

        /*
         * Detect chromosome switches via the per-line chromosome index that
         * compressChrName (field 2) already recorded in mappedChr. No reset
         * index list is stored in the metadata: the decoder replays the same
         * comparison against its own decoded mappedChr, which keeps the
         * bitstream independent of RNAME ordering (CRAM-style).
         */
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
        /*
         * Delta is encoded as an unsigned varint (LEB128). After the baseline
         * is reset at every RNAME change, coordinate-sorted deltas are
         * non-negative, so the sign bit never needs to be stored.
         */
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
            /*
             * RNEXT: the mate's reference is "=" (same as RNAME) in ~99.6% of
             * paired reads, so a fixed uint16 per line wastes a byte of zeros on
             * nearly every row. Store a one-byte code instead:
             *   0x00 = "=" (0xFFFE)
             *   0x01 = "*" (0xFFFF)
             *   0x02..0xFD = chrIndex + 2   (real contig index 0..251)
             *   0xFF = escape, followed by the raw uint16 (index >= 252,
             *          extremely rare; requires > 250 reference contigs)
             * 0xFE is unused.
             */
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
            uint8_t ch = '\n';
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

            /*
             * When FLAG/POS/PNEXT are compressed in textual form (affix), the
             * tracking maps must still be filled here: compressBaseWithRef
             * relies on mappedPos/mappedChr/mappedFlag to restore reference
             * positions in the SEQ phase, and the decompression side does too.
             * The binary form fills them in compressNumber; the textual form
             * fills them here, keeping both sides consistent.
             */
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

/*
 * The OPTION field (all tags starting at column 12) is compressed by CRAM-style
 * tag columnization.
 *
 * Idea: OPTION text = constant tag structure (name/colon/type) + varying
 * values. Generic byte-stream compression (affix/bwt) absorbs the tag structure
 * via prefix matching, but values are still stored as decimal text, wasting
 * roughly half. Here the structure and the values are separated:
 *   1. Build a per-block tag dictionary (each name+type stored once).
 *   2. Each record stores only its sequence of tag ids (in order).
 *   3. Each tag's values form a separate column: integers as delta +
 *      fixed-width binary, the rest as length-prefixed byte streams.
 *   4. Each column is compressed independently with bwt_cm.
 */
/*
 * Parse the OPTION column line by line: build the tag dictionary (name and type, in first-seen
 * order), each line's tag sequence, and each tag's own column of values. Tags are the
 * "TAG:TYPE:VALUE" triples of a SAM line's trailing fields, and a line that carries none
 * contributes an empty sequence.
 *
 * Returns 0, or -1 when a line has no OPTION field at all.
 */
int32_t SamCodecActuator::collectOptionTags(OptionParseState& st)
{
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        uint32_t lineEnd = npos[lineIdx] - lineStart;
        uint8_t* line = buffer + lineStart;
        if (*line == '@') {
            continue;
        }
        uint32_t contentIdx = lineIdx - headEndLine;
        uint32_t fieldStartPos = contentPos[contentIdx][10] + 1;
        uint32_t fieldLen = lineEnd - fieldStartPos;
        if (fieldLen == 0) {
            st.recIds.push_back({});
            continue;
        }
        st.srcLen += fieldLen;

        std::vector<uint8_t> ids;
        const uint8_t* p = line + fieldStartPos;
        uint32_t segStart = 0;
        for (uint32_t i = 0; i <= fieldLen; ++i) {
            bool isEnd = (i == fieldLen) || (p[i] == '\t');
            if (isEnd) {
                uint32_t segLen = i - segStart;
                if (segLen > 0 && p[segStart] != '\n') {
                    const uint8_t* ts = p + segStart;
                    const uint8_t* colon1 = nullptr, *colon2 = nullptr;
                    for (uint32_t j = 1; j < segLen; ++j) {
                        if (ts[j] == ':' && colon1 == nullptr) colon1 = ts + j;
                        else if (ts[j] == ':' && colon2 == nullptr) colon2 = ts + j;
                    }
                    if (colon2 != nullptr) {
                        std::string name((char*)ts, colon1 - ts);
                        std::string type((char*)(colon1 + 1), colon2 - colon1 - 1);
                        std::string value((char*)(colon2 + 1), ts + segLen - colon2 - 1);
                        int id;
                        auto it = st.tagId.find(name);
                        if (it == st.tagId.end()) {
                            id = (int)st.tagDict.size();
                            st.tagDict.push_back({name, type});
                            st.tagId[name] = id;
                            st.tagVals.push_back({});
                        } else {
                            id = it->second;
                        }
                        ids.push_back((uint8_t)id);
                        st.tagVals[id].push_back(value);
                    }
                }
                segStart = i + 1;
            }
        }
        st.recIds.push_back(ids);
    }
    return 0;
}

int32_t SamCodecActuator::compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {

    /* Parse OPTION line by line: see collectOptionTags. */
    OptionParseState opt;
    if (collectOptionTags(opt) != 0) {
        return -1;
    }
    fieldSrcLen = opt.srcLen;

    const uint32_t nTag = (uint32_t)opt.tagDict.size();
    if (nTag == 0) {
        fieldMeta["mode"] = "tag_split";
        fieldMeta["srclen"] = 0;
        fieldMeta["dstlen"] = 0;
        fieldMeta["tags"] = Json::Value(Json::arrayValue);
        return 0;
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

    /*
     * Every column is encoded as "line by line + delimiter" and compressed with
     * bwt_cm's line-wise mode. A single large blob must not be used: bwt_cm
     * fails to decode near-constant large blobs (id sequences almost constant,
     * integer columns with deltas almost all 0), whereas line-wise encoding with
     * '\n' delimiters is a verified working usage.
     */

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
        std::shared_ptr<coder_bwt_cm> c = std::make_shared<coder_bwt_cm>(io.get());
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

    /* Each tag value column: one value per line, terminated by '\n' (SAM values never contain '\n'). Integer columns are stored as text to avoid the bwt large-blob defect. */
    for (uint32_t t = 0; t < nTag; ++t) {
        const std::vector<std::string>& vals = opt.tagVals[t];
        std::vector<uint8_t> col;
        for (size_t k = 0; k < vals.size(); ++k) {
            col.insert(col.end(), vals[k].begin(), vals[k].end());
            col.push_back('\n');
        }
        std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "OPTION tag");
        std::shared_ptr<coder_bwt_cm> c = std::make_shared<coder_bwt_cm>(io.get());
        size_t pos = 0;
        for (size_t k = 0; k < vals.size(); ++k) {
            size_t start = pos;
            while (pos < col.size() && col[pos] != '\n') pos++;
            pos++;
            if (pos > start) c->encode_line(col.data() + start, (uint32_t)(pos - start));
        }
        c->encode_flush();
        if (io->err != coder_io::IO_OK) return -1;
        Json::Value s;
        s["sname"] = "tag";
        s["tag"] = opt.tagDict[t].first;
        s["type"] = opt.tagDict[t].second;
        s["srclen"] = (Json::Value::UInt)col.size();
        s["dstlen"] = (Json::Value::UInt)io->data_len;
        s["coder"] = io->meta;
        streams.append(s);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + io->data_len);
        totalDst += io->data_len;
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
            uint32_t sequeceLength = parseCigar(fieldStart, fieldLength);
            baseLengthBuffer[lineIdx - headEndLine] = sequeceLength;
            /* Reference span for TLEN reconstruction (field 8 is compressed after this field). */
            cigarReadLen[lineIdx] = parseCigarRefConsumed(fieldStart, fieldLength);
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

uint32_t SamCodecActuator::parseCigarRefConsumed(uint8_t* cigarString, uint32_t cigarLength) {
    /* The reference span the CIGAR consumes: only M/D/N/=/X count. The parse itself is the shared
       one (see sam_seq_payload.h), so the span and the operation list cannot disagree. */
    return cigarRefConsumed(cigarString, cigarLength);
}

/*
 * Reconstruct TLEN per the SAM spec:
 *   |TLEN| = rightmost mapped base - leftmost mapped base + 1
 * The left fragment is positive and the right one negative. When pos < pnext,
 * this read is on the left and the right end = pnext + mate reference span - 1;
 * when pos > pnext, this read is on the right, the left end = pos, and the
 * template length is negative.
 */
int32_t SamCodecActuator::computeTLEN(uint32_t lineIdx, bool minusOne) {
    /* Not paired (FLAG bit 0) or either end unmapped (bits 2/3): TLEN is set to 0. */
    auto flagIt = mappedFlag.find(lineIdx);
    if (flagIt == mappedFlag.end() || !(flagIt->second & 0x1) ||
        (flagIt->second & 0x4) || (flagIt->second & 0x8)) {
        return 0;
    }

    /* Reference sequence unavailable or mate on a different reference: TLEN is set to 0. */
    auto chrIt = mappedChr.find(lineIdx);
    if (chrIt == mappedChr.end() || chrIt->second == 0xFFFF) {
        return 0;
    }
    auto nextChrIt = nextMappedChr.find(lineIdx);
    if (nextChrIt == nextMappedChr.end() || nextChrIt->second == 0xFFFF) {
        return 0;
    }
    if (nextChrIt->second != 0xFFFE && nextChrIt->second != chrIt->second) {
        return 0;
    }

    auto posIt = mappedPos.find(lineIdx);
    if (posIt == mappedPos.end()) {
        return 0;
    }
    int64_t pos = posIt->second;

    auto pnextIt = nextMappedPos.find(lineIdx);
    if (pnextIt == nextMappedPos.end()) {
        return 0;
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

    int64_t templateLen;
    /* The rule and the two conventions it covers live in sam_field_rules.h, next to the BAM
       actuator that applies the same one to the same column. */
    templateLen = samTemplateLen(pos, pnext, refSpan, mateRefSpan, minusOne);
    return (int32_t)templateLen;
}

template<typename CoderType>
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

    /*
     * First pass: count how many records each template-length convention
     * matches (bwa minus 1 / minimap2 plain); the one with more matches is
     * chosen as this block's convention. The convention is stored in meta for
     * the decompression side.
     */
    uint64_t convMinus1 = 0, convPlain = 0;
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
            if (computeTLEN(lineIdx, true) == currentTLEN) convMinus1++;
            if (computeTLEN(lineIdx, false) == currentTLEN) convPlain++;
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
            int32_t computedTLEN = computeTLEN(lineIdx, minusOne);
            if (computedTLEN != currentTLEN) {
                tlenExceptions.push_back(std::make_pair(contentIdx, currentTLEN));
            }
        }
    }

    if (!tlenExceptions.empty()) {
        std::shared_ptr<coder_io> tlenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "TLEN exceptions");
        std::shared_ptr<CoderType> tlenCoder = std::make_shared<CoderType>(tlenIo.get());

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
            tlenExcSrcLen += tlenPutVarint(tlenExcBuffer + tlenExcSrcLen, lineIdx - prevLine);
            prevLine = lineIdx;
            tlenExcSrcLen += tlenPutVarint(tlenExcBuffer + tlenExcSrcLen,
                                           tlenZigzag32(tlenExceptions[i].second));
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
        /* Which layout the stream is in (see tlenDecodeVarints). Set here, where a stream exists,
           so a block without exceptions carries no marker. */
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
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    fieldSrcLen = 0;
    /* Pre-allocate by the actual block data length (SAM blocks may be split by read count and exceed the byte block_size). */
    const auto seqT0 = pbgzprof::nowOrZero();
    std::unique_ptr<uint8_t[]> tmpBuffer = std::make_unique<uint8_t[]>(inBlockPtr->getDataLen());
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
        const SamFieldRange field = samFieldRange(line, contentPos[contentIdx], fieldIdx, lineEnd);
        uint8_t* fieldStart = field.start;
        uint32_t fieldLength = field.lengthWithTab();

        if (minBaseLength == maxBaseLength) {
            // Remove trailing \t
            memcpy(tmpBuffer.get() + fieldSrcLen, fieldStart, fieldLength - 1);
            fieldSrcLen += fieldLength - 1;
        } else {
            // Encode the field data
            memcpy(tmpBuffer.get() + fieldSrcLen, fieldStart, fieldLength);
            fieldSrcLen += fieldLength;
        }
    }

    // Create encoder for regular field compression
    pbgzprof::addSince(pbgzprof::SEQ_PREP, seqT0);
    const auto seqT1 = pbgzprof::nowOrZero();
    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ");
    std::shared_ptr<coder> fieldCoder = makeFieldEncoder(fieldIdx, samFieldDefaultCoder(fieldIdx, CoderType::FC), fieldIo.get(), false);

    pbgzprof::addSince(pbgzprof::SEQ_CTOR, seqT1);
    const auto seqT2 = pbgzprof::nowOrZero();
    fieldCoder->encode_line(tmpBuffer.get(), fieldSrcLen);
    fieldCoder->encode_flush();
    pbgzprof::addSince(pbgzprof::SEQ_CODEC, seqT2);
    // Smart pointer automatically cleans up
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode SEQ field overflow: output buffer too small");
        return -1;
    }

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    // Set field metadata
    fieldMeta["minlen"] = minBaseLength;
    fieldMeta["maxlen"] = maxBaseLength;
    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_INFO("SAM base field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

/*
 * The SEQ match stream split into (run lengths, surviving values), or left whole when the
 * split would not pay for itself: a split with useRle false carries no buffers at all, and
 * the caller codes the match stream as it is.
 */
/*
 * Where one field of one content line sits.
 *
 * `line` is the start of the line, `tabs` the tab positions of that line (contentPos for this
 * record) and `fieldIdx` the 1-based SAM column; `lineEnd` is the line's length, used for the
 * last field, which no tab follows. The range excludes both surrounding tabs.
 *
 * Every column walk in this file asks this question - SEQ, QUAL, Nodifier, the fields
 * compressed as a whole - and it used to be spelled out inline at each of them.
 */
/* SamFieldRange and samFieldRange are defined near the top of this file, before their first
   user: every column walk below needs them. */

/*
 * Record that `ch` occurred at block position `pos`: the first occurrence registers the
 * character (and reserves its list, unless it is the 'N' class the constructor sized), every
 * occurrence appends the position. The other half of the SEQ column's exception handling is
 * the emitter, writeSeqExceptionStreams.
 */
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

/*
 * Append one record's payload to the block's match stream, refusing to run past `capacity`
 * (the caller reports that as a block error, with the lengths it already knows). What the
 * stream is then coded as - whole, or split into runs and values - is splitSeqMatchStream's
 * business.
 */
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

/*
 * Write the base-length auxiliary stream, for the layout where base lengths vary: the reads whose
 * length is not already implied by their CIGAR, as (ordinal, length) pairs.
 *
 * Both columns go out as forward deltas in varint, with the ordinals written as one run before the
 * lengths rather than interleaved with them: a file that is entirely unmapped - or entirely mapped
 * - collapses its whole ordinal column into the same repeated delta, and the lengths keep
 * statistics of their own. The absolute 4-byte pairs this replaced cost 8 B per read; on
 * ERR14949932 (all 3,343,586 reads unmapped) 26,748,688 B of source became 6,171,912 B, where
 * CRAM's read-length series needs 1,302,132 B for the same reads.
 *
 * Does nothing when every read has the same length.
 */
int32_t SamCodecActuator::writeSeqBaseLengthStream(Json::Value& metaSubs, Json::Value& metaStreams,
                                                  uint32_t& totalSrcLen, uint32_t& totalDstLen)
{
    if (minBaseLength != maxBaseLength) {
        std::shared_ptr<coder_io> lenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ length");
        /*
         * The reads whose length is not already implied by their CIGAR, as (ordinal, length)
         * pairs. Both columns go out as forward deltas in varint, with the ordinals written as
         * one run before the lengths rather than interleaved with them.
         */
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

        CoderFactory::applyLevel(lenIo.get(), CoderType::BWT_CM, engineCompressLevel());
        std::shared_ptr<coder_bwt_cm> lenCoder = std::make_shared<coder_bwt_cm>(lenIo.get());
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
        /* The delta layout, and how many entries each of the two runs holds; a stream without
           these carries the absolute pairs this replaced. */
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

/*
 * Write the SEQ match stream: sub-stream "m" (the run-length segment under RLE, the whole match
 * stream otherwise) and, under RLE, sub-stream "mval" (the surviving non-zero values).
 *
 * A pinned verdict, or one measured on a sample coder_fc could carry, can still name a segment too
 * short for it, so an FC verdict on such a segment falls back to BWT_CM. "mval" stays pinned to
 * BWT_CM either way; the decoder dispatches on the magic each sub-stream records.
 */
int32_t SamCodecActuator::writeSeqMatchStreams(const SeqRleSplit& rle, const uint8_t* matchBuffer,
                                               uint32_t matchLen, uint32_t srcLen, coder_io* matchIo,
                                               Json::Value& metaSubs, Json::Value& metaStreams,
                                               uint32_t& totalDstLen)
{
    {
        const uint32_t payLen = rle.useRle ? rle.runLength : matchLen;
        const uint8_t* payBuf = rle.useRle ? rle.run.get() : matchBuffer;

        CoderType payType = seqMatchCoderFor();
        /* A pinned verdict, or one measured on a sample coder_fc could carry, can still
           name a segment too short for it (FC_MIN_LEN bytes or less). */
        if (payType == CoderType::FC &&
            (payLen <= FC_MIN_LEN || payLen >= (uint32_t)FC_MAX_LEN)) {
            payType = CoderType::BWT_CM;
        }
        CoderFactory::applyLevel(matchIo, payType, engineCompressLevel());
        std::shared_ptr<coder> payCoder = CoderFactory::makeEncoder(payType, matchIo);
        payCoder->encode_line(const_cast<uint8_t*>(payBuf), payLen);
        payCoder->encode_flush();
        LOG_DEBUG("SEQ match stream: rle=%d coder=%s src=%u dst=%u", (int)rle.useRle,
                  coderTypeToMagic(payType), payLen, (uint32_t)matchIo->data_len);
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
        /*
         * orgrawlen = 展开后的原始长度（= 总碱基数），解压端据此分配并还原；
         * rlelen    = 本子流（游程段）的解码输出长度，等于上面写入的 srclen。
         */
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
        /* This half stays pinned to BWT_CM: only the run-length stream above follows
           the field selection. Decoding dispatches on the recorded magic either way. */
        {
            CoderFactory::applyLevel(valIo.get(), CoderType::BWT_CM, engineCompressLevel());
            std::shared_ptr<coder_bwt_cm> bwt = std::make_shared<coder_bwt_cm>(valIo.get());
            bwt->encode_line(rle.val.get(), rle.valLength);
            bwt->encode_flush();
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

int32_t SamCodecActuator::compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    if (pRefeGene == nullptr) {
        LOG_ERROR("Reference genome is not available for base compression with reference");
        return -1;
    }

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    /*
     * SEQ characters the 2-bit path cannot carry, one accumulator per character (see
     * SeqExceptionClass), with excPresent listing the ones that occurred. The positions
     * increase strictly within the block, so the deltas are small - a file with ~11% Ns
     * averages ~5, i.e. one byte each - and unlike the absolute form their statistics do not
     * widen with the block: on ERR14949932 the absolute form of the 'N' class compresses
     * 14.7% worse at -l8 (100k reads per block) than at -l5 (10k).
     *
     * 'N' is the only character whose absolute form is a candidate, so it is the only one
     * given the size estimate up front; the others are small and grow as they are filled.
     */
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
    /* Scratch for the reference 2-bit-per-base sequence used by the CIGAR
       M/=/X segment rebuild (a single op never exceeds baseMaxLength bases).
       getStretch2Bits1Char can write up to outLen + 3 bytes because of its
       unaligned 4-byte writes, so size the buffer with extra slack. */
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
    /*
     * SEQ match 流（逐碱基 2bit：M/=/X 段为与参考的 XOR，I/S 段为原值）。
     *
     * 对齐 FastqCodecActuator::compressBaseWithRef 的做法（src/fastq_actuator.cpp:954）：
     *   1) 先把整块的 2bit 结果累积到 matchBuffer，拿到总长度 matchLen；
     *   2) 先攒够再整体编码（见下面对 matchCoderType 的说明），因为 coder_fc 的
     *      decode_line 只支持整块解压（"only support block decompress"），不能逐行。
     */
    std::unique_ptr<uint8_t[]> matchBuffer = std::make_unique<uint8_t[]>(inBlockPtr->getDataLen() + 1);
    uint32_t matchLen = 0;
    int64_t srcLen = 0;
    uint32_t totalBaseLength = 0;

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

        /*
         * CIGAR-segment-based reference rebuild: for a mapped read the SEQ is encoded per-base as a
         * 2-bit XOR against the reference, but only on the operations that consume reference
         * (M/=/X); operations that consume SEQ but not reference (I/S) are stored directly, and
         * reference-only operations (D/N) only advance the reference position without producing
         * any SEQ bytes. This is required so that indels do not cause the read to be compared
         * against a contiguous reference window it does not actually align to, and the
         * decompression side mirrors the same walk (see buildSeqReferenceCodedBases).
         *
         * The branch itself is the shared one (buildSeqRecordPayload), because the codec
         * pre-selection builds the same payload to trial the match coder on it; what stays here is
         * what only the block pass does with it - the reference's own match statistics and the
         * exception characters.
         */
        static const std::vector<CigarOp> emptyCigarOps;
        /* The reference span this record's CIGAR consumes, from the parse the caller did. */
        const auto crlIt = cigarReadLen.find(lineIdx);
        const uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
        const std::vector<CigarOp>& ops =
            (contentIdx < cigarOpList.size()) ? cigarOpList[contentIdx] : emptyCigarOps;

        const SeqRecordPayload payload =
            buildSeqRecordPayload(chrId, flag, startPos, refConsumed, ops, seqStart, seqLength,
                                  pRefeGene, baseMappedBuffer.get(), ref2bitBuf.get());
        const uint32_t outLen = payload.length;
        const bool rawBases = payload.rawBases;
        if (payload.usedReference) {
            pRefeGene->updateMatchedGene((uint64_t)payload.refPos,
                (crlIt != cigarReadLen.end()) ? crlIt->second : seqLength);
        }
        /*
         * Record the characters the payload cannot carry, one list per character (see
         * SeqExceptionClass).
         *
         * A record coded against the reference has 2-bit codes in the payload, so A/C/G/T ride
         * there and every other spelling - N, lower case, IUPAC codes - is recorded instead of
         * being silently folded, the case-insensitivity of that code being exactly the problem.
         * A record that writes its own characters carries the plain bases and N itself; the
         * payload is still a stream of a DNA-oriented coder, though, and it does not hand back
         * every byte it is given, so anything outside A/C/G/T/N is recorded here as well.
         */
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
    /*
     * RLE preprocessing, applied only when the zero runs are long enough to pay for it.
     *
     * The match stream is a sparse 0..3 byte stream in which a block whose reads all take their
     * bases from the reference is ~99.28% zeros. Splitting it into two independent sub-streams
     *   "m"    - varint-encoded run lengths of the zero runs
     *   "mval" - the surviving non-zero values (1 byte each)
     * lets each be modelled on its own instead of forcing one context model to cope with two very
     * different distributions.
     *
     * That split is a loss whenever the zeros do not come in long runs, which is exactly the block
     * where no read is aligned to the reference: every base is then coded from its own value, the
     * zeros are just the 'A's (A codes as 0), and their runs are short. Measured on two all-unmapped
     * datasets - an ONT block (3995 reads, 56,570,656 bases, 30.6% 'A', average zero run 0.44) and a
     * short-read file (64 bp reads, 43.1% 'A', average run 0.76) - the split needs 3.4% and 19.2%
     * more for the SEQ column than the same stream coded whole. The reason it is worse rather than
     * merely useless is that the split destroys what the modeller exploits: each non-zero base
     * becomes one byte of "mval" plus about one byte of "m" (so the input the coders see grows by
     * 14%-39%), and the two streams are no longer adjacent, which costs most where the data has
     * cross-read structure to begin with (the short-read file, whose reads repeat each other).
     *
     * So the split is conditional, on how much of the block actually came from the reference: it is
     * written only when the match stream is at least 98% zeros. Where the crossover sits was
     * measured on 90 bp blocks built from con_sorted with a known fraction of reads marked unmapped
     * (same block, same coder, only the layout differing):
     *
     *   zeros  reads unmapped   SEQ split   SEQ whole
     *   99.0%      0%             27,018      27,649    the split wins, by 2.3%
     *   95.6%      5%             94,024      83,938    it loses, by 12%
     *   92.1%     10%            134,213     111,138
     *   85.2%     20%            198,046     150,558
     *   78.3%     30%            252,684     184,408
     *
     * The crossover is sharp and sits between 99.0% and 95.6%; 98% is the middle of that window,
     * with the two kinds of block production actually hands it well away from the line - a block
     * whose reads all take their bases from the reference is 99.2% zeros, one where none do is
     * 30-36%. So the rule is a threshold rather than a measurement: the two regimes are three
     * orders of magnitude apart in this statistic, and there is nothing in between to compare.
     *
     * The decision is per block, which is what lets one file hold both kinds of block - a reference
     * used for part of a file and not the rest - and each be written the better way.
     *
     * The decoding side needs no change and no new format version: it dispatches per sub-stream on
     * whether the block's SEQ meta carries the "rle" member (see initDecoder / decompressBase), and
     * a stream written without the split is byte-for-byte the form used before the split existed.
     */
    SeqRleSplit rle = splitSeqMatchStream(matchBuffer.get(), matchLen);

    /*
     * Sub-stream "m": the run-length segment under RLE, otherwise the whole match
     * stream. Which coder writes it was decided once per file, before any block ran, by
     * the codec pre-selection (see CodecSelector::selectSeqReferenceCoder / seqMatchCoderFor);
     * here that verdict is applied, with the FC length guard above. The decoder dispatches on
     * the magic the sub-stream's own meta records, so an archive written under a different
     * verdict still decodes.
     *
     * What each coder is worth on this transformed layout, measured on
     * con_sorted.sam (1M reads, -l8, dual-stream RLE):
     *   coder_bwt_cm : SEQ 362,880 B
     *   coder_fc     : SEQ 440,989 B   (+17%; generic LZP+BWT+MTF text path)
     *   coder_arith  : SEQ 580,606 B   (+60%, and its whole-block decode fails)
     * so arith is not among the candidates. The stream is sparse (~99.45% zero)
     * with long runs and cross-read patterns, which needs context modeling.
     */
    /* Sub-stream "m": run-length segment under RLE, otherwise the whole match stream. */
    {
    /* Sub-stream "m": run-length segment under RLE, otherwise the whole match stream. */
    if (writeSeqMatchStreams(rle, matchBuffer.get(), matchLen, srcLen, matchIo.get(), metaSubs,
                             metaStreams, totalDstLen) != 0) {
        return -1;
    }
    totalSrcLen += srcLen;
    }


    /*
     * SEQ exception sub-streams: one per character that is actually present, each carrying
     * the positions of that character (see SeqExceptionClass). A file of plain A/C/G/T writes
     * none of them, and only the 'N' character has two forms to be measured against each
     * other. They are emitted in ascending character order, so the block layout follows the
     * set of characters rather than the order in which they happened to appear.
     *
     * Both halves live in writeSeqExceptionStreams / writeSeqExceptionStream.
     */
    if (writeSeqExceptionStreams(excClasses, excPresent, metaStreams, totalSrcLen, totalDstLen) != 0) {
        return -1;
    }

    if (writeSeqBaseLengthStream(metaSubs, metaStreams, totalSrcLen, totalDstLen) != 0) {
        return -1;
    }

    // Set metadata
    /* Count of the 'N' stream, which older archives read as their whole N list; every stream
       carries its own count in its stream meta. */
    fieldMeta["ncount"] = (Json::Value::UInt)excClasses[(uint8_t)'N'].count;
    fieldMeta["litbases"] = (Json::Value::UInt)1;
    /*
     * How to read a payload byte that is above the 2-bit range: as one of the record's own
     * characters, because the record does not use the reference (see the fallback in the record
     * loop). Written unconditionally - it describes the layout this build writes, not whether this
     * particular block held such a record - and a stream without it is read the old way, every byte
     * a 2-bit code with the exception streams applied over it.
     */
    fieldMeta["minlen"] = minBaseLength;
    fieldMeta["maxlen"] = maxBaseLength;
    fieldMeta["totalsrclen"] = totalSrcLen;
    fieldMeta["totaldstlen"] = totalDstLen;
    fieldMeta["streams"] = metaStreams;
    fieldMeta["field"] = fieldIdx;
    /*
     * Only fieldSrcLen is corrected: the SEQ column's raw size = the match
     * stream's source length (srcLen, i.e. the total SEQ text length across
     * records, matching the QUAL column's accounting); the npos/baselen
     * auxiliary streams are excluded. meta's totalsrclen is kept as is (it
     * includes the auxiliary streams).
     */
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
    CoderFactory::applyLevel(excIo.get(), CoderType::BWT_CM, engineCompressLevel());
    std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(excIo.get());
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

        /*
         * Only the 'N' class keeps more than one form, and the form arrives as a verdict from the
         * codec pre-selection, which measured the three candidates on the first block before any
         * block was dispatched (see CodecSelector::selectSeqReferenceCoder). Nothing is measured
         * here: writing a block is not a place to be comparing coders, and every block of the file
         * writes the one form the file's first block settled. A file the pre-selection could not
         * judge - its first block carried no N at all - keeps the absolute form, which is what this
         * stream held before the trial existed. The other characters are always in the delta form.
         */
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

/*
 * Encode the QUAL column record by record, through whichever coder the pre-selection's trial
 * picked (see qual_record_encoder) - all this step knows is how to hand one record over.
 *
 * Three details it does carry: a record whose quality is the single '*' for "missing" is expanded
 * into seqLen '*' first, or the length the decoder fetches by SEQ/CIGAR would not match the one
 * byte written here and the whole column would shift out of alignment (the decoder folds it back,
 * see decompressQuality); the SEQ column is passed along as context; and the strand direction is
 * read out of FLAG by hand (bit 0x10), since only coder_fcv2 wants it.
 *
 * Returns the column's source length - the quality text itself, without the auxiliary stream.
 */
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

        /*
         * Reads with missing quality (a single '*'): expand into seqLen '*'
         * before entering the stream; otherwise the per-record stream length (1)
         * would not match what the decompression side fetches by SEQ/CIGAR
         * length (seqLen), and the whole quality column would shift out of
         * alignment. The decode side folds it back into a single '*' (see
         * decompressQuality).
         */
        if (qualLength == 1 && qualStart[0] == '*' && seqLen > 0) {
            qualMissingBuf.assign(seqLen, '*');
            qualStart = qualMissingBuf.data();
            qualLength = seqLen;
        }

        /*
         * Take the strand direction from the FLAG field. Per the SAM spec, bit 0x10 being
         * set means SEQ and QUAL are stored relative to the reference's forward strand, i.e.
         * reversed relative to the order read out by the sequencer, and a coder that models
         * sequencing cycles must use this to restore the real cycle number. The field is
         * parsed by hand rather than with strtol because it is not NUL-terminated and this
         * is a hot path taken for every record.
         */
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

/*
 * Write the quality column's frequency table as an auxiliary stream: the alphabet in
 * frequency-descending order, two uint16 per entry (byte and weight), coded with BWT_CM.
 *
 * The decoder needs it to rebuild the coder's alphabet before the first record, and the stream
 * counts toward the field's totalsrclen (the original accounting) although the column's own source
 * length excludes it.
 */
int32_t SamCodecActuator::writeQualFreqStream(Json::Value& subMeta, Json::Value& streamMeta,
                                               uint32_t& totalSrcLength, uint32_t& totalDstLength)
{
    std::shared_ptr<coder_io> qualityFreqIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL freq table");
    std::shared_ptr<coder_bwt_cm> qualityFreqCoder = std::make_shared<coder_bwt_cm>(qualityFreqIo.get());
    CoderFactory::applyLevel(qualityFreqIo.get(), CoderType::BWT_CM, engineCompressLevel());
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

    /*
     * Quality values can use two coders. coder_qual is the original one and uses
     * SEQ as context; fcv2 is a context-mixing coder that uses the previous and
     * second-previous quality values, the in-record sequencing-cycle number, and
     * the strand direction as context.
     *
     * Which one is used is decided by the preprocessing trial results. The
     * quality column goes through a dedicated evaluation path (QualSelector),
     * whose two candidates are these coders; samples are collected per record
     * and therefore preserve record boundaries and strand direction.
     *
     * When preprocessing did not run, was skipped because the sample was too
     * small, or the current engine provides no preprocessing information,
     * coderFor returns the passed-in fallback QUAL, i.e. the original
     * coder_qual, behaving exactly as before selection was wired in.
     */
    int64_t qualPriorAddress = -1;
    CoderType pickedQualCoder = CoderType::QUAL;
    const PreprocessInfo* qualPreInfo =
        (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfo() : nullptr;
    if (qualPreInfo != nullptr) {
        pickedQualCoder = qualPreInfo->coderFor(SAM_QUAL, samFieldDefaultCoder(SAM_QUAL, CoderType::QUAL));
    }
    /*
     * Turning that verdict into an encoder - the quality alphabet, the fcv2 context tier
     * the trial settled on, the trained prior - is the factory's business, not this
     * actuator's (see CoderFactory::makeQualEncoder). All this side has to know is how to
     * feed a record, which is what qual_record_encoder offers whichever coder won.
     *
     * The prior's absolute address must be written into the block meta: the decoding side
     * can only match up if it obtains the same snapshot, and in random-access scenarios it
     * cannot infer the address by walking the sequential stream. It is registered only when
     * the snapshot was actually loaded, so the decoder is never promised one it cannot get.
     */
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

    /*
     * Only fieldSrcLen is corrected: the QUAL column's raw size = total length
     * of the quality text (streamSrcLen); the frequency-table auxiliary stream
     * is excluded. meta's totalsrclen is kept as is (it includes freqSrcLen).
     */
    fieldSrcLen = streamSrcLen;

    /*
     * In verbose mode, print per block which coder this block's QUAL actually
     * took, along with this block's own compression ratio.
     *
     * Why it must be done here: QUAL selection follows a "decided once in
     * preprocessing, reused by all subsequent blocks" pattern, but a single
     * compression run has ~10 threads compressing different blocks
     * concurrently, and the final log is a single summary line—looking at just
     * that line you cannot confirm each block really used the expected coder,
     * nor see the ratio jitter between blocks.
     *
     * Three constraints:
     *   1. Only output when verbose=true; otherwise keep the default output
     *      byte-identical;
     *   2. Each line must be emitted with a single fprintf call, because
     *      multiple threads write to stderr concurrently and splitting into
     *      several calls would interleave and become unreadable; a single
     *      fprintf write into the kernel buffer is usually atomic;
     *   3. When totalSrcLength is 0, set the ratio directly to 0 to avoid
     *      division by zero.
     */
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
    optionRecLines.clear();
    // Parse meta information
    initMetaInfo();

    /*
     * -m fast blocks carry a structured "bam" meta (columns instead of SAM
     * text). They are routed here because their block type is BAM like the
     * textual path's; hand them to the column decoder.
     */
    if (meta.isMember("bam")) {
        BamCodecActuator bamActuator(inBlockPtr, outBlockPtr, pbgzEngine, pRefeGene);
        int32_t ret = bamActuator.decompress();
        return ret;
    }

    /*
     * Pre-allocate at the block entry: the file header's block_size (the upper
     * bound fixed at compression time) x2 guarantees a large enough output
     * buffer—a deterministic value, not an estimate, unaffected by
     * fieldcount/read length. This is the primary defense against
     * out-of-bounds writes across all fields. block_size is read back from
     * baseFileMeta (DecompressEngine::createBlockReader) and is already
     * available here; when it is 0 (old files that did not write it), fall back
     * to the default getBlockSize(). coder_io's putc checks and the decode
     * error-return chain act as a backstop (see decompressQuality et al.).
     */
    size_t bs = pbgzEngine->getFileBlockSize();
    if (bs == 0) {
        bs = ConfigManager::getInstance().getBlockSizeByCompressLevel(pbgzEngine->getParameter().compressLevel);
    }
    /*
     * Pre-allocate the output buffer as the larger of "block-size upper bound
     * x2" and "this block's actual data length": SAM blocks may be split by
     * read count and exceed the byte block_size (e.g. with -l 1, 10000 reads
     * ~= 1.8MB > the 512KB block upper bound). inBlockPtr->getDataLen() comes
     * from the block meta's datalen, i.e. the original pre-compression length,
     * which is the output bound.
     */
    size_t outCapacity = bs * 2;
    if ((size_t)inBlockPtr->getDataLen() > outCapacity) {
        outCapacity = (size_t)inBlockPtr->getDataLen();
    }
    /* The block records its original decoded text length (textlen, written by
       the compressor): that is the exact output bound, independent of -l and of
       how large the compressed payload happens to be. */
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
        // Verify checksum of decompressed content
        std::string md5;
        calcMd5sum(md5, outBlockPtr->getBuffer(), outBlockPtr->getDataLen());
        if (md5 != meta["md5"].asString()) {
            LOG_ERROR("MD5 check failed for SAM data, blockid = %d, expected: %s, got: %s", outBlockPtr->getBlockId(),
                meta["md5"].asString().c_str(), md5.c_str());
            return -1;
        }
    }
    if (targeBlock != outBlockPtr) {
        MemoryUtil::safeDeleteClass(targeBlock);
    }

    return 0;
}

/*
 * Decode one line's fields, in order, appending each to the output block.
 *
 * The dispatch is per field, because the layout of a column is the stream's own business: ID has
 * its own split/coder machinery, SEQ and QUAL need the coder that wrote them, POS/PNEXT/CIGAR/TLEN
 * come from the pre-decode cache when preDecodeForTLEN filled it (see copyPreDecodedField below)
 * and fall back to their own decoders otherwise, and the trailing optional fields are ordinary
 * text with the block's final tab turned into a newline.
 *
 * SEQ's landing spot and length travel back in `st`, because QUAL is decoded against them.
 */

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
    /* +8 slack: getStretch2Bits1Char may write up to outLen+3 bytes due to its
       unaligned 4-byte writes, and outLen (= actualBaseLen) can reach maxBaseLength. */
    refeStrecchBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength + 8);
    uint32_t totalBaseLen = 0;
    /* The SEQ exception streams and their cursors live in seqExc (see SeqExcStream). */

    uint8_t* pBaseOut = nullptr;
    if (streams[9]["coder"]["magic"].asString() == "coder_fc") {
        /*
         * coder_fc is a "whole-block" coder: SEQ must be fully decoded in one
         * go, but the final SAM output is interleaved line by line
         * (ID\tFLAG\t...\tSEQ\tQUAL\n), so SEQ can only land somewhere else
         * first and be moved line by line afterwards. It is staged at the
         * **tail** of the outputBlock buffer (same landing spot as in
         * initDecoder), while the head appends line content normally; the two
         * grow toward each other without overlapping:
         *   head output <= block_size, tail SEQ <= block_size, and the block
         *   entry has already done a one-shot ensureCapacity at block_size*2,
         *   giving a deterministic capacity bound.
         * Why not a separate malloc'd buffer: the number of memcpy calls is
         * identical (lines must be moved either way); a separate buffer only
         * adds one malloc/free per block, page faults on first write, and an
         * extra memory peak of "threads x one full SEQ".
         *
         * Invariant: never realloc outputBlock within a block, or both the tail
         * pointer here and basePtr below would dangle.
         */
        pBaseOut = pBaseEnd - streams[9]["totalsrclen"].asUInt();
    }

    /*
     * First pre-decode POS/CIGAR/PNEXT: TLEN reconstruction needs the full mate
     * index, and during line-by-line decoding a mate may lie in the second half
     * of the block, not yet decoded. The pre-decoded results are cached per
     * line and the main loop copies them directly, avoiding a second decode.
     * POS delta baseline resets are detected from the decoded RNAME indices
     * (mappedChr), so no reset list is needed.
     */
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

    if (headerMeta["coder"]["magic"].asString() != "coder_bwt_cm") {
        return -1;
    }

    uint32_t dstLen = headerMeta["dstlen"].asUInt();
    if (refPosChrIndex == 65535) {  // not set position paramter
        // Create SAM file header decompressor
        std::shared_ptr<coder_io> headerIo = makeCoderIo(inBlockPtr->getBuffer(), dstLen, "SAM header");
        std::shared_ptr<coder_bwt_cm> headerDecoder = std::make_shared<coder_bwt_cm>(headerIo.get());

        // Set decoder level
        if (headerMeta["coder"].isMember("level")) {
            headerDecoder->set_level(headerMeta["coder"]["level"].asInt());
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

/*
 * A "cnt" sub-stream's deltas (see id_int::Counter) expanded into the zigzag varint form the numeric
 * layout carries, so that everything downstream of the buffer sees one delta per line either way.
 * Both initDecoder - which fills idNumericBufs - and preDecodeForTLEN, which rebuilds its own copy
 * because it must not consume the buffer the main loop reads from, need exactly this.
 */
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
        /* The numeric layout's own convention, which is what the per-line accumulation expects:
           the magnitude in the high bits, the sign in the low one. */
        const uint64_t mag = (dv >= 0) ? (uint64_t)dv : (uint64_t)(-dv);
        const uint32_t zz = (uint32_t)((dv >= 0) ? (mag << 1) : ((mag << 1) | 1u));
        if (wp + 5 > outCap) {
            return -1;
        }
        wp += tlenPutVarint(out + wp, zz);
    }
    return (int32_t)wp;
}

/*
 * Build the ID/QNAME column: one decoder per split sub-stream, plus the payload state each
 * segment layout needs (the reading side of collectIdSplitSegments, buildIdDictLayout,
 * buildIdHexLayout and buildIdValueLayouts).
 *
 * A decoder is only built for the layouts whose stream carries a coder of its own. A constant,
 * a dictionary, a fixed alphabet, a counter and the value-domain layouts all put their bytes
 * in the block raw and are walked per line by decompressIdField, so their slot holds no
 * decoder - which is what the nulls in idDecoders mean.
 */
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
        } else if (coderName == "coder_affix_match" || coderName == "coder_bwt_cm") {
            /*
             * Which decoder the sub-stream needs is the stream's own business (see
             * CoderFactory::makeFieldDecoder), and the level it was written with
             * travels in its own meta.
             */
            FieldDecoderArgs idArgs;
            idArgs.level = idStreamMeta[i]["coder"]["level"].asInt();
            std::shared_ptr<coder> idDec =
                CoderFactory::makeFieldDecoder(coderName, io.get(), idArgs);
            if (idDec == nullptr) {
                LOG_ERROR("Unsupport coder name:%s", coderName.c_str());
                return -1;
            }
            idDecoders.push_back(idDec);
        } else if (coderName == "coder_qname") {
            idUsesQnameCoder = true;
            idDecoders.push_back(CoderFactory::makeFieldDecoder(coderName, io.get(),
                                                                FieldDecoderArgs()));
        } else {
            LOG_ERROR("Unsupport coder name:%s", coderName.c_str());
            return -1;
        }

        /*
         * Numeric sub-streams carry no per-line terminator, so they cannot
         * be decoded one segment per line the way the textual path does.
         * Decode the whole varint stream once here and serve one value per
         * line from idNumericPos (see decompressIdField). The value-domain
         * layout ("intd") works the same way, except that its payload needs
         * no undoing: the range decoder walks it one value per line.
         */
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
            /*
             * A counter's deltas, run-length coded (see id_int::Counter), expanded into the
             * very zigzag varint form the numeric layout carries - the layout this one
             * competes with - so nothing downstream has to know which of the two won: one
             * delta per line, same buffer, same per-line accumulation.
             */
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
            /*
             * A low-cardinality segment: the block's few distinct texts sit at the head of the
             * payload and every line is an index into them, coded through id_int::Model - the
             * same value layout the "intd" case uses, so the per-sub-stream state already
             * kept here serves both. What the index means differs, and the reconstruction
             * below knows that.
             */
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
            /*
             * The index stream is a counter layout (see id_int::Counter): its deltas are
             * expanded here into one absolute index per line, which is what the per-line
             * reconstruction reads.
             */
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
                wp += tlenPutVarint(nbuf + wp, (uint32_t)acc);
            }
            nlen = wp;
        } else if (isHex) {
            /*
             * Fixed-alphabet characters (see the encoder): the payload holds one uniform
             * code per character over the block's alphabet, and a length per line when the
             * lengths vary. Nothing is decoded here; the reconstruction walks it per line.
             */
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

/*
 * Build the SEQ column: with a reference, the match / whole-match sub-stream, the exception
 * streams and (when base lengths vary) the base-length stream; without one, the single stream
 * for the whole column, which initSeqWholeBlockDecoder either prepares or decodes outright.
 * readOffset walks the block in both cases.
 */
int32_t SamCodecActuator::initSeqFieldDecoders(Json::Value& baseMeta, uint32_t idx,
                                              RoughIOBlock* outputBlock)
{
    maxBaseLength = baseMeta["maxlen"].asUInt();
    minBaseLength = baseMeta["minlen"].asUInt();
    LOG_DEBUG("maxBaseLen = %d, minBaseLen = %d", maxBaseLength, minBaseLength);
    bool isUseReference = pRefeGene != nullptr && baseMeta.isMember("streams");
    if (!isUseReference) {
        /*
         * One stream for the whole column: initSeqWholeBlockDecoder builds its decoder
         * (or decodes it outright, for the whole-block coder) and answers what it
         * consumed.
         */
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
        /* The sub-stream's meta travels with the view: a whole-block coder
           needs it to decode (coder_fc reads bi/bn/lv and its trailing index
           offset from there, while coder_bwt_cm carries the block size in
           the stream itself and ignores it). */
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
        /*
         * SEQ exception streams (see SeqExceptionClass): one per character present,
         * each a strictly increasing list of block offsets naming its character in
         * "ch" - read by predecodeSeqExceptionStreams, which is also where the older
         * absolute "npos" layout is handled.
         */
        if (predecodeSeqExceptionStreams(baseMeta, baseMetaStreams, id) != 0) {
            return -1;
        }

        if (minBaseLength != maxBaseLength) {
            id++;
            if (predecodeSeqBaseLengthStream(baseMetaStreams, id) != 0) {
                return -1;
            }
        }
    }
    return 0;
}

/*
 * Record the stream offset of a field whose column is decoded lazily, and answer whether this
 * was such a field (in which case the caller moves on to the next one).
 *
 * OPTION tag-split builds no decoder here and does not advance readOffset: the start of the
 * field's stream is recorded, and decompressOptionField decodes all of its columns from the
 * streams array on the first OPTION line. readOffset cannot be used for that at the time -
 * line-by-line decoders such as decompressIdField advance it, so by the time the OPTION line
 * is reached it has already shifted (measured as off by the ID field's dstlen).
 *
 * PNEXT qname-rebuild keeps only an exception stream in streams[7]["streams"][0]: no
 * line-by-line decoder either, just its offset, with readOffset advanced past it so the later
 * fields stay aligned.
 */
bool SamCodecActuator::recordDeferredFieldOffset(uint32_t idx, Json::Value& fieldMeta)
{
    if (idx == 11 && fieldMeta.isMember("mode") && fieldMeta["mode"].asString() == "tag_split") {
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

    /*
     * When FLAG/POS/PNEXT are compressed in textual form (affix), the tracking
     * maps must be repopulated here: the later QUAL and reference-sequence
     * stages all read mappedPos/mappedChr/mappedFlag. This must stay consistent
     * with what decompressNumber fills for the binary form.
     */
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

/*
 * Decompression side of the OPTION field: before decoding the whole block, all
 * id-sequence columns and tag-value columns are decoded and cached; when
 * emitting lines, `NAME:TYPE:VALUE` is reassembled from that line's id
 * sequence, fully symmetric with the compression side.
 *
 * On compression the entire OPTION field is one field (field 11), and on
 * decompression all lines are decoded at once with output advancing line by
 * line from the block start. Everything is lazily decoded on the first call
 * (lineNo == 0 with an empty cache).
 */
int32_t SamCodecActuator::decompressOptionField(uint32_t lineNo, uint8_t splitFlag,
                                                RoughIOBlock* outputBlock,
                                                const Json::Value& fieldMeta) {
    if (!fieldMeta.isMember("tags") || fieldMeta["mode"].asString() != "tag_split") {
        /*
         * affix form: the OPTION of all lines in the block is one column. When a
         * line has no OPTION, the decoded result is empty and only a delimiter
         * was appended (fieldLen == 1); the '\t' just appended after QUAL must
         * be removed (same handling as empty lines in tag_split). The per-block
         * field count is the block-wide maximum, so lines without OPTION also
         * reach this path.
         */
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

    /* Lazy: the entire OPTION column is decoded at once. */
    if (optionCacheEmpty) {
        if (0 != decodeOptionColumn(fieldMeta)) {
            return -1;
        }
    }
    if (lineNo < optionRecLines.size()) {
        const std::string& content = optionRecLines[lineNo];
        if (!content.empty()) {
            memcpy(outputBlock->getCurrent(), content.data(), content.size());
            outputBlock->setDataLen(outputBlock->getDataLen() + (uint32_t)content.size());
            *(outputBlock->getCurrent()) = splitFlag;
            outputBlock->setDataLen(outputBlock->getDataLen() + 1);
            return (int32_t)content.size() + 1;
        }
    }
    /* This line has no OPTION: turn the '\t' just appended after QUAL back into '\n' (consistent with the old logic). */
    uint8_t* pEnd = outputBlock->getCurrent();
    if (pEnd > outputBlock->getBuffer()) {
        *(pEnd - 1) = '\n';
    }
    return 0;
}

int32_t SamCodecActuator::decodeOptionColumn(const Json::Value& fieldMeta) {
    optionRecLines.clear();
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

    /* Decode the id-sequence column (line-wise text, comma-separated). */
    std::vector<std::vector<uint8_t>> recIds;
    uint32_t idStreamDst = 0;
    for (uint32_t i = 0; i < streams.size(); ++i) {
        if (streams[i]["sname"].asString() == "ids") {
            idStreamDst = streams[i]["dstlen"].asUInt();
            break;
        }
    }
    {
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + optBase, idStreamDst, "OPTION ids");
        ioVector.push_back(io);
        std::shared_ptr<coder_bwt_cm> c = std::make_shared<coder_bwt_cm>(io.get());
        uint8_t tmp[1 << 16];
        int32_t n;
        std::string line;
        std::vector<std::vector<uint8_t>> lines;
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
            lines.push_back(ids);
        }
        recIds.swap(lines);
        optBase += idStreamDst;
    }

    /* Decode each tag-value column (line-wise text). */
    std::vector<std::vector<std::string>> tagVals(nTag);
    for (uint32_t i = 0; i < streams.size(); ++i) {
        const Json::Value& s = streams[i];
        if (s["sname"].asString() != "tag") continue;
        uint32_t tagIdx = 0;
        for (uint32_t t = 0; t < nTag; ++t) {
            if (tags[t][0].asString() == s["tag"].asString()) { tagIdx = t; break; }
        }
        uint32_t dstlen = s["dstlen"].asUInt();
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + optBase, dstlen, "OPTION tag");
        ioVector.push_back(io);
        std::shared_ptr<coder_bwt_cm> c = std::make_shared<coder_bwt_cm>(io.get());
        uint8_t tmp[1 << 16];
        int32_t n;
        while ((n = c->decode_line(tmp, sizeof(tmp), '\n', false)) > 0) {
            size_t tlen = (size_t)n;
            if (tlen > 0 && tmp[tlen-1] == '\n') tlen--;
            tagVals[tagIdx].emplace_back((char*)tmp, tlen);
        }
        optBase += dstlen;
    }

    /* Reassemble the OPTION text line by line. */
    const uint32_t lines = (uint32_t)recIds.size();
    optionRecLines.resize(lines);
    std::vector<size_t> colPos(nTag, 0);
    for (uint32_t r = 0; r < lines; ++r) {
        std::string out;
        const auto& ids = recIds[r];
        for (size_t k = 0; k < ids.size(); ++k) {
            uint32_t tid = ids[k];
            if (tid >= nTag || colPos[tid] >= tagVals[tid].size()) continue;
            const std::string& v = tagVals[tid][colPos[tid]++];
            if (k) out += '\t';
            out += tags[tid][0].asString();
            out += ':';
            out += tags[tid][1].asString();
            out += ':';
            out += v;
        }
        optionRecLines[r] = out;
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
        /* PNEXT is valid only if:
         * - FLAG bit 0x1 is set (paired-end sequencing)
         * - FLAG bit 0x8 is not set (mate is mapped)
         */
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
    /*
     * Decode an unsigned varint (LEB128), mirroring compressPosFieldDelta.
     * Bytes are read one at a time (fixed-length decode_line) until the
     * continuation bit is clear. The accumulator is cleared whenever the
     * chromosome (RNAME) changes, detected from the per-line chromosome
     * index decoded by decompressChrName — no reset list is stored in the
     * metadata.
     */
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

/*
 * Pre-decode POS/CIGAR/PNEXT and cache the results per line. TLEN
 * reconstruction needs the full mate index and reference spans; while the main
 * loop decodes line by line, a mate may lie in the second half of the block and
 * be unavailable. Here all lines are decoded once up front; the main loop then
 * copies the cache directly without decoding again.
 */
/*
 * The SEQ column with no reference: it is one stream, and this builds its decoder (or decodes
 * it outright, for the whole-block coder). Answers the bytes it consumed, or -1 on an error
 * already reported.
 */
int32_t SamCodecActuator::initSeqWholeBlockDecoder(const Json::Value& baseMeta, uint32_t idx,
                                                   RoughIOBlock* outputBlock)
{
    const std::string coderName = baseMeta["coder"]["magic"].asString();
    const uint32_t dstLength = baseMeta["totaldstlen"].asUInt();
    const uint32_t srcLength = baseMeta["totalsrclen"].asUInt();
    LOG_DEBUG("srclen = %d, dstlen = %d", srcLength, dstLength);

    if (coderName == "coder_fc") {
        /*
         * The whole-block SEQ is decoded into the tail of the outputBlock buffer as staging,
         * then decompressBase moves it to the head line by line. No separate buffer is
         * allocated: that would add one malloc/free of the whole SEQ per block plus page
         * faults on first touch, and raise peak memory by threads x SEQ size, while the number
         * of copies stays the same.
         *
         * Invariant: never realloc outputBlock within a block, or this tail pointer would
         * dangle. Capacity is guaranteed by the one-shot block_size*2 pre-allocation at the
         * block entry (see decompress())-head output <= block_size and tail staging <=
         * block_size, exactly 2x.
         */
        coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstLength, &ioErrSink, "SEQ");
        baseIo.meta = baseMeta;
        baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
        coder_fc baseDecoder = coder_fc(&baseIo);
        if (baseDecoder.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcLength,
                                    srcLength, UINT8_MAX, false) < 0) {
            LOG_ERROR("Decode SEQ by coder_fc failed, srclen = %u", srcLength);
            return -1;
        }
    } else if (coderName == "coder_bwt_cm") {
        std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "SEQ");
        ioVector.push_back(io);
        fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
    } else {
        LOG_ERROR("Unsupported coder name:%s", coderName.c_str());
        return -1;
    }
    return (int32_t)dstLength;
}

/*
 * The reference layout's match stream. Under the RLE split (see compressBaseWithRef) the two
 * sub-streams
 *   "m"    -> varint run lengths of the zero runs
 *   "mval" -> the surviving non-zero values (1 byte each)
 * are decoded here and expanded back into the original sparse byte stream (orgrawlen bytes),
 * which decompressBase then serves per record through matchBlockDecode. Without the split the
 * stream is simply kept as the record's decoder.
 *
 * `streamId` points at "m" on entry and at the last sub-stream the layout used on return, so
 * the caller can carry on from there; `extraOffset` reports the bytes it must skip past (the
 * "mval" stream, when there is one).
 */
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

    /*
     * The value sub-stream (if any) sits immediately after the run sub-stream. A fully
     * reference-matching block has no surviving non-zero values, so compressBaseWithRef writes
     * no "mval" sub-stream at all; the caller then carries on from the run stream, so the npos
     * / baselen sub-streams that follow are addressed as if the value stream never existed.
     */
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
        /* Each sub-stream carries its own coder name: the two halves of the RLE layout are
           written independently, so the value stream must not be decoded with the run
           stream's coder. */
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

    /* Expand: each run is followed by one surviving value. calloc zero-fills, so the expanded
       zeros need no explicit write. */
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

/*
 * The match stream of the reference layout when it was written whole - the legacy layout with
 * no "rle" member. Mirrors predecodeSeqMatchStream, which is the split case.
 */
int32_t SamCodecActuator::predecodeSeqWholeMatchStream(const Json::Value& streams, uint32_t streamId,
                                                       coder_io* matchIo, uint32_t idx)
{
    const std::string matchCoderName = streams[streamId]["coder"]["magic"].asString();
    /*
     * Which decoder the sub-stream needs is the stream's own business (see
     * CoderFactory::makeFieldDecoder), and whether it hands back a whole block or one record at
     * a time is the coder's own trait rather than a name this side knows (see
     * decoderIsWholeBlock): a whole-block coder's stream is decoded here and decompressBase
     * slices it per record, while every other coder - including the legacy layout with no "rle"
     * member - is read line by line.
     */
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

/*
 * The baselen sub-stream: the reads whose length is not already implied by their CIGAR, which
 * the SEQ column needs before it can slice the payload per record (see decompressBase). The
 * delta layout holds the entry ordinals as one run of forward deltas followed by the lengths;
 * an archive written before it holds absolute 4-byte (ordinal, length) pairs instead.
 */
int32_t SamCodecActuator::predecodeSeqBaseLengthStream(const Json::Value& streams, uint32_t streamId)
{
    if (streams[streamId]["sname"].asString() != "baselen") {
        LOG_ERROR("check sub stream failed. sname not match");
        return -1;
    }
    const uint32_t dstlen = streams[streamId]["dstlen"].asUInt();
    const uint32_t srclen = streams[streamId]["srclen"].asUInt();
    uint8_t* baseLenBuffer = MemoryUtil::safeAlloc<uint8_t>(srclen);

    if (streams[streamId]["coder"]["magic"].asString() == "coder_bwt_cm") {
        coder_io baseLenIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ length");
        auto baseLenCoder = std::make_unique<coder_bwt_cm>(&baseLenIo);
        if (baseLenCoder->decode_line((uint8_t*)baseLenBuffer, srclen, UINT8_MAX, false) < 0) {
            MemoryUtil::safeFree(baseLenBuffer);
            LOG_ERROR("Decode base lengths failed");
            return -1;
        }
    } else {
        MemoryUtil::safeFree(baseLenBuffer);
        LOG_ERROR("check sub stream failed. coder not match. coder = %s.",
            streams[streamId]["coder"]["magic"].asString().c_str());
        return -1;
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

/*
 * The SEQ column's exception streams (see SeqExceptionClass): one per character that occurred,
 * each a strictly increasing list of block offsets with its character in the stream meta's
 * "ch". A character that never occurs has no stream, so this reads sub-streams until one turns
 * up that is not an exception stream (see seqExcFormOf / seqExcLayoutOf); the per-record refill
 * then walks the lists collected here.
 *
 * The one older archive still read is the first layout's lone "npos": a single list of absolute
 * offsets holding every position that layout knew of, 'N' and 'n' alike, with no "ch" and with
 * its count in the field-level "ncount". It wrote an 'N' back at each of those positions, and
 * so does this.
 */
int32_t SamCodecActuator::predecodeSeqExceptionStreams(const Json::Value& baseMeta,
                                                       const Json::Value& streams,
                                                       uint32_t& streamId)
{
    seqExc.clear();
    for (;;) {
        const std::string name = (streams.size() > (size_t)streamId + 1)
                                     ? streams[streamId + 1]["sname"].asString()
                                     : std::string();
        /* Every form name identifies an exception stream, so the block's exception streams are
           the ones that come next and this stops at the first other name. */
        if (seqExcFormOf(name) == SeqExcForm::None) {
            break;   /* this block has no further exception stream */
        }
        streamId++;
        const SeqExcLayout layout = seqExcLayoutOf(baseMeta, streams[streamId]);
        const uint8_t ch = layout.ch;
        const uint32_t dstlen = streams[streamId]["dstlen"].asUInt();
        const uint32_t srclen = streams[streamId]["srclen"].asUInt();
        if (streams[streamId]["coder"]["magic"].asString() != "coder_bwt_cm") {
            LOG_ERROR("check sub stream failed. coder not match. coder = %s.",
                streams[streamId]["coder"]["magic"].asString().c_str());
            return -1;
        }

        coder_io excIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ npos");
        auto excCoder = std::make_unique<coder_bwt_cm>(&excIo);
        /* The payload carries the positions alone - the character travels in the stream's meta
           ("ch") - and every form is expanded by the same reader. */
        std::vector<uint8_t> raw(srclen);
        std::vector<uint32_t> pos;
        if (excCoder->decode_line(raw.data(), srclen, UINT8_MAX, false) < 0 ||
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

/*
 * The QUAL column. Its frequency-table sub-stream is decoded first, because the alphabet it
 * carries is what the value coder is built with; the value stream's decoder is then built by
 * the magic its meta carries (see CoderFactory::makeQualDecoder, the only place that knows
 * which coder takes which parameters).
 */
int32_t SamCodecActuator::initQualFieldDecoders(const Json::Value& qualMeta)
{
    const Json::Value& qualStreamMeta = qualMeta["streams"];
    if (qualStreamMeta.size() != 2) {
        LOG_ERROR("quality streams check failded, size = %d", qualStreamMeta.size());
        return -1;
    }
    if (qualStreamMeta[1]["coder"]["magic"] != "coder_bwt_cm") {
        LOG_ERROR("Unsupport coder type: %s", qualStreamMeta[1]["coder"]["magic"].asString().c_str());
        return -1;
    }

    const uint32_t qualDstLength = qualStreamMeta[0]["dstlen"].asUInt();
    const uint32_t freqDstLength = qualStreamMeta[1]["dstlen"].asUInt();
    coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualDstLength, freqDstLength, &ioErrSink, "QUAL freq table");
    auto qualFreqCoder = std::make_unique<coder_bwt_cm>(&qualFreqIo);
    const uint32_t qualFreqSrcLength = qualStreamMeta[1]["srclen"].asUInt();
    /* Same as fastq_actuator: counting with uint8_t wraps when the alphabet exceeds 127
       symbols, causing a heap out-of-bounds write. */
    const uint32_t qualFreqArrLength = qualFreqSrcLength / sizeof(uint16_t);
    uint16_t* qualFreqArr = new uint16_t[qualFreqArrLength];
    const int32_t qualFreq = qualFreqCoder->decode_line((uint8_t*)qualFreqArr, qualFreqSrcLength, UINT8_MAX, false);
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
    /*
     * Turning the magic the stream carries into a decoder is the factory's business (see
     * CoderFactory::makeQualDecoder); what this side supplies is the alphabet it has just
     * decoded and the prior the engine holds. The stream's meta says whether the encoder
     * started from one, and the package index is the key the prior is stored under: under
     * piped input the absolute offset degenerates to 0, so looking it up by offset always
     * misses.
     */
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

/*
 * Any ordinary column: one stream, whose decoder is the one its magic names. The stream's
 * compression level and - for coder_arith - the file-level POS delta prior travel with the
 * meta; a missing prior must fall back to the coder's uniform table, exactly as the encoder
 * cold-starts (see the encoder's note on that fallback).
 *
 * The position of the stream is recorded as well: preDecodeForTLEN rebuilds decoders for
 * POS(3)/CIGAR(5)/PNEXT(7) from these, which is what lets TLEN reconstruction obtain the full
 * mate index before the column is decoded line by line.
 */
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

    /*
     * The level is handed over for the coders whose stream records one: the textual coder
     * always reads the field meta's value, while coder_arith takes it only when the meta
     * carries one; a bwt_cm stream keeps its own settings instead.
     */
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

/*
 * The TLEN column's exception stream: the records whose TLEN the pair walk cannot produce.
 * Reading it fills tlenCache before the column is decoded, and its layout is read back from the
 * field meta (see the marker note on tlenDecodeVarints): the marker is authoritative where it
 * is there, and without one the archive is older - the fixed layout, unless the stream's size
 * rules that out, in which case the archive was written between the change to varint and the
 * marker and holds varints.
 */
int32_t SamCodecActuator::decodeTlenExceptionStream(const Json::Value& stream,
                                                    const Json::Value& tlenMeta)
{
    const uint32_t srclen = stream["srclen"].asUInt();
    const uint32_t dstlen = stream["dstlen"].asUInt();
    if (stream["coder"]["magic"].asString() != "coder_bwt_cm") {
        LOG_ERROR("Unsupported TLEN exception coder type: %s", stream["coder"]["magic"].asString().c_str());
        return -1;
    }
    uint8_t* excBuffer = MemoryUtil::safeAlloc<uint8_t>(srclen);
    if (excBuffer == nullptr) {
        return -1;
    }
    coder_io tlenIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "TLEN exceptions");
    if (stream["coder"].isMember("level")) {
        tlenIo.meta["level"] = stream["coder"]["level"].asInt();
    }
    coder_bwt_cm tlenDecoder(&tlenIo);
    if (tlenDecoder.decode_line(excBuffer, srclen, UINT8_MAX, false) < 0) {
        MemoryUtil::safeFree(excBuffer);
        LOG_ERROR("Decode TLEN exceptions failed");
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

/*
 * The TLEN column. Under the reconstruction optimization it carries nothing but the exception
 * stream above; otherwise it is an ordinary column whose decoder is built here.
 */
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
    if (coderName != "coder_bwt_cm") {
        LOG_ERROR("Unsupported TLEN coder type: %s", coderName.c_str());
        return -1;
    }
    std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "TLEN");
    ioVector.push_back(io);
    fieldDecoders[fieldIdx] = std::make_shared<coder_bwt_cm>(io.get());
    readOffset += dstLen;
    return 0;
}

/*
 * The QNAME column's own decoders, rebuilt from the recorded sub-stream offsets so the main
 * loop's idDecoders are not consumed by this pre-decode (mirrors decompressIdField).
 *
 * Two kinds of sub-stream are read up front: the numeric layout, which has no per-line
 * terminator (one varint per line, served from numericPos), and the counter layout, which
 * carries the same one-delta-per-line content and rebuilds it from its own payload without
 * touching state the main loop decodes through. The value-domain and constant layouts stay
 * with the main loop and get a null decoder, which the rebuild loop reads as "handled
 * elsewhere".
 */
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
        if (idStreamCoders[si] == "coder_affix_match" || idStreamCoders[si] == "coder_bwt_cm") {
            FieldDecoderArgs idPreArgs;
            idPreArgs.level = lvl;
            state.decoders.push_back(
                CoderFactory::makeFieldDecoder(idStreamCoders[si], io.get(), idPreArgs));
        } else if (idStreamCoders[si] == "coder_qname") {
            state.decoders.push_back(CoderFactory::makeFieldDecoder(idStreamCoders[si], io.get(),
                                                                    FieldDecoderArgs()));
        } else {
            state.decoders.push_back(nullptr);
        }
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

/*
 * Pre-decode the FLAG column (field 1).
 *
 * PNEXT is stored as a delta against POS only when the mate position is meaningful (FLAG bit
 * 0x1 set and bit 0x8 clear) and holds the original value otherwise, so this column has to be
 * known before PNEXT is rebuilt; it mirrors compressPNextFieldDelta.
 *
 * The column is read through a decoder of its own, so the one the main loop keeps is not
 * consumed. Which coder that decoder is comes from the stream's magic (see
 * CoderFactory::makeFieldDecoder); a null result means this build cannot read that coder, and
 * the flag column is simply left unreconstructed.
 */
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

/*
 * Pre-decode RNAME (field 2) so the POS delta chain below can detect chromosome switches from
 * mappedChr, mirroring the encoder. Another decoder of its own, so the main loop's RNAME
 * decoder state stays untouched.
 */
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

/*
 * Build the pre-decoder of one of the columns TLEN reconstruction needs. The decoder is the
 * one the stream's magic names (see CoderFactory::makeFieldDecoder), with the level the field
 * meta records and - for coder_arith - the file-level POS prior the encoder used, or the
 * arithmetic decode diverges on the first symbol (the factory applies it).
 *
 * Answers false when there is nothing to pre-decode here: the field has no stream of its own,
 * or this build cannot read the coder that wrote it.
 */
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

/*
 * Rebuild one column from its pre-decoder, line by line, and leave it in the cache the TLEN
 * walk reads (tlenPreDecodedFields) together with the maps the reconstruction needs: POS
 * fills mappedPos, CIGAR fills baseLengthBuffer and the parsed op lists, PNEXT fills
 * nextMappedPos.
 *
 * Three field forms travel here. CIGAR is always textual. POS and PNEXT are fixed-width binary
 * by default and textual when their mode says so; POS in "pos_delta" mode is a zigzag varint
 * stream read byte by byte until the continuation bit clears (mirroring
 * compressPosFieldDelta), whose values are the deltas from the previous line - reset to 0 at a
 * chromosome switch, which is what the RNAME pre-decode above is for. PNEXT in "pnext_delta"
 * mode holds a delta against POS only for records whose mate is mapped (FLAG 0x1 set, 0x8
 * clear) and the absolute value otherwise.
 */
int32_t SamCodecActuator::rebuildTlenColumn(uint32_t fieldIdx, const TlenColumnPredecoder& pre)
{
    const std::string& coderName = pre.coderName;
    const std::string& mode = pre.mode;
    std::vector<std::string>& fieldCache = tlenPreDecodedFields[fieldIdx];
    fieldCache.clear();
    fieldCache.reserve(samLine);

    uint8_t tmpBuf[1024];
    /* Delta chain of pos_delta mode: the absolute POS reconstructed from the previous line. */
    int64_t posPrev = 0;
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        /*
         * coder_affix_match points last at the caller's output buffer unless need2hold is set;
         * it must be set here, otherwise the cross-line context would be overwritten by the
         * next decode. bwt_cm buffers internally on its own.
         */
        const bool need2hold = CoderFactory::decoderHoldsCallerBuffer(coderName);
        /*
         * Field form: CIGAR is always textual; POS/PNEXT default to fixed-width binary and are
         * decoded as text only when mode is textual. POS uses a zigzag varint stream: the value
         * is read byte-by-byte (fixed-length decode_line) until the continuation bit is clear.
         * Textual data must never be decoded as fixed-length, otherwise bwt_cm's fixed-length
         * branch would spin past the block boundary.
         */
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
            len = pre.decoder->decode_line(tmpBuf, (uint32_t)sizeof(tmpBuf), '\t', need2hold);
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
                cigarReadLen[lineNo] = parseCigarRefConsumed(tmpBuf, len);
                baseLengthBuffer[lineNo] = parseCigar(tmpBuf, len);
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
                        /* PNEXT is valid only if:
                         * - FLAG bit 0x1 is set (paired-end sequencing)
                         * - FLAG bit 0x8 is not set (mate is mapped)
                         */
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

/*
 * Pre-decode QNAME (field 0) into decodedQnames.
 *
 * PNEXT in pnext_qname_rebuild mode is rebuilt by pairing records that share a QNAME, so every
 * line's QNAME must be known before PNEXT is reconstructed. The column is decoded with the
 * per-sub-stream decoders initDecoder already built (see predecodeIdStreams), mirroring
 * decompressIdField.
 */
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
            // Reconstruct QNAME from split sub-streams. Each sub-stream
            // encodes its segment including the trailing separator (the
            // split symbol is part of the segment bytes), mirroring
            // decompressIdField; so we simply append each decoded segment.
            for (uint32_t si = 0; si < preIdDec.size(); ++si) {
                if (!preIdDec[si] && numBufs[si] == nullptr) continue;
                if (numBufs[si] != nullptr) {
                    /* Numeric sub-stream: one varint per line, no separator byte
                       (the textual path strips a trailing '\t' here anyway). */
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
                    /* Value-domain segment: the payload holds the values themselves, one per
                       line, and carries no separator (the split symbol belongs to the segment
                       that ends with it, exactly as in the numeric layout). */
                    char ntmp[24];
                    const uint64_t v = idIntModes[si].decode(idIntCoders[si]);
                    const int nn = snprintf(ntmp, sizeof(ntmp), "%llu", (unsigned long long)v);
                    if (nn > 0) {
                        qname.append(ntmp, (size_t)nn);
                    }
                    continue;
                }
                if (!preIdDec[si]) {
                    /* Every layout without a coder of its own is handled above; getting here
                       means a mode this build does not know, which must not be dereferenced. */
                    LOG_ERROR("Predecode id sub-stream(%u) has no decoder", si);
                    return -1;
                }
                uint8_t sep = (si < idAnalysis.symbols.size()) ? idAnalysis.symbols[si] : UINT8_MAX;
                // coder_affix_match keeps cross-line context in an internal
                // buffer only when need2hold is set; without it, it points
                // last at our fixed qnameBuf and the next line would corrupt
                // the prefix reference.
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

/*
 * PNEXT (field 7) in pnext_qname_rebuild mode.
 *
 * The stream carries only the exception (line, delta) pairs in streams[7]["streams"][0]; every
 * other line is rebuilt here by pairing the records that share a QNAME and are mutually mapped.
 * mappedPos is already populated by the POS (field 3) pre-decode, and the flags and QNAMEs come
 * from the two pre-decodes above. The result fills tlenPreDecodedFields[7] and nextMappedPos
 * for every line.
 */
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
            /* One buffer size for both varint layouts, so the parse loop below is shared:
               an entry is at most a u32 varint (<= 5 B) plus a zigzag i64 varint (<= 10 B). */
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
                /* Every ordinal as a forward delta, then every delta as a zigzag varint;
                   where all the deltas are the same value that value is in the meta and the
                   second run is not there at all (see compressPNextFieldDelta). */
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
                /* Absolute contentIdx then a zigzag-LEB128 delta per pair, the layout this
                   replaced; archives written by it still decode here. */
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
                    /* The exception stores the 0-based data-line index (contentIdx on
                       the compression side, i.e. lineIdx - headEndLine). On this
                       decompression side lineNo is also 0-based data-line indexed, so
                       no headEndLine offset is added here. */
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

    /*
     * FLAG (field 1) first: PNEXT is stored as a delta against POS only when the mate position
     * is meaningful, and that is decided by the flag bits. See predecodeFlagColumn, which is
     * also what mirrors compressPNextFieldDelta.
     */
    if (predecodeFlagColumn() != 0) {
        return -1;
    }

    /*
     * Pre-decode RNAME (field 2) so the POS delta chain below can detect
     * chromosome switches from mappedChr, mirroring the encoder. A separate
     * decoder is used so the main loop's RNAME decoder state is untouched.
     */
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

    /*
     * PNEXT (field 7): in pnext_qname_rebuild mode the stream carries only the exceptions and
     * the rest of the column is rebuilt from the QNAME groups (see rebuildPnextByQname).
     */
    std::string pnextMode = streams.isValidIndex(7) && streams[7].isMember("mode")
        ? streams[7]["mode"].asString() : "";
    if (pnextMode == "pnext_qname_rebuild" && rebuildPnextByQname() != 0) {
        return -1;
    }

    buildTlenMateIndex();
    return 0;
}


/*
 * Reconstruct one split segment of this line's ID, whatever layout it was written in, and append
 * it to the output block.
 *
 * The layouts differ in where the bytes come from and whether the split symbol travels with them:
 * "numeric" (and "cnt", which initDecoder expanded into exactly that varint buffer) keeps one
 * value per line in a buffer decoded once, so the value is rendered back to decimal text and the
 * split symbol re-appended here; "dict" carries one index per line into the block's few texts;
 * "hexd" is a uniform code per character over a fixed alphabet; a constant segment is the text
 * itself, stored once; the value-domain layout ("intd"/"intu") decodes one value per line from its
 * model; and the legacy textual layout is one decode_line per line, split symbol included.
 *
 * Returns 0 when the segment is done, -1 on a corrupt or truncated stream.
 */
int32_t SamCodecActuator::reconstructIdSegment(uint32_t splitIdx, Json::Value& splitMeta,
                                                RoughIOBlock* outputBlock, uint32_t& idLength)
{
    const uint32_t splitDstLen = splitMeta["dstlen"].asUInt();
    std::string coderName = splitMeta["coder"]["magic"].asString();

    /* "cnt" is the same layout from here on: initDecoder already expanded the counter's deltas
       into the varint buffer that the numeric layout carries (see id_int::Counter). */
    const std::string idMode = splitMeta.isMember("mode") ? splitMeta["mode"].asString()
                                                         : std::string();
    const bool numericMode = (idMode == "numeric" || idMode == "cnt");
    if (numericMode) {
        /*
         * Numeric layout: the whole varint stream was decoded once by
         * initDecoder (no per-line terminator exists). Serve exactly one value
         * per line here, render it back to decimal text and re-append the split
         * symbol, so the caller sees the same output shape as the textual path.
         */
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
        /* Dictionary segment: one index per line, and the text those few index values stand for
           came with the payload (the trailing split symbol included, like the textual layout).
           initDecoder already turned the index stream into one absolute index per line. */
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
        /* Fixed-alphabet segment: the length when they vary, then one uniform character code
           each, then the split symbol this segment ends with (the encoder stored neither). */
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
        /* Constant segment: stored once for the block, so nothing per line was coded for it.
           The text carries its own trailing split symbol, exactly as the textual path emits it. */
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
        /* Value-domain segment: one value per line, then the split symbol this segment ends
           with (the encoder stored neither as text). */
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
    /*
     * Record the chromosome index for every line (including the special
     * "*" = 0xFFFF and "=" = 0xFFFE values), mirroring compressChrName.
     * The POS delta decoder uses consecutive lines' indices to detect
     * chromosome switches instead of a reset-index list in the metadata.
     */
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
    /*
     * ensureCapacity **must not** be called here: the caller
     * decompressSamByFields captured basePtr = outputBlock->getCurrent() in the
     * SEQ phase for QUAL to use, and any realloc within the block would leave it
     * dangling and produce random garbage. Sufficient buffer space is guaranteed
     * by the one-shot pre-allocation at the block entry (the deterministic
     * block_size upper bound x 2, see decompress()); no further growth happens
     * line by line.
     */
    bool isUserReference = pRefeGene != nullptr && fieldMeta.isMember("streams");
    uint32_t actualBaseLen = 0;
    if (!isUserReference) {
        if (minBaseLength == maxBaseLength) {
            actualBaseLen = baseLengthBuffer[lineNo] == 0 ? maxBaseLength : baseLengthBuffer[lineNo];
            if (fieldMeta["coder"]["magic"].asString() == "coder_fc") {
                memcpy(outputBlock->getCurrent(), pBaseOut, actualBaseLen);
                pBaseOut += actualBaseLen;
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
            } else if (fieldMeta["coder"]["magic"].asString() == "coder_bwt_cm") {
                int32_t decLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), actualBaseLen, UINT8_MAX, false);
                if (decLen < 0 || (uint32_t)decLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %lld, line %d, expect len %d, actual len %d", (long long)inBlockPtr->getBlockId(), lineNo, actualBaseLen, decLen);
                    return -1;
                }
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
            } else {
                LOG_ERROR("Not supported coder name:%s",fieldMeta["coder"]["magic"].asString().c_str());
                return -1;
            }

            *(outputBlock->getCurrent()) = '\t';
            outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        } else {
            if (CoderFactory::decoderIsWholeBlock(fieldMeta["coder"]["magic"].asString())) {
                /* A whole-block coder: its stream is already decoded into the staging buffer
                   (see initDecoder), so this record is copied out of it up to its tab. */
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
            } else if (fieldMeta["coder"]["magic"].asString() == "coder_bwt_cm") {
                /* decode_line's split mode reports BUF_SMALL as soon as it has
                   filled out_len characters without yet seeing the split char.
                   A row whose SEQ is exactly maxBaseLength long still needs one
                   more character for its trailing '\t', so the per-row limit
                   must be maxBaseLength + 1 (the row text is stored with its
                   terminating tab). Without this, a max-length row always
                   fails with CODER_ERR_BUF_SMALL. */
                int32_t decLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), maxBaseLength + 1, '\t', false);
                if (decLen <= 0) {
                    LOG_ERROR("base decode failed in block %lld, line %d: %d", (long long)inBlockPtr->getBlockId(), lineNo, decLen);
                    return -1;
                }
                actualBaseLen = decLen;
                outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
                actualBaseLen -= 1; // Remove \t length
            } else {
                LOG_ERROR("Not supported coder name:%s", fieldMeta["coder"]["magic"].asString().c_str());
                return -1;
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
        /*
         * A payload byte above the 2-bit range is one of the record's own characters, put there
         * because the record does not use the reference (see the encoder). The two encodings do
         * not overlap and the decoder knows which record is which by the same test the encoder
         * applied, so no per-byte tag is needed. Without this member every byte is a 2-bit code,
         * which is what an archive written before it carries.
         */
        const bool litBases = fieldMeta.isMember("litbases");
        /* `src` is the buffer this record's payload was read into: the two entry points below use
           different staging buffers, and both call here. */
        const auto writeDirectBases = [&](const uint8_t* src, int32_t len) {
            if (litBases) {
                memcpy(outputBlock->getCurrent(), src, (size_t)len);
            } else {
                for (int32_t o = 0; o < len; ++o) {
                    outputBlock->getCurrent()[o] = atcg4[src[o]];
                }
            }
        };
        /* The layout is identified by its sub-stream names, not by which coder
           produced them: the run stream's coder follows the field selection
           (BWT_CM or FC), and initDecoder has already expanded both halves into
           buffers before this point. */
        if (baseStream[0]["sname"].asString() == "m") {
            int32_t decoderLen = 0;
            uint16_t mapFlag = mappedFlag.find(lineNo) == mappedFlag.end() ? 4 : mappedFlag[lineNo];
            // Not matched
            if (mapFlag & 0x04) {
                decoderLen = readMatchLine(fieldIdx, baseSquashBuffer, actualBaseLen, lineNo);
                if ((uint32_t)decoderLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                    return -1;
                }
                writeDirectBases(baseSquashBuffer, decoderLen);
            } else {
                /*
                 * CIGAR-segment-based reference rebuild, mirroring the
                 * compression side exactly: M/=/X segments are restored from
                 * the reference (decode_line yields the per-base 2-bit XOR
                 * against the reference), I/S segments are direct 2-bit codes,
                 * and D/N only advance the reference position. Same refPos
                 * progression, so the round-trip is lossless.
                 */
                bool findMappedPos = false;
                int64_t refeMappedPos = 0;
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
                    /*
                     * No op list, or an empty one: there is nothing to walk, so this record took the
                     * encoder's fallback and its bases are not reference-coded - the payload holds
                     * either its characters (this build, see "litbases") or their 2-bit codes (an
                     * older one), which is what writeDirectBases resolves. An empty list is
                     * reachable with a valid mapping: a CIGAR of "*" on a record whose FLAG does not
                     * say unmapped, so the empty case has to be spelled out rather than left to the
                     * walk below, which would read the payload as XOR values.
                     */
                    const bool hasOps = (lineNo < cigarOpList.size()) && !cigarOpList[lineNo].empty();
                    if (!hasOps) {
                        writeDirectBases(baseDiffSquashBuffer, (int32_t)actualBaseLen);
                    } else {
                        const std::vector<CigarOp>& ops = cigarOpList[lineNo];
                        uint32_t readPos = 0;
                        int64_t refPosLocal = refeMappedPos;
                        for (size_t oi = 0; oi < ops.size(); ++oi) {
                            const CigarOp& op = ops[oi];
                            switch (op.op) {
                                case 'M': case '=': case 'X':
                                    if (readPos + op.len > actualBaseLen) { readPos = actualBaseLen; break; }
                                    pRefeGene->getStretch2Bits1Char(refeStrecchBuffer, op.len, refPosLocal);
                                    for (uint32_t i = 0; i < op.len; ++i) {
                                        baseSquashBuffer[i] = refeStrecchBuffer[i] ^ baseDiffSquashBuffer[readPos + i];
                                    }
                                    pRefeGene->getActgFrom2Bits(baseSquashBuffer, op.len, out + readPos);
                                    readPos += op.len;
                                    refPosLocal += op.len;
                                    break;
                                case 'I': case 'S':
                                    if (readPos + op.len > actualBaseLen) { readPos = actualBaseLen; break; }
                                    for (uint32_t i = 0; i < op.len; ++i) {
                                        out[readPos + i] = atcg4[baseDiffSquashBuffer[readPos + i] & 0x3];
                                    }
                                    readPos += op.len;
                                    break;
                                case 'D': case 'N':
                                    refPosLocal += op.len;
                                    break;
                                case 'H': case 'P':
                                default:
                                    break; // consume neither SEQ nor reference
                            }
                        }
                        /* Any residual positions (e.g. inconsistent CIGAR) fall back to direct. */
                        for (; readPos < actualBaseLen; ++readPos) {
                            out[readPos] = atcg4[baseDiffSquashBuffer[readPos] & 0x3];
                        }
                    }
                }
            }

            /*
             * Write the exception characters back (see SeqExceptionClass): one cursor per
             * stream, each a strictly increasing list of block offsets, and a position belongs
             * to at most one stream because one character is not another.
             */
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
    /*
     * One call, whichever coder the stream says wrote this column: the record goes to dst,
     * and the record's already-decoded SEQ (basePtr, corresponding to seqStart on the
     * compression side) is the context. The strand direction travels with the stream itself,
     * and the byte-stream coders need neither: they fetch as many values as this record asks
     * for - the encoding side fed the same records one by one, with no delimiter in between.
     */
    if (qualDecoder->decode_record(dst, actualBaseLen, basePtr, actualBaseLen) < 0) {
        LOG_ERROR("Decode quality failed, len = %u", actualBaseLen);
        return -1;
    }

    /*
     * Reads with missing quality (a single '*' in the original file) were
     * expanded into seqLen '*' on the compression side; here they are folded
     * back into a single '*', ensuring the reconstructed output matches the
     * original.
     */
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
    /* cigarOpList is indexed by the 0-based data line (lineIdx - headEndLine),
       which is always non-negative here, so use an unsigned size_t index to
       avoid signed/unsigned comparison warnings and resize clutter. */
    const size_t contentIdx = static_cast<size_t>(lineIdx - headEndLine);
    if (fieldLen > 1) {
        uint32_t seqLength = parseCigar(outputBlock->getCurrent(), fieldLen);
        baseLengthBuffer[lineIdx] = seqLength;
        cigarReadLen[lineIdx] = parseCigarRefConsumed(outputBlock->getCurrent(), fieldLen);
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

uint32_t SamCodecActuator::parseCigar(uint8_t* cigarString, uint32_t cigarLength) {
    // CIGAR format like 6S30M1I114S, M/I/S/=/X: consume SEQ, D/N/H/P don't consume SEQ, so actual SEQ length is the sum of operations that consume SEQ
    if (cigarString == nullptr || cigarLength == 0) {
        return 0;
    }

    uint32_t seqLength = 0;
    uint32_t currentNumber = 0;
    for (uint32_t i = 0; i < cigarLength; ++i) {
        char ch = cigarString[i];
        if (ch >= '0' && ch <= '9') {
            // Accumulate numbers
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            // When encountering operator, determine if it consumes SEQ
            if (currentNumber > 0) {
                switch (ch) {
                    case 'M':  // Match or mismatch
                    case 'I':  // Insertion to reference sequence
                    case 'S':  // Soft clipping at sequence start
                    case '=':  // Match
                    case 'X':  // Mismatch
                    case 'm':  // Lowercase version
                    case 'i':  // Lowercase version
                    case 's':  // Lowercase version
                    case 'x':  // Lowercase version
                        // These operations consume SEQ length
                        seqLength += currentNumber;
                        break;
                    case 'D':  // Deletion from reference sequence
                    case 'N':  // Skip from reference sequence
                    case 'H':  // Hard clipping at sequence start
                    case 'P':  // Padding (silent deletion)
                    case 'd':  // Lowercase version
                    case 'n':  // Lowercase version
                    case 'h':  // Lowercase version
                    case 'p':  // Lowercase version
                        // These operations don't consume SEQ length
                        break;
                    default:
                        // Unknown operator, ignore
                        break;
                }
                currentNumber = 0;
            }
        }
    }
    return seqLength;
}

void SamCodecActuator::parseCigarOps(uint8_t* cigarString, uint32_t cigarLength, std::vector<CigarOp>& ops) {
    /* The shared parse (see sam_seq_payload.h); the span it answers is parseCigarRefConsumed's.
       Qualified, or the name in this scope would be the member itself. */
    ::parseCigarOps(cigarString, cigarLength, ops);
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