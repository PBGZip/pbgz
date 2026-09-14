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
#include "coder/coder_affix_match.h"
#include "coder/coder_qname.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "utils/md5_util.h"
#include "coder/coder_json.h"

#include <cstdlib>
#include <cstring>
#include "log/logger.h"
#include "sam_info.h"
#include "actg.h"
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
    idPosLength = 0;
    headerSrcLen = 0;
    headerDstLen = 0;
    readOffset = 0;
    baseNCount = 0;
    qualCoder = nullptr;
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

    // Release qualCoder
    qualCoder.reset();
    qualFcv2Decoder.reset();
    qualCmDecoder.reset();

    ioVector.clear();
}

int32_t SamCodecActuator::preAnalysisIdFirstLine(uint8_t* pBuffer, uint32_t bufLen) {
    if (pBuffer == nullptr || bufLen == 0) {
        return -1;
    }

    std::vector<int32_t> currentLinePos;
    for (uint32_t i = 0; i < bufLen; ++i) {
        char ch = pBuffer[i];
        if (idSplitDefault.find(ch) != std::string::npos) {
            idSplitSymbols.push_back(ch);
            currentLinePos.push_back(static_cast<int32_t>(i));
        }
    }

    // Initialize max and min length for each separator
    for (size_t i = 0; i < idSplitSymbols.size(); ++i) {
        idSplitMinLen.push_back(UINT32_MAX);
        idSplitMaxLen.push_back(0);
    }

    // Store first line ID analysis information to idSplitPos
    uint32_t lastPos = 0;
    for (size_t idx = 0; idx < currentLinePos.size(); ++idx) {
        uint32_t pos = currentLinePos[idx];
        uint32_t curLen = pos - lastPos + (0 == idx ? 1 : 0);   // First line no needs offset
        if (curLen < idSplitMinLen[idx]) {
            idSplitMinLen[idx] = curLen;
        }
        if (curLen > idSplitMaxLen[idx]) {
            idSplitMaxLen[idx] = curLen;
        }
        idPosLength++;
        lastPos = pos;
    }

    // Add the current line positions to idSplitPos
    idSplitPos.push_back(currentLinePos);
    return 0;
}

int32_t SamCodecActuator::preAnalysisIdLine(uint8_t* pBuffer, uint32_t bufferLen) {
    std::vector<int32_t> currentLinePos;
    uint32_t lastPos = 0;
    uint32_t lastFindPos = 0;
    for (uint32_t idx = 0; idx < idSplitSymbols.size(); ++idx) {
        uint8_t symbol = idSplitSymbols[idx];
        uint32_t pos = lastPos;
        for (; pos < bufferLen; ++pos) {
            if (*(pBuffer + pos) == symbol) {
                uint32_t curLen = pos - lastFindPos + (0 == idx ? 1 : 0);
                if (curLen < idSplitMinLen[idx]) {
                    idSplitMinLen[idx] = curLen;
                }
                if (curLen > idSplitMaxLen[idx]) {
                    idSplitMaxLen[idx] = curLen;
                }
                currentLinePos.push_back(static_cast<int32_t>(pos));
                idPosLength++;
                lastPos = pos + 1;
                lastFindPos = pos;
                break;
            }
            lastPos = pos + 1;
        }

        if (pos >= bufferLen) {
            idPosLength = UINT32_MAX; // Mark as unavailable
            break;
        }
    }

    // Add the current line positions to idSplitPos
    idSplitPos.push_back(currentLinePos);
    return 0;
}

int32_t SamCodecActuator::preAnalysis() {
    if (inBlockPtr == nullptr) {
        LOG_ERROR("Input block pointer is null");
        return -1;
    }

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
            if (line.length() < 3) {
                return -1;
            }
            headEndLine++;
            if (line.substr(0, 3) == "@HD" || line.substr(0, 3) == "@RG" || line.substr(0, 3) == "@CO") {
                continue;
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
                continue;
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
        } else {
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
    if (idPosLength != UINT32_MAX) {
        for (uint32_t idx = 0; idx < idSplitSymbols.size(); ++idx) {
            if (idSplitMinLen[idx] == 0) {
                idPosLength = UINT32_MAX;
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
                // ID field: compress based on analysis result
                if (idPosLength == UINT32_MAX) {
                    LOG_DEBUG("Id will compress in all.");
                    fieldDstLen = compressIdFieldInAll(fieldSrcLen, fieldMeta);
                } else {
                    /*
                     * QNAME trial-based selection: affix segmentation vs
                     * coder_qname (cross-line deduplication). For concatenated
                     * FASTQs (alternating prefixes, globally random numbers)
                     * qname wins; for a single FASTQ (constant prefix, locally
                     * ordered fragment numbers) affix's shared adjacent prefixes
                     * win. Both coders run a trial on the first QNAME_TRIAL_LINES
                     * lines and are rolled back; the one with the smaller
                     * measured output is used for the full encoding.
                     */
                    const uint32_t QNAME_TRIAL_LINES = 20000;
                    int32_t lenAffix = 0;
                    int32_t lenQname = 0;
                    const int64_t startLen = outBlockPtr->getDataLen();
                    Json::Value metaAffix, metaQname;
                    uint32_t srcAffix = 0, srcQname = 0;
                    /*
                     * The trial (affix split vs coder_qname, each over the first
                     * QNAME_TRIAL_LINES lines) measures a property of the file's naming
                     * scheme, not of one block: a single FASTQ's constant prefix versus a
                     * concatenated file's alternating prefixes does not change from block to
                     * block, and for long reads the "first 20000 lines" usually *is* the whole
                     * block, so the trial was re-encoding every line twice per block
                     * (measured: over half of the QNAME cost). The verdict is therefore
                     * published on the engine once and reused by every later block.
                     */
                    bool useAffix = false;
                    const int qnameDecision = (pbgzEngine != nullptr)
                        ? pbgzEngine->samQnameUseAffix.load(std::memory_order_relaxed) : -1;
                    if (qnameDecision >= 0) {
                        useAffix = (qnameDecision == 1);
                    } else {
                        {
                            PBGZ_PROF_SCOPE(pbgzprof::FIELD_QNAME_TRIAL);
                            int32_t lAffix = compressIdFieldSplit(srcAffix, metaAffix, QNAME_TRIAL_LINES);
                            outBlockPtr->setDataLen(startLen);
                            int32_t lQname = compressIdFieldQname(srcQname, metaQname, QNAME_TRIAL_LINES);
                            outBlockPtr->setDataLen(startLen);
                            lenAffix = lAffix;
                            lenQname = lQname;
                        }
                        useAffix = (lenQname < 0) || (lenAffix >= 0 && lenAffix <= lenQname);
                        if (pbgzEngine != nullptr) {
                            pbgzEngine->samQnameUseAffix.store(useAffix ? 1 : 0, std::memory_order_relaxed);
                        }
                    }
                    if (useAffix) {
                        LOG_DEBUG("QNAME: affix=%d qname=%d -> affix", lenAffix, lenQname);
                        fieldDstLen = compressIdFieldSplit(fieldSrcLen, fieldMeta);
                    } else {
                        LOG_DEBUG("QNAME: affix=%d qname=%d -> qname", lenAffix, lenQname);
                        fieldDstLen = compressIdFieldQname(fieldSrcLen, fieldMeta);
                    }
                }
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

        // LOG_INFO("Compress Rate for block(%d fieldId=%d): src = %d, dst = %d, ratio = %.2f%%.",
        //     outBlockPtr->getBlockId(), fieldIdx, fieldSrcLen, fieldDstLen, ((fieldDstLen * 1.0) * 100) / fieldSrcLen);

        streamMeta.append(fieldMeta);
        totalSrcLen += fieldSrcLen;
        totalDstLen += fieldDstLen;
    }

    // Set SAM metadata
    samMeta["lines"] = lineNum;
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
 * The lowest block allowed to make a per-file choice. Block 0 carries the SAM header
 * and no records, so it never reaches the code that could measure anything - block 1 is
 * the decider in practice, and a file small enough to fit in one block decides on
 * block 0.
 */
const int64_t kFileChoiceDeciderBlockId = 1;

/*
 * Compressed size of one form of the N positions, or 0 when it cannot be coded at all.
 * Only the size is wanted, so the bytes land in a scratch buffer and are discarded; the
 * caller writes the chosen form into the block afterwards.
 */
static uint32_t nposFormSize(const uint8_t* data, uint32_t len, uint8_t level)
{
    if (len == 0) {
        return 0;
    }
    std::vector<uint8_t> scratch((size_t)len * 2 + 65536);
    coder_io trialIo(scratch.data(), (int32_t)scratch.size());
    CoderFactory::applyLevel(&trialIo, CoderType::BWT_CM, level);
    std::shared_ptr<coder> trialCoder = CoderFactory::makeEncoder(CoderType::BWT_CM, &trialIo);
    trialCoder->encode_line(data, len);
    trialCoder->encode_flush();
    return (trialIo.err == coder_io::IO_OK && trialIo.data_len > 0) ? (uint32_t)trialIo.data_len : 0;
}

/*
 * The SEQ characters the 2-bit path cannot carry, one sub-stream per character that is
 * actually present in the block.
 *
 * One stream per character keeps their statistics apart - a run of N and a run of n are
 * different distributions - and it scales to any alphabet: SAMv1 defines SEQ as [A-Za-z=.]+
 * and says "No assumptions can be made on the letter cases", while the 2-bit tables are
 * deliberately case-insensitive (they index on bits 1-2 of the byte, so the case bit is
 * dropped). Lower-case a/c/g/t, IUPAC codes and the like would be folded or mangled by that
 * path, and single-cell and consensus data do contain them.
 *
 * A character that does not occur writes nothing - no stream and no meta - and the character
 * is not in the payload either: the stream's own meta carries it ("ch"), so the payload is
 * only the positions.
 *
 * The positions are written in one of three forms, and the stream name says which one:
 *   "nposd" - forward deltas, varint (every character; the default)
 *   "npos"  - absolute 4-byte block offsets
 *   "nposr" - runs: one gap and one length per run of consecutive positions
 * Only the 'N' class keeps more than one of them as a candidate, and the file-level trial in
 * PreprocessInfo::nposForm picks between them (see the emit loop). A lone "npos" with no "ch"
 * is the first layout's single list, which held 'N' and 'n' together and is still read; see
 * the decoder. The names themselves live in sam_field_layout.h, next to the readers.
 */

/*
 * Encoder-side accumulator for one exception character: a strictly increasing list of its
 * positions in the variable-length delta form above, plus (for the class the trial measures)
 * the absolute form and the runs of consecutive positions.
 */
struct SeqExceptionClass {
    /* Set on the character whose alternative forms the trial compares; positions and runs are
       kept for it and nothing else. */
    bool keepForms = false;
    std::vector<uint32_t> abs;
    std::vector<uint8_t> varint;
    /* Runs of consecutive positions: one gap (from the end of the previous run) and one length
       each. The N positions of ERR14949932 arrive in runs of 34.3 on average, 98% of them
       exactly 35 long, so describing them as runs makes the source ~1% of the size of the
       one-varint-per-position form before the coder sees either. */
    std::vector<uint32_t> runGaps;
    std::vector<uint32_t> runLens;
    uint32_t runStart = 0;
    uint32_t runLen = 0;
    uint32_t prevRunEnd = 0;
    uint32_t count = 0;
    uint32_t last = 0;

    void add(uint32_t pos)
    {
        if (keepForms) {
            abs.push_back(pos);
        }
        uint32_t delta = pos - last;   /* strictly increasing within the block */
        if (keepForms) {
            if (delta == 1 && runLen != 0) {
                runLen++;
            } else {
                closeRun();
                runStart = pos;
                runLen = 1;
            }
        }
        last = pos;
        while (delta >= 0x80) {
            varint.push_back((uint8_t)(delta | 0x80));
            delta >>= 7;
        }
        varint.push_back((uint8_t)delta);
        count++;
    }

    /* Ends the run in progress, if any, so runGaps/runLens describe every position added. */
    void closeRun()
    {
        if (runLen != 0) {
            runGaps.push_back(runStart - prevRunEnd);
            runLens.push_back(runLen);
            prevRunEnd = last + 1;
            runLen = 0;
        }
    }
};

/* The run form's payload: every gap, then every length, each as a varint (see above). */
static void buildRunForm(const SeqExceptionClass& exc, std::vector<uint8_t>& out)
{
    for (uint32_t gap : exc.runGaps) {
        appendVarint(out, gap);
    }
    for (uint32_t len : exc.runLens) {
        appendVarint(out, len);
    }
}

/*
 * Run measure() once for the whole file and cache the answer in slot.
 *
 * The verdict belongs to the lowest block that can measure it, not to whichever block
 * reaches this code first: after block 0 completes, the other coder threads are released
 * together, so "first to get here" is a property of the scheduler while the measurement
 * is a property of the file, and letting the former decide would cost the reproducibility
 * of the output. The blocks released alongside the decider wait on decideCond for its
 * verdict instead of repeating the same trial; a wait that times out means the deciding
 * block had nothing to measure, in which case the first waiter to wake decides.
 */
int32_t SamCodecActuator::decideOncePerFile(std::atomic<int32_t>* slot, int32_t fallback,
                                           const std::function<int32_t()>& measure)
{
    PreprocessInfo* preInfo = preprocessInfoMut();
    if (preInfo == nullptr) {
        /* Nowhere to keep a per-file verdict (decompression, or a caller driving the
           actuator without an engine), so the caller's fixed choice is used. */
        return fallback;
    }

    const int32_t cached = slot->load(std::memory_order_relaxed);
    if (cached >= 0) {
        return cached;
    }

    std::unique_lock<std::mutex> lock(preInfo->decideMutex);
    const int32_t decided = slot->load(std::memory_order_relaxed);
    if (decided >= 0) {
        return decided;   /* decided while this block waited for the lock */
    }

    if (inBlockPtr->getBlockId() > kFileChoiceDeciderBlockId) {
        const bool arrived = preInfo->decideCond.wait_for(lock, std::chrono::seconds(1),
                                                         [slot]() {
                                                             return slot->load(
                                                                        std::memory_order_relaxed) >= 0;
                                                         });
        if (arrived) {
            return slot->load(std::memory_order_relaxed);
        }
        /* The deciding block never measured anything; waiting longer cannot help. */
    }

    const int32_t verdict = measure();
    slot->store(verdict, std::memory_order_relaxed);
    preInfo->decideCond.notify_all();
    return verdict;
}

PreprocessInfo* SamCodecActuator::preprocessInfoMut() const
{
    return (pbgzEngine != nullptr) ? pbgzEngine->getPreprocessInfoMut() : nullptr;
}

/*
 * The coder for the SEQ match stream: BWT_CM and FC are each tried on one block's
 * segment and the smaller wins, after which every later block reuses that verdict.
 *
 * Deciding per block instead, measured on ERR14949932 (1M reads): 10.8% of the
 * compression time at -l5 and 3.3% at -l8 for a byte-identical output, because FC never
 * won a block, 0 of 335 at -l5 and 0 of 34 at -l8.
 */
CoderType SamCodecActuator::seqMatchCoderFor(const uint8_t* payBuf, uint32_t payLen)
{
    if (payLen == 0) {
        /* No such stream in this block, so it is not a basis for the decision. */
        return CoderType::BWT_CM;
    }
    PreprocessInfo* preInfo = preprocessInfoMut();
    if (preInfo == nullptr) {
        return CoderType::BWT_CM;
    }

    const int32_t verdict = decideOncePerFile(&preInfo->seqMatchCoder, (int32_t)CoderType::BWT_CM,
                                             [&]() -> int32_t {
        const CoderType kCandidates[2] = {CoderType::BWT_CM, CoderType::FC};
        CoderType best = CoderType::BWT_CM;
        uint32_t bestLen = UINT32_MAX;
        std::vector<uint8_t> scratch((payLen << 1) + 65536);
        for (int i = 0; i < 2; ++i) {
            const CoderType type = kCandidates[i];
            /* coder_fc refuses inputs of FC_MIN_LEN bytes or less. */
            if (type == CoderType::FC &&
                (payLen <= FC_MIN_LEN || payLen >= (uint32_t)FC_MAX_LEN)) {
                continue;
            }
            coder_io trialIo(scratch.data(), (int32_t)scratch.size());
            CoderFactory::applyLevel(&trialIo, type, engineCompressLevel());
            std::shared_ptr<coder> trialCoder = CoderFactory::makeEncoder(type, &trialIo);
            trialCoder->encode_line(payBuf, payLen);
            trialCoder->encode_flush();
            if (trialIo.err != coder_io::IO_OK || trialIo.data_len <= 0) {
                continue;   /* this coder cannot carry the segment at all */
            }
            if ((uint32_t)trialIo.data_len < bestLen) {
                best = type;
                bestLen = (uint32_t)trialIo.data_len;
            }
        }
        return (int32_t)best;
    });
    return (CoderType)verdict;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;
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
    Json::Value metaSubs;
    Json::Value metaStreams;
    uint32_t totalDstLen = 0;

    if (!pnextExceptions.empty()) {
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
        excStream.reserve(pnextExceptions.size() * 2 + 16);
        uint32_t prevOrdinal = 0;
        bool sameDelta = true;
        const int64_t firstDelta = pnextExceptions[0].second;
        for (const auto& e : pnextExceptions) {
            appendVarint(excStream, e.first - prevOrdinal);
            prevOrdinal = e.first;
            if (e.second != firstDelta) {
                sameDelta = false;
            }
        }
        if (sameDelta) {
            fieldMeta["exc_delta"] = (Json::Value::UInt64)zigzag64(firstDelta);
        } else {
            for (const auto& e : pnextExceptions) {
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

        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;

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
}

int32_t SamCodecActuator::compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                                               uint32_t trialLines) {
    // Similar to FastqActuator::compressIdInSplit, create multiple streams for each split symbol
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    // Process each split symbol (similar to FastqActuator)
    for (uint32_t i = 0; i < idSplitSymbols.size(); ++i) {
        std::vector<size_t>& npos = inBlockPtr->getNpos();
        uint32_t lineNum = npos.size();
        uint8_t* buffer = inBlockPtr->getBuffer();

        std::shared_ptr<coder_io> idIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QNAME sub-stream");
        uint32_t srcLength = 0;

        /*
         * Numeric sub-stream mode: when every occurrence of a segment is a plain
         * decimal integer (typical for FASTQ/SEQ read ordinals such as
         * "SRR2769247.<n>"), the segment is stored as a binary stream of
         * zigzag(delta) LEB128 varints instead of decimal text. Adjacent ordinals
         * are strongly correlated (delta == 1 dominates), and the decimal text
         * wastes that structure across digit boundaries. Non-numeric segments keep
         * the original textual coder_affix_match path unchanged.
         */
        std::vector<std::string> segTexts;      /* full segment, trailing split symbol included */
        std::vector<std::string> segPayloads;   /* segment without the trailing split symbol */
        bool allNumeric = true;
        for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
            if (trialLines != 0 && lineIdx - headEndLine >= trialLines) {
                break;
            }
            uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
            uint8_t* line = buffer + lineStart;

            if (*line == '@') {
                continue;
            }

            uint8_t* segmentStart = nullptr;
            uint32_t segmentLength = 0;
            uint32_t contentId = lineIdx - headEndLine;

            if (i == 0) {
                segmentStart = line;
                segmentLength = idSplitPos[contentId][0] + 1;
            } else {
                uint32_t prevPos = idSplitPos[contentId][i-1];
                uint32_t currPos = idSplitPos[contentId][i];
                segmentStart = line + prevPos + 1;
                segmentLength = currPos - prevPos;
            }

            /*
             * Keep the full segment (trailing split symbol included) for the text
             * path, and the bare payload for the numeric path / the all-digits test.
             */
            const bool hasSep = (segmentLength > 0 && segmentStart[segmentLength - 1] == idSplitSymbols[i]);
            const uint32_t payloadLen = hasSep ? segmentLength - 1 : segmentLength;
            std::string seg((char*)segmentStart, segmentLength);
            std::string payload((char*)segmentStart, payloadLen);
            if (!payload.empty()) {
                if (payload.find_first_not_of("0123456789") != std::string::npos) {
                    allNumeric = false;
                } else if (payload.size() > 1 && payload[0] == '0') {
                    /* A numeric segment with a leading zero (e.g. a zero-padded
                       date/month "01" inside a Nanopore QNAME) cannot round-trip
                       through the zigzag-delta varint stream: the decoder re-emits
                       plain decimal text and would silently drop the "0". Force the
                       lossless textual path for such segments. */
                    allNumeric = false;
                }
            }
            segTexts.push_back(seg);
            segPayloads.push_back(payload);
            srcLength += segmentLength;
        }

        bool numericMode = false;
        uint32_t numericSrcLen = 0;
        if (allNumeric && !segTexts.empty()) {
            /*
             * Worst case 10 bytes per varint (64-bit value, 7 payload bits each).
             * The encoder falls back to the text path if it would not be shorter.
             */
            const uint32_t cap = (uint32_t)(segTexts.size() * 10);
            uint8_t* varBuf = MemoryUtil::safeAlloc<uint8_t>(cap);
            if (varBuf != nullptr) {
                uint32_t pos = 0;
                uint64_t prevVal = 0;
                uint32_t deltaOneCnt = 0;    /* how often the ordinal just increments */
                uint32_t deltaTotal = 0;
                bool ok = true;
                for (size_t k = 0; k < segPayloads.size() && ok; ++k) {
                    if (segPayloads[k].empty()) {
                        continue;               /* keep the row count aligned */
                    }
                    errno = 0;
                    char* endPtr = nullptr;
                    unsigned long long v = std::strtoull(segPayloads[k].c_str(), &endPtr, 10);
                    if (errno == ERANGE || endPtr != segPayloads[k].c_str() + segPayloads[k].size()) {
                        ok = false;
                        break;
                    }
                    const uint64_t cur = (uint64_t)v;
                    const uint64_t delta = (cur >= prevVal) ? (cur - prevVal) : (prevVal - cur);
                    const uint64_t zz = (cur >= prevVal) ? (delta << 1) : ((delta << 1) | 1u);
                    if (zz > 0xFFFFFFFFull) {   /* keep it 32-bit decodable */
                        ok = false;
                        break;
                    }
                    if (pos + 10 > cap) {
                        ok = false;
                        break;
                    }
                    pos += tlenPutVarint(varBuf + pos, (uint32_t)zz);
                    if (cur == prevVal + 1) {
                        deltaOneCnt++;
                    }
                    deltaTotal++;
                    prevVal = cur;
                }
                /*
                 * Varints only pay off when consecutive ordinals are correlated
                 * (successive reads, or mate pairs stored adjacently). When the
                 * file is ordered by alignment position and the ordinals scatter
                 * (e.g. "ERR031968.<random id>"), the deltas are large and nearly
                 * incompressible, and the decimal text wins because affix matching
                 * still finds repeated digit prefixes. Require a meaningful share
                 * of pure +1 steps before switching away from the text path.
                 */
                const bool correlated = (deltaTotal > 0) && (deltaOneCnt * 10 >= deltaTotal);
                if (ok && pos > 0 && correlated) {
                    std::shared_ptr<coder_bwt_cm> numCoder = std::make_shared<coder_bwt_cm>(idIo.get());
                    numCoder->encode_line(varBuf, pos);
                    numCoder->encode_flush();
                    if (idIo->err == coder_io::IO_OK && idIo->data_len > 0) {
                        numericMode = true;
                        numericSrcLen = pos;
                    }
                }
                MemoryUtil::safeFree(varBuf);
            }
        }

        if (!numericMode) {
            std::shared_ptr<coder_affix_match> idCoder = std::make_shared<coder_affix_match>(idIo.get());
            for (size_t k = 0; k < segTexts.size(); ++k) {
                if (segTexts[k].empty()) {
                    continue;
                }
                idCoder->encode_line((uint8_t*)segTexts[k].data(), (uint32_t)segTexts[k].size());
            }
            idCoder->encode_flush();
        }
        if (idIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode id segment overflow: output buffer too small");
            return -1;
        }

        // Update output block data length
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idIo->data_len);

        // Create metadata for this stream (similar to FastqActuator)
        Json::Value tmpMeta;
        tmpMeta["srclen"] = numericMode ? numericSrcLen : srcLength;
        tmpMeta["dstlen"] = idIo->data_len;
        tmpMeta["coder"] = idIo->meta;
        tmpMeta["splitidx"] = i; // Index of split symbol
        /* Absent (or "text") = legacy textual layout; "numeric" = zigzag delta varints. */
        if (numericMode) {
            tmpMeta["mode"] = "numeric";
        }

        streamMeta.append(tmpMeta);
        totalSrcLength += srcLength;
        totalDstLength += idIo->data_len;
    }

    // Set field metadata with streams (similar to FastqActuator)
    fieldMeta["totalsrclen"] = totalSrcLength;
    fieldMeta["totaldstlen"] = totalDstLength;
    fieldMeta["splitsym"] = std::string((char*)idSplitSymbols.data(), idSplitSymbols.size());
    fieldMeta["streams"] = streamMeta;
    fieldMeta["field"] = 0;

    fieldSrcLen = totalSrcLength;

    LOG_DEBUG("SAM ID compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            totalSrcLength, totalDstLength, (double)(totalDstLength * 100)/(double)totalSrcLength);

    return totalDstLength;
}

int32_t SamCodecActuator::compressIdFieldQname(uint32_t& fieldSrcLen, Json::Value& fieldMeta,
                                               uint32_t trialLines) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QNAME");
    std::shared_ptr<coder_qname> qnameCoder = std::make_shared<coder_qname>(fieldIo.get());

    fieldSrcLen = 0;
    for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
        if (trialLines != 0 && lineIdx - headEndLine >= trialLines) {
            break;
        }
        uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
        if (buffer[lineStart] == '@') {
            continue;
        }
        uint32_t contentId = lineIdx - headEndLine;
        uint8_t* idStart = buffer + lineStart;
        /* The QNAME field runs up to the first tab (inclusive, consistent with compressIdFieldInAll). */
        uint32_t idLength = contentPos[contentId][0] + 1;
        qnameCoder->encode_line(idStart, idLength);
        fieldSrcLen += idLength;
    }

    qnameCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode QNAME (coder_qname) overflow: output buffer too small");
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    Json::Value streamMeta;
    Json::Value tmpMeta;
    tmpMeta["srclen"] = fieldSrcLen;
    tmpMeta["dstlen"] = fieldIo->data_len;
    tmpMeta["coder"] = fieldIo->meta;
    tmpMeta["splitidx"] = 0;
    streamMeta.append(tmpMeta);

    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["splitsym"] = "\t";
    fieldMeta["field"] = 0;

    LOG_INFO("SAM QNAME (coder_qname) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
            fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
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
        uint32_t prevTabPos = contentPos[contentId][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentId].size()) ? contentPos[contentId][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos - 1;

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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t  fieldLength = currTabPos - prevTabPos;

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
int32_t SamCodecActuator::compressOptionField(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    /* Parse OPTION line by line: build the tag dictionary, each line's id sequence, and each tag's value column. */
    std::vector<std::pair<std::string, std::string>> tagDict;
    std::map<std::string, int> tagId;
    std::vector<std::vector<uint8_t>> recIds;
    std::vector<std::vector<std::string>> tagVals;
    fieldSrcLen = 0;

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
            recIds.push_back({});
            continue;
        }
        fieldSrcLen += fieldLen;

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
                        auto it = tagId.find(name);
                        if (it == tagId.end()) {
                            id = (int)tagDict.size();
                            tagDict.push_back({name, type});
                            tagId[name] = id;
                            tagVals.push_back({});
                        } else {
                            id = it->second;
                        }
                        ids.push_back((uint8_t)id);
                        tagVals[id].push_back(value);
                    }
                }
                segStart = i + 1;
            }
        }
        recIds.push_back(ids);
    }

    const uint32_t nTag = (uint32_t)tagDict.size();
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
        e.append(tagDict[t].first);
        e.append(tagDict[t].second);
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
        for (size_t r = 0; r < recIds.size(); ++r) {
            for (size_t k = 0; k < recIds[r].size(); ++k) {
                if (k) idStream.push_back(',');
                uint32_t v = recIds[r][k];
                char buf[4];
                int n = snprintf(buf, sizeof(buf), "%u", v);
                for (int b = 0; b < n; ++b) idStream.push_back((uint8_t)buf[b]);
            }
            idStream.push_back('\n');
        }
        std::shared_ptr<coder_io> io = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "OPTION ids");
        std::shared_ptr<coder_bwt_cm> c = std::make_shared<coder_bwt_cm>(io.get());
        size_t pos = 0;
        for (size_t r = 0; r < recIds.size(); ++r) {
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
        const std::vector<std::string>& vals = tagVals[t];
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
        s["tag"] = tagDict[t].first;
        s["type"] = tagDict[t].second;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t  fieldLength = currTabPos - prevTabPos;

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
    /* Counts only CIGAR operations that consume reference sequence: M/D/N/=/X (including lowercase). */
    if (cigarString == nullptr || cigarLength == 0) {
        return 0;
    }

    uint32_t refConsumed = 0;
    uint32_t currentNumber = 0;
    for (uint32_t i = 0; i < cigarLength; ++i) {
        char ch = cigarString[i];
        if (ch >= '0' && ch <= '9') {
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            if (currentNumber > 0) {
                switch (ch) {
                    case 'M': case 'D': case 'N': case '=': case 'X':
                    case 'm': case 'd': case 'n':
                        refConsumed += currentNumber;
                        break;
                    default:
                        break;
                }
                currentNumber = 0;
            }
        }
    }
    return refConsumed;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;

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

int32_t SamCodecActuator::compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    if (pRefeGene == nullptr) {
        LOG_ERROR("Reference genome is not available for base compression with reference");
        return -1;
    }

    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    uint64_t nOffset = 0;
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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* seqStart = line + prevTabPos + 1;
        uint32_t seqLength = currTabPos - prevTabPos - 1;

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
         * Process sequence: record the characters the 2-bit path cannot carry, one list per
         * character (see SeqExceptionClass). A/C/G/T need no entry: that path carries them,
         * and it is case-insensitive, which is exactly why every other spelling - lower
         * case, IUPAC codes - has to be recorded instead of being silently folded.
         */
        uint32_t outLen = 0;
        for (uint32_t n = 0; n < seqLength; n++) {
            const uint8_t ch = seqStart[n];
            if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') {
                SeqExceptionClass& exc = excClasses[ch];
                if (exc.count == 0) {
                    excPresent.push_back(ch);
                    if (ch != (uint8_t)'N') {
                        exc.varint.reserve(64);
                    }
                }
                exc.add(totalBaseLength + n);
            }
            nOffset++;
        }

        totalBaseLength += seqLength;

        /*
         * CIGAR-segment-based reference rebuild.
         *
         * For a mapped read the SEQ is encoded per-base as a 2-bit XOR against
         * the reference, but only on the operations that consume reference
         * (M/=/X). Operations that consume SEQ but not reference (I/S) are
         * stored directly, and reference-only operations (D/N) only advance the
         * reference position without producing any SEQ bytes. This is required
         * so that indels do not cause the read to be compared against a
         * contiguous reference window it does not actually align to.
         *
         * The decompression side mirrors this exact CIGAR walk (same refPos
         * progression and same per-segment encoding), guaranteeing a lossless
         * round-trip.
         */
        bool useReference = (chrId != 0xFFFF && chrId != 0xFFFE && !(flag & 0x04));
        int64_t refPos = 0;
        if (useReference) {
            int64_t chrStartPos = SamInfo::getInstance().getPositionByIndex(chrId);
            if (chrStartPos == -1) {
                useReference = false;
            } else {
                refPos = chrStartPos + startPos - 1; // SAM is 1-based
                /* Reference-consumed span (M/D/N/=/X); guard against reads whose
                   alignment runs past the reference end. */
                auto crlIt = cigarReadLen.find(lineIdx);
                uint32_t refConsumed = (crlIt != cigarReadLen.end()) ? crlIt->second : 0;
                uint64_t needSquash = (uint64_t)(refConsumed >> 2) + !!(refConsumed & 0x3) + 1;
                if (refPos < 0 || ((uint64_t)(refPos >> 2)) + needSquash > (uint64_t)pRefeGene->getSquashLength()) {
                    useReference = false;
                }
            }
        }

        const std::vector<CigarOp>& ops = (contentIdx < cigarOpList.size()) ? cigarOpList[contentIdx] : emptyCigarOps;
        if (useReference && !ops.empty()) {
            // Build the per-base 2-bit stream: XOR with reference for M/=/X,
            // direct 2-bit for I/S.
            uint32_t readPos = 0;
            int64_t refPosLocal = refPos;
            bool consistent = true;
            for (size_t oi = 0; oi < ops.size(); ++oi) {
                const CigarOp& op = ops[oi];
                switch (op.op) {
                    case 'M': case '=': case 'X':
                        if (readPos + op.len > seqLength) { consistent = false; break; }
                        pRefeGene->getStretch2Bits1Char(ref2bitBuf.get(), op.len, refPosLocal);
                        for (uint32_t i = 0; i < op.len; ++i) {
                            uint8_t read2 = (seqStart[readPos + i] >> 1) & 0x3;
                            baseMappedBuffer[readPos + i] = read2 ^ ref2bitBuf[i];
                        }
                        readPos += op.len;
                        refPosLocal += op.len;
                        break;
                    case 'I': case 'S':
                        if (readPos + op.len > seqLength) { consistent = false; break; }
                        for (uint32_t i = 0; i < op.len; ++i) {
                            baseMappedBuffer[readPos + i] = (seqStart[readPos + i] >> 1) & 0x3;
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
                if (!consistent) break;
            }
            if (!consistent || readPos != seqLength) {
                // CIGAR/SEQ length mismatch: fall back to direct encoding.
                actgEncode(seqStart, baseMappedBuffer.get(), seqLength);
                outLen = seqLength;
            } else {
                outLen = seqLength;
                auto crlIt = cigarReadLen.find(lineIdx);
                pRefeGene->updateMatchedGene((uint64_t)refPos,
                    (crlIt != cigarReadLen.end()) ? crlIt->second : seqLength);
            }
        } else {
            // No valid mapping or no CIGAR ops, encode directly
            actgEncode(seqStart, baseMappedBuffer.get(), seqLength);
            outLen = seqLength;
        }
        // Encode the mapped data
        if (outLen > 0) {
            /* Accumulate instead of encoding line-by-line: the coder is picked by
               the total length after the loop (coder_fc is block-only). */
            if (matchLen + seqLength > inBlockPtr->getDataLen()) {
                LOG_ERROR("SEQ match buffer overflow in block %llu: %u + %u > %u",
                          inBlockPtr->getBlockId(), matchLen, seqLength, inBlockPtr->getDataLen());
                return -1;
            }
            std::memcpy(matchBuffer.get() + matchLen, baseMappedBuffer.get(), seqLength);
            matchLen += seqLength;
            srcLen += seqLength;
        }
    }
    /*
     * RLE preprocessing (on by default).
     * The match stream is a sparse 0..3 byte stream in which ~99.28% of the symbols
     * are 0 (reference-matching bases). Split it into two independent sub-streams:
     *   "m"    - varint-encoded run lengths of the zero runs
     *   "mval" - the surviving non-zero values (1 byte each)
     * Encoding them separately lets each stream be modelled on its own instead of
     * forcing one context model to cope with two very different distributions.
     * Backward compatibility: old archives carry no "rle" member and keep taking the
     * original single-stream, line-by-line path, so a new binary still reads them.
     */
    uint32_t rleRunLen = 0;
    uint32_t rleValLen = 0;
    std::unique_ptr<uint8_t[]> runBuffer;
    std::unique_ptr<uint8_t[]> valBuffer;
    const bool useRle = (matchLen > 0);
    if (useRle) {
        uint32_t nNonZero = 0;
        for (uint32_t i = 0; i < matchLen; i++) {
            if (matchBuffer[i] != 0) {
                nNonZero++;
            }
        }
        runBuffer = std::make_unique<uint8_t[]>((nNonZero + 1) * 5 + 16);
        valBuffer = std::make_unique<uint8_t[]>(nNonZero + 16);
        uint32_t rp = 0;
        uint32_t vp = 0;
        uint32_t run = 0;
        for (uint32_t i = 0; i < matchLen; i++) {
            if (matchBuffer[i] == 0) {
                run++;
            } else {
                uint32_t v = run;
                while (v >= 0x80) { runBuffer[rp++] = (uint8_t)(v | 0x80); v >>= 7; }
                runBuffer[rp++] = (uint8_t)v;
                valBuffer[vp++] = matchBuffer[i];
                run = 0;
            }
        }
        uint32_t tail = run;                /* trailing run of zeros */
        while (tail >= 0x80) { runBuffer[rp++] = (uint8_t)(tail | 0x80); tail >>= 7; }
        runBuffer[rp++] = (uint8_t)tail;

        rleRunLen = rp;                     /* length of the run-length segment */
        rleValLen = vp;                     /* length of the value segment */
    }

    /*
     * Sub-stream "m": the run-length segment under RLE, otherwise the whole match
     * stream. BWT_CM and FC both code it and the smaller wins - which of the two is
     * smaller is a property of the data rather than something a rule predicts - and
     * the decoder dispatches on the magic the winner recorded, so the choice needs no
     * other bookkeeping. The trial runs once per file, not per block; see
     * PreprocessInfo::seqMatchCoder.
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
        const uint32_t payLen = useRle ? rleRunLen : matchLen;
        uint8_t* payBuf = useRle ? runBuffer.get() : matchBuffer.get();

        CoderType payType = seqMatchCoderFor(payBuf, payLen);
        /* The verdict came from another block's segment, so a segment too short for
           coder_fc (it refuses FC_MIN_LEN bytes or less) still falls back. */
        if (payType == CoderType::FC &&
            (payLen <= FC_MIN_LEN || payLen >= (uint32_t)FC_MAX_LEN)) {
            payType = CoderType::BWT_CM;
        }
        CoderFactory::applyLevel(matchIo.get(), payType, engineCompressLevel());
        std::shared_ptr<coder> payCoder = CoderFactory::makeEncoder(payType, matchIo.get());
        payCoder->encode_line(payBuf, payLen);
        payCoder->encode_flush();
    }
    if (matchIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode base match stream overflow: output buffer too small");
        return -1;
    }
    // First sub-stream: run lengths (RLE) / whole match stream (legacy layout)
    metaSubs.clear();
    metaSubs["srclen"] = useRle ? (Json::Value::UInt)rleRunLen : (Json::Value::UInt)srcLen;
    metaSubs["dstlen"] = matchIo->data_len;
    metaSubs["coder"] = matchIo->meta;
    metaSubs["sname"] = "m";
    if (useRle) {
        /*
         * orgrawlen = 展开后的原始长度（= 总碱基数），解压端据此分配并还原；
         * rlelen    = 本子流（游程段）的解码输出长度，等于上面写入的 srclen。
         */
        metaSubs["rle"] = (Json::Value::UInt)1;
        metaSubs["orgrawlen"] = (Json::Value::UInt)srcLen;
        metaSubs["rlelen"] = (Json::Value::UInt)rleRunLen;
    }
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + matchIo->data_len);
    totalDstLen += matchIo->data_len;

    /* Second sub-stream: the surviving non-zero values (RLE only). */
    if (useRle && rleValLen > 0) {
        std::shared_ptr<coder_io> valIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match val");
        /* This half stays pinned to BWT_CM: only the run-length stream above follows
           the field selection. Decoding dispatches on the recorded magic either way. */
        {
            CoderFactory::applyLevel(valIo.get(), CoderType::BWT_CM, engineCompressLevel());
            std::shared_ptr<coder_bwt_cm> bwt = std::make_shared<coder_bwt_cm>(valIo.get());
            bwt->encode_line(valBuffer.get(), rleValLen);
            bwt->encode_flush();
        }
        if (valIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ match value stream overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = (Json::Value::UInt)rleValLen;
        metaSubs["dstlen"] = valIo->data_len;
        metaSubs["coder"] = valIo->meta;
        metaSubs["sname"] = "mval";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + valIo->data_len);
        totalDstLen += valIo->data_len;
    }
    totalSrcLen += srcLen;

    /*
     * SEQ exception sub-streams: one per character that is actually present, each carrying
     * the positions of that character (see SeqExceptionClass). A file of plain A/C/G/T writes
     * none of them, and only the 'N' character has two forms to be measured against each
     * other. They are emitted in ascending character order, so the block layout follows the
     * set of characters rather than the order in which they happened to appear.
     */
    std::sort(excPresent.begin(), excPresent.end());
    auto emitExceptions = [&](SeqExceptionClass& exc, uint8_t ch, const char* sname,
                              uint32_t srcLen, const uint8_t* src, uint32_t runCount) -> int32_t {
        std::shared_ptr<coder_io> excIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ npos");
        CoderFactory::applyLevel(excIo.get(), CoderType::BWT_CM, engineCompressLevel());
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(excIo.get());
        subCoder->encode_line(src, srcLen);
        subCoder->encode_flush();
        if (excIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode SEQ exception positions overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = srcLen;
        metaSubs["dstlen"] = excIo->data_len;
        metaSubs["coder"] = excIo->meta;
        metaSubs["sname"] = sname;
        metaSubs["ch"] = (Json::Value::UInt)ch;
        metaSubs["count"] = exc.count;
        if (runCount > 0) {
            /* The run form: how many runs the two varint sections below hold. */
            metaSubs["runs"] = (Json::Value::UInt)runCount;
        }
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + excIo->data_len);
        totalSrcLen += srcLen;
        totalDstLen += excIo->data_len;
        return 0;
    };

    for (size_t ei = 0; ei < excPresent.size(); ++ei) {
        const uint8_t ch = excPresent[ei];
        SeqExceptionClass& exc = excClasses[ch];

        /* Only the 'N' class keeps more than one form; the file-level trial (see
           PreprocessInfo::nposForm) measures the candidates on the deciding block and every
           block is written the way it says. Other characters are always in the delta form. */
        int32_t form = 1;
        if (ch == (uint8_t)'N') {
            PreprocessInfo* preInfo = preprocessInfoMut();
            form = (preInfo != nullptr)
                ? decideOncePerFile(&preInfo->nposForm, 0, [&]() -> int32_t {
                      const uint32_t absSize =
                          nposFormSize((const uint8_t*)exc.abs.data(),
                                       (uint32_t)(exc.abs.size() * sizeof(uint32_t)),
                                       engineCompressLevel());
                      const uint32_t varSize = nposFormSize(exc.varint.data(),
                                                            (uint32_t)exc.varint.size(),
                                                            engineCompressLevel());
                      exc.closeRun();
                      std::vector<uint8_t> runBuf;
                      runBuf.reserve((exc.runGaps.size() << 1) + 16);
                      buildRunForm(exc, runBuf);
                      const uint32_t runSize = nposFormSize(runBuf.data(),
                                                            (uint32_t)runBuf.size(),
                                                            engineCompressLevel());
                      /* 0 = absolute ("npos"), 1 = deltas ("nposd"), 2 = runs ("nposr"). The
                         smallest wins; a tie keeps the earlier form, which is what an archive
                         written before that form existed carries. */
                      int32_t best = 0;
                      uint32_t bestSize = absSize;
                      if (varSize > 0 && (bestSize == 0 || varSize < bestSize)) {
                          best = 1;
                          bestSize = varSize;
                      }
                      if (runSize > 0 && (bestSize == 0 || runSize < bestSize)) {
                          best = 2;
                      }
                      return best;
                  })
                : 0;   /* nowhere to keep a per-file verdict: the absolute form, as before */
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
        if (emitExceptions(exc, ch, sname, srcLen, src, runCount) != 0) {
            return -1;
        }
    }

    if (minBaseLength != maxBaseLength) {
        std::shared_ptr<coder_io> lenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ length");
        /*
         * The reads whose length is not already implied by their CIGAR, as (ordinal, length)
         * pairs. Both columns go out as forward deltas in varint, with the ordinals written as
         * one run before the lengths rather than interleaved with them: a file that is entirely
         * unmapped - or entirely mapped - collapses its whole ordinal column into the same
         * repeated delta, and the lengths keep statistics of their own.
         *
         * The absolute 4-byte pairs this replaces cost 8 B per read: on ERR14949932, whose
         * 3,343,586 reads are all unmapped, 26,748,688 B of source compressed to 6,171,912 B,
         * where CRAM's read-length series needs 1,302,132 B for the same reads.
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
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + lenIo->data_len);
        totalSrcLen += (uint32_t)lenBuf.size();
        totalDstLen += lenIo->data_len;
    }

    // Set metadata
    /* Count of the 'N' stream, which older archives read as their whole N list; every stream
       carries its own count in its stream meta. */
    fieldMeta["ncount"] = (Json::Value::UInt)excClasses[(uint8_t)'N'].count;
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

int32_t SamCodecActuator::compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::chrono::steady_clock::time_point qualT0, qualT1, qualT2;
    qualT0 = pbgzprof::nowOrZero();
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

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
    const bool useFcv2 = (pickedQualCoder == CoderType::FCV2);
    const bool useBwtCm = (pickedQualCoder == CoderType::BWT_CM);
    std::shared_ptr<coder_qual> qualityCoder;
    std::shared_ptr<coder_fcv2> fcv2Coder;
    std::shared_ptr<coder_bwt_cm> qualCmCoder;
    if (useFcv2) {
        /*
         * If preprocessing picks fcv2, it also hands over the context parameter
         * tier; the encoding side creates the coder with the same tier used for
         * prior training, and the tier is written into the stream header for the
         * decoding side to read back. Without preprocessing, it falls back to
         * the default tier, consistent with cold-start behavior.
         */
        Fcv2Cfg fcv2Cfg;
        const FieldCodecSelection* qualSel =
            (qualPreInfo != nullptr) ? qualPreInfo->getField(SAM_QUAL) : nullptr;
        if (qualSel != nullptr && qualSel->selectedCoder == CoderType::FCV2) {
            fcv2Cfg.cycleMax = qualSel->fcv2Params.cycleMax;
            fcv2Cfg.cycleBucket = qualSel->fcv2Params.cycleBucket;
            fcv2Cfg.deltaMax = qualSel->fcv2Params.deltaMax;
            fcv2Cfg.deltaBucket = qualSel->fcv2Params.deltaBucket;
            fcv2Cfg.prevShift = qualSel->fcv2Params.prevShift;
            fcv2Cfg.useDelta = qualSel->fcv2Params.useDelta;
            fcv2Cfg.useDedup = qualSel->fcv2Params.useDedup;
            fcv2Cfg.useQa = qualSel->fcv2Params.useQa;
            fcv2Cfg.modelCount = qualSel->fcv2Params.modelCount;
        }
        std::vector<uint32_t> freqByByte(256, 0);
        for (uint32_t i = 0; i < qualFreqTable.size(); ++i) {
            uint32_t b = (uint32_t)qualFreqTable[i].first + (uint32_t)'!';
            if (b < 256) {
                /* qualFreqTable keeps only the alphabet in frequency-descending
                   order; the actual counts were not retained. A monotonically
                   decreasing weight is derived from the rank here, used only to
                   shape the Huffman tree. */
                freqByByte[b] = (uint32_t)(qualFreqTable.size() - i);
            }
        }
        /*
         * When a prior exists, start from it: the model need not relearn from
         * fixed initial values, and the gain is larger for smaller blocks. The
         * prior's absolute address must be written into the block meta—the
         * decoding side can only match up if it obtains the same snapshot, and
         * in random-access scenarios it cannot infer this address by walking the
         * sequential stream.
         *
         * A load failure on the compression side is a benign fallback (the
         * fixed-initial model is kept, only the ratio suffers), so the address
         * is registered only when loaded is true, avoiding a promise the
         * decoding side cannot honor.
         */
        qualPriorAddress = (pbgzEngine != nullptr) ? pbgzEngine->getQualPriorAddress() : -1;
        AuxPayloadPtr priorBlob =
            (pbgzEngine != nullptr) ? pbgzEngine->getQualPrior(0) : AuxPayloadPtr();
        bool priorLoaded = false;
        pbgzprof::addSince(pbgzprof::QUAL_PREP, qualT0);
        qualT1 = pbgzprof::nowOrZero();
        if (priorBlob && !priorBlob->empty()) {
            fcv2Coder = std::make_shared<coder_fcv2>(qualityIo.get(), freqByByte, fcv2Cfg,
                                                     *priorBlob, &priorLoaded);
        } else {
            fcv2Coder = std::make_shared<coder_fcv2>(qualityIo.get(), freqByByte, fcv2Cfg);
        }
        if (!priorLoaded) {
            qualPriorAddress = -1;
        }
    } else if (useBwtCm) {
        /*
         * bwt_cm needs no alphabet and no SEQ or strand direction; feed quality
         * values record by record. It is the second choice when fcv2 does not
         * apply: empirically 7.4 percentage points better than coder_qual.
         */
        qualT1 = pbgzprof::nowOrZero();   /* nothing to prepare on this path */
        qualCmCoder = std::make_shared<coder_bwt_cm>(qualityIo.get());
        CoderFactory::applyLevel(qualityIo.get(), CoderType::BWT_CM, engineCompressLevel());
    } else {
        qualT1 = pbgzprof::nowOrZero();
        qualityCoder = std::make_shared<coder_qual>(qualityIo.get(), true, qualFreqTable);
    }
    /* Every path above starts its measurement, so this is the constructor's cost on all of them. */
    pbgzprof::addSince(pbgzprof::QUAL_CTOR, qualT1);
    qualT2 = pbgzprof::nowOrZero();

    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;
    Json::Value streamMeta;

    // Encode quality data
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

        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* qualStart = line + prevTabPos + 1;
        uint32_t qualLength = currTabPos - prevTabPos - 1;

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

        if (useFcv2) {
            /*
             * Take the strand direction from the FLAG field. Per the SAM spec,
             * bit 0x10 being set means SEQ and QUAL are stored relative to the
             * reference's forward strand, i.e. reversed relative to the order
             * read out by the sequencer, and the coder must use this to restore
             * the real cycle number. The field is parsed by hand rather than
             * with strtol because it is not NUL-terminated and this is a hot
             * path taken for every record.
             */
            bool rev = false;
            if (contentPos[contentIdx].size() > 1) {
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
            fcv2Coder->encode_record(qualStart, qualLength, rev, seqStart, seqLen);
        } else if (useBwtCm) {
            qualCmCoder->encode_line(qualStart, qualLength);
        } else {
            // Encode quality with sequence context
            qualityCoder->encode_qual_gen2(seqStart, qualStart, qualLength);
        }
        streamSrcLen += qualLength;
    }

    if (useFcv2) {
        fcv2Coder->encode_flush();
    } else if (useBwtCm) {
        qualCmCoder->encode_flush();
    } else {
        qualityCoder->encode_flush();
    }
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

    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        uint8_t* basePtr = nullptr;
        uint32_t actualBaseLen = 0;
        // Decode each field for this line
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
            Json::Value& idMeta = streamMeta[idx];
            std::string idSplit = idMeta["splitsym"].asString();
            for (uint32_t i = 0; i < idSplit.length(); ++i) {
                idSplitSymbols.push_back(idSplit.c_str()[i]);
            }
            Json::Value& idStreamMeta = idMeta["streams"];
            if (idStreamMeta.size() != idSplitSymbols.size()) {
                LOG_ERROR("id streams not match id split, expect:%u, actual:%u", idSplitSymbols.size(), idStreamMeta.size());
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
                if (coderName == "coder_affix_match") {
                    idDecoders.push_back(std::make_shared<coder_affix_match>(io.get()));
                    idDecoders.back()->set_level(idStreamMeta[i]["coder"]["level"].asInt());
                } else if (coderName == "coder_bwt_cm") {
                    idDecoders.push_back(std::make_shared<coder_bwt_cm>(io.get()));
                    idDecoders.back()->set_level(idStreamMeta[i]["coder"]["level"].asInt());
                } else if (coderName == "coder_qname") {
                    idUsesQnameCoder = true;
                    idDecoders.push_back(std::make_shared<coder_qname>(io.get()));
                } else {
                    LOG_ERROR("Unsupport coder name:%s", coderName.c_str());
                    return -1;
                }

                /*
                 * Numeric sub-streams carry no per-line terminator, so they cannot
                 * be decoded one segment per line the way the textual path does.
                 * Decode the whole varint stream once here and serve one value per
                 * line from idNumericPos (see decompressIdField).
                 */
                const bool isNumeric = idStreamMeta[i].isMember("mode") &&
                                       idStreamMeta[i]["mode"].asString() == "numeric";
                uint8_t* nbuf = nullptr;
                uint32_t nlen = 0;
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
                }
                idNumericBufs.push_back(nbuf);
                idNumericLens.push_back(nlen);
                idNumericPos.push_back(0);
                idNumericAcc.push_back(0);

                readOffset += dstLength;
            }
        } else if (idx == 9) {
            Json::Value& baseMeta = streamMeta[idx];
            maxBaseLength = baseMeta["maxlen"].asUInt();
            minBaseLength = baseMeta["minlen"].asUInt();
            LOG_DEBUG("maxBaseLen = %d, minBaseLen = %d", maxBaseLength, minBaseLength);
            bool isUseReference = pRefeGene != nullptr && baseMeta.isMember("streams");
            if (!isUseReference) {
                // Scenario without reference genome
                std::string coderName = streamMeta[idx]["coder"]["magic"].asString();
                uint32_t dstLength = streamMeta[idx]["totaldstlen"].asUInt();
                uint32_t srcLength = streamMeta[idx]["totalsrclen"].asUInt();
                LOG_DEBUG("srclen = %d, dstlen = %d", srcLength, dstLength);
                if (coderName == "coder_fc") {
                    /*
                     * The whole-block SEQ is decoded into the tail of the
                     * outputBlock buffer as staging, then decompressBase moves
                     * it to the head line by line. No separate buffer is
                     * allocated: that would add one malloc/free of the whole SEQ
                     * per block plus page faults on first touch, and raise peak
                     * memory by threads x SEQ size, while the number of copies
                     * stays the same.
                     *
                     * Invariant: never realloc outputBlock within a block, or
                     * this tail pointer would dangle. Capacity is guaranteed by
                     * the one-shot block_size*2 pre-allocation at the block
                     * entry (see decompress())—head output <= block_size and
                     * tail staging <= block_size, exactly 2x.
                     */
                    coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstLength, &ioErrSink, "SEQ");
                    baseIo.meta = baseMeta;
                    baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
                    coder_fc baseDecoder = coder_fc(&baseIo);
                    if (baseDecoder.decode_line(outputBlock->getBuffer() + outputBlock->getBufferSize() - srcLength, srcLength, UINT8_MAX, false) < 0) {
                        LOG_ERROR("Decode SEQ by coder_fc failed, srclen = %u", srcLength);
                        return -1;
                    }
                } else if(coderName == "coder_bwt_cm") {
                    std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "SEQ");
                    ioVector.push_back(io);
                    fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
                } else {
                    LOG_ERROR("Unsupported coder name:%s", coderName.c_str());
                    return -1;
                }
                readOffset += dstLength;
            } else {
                uint32_t id = 0;
                // Scenario using reference genome
                Json::Value& baseMetaStreams = baseMeta["streams"];
                if (baseMetaStreams[id]["sname"] != "m") {
                    LOG_ERROR("check sub stream failed:%s", baseMetaStreams[id]["sname"].asString().c_str());
                    return -1;
                }

                uint32_t dstLength = baseMetaStreams[id]["dstlen"].asUInt();
                const std::string matchCoderName = baseMetaStreams[id]["coder"]["magic"].asString();
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
                    /*
                     * RLE layout written by compressBaseWithRef (two sub-streams):
                     *   "m"    -> varint run lengths of the zero runs
                     *   "mval" -> the surviving non-zero values (1 byte each)
                     * Expand them back into the original sparse byte stream
                     * (orgrawlen bytes) and serve it per record via matchBlockDecode.
                     */
                    uint32_t runLen = baseMetaStreams[id]["srclen"].asUInt();
                    uint32_t srcLen = baseMetaStreams[id]["orgrawlen"].asUInt();
                    std::shared_ptr<coder> runDecoder;
                    if (matchCoderName == "coder_bwt_cm") {
                        runDecoder = std::make_shared<coder_bwt_cm>(io.get());
                    } else if (matchCoderName == "coder_fc") {
                        runDecoder = std::make_shared<coder_fc>(io.get());
                    } else {
                        LOG_ERROR("check sub stream failed, coder name not match: %s", matchCoderName.c_str());
                        return -1;
                    }
                    std::unique_ptr<uint8_t[]> runBuf = std::make_unique<uint8_t[]>(runLen + 1);
                    if (runDecoder->decode_line(runBuf.get(), runLen, UINT8_MAX, false) < 0) {
                        LOG_ERROR("Decode SEQ match run stream failed");
                        return -1;
                    }

                    /*
                     * The value sub-stream (if any) sits immediately after the run
                     * sub-stream. A fully reference-matching block has no surviving
                     * non-zero values, so compressBaseWithRef writes no "mval"
                     * sub-stream at all; id then stays on the run stream so the
                     * npos / baselen sub-streams that follow are addressed as if the
                     * value stream never existed.
                     */
                    uint32_t valDstLen = 0;
                    uint32_t valLen = 0;
                    std::unique_ptr<uint8_t[]> valBuf;
                    if (id + 1 < baseMetaStreams.size() &&
                        baseMetaStreams[id + 1]["sname"].asString() == "mval") {
                        id++;                       /* step past "m" onto the "mval" sub-stream */
                        const uint32_t valId = id;
                        valDstLen = baseMetaStreams[valId]["dstlen"].asUInt();
                        valLen = baseMetaStreams[valId]["srclen"].asUInt();
                        std::shared_ptr<coder_io> valIo = makeCoderIo(
                            inBlockPtr->getBuffer() + readOffset + dstLength, valDstLen, "SEQ match val");
                        valIo->meta = baseMetaStreams[valId];
                        ioVector.push_back(valIo);
                        /* Each sub-stream carries its own coder name: the two halves of
                           the RLE layout are written independently, so the value stream
                           must not be decoded with the run stream's coder. */
                        const std::string valCoderName = baseMetaStreams[valId]["coder"]["magic"].asString();
                        std::shared_ptr<coder> valDecoder;
                        if (valCoderName == "coder_bwt_cm") {
                            valDecoder = std::make_shared<coder_bwt_cm>(valIo.get());
                        } else if (valCoderName == "coder_fc") {
                            valDecoder = std::make_shared<coder_fc>(valIo.get());
                        } else {
                            LOG_ERROR("check sub stream failed, coder name not match: %s", valCoderName.c_str());
                            return -1;
                        }
                        valBuf = std::make_unique<uint8_t[]>(valLen + 1);
                        if (valDecoder->decode_line(valBuf.get(), valLen, UINT8_MAX, false) < 0) {
                            LOG_ERROR("Decode SEQ match value stream failed");
                            return -1;
                        }
                        matchExtraOffset = valDstLen;   /* skipped past by the caller below */
                    } else {
                        /* No surviving values in this fully-matching block. */
                        matchExtraOffset = 0;
                    }

                    /* Expand: each run is followed by one surviving value. calloc zero-fills,
                       so the expanded zeros need no explicit write. */
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
                } else if (matchCoderName == "coder_bwt_cm") {
                    /* Legacy layout (no rle member): decode line by line, as before. */
                    fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
                    matchBlockDecode = false;
                } else if (matchCoderName == "coder_fc") {
                    /* coder_fc only supports block decompression, so decode the whole
                       match stream here and slice it per record in decompressBase. */
                    fieldDecoders[idx] = std::make_shared<coder_fc>(io.get());
                    matchBlockDecode = true;
                    matchBlockLength = baseMetaStreams[id]["srclen"].asUInt();
                    matchBlockOffset = 0;
                    MemoryUtil::safeFree(matchBlockBuffer);
                    matchBlockBuffer = MemoryUtil::safeAlloc<uint8_t>(matchBlockLength + 1);
                    if (matchBlockBuffer == nullptr) {
                        LOG_ERROR("Alloc SEQ match block buffer failed");
                        return -1;
                    }
                    if (fieldDecoders[idx]->decode_line(matchBlockBuffer, matchBlockLength,
                                                        UINT8_MAX, false) < 0) {
                        LOG_ERROR("Decode SEQ match stream (coder_fc) failed");
                        return -1;
                    }
                } else {
                    LOG_ERROR("check sub stream failed, coder name not match: %s", matchCoderName.c_str());
                    return -1;
                }

                readOffset += dstLength + matchExtraOffset;   /* +mval sub-stream under RLE */
                matchExtraOffset = 0;
                /*
                 * SEQ exception streams (see SeqExceptionClass): one per character present,
                 * each a strictly increasing list of block offsets, each naming its character
                 * in "ch". A character that never occurs has no stream, so this reads streams
                 * until one turns up that is not an exception stream; the per-record refill
                 * then walks the lists it has collected.
                 *
                 * The one older archive still read is the first layout's lone "npos": a single
                 * list of absolute offsets holding every position that layout knew of, 'N' and
                 * 'n' alike, with no "ch" and with its count in the field-level "ncount". It
                 * wrote an 'N' back at each of those positions, and so does this.
                 */
                seqExc.clear();
                for (;;) {
                    const std::string name = (baseMetaStreams.size() > (size_t)id + 1)
                                                 ? baseMetaStreams[id + 1]["sname"].asString()
                                                 : std::string();
                    /* Every form name identifies an exception stream, so the block's exception
                       streams are the ones that come next and this stops at the first other name.
                       What the name and the rest of the meta say is read in one place (see
                       seqExcLayoutOf), which is also where the older layouts are handled. */
                    if (seqExcFormOf(name) == SeqExcForm::None) {
                        break;   /* this block has no further exception stream */
                    }
                    id++;
                    const SeqExcLayout layout = seqExcLayoutOf(baseMeta, baseMetaStreams[id]);
                    const uint8_t ch = layout.ch;
                    uint32_t dstlen = baseMetaStreams[id]["dstlen"].asUInt();
                    uint32_t srclen = baseMetaStreams[id]["srclen"].asUInt();
                    if (baseMetaStreams[id]["coder"]["magic"].asString() != "coder_bwt_cm") {
                        LOG_ERROR("check sub stream failed. coder not match. coder = %s.",
                            baseMetaStreams[id]["coder"]["magic"].asString().c_str());
                        return -1;
                    }

                    coder_io excIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ npos");
                    auto excCoder = std::make_unique<coder_bwt_cm>(&excIo);
                    /* The payload carries the positions alone - the character travels in the
                       stream's meta ("ch") - and every form is expanded by the same reader. */
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
                }

                if (minBaseLength != maxBaseLength) {
                    id++;
                    if (baseMetaStreams[id]["sname"].asString() != "baselen") {
                        LOG_ERROR("check sub stream failed. sname not match");
                        return -1;
                    }

                    uint32_t dstlen = baseMetaStreams[id]["dstlen"].asUInt();
                    uint32_t srclen = baseMetaStreams[id]["srclen"].asUInt();
                    uint8_t* baseLenBuffer = MemoryUtil::safeAlloc<uint8_t>(srclen);

                    if (baseMetaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
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
                            baseMetaStreams[id]["coder"]["magic"].asString().c_str());
                        return -1;
                    }
                    /* The delta layout: the entry ordinals as one run of forward deltas, then
                       the lengths (see the encoder). An archive written before it holds
                       absolute 4-byte (ordinal, length) pairs instead. */
                    if (baseMetaStreams[id]["delta"].isUInt()) {
                        const uint32_t count = baseMetaStreams[id]["count"].asUInt();
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
                }
            }
        } else if (idx == 10) {
            Json::Value& qualMeta = streamMeta[idx];
            Json::Value& qualStreamMeta = qualMeta["streams"];
            if (qualStreamMeta.size() != 2) {
                LOG_ERROR("quality streams check failded, size = %d", qualStreamMeta.size());
                return -1;
            }
            if (qualStreamMeta[1]["coder"]["magic"] == "coder_bwt_cm") {
                uint32_t qualDstLength = qualStreamMeta[0]["dstlen"].asUInt();
                uint32_t freqDstLength = qualStreamMeta[1]["dstlen"].asUInt();

                coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualDstLength, freqDstLength, &ioErrSink, "QUAL freq table");
                auto qualFreqCoder = std::make_unique<coder_bwt_cm>(&qualFreqIo);
                uint32_t qualFreqSrcLength = qualStreamMeta[1]["srclen"].asUInt();
                /* Same as fastq_actuator: counting with uint8_t wraps when the alphabet exceeds 127 symbols, causing a heap out-of-bounds write. */
                uint32_t qualFreqArrLength = qualFreqSrcLength / sizeof(uint16_t);
                uint16_t* qualFreqArr = new uint16_t[qualFreqArrLength];
                int32_t qualFreq = qualFreqCoder->decode_line((uint8_t*)qualFreqArr,qualFreqSrcLength, UINT8_MAX, false);
                if (qualFreq < 0 || (uint32_t)qualFreq != qualFreqSrcLength) {
                    LOG_ERROR("Decode quality frequncy failed");
                    delete [] qualFreqArr;
                    return -1;
                }
                for (uint32_t i = 0; i < qualFreqArrLength; i += 2) {
                    qualFreqTable.push_back(std::make_pair(qualFreqArr[i], qualFreqArr[i + 1]));
                }
                delete [] qualFreqArr;

                if (qualStreamMeta[0]["coder"]["magic"].asString() == "coder_qual") {
                    std::shared_ptr<coder_io> qualIo = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualDstLength, "QUAL");
                    ioVector.push_back(qualIo);
                    qualCoder = std::make_shared<coder_qual>(qualIo.get(), true, qualFreqTable);
                } else if (qualStreamMeta[0]["coder"]["magic"].asString() == "coder_fcv2") {
                    /*
                     * fcv2's alphabet and per-symbol frequencies are written in
                     * its own stream header; begin_decode reads them back to
                     * rebuild the Huffman tree, so passing an empty frequency
                     * table here at construction is enough.
                     */
                    std::shared_ptr<coder_io> qualIo = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualDstLength, "QUAL");
                    ioVector.push_back(qualIo);
                    std::vector<uint32_t> emptyFreq(256, 0);
                    /*
                     * If the encoding side started from a prior, the decoding
                     * side must start from the same snapshot, or the model
                     * diverges immediately. Unlike the benign fallback on the
                     * compression side, a load failure here must be fatal:
                     * silently falling back to fixed initial values would decode
                     * a stream of seemingly valid wrong data, far more dangerous
                     * than failing outright.
                     *
                     * The prior's address comes from the block meta rather than
                     * sequential inference, so random access works as well.
                     */
                    if (qualStreamMeta[0].isMember("prior")) {
                        /*
                         * The offset in meta is only a seek handle and
                         * validation value; the index key is the package index:
                         * under piped input the absolute offset degenerates to
                         * 0, and looking it up by offset would always miss.
                         */
                        const int64_t priorAddress = qualStreamMeta[0]["prior"].asInt64();
                        AuxPayloadPtr priorBlob =
                            (pbgzEngine != nullptr)
                                ? pbgzEngine->getQualPrior(inBlockPtr->getPackageIndex())
                                : AuxPayloadPtr();
                        if (!priorBlob || priorBlob->empty()) {
                            LOG_ERROR("Qual prior at offset %lld required but unavailable",
                                      (long long)priorAddress);
                            return -1;
                        }
                        bool priorLoaded = false;
                        qualFcv2Decoder = std::make_shared<coder_fcv2>(qualIo.get(), emptyFreq,
                                                                       *priorBlob, &priorLoaded);
                        if (!priorLoaded) {
                            LOG_ERROR("Qual prior at offset %lld failed to load",
                                      (long long)priorAddress);
                            return -1;
                        }
                    } else {
                        qualFcv2Decoder = std::make_shared<coder_fcv2>(qualIo.get(), emptyFreq);
                    }
                    const int32_t fcv2Ret = qualFcv2Decoder->begin_decode();
                    if (fcv2Ret == coder_ns::CODER_ERR_UNSUPPORTED_VERSION) {
                        /* Neither the current header layout nor the one from before the version
                           byte fits this stream, so its layout is unknown - naming that is more
                           useful than reporting it as corruption (see FCV2_STREAM_VERSION). */
                        LOG_ERROR("QUAL stream fits neither the current fcv2 header layout (v%u) nor "
                                  "the layout from before it: the archive was written by an "
                                  "incompatible version and cannot be decoded",
                                  (unsigned)FCV2_STREAM_VERSION);
                        return -1;
                    }
                    if (fcv2Ret != 0) {
                        LOG_ERROR("fcv2 begin_decode failed");
                        return -1;
                    }
                } else if (qualStreamMeta[0]["coder"]["magic"].asString() == "coder_bwt_cm") {
                    std::shared_ptr<coder_io> qualIo = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualDstLength, "QUAL");
                    ioVector.push_back(qualIo);
                    qualCmDecoder = std::make_shared<coder_bwt_cm>(qualIo.get());
                } else {
                    LOG_ERROR("Unsupport coder type: %s", streamMeta[0]["coder"]["magic"].asString().c_str());
                    return -1;
                }
                readOffset += (qualDstLength + freqDstLength);
            } else {
                LOG_ERROR("Unsupport coder type: %s", streamMeta[1]["coder"]["magic"].asString().c_str());
                return -1;
            }
        } else if (idx == 8) {
            /* TLEN field: prefer the reconstruction optimization; the exception stream is decoded separately. */
            Json::Value& tlenMeta = streamMeta[idx];
            bool isOptimized = tlenMeta.isMember("optimized") && tlenMeta["optimized"].asBool();

            if (isOptimized) {
                if (tlenMeta.isMember("streams") && tlenMeta["streams"].isArray()) {
                    Json::Value& tlenStreams = tlenMeta["streams"];
                    for (Json::Value::iterator it = tlenStreams.begin(); it != tlenStreams.end(); ++it) {
                        Json::Value& stream = *it;
                        if (stream["sname"].asString() != "tlenexc") {
                            continue;
                        }
                        uint32_t srclen = stream["srclen"].asUInt();
                        uint32_t dstlen = stream["dstlen"].asUInt();
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
                        /*
                         * The layout the stream is in (see the marker note on tlenDecodeVarints).
                         * The marker is authoritative where it is there. Without one the archive
                         * is older than the marker, which means the fixed layout - unless the
                         * stream's size rules that out, in which case the archive was written
                         * between the change to varint and the marker and holds varints.
                         */
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
                        if (decodedOk) {
                            tlenCache.insert(decoded.begin(), decoded.end());
                        } else {
                            MemoryUtil::safeFree(excBuffer);
                            LOG_ERROR("Corrupt TLEN exception stream: neither the varint layout nor the "
                                      "fixed-pair layout fits (%u bytes, %u exceptions)", srclen, declared);
                            return -1;
                        }
                        MemoryUtil::safeFree(excBuffer);
                        readOffset += dstlen;
                    }
                }
            } else {
                std::string coderName = tlenMeta["coder"]["magic"].asString();
                uint32_t dstLen = tlenMeta["dstlen"].asUInt();
                if (coderName != "coder_bwt_cm") {
                    LOG_ERROR("Unsupported TLEN coder type: %s", coderName.c_str());
                    return -1;
                }
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "TLEN");
                ioVector.push_back(io);
                fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
                readOffset += dstLen;
            }
        } else if (idx == 11 && streamMeta[idx].isMember("mode") &&
                   streamMeta[idx]["mode"].asString() == "tag_split") {
            /*
             * OPTION tag-split: no decoder is built in initDecoder and
             * readOffset is not advanced—the start of this field's stream is
             * recorded, and decompressOptionField lazily decodes all columns
             * from the streams array when the first OPTION line is reached.
             * readOffset cannot be used directly: line-by-line decoders such as
             * decompressIdField advance readOffset, so it has already shifted by
             * the time the OPTION line is reached (measured as off by the ID
             * field's dstlen).
             */
            fieldIoStart[idx] = readOffset;
            continue;
        } else if (idx == 7 && streamMeta[idx].isMember("mode") &&
                   streamMeta[idx]["mode"].asString() == "pnext_qname_rebuild") {
            /*
             * PNEXT qname-rebuild mode stores only an exception stream in
             * streams[7]["streams"][0]. No line-by-line decoder is built; the
             * exception stream offset is recorded so preDecodeForTLEN can decode
             * it, and readOffset is advanced past it so later fields align.
             */
            fieldIoStart[idx] = readOffset;
            if (streamMeta[idx].isMember("streams") && streamMeta[idx]["streams"].size() > 0) {
                uint32_t excDstLen = streamMeta[idx]["streams"][0]["dstlen"].asUInt();
                readOffset += excDstLen;
            }
            continue;
        } else {
            std::string coderName = streamMeta[idx]["coder"]["magic"].asString();
            uint32_t dstLen = streamMeta[idx]["dstlen"].asUInt();
            /*
             * Record the compressed stream positions of POS(3)/CIGAR(5)/PNEXT(7)
             * so that preDecodeForTLEN can rebuild decoders and pre-decode the
             * whole block; only then can TLEN reconstruction obtain the full
             * mate index.
             */
            fieldIoStart[idx] = readOffset;
            fieldIoDstLen[idx] = dstLen;
            if (streamMeta[idx]["coder"].isMember("level")) {
                fieldIoLevel[idx] = streamMeta[idx]["coder"]["level"].asInt();
            }
            if (coderName == "coder_bwt_cm") {
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "SAM field");
                ioVector.push_back(io);
                fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
            } else if (coderName == "coder_affix_match") {
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "SAM field");
                ioVector.push_back(io);
                fieldDecoders[idx] = std::make_shared<coder_affix_match>(io.get());
                fieldDecoders[idx]->set_level(streamMeta[idx]["coder"]["level"].asInt());
            } else if (coderName == "coder_arith") {
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLen, "SAM field");
                ioVector.push_back(io);
                fieldDecoders[idx] = std::make_shared<coder_arith>(io.get());
                if (streamMeta[idx]["coder"].isMember("level")) {
                    fieldDecoders[idx]->set_level(streamMeta[idx]["coder"]["level"].asInt());
                }
                /* Apply the file-level POS delta prior. The encoder writes a
                 * prior only when arith won the POS trial and the file is
                 * large enough to make it pay; on smaller files it cold-starts
                 * from the uniform table, so a missing prior here must fall
                 * back to uniform too, matching the encoder exactly. */
                AuxPayloadPtr posPrior =
                    (pbgzEngine != nullptr) ? pbgzEngine->getPosPrior() : AuxPayloadPtr();
                if (posPrior && !posPrior->empty()) {
                    static_cast<coder_arith*>(fieldDecoders[idx].get())
                        ->set_prior(posPrior->data(), (uint32_t)posPrior->size());
                }
            } else {
                LOG_ERROR("Unsupport coder type: %s", streamMeta[0]["coder"]["magic"].asString().c_str());
                return -1;
            }
            readOffset += dstLen;
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
int32_t SamCodecActuator::preDecodeForTLEN() {
    if (samLine == 0 || !meta.isMember("sam")) {
        return 0;
    }

    uint8_t* buffer = inBlockPtr->getBuffer();
    uint8_t tmpBuf[1024];
    Json::Value& streams = meta["sam"]["streams"];

    /*
     * Pre-decode QNAME (field 0) into decodedQnames. PNEXT (pnext_qname_rebuild
     * mode) is rebuilt by pairing records with the same QNAME, so every line's
     * QNAME must be known before PNEXT is reconstructed. QNAME is decoded with
     * the per-sub-stream decoders already built by initDecoder (idDecoders),
     * mirroring decompressIdField.
     */
    decodedQnames.clear();
    decodedQnames.resize(samLine);
    if (!idStreamOffsets.empty()) {
        // Rebuild independent decoders from the recorded sub-stream offsets so we
        // do not consume the idDecoders used by the main loop.
        std::vector<std::shared_ptr<coder_io>> preIdIo;
        std::vector<std::shared_ptr<coder>> preIdDec;
        for (uint32_t si = 0; si < idStreamOffsets.size(); ++si) {
            std::shared_ptr<coder_io> io = makeCoderIo(buffer + idStreamOffsets[si], idStreamDstLens[si], "QNAME predecode");
            preIdIo.push_back(io);
            int32_t lvl = -1;
            if (streams[0].isMember("streams") && streams[0]["streams"].isValidIndex(si) &&
                streams[0]["streams"][si]["coder"].isMember("level")) {
                lvl = streams[0]["streams"][si]["coder"]["level"].asInt();
            }
            if (idStreamCoders[si] == "coder_affix_match") {
                auto dec = std::make_shared<coder_affix_match>(io.get());
                if (lvl >= 0) dec->set_level(lvl);
                preIdDec.push_back(dec);
            } else if (idStreamCoders[si] == "coder_bwt_cm") {
                auto dec = std::make_shared<coder_bwt_cm>(io.get());
                if (lvl >= 0) dec->set_level(lvl);
                preIdDec.push_back(dec);
            } else if (idStreamCoders[si] == "coder_qname") {
                preIdDec.push_back(std::make_shared<coder_qname>(io.get()));
            } else {
                preIdDec.push_back(nullptr);
            }
        }

        /* Numeric sub-streams have no per-line terminator: decode the whole varint
           stream up front and serve one value per line, mirroring decompressIdField. */
        std::vector<uint8_t*> numBufs(idStreamOffsets.size(), nullptr);
        std::vector<uint32_t> numLens(idStreamOffsets.size(), 0);
        std::vector<uint32_t> numPos(idStreamOffsets.size(), 0);
        std::vector<uint64_t> numAcc(idStreamOffsets.size(), 0);
        for (uint32_t si = 0; si < idStreamOffsets.size(); ++si) {
            if (si >= preIdDec.size() || !preIdDec[si]) {
                continue;
            }
            if (!(streams[0].isMember("streams") && streams[0]["streams"].isValidIndex(si) &&
                  streams[0]["streams"][si].isMember("mode") &&
                  streams[0]["streams"][si]["mode"].asString() == "numeric")) {
                continue;
            }
            const uint32_t srclen = streams[0]["streams"][si]["srclen"].asUInt();
            uint8_t* b = MemoryUtil::safeAlloc<uint8_t>(srclen + 1);
            if (b == nullptr) {
                return -1;
            }
            int32_t dl = preIdDec[si]->decode_line(b, srclen, UINT8_MAX, false);
            if (dl < 0) {
                MemoryUtil::safeFree(b);
                LOG_ERROR("Predecode id numeric sub-stream(%u) failed: %d", si, dl);
                return -1;
            }
            numBufs[si] = b;
            numLens[si] = (uint32_t)dl;
        }

        char qnameBuf[512];
        for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
            std::string qname;
            if (idUsesQnameCoder && !preIdDec.empty()) {
                // Single-stream coder_qname: one QNAME per line (with trailing '\t').
                bool hold = (idStreamCoders[0] == "coder_affix_match");
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
                    if (!preIdDec[si]) continue;
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
                    uint8_t sep = (si < idSplitSymbols.size()) ? idSplitSymbols[si] : UINT8_MAX;
                    // coder_affix_match keeps cross-line context in an internal
                    // buffer only when need2hold is set; without it, it points
                    // last at our fixed qnameBuf and the next line would corrupt
                    // the prefix reference.
                    bool hold = (idStreamCoders[si] == "coder_affix_match");
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
    }

    /*
     * Decode FLAG (field 1) first so that PNEXT (field 7) can be decoded
     * correctly: PNEXT is stored as delta against POS only when the mate
     * position is meaningful (FLAG bit 0x1 set and bit 0x8 clear), otherwise
     * it stores the original value. This mirrors compressPNextFieldDelta.
     */
    auto flagStartIt = fieldIoStart.find(1);
    if (flagStartIt != fieldIoStart.end() && streams.isValidIndex(1) && streams[1].isMember("coder")) {
        uint32_t off = flagStartIt->second;
        uint32_t dstlen = fieldIoDstLen[1];
        std::string coderName = streams[1]["coder"]["magic"].asString();
        std::string mode = streams[1].isMember("mode") ? streams[1]["mode"].asString() : "";

        std::shared_ptr<coder_io> tmpIo = makeCoderIo(buffer + off, dstlen, "FLAG predecode");
        std::shared_ptr<coder> tmpDec;
        if (coderName == "coder_bwt_cm") {
            tmpDec = std::make_shared<coder_bwt_cm>(tmpIo.get());
        } else if (coderName == "coder_affix_match") {
            tmpDec = std::make_shared<coder_affix_match>(tmpIo.get());
        } else {
            tmpDec = nullptr;
        }
        if (tmpDec) {
            auto lvIt = fieldIoLevel.find(1);
            if (lvIt != fieldIoLevel.end()) {
                tmpDec->set_level(lvIt->second);
            }

            bool binary = (mode == "number" || mode == "");
            for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
                bool need2hold = (coderName == "coder_affix_match");
                int32_t len;
                if (binary) {
                    len = tmpDec->decode_line(tmpBuf, sizeof(uint16_t), UINT8_MAX, need2hold);
                    if (len >= (int32_t)sizeof(uint16_t)) {
                        mappedFlag[lineNo] = *(uint16_t*)tmpBuf;
                    } else {
                        mappedFlag[lineNo] = 0;
                    }
                } else {
                    len = tmpDec->decode_line(tmpBuf, (uint32_t)sizeof(tmpBuf), '\t', need2hold);
                    if (len > 1) {
                        std::string flagStr((char*)tmpBuf, (size_t)len - 1);
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
        }
    }

    /*
     * Pre-decode RNAME (field 2) so the POS delta chain below can detect
     * chromosome switches from mappedChr, mirroring the encoder. A separate
     * decoder is used so the main loop's RNAME decoder state is untouched.
     */
    auto rnameStartIt = fieldIoStart.find(2);
    if (rnameStartIt != fieldIoStart.end() && streams.isValidIndex(2) && streams[2].isMember("coder")) {
        uint32_t rnOff = rnameStartIt->second;
        uint32_t rnDst = fieldIoDstLen[2];
        std::string rnCoderName = streams[2]["coder"]["magic"].asString();
        std::shared_ptr<coder_io> rnIo = makeCoderIo(buffer + rnOff, rnDst, "RNAME predecode");
        std::shared_ptr<coder> rnDec;
        if (rnCoderName == "coder_bwt_cm") {
            rnDec = std::make_shared<coder_bwt_cm>(rnIo.get());
        } else if (rnCoderName == "coder_affix_match") {
            rnDec = std::make_shared<coder_affix_match>(rnIo.get());
        }
        if (rnDec) {
            auto rnLvIt = fieldIoLevel.find(2);
            if (rnLvIt != fieldIoLevel.end()) {
                rnDec->set_level(rnLvIt->second);
            }
            bool rnHold = (rnCoderName == "coder_affix_match");
            for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
                uint16_t chrIndex = 0;
                int32_t len = rnDec->decode_line((uint8_t*)&chrIndex, sizeof(chrIndex), UINT8_MAX, rnHold);
                mappedChr[lineNo] = (len >= (int32_t)sizeof(chrIndex)) ? chrIndex : 0xFFFF;
            }
        }
    }

    /* Fields needed to reconstruct TLEN: POS(3), CIGAR(5), PNEXT(7). */
    static const uint32_t tlenFields[] = {3, 5, 7};

    for (uint32_t f : tlenFields) {
        auto startIt = fieldIoStart.find(f);
        if (startIt == fieldIoStart.end() || !streams.isValidIndex(f) || !streams[f].isMember("coder")) {
            continue;
        }

        uint32_t off = startIt->second;
        uint32_t dstlen = fieldIoDstLen[f];
        std::string coderName = streams[f]["coder"]["magic"].asString();
        std::string mode = streams[f].isMember("mode") ? streams[f]["mode"].asString() : "";

        std::shared_ptr<coder_io> tmpIo = makeCoderIo(buffer + off, dstlen, "TLEN predecode");
        std::shared_ptr<coder> tmpDec;
        if (coderName == "coder_bwt_cm") {
            tmpDec = std::make_shared<coder_bwt_cm>(tmpIo.get());
        } else if (coderName == "coder_affix_match") {
            tmpDec = std::make_shared<coder_affix_match>(tmpIo.get());
        } else if (coderName == "coder_arith") {
            tmpDec = std::make_shared<coder_arith>(tmpIo.get());
        } else {
            continue;
        }
        auto lvIt = fieldIoLevel.find(f);
        if (lvIt != fieldIoLevel.end()) {
            tmpDec->set_level(lvIt->second);
        }
        /* A coder_arith decoder must use the same file-level prior the encoder
         * did, or the arithmetic decode diverges on the first symbol. */
        if (coderName == "coder_arith") {
            AuxPayloadPtr posPrior =
                (pbgzEngine != nullptr) ? pbgzEngine->getPosPrior() : AuxPayloadPtr();
            if (posPrior && !posPrior->empty()) {
                static_cast<coder_arith*>(tmpDec.get())
                    ->set_prior(posPrior->data(), (uint32_t)posPrior->size());
            }
        }

        std::vector<std::string>& fieldCache = tlenPreDecodedFields[f];
        fieldCache.clear();
        fieldCache.reserve(samLine);

        /* Delta chain of pos_delta mode: the absolute POS reconstructed from the previous line. */
        int64_t posPrev = 0;
        for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
            /*
             * coder_affix_match points last at the caller's output buffer
             * unless need2hold is set; it must be set here, otherwise the
             * cross-line context would be overwritten by the next decode.
             * bwt_cm buffers internally on its own.
             */
            bool need2hold = (coderName == "coder_affix_match");
            /*
             * Field form: CIGAR is always textual (encoded line-wise by
             * compressCigar); POS/PNEXT default to fixed-width binary and are
             * decoded as text only when mode is textual. POS uses a zigzag
             * varint stream: the value is read byte-by-byte (fixed-length
             * decode_line) until the continuation bit is clear. Textual data
             * must never be decoded as fixed-length, otherwise bwt_cm's
             * fixed-length branch would spin past the block boundary.
             */
            if (f == 3 && mode == "pos_delta") {
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
                    vlen = tmpDec->decode_line(&b, 1, UINT8_MAX, need2hold);
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

            bool binary = (f != 5) && (mode == "number" || mode == "");
            int32_t len;
            if (binary) {
                len = tmpDec->decode_line(tmpBuf, 4, UINT8_MAX, need2hold);
            } else {
                len = tmpDec->decode_line(tmpBuf, (uint32_t)sizeof(tmpBuf), '\t', need2hold);
            }
            if (len <= 1) {
                fieldCache.emplace_back();
                continue;
            }

            switch (f) {
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
    }

    /*
     * PNEXT (field 7) in pnext_qname_rebuild mode stores only exception
     * (line, delta) pairs in streams[7]["streams"][0]. Non-exception lines are
     * rebuilt here by pairing records that share the same QNAME and are
     * mutually mapped. mappedPos is already populated by the POS (field 3)
     * pre-decode above. The result fills tlenPreDecodedFields[7] and
     * nextMappedPos for every line.
     */
    std::string pnextMode = streams.isValidIndex(7) && streams[7].isMember("mode")
        ? streams[7]["mode"].asString() : "";
    if (pnextMode == "pnext_qname_rebuild") {
        pnextCache.clear();
        // Decode the exception stream (pairs of contentIdx, delta).
        if (streams[7].isMember("streams") && streams[7]["streams"].size() > 0) {
            Json::Value& excMeta = streams[7]["streams"][0];
            uint32_t excOff = fieldIoStart[7];
            uint32_t excDstLen = excMeta["dstlen"].asUInt();
            std::string excCoderName = excMeta["coder"]["magic"].asString();
            std::shared_ptr<coder_io> excIo = makeCoderIo(buffer + excOff, excDstLen, "PNEXT exc predecode");
            std::shared_ptr<coder> excDec;
            if (excCoderName == "coder_bwt_cm") {
                excDec = std::make_shared<coder_bwt_cm>(excIo.get());
            } else if (excCoderName == "coder_affix_match") {
                excDec = std::make_shared<coder_affix_match>(excIo.get());
            }
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
                size_t excOff = 0;
                auto nextVarint = [&raw, &excOff](uint64_t& out) -> bool {
                    out = 0;
                    int shift = 0;
                    while (excOff < raw.size()) {
                        const uint8_t b = raw[excOff++];
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
    }

    /* Full mate index: built at once once all (pos, pnext) pairs are known. */
    tlenMateIndex.clear();
    for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
        auto posIt = mappedPos.find(lineNo);
        auto pnextIt = nextMappedPos.find(lineNo);
        if (posIt != mappedPos.end() && pnextIt != nextMappedPos.end()) {
            tlenMateIndex[std::make_pair(posIt->second, pnextIt->second)] = lineNo;
        }
    }
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
        Json::Value& splitMeta = idStreams[splitIdx];
        uint32_t splitDstLen = splitMeta["dstlen"].asUInt();
        std::string coderName = splitMeta["coder"]["magic"].asString();

        const bool numericMode = splitMeta.isMember("mode") &&
                                 splitMeta["mode"].asString() == "numeric";
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

            const uint8_t sep = (splitIdx < idSplitSymbols.size())
                                ? (uint8_t)idSplitSymbols[splitIdx] : (uint8_t)'\t';
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
            continue;
        }

        // Decode segment (legacy textual layout)
        int32_t segmentLen = idDecoders[splitIdx]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(),
            (splitIdx < idSplitSymbols.size()) ? idSplitSymbols[splitIdx] : UINT8_MAX, false);
        if (segmentLen < 0) {
            LOG_ERROR("Decode id field segment(%u) failed: %d", splitIdx, segmentLen);
            return -1;
        }
        readOffset += splitDstLen;
        idLength += splitDstLen;
        outputBlock->setDataLen(outputBlock->getDataLen() + segmentLen);
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
            if (fieldMeta["coder"]["magic"].asString() == "coder_fc") {
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
                for (int32_t o = 0; o < decoderLen; ++o) {
                    outputBlock->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
                }
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
                    for (int32_t o = 0; o < decoderLen; ++o) {
                        outputBlock->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
                    }
                } else {
                    decoderLen = readMatchLine(fieldIdx, baseDiffSquashBuffer, actualBaseLen, lineNo);
                    if ((uint32_t)decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d,expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    uint8_t* out = outputBlock->getCurrent();
                    if (lineNo >= cigarOpList.size()) {
                        // No CIGAR op list available: treat as direct 2-bit codes.
                        for (uint32_t o = 0; o < actualBaseLen; ++o) {
                            out[o] = atcg4[baseDiffSquashBuffer[o]];
                        }
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
    if (qualFcv2Decoder != nullptr) {
        /* The strand direction comes with the stream; the base context uses this
           record's already-decoded SEQ (basePtr), corresponding to seqStart on
           the compression side. */
        if (qualFcv2Decoder->decode_record(dst, actualBaseLen, basePtr, actualBaseLen) < 0) {
            LOG_ERROR("Decode quality by fcv2 failed, len = %u", actualBaseLen);
            return -1;
        }
    } else if (qualCmDecoder != nullptr) {
        /*
         * bwt_cm likewise needs no SEQ context. The compression side fed
         * records one by one via encode_line; here they are fetched one by one
         * with the same length, without a delimiter; the length is given by the
         * caller.
         */
        if (qualCmDecoder->decode_line(dst, actualBaseLen) < 0) {
            LOG_ERROR("Decode quality by bwt_cm failed, len = %u", actualBaseLen);
            return -1;
        }
    } else {
        qualCoder->decode_qual_gen2(basePtr, dst, actualBaseLen);
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
    ops.clear();
    if (cigarString == nullptr || cigarLength == 0) {
        return;
    }
    uint32_t currentNumber = 0;
    for (uint32_t i = 0; i < cigarLength; ++i) {
        char ch = cigarString[i];
        if (ch >= '0' && ch <= '9') {
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            if (currentNumber > 0 && (ch == 'M' || ch == 'I' || ch == 'D' || ch == 'N' ||
                                      ch == 'S' || ch == 'H' || ch == 'P' || ch == '=' || ch == 'X')) {
                ops.push_back(CigarOp{ch, currentNumber});
            }
            currentNumber = 0;
        }
    }
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