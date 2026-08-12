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
#include "coder/coder_io.h"
#include "coder/coder_fc.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_qual.h"
#include "coder/coder_fcv2.h"
#include "utils/md5_util.h"
#include "coder/coder_json.h"
#include "log/logger.h"
#include "sam_info.h"
#include "actg.h"
#include "pbgz_manager.h"
#include "utils/path_util.h"
#include "pbgz_index.h"
#include "decompress_engine.h"
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
    baseNPosBuffer = nullptr;
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
}

SamCodecActuator::~SamCodecActuator() {
    MemoryUtil::safeFree(baseNPosBuffer);
    MemoryUtil::safeFree(baseLengthBuffer);
    MemoryUtil::safeFree(baseSquashBuffer);
    MemoryUtil::safeFree(baseDiffSquashBuffer);
    MemoryUtil::safeFree(refeStrecchBuffer);

    // Release idDecoders
    idDecoders.clear();

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
                                fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fastq is %s. \n", refFileName.c_str(), inputFastq.c_str());
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
                                        fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fastq is %s \n", refGeneName.c_str(), inputFastq.c_str());
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
                            LOG_ERROR("Not a valid sam data, baselen = %u, qualityLen = %u ", baseLen, qualityLen);
                            return -1;
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

    if (headEndLine > 0) {
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
        calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
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
        switch (fieldIdx) {
            case 0: // QNAME as FQ:ID
                // ID field: compress based on analysis result
                if (idPosLength == UINT32_MAX) {
                    LOG_DEBUG("Id will compress in all.");
                    fieldDstLen = compressIdFieldInAll(fieldSrcLen, fieldMeta);
                } else {
                    LOG_DEBUG("Id will compress in split.");
                    fieldDstLen = compressIdFieldSplit(fieldSrcLen, fieldMeta);
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
                 * POS 按与上一行 POS 的差值压缩。实测差分约 0.34 B/行，优于定宽二进制
                 * （1.00 B/行）和文本 affix（0.53 B/行）。固定走差分路径，不参与编码器选择。
                 */
                fieldDstLen = compressPosFieldDelta<coder_bwt_cm>(fieldIdx, fieldSrcLen, fieldMeta);
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
                 * PNEXT 按与 POS 的差值压缩。连续行的差值远小于原始坐标，实测
                 * bwt_cm 压差值文本约 0.86 bytes/line，远好于压原始 PNEXT 或 affix。
                 * 试压选择基于的是原始 PNEXT 文本，和这里实际采用的差值编码不一致，
                 * 所以该字段固定走差值路径，不参与编码器选择。
                 */
                fieldDstLen = compressPNextFieldDelta<coder_bwt_cm>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 8: // TLEN
                /*
                 * TLEN 不存原值：解压时按 POS/PNEXT/CIGAR 推算，只对推算不上的行
                 * 存异常。CIGAR/FLAG/RNAME/RNEXT/POS/PNEXT 都在本字段之前压缩，
                 * 推算所需的跟踪 map 已就绪。
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
            case 11: // Optional fields
                 fieldDstLen = compressRegularField(fieldIdx, fieldSrcLen, fieldMeta);
                break;
        }

        // Record statistics for this field
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
    samMeta["streams"] = streamMeta;
    meta["sam"] = samMeta;

    // Calculate MD5 of original data
    std::string md5;
    calcMd5sum(md5, inBlockPtr->getBuffer(), inBlockPtr->getDataLen());
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

template<typename CoderType>
int32_t SamCodecActuator::compressPNextFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "PNEXT delta");
    std::shared_ptr<CoderType> fieldCoder = std::make_shared<CoderType>(fieldIo.get());

    fieldSrcLen = 0;
    uint32_t deltaLength = 0;

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

        uint8_t buff[32] = {0};
        int length = 0;
        if (fieldLength > 1) {
            std::string pNextStr = std::string((char*)fieldStart, fieldLength - 1);
            int64_t pNextValue = (int64_t)std::stoll(pNextStr);
            nextMappedPos[lineIdx] = pNextValue;

            int64_t pos = mappedPos.find(lineIdx) == mappedPos.end() ? 0 : mappedPos[lineIdx];
            int64_t pNextDelta = pNextValue - pos;
            length = snprintf((char*)buff, sizeof(buff), "%" PRId64 "\t", pNextDelta);
            fieldSrcLen += fieldLength;
        } else {
            /* 空/异常值按差值为 0 处理。 */
            length = snprintf((char*)buff, sizeof(buff), "0\t");
            fieldSrcLen += 2;
        }
        fieldCoder->encode_line(buff, (uint32_t)length);
        deltaLength += (uint32_t)length;
    }

    fieldCoder->encode_flush();
    if (fieldIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode PNEXT delta overflow: output buffer too small");
        return -1;
    }
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    fieldMeta["srclen"] = deltaLength;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;
    fieldMeta["mode"] = "pnext_delta";

    LOG_INFO("SAM field(%d) (PNEXT) delta compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

template<typename CoderType>
int32_t SamCodecActuator::compressPosFieldDelta(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "POS delta");
    std::shared_ptr<CoderType> fieldCoder = std::make_shared<CoderType>(fieldIo.get());

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
        uint32_t prevTabPos = contentPos[contentIdx][fieldIdx - 1];
        uint32_t currTabPos = (fieldIdx < contentPos[contentIdx].size()) ? contentPos[contentIdx][fieldIdx] : lineEnd;
        uint8_t* fieldStart = line + prevTabPos + 1;
        uint32_t fieldLength = currTabPos - prevTabPos;

        uint8_t buff[32] = {0};
        int length = 0;
        if (fieldLength > 1) {
            std::string posStr = std::string((char*)fieldStart, fieldLength - 1);
            int64_t posValue = (int64_t)std::stoll(posStr);
            mappedPos[lineIdx] = posValue;
            int64_t delta = posValue - prevPos;
            length = snprintf((char*)buff, sizeof(buff), "%" PRId64 "\t", delta);
            fieldSrcLen += fieldLength;
            prevPos = posValue;
        } else {
            /* 空/异常值按差值为 0 处理。 */
            mappedPos[lineIdx] = 0;
            length = snprintf((char*)buff, sizeof(buff), "0\t");
            fieldSrcLen += 2;
        }
        fieldCoder->encode_line(buff, (uint32_t)length);
        deltaLength += (uint32_t)length;
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

int32_t SamCodecActuator::compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    // Similar to FastqActuator::compressIdInSplit, create multiple streams for each split symbol
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    // Process each split symbol (similar to FastqActuator)
    for (uint32_t i = 0; i < idSplitSymbols.size(); ++i) {
        std::shared_ptr<coder_io> idIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QNAME sub-stream");
        std::shared_ptr<coder_affix_match> idCoder = std::make_shared<coder_affix_match>(idIo.get());

        /*
         * 纯数字段检测（取首行）：排除行尾分隔符，否则像 "265269\t" 这样的纯数字列
         * 永远判不出来。检测结果只用于日志——实测 affix 压这类短列段比 bwt_cm 文本
         * 和逐行 delta 都更好，所以所有段统一用 affix，不做分支。
         */
        std::vector<size_t>& npos = inBlockPtr->getNpos();
        uint32_t lineNum = npos.size();
        uint8_t* buffer = inBlockPtr->getBuffer();
        bool idDigit = false;
        if (idSplitMaxLen[i] != idSplitMinLen[i]) {
            uint32_t lineStart = (headEndLine == 0) ? 0 : npos[headEndLine - 1] + 1;
            uint8_t* line = buffer + lineStart;
            uint8_t* segPtr = nullptr;
            uint32_t segLen = 0;
            if (i == 0) {
                segPtr = line;
                segLen = idSplitPos[0][0] + 1;
            } else {
                uint32_t prevPos = idSplitPos[0][i - 1];
                uint32_t currPos = idSplitPos[0][i];
                segPtr = line + prevPos + 1;
                segLen = currPos - prevPos;
            }
            idDigit = (segLen > 1);
            for (uint32_t j = 0; idDigit && j + 1 < segLen; ++j) {
                if ((segPtr[j] & 0xF0) != 0x30) {
                    idDigit = false;
                    break;
                }
            }
        }
        LOG_DEBUG("SAM QNAME split column(%u) uses coder_affix_match (idDigit=%d)", i, (int)idDigit);

        uint32_t srcLength = 0;
        // Process each line and compress the specific split segment
        for (uint32_t lineIdx = headEndLine; lineIdx < lineNum; ++lineIdx) {
            uint32_t lineStart = (lineIdx == 0) ? 0 : npos[lineIdx - 1] + 1;
            uint8_t* line = buffer + lineStart;

            // Skip header lines (starting with @)
            if (*line == '@') {
                continue;
            }

            // Extract the specific split segment for this symbol
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

            // Encode the segment data
            if (segmentLength > 0) {
                idCoder->encode_line(segmentStart, segmentLength);
                srcLength += segmentLength;
            }
        }

        // Flush the encoder for this segment
        idCoder->encode_flush();
        if (idIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode id segment overflow: output buffer too small");
            return -1;
        }

        // Update output block data length
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + idIo->data_len);

        // Create metadata for this stream (similar to FastqActuator)
        Json::Value tmpMeta;
        tmpMeta["srclen"] = srcLength;
        tmpMeta["dstlen"] = idIo->data_len;
        tmpMeta["coder"] = idIo->meta;
        tmpMeta["splitidx"] = i; // Index of split symbol

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
        // Encode the chromosome index
        chrCoder->encode_line(reinterpret_cast<const uint8_t*>(&chrIndex), sizeof(chrIndex));
        srcLen += sizeof(chrIndex);
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
             * FLAG/POS/PNEXT 被以文本形态（affix）压缩时，这里的跟踪 map 必须照常
             * 填好：compressBaseWithRef 在 SEQ 阶段要靠 mappedPos/mappedChr/mappedFlag
             * 还原参考位置，解压侧也要靠它们还原。二进制形态在 compressNumber 里填，
             * 文本形态在这里填，两边保持一致。
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
        } else {
            uint32_t sequeceLength = parseCigar(fieldStart, fieldLength);
            baseLengthBuffer[lineIdx - headEndLine] = sequeceLength;
            /* 参考跨度供 TLEN 推算（field 8 在本字段之后压缩）。 */
            cigarReadLen[lineIdx] = parseCigarRefConsumed(fieldStart, fieldLength);
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
    /* 只统计消耗参考序列的 CIGAR 操作：M/D/N/=/X（含小写）。 */
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
 * 按 SAM 规范推算 TLEN：
 *   |TLEN| = 最右 mapped base - 最左 mapped base + 1
 * 左侧片段为正、右侧片段为负。pos < pnext 时本读在左，右端 = pnext + 伙伴参考跨度 - 1；
 * pos > pnext 时本读在右，左端 = pos，模板长为负。
 */
int32_t SamCodecActuator::computeTLEN(uint32_t lineIdx) {
    /* 非配对（FLAG bit0）或任一端未比对（bit2/bit3），TLEN 定为 0。 */
    auto flagIt = mappedFlag.find(lineIdx);
    if (flagIt == mappedFlag.end() || !(flagIt->second & 0x1) ||
        (flagIt->second & 0x4) || (flagIt->second & 0x8)) {
        return 0;
    }

    /* 参照序列不可用或伙伴在不同参照上，TLEN 定为 0。 */
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

    /* 本读的参考跨度来自 CIGAR。 */
    auto readLenIt = cigarReadLen.find(lineIdx);
    uint32_t refSpan = (readLenIt != cigarReadLen.end()) ? readLenIt->second : 0;

    /* 伙伴跨度通过 (pnext, pos) 反查完整索引；找不到按 0 处理（跨度缺失只损失压缩率）。 */
    uint32_t mateRefSpan = 0;
    auto mateIt = tlenMateIndex.find(std::make_pair(pnext, pos));
    if (mateIt != tlenMateIndex.end()) {
        auto spanIt = cigarReadLen.find(mateIt->second);
        if (spanIt != cigarReadLen.end()) {
            mateRefSpan = spanIt->second;
        }
    }

    int64_t templateLen;
    if (pos < pnext) {
        templateLen = pnext + (int64_t)mateRefSpan - pos;
    } else if (pos > pnext) {
        templateLen = -(pos + (int64_t)refSpan - pnext);
    } else {
        templateLen = pnext + (int64_t)mateRefSpan - pos;
    }
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

    /* 推算不上的行存 (相对行号, 实际值)。 */
    std::vector<std::pair<uint32_t, int32_t>> tlenExceptions;

    /* POS/PNEXT 已在前面字段压好，这里建完整伙伴索引。 */
    tlenMateIndex.clear();
    for (const auto& entry : nextMappedPos) {
        auto posIt2 = mappedPos.find(entry.first);
        if (posIt2 == mappedPos.end()) {
            continue;
        }
        tlenMateIndex[std::make_pair(posIt2->second, entry.second)] = entry.first;
    }

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
            int32_t computedTLEN = computeTLEN(lineIdx);
            if (computedTLEN != currentTLEN) {
                tlenExceptions.push_back(std::make_pair(contentIdx, currentTLEN));
            }
        }
    }

    if (!tlenExceptions.empty()) {
        std::shared_ptr<coder_io> tlenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "TLEN exceptions");
        std::shared_ptr<CoderType> tlenCoder = std::make_shared<CoderType>(tlenIo.get());

        uint32_t tlenExcSrcLen = (uint32_t)(tlenExceptions.size() * (2 * sizeof(int32_t)));
        int32_t* tlenExcBuffer = MemoryUtil::safeAlloc<int32_t>(tlenExcSrcLen);
        if (tlenExcBuffer == nullptr) {
            return -1;
        }
        for (uint32_t i = 0; i < tlenExceptions.size(); ++i) {
            tlenExcBuffer[2 * i] = (int32_t)tlenExceptions[i].first;
            tlenExcBuffer[(2 * i) + 1] = tlenExceptions[i].second;
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
    std::unique_ptr<uint8_t[]> tmpBuffer = std::make_unique<uint8_t[]>(outBlockPtr->getBlockSize());
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
    std::shared_ptr<coder_io> fieldIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ");
    std::shared_ptr<coder> fieldCoder = makeFieldEncoder(fieldIdx, samFieldDefaultCoder(fieldIdx, CoderType::FC), fieldIo.get(), false);

    fieldCoder->encode_line(tmpBuffer.get(), fieldSrcLen);
    fieldCoder->encode_flush();
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

    uint32_t offset = 0;
    uint64_t nOffset = 0;

    // Initialize mapping buffers similar to FastqActuator
    const uint32_t baseMaxLength = inBlockPtr->getMaxLineLen() + 4;
    const uint32_t lsquash = (baseMaxLength >> 2) + !!(baseMaxLength & 0x3);

    uint32_t baseMappedLength = (baseMaxLength << 1);

    std::unique_ptr<uint8_t[]> basePairBuffer = std::make_unique<uint8_t[]>(baseMaxLength);
    baseSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(lsquash);
    std::unique_ptr<uint8_t[]> baseMappedBuffer = std::make_unique<uint8_t[]>(baseMappedLength);
    baseNPosBuffer = MemoryUtil::safeAlloc<uint32_t>(baseNCount);

    // Create metadata structure
    Json::Value metaSubs;
    Json::Value metaStreams;
    uint32_t totalSrcLen = 0;
    uint32_t totalDstLen = 0;

    // Second pass: compress with reference
    std::shared_ptr<coder_io> matchIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ match");
    std::shared_ptr<coder_bwt_cm> matchCm = std::make_shared<coder_bwt_cm>(matchIo.get());
    CoderFactory::applyLevel(matchIo.get(), CoderType::BWT_CM, engineCompressLevel());
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

        // Process sequence: remove N's and record positions
        uint32_t outLen = 0;
        for (uint32_t n = 0; n < seqLength; n++) {
            char ch = seqStart[n];
            if (ch == 'N' || ch == 'n') {
                baseNPosBuffer[nOffset] = totalBaseLength + n;
                nOffset++;
            }
        }

        totalBaseLength += seqLength;
        if (chrId != 0xFFFF && chrId != 0xFFFE && !(flag & 0x04)) {
            // Get chromosome start position from SamInfo
            int64_t chrStartPos =  SamInfo::getInstance().getPositionByIndex(chrId);
            if (chrStartPos == -1) {
                // No valid mapping, encode directly
                // LOG_DEBUG("chrStartPos not found, line = %d", contentIdx);
                actgEncode(seqStart, baseMappedBuffer.get(), seqLength);
                outLen = seqLength;
            } else {
                // Calculate actual reference position
                int64_t refPos = chrStartPos + startPos - 1; // SAM is 1-based
                // Determine strand direction from FLAG bit 4
                uint32_t squashBufferLength = actgSquash(seqStart, seqLength, baseSquashBuffer);
                uint8_t shiftBitLength = refPos % 4;
                int64_t refSquashPos = refPos / 4;
                uint32_t baseSquashLength = (seqLength >> 2) + !!(seqLength & 0x3) + 1;
                if (refSquashPos + baseSquashLength > pRefeGene->getSquashLength()) {
                    // LOG_DEBUG("Mapped pos is out of bound. line = %d", contentIdx);
                    actgEncode(seqStart, baseMappedBuffer.get(), seqLength);
                    outLen = seqLength;
                } else {
                    const uint8_t* beginRefPos = pRefeGene->getSquash() + refSquashPos;
                    uint8_t* refeMappedPos = MemoryUtil::safeAlloc<uint8_t>(squashBufferLength);;
                    if (shiftBitLength == 0) {
                        memcpy(refeMappedPos, beginRefPos, squashBufferLength);
                    } else if (shiftBitLength == 1) {
                        for (uint32_t i = 0; i < squashBufferLength; ++i) {
                            refeMappedPos[i] = ((beginRefPos[i] << 2) & 0xFC) + ((beginRefPos[i + 1] >> 6) & 0x03);
                        }
                    } else if (shiftBitLength == 2) {
                        for (uint32_t i = 0; i < squashBufferLength; ++i) {
                            refeMappedPos[i] = ((beginRefPos[i] << 4) & 0xF0) + ((beginRefPos[i + 1] >> 4) & 0x0F);
                        }
                    } else if (shiftBitLength == 3) {
                        for (uint32_t i = 0; i < squashBufferLength; ++i) {
                            refeMappedPos[i] = ((beginRefPos[i] << 6) & 0xC0) + ((beginRefPos[i + 1] >> 2) & 0x3F);
                        }
                    }
                    outLen = actgStretchMappingXor(baseSquashBuffer, refeMappedPos, squashBufferLength, baseMappedBuffer.get());
                    MemoryUtil::safeFree(refeMappedPos);
                    pRefeGene->updateMatchedGene(refPos, baseSquashLength << 2);
                }
            }
        } else {
            // No valid mapping, encode directly
            // LOG_DEBUG("Not mapping, line = %d", contentIdx);
            actgEncode(seqStart, baseMappedBuffer.get(), seqLength);
            outLen = seqLength;
        }
        // Encode the mapped data
        if (outLen > 0) {
            matchCm->encode_line(baseMappedBuffer.get(), seqLength);
            srcLen += seqLength;
        }
        offset++;
    }
    // Flush match encoder
    matchCm->encode_flush();
    if (matchIo->err != coder_io::IO_OK) {
        LOG_ERROR("Encode base match stream overflow: output buffer too small");
        return -1;
    }
    // First sub-stream: match stream between reads and reference
    metaSubs.clear();
    metaSubs["srclen"] = srcLen;
    metaSubs["dstlen"] = matchIo->data_len;
    metaSubs["coder"] = matchIo->meta;
    metaSubs["sname"] = "m";
    metaStreams.append(metaSubs);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + matchIo->data_len);
    totalDstLen += matchIo->data_len;
    totalSrcLen += srcLen;

    // Second sub-stream: positions of N's
    if (baseNCount > 0) {
        std::shared_ptr<coder_io> nposIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ npos");
        uint32_t nposSrcLen = (baseNCount << 2); // 4 bytes per N position
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(nposIo.get());
        CoderFactory::applyLevel(nposIo.get(), CoderType::BWT_CM, engineCompressLevel());
        subCoder->encode_line((uint8_t*)baseNPosBuffer, nposSrcLen);
        subCoder->encode_flush();
        if (nposIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode base N positions overflow: output buffer too small");
            return -1;
        }
        metaSubs.clear();
        metaSubs["srclen"] = nposSrcLen;
        metaSubs["dstlen"] = nposIo->data_len;
        metaSubs["coder"] = nposIo->meta;
        metaSubs["sname"] = "npos";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + nposIo->data_len);
        totalSrcLen += nposSrcLen;
        totalDstLen += nposIo->data_len;
    }

    if (minBaseLength != maxBaseLength) {
        std::shared_ptr<coder_io> lenIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "SEQ length");
        uint32_t baseLenSrcLen = unmapedReadLength.size() << 1;
        uint32_t* baseLenBuffer = MemoryUtil::safeAlloc<uint32_t>(baseLenSrcLen);
        if (baseLenBuffer == nullptr) {
            return -1;
        }
        for (uint32_t i = 0; i < unmapedReadLength.size(); i++) {
            std::pair baseLenPos = unmapedReadLength[i];
            baseLenBuffer[2 * i] = baseLenPos.first;
            baseLenBuffer[(2 * i) + 1] = baseLenPos.second;
        }
        std::shared_ptr<coder_bwt_cm> lenCoder = std::make_shared<coder_bwt_cm>(lenIo.get());
        CoderFactory::applyLevel(lenIo.get(), CoderType::BWT_CM, engineCompressLevel());
        lenCoder->encode_line((uint8_t*)baseLenBuffer, baseLenSrcLen<<2);
        lenCoder->encode_flush();
        MemoryUtil::safeFree(baseLenBuffer);
        if (lenIo->err != coder_io::IO_OK) {
            LOG_ERROR("Encode base length stream overflow: output buffer too small");
            return -1;
        }

        metaSubs.clear();
        metaSubs["srclen"] = baseLenSrcLen<<2;
        metaSubs["dstlen"] = lenIo->data_len;
        metaSubs["coder"] = lenIo->meta;
        metaSubs["sname"] = "baselen";
        metaStreams.append(metaSubs);
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + lenIo->data_len);
        totalSrcLen += baseLenSrcLen<<2;
        totalDstLen += lenIo->data_len;
    }

    // Set metadata
    fieldMeta["ncount"] = baseNCount;
    fieldMeta["minlen"] = minBaseLength;
    fieldMeta["maxlen"] = maxBaseLength;
    fieldMeta["totalsrclen"] = totalSrcLen;
    fieldMeta["totaldstlen"] = totalDstLen;
    fieldMeta["streams"] = metaStreams;
    fieldMeta["field"] = fieldIdx;
    /*
     * 只修正 fieldSrcLen：SEQ 列原始大小 = 匹配流源长（srcLen，即各记录 SEQ 文本
     * 总长，与 QUAL 列口径一致）；npos/baselen 辅助流不计入。meta 的 totalsrclen
     * 保持原样（含辅助流）。
     */
    fieldSrcLen = (uint32_t)srcLen;

    LOG_INFO("SAM base field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%",
        totalSrcLen, totalDstLen, (double)(totalDstLen * 100)/(double)totalSrcLen);

    return totalDstLen;
}

int32_t SamCodecActuator::compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<size_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create quality encoder similar to FastqActuator
    std::shared_ptr<coder_io> qualityIo = makeCoderIo(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), "QUAL");

    /*
     * 质量值可以用两种编码器。coder_qual 是原有的，以 SEQ 为上下文；fcv2 是上下文
     * 混合编码器，以前一个和前两个质量值、记录内的测序循环序号、链方向为上下文。
     *
     * 用哪一个由预处理阶段的试压结果决定。质量值列走的是专用评估路径
     * （QualSelector），候选就是这两个，样本按记录采集因而保留了记录边界和链方向。
     *
     * 预处理没跑成、样本太小被跳过、或者当前引擎不提供预处理信息时，coderFor 返回
     * 传入的兜底类型 QUAL，也就是沿用原有的 coder_qual，行为与接通选择之前一致。
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
        std::vector<uint32_t> freqByByte(256, 0);
        for (uint32_t i = 0; i < qualFreqTable.size(); ++i) {
            uint32_t b = (uint32_t)qualFreqTable[i].first + (uint32_t)'!';
            if (b < 256) {
                /* qualFreqTable 只保留了按频率降序的字母表，实际计数没有留下来，
                   这里用排名折算出一个单调递减的权重，仅用于确定哈夫曼树形状。 */
                freqByByte[b] = (uint32_t)(qualFreqTable.size() - i);
            }
        }
        /*
         * 有先验就从先验起步：模型不必从固定初值重新学，块越小收益越明显。
         * 先验的绝对地址要一并写进块 meta——解码方只有拿到同一份快照才能对上，
         * 而它在随机读场景下无法沿着顺序流推断出这个地址。
         *
         * 加载失败在压缩侧是良性回退（保留固定初值模型，仅损失压缩率），
         * 所以只有 loaded 为真时才登记地址，避免给解码方一个它无法兑现的承诺。
         */
        qualPriorAddress = (pbgzEngine != nullptr) ? pbgzEngine->getQualPriorAddress() : -1;
        AuxPayloadPtr priorBlob =
            (pbgzEngine != nullptr) ? pbgzEngine->getQualPrior(0) : AuxPayloadPtr();
        bool priorLoaded = false;
        if (priorBlob && !priorBlob->empty()) {
            fcv2Coder = std::make_shared<coder_fcv2>(qualityIo.get(), freqByByte,
                                                     *priorBlob, &priorLoaded);
        } else {
            fcv2Coder = std::make_shared<coder_fcv2>(qualityIo.get(), freqByByte);
        }
        if (!priorLoaded) {
            qualPriorAddress = -1;
        }
    } else if (useBwtCm) {
        /*
         * bwt_cm 不需要字母表，也不需要 SEQ 或链方向，逐条喂质量值即可。
         * 它是 fcv2 不适用时的第二顺位：实测比 coder_qual 好 7.4 个百分点。
         */
        qualCmCoder = std::make_shared<coder_bwt_cm>(qualityIo.get());
        CoderFactory::applyLevel(qualityIo.get(), CoderType::BWT_CM, engineCompressLevel());
    } else {
        qualityCoder = std::make_shared<coder_qual>(qualityIo.get(), true, qualFreqTable);
    }

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
        if (fieldIdx >= 1 && contentPos[contentIdx].size() >= fieldIdx) {
            seqStart = line + contentPos[contentIdx][fieldIdx - 2] + 1;
        }

        if (useFcv2) {
            /*
             * 从 FLAG 字段取链方向。按 SAM 规范，0x10 置位表示 SEQ 和 QUAL 在文件里
             * 是相对参考正链存放的，也就是相对测序仪读出顺序已经反转，编码器需要据此
             * 还原真实的循环序号。这里手工解析而不是用 strtol，是因为字段没有以 \0
             * 结尾，且这是每条记录都会走到的热路径。
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
            fcv2Coder->encode_record(qualStart, qualLength, rev);
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
        /* 先验块容器头的绝对文件偏移，解码方凭它取回同一份快照。 */
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

    /* 频率表辅助流仍计入 meta 的 totalsrclen（保持原有口径）。 */
    totalSrcLength += freqSrcLen;
    totalDstLength += qualityFreqIo->data_len;

    // Set field metadata
    fieldMeta["totalsrclen"] = totalSrcLength;
    fieldMeta["totaldstlen"] = totalDstLength;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["field"] = fieldIdx;

    /*
     * 只修正 fieldSrcLen：QUAL 列原始大小 = 质量值文本总长（streamSrcLen），
     * 频率表辅助流不计入；meta 的 totalsrclen 保持原样（含 freqSrcLen）。
     */
    fieldSrcLen = streamSrcLen;

    /*
     * verbose 模式下逐块打印本块 QUAL 到底走了哪一路编码器、以及本块自己的压缩率。
     *
     * 为什么必须在这里做：QUAL 的选择是"预处理阶段一次决定、后续所有块沿用"的模式，
     * 但一次压缩里有 ~10 个线程并发压不同的块，最终日志只有一行汇总，光看那一行没
     * 办法确认每个块真的都走了预期的编码器，也没法看到块与块之间压缩率的抖动。
     *
     * 三点约束：
     *   1. 只在 verbose=true 时输出，非 verbose 时保持默认输出字节级一致；
     *   2. 必须整行一次 fprintf 打完，因为多线程并发写 stderr，拆成多次调用会互相
     *      穿插，肉眼没法读；单次 fprintf 写入内核缓冲区通常是原子的；
     *   3. totalSrcLength 为 0 时直接把比例置 0，避免除零。
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
    // Parse meta information
    initMetaInfo();

    /*
     * 块入口预分配：按文件头的 block_size（压缩时确定的上界）×2 保证输出缓冲够——
     * 确定值，不是估算，不受 fieldcount/读长影响。堵所有字段越界的主防线。
     * block_size 从 baseFileMeta 读回（DecompressEngine::createBlockReader），此刻
     * 已就绪；为 0（老文件没写）时落回默认 getBlockSize()。
     * coder_io 的 putc 检查与 decode 错误返回链作兜底（见 decompressQuality 等）。
     */
    size_t bs = pbgzEngine->getFileBlockSize();
    if (bs == 0) {
        bs = ConfigManager::getInstance().getBlockSizeByCompressLevel(pbgzEngine->getParameter().compressLevel);
    }
    if (outBlockPtr->ensureCapacity(bs * 2) != 0) {
        LOG_ERROR("preallocate output buffer failed, block_size=%zu", bs);
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

    // ensureCapacity 已在 decompress() 开头调用，此处不重复

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
    refeStrecchBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    uint32_t totalBaseLen = 0;
    uint32_t nposOffset = 0;

    uint8_t* pBaseOut = nullptr;
    if (streams[9]["coder"]["magic"].asString() == "coder_fc") {
        /*
         * coder_fc 是"整块"编码器：SEQ 必须一次性全解出来，但最终 SAM 输出是按行
         * 交错的（ID\tFLAG\t...\tSEQ\tQUAL\n），所以 SEQ 只能先落到别处、再逐行搬。
         * 这里把它暂存在 outputBlock 缓冲的**尾部**（见 initDecoder 同样的落点），
         * 头部正常向前追加行内容，两者对向生长、互不覆盖：
         *   头部输出 <= block_size，尾部 SEQ <= block_size，块入口已按 block_size*2
         *   一次性 ensureCapacity，容量有确定上界。
         * 不用独立 malloc 缓冲的原因：memcpy 次数完全一样（都得逐行搬），独立缓冲只
         * 会多出每块一次 malloc/free、首次写的缺页开销，以及"线程数 × 整条 SEQ"的
         * 额外内存峰值。
         *
         * 不变量：块内绝对不能对 outputBlock 做 realloc，否则这里的尾部指针以及
         * 下面 basePtr 都会悬空。
         */
        pBaseOut = pBaseEnd - streams[9]["totalsrclen"].asUInt();
    }

    /*
     * 先把 POS/CIGAR/PNEXT 预解码出来：TLEN 推算需要完整伙伴索引，逐行解码时伙伴
     * 可能在块的后半段还没解到。预解码结果按行缓存，主循环直接拷贝，避免二次解码。
     */
    posDeltaPrev = 0;
    if (0 != preDecodeForTLEN()) {
        LOG_ERROR("Pre-decode for TLEN failed.");
        return -1;
    }

    /* 拷贝预解码字段的字节，避免再次解码；未缓存时返回 -1。 */
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
                decoderLen = decompressBase(fieldIdx, streams[fieldIdx], pBaseOut, lineNo, nposOffset, totalBaseLen, outputBlock);
                actualBaseLen = decoderLen;
            } else if (fieldIdx == 10 ) {  /// QUAL
                decoderLen = decompressQuality(basePtr, actualBaseLen, outputBlock);
                // No optional fields scenario, replace appended \t with \n
                if (fieldIdx + 1 == fieldCount) {
                    uint8_t* pEnd = outputBlock->getCurrent();
                    *(pEnd - 1) = '\n';
                }
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
            for (uint32_t i = 0; i < idStreamMeta.size(); ++i) {
                std::string coderName = idStreamMeta[i]["coder"]["magic"].asString();
                uint32_t dstLength = idStreamMeta[i]["dstlen"].asUInt();
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "QNAME sub-stream");
                ioVector.push_back(io);
                if (coderName == "coder_affix_match") {
                    idDecoders.push_back(std::make_shared<coder_affix_match>(io.get()));
                    idDecoders.back()->set_level(idStreamMeta[i]["coder"]["level"].asInt());
                } else if (coderName == "coder_bwt_cm") {
                    idDecoders.push_back(std::make_shared<coder_bwt_cm>(io.get()));
                    idDecoders.back()->set_level(idStreamMeta[i]["coder"]["level"].asInt());
                } else {
                    LOG_ERROR("Unsupport coder name:%s", coderName.c_str());
                    return -1;
                }
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
                     * 整块 SEQ 解到 outputBlock 缓冲的尾部作暂存，decompressBase 再逐行
                     * 搬到头部。不另开独立缓冲：那样每块要多一次整条 SEQ 的 malloc/free
                     * 和首次触碰缺页，峰值内存还要多出 线程数×SEQ大小，而拷贝次数一样。
                     *
                     * 不变量：块内绝不能对 outputBlock 做 realloc，否则这个尾部指针悬空。
                     * 容量由块入口一次性按 block_size*2 预分配保证（见 decompress()）——
                     * 头部输出 ≤ block_size、尾部暂存 ≤ block_size，正好 2 倍。
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
                std::shared_ptr<coder_io> io = makeCoderIo(inBlockPtr->getBuffer() + readOffset, dstLength, "SEQ match");
                ioVector.push_back(io);
                if (baseMetaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
                    fieldDecoders[idx] =  std::make_shared<coder_bwt_cm>(io.get());
                } else {
                    LOG_ERROR("check sub stream failed, coder name not match");
                    return -1;
                }

                readOffset += dstLength;
                baseNCount = baseMeta["ncount"].asUInt();
                if (baseNCount > 0) {
                    id++;
                    if (baseMetaStreams[id]["sname"].asString() != "npos") {
                        LOG_ERROR("check sub stream failed. sname not match");
                        return -1;
                    }

                    uint32_t dstlen = baseMetaStreams[id]["dstlen"].asUInt();
                    uint32_t srclen = baseMetaStreams[id]["srclen"].asUInt();
                    baseNPosBuffer = MemoryUtil::safeAlloc<uint32_t>(baseNCount);
                    if (baseMetaStreams[id]["coder"]["magic"].asString() == "coder_bwt_cm") {
                        coder_io nposIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "SEQ npos");
                        auto nposCoder = std::make_unique<coder_bwt_cm>(&nposIo);
                        if (nposCoder->decode_line((uint8_t*)baseNPosBuffer, srclen, UINT8_MAX, false) < 0) {
                            LOG_ERROR("Decode base N positions failed");
                            return -1;
                        }
                    } else {
                        LOG_ERROR("check sub stream failed. coder not match. coder = %s.",
                            baseMetaStreams[id]["coder"]["magic"].asString().c_str());
                        return -1;
                    }
                    readOffset += dstlen;
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
                    uint32_t* baseLenPtr = (uint32_t*)baseLenBuffer;
                    uint32_t baseLenCount = srclen >> 2;
                    for (uint32_t i = 0; i < baseLenCount; i += 2) {
                        unmapedReadLength.push_back(std::make_pair(baseLenPtr[i], baseLenPtr[i+1]));
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
                /* 同 fastq_actuator: 计数用 uint8_t 会在字母表超过 127 个符号时回绕, 导致堆越界。 */
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
                     * fcv2 的字母表和各符号频率都写在它自己的码流头部，begin_decode
                     * 会读回来重建哈夫曼树，所以这里构造时传空频率表即可。
                     */
                    std::shared_ptr<coder_io> qualIo = makeCoderIo(inBlockPtr->getBuffer() + readOffset, qualDstLength, "QUAL");
                    ioVector.push_back(qualIo);
                    std::vector<uint32_t> emptyFreq(256, 0);
                    /*
                     * 编码方若从先验起步，解码方必须用同一份快照起步，否则模型立刻发散。
                     * 与压缩侧的良性回退不同，这里加载失败必须致命：静默回退到固定初值
                     * 会解出一片看似合法的错误数据，比直接失败危险得多。
                     *
                     * 先验的地址来自块 meta 而非顺序推断，随机读因此同样成立。
                     */
                    if (qualStreamMeta[0].isMember("prior")) {
                        /*
                         * meta 里的偏移只是 seek 手段与校验值，索引键用包序号：
                         * 绝对偏移在管道输入下退化为 0，按它查表会全部落空。
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
                    if (qualFcv2Decoder->begin_decode() != 0) {
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
            /* TLEN 字段：优先走推算优化，异常流单独解码。 */
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
                        int32_t* excBuffer = MemoryUtil::safeAlloc<int32_t>(srclen >> 2);
                        if (excBuffer == nullptr) {
                            return -1;
                        }
                        coder_io tlenIo(inBlockPtr->getBuffer() + readOffset, dstlen, &ioErrSink, "TLEN exceptions");
                        if (stream["coder"].isMember("level")) {
                            tlenIo.meta["level"] = stream["coder"]["level"].asInt();
                        }
                        coder_bwt_cm tlenDecoder(&tlenIo);
                        if (tlenDecoder.decode_line((uint8_t*)excBuffer, srclen, UINT8_MAX, false) < 0) {
                            MemoryUtil::safeFree(excBuffer);
                            LOG_ERROR("Decode TLEN exceptions failed");
                            return -1;
                        }
                        uint32_t excCount = srclen / (2 * sizeof(int32_t));
                        for (uint32_t i = 0; i < excCount; ++i) {
                            tlenCache[(uint32_t)excBuffer[2 * i]] = excBuffer[2 * i + 1];
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
        } else {
            std::string coderName = streamMeta[idx]["coder"]["magic"].asString();
            uint32_t dstLen = streamMeta[idx]["dstlen"].asUInt();
            /*
             * 记下 POS(3)/CIGAR(5)/PNEXT(7) 的压缩流位置，preDecodeForTLEN 据此重建
             * 解码器把全块预解码，TLEN 推算才能拿到完整伙伴索引。
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
    /* decode_line 负值是错误码（码流损坏/缓冲不足），原样上抛，不能当长度用。 */
    int32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), outputBlock->getRemain(), splitFlag, false);
    if (fieldLen < 0) {
        LOG_ERROR("Decode regular field(%u) failed: %d", fieldIdx, fieldLen);
        return -1;
    }
    outputBlock->setDataLen(outputBlock->getDataLen() + fieldLen);

    /*
     * FLAG/POS/PNEXT 以文本形态（affix）压缩时，这里要把跟踪 map 填回来：
     * 后续 QUAL 阶段与参考序列阶段都要读 mappedPos/mappedChr/mappedFlag。
     * 与二进制形态 decompressNumber 里填的保持一致。
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

int32_t SamCodecActuator::decompressPNextFieldDelta(uint32_t fieldIdx, uint32_t lineNo, uint8_t splitFlag, RoughIOBlock* outputBlock) {
    uint8_t deltaBuffer[32] = {0};
    int32_t deltaLen = fieldDecoders[fieldIdx]->decode_line(deltaBuffer, sizeof(deltaBuffer), splitFlag, true);
    if (deltaLen < 0) {
        LOG_ERROR("Decode PNEXT delta failed at line %u", lineNo);
        return -1;
    }

    int64_t pos = mappedPos.find(lineNo) == mappedPos.end() ? 0 : mappedPos[lineNo];
    int64_t pNext = pos;
    if ((uint32_t)deltaLen > 1) {
        std::string pNextDeltaStr = std::string((char*)deltaBuffer, (size_t)deltaLen - 1);
        try {
            pNext = (int64_t)std::stoll(pNextDeltaStr) + pos;
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to parse delta value '%s' for line %d: %s", pNextDeltaStr.c_str(), lineNo, e.what());
            return -1;
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
    uint8_t deltaBuffer[32] = {0};
    int32_t deltaLen = fieldDecoders[fieldIdx]->decode_line(deltaBuffer, sizeof(deltaBuffer), splitFlag, true);
    if (deltaLen < 0) {
        LOG_ERROR("Decode POS delta failed at line %u", lineNo);
        return -1;
    }

    int64_t pos = posDeltaPrev;
    if ((uint32_t)deltaLen > 1) {
        std::string deltaStr((char*)deltaBuffer, (size_t)deltaLen - 1);
        try {
            pos += (int64_t)std::stoll(deltaStr);
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to parse POS delta '%s' for line %d: %s", deltaStr.c_str(), lineNo, e.what());
            return -1;
        }
    }
    posDeltaPrev = pos;
    mappedPos[lineNo] = pos;
    char buff[32];
    int posLen = snprintf(buff, sizeof(buff), "%" PRId64 "\t", pos);
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
            val = computeTLEN(lineNo);
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
 * 预解码 POS/CIGAR/PNEXT 并把结果按行缓存。TLEN 推算需要完整伙伴索引与参考跨度，
 * 主循环逐行解码时伙伴可能落在块的后半段，拿不到；这里先解一遍全部行，主循环随后
 * 直接拷贝缓存，不再二次解码。
 */
int32_t SamCodecActuator::preDecodeForTLEN() {
    if (samLine == 0 || !meta.isMember("sam")) {
        return 0;
    }

    uint8_t* buffer = inBlockPtr->getBuffer();
    uint8_t tmpBuf[1024];
    Json::Value& streams = meta["sam"]["streams"];

    /* 推算 TLEN 需要的字段：POS(3)、CIGAR(5)、PNEXT(7)。 */
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
        } else {
            continue;
        }
        auto lvIt = fieldIoLevel.find(f);
        if (lvIt != fieldIoLevel.end()) {
            tmpDec->set_level(lvIt->second);
        }

        std::vector<std::string>& fieldCache = tlenPreDecodedFields[f];
        fieldCache.clear();
        fieldCache.reserve(samLine);

        /* pos_delta 模式的差值链：上一行还原出的绝对 POS。 */
        int64_t posPrev = 0;
        for (uint32_t lineNo = 0; lineNo < samLine; ++lineNo) {
            /*
             * coder_affix_match 把 last 指向调用方输出缓冲，除非 need2hold 置位；
             * 这里必须置位，否则跨行上下文会被下一次解码覆盖。bwt_cm 内部自行缓冲。
             */
            bool need2hold = (coderName == "coder_affix_match");
            /*
             * 字段形态：CIGAR 恒为文本（compressCigar 逐行编码）；POS/PNEXT 缺省按
             * 二进制定宽处理，mode 为文本形态时才按文本解。文本数据绝不能按定长解，
             * 否则 bwt_cm 的定长分支会越过块边界空转。
             */
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
                    if (mode == "pos_delta") {
                        int64_t delta = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                        int64_t pos = posPrev + delta;
                        posPrev = pos;
                        mappedPos[lineNo] = pos;
                        fieldCache.emplace_back(std::to_string(pos) + '\t');
                    } else if (binary) {
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
                    fieldCache.emplace_back((const char*)tmpBuf, (size_t)len);
                    break;
                case 7:
                    if (binary) {
                        int64_t pnextBin = (int64_t)(uint32_t)(*(uint32_t*)tmpBuf);
                        nextMappedPos[lineNo] = pnextBin;
                        fieldCache.emplace_back(std::to_string(pnextBin) + '\t');
                    } else if (mode == "pnext_delta") {
                        int64_t delta = (int64_t)std::stoll(std::string((char*)tmpBuf, len - 1));
                        int64_t pos = mappedPos.count(lineNo) ? mappedPos[lineNo] : 0;
                        int64_t pnext = delta + pos;
                        nextMappedPos[lineNo] = pnext;
                        fieldCache.emplace_back(std::to_string(pnext) + '\t');
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

    /* 完整伙伴索引：所有 (pos, pnext) 已知后一次性建立。 */
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
    // Reconstruct ID from split segments
    for (uint32_t splitIdx = 0; splitIdx < idStreams.size(); ++splitIdx) {
        Json::Value& splitMeta = idStreams[splitIdx];
        uint32_t splitDstLen = splitMeta["dstlen"].asUInt();
        std::string coderName = splitMeta["coder"]["magic"].asString();

        // Decode segment
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
    if (fieldDecoders[fieldIdx]->decode_line((uint8_t*)&chrIndex, sizeof(chrIndex), UINT8_MAX, false) < 0) {
        LOG_ERROR("Decode chr name field(%u) failed, lineNo = %u", fieldIdx, lineNo);
        return -1;
    }
    if (chrIndex == 0xFFFF) {
        if (fieldIdx == 6) {
            nextMappedChr[lineNo] = 0xFFFF;
        }
        *outputBlock->getCurrent() = '*';
        *(outputBlock->getCurrent() + 1) = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 2);
        return 2;
    } else if (chrIndex == 0xFFFE) {
        if (fieldIdx == 6) {
            nextMappedChr[lineNo] = 0xFFFE;
        }
        *outputBlock->getCurrent() = '=';
        *(outputBlock->getCurrent() + 1) = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 2);
        return 2;
    } else {
        if (fieldIdx == 2) {
            mappedChr[lineNo] = chrIndex;
        } else if (fieldIdx == 6) {
            nextMappedChr[lineNo] = chrIndex;
        }
        std::string chrName = SamInfo::getInstance().getChromosomeInfo(chrIndex).name;
        memcpy(outputBlock->getCurrent(), chrName.c_str(), chrName.length());
        outputBlock->setDataLen(outputBlock->getDataLen() + chrName.length());
        *outputBlock->getCurrent() = '\t';
        outputBlock->setDataLen(outputBlock->getDataLen() + 1);
        return  chrName.length() + 1;
    }
}

int32_t SamCodecActuator::decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& nposOffset, uint32_t& totalBaseLen, RoughIOBlock* outputBlock) {
    /*
     * 这里**不能**再 ensureCapacity：调用方 decompressSamByFields 在 SEQ 阶段抓了
     * basePtr = outputBlock->getCurrent() 留给 QUAL 用，块内任何 realloc 都会让它
     * 悬空，产出随机错误内容。缓冲够用由块入口一次性预分配保证（按文件头 block_size
     * 的确定上界 ×2，见 decompress()），逐行不再扩容。
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
                int32_t decLen = fieldDecoders[fieldIdx]->decode_line(outputBlock->getCurrent(), maxBaseLength, '\t', false);
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
        if (baseStream[0]["coder"]["magic"].asString() == "coder_bwt_cm") {
            int32_t decoderLen = 0;
            uint16_t mapFlag = mappedFlag.find(lineNo) == mappedFlag.end() ? 4 : mappedFlag[lineNo];
            // Not matched
            if (mapFlag & 0x04) {
                decoderLen = fieldDecoders[fieldIdx]->decode_line(baseSquashBuffer, actualBaseLen, UINT8_MAX, false);
                if (decoderLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                    return -1;
                }
                for (uint32_t o = 0; o < decoderLen; ++o) {
                    outputBlock->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
                }
            } else {
                // Get position in reference genome
                bool findMappedPos = false;
                int64_t refeMappedPos = 0;
                do {
                    uint16_t chrIdx = mappedChr.find(lineNo) == mappedChr.end() ? 0xFFFF : mappedChr[lineNo];
                    if (chrIdx ==  0xFFFF) {
                        break;
                    }
                    int64_t refeChrPos = SamInfo::getInstance().getPositionByIndex(chrIdx);
                    if (refeChrPos == -1) {
                        break;
                    }
                    refeMappedPos = refeChrPos + mappedPos[lineNo] - 1;
                    uint32_t baseSquashLength = actualBaseLen / 4 + !!(actualBaseLen & 0x3) + 1;
                    if ((refeMappedPos / 4) + baseSquashLength > pRefeGene->getSquashLength()) {
                        break;
                    }
                    findMappedPos = true;
                } while(0);

                if (!findMappedPos) {
                    decoderLen = fieldDecoders[fieldIdx]->decode_line(baseSquashBuffer, actualBaseLen, UINT8_MAX, false);
                    if (decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    for (uint32_t o = 0; o < decoderLen; ++o) {
                        outputBlock->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
                    }
                } else {
                    decoderLen = fieldDecoders[fieldIdx]->decode_line(baseDiffSquashBuffer, actualBaseLen, UINT8_MAX, false);
                    if (decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d,expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    pRefeGene->getStretch2Bits1Char(refeStrecchBuffer, actualBaseLen, refeMappedPos);
                    actgXor(refeStrecchBuffer, baseDiffSquashBuffer, baseSquashBuffer, actualBaseLen);
                    pRefeGene->getActgFrom2Bits(baseSquashBuffer, actualBaseLen, outputBlock->getCurrent());
                }
            }

            // Fill N back
            for (uint32_t n = 0; n < actualBaseLen; ++n) {
                if (nposOffset < baseNCount && baseNPosBuffer[nposOffset] == totalBaseLen + n) {
                    *(outputBlock->getCurrent() + n) = 'N';
                    nposOffset++;
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
    if (qualFcv2Decoder != nullptr) {
        /* fcv2 不需要 SEQ 作为上下文，链方向也由它自己的码流带着，只要长度。 */
        if (qualFcv2Decoder->decode_record(outputBlock->getCurrent(), actualBaseLen) < 0) {
            LOG_ERROR("Decode quality by fcv2 failed, len = %u", actualBaseLen);
            return -1;
        }
    } else if (qualCmDecoder != nullptr) {
        /*
         * bwt_cm 同样不需要 SEQ 上下文。压缩侧是按记录逐条 encode_line 喂进去的，
         * 这里也按同样的长度逐条取回；不传分隔符，长度由调用方给定。
         */
        if (qualCmDecoder->decode_line(outputBlock->getCurrent(), actualBaseLen) < 0) {
            LOG_ERROR("Decode quality by bwt_cm failed, len = %u", actualBaseLen);
            return -1;
        }
    } else {
        qualCoder->decode_qual_gen2(basePtr, outputBlock->getCurrent(), actualBaseLen);
    }
    outputBlock->setDataLen(outputBlock->getDataLen() + actualBaseLen);
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
    if (fieldLen > 1) {
        uint32_t seqLength = parseCigar(outputBlock->getCurrent(), fieldLen);
        baseLengthBuffer[lineIdx] = seqLength;
        cigarReadLen[lineIdx] = parseCigarRefConsumed(outputBlock->getCurrent(), fieldLen);
    } else {
        baseLengthBuffer[lineIdx] = 0;
        cigarReadLen[lineIdx] = 0;
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