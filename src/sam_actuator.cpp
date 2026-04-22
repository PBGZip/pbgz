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
#include "utils/md5_util.h"
#include "coder/coder_json.h"
#include "log/logger.h"
#include "sam_info.h"
#include "actg.h"
#include "pbgz_manager.h"
#include "utils/path_util.h"

SamCodecActuator::SamCodecActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, Reference* pReferene): CodecActuator(inPtr, outPtr) {
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
    notifyFlag = false;
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

    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
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
                size_t pos = line.find("UR:");
                if (pos != std::string::npos) {
                    pos += 3;
                    size_t endPos = line.find("\t", pos);
                    if (endPos == std::string::npos) {
                        endPos = line.length();
                    }
                    std::string refFilePath = line.substr(pos, endPos - pos);
                    LOG_INFO("Reference fasta name: %s.", refFilePath.c_str());
                    if (pRefeGene != nullptr) {
                        std::string inputFastq = PathUtil::getFileName(pRefeGene->getFastaFileName());
                        if (refFilePath != inputFastq) {
                            fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fastq is %s. \n", refFilePath.c_str(), inputFastq.c_str());
                        }
                    }
                }
                continue;
            } else if (line.substr(0, 3) == "@PG") {
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
                        std::string refGeneName = tokens[2];
                        LOG_INFO("Reference gene name extracted from @PG CL: %s", refGeneName.c_str());
                        // Reference gene name can be saved or used here as needed
                        if (pRefeGene != nullptr) {
                            std::string inputFastq = PathUtil::getFileName(pRefeGene->getFastaFileName());
                            if (PathUtil::isGzFile(pRefeGene->getFastaFileName())) {
                                inputFastq = PathUtil::getFileNameFromGz(pRefeGene->getFastaFileName());
                            }
                            if (refGeneName != inputFastq) {
                                fprintf(stderr, "Warning: fasta file not match, SAM fasta is %s, input fastq is %s \n", refGeneName.c_str(), inputFastq.c_str());
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
                    } 
                    // CIGAR is the 6th field
                    if (linePos.size() == 5) {
                        uint32_t cigarLen = baseLen = i - linePos.at(4) - 1;
                        if (cigarLen > 1) {
                            lineCigarMatchFlag = true;
                        } else {
                            lineCigarMatchFlag = false;
                        }
                    }
                    // Base is the 10th field
                    if (linePos.size() == 9) {
                        baseLen = i - linePos.at(8) - 1;
                        if (baseLen > maxBaseLength) {
                            maxBaseLength = baseLen;
                        } 
                        if (!lineCigarMatchFlag) {
                            if (baseLen < minBaseLength) {
                                minBaseLength =  baseLen;
                            }
                        }
                    }
                    // Quality value is the 11th field
                    if (linePos.size() == 10) {
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
        return 0;
    }

    notifyFlag = true;

    // Check if preAnalysis was successful
    if (contentPos.empty()) {
        LOG_ERROR("preAnalysis() must be called before compress()");
        return -1;
    }

    return compressSamByFields();
}

int32_t SamCodecActuator::compressSamHeader() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter for SAM header compression");
        return -1;
    }

    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint8_t* buffer = inBlockPtr->getBuffer();
    
    // Clear previous chromosome information but keep chromosome ID counter for multi-block compatibility
    // Note: Do not clear chromosome information as multiple blocks may need to process chromosome info
    headerSrcLen = 0;
    // Create SAM file header compressor
    std::shared_ptr<coder_io> headerIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> headerCoder = std::make_shared<coder_bwt_cm>(headerIo.get());
    
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
    
    LOG_DEBUG("SAM header compression completed: %u lines, %u bytes -> %u bytes, compress ratio = %.2f%%", 
             headEndLine, headerSrcLen, headerDstLen, (double)(headerDstLen* 100)/(double)headerSrcLen);
    
    return 0;
}

int32_t SamCodecActuator::compressSamByFields() {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
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
                fieldDstLen = compressNumber<uint16_t>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 2: // RNAME
                fieldDstLen = compressChrName(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 3: // POS
                fieldDstLen = compressNumber<uint32_t>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 4: // MAPQ
                fieldDstLen = compressNumber<uint8_t>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 5: // CIGAR 
                fieldDstLen = compressCigar(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 6: // RNEXT
                fieldDstLen = compressChrName(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 7: // PNEXT
                fieldDstLen = compressNumber<uint32_t>(fieldIdx, fieldSrcLen, fieldMeta);
                break;
            case 8: // TLEN 
                fieldDstLen = compressNumber<int32_t>(fieldIdx, fieldSrcLen, fieldMeta);
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

        // LOG_INFO("Compress Rate for block(%d fieldId=%d): src = %d, dst = %d, ratio =  %0.2f%%.",
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

int32_t SamCodecActuator::compressIdFieldSplit(uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    // Similar to FastqActuator::compressIdInSplit, create multiple streams for each split symbol
    Json::Value streamMeta;
    uint32_t totalSrcLength = 0;
    uint32_t totalDstLength = 0;

    // Process each split symbol (similar to FastqActuator)
    for (uint32_t i = 0; i < idSplitSymbols.size(); ++i) {
        std::shared_ptr<coder_io> idIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        std::shared_ptr<coder_bwt_cm> idCoder = std::make_shared<coder_bwt_cm>(idIo.get());
        uint32_t srcLength = 0;
        // Process each line and compress the specific split segment
        std::vector<uint32_t>& npos = inBlockPtr->getNpos();
        uint32_t lineNum = npos.size();
        uint8_t* buffer = inBlockPtr->getBuffer();
        
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
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    
    // Create encoder for ID field whole compression
    std::shared_ptr<coder_io> fieldIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> fieldCoder = std::make_shared<coder_bwt_cm>(fieldIo.get());
    
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
    
    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);

    Json::Value streamMeta;
    streamMeta["srclen"] = fieldSrcLen;
    streamMeta["dstlen"] = fieldIo->data_len;
    streamMeta["coder"] = fieldIo->meta;
    
    // Set field metadata
    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["splitsym"] = "\t";
    fieldMeta["field"] = 0;

    LOG_DEBUG("SAM ID compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
            fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);
    
    return fieldIo->data_len;
}

int32_t SamCodecActuator::compressChrName(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();

    // Create encoder for regular field compression
    std::shared_ptr<coder_io> chrIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> chrCoder = std::make_shared<coder_bwt_cm>(chrIo.get());

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
        fieldSrcLen += sizeof(chrIndex);
    }
    
    // Flush the encoder
    chrCoder->encode_flush();
    
    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + chrIo->data_len);
    
    // Set field metadata
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = chrIo->data_len;
    fieldMeta["coder"] = chrIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_DEBUG("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
            fieldIdx, fieldSrcLen, chrIo->data_len, (double)(chrIo->data_len * 100)/(double)fieldSrcLen);

    return chrIo->data_len;
 }

int32_t SamCodecActuator::compressRegularField(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    
    // Create encoder for regular field compression
    std::shared_ptr<coder_io> fieldIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> fieldCoder = std::make_shared<coder_bwt_cm>(fieldIo.get());
    
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
        }
    }
    
    // Flush the encoder for this field
    fieldCoder->encode_flush();
    
    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);
    
    // Set field metadata
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_DEBUG("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

int32_t SamCodecActuator::compressCigar(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    
    // Create encoder for regular field compression
    std::shared_ptr<coder_io> fieldIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> fieldCoder = std::make_shared<coder_bwt_cm>(fieldIo.get());
    
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
        } else {
            uint32_t sequeceLength = parseCigar(fieldStart, fieldLength);
            baseLengthBuffer[lineIdx - headEndLine] = sequeceLength;
        }
        
        // Encode the field data
        if (fieldLength > 0) {
            fieldCoder->encode_line(fieldStart, fieldLength);
            fieldSrcLen += fieldLength;
        }
    }
    
    // Flush the encoder for this field
    fieldCoder->encode_flush();
    
    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);
    
    // Set field metadata
    fieldMeta["srclen"] = fieldSrcLen;
    fieldMeta["dstlen"] = fieldIo->data_len;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_DEBUG("SAM field(%d) compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
        fieldIdx, fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);

    return fieldIo->data_len;
}

int32_t SamCodecActuator::compressBaseWithoutRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
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
    std::shared_ptr<coder_io> fieldIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_fc> fieldCoder = std::make_shared<coder_fc>(fieldIo.get());

    fieldCoder->encode_line(tmpBuffer.get(), fieldSrcLen);
    fieldCoder->encode_flush();
    // Smart pointer automatically cleans up

    // Update output block data length
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldIo->data_len);
    
    // Set field metadata
    fieldMeta["minlen"] = minBaseLength;
    fieldMeta["maxlen"] = maxBaseLength;
    fieldMeta["totalsrclen"] = fieldSrcLen;
    fieldMeta["totaldstlen"] = fieldIo->data_len;;
    fieldMeta["coder"] = fieldIo->meta;
    fieldMeta["field"] = fieldIdx;

    LOG_DEBUG("SAM base field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
        fieldSrcLen, fieldIo->data_len, (double)(fieldIo->data_len * 100)/(double)fieldSrcLen);
    
    return fieldIo->data_len;
}

int32_t SamCodecActuator::compressBaseWithRef(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    if (pRefeGene == nullptr) {
        LOG_ERROR("Reference genome is not available for base compression with reference");
        return -1;
    }

    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
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
    std::shared_ptr<coder_io> matchIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> matchCm = std::make_shared<coder_bwt_cm>(matchIo.get());
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
        std::shared_ptr<coder_io> nposIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        uint32_t nposSrcLen = (baseNCount << 2); // 4 bytes per N position
        std::shared_ptr<coder_bwt_cm> subCoder = std::make_shared<coder_bwt_cm>(nposIo.get());
        subCoder->encode_line((uint8_t*)baseNPosBuffer, nposSrcLen);
        subCoder->encode_flush();
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
        std::shared_ptr<coder_io> lenIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
        uint32_t baseLenSrcLen = unmapedReadLength.size() << 1;
        uint32_t* baseLenBuffer = MemoryUtil::safeAlloc<uint32_t>(baseLenSrcLen);
        if (baseLenBuffer == nullptr) {
            return -1;
        }
        for (uint32_t i = 0; i < unmapedReadLength.size(); i += 2) {
            std::pair baseLenPos = unmapedReadLength[i];
            baseLenBuffer[i] = baseLenPos.first;
            baseLenBuffer[i + 1] = baseLenPos.second;
        }
        std::shared_ptr<coder_bwt_cm> lenCoder = std::make_shared<coder_bwt_cm>(lenIo.get());
        lenCoder->decode_line((uint8_t*)baseLenBuffer, baseLenSrcLen<<2);
        lenCoder->encode_flush();

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
    fieldSrcLen = totalSrcLen;

    LOG_DEBUG("SAM base field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
        totalSrcLen, totalDstLen, (double)(totalDstLen * 100)/(double)totalSrcLen);
    
    return totalDstLen;
}

int32_t SamCodecActuator::compressQuality(uint32_t fieldIdx, uint32_t& fieldSrcLen, Json::Value& fieldMeta) {
    std::vector<uint32_t>& npos = inBlockPtr->getNpos();
    uint32_t lineNum = npos.size();
    uint8_t* buffer = inBlockPtr->getBuffer();
    
    // Create quality encoder similar to FastqActuator
    std::shared_ptr<coder_io> qualityIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_qual> qualityCoder = std::make_shared<coder_qual>(qualityIo.get(), true, qualFreqTable);

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
        
        // Encode quality with sequence context
        qualityCoder->encode_qual_gen2(seqStart, qualStart, qualLength);
        streamSrcLen += qualLength;
    }
    
    qualityCoder->encode_flush();
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityIo->data_len);
    
    Json::Value subMeta;
    subMeta["srclen"] = streamSrcLen;
    subMeta["dstlen"] = qualityIo->data_len;
    subMeta["coder"] = qualityIo->meta;
    streamMeta.append(subMeta);

    totalSrcLength += streamSrcLen;
    totalDstLength += qualityIo->data_len;

    // Encode quality frequency table
    std::shared_ptr<coder_io> qualityFreqIo = std::make_shared<coder_io>(outBlockPtr->getCurrent(), outBlockPtr->getRemain());
    std::shared_ptr<coder_bwt_cm> qualityFreqCoder = std::make_shared<coder_bwt_cm>(qualityFreqIo.get());
    std::shared_ptr<uint16_t[]> qualityFreqArray(new uint16_t[qualFreqTable.size() << 1]);
    for (uint32_t i = 0; i < qualFreqTable.size(); ++i) {
        int idx = i << 1;
        qualityFreqArray[idx] = qualFreqTable[i].first;
        qualityFreqArray[idx + 1] = qualFreqTable[i].second;
    }

    uint32_t freqSrcLen = (qualFreqTable.size() << 1) * sizeof(uint16_t);
    qualityFreqCoder->encode_line((uint8_t*)qualityFreqArray.get(), freqSrcLen);
    qualityFreqCoder->encode_flush();
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + qualityFreqIo->data_len);

    subMeta.clear();
    subMeta["srclen"] = freqSrcLen;
    subMeta["dstlen"] = qualityFreqIo->data_len;
    subMeta["coder"] = qualityFreqIo->meta;
    subMeta["streamname"] = "qualityfreq";
    streamMeta.append(subMeta);

    totalSrcLength += freqSrcLen;
    totalDstLength += qualityFreqIo->data_len;

    // Set field metadata
    fieldMeta["totalsrclen"] = totalSrcLength;
    fieldMeta["totaldstlen"] = totalDstLength;
    fieldMeta["streams"] = streamMeta;
    fieldMeta["field"] = fieldIdx;
    
    fieldSrcLen = totalSrcLength;

    LOG_DEBUG("SAM quality field compression completed: %u bytes -> %u bytes, compress ratio = %.2f%%", 
        totalSrcLength, totalDstLength, (double)(totalDstLength * 100)/(double)totalSrcLength);
    
    return totalDstLength;
}

int32_t SamCodecActuator::decompress() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter, inBlockPtr or outBlockPtr is nullptr for SAM decompression");
        return -1;
    }
    // Parse meta information
    coder_json metaCoder;
    metaCoder.decoder(inBlockPtr->getMetaBuffer(), inBlockPtr->getMetaLen(), meta);
    if (meta.isMember("header")) {
        if (0 != decompressHeader()) {
            LOG_ERROR("Decompress header failed. block id = %d", inBlockPtr->getBlockId());
            return -1;
        }
    } else {
        LOG_DEBUG("No header info for block: %d", inBlockPtr->getBlockId());
        notifyFlag = true;
    }

    if (0 != decompressSamByFields()) {
        LOG_ERROR("Decompress fields failed. block id = %d", inBlockPtr->getBlockId());
        return -1;
    }

    // Verify checksum of decompressed content
    std::string md5;
    calcMd5sum(md5, outBlockPtr->getBuffer(), outBlockPtr->getDataLen());
    if (md5 != meta["md5"].asString()) {
        LOG_ERROR("MD5 check failed for SAM data, blockid = %d", outBlockPtr->getBlockId());
        return -1;
    }

    return 0;
}

int32_t SamCodecActuator::decompressSamByFields() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
        LOG_ERROR("Invalid parameter, inBlockPtr or outBlockPtr is nullptr for SAM field-by-field decompression");
        return -1;
    }

    if (!meta.isMember("sam")) {
        LOG_INFO("No SAM info for field-by-field decompression");
        return 0;
    }

    // Initialize decoders based on compression metadata
    if (0 != initDecoder(outBlockPtr)) {
        LOG_ERROR("Init decoder failed.");
        return -1;
    }

    // Set block information
    outBlockPtr->setBlockId(inBlockPtr->getBlockId());
    outBlockPtr->setBlockType(inBlockPtr->getBlockType());
    
    Json::Value& samMeta = meta["sam"];
    Json::Value& streams = samMeta["streams"];
    uint32_t fieldCount = samMeta["fieldcount"].asUInt();
    uint32_t lineNum = samMeta["lines"].asUInt();
    uint8_t* pBaseEnd = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize();
    baseSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    baseDiffSquashBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    refeStrecchBuffer = MemoryUtil::safeAlloc<uint8_t>(maxBaseLength);
    uint32_t totalBaseLen = 0;
    uint32_t nposOffset = 0;

    uint8_t* pBaseOut = nullptr;
    if (streams[9]["coder"]["magic"].asString() == "coder_fc") {
        pBaseOut = pBaseEnd - streams[9]["totalsrclen"].asUInt();
    }

    for (uint32_t lineNo = 0; lineNo < lineNum; ++lineNo) {
        uint8_t* basePtr = nullptr;
        uint32_t actualBaseLen = 0;
        // Decode each field for this line
        for (uint32_t fieldIdx = 0; fieldIdx < fieldCount; ++fieldIdx) {
            int32_t decoderLen = 0;
            if (fieldIdx == 0) {    /// ID 
                decoderLen = decompressIdField(fieldIdx, streams[fieldIdx]);
            } else if (fieldIdx == 1) {  /// FLAG
                decoderLen = decompressNumber<uint16_t>(fieldIdx, lineNo);
            } else if (fieldIdx == 2) {  /// RNAME
                decoderLen = decompressChrName(fieldIdx, lineNo);
            } else if (fieldIdx == 3) {  /// POS
                decoderLen = decompressNumber<uint32_t>(fieldIdx, lineNo);
            } else if (fieldIdx == 4) {  /// MAPQ
                decoderLen = decompressNumber<uint8_t>(fieldIdx, lineNo);
            } else if (fieldIdx == 5) {  /// CIGAR
                decoderLen = decompressCigar(fieldIdx, '\t', lineNo);
            } else if (fieldIdx == 6) {  /// RNEXT
                decoderLen = decompressChrName(fieldIdx, lineNo);
            } else if (fieldIdx == 7) {  /// PNEXT
                decoderLen = decompressNumber<uint32_t>(fieldIdx, lineNo);
            } else if (fieldIdx == 8) {  /// TLEN
                decoderLen = decompressNumber<int32_t>(fieldIdx, lineNo);
            } else if (fieldIdx == 9) {  /// SEQ
                basePtr = outBlockPtr->getCurrent();
                decoderLen = decompressBase(fieldIdx, streams[fieldIdx], pBaseOut, lineNo, nposOffset, totalBaseLen);
                actualBaseLen = decoderLen;
            } else if (fieldIdx == 10 ) {  /// QUAL  

                // LOG_DEBUG("lineNo = %d, actualBaseLen = %d", lineNo, actualBaseLen); 

                // PbgzManager::getInstance().printBufferContent(basePtr, actualBaseLen);
                
                decompressQuality(basePtr, actualBaseLen);
                // No optional fields scenario, replace appended \t with \n
                if (fieldIdx + 1 == fieldCount) {
                    uint8_t* pEnd = outBlockPtr->getCurrent();
                    *(pEnd - 1) = '\n';
                }
            } else {   /// Optional fields
                // Decode field until tab or end
                if (fieldIdx + 1 == fieldCount) {
                    decoderLen = decompressRegularField(fieldIdx, '\n');
                    // Single newline means it's appended, need to change \t after quality to \n and remove appended \n
                    if (decoderLen == 1) {
                        uint8_t* pEnd = outBlockPtr->getCurrent();
                        *(pEnd - 2) = '\n';
                        outBlockPtr->setDataLen(outBlockPtr->getDataLen() - 1);
                    }
                } else {
                    decoderLen = decompressRegularField(fieldIdx, '\t');
                }
            }

            if (decoderLen < 0) {
                LOG_ERROR("Decode field(%d) failed. lineNo = %d", fieldIdx, lineNo);
                return -1;
            }
        }
    }

    return 0;
}

int32_t SamCodecActuator::decompressHeader() {
    if (inBlockPtr == nullptr || outBlockPtr == nullptr) {
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
    
    headEndLine = headerMeta["lines"].asInt64();
    uint32_t dstLen = headerMeta["dstlen"].asUInt();

    // Create SAM file header decompressor
    std::shared_ptr<coder_io> headerIo = std::make_shared<coder_io>(inBlockPtr->getBuffer(), dstLen);
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
        uint32_t decodedLen = headerDecoder->decode_line(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), '\n', false);
        if (decodedLen == 0) {
            break; // No more data
        }
        std::string headStr = std::string((char*)outBlockPtr->getCurrent(), decodedLen);
        if (headStr.substr(0, 3) == "@SQ") {
            SamUtil::parseChromosomeInfo(headStr);
        }

        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + decodedLen);
        lineCount++;
        decoderTotalLen += decodedLen;
    }
    readOffset += dstLen;
    LOG_DEBUG("SAM header decompression completed: %u lines, %u bytes - > %u bytes.", headEndLine, dstLen, decoderTotalLen);
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
                std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLength);
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
                    coder_io baseIo(inBlockPtr->getBuffer() + readOffset, dstLength);
                    baseIo.meta = baseMeta;
                    baseIo.meta["dstlen"] = baseMeta["totaldstlen"].asUInt();
                    coder_fc baseDecoder = coder_fc(&baseIo);
                    baseDecoder.decode_line(outBlockPtr->getBuffer() + outBlockPtr->getBufferSize() - srcLength, srcLength, UINT8_MAX, false);
                } else if(coderName == "coder_bwt_cm") {
                    std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLength);
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
                std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLength);
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
                        coder_io nposIo(inBlockPtr->getBuffer() + readOffset, dstlen);
                        coder_bwt_cm nposCoder(&nposIo);
                        nposCoder.decode_line((uint8_t*)baseNPosBuffer, srclen, UINT8_MAX, false);
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
                        coder_io baseLenIo(inBlockPtr->getBuffer() + readOffset, dstlen);
                        coder_bwt_cm baseLenCoder(&baseLenIo);
                        baseLenCoder.decode_line((uint8_t*)baseLenBuffer, srclen, UINT8_MAX, false);
                    } else {
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

                coder_io qualFreqIo(inBlockPtr->getBuffer() + readOffset + qualDstLength, freqDstLength);
                coder_bwt_cm qualFreqCoder(&qualFreqIo);
                uint32_t qualFreqSrcLength = qualStreamMeta[1]["srclen"].asUInt();
                uint8_t qualFreqArrLength = qualFreqSrcLength / sizeof(uint16_t);
                uint16_t* qualFreqArr = new uint16_t[qualFreqArrLength];
                uint32_t qualFreq = qualFreqCoder.decode_line((uint8_t*)qualFreqArr,qualFreqSrcLength, UINT8_MAX, false);
                if (qualFreq != qualFreqSrcLength) {
                    LOG_ERROR("Decode quality frequncy failed");
                    delete [] qualFreqArr;
                    return -1;
                }
                for (uint32_t i = 0; i < qualFreqArrLength; i += 2) {
                    qualFreqTable.push_back(std::make_pair(qualFreqArr[i], qualFreqArr[i + 1]));
                }
                delete [] qualFreqArr;

                if (qualStreamMeta[0]["coder"]["magic"].asString() == "coder_qual") {
                    std::shared_ptr<coder_io> qualIo = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, qualDstLength);
                    ioVector.push_back(qualIo);
                    qualCoder = std::make_shared<coder_qual>(qualIo.get(), true, qualFreqTable);
                } else {
                    LOG_ERROR("Unsupport coder type: %s", streamMeta[0]["coder"]["magic"].asString().c_str());
                    return -1;
                }
                readOffset += (qualDstLength + freqDstLength);
            } else {
                LOG_ERROR("Unsupport coder type: %s", streamMeta[1]["coder"]["magic"].asString().c_str());
                return -1;
            }
        } else {
            std::string coderName = streamMeta[idx]["coder"]["magic"].asString();
            uint32_t dstLen = streamMeta[idx]["dstlen"].asUInt();
            if (coderName == "coder_bwt_cm") {
                std::shared_ptr<coder_io> io = std::make_shared<coder_io>(inBlockPtr->getBuffer() + readOffset, dstLen);
                ioVector.push_back(io);
                fieldDecoders[idx] = std::make_shared<coder_bwt_cm>(io.get());
            } else {
                LOG_ERROR("Unsupport coder type: %s", streamMeta[0]["coder"]["magic"].asString().c_str());
                return -1;
            }
            readOffset += dstLen;
        }
    }
    return 0;
}

int32_t SamCodecActuator::decompressRegularField(uint32_t fieldIdx, uint8_t splitFlag) {
    uint32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), splitFlag, false);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldLen);
    return fieldLen;
}

int32_t SamCodecActuator::decompressIdField(uint32_t fieldIdx, Json::Value& fieldMeta) {
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
        uint32_t segmentLen = idDecoders[splitIdx]->decode_line(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), 
            (splitIdx < idSplitSymbols.size()) ? idSplitSymbols[splitIdx] : UINT8_MAX, false);
        readOffset += splitDstLen;
        idLength += splitDstLen;
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + segmentLen);
    }
    return idLength;
}

int32_t SamCodecActuator::decompressChrName(uint32_t fieldIdx, uint32_t lineNo) {
    uint16_t chrIndex = 0;
    fieldDecoders[fieldIdx]->decode_line((uint8_t*)&chrIndex, sizeof(chrIndex), UINT8_MAX, false);
    if (chrIndex == 0xFFFF) {
        *outBlockPtr->getCurrent() = '*';
        *(outBlockPtr->getCurrent() + 1) = '\t';
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 2);
        return 2;
    } else if (chrIndex == 0xFFFE) {
        *outBlockPtr->getCurrent() = '=';
        *(outBlockPtr->getCurrent() + 1) = '\t';
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 2);
        return 2;
    } else {
        if (fieldIdx == 2) {
            mappedChr[lineNo] = chrIndex;
        } else if (fieldIdx == 6) {
            nextMappedChr[lineNo] = chrIndex;
        }
        std::string chrName = SamInfo::getInstance().getChromosomeInfo(chrIndex).name;
        memcpy(outBlockPtr->getCurrent(), chrName.c_str(), chrName.length());
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + chrName.length());
        *outBlockPtr->getCurrent() = '\t';
        outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
        return  chrName.length() + 1;
    }
}

int32_t SamCodecActuator::decompressBase(uint32_t fieldIdx, Json::Value& fieldMeta, uint8_t*& pBaseOut, uint32_t lineNo,
                                    uint32_t& nposOffset, uint32_t& totalBaseLen) {
    uint8_t* pBaseEnd = outBlockPtr->getBuffer() + outBlockPtr->getBufferSize();
    bool isUserReference = pRefeGene != nullptr && fieldMeta.isMember("streams");
    uint32_t actualBaseLen = 0;
    if (!isUserReference) {
        if (minBaseLength == maxBaseLength) {
            actualBaseLen = baseLengthBuffer[lineNo] == 0 ? maxBaseLength : baseLengthBuffer[lineNo];
            LOG_DEBUG("actualBaseLen = %d", actualBaseLen);
            if (fieldMeta["coder"]["magic"].asString() == "coder_fc") {
                memcpy(outBlockPtr->getCurrent(), pBaseOut, actualBaseLen);
                pBaseOut += actualBaseLen;
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
            } else if (fieldMeta["coder"]["magic"].asString() == "coder_bwt_cm") {
                fieldDecoders[fieldIdx]->decode_line(outBlockPtr->getCurrent(), actualBaseLen, UINT8_MAX, false);
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
            } else {
                LOG_ERROR("Not supported coder name:%s",fieldMeta["coder"]["magic"].asString().c_str());
                return -1;
            }

            *(outBlockPtr->getCurrent()) = '\t';
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
        } else {
            if (fieldMeta["coder"]["magic"].asString() == "coder_fc") {
                uint8_t* pBaseTmp = outBlockPtr->getCurrent();
                uint8_t* ptr = pBaseOut;
                for (; ptr < pBaseEnd; ++ptr) {
                    *pBaseTmp++ = *ptr;
                    if (*ptr == '\t') {
                        break;
                    }
                }
                actualBaseLen = ptr - pBaseOut + 1;
                pBaseOut += actualBaseLen;
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
                actualBaseLen -= 1; // Remove \t length
            } else if (fieldMeta["coder"]["magic"].asString() == "coder_bwt_cm") {
                actualBaseLen = fieldDecoders[fieldIdx]->decode_line(outBlockPtr->getCurrent(), maxBaseLength, '\t', false);
                outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
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
            uint32_t decoderLen = 0;
            uint16_t mapFlag = mappedFlag.find(lineNo) == mappedFlag.end() ? 4 : mappedFlag[lineNo];
            // Not matched
            if (mapFlag & 0x04) {
                decoderLen = fieldDecoders[fieldIdx]->decode_line(baseSquashBuffer, actualBaseLen, UINT8_MAX, false);
                if (decoderLen != actualBaseLen) {
                    LOG_ERROR("base decode failed in block %llu, line %d, expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                    return -1;
                }
                for (uint32_t o = 0; o < decoderLen; ++o) {
                    outBlockPtr->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
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
                        outBlockPtr->getCurrent()[o] = atcg4[baseSquashBuffer[o]];
                    }
                } else {
                    decoderLen = fieldDecoders[fieldIdx]->decode_line(baseDiffSquashBuffer, actualBaseLen, UINT8_MAX, false);
                    if (decoderLen != actualBaseLen) {
                        LOG_ERROR("base decode failed in block %llu, line %d,expect len %d, actural len %d", inBlockPtr->getBlockId(), lineNo, actualBaseLen, decoderLen);
                        return -1;
                    }
                    pRefeGene->getStretch2Bits1Char(refeStrecchBuffer, actualBaseLen, refeMappedPos);
                    actgXor(refeStrecchBuffer, baseDiffSquashBuffer, baseSquashBuffer, actualBaseLen);
                    pRefeGene->getActgFrom2Bits(baseSquashBuffer, actualBaseLen, outBlockPtr->getCurrent());   
                }
            }
            
            // Fill N back
            for (uint32_t n = 0; n < actualBaseLen; ++n) {
                if (nposOffset < baseNCount && baseNPosBuffer[nposOffset] == totalBaseLen + n) {
                    *(outBlockPtr->getCurrent() + n) = 'N';
                    nposOffset++;
                } 
            }
            totalBaseLen += actualBaseLen;
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);

            *(outBlockPtr->getCurrent()) = '\t';
            outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
        } else {
            LOG_ERROR("Not supported coder name:%s", fieldMeta["coder"]["magic"].asString().c_str());
            return -1;
        }
    }
    return actualBaseLen;
}

int32_t SamCodecActuator::decompressQuality(uint8_t* basePtr, uint32_t actualBaseLen) {
    qualCoder->decode_qual_gen2(basePtr, outBlockPtr->getCurrent(), actualBaseLen);
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + actualBaseLen);
    *(outBlockPtr->getCurrent()) = '\t';
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + 1);
    return actualBaseLen;
}

int32_t SamCodecActuator::decompressCigar(uint32_t fieldIdx, uint8_t splitFlag, uint32_t lineIdx) {
    uint32_t fieldLen = fieldDecoders[fieldIdx]->decode_line(outBlockPtr->getCurrent(), outBlockPtr->getRemain(), splitFlag, false);
    if (fieldLen > 1) {
        uint32_t seqLength = parseCigar(outBlockPtr->getCurrent(), fieldLen);
        baseLengthBuffer[lineIdx] = seqLength;
    } else {
        baseLengthBuffer[lineIdx] = 0;
    }
    
    outBlockPtr->setDataLen(outBlockPtr->getDataLen() + fieldLen);
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

