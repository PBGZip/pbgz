/*
 * compress_selector.cpp - Implementation for compression selector
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

#include "compress_selector.h"
#include <stdexcept>
#include <algorithm>
#include <vector>
#include "log/logger.h"
#include "io_block.h"
#include "coder/coder_io.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_affix_match.h"
#include "coder/coder_qual.h"
#include "utils/memory_util.h"


void CompressionSelectorManager::clearAllSelections() {
    std::lock_guard<std::mutex> lock(mutex_);
    fieldStats.clear();
}

std::vector<CompressionSelectorManager::CompressionStats> CompressionSelectorManager::getStatsForField(uint32_t fieldIdx) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = fieldStats.find(fieldIdx);
    if (it != fieldStats.end()) {
        return std::vector<CompressionStats>(it->second.begin(), it->second.end());
    }
    // Return empty vector if not found
    return std::vector<CompressionStats>();
}

void CompressionSelectorManager::addStatsForField(uint32_t fieldIdx, const CompressionStats& stats) {
    std::lock_guard<std::mutex> lock(mutex_);
    fieldStats[fieldIdx].insert(stats);
}

CompressionModel CompressionSelectorManager::getBestModelForField(uint32_t fieldIdx) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = fieldStats.find(fieldIdx);
    if (it == fieldStats.end() || it->second.empty()) {
        return CompressionModel::MODEL_BWT_CM;
    }
    
    // Best is the first element (set is sorted by compressedSize)
    return static_cast<CompressionModel>(it->second.begin()->modelType);
}

void CompressionSelectorManager::printCompressionStats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    LOG_INFO("=== Compression Statistics by Field ===");
    
    for (const auto& fieldEntry : fieldStats) {
        uint32_t fieldIdx = fieldEntry.first;
        
        if (fieldEntry.second.empty()) {
            continue;
        }
        
        auto& bestStats = *fieldEntry.second.begin();
        LOG_INFO("Field %u, Model %s, Ratio: %.2f%%, Original size: %lu bytes, Compressed size: %lu bytes",
                 fieldIdx, modelName(static_cast<CompressionModel>(bestStats.modelType)),
                 bestStats.compressionRatio, bestStats.sourceSize, bestStats.compressedSize);
    }
    
    LOG_INFO("=======================================");
}

CompressionField SamCompressionSlector::getSelctorField(uint32_t samFieldIdx) {
     switch (samFieldIdx) {
        case 0: return CompressionField::SAM_QNAME;
        case 1: return CompressionField::SAM_FLAG;
        case 2: return CompressionField::SAM_RNAME;
        case 3: return CompressionField::SAM_POS;
        case 4: return CompressionField::SAM_MAPQ;
        case 5: return CompressionField::SAM_CIGAR;
        case 6: return CompressionField::SAM_RNEXT;
        case 7: return CompressionField::SAM_PNEXT;
        case 8: return CompressionField::SAM_TLEN;
        case 9: return CompressionField::SAM_SEQ;
        case 10: return CompressionField::SAM_QUAL;
        default: return CompressionField::SAM_OPTION;
    }
}

void SamCompressionSlector::testSamRegularFiled(uint32_t samFieldIdx, RoughIOBlock* blockPtr) {
    CompressionField fieldEnum = getSelctorField(samFieldIdx);
    auto configIt = compressFieldConfg.find(fieldEnum);
    if (configIt == compressFieldConfg.end()) {
        return;
    }
    std::vector<CompressionModel>& encoderList = configIt->second;
    if (encoderList.empty()) {
        return;
    }
    
    CompressionSelectorManager& selector = CompressionSelectorManager::getInstance();
    
    // Get SAM data from block
    uint8_t* buffer = blockPtr->getBuffer();
    std::vector<uint32_t>& npos = blockPtr->getNpos();
    uint32_t lineNum = npos.size();
    
    // Extract field data from SAM lines
    std::vector<std::vector<uint8_t>> fieldDataList;
    std::vector<std::vector<uint8_t>> seqDataList; // For QUAL_GEN2 context
    uint32_t totalSrcLen = 0;
    
    // Build quality frequency table for QUAL_GEN2 (same as SamCodecActuator::preAnalysis)
    std::vector<std::pair<uint16_t, uint16_t>> qualFreqTable;
    std::pair<uint8_t, uint32_t> qualityFrequnce[256];
    for (int32_t k = 0; k < 256; ++k) {
        qualityFrequnce[k].first = k;
        qualityFrequnce[k].second = 0;
    }
    
    for (uint32_t i = 0; i < lineNum; ++i) {
        uint32_t lineStart = (i == 0) ? 0 : npos[i - 1] + 1;
        uint32_t lineEnd = npos[i] - lineStart;
        uint8_t* line = buffer + lineStart;
        
        // Skip header lines
        if (line[0] == '@') {
            continue;
        }
        
        // Find tab positions in this line
        std::vector<uint32_t> tabPos;
        for (uint32_t j = 0; j < lineEnd; ++j) {
            if (line[j] == '\t') {
                tabPos.push_back(j);
            }
        }
        
        // Extract field data based on field index
        if (samFieldIdx > tabPos.size()) {
            continue; // Field doesn't exist in this line
        }
        
        uint32_t fieldStart = (samFieldIdx == 0) ? 0 : tabPos[samFieldIdx - 1] + 1;
        uint32_t fieldEnd = (samFieldIdx < tabPos.size()) ? tabPos[samFieldIdx] : lineEnd;
        
        if (fieldEnd > fieldStart) {
            uint32_t fieldLen = fieldEnd - fieldStart;
            std::vector<uint8_t> fieldData(line + fieldStart, line + fieldStart + fieldLen);
            fieldDataList.push_back(fieldData);
            totalSrcLen += fieldLen;
            
            // Extract SEQ field (field 9) for QUAL_GEN2 context
            bool skipQUAL = false;
            if (samFieldIdx == 10) {
                // Collect quality frequency from actual data
                for (uint32_t p = 0; p < fieldLen; ++p) {
                    uint8_t ch = fieldData[p];
                    qualityFrequnce[ch].second++;
                }
                
                if (tabPos.size() >= 10) {
                    uint32_t seqStart = (tabPos[8] + 1 < lineEnd) ? tabPos[8] + 1 : 0;
                    uint32_t seqEnd = (tabPos[9] < lineEnd) ? tabPos[9] : lineEnd;
                    if (seqEnd > seqStart) {
                        uint32_t seqLen = seqEnd - seqStart;
                        std::vector<uint8_t> seqData(line + seqStart, line + seqStart + seqLen);
                        // Check if SEQ length matches QUAL length
                        if (seqLen != fieldLen) {
                            skipQUAL = true;
                        } else {
                            seqDataList.push_back(seqData);
                        }
                    } else {
                        skipQUAL = true;
                    }
                } else {
                    skipQUAL = true;
                }
                
                if (skipQUAL) {
                    // Remove the already added QUAL data from fieldDataList
                    fieldDataList.pop_back();
                    totalSrcLen -= fieldLen;
                    continue;
                }
            }
        }
    }
    
    if (fieldDataList.empty()) {
        LOG_WARNING("No data found for SAM field %d", samFieldIdx);
        return;
    }
    
    // Build qualFreqTable from collected quality frequency (same as SamCodecActuator::preAnalysis)
    if (samFieldIdx == 10) {
        std::sort(qualityFrequnce, qualityFrequnce + 256, 
            [](const std::pair<uint8_t, uint32_t> &a, const std::pair<uint8_t, uint32_t> &b) { return a.second > b.second;});
        for (int k = 0; k < 256; k++) {
            if (qualityFrequnce[k].second == 0) {
                continue;
            }
            qualFreqTable.push_back(std::make_pair(qualityFrequnce[k].first - '!', 1));
        }
    }
    
    // Test each encoder
    for (const auto& model : encoderList) {
        // Allocate temporary output block for compression
        RoughIOBlock* outBlock = new RoughIOBlock(totalSrcLen * 2 + 1024);
        
        // Create output IO
        std::shared_ptr<coder_io> outputIo = std::make_shared<coder_io>(
            outBlock->getCurrent(), 
            outBlock->getRemain()
        );
        
        try {
            // Perform compression based on model type
            uint32_t totalDstLen = 0;
            
            if (model == CompressionModel::MODEL_BWT_CM) {
                std::shared_ptr<coder_bwt_cm> coder = std::make_shared<coder_bwt_cm>(outputIo.get());
                for (const auto& fieldData : fieldDataList) {
                    coder->encode_line(fieldData.data(), static_cast<uint32_t>(fieldData.size()));
                }
                coder->encode_flush();
            } else if (model == CompressionModel::MODEL_AFFIX_MATCH) {
                std::shared_ptr<coder_affix_match> coder = std::make_shared<coder_affix_match>(outputIo.get());
                for (const auto& fieldData : fieldDataList) {
                    coder->encode_line(fieldData.data(), static_cast<uint32_t>(fieldData.size()));
                }
                coder->encode_flush();
            } else if (model == CompressionModel::MODEL_QUAL_GEN2) {
                std::shared_ptr<coder_qual> coder = std::make_shared<coder_qual>(outputIo.get(), true, qualFreqTable);
                // For QUAL_GEN2, use SEQ as context - lined up with QUAL data
                size_t seqIdx = 0;
                for (const auto& fieldData : fieldDataList) {
                    uint8_t* seqPtr = nullptr;
                    if (seqIdx < seqDataList.size() && !seqDataList[seqIdx].empty()) {
                        seqPtr = const_cast<uint8_t*>(seqDataList[seqIdx].data());
                    }
                    coder->encode_qual_gen2(seqPtr, const_cast<uint8_t*>(fieldData.data()), static_cast<uint32_t>(fieldData.size()));
                    seqIdx++;
                }
                coder->encode_flush();
            }
        
            totalDstLen = outputIo->data_len;
            outBlock->setDataLen(totalDstLen);
            
            // Calculate compression statistics
            CompressionSelectorManager::CompressionStats stats;
            stats.modelType = static_cast<uint32_t>(model);
            stats.sourceSize = totalSrcLen;
            stats.compressedSize = totalDstLen;
            stats.compressionRatio = (totalSrcLen > 0) ? (double)(totalDstLen * 100.0) / totalSrcLen : 0.0;
            
            selector.addStatsForField(samFieldIdx, stats);
            
            LOG_DEBUG("Field %d with model %s: %u -> %u bytes, ratio %.2f%%", 
                samFieldIdx, modelName(model), totalSrcLen, totalDstLen, stats.compressionRatio);
            
        } catch (const std::exception& e) {
            LOG_ERROR("Compression test failed for field %d with model %s: %s", 
                samFieldIdx, modelName(model), e.what());
        }
        delete outBlock;
    }
}

void SamCompressionSlector::testSamRegularFiled(RoughIOBlock* blockPtr){ 
    if (!blockPtr) {
        LOG_ERROR("Block pointer is null compression testing");
        return;
    }

    static bool isNeedTest = true;
    if (!isNeedTest) {
        return;
    }
    isNeedTest = false;

    for (uint32_t idx = 0; idx < 12; ++idx) {
        testSamRegularFiled(idx, blockPtr);
    }

    return;
}
