/*
 * compress_selector.h - Header file for compression selector
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

#pragma once

#include <map>
#include <set>
#include <vector>
#include <cstdint>
#include <mutex>
#include <string>
#include "io_block.h"

// Compression model enumeration
enum class CompressionModel : uint32_t {
    MODEL_QUAL_GEN2 = 0,     // Original QUAL compression with frequency table
    MODEL_BWT_CM = 1,        // BWT compression with context modeling
    MODEL_AFFIX_MATCH = 2    // Affix compression with match optimization
};

enum class CompressionField : uint32_t {
    // FASTQ fields
    FASTQ_ID = 0,
    FASTQ_BASE = 1,
    FASTQ_COMMENT = 2,
    FASTQ_QUALITY = 3,
    
    // SAM fields
    SAM_QNAME = 4,
    SAM_FLAG = 5,
    SAM_RNAME = 6,
    SAM_POS = 7,
    SAM_MAPQ = 8,
    SAM_CIGAR = 9,
    SAM_RNEXT = 10,
    SAM_PNEXT = 11,
    SAM_TLEN = 12,
    SAM_SEQ = 13,
    SAM_QUAL = 14,
    SAM_OPTION = 15,
};

// Compression selector for SAM fields
class CompressionSelectorManager {
public:
    static CompressionSelectorManager& getInstance() {
        static CompressionSelectorManager instance;
        return instance;
    }
    
    // Prevent copy and assignment
    CompressionSelectorManager(const CompressionSelectorManager&) = delete;
    CompressionSelectorManager& operator=(const CompressionSelectorManager&) = delete;
    
    // Clear all statistics
    void clearAllSelections();
    
    // Get statistics about compression performance
    struct CompressionStats {
        double compressionRatio;
        uint64_t sourceSize;
        uint64_t compressedSize;
        uint32_t modelType;
        
        // Comparison operator for strict weak ordering - smaller compressedSize is better
        bool operator<(const CompressionStats& other) const {
            if (compressedSize != other.compressedSize)
                return compressedSize < other.compressedSize;
            return modelType < other.modelType;
        }
    };
    
    std::vector<CompressionStats> getStatsForField(uint32_t fieldIdx) const;
    void addStatsForField(uint32_t fieldIdx, const CompressionStats& stats);
    
    // Get the best compression model for a field based on compression ratio
    CompressionModel getBestModelForField(uint32_t fieldIdx) const;
    
    // Print compression statistics for all fields with selected encoders
    void printCompressionStats() const;

private:
    CompressionSelectorManager() = default;
    ~CompressionSelectorManager() = default;
    
    mutable std::mutex mutex_;
    std::map<uint32_t, std::set<CompressionStats>> fieldStats;
};

inline const char* modelName(CompressionModel model) {
    switch (model) {
        case CompressionModel::MODEL_QUAL_GEN2: return "QUAL_GEN2";
        case CompressionModel::MODEL_BWT_CM:    return "BWT_CM";
        case CompressionModel::MODEL_AFFIX_MATCH: return "AFFIX_MATCH";
        default: return "UNKNOWN";
    }
}

inline std::map<CompressionField, std::vector<CompressionModel>> compressFieldConfg = {
    {CompressionField::SAM_FLAG, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_POS, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_MAPQ, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_CIGAR, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_PNEXT, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_TLEN, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}},
    {CompressionField::SAM_QUAL, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_QUAL_GEN2}},
    {CompressionField::SAM_OPTION, {CompressionModel::MODEL_BWT_CM, CompressionModel::MODEL_AFFIX_MATCH}}
};

class SamCompressionSlector {
public:
    // Get CompressionField enum from SAM field index (0-based)
    static CompressionField getSelctorField(uint32_t samFieldIdx);

    static void testSamRegularFiled(uint32_t samFieldIdx, RoughIOBlock* blockPtr);

    static void testSamRegularFiled(RoughIOBlock* inBlockPtr);
};
