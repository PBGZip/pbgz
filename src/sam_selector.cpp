#include "sam_selector.h"
#include <stdexcept>

namespace sam_compressor {

bool SamSelectorManager::hasSelectionForField(uint32_t fieldIdx) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return fieldSelections_.find(fieldIdx) != fieldSelections_.end();
}

CompressionModel SamSelectorManager::getSelectionForField(uint32_t fieldIdx) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = fieldSelections_.find(fieldIdx);
    if (it != fieldSelections_.end()) {
        return it->second;
    }
    // Default to BWT_CM if no selection exists
    return CompressionModel::MODEL_BWT_CM;
}

void SamSelectorManager::setSelectionForField(uint32_t fieldIdx, CompressionModel model) {
    std::lock_guard<std::mutex> lock(mutex_);
    fieldSelections_[fieldIdx] = model;
}

void SamSelectorManager::clearAllSelections() {
    std::lock_guard<std::mutex> lock(mutex_);
    fieldSelections_.clear();
    fieldStats_.clear();
}

SamSelectorManager::CompressionStats SamSelectorManager::getStatsForField(uint32_t fieldIdx) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = fieldStats_.find(fieldIdx);
    if (it != fieldStats_.end()) {
        return it->second;
    }
    // Return default stats if not found
    return CompressionStats{1.0, 0, 0, 0};
}

void SamSelectorManager::updateStatsForField(uint32_t fieldIdx, const CompressionStats& stats) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = fieldStats_.find(fieldIdx);
    if (it != fieldStats_.end()) {
        // Update existing stats with weighted average
        CompressionStats& existing = it->second;
        uint64_t totalTrials = existing.trialCount + stats.trialCount;
        existing.compressionRatio = (existing.compressionRatio * existing.trialCount + 
                                    stats.compressionRatio * stats.trialCount) / totalTrials;
        existing.sourceSize += stats.sourceSize;
        existing.compressedSize += stats.compressedSize;
        existing.trialCount = totalTrials;
    } else {
        // Store new stats
        fieldStats_[fieldIdx] = stats;
    }
}

} // namespace sam_compressor