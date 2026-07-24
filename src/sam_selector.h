#ifndef SAM_SELECTOR_H
#define SAM_SELECTOR_H

#include <map>
#include <cstdint>
#include <mutex>

namespace sam_compressor {

// Compression model enumeration
enum class CompressionModel : uint32_t {
    MODEL_QUAL_GEN2 = 0,     // Original QUAL compression with frequency table
    MODEL_BWT_CM = 1,        // BWT compression with context modeling
    MODEL_AFFIX_MATCH = 2    // Affix compression with match optimization
};

// Compression selector for SAM fields
class SamSelectorManager {
public:
    static SamSelectorManager& getInstance() {
        static SamSelectorManager instance;
        return instance;
    }
    
    // Prevent copy and assignment
    SamSelectorManager(const SamSelectorManager&) = delete;
    SamSelectorManager& operator=(const SamSelectorManager&) = delete;
    
    // Check if a selection exists for a specific field
    bool hasSelectionForField(uint32_t fieldIdx) const;
    
    // Get the selected compression model for a field
    CompressionModel getSelectionForField(uint32_t fieldIdx) const;
    
    // Set the compression model for a specific field
    void setSelectionForField(uint32_t fieldIdx, CompressionModel model);
    
    // Clear all selections
    void clearAllSelections();
    
    // Get statistics about compression performance
    struct CompressionStats {
        double compressionRatio;
        uint64_t sourceSize;
        uint64_t compressedSize;
        uint32_t trialCount;
    };
    
 CompressionStats getStatsForField(uint32_t fieldIdx) const;
    void updateStatsForField(uint32_t fieldIdx, const CompressionStats& stats);

private:
    SamSelectorManager() = default;
    ~SamSelectorManager() = default;
    
    mutable std::mutex mutex_;
    std::map<uint32_t, CompressionModel> fieldSelections_;
    std::map<uint32_t, CompressionStats> fieldStats_;
};

} // namespace sam_compressor

#endif // SAM_SELECTOR_H