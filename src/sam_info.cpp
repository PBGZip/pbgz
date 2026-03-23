#include "sam_info.h"

void SamInfo::addChrNameIndex(const std::string& chrName) {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    chrNameIndex[chrName] = chrNameIndex.size();
}

uint16_t SamInfo::getChrNameIndex(const std::string& chrName) const {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    auto it = chrNameIndex.find(chrName);
    if (it != chrNameIndex.end()) {
        return it->second;
    }
    // 返回一个表示未找到的值，这里使用65535（uint16_t的最大值）
    return 65535;
}

bool SamInfo::hasChrName(const std::string& chrName) const {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    return chrNameIndex.find(chrName) != chrNameIndex.end();
}

void SamInfo::addChromosomeInfo(const std::string& name, uint32_t length) {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    
    // Check if chromosome with this name already exists
    auto it = chrNameIndex.find(name);
    if (it != chrNameIndex.end()) {
        // Chromosome already exists, update its info if needed
        for (auto& chrInfo : chromosomeInfoList) {
            if (chrInfo.name == name) {
                chrInfo.length = length; // Update length if different
                // 重新计算所有染色体的position偏移量
                calculateChromosomePositions();
                return;
            }
        }
    } else {
        // New chromosome, get next available ID and add it
        uint16_t newId = chrIdCounter++;
        ChromosomeInfo chrInfo;
        chrInfo.id = newId;
        chrInfo.name = name;
        chrInfo.length = length;
        chromosomeInfoList.push_back(chrInfo);
        
        // Also add to name-to-id map for quick lookup
        chrNameIndex[name] = newId;
        
        // 计算所有染色体的position偏移量
        calculateChromosomePositions();
    }
}

const ChromosomeInfo& SamInfo::getChromosomeInfo(uint16_t id) const {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    static ChromosomeInfo emptyInfo; // 返回空信息的默认对象
    if (id < chromosomeInfoList.size()) {
        return chromosomeInfoList[id];
    }
    return emptyInfo;
}

const std::vector<ChromosomeInfo>& SamInfo::getAllChromosomeInfo() const {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    return chromosomeInfoList;
}

void SamInfo::clearChromosomeInfo() {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    chromosomeInfoList.clear();
    chrNameIndex.clear();
}

uint16_t SamInfo::getNextChrId() {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    return chrIdCounter++;
}

void SamInfo::resetChrIdCounter() {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    chrIdCounter = 0;
}

void SamInfo::calculateChromosomePositions() {
    // 按照ID从小到大排序染色体信息列表
    std::sort(chromosomeInfoList.begin(), chromosomeInfoList.end(), 
              [](const ChromosomeInfo& a, const ChromosomeInfo& b) {
                  return a.id < b.id;
              });
    
    // 计算每个染色体的position偏移量
    uint64_t accumulatedPosition = 0;
    for (auto& chrInfo : chromosomeInfoList) {
        chrInfo.position = accumulatedPosition;
        accumulatedPosition += chrInfo.length;
    }

    return;
}

int64_t SamInfo::getPosistionByIndex(uint32_t index) {
    if (chromosomeInfoList.size() < index) {
        return -1;
    }

    return chromosomeInfoList[index].position;
}

