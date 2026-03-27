#pragma once

#include <string>
#include <stdint.h>
#include <map>
#include <mutex>
#include <vector>
#include <algorithm>

// 染色体信息结构体
struct ChromosomeInfo {
    uint16_t id;        // 染色体ID
    std::string name;   // 染色体名称
    uint32_t length;    // 染色体长度
    uint64_t position;  // 染色体相对于起始位置的偏移
    
    ChromosomeInfo() : id(0), length(0), position(0) {}
    ChromosomeInfo(uint16_t id, const std::string& name, uint32_t length) 
        : id(id), name(name), length(length), position(0) {}
};

class SamInfo {

public:
    SamInfo() {
        chrIdCounter = 0; 
    }

    static SamInfo& getInstance() {
        static SamInfo instance;
        return instance;
    }

    // 添加染色体名称和索引的映射
    void addChrNameIndex(const std::string& chrName);
    
    // 获取染色体名称对应的索引
    uint16_t getChrNameIndex(const std::string& chrName) const;
    
    // 检查染色体名称是否存在
    bool hasChrName(const std::string& chrName) const;
    
    // 添加染色体信息（自动获取ID，只需要名称和长度）
    void addChromosomeInfo(const std::string& name, uint32_t length);
    
    // 获取染色体信息
    const ChromosomeInfo& getChromosomeInfo(uint16_t id) const;
    
    // 获取所有染色体信息
    const std::vector<ChromosomeInfo>& getAllChromosomeInfo() const;
    
    // 清空染色体信息
    void clearChromosomeInfo();
    
    // 获取下一个染色体ID（全局计数器）
    uint16_t getNextChrId();
    
    // 重置染色体ID计数器
    void resetChrIdCounter();
    
    // 计算所有染色体的position偏移量
    void calculateChromosomePositions();

    int64_t getPosistionByIndex(uint32_t index);

private:
    mutable std::mutex chrNameMutex;
    std::map<std::string, uint16_t> chrNameIndex;
    std::vector<ChromosomeInfo> chromosomeInfoList; // 染色体信息列表
    uint16_t chrIdCounter; // 全局染色体ID计数器
};
