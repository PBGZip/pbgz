/*
 * pbgz_stat.h - Header file for pbgz project
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

#include <cstdint>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <memory>
#include <functional>

enum class MetricType {
    VALUE,           // 设置值
    ACCUMULATED,     // 累加值
    CALCULATED       // 计算值 (基于其他指标)
};

namespace StatMetricIds {
    const uint16_t ORIGINAL_SIZE = 1;
    const uint16_t COMPRESSED_SIZE = 2;
    const uint16_t COMPRESSION_RATIO = 3;
}

struct StatMetric {
    uint16_t metricId;
    std::string metricName;
    MetricType metricType;
    std::function<uint64_t(const std::map<uint16_t, uint64_t>&)> calcFunc;
    
    StatMetric() : metricId(0), metricType(MetricType::VALUE) {}
    StatMetric(uint16_t id, const std::string& name, MetricType type, 
               std::function<uint64_t(const std::map<uint16_t, uint64_t>&)> func = nullptr)
        : metricId(id), metricName(name), metricType(type), calcFunc(func) {}
};

// 全局统计指标定义：指标ID对应的名称、类型、计算公式
const std::map<uint16_t, StatMetric> gStatMetrics = {
    {StatMetricIds::ORIGINAL_SIZE, StatMetric(StatMetricIds::ORIGINAL_SIZE, "Original Size", MetricType::ACCUMULATED)},
    {StatMetricIds::COMPRESSED_SIZE, StatMetric(StatMetricIds::COMPRESSED_SIZE, "Compressed Size", MetricType::ACCUMULATED)},
    {StatMetricIds::COMPRESSION_RATIO, StatMetric(StatMetricIds::COMPRESSION_RATIO, "Compression Ratio", MetricType::CALCULATED,
        [](const std::map<uint16_t, uint64_t>& metrics) -> uint64_t {
            auto original = metrics.find(StatMetricIds::ORIGINAL_SIZE);
            auto compressed = metrics.find(StatMetricIds::COMPRESSED_SIZE);
            if (original != metrics.end() && compressed != metrics.end() && original->second > 0) {
                double ratio = (double)compressed->second / original->second;
                return (uint64_t)(ratio * 10000);
            }
            return 0;
        })}
};

class StatObject {
public:
    uint16_t objectId;
    std::string objectName;
    std::map<uint16_t, uint64_t> metricValues;  // metricId -> 统计值
    
    StatObject(uint16_t id = 0, const std::string& name = "") : objectId(id), objectName(name) {}
    
    void setMetricValue(uint16_t metricId, uint64_t value) {
        metricValues[metricId] = value;
    }
    
    void addMetricValue(uint16_t metricId, uint64_t value) {
        metricValues[metricId] += value;
    }
    
    uint64_t getMetricValue(uint16_t metricId) const {
        auto it = metricValues.find(metricId);
        return it != metricValues.end() ? it->second : 0;
    }
};

class StatUnit {
public:
    uint16_t unitId;
    std::string unitName;
    std::vector<uint16_t> metrics;  // 该统计单元包含的统计指标ID列表
    std::map<uint16_t, StatObject> objects;  // objectId -> StatObject
    
    StatUnit(uint16_t id = 0, const std::string& name = "") : unitId(id), unitName(name) {}
    
    void addMetricId(uint16_t metricId) {
        metrics.push_back(metricId);
    }
    
    void addStatObject(uint16_t objectId, const std::string& objectName) {
        if (objects.find(objectId) == objects.end()) {
            objects[objectId] = StatObject(objectId, objectName);
        }
    }
    
    StatObject* getStatObject(uint16_t objectId) {
        auto it = objects.find(objectId);
        if (it != objects.end()) {
            return &it->second;
        }
        return nullptr;
    }
    
    const std::vector<uint16_t>& getMetricIds() const {
        return metrics;
    }
    
    const std::map<uint16_t, StatObject>& getObjects() const {
        return objects;
    }
};

class PbgzStat {
public:
    PbgzStat();
    virtual ~PbgzStat();
    
    virtual int init();
    
    // 统计单元管理
    void registerStatUnit(uint16_t unitId, const std::string& unitName);
    void addStatUnitMetric(uint16_t unitId, uint16_t metricId);
    
    // 统计对象管理
    void addStatObject(uint16_t unitId, uint16_t objectId, const std::string& objectName);
    
    // 统计值操作
    void setMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId, uint64_t value);
    void addMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId, uint64_t value);
    uint64_t getMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId) const;
    
    // 获取指定统计单元中所有对象的某个指标值
    std::map<uint16_t, uint64_t> getUnitMetricValues(uint16_t unitId, uint16_t metricId) const;
    
    // 打印统计结果
    void printStats();
    
protected:
    // 三层结构：
    // 第一层：统计单元（StatUnit），包含指标和对象
    std::map<uint16_t, StatUnit> statUnits;
};

namespace StatUnitIds{
    const uint16_t COMPRESSION_RATIO = 1;
}; // namespace StatUnitIds

namespace StatObjectId {
    const uint16_t FASTQ_ID = 1;
    const uint16_t FASTQ_BASE = 2;
    const uint16_t FASTQ_COMMENT = 3;
    const uint16_t FASTQ_QUALITY = 4;
    
    const uint16_t SAM_QNAME = 5;
    const uint16_t SAM_FLAG = 6;
    const uint16_t SAM_RNAME = 7;
    const uint16_t SAM_POS = 8;
    const uint16_t SAM_MAPQ = 9;
    const uint16_t SAM_CIGAR = 10;
    const uint16_t SAM_RNEXT = 11;
    const uint16_t SAM_PNEXT = 12;
    const uint16_t SAM_TLEN = 13;
    const uint16_t SAM_SEQ = 14;
    const uint16_t SAM_QUAL = 15;
}

namespace StatObjectNames {
    namespace Fastq {
        const std::string ID = "ID";
        const std::string BASE = "Base";
        const std::string COMMENT = "Comment";
        const std::string QUALITY = "Quality";
    }
    
    namespace Sam {
        const std::string QNAME = "QNAME";
        const std::string FLAG = "FLAG";
        const std::string RNAME = "RNAME";
        const std::string POS = "POS";
        const std::string MAPQ = "MAPQ";
        const std::string CIGAR = "CIGAR";
        const std::string RNEXT = "RNEXT";
        const std::string PNEXT = "PNEXT";
        const std::string TLEN = "TLEN";
        const std::string SEQ = "SEQ";
        const std::string QUAL = "QUAL";
    }
}

class FastqStat : public PbgzStat {
public:
    FastqStat();
    ~FastqStat() override;
    
    int init() override;
    
private:
};

class SamStat : public PbgzStat {
public:
     SamStat();
    ~SamStat() override;
    
    int init() override;
    
private:
};