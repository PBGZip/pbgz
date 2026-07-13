/*
 * pbgz_stat.cpp - Source file for pbgz project
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

#include "pbgz_stat.h"
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <cstdio>

PbgzStat::PbgzStat() {
}

PbgzStat::~PbgzStat() {
    statUnits.clear();
}

int PbgzStat::init() {
    statUnits.clear();
    return 0;
}

void PbgzStat::registerStatUnit(uint16_t unitId, const std::string& unitName) {
    if (statUnits.find(unitId) == statUnits.end()) {
        statUnits[unitId] = StatUnit(unitId, unitName);
    }
}

void PbgzStat::addStatUnitMetric(uint16_t unitId, uint16_t metricId) {
    auto it = statUnits.find(unitId);
    if (it != statUnits.end()) {
        // 从全局变量中查找指标定义
        auto metricIt = gStatMetrics.find(metricId);
        if (metricIt != gStatMetrics.end()) {
            // 只保存metricId
            it->second.addMetricId(metricId);
        }
    }
}

void PbgzStat::setMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId, uint64_t value) {
    auto unitIt = statUnits.find(unitId);
    if (unitIt != statUnits.end()) {
        StatObject* obj = unitIt->second.getStatObject(objectId);
        if (obj) {
            obj->setMetricValue(metricId, value);
        }
    }
}

void PbgzStat::addMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId, uint64_t value) {
    auto unitIt = statUnits.find(unitId);
    if (unitIt != statUnits.end()) {
        StatObject* obj = unitIt->second.getStatObject(objectId);
        if (obj) {
            obj->addMetricValue(metricId, value);
        }
    }
}

uint64_t PbgzStat::getMetricValue(uint16_t unitId, uint16_t objectId, uint16_t metricId) const {
    auto unitIt = statUnits.find(unitId);
    if (unitIt != statUnits.end()) {
        const StatUnit& unit = unitIt->second;
        auto objIt = unit.getObjects().find(objectId);
        if (objIt != unit.getObjects().end()) {
            const StatObject& obj = objIt->second;
            
            // 从全局变量中查找指标定义
            auto metricIt = gStatMetrics.find(metricId);
            if (metricIt != gStatMetrics.end()) {
                const StatMetric& metric = metricIt->second;
                
                // 如果是计算类型的指标，需要动态计算
                if (metric.metricType == MetricType::CALCULATED && metric.calcFunc) {
                    std::map<uint16_t, uint64_t> allValues;
                    // 收集该对象所有指标的值
                    const auto& metricIds = unit.getMetricIds();
                    for (const auto& id : metricIds) {
                        allValues[id] = obj.getMetricValue(id);
                    }
                    return metric.calcFunc(allValues);
                } else {
                    // 直接返回存储的值
                    return obj.getMetricValue(metricId);
                }
            }
            return obj.getMetricValue(metricId);
        }
    }
    return 0;
}

void PbgzStat::addStatObject(uint16_t unitId, uint16_t objectId, const std::string& objectName) {
    auto it = statUnits.find(unitId);
    if (it != statUnits.end()) {
        it->second.addStatObject(objectId, objectName);
    }
}

std::map<uint16_t, uint64_t> PbgzStat::getUnitMetricValues(uint16_t unitId, uint16_t metricId) const {
    std::map<uint16_t, uint64_t> result;
    auto unitIt = statUnits.find(unitId);
    if (unitIt != statUnits.end()) {
        const StatUnit& unit = unitIt->second;
        for (const auto& objPair : unit.getObjects()) {
            // 对于计算类型的指标，需要通过getMetricValue获取计算后的值
            result[objPair.first] = getMetricValue(unitId, objPair.first, metricId);
        }
    }
    return result;
}

void PbgzStat::printStats() {
    fprintf(stderr, "=== Statistics ===\n");
    
    if (statUnits.empty()) {
        fprintf(stderr, "No statistics data available.\n");
        return;
    }
    
    // 遍历每个统计单元，按表格打印
    for (const auto& unitPair : statUnits) {
        uint16_t unitId = unitPair.first;
        const StatUnit& unit = unitPair.second;
        
        // 获取该单元的指标ID列表
        const std::vector<uint16_t>& metricIds = unit.getMetricIds();
        if (metricIds.empty()) {
            fprintf(stderr, "No metrics defined for unit %s.\n", unit.unitName.c_str());
            continue;
        }
        
        // 获取该单元的所有统计对象
        const auto& objects = unit.getObjects();
        if (objects.empty()) {
            fprintf(stderr, "No statistics objects for unit %s.\n", unit.unitName.c_str());
            continue;
        }
        
        // 计算列宽
        std::vector<int> colWidths(metricIds.size());
        int maxObjectNameWidth = 15; // 默认对象名称列宽
        
        // 计算对象名称的最大宽度
        for (const auto& objPair : objects) {
            const std::string& objName = objPair.second.objectName;
            if (objName.length() > static_cast<size_t>(maxObjectNameWidth)) {
                maxObjectNameWidth = static_cast<int>(objName.length());
            }
        }
        
        // 计算每列的宽度
        for (size_t i = 0; i < metricIds.size(); i++) {
            auto metricIt = gStatMetrics.find(metricIds[i]);
            if (metricIt != gStatMetrics.end()) {
                const std::string& metricName = metricIt->second.metricName;
                // 为压缩比添加额外宽度以容纳%符号
                if (metricIt->second.metricType == MetricType::CALCULATED) {
                    colWidths[i] = std::max(26, (int)metricName.length());
                } else {
                    colWidths[i] = std::max(20, (int)metricName.length());
                }
            }
        }
        
        // 对象名称列宽度
        const int nameWidth = maxObjectNameWidth;
        
        // 打印统计单元标题
        fprintf(stderr, "\n=== %s ===\n", unit.unitName.c_str());
        
        // 打印表头
        fprintf(stderr, "%-*s", nameWidth, "Object");
        for (size_t i = 0; i < metricIds.size(); i++) {
            auto metricIt = gStatMetrics.find(metricIds[i]);
            if (metricIt != gStatMetrics.end()) {
                fprintf(stderr, "%-*s", colWidths[i], metricIt->second.metricName.c_str());
            }
        }
        fprintf(stderr, "\n");
        
        // 打印分隔线
        fprintf(stderr, "%.*s", nameWidth, "----------------------------------------------------------------");
        for (int width : colWidths) {
            fprintf(stderr, "%.*s", width, "----------------------------------------------------------------");
        }
        fprintf(stderr, "\n");
        
        // 打印数据行
        for (const auto& objPair : objects) {
            const StatObject& obj = objPair.second;
            
            // 打印对象名称
            fprintf(stderr, "%-*s", nameWidth, obj.objectName.c_str());
            
            // 打印每个指标的值
            for (size_t i = 0; i < metricIds.size(); i++) {
                uint16_t metricId = metricIds[i];
                auto metricIt = gStatMetrics.find(metricId);
                if (metricIt == gStatMetrics.end()) continue;
                
                const StatMetric& metric = metricIt->second;
                uint64_t value = getMetricValue(unitId, objPair.first, metricId);
                
                // 格式化压缩比作为百分比
                if (metric.metricType == MetricType::CALCULATED) {
                    char buffer[32];
                    snprintf(buffer, sizeof(buffer), "%.2f%%", (double)value / 100.0);
                    fprintf(stderr, "%-*s", colWidths[i], buffer);
                } else {
                    fprintf(stderr, "%-*lu", colWidths[i], value);
                }
            }
            fprintf(stderr, "\n");
        }
        
        // 打印底部线
        fprintf(stderr, "%.*s", nameWidth, "================================================================");
        for (int width : colWidths) {
            fprintf(stderr, "%.*s", width, "================================================================");
        }
        fprintf(stderr, "\n");
    }
}

FastqStat::FastqStat() : PbgzStat() {
}

FastqStat::~FastqStat() {
}

int FastqStat::init() {
    PbgzStat::init();
    
    // 注册统计单元（Compressed Size）
    registerStatUnit(StatUnitIds::COMPRESSION_RATIO, "Compressed Size");
    
    // 添加统计指标（通过指标ID）
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::ORIGINAL_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSED_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSION_RATIO);
    
    // 添加统计对象
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_ID, StatObjectNames::Fastq::ID);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_BASE, StatObjectNames::Fastq::BASE);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_COMMENT, StatObjectNames::Fastq::COMMENT);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::FASTQ_QUALITY, StatObjectNames::Fastq::QUALITY);
    
    return 0;
}

SamStat::SamStat() : PbgzStat() {
}

SamStat::~SamStat() {
}

int SamStat::init() {
    PbgzStat::init();
    
    // 注册统计单元（Compressed Size）
    registerStatUnit(StatUnitIds::COMPRESSION_RATIO, "Compressed Size");
    
    // 添加统计指标（通过指标ID）
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::ORIGINAL_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSED_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSION_RATIO);
    
    // 添加统计对象
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QNAME, StatObjectNames::Sam::QNAME);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_FLAG, StatObjectNames::Sam::FLAG);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNAME, StatObjectNames::Sam::RNAME);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_POS, StatObjectNames::Sam::POS);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_MAPQ, StatObjectNames::Sam::MAPQ);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_CIGAR, StatObjectNames::Sam::CIGAR);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_RNEXT, StatObjectNames::Sam::RNEXT);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_PNEXT, StatObjectNames::Sam::PNEXT);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_TLEN, StatObjectNames::Sam::TLEN);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_SEQ, StatObjectNames::Sam::SEQ);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_QUAL, StatObjectNames::Sam::QUAL);
    addStatObject(StatUnitIds::COMPRESSION_RATIO, StatObjectId::SAM_OPTION, StatObjectNames::Sam::OPTION);
    
    return 0;
}