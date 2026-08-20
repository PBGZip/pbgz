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
        // Look up the metric definition in the global table
        auto metricIt = gStatMetrics.find(metricId);
        if (metricIt != gStatMetrics.end()) {
            // Only store the metricId
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
            
            // Look up the metric definition in the global table
            auto metricIt = gStatMetrics.find(metricId);
            if (metricIt != gStatMetrics.end()) {
                const StatMetric& metric = metricIt->second;
                
                // Calculated-type metrics need to be computed dynamically
                if (metric.metricType == MetricType::CALCULATED && metric.calcFunc) {
                    std::map<uint16_t, uint64_t> allValues;
                    // Collect the values of all metrics of this object
                    const auto& metricIds = unit.getMetricIds();
                    for (const auto& id : metricIds) {
                        allValues[id] = obj.getMetricValue(id);
                    }
                    return metric.calcFunc(allValues);
                } else {
                    // Return the stored value directly
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
            // For calculated-type metrics, retrieve the computed value through getMetricValue
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
    
    // Iterate over each statistic unit and print it as a table
    for (const auto& unitPair : statUnits) {
        uint16_t unitId = unitPair.first;
        const StatUnit& unit = unitPair.second;
        
        // Get the list of metric IDs of this unit
        const std::vector<uint16_t>& metricIds = unit.getMetricIds();
        if (metricIds.empty()) {
            fprintf(stderr, "No metrics defined for unit %s.\n", unit.unitName.c_str());
            continue;
        }
        
        // Get all statistic objects of this unit
        const auto& objects = unit.getObjects();
        if (objects.empty()) {
            fprintf(stderr, "No statistics objects for unit %s.\n", unit.unitName.c_str());
            continue;
        }
        
        // Compute column widths
        std::vector<int> colWidths(metricIds.size());
        int maxObjectNameWidth = 15; // Default object name column width
        
        // Compute the maximum width of the object names
        for (const auto& objPair : objects) {
            const std::string& objName = objPair.second.objectName;
            if (objName.length() > static_cast<size_t>(maxObjectNameWidth)) {
                maxObjectNameWidth = static_cast<int>(objName.length());
            }
        }
        
        // Compute the width of each column
        for (size_t i = 0; i < metricIds.size(); i++) {
            auto metricIt = gStatMetrics.find(metricIds[i]);
            if (metricIt != gStatMetrics.end()) {
                const std::string& metricName = metricIt->second.metricName;
                // Add extra width for the compression ratio to accommodate the % sign
                if (metricIt->second.metricType == MetricType::CALCULATED) {
                    colWidths[i] = std::max(26, (int)metricName.length());
                } else {
                    colWidths[i] = std::max(20, (int)metricName.length());
                }
            }
        }
        
        // Object name column width
        const int nameWidth = maxObjectNameWidth;
        
        // Print the statistic unit title
        fprintf(stderr, "\n=== %s ===\n", unit.unitName.c_str());
        
        // Print the header
        fprintf(stderr, "%-*s", nameWidth, "Object");
        for (size_t i = 0; i < metricIds.size(); i++) {
            auto metricIt = gStatMetrics.find(metricIds[i]);
            if (metricIt != gStatMetrics.end()) {
                fprintf(stderr, "%-*s", colWidths[i], metricIt->second.metricName.c_str());
            }
        }
        fprintf(stderr, "\n");
        
        // Print the separator line
        fprintf(stderr, "%.*s", nameWidth, "----------------------------------------------------------------");
        for (int width : colWidths) {
            fprintf(stderr, "%.*s", width, "----------------------------------------------------------------");
        }
        fprintf(stderr, "\n");
        
        // Print the data rows
        for (const auto& objPair : objects) {
            const StatObject& obj = objPair.second;
            
            // Print the object name
            fprintf(stderr, "%-*s", nameWidth, obj.objectName.c_str());
            
            // Print the value of each metric
            for (size_t i = 0; i < metricIds.size(); i++) {
                uint16_t metricId = metricIds[i];
                auto metricIt = gStatMetrics.find(metricId);
                if (metricIt == gStatMetrics.end()) continue;
                
                const StatMetric& metric = metricIt->second;
                uint64_t value = getMetricValue(unitId, objPair.first, metricId);
                
                // Format the compression ratio as a percentage
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
        
        // Print the bottom line
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
    
    // Register the statistic unit (Compressed Size)
    registerStatUnit(StatUnitIds::COMPRESSION_RATIO, "Compressed Size");
    
    // Add the statistic metrics (by metric ID)
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::ORIGINAL_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSED_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSION_RATIO);
    
    // Add the statistic objects
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
    
    // Register the statistic unit (Compressed Size)
    registerStatUnit(StatUnitIds::COMPRESSION_RATIO, "Compressed Size");
    
    // Add the statistic metrics (by metric ID)
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::ORIGINAL_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSED_SIZE);
    addStatUnitMetric(StatUnitIds::COMPRESSION_RATIO, StatMetricIds::COMPRESSION_RATIO);
    
    // Add the statistic objects
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