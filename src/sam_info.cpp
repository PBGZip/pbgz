/*
 * sam_info.cpp - Implementation of sam_info
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

#include "sam_info.h"
#include "log/logger.h"

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
    // Return a value indicating not found, using 65535 (maximum value of uint16_t)
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
                // Recalculate position offsets for all chromosomes
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
        
        // Calculate position offsets for all chromosomes
        calculateChromosomePositions();
    }
}

const ChromosomeInfo& SamInfo::getChromosomeInfo(uint16_t id) const {
    std::lock_guard<std::mutex> lock(chrNameMutex);
    static ChromosomeInfo emptyInfo; // Default object returning empty information
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
    // Sort chromosome information list by ID in ascending order
    std::sort(chromosomeInfoList.begin(), chromosomeInfoList.end(), 
              [](const ChromosomeInfo& a, const ChromosomeInfo& b) {
                  return a.id < b.id;
              });
    
    // Calculate position offset for each chromosome
    uint64_t accumulatedPosition = 0;
    for (auto& chrInfo : chromosomeInfoList) {
        chrInfo.position = accumulatedPosition;
        accumulatedPosition += chrInfo.length;
    }
    return;
}

int64_t SamInfo::getPositionByIndex(uint32_t index) {
    if (chromosomeInfoList.size() < index) {
        return -1;
    }
    return chromosomeInfoList[index].position;
}


int32_t SamUtil::parseChromosomeInfo(const std::string& sqLine) {
    // Parse @SQ line, format: @SQ SN:chr_name LN:length [other optional fields]
    std::string chrName;
    uint32_t chrLength = 0;
    
    // Parse SN field (chromosome name)
    size_t snPos = sqLine.find("SN:");
    if (snPos != std::string::npos) {
        snPos += 3;
        size_t snEnd = sqLine.find("\t", snPos);
        if (snEnd == std::string::npos) {
            snEnd = sqLine.length();
        }
        chrName = sqLine.substr(snPos, snEnd - snPos);
    } else {
        LOG_ERROR("Cannot find SN field in @SQ line: %s", sqLine.c_str());
        return -1;
    }
    
    // Parse LN field (chromosome length)
    size_t lnPos = sqLine.find("LN:");
    if (lnPos != std::string::npos) {
        lnPos += 3;
        size_t lnEnd = sqLine.find("\t", lnPos);
        if (lnEnd == std::string::npos) {
            lnEnd = sqLine.length();
        }
        std::string lengthStr = sqLine.substr(lnPos, lnEnd - lnPos);
        try {
            chrLength = std::stoul(lengthStr);
        } catch (const std::exception& e) {
            LOG_ERROR("Invalid chromosome length in @SQ line: %s", sqLine.c_str());
            return -1;
        }
    } else {
        LOG_ERROR("Cannot find LN field in @SQ line: %s", sqLine.c_str());
        return -1;
    }
    
    // Add chromosome information to SamInfo (automatically gets ID internally)
    SamInfo::getInstance().addChromosomeInfo(chrName, chrLength);
    LOG_INFO("Parsed chromosome info: Name=%s, Length=%u", chrName.c_str(), chrLength);
    return 0;
}
