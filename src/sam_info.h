/*
 * sam_info.h - Header file for sam_info
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

#include <string>
#include <stdint.h>
#include <map>
#include <mutex>
#include <vector>
#include <algorithm>

// Chromosome information structure
struct ChromosomeInfo {
    uint16_t id;        // Chromosome ID
    std::string name;   // Chromosome name
    uint32_t length;    // Chromosome length
    uint64_t position;  // Position offset of chromosome relative to start position
    
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

     // Add mapping of chromosome name and index
    void addChrNameIndex(const std::string& chrName);
    
     // Get index corresponding to chromosome name
    uint16_t getChrNameIndex(const std::string& chrName) const;
    
     // Check if chromosome name exists
    bool hasChrName(const std::string& chrName) const;
    
     // Add chromosome information (automatically gets ID, only needs name and length)
    void addChromosomeInfo(const std::string& name, uint32_t length);
    
     // Get chromosome information
    const ChromosomeInfo& getChromosomeInfo(uint16_t id) const;
    
     // Get all chromosome information
    const std::vector<ChromosomeInfo>& getAllChromosomeInfo() const;
    
     // Clear chromosome information
    void clearChromosomeInfo();
    
     // Get next chromosome ID (global counter)
    uint16_t getNextChrId();
    
     // Reset chromosome ID counter
    void resetChrIdCounter();
    
     // Calculate position offsets for all chromosomes
    void calculateChromosomePositions();

    int64_t getPositionByIndex(uint32_t index);

private:
    mutable std::mutex chrNameMutex;
    std::map<std::string, uint16_t> chrNameIndex;
    std::vector<ChromosomeInfo> chromosomeInfoList; // Chromosome information list
    uint16_t chrIdCounter; // Global chromosome ID counter
};

namespace SamUtil {
    int32_t parseChromosomeInfo(const std::string& sqLine);

    /*
     * Pre-parse the header at the head of a block into the chromosome table.
     *
     * The header is the one part of a SAM file that has to be understood before any block can be
     * classified by its RNAME, since the names come from it (SamInfo::getChrNameIndex, used by
     * SamCodecActuator::scanDataLine); a block whose RNAME cannot be resolved is downgraded to
     * BINARY. Until now the table was only ever built by the first block's pre-analysis, on a
     * coder thread, which is why the pipeline has to serialise that block against all the others
     * (see the first-block serialization in PbgzEngine::startWorkTask).
     *
     * The reader keeps the header lines in a block of their own (SamBlockReader::readBlock with
     * splitSamHeader; BamBlockReader emits the BAM header the same way) and the reader thread sees
     * that block before anything is queued, so the table can be built there, ahead of every coder
     * thread - and no data block then races an unbuilt table. Blocks that carry their header merged
     * into the first data block go through here too.
     *
     * Only @SQ is consumed: it is the part every later stage reads. The other header kinds have no
     * engine-level effect - the UR: / CL: reference checks live in the actuator, which is the side
     * that holds the loaded reference. The pass is idempotent against the actuator's own one
     * (SamInfo keys by name), so whichever runs first fills the table and the other only confirms.
     *
     * `npos` holds the newline positions of the block, as the reader recorded them. Returns the
     * number of header lines read, 0 when the block does not open with one.
     */
    int32_t prefetchHeaderChromosomes(const uint8_t* buffer, int64_t dataLen,
                                      const std::vector<size_t>& npos);
}