/*
 * pbgz_types.h - Header file for pbgz project
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
#include <unistd.h>
#include <chrono>

typedef struct PbgzParameter{
    std::string inputFile;
    bool isDecToGZ;   // Whether to decompress to GZ format
    bool isDecToBam;  // Whether to decompress SAM to BAM format
    std::string outputFile;  // 
    std::string outputDir;   // 
    bool isOverwriteOutFile;         // Whether to force overwrite output file
    std::string referenceGenic;      // Reference genome
    std::string niIndexFile;         // Prebuilt index; it is only an acceleration aid, the reference genome is always authoritative via referenceGenic
    bool isUnpackRef;  
    uint32_t threadNum; 
    uint8_t compressLevel;
    bool isRemoveOriginFile;    // Whether to remove source file
    uint8_t logLevel;
    std::string logFile;
    std::string refeGenePos;
    bool isMakeIndex;
    bool showHelp;
    bool showStat;
    bool verbose;

    PbgzParameter() {
        isDecToGZ = false;
        isDecToBam = false;
        isOverwriteOutFile = false;
        isUnpackRef = false;
        threadNum = sysconf(_SC_NPROCESSORS_ONLN);
        compressLevel = 5;
        isRemoveOriginFile = false;
        logLevel = 6;
        isMakeIndex = false;
        showHelp = false;
        showStat = false;
        verbose = false;
    }
} PbgzParameter;


const std::string STDIN = "/dev/stdin";
const std::string STDOUT = "/dev/stdout";

const uint8_t PBGZ_VERSION_MAJOR = 2; // Major version of PBGZ file format
const uint8_t PBGZ_VERSION_MINOR = 2; // Minor version of PBGZ file format
const uint8_t PBGZ_VERSION_PATCH = 0; // Patch version

const int32_t GENE3_MAX_BASE = 1048760; // Maximum bases per FASTQ generation 3 record
const int32_t GENE2_MAX_BASE = 2048;  // Maximum bases per FASTQ generation 2 record
