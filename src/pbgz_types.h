#pragma once

#include <string>
#include <unistd.h>
#include <chrono>

typedef struct PbgzParameter{
    std::string inputFile;
    bool isDecToGZ;   // Whether to decompress to GZ format
    std::string outputFile;  // 
    std::string outputDir;   // 
    bool isOverwriteOutFile;         // Whether to force overwrite output file
    std::string referenceGenic;      // Reference genome
    bool isUnpackRef;  
    uint32_t threadNum; 
    uint8_t compressLevel;
    bool isRemoveOriginFile;    // Whether to remove source file
    uint8_t logLevel;
    std::string logFile;
    std::string refeGenePos;
    bool isMakeIndex;
    bool showHelp;

    PbgzParameter() {
        isDecToGZ = false;
        isOverwriteOutFile = false;
        isUnpackRef = false;
        threadNum = sysconf(_SC_NPROCESSORS_ONLN);
        compressLevel = 5;
        isRemoveOriginFile = false;
        logLevel = 6;
        isMakeIndex = false;
        showHelp = false;
    }
} PbgzParameter;


const std::string STDIN = "/dev/stdin";
const std::string STDOUT = "/dev/stdout";

const uint8_t PBGZ_VERSION_MAJOR = 2; // Major version of PBGZ file format
const uint8_t PBGZ_VERSION_MINOR = 1; // Minor version of PBGZ file format
const uint8_t PBGZ_VERSION_PATCH = 0; // Patch version

const int32_t GENE3_MAX_BASE = 1048760; // Maximum bases per FASTQ generation 3 record
const int32_t GENE2_MAX_BASE = 2048;  // Maximum bases per FASTQ generation 2 record
