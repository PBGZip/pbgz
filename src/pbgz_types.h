#pragma once

#include <string>
#include <unistd.h>
#include <chrono>

typedef struct PbgzParameter{
    std::string inputFile;
    bool showHelp;
    bool showVersion;
    bool isDecompress;   // 标记是解压还是压缩， 默认时压缩
    bool isDecToGZ;   // 是否解压到GZ格式
    std::string outputFile;  // 
    std::string outputDir;   // 
    bool isOverwriteOutFile;         // 是否强制覆盖输出文件
    std::string referenceGenic;      // 参考基因
    bool isUnpackRef;  
    uint32_t threadNum; 
    uint8_t compressLevel;
    bool isRemoveOriginFile;    //是否删除源文件
    uint8_t logLevel;
    std::string logFile;

    PbgzParameter() {
        showHelp = false;
        showVersion = false;
        isDecompress = false;
        isDecToGZ = false;
        isOverwriteOutFile = false;
        isUnpackRef = false;
        threadNum = sysconf(_SC_NPROCESSORS_ONLN);
        compressLevel = 2;
        isRemoveOriginFile = false;
        logLevel = 6;
    }
} PbgzParameter;


const std::string STDIN = "/dev/stdin";
const std::string STDOUT = "/dev/stdout";

const uint8_t PBGZ_VERSION_MAJOR = 2; // Major version of PBGZ file format
const uint8_t PBGZ_VERSION_MINOR = 0; // Minor version of PBGZ file format
const uint8_t PBGZ_VERSION_PATCH = 0; // Patch version

const int32_t GENE3_MAX_BASE = 1048760; // Maximum bases per FASTQ generation 3 record
const int32_t GENE2_MAX_BASE = 2048;  // Maximum bases per FASTQ generation 2 record

