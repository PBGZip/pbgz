#pragma once

#include<vector>
#include<string>

#include "coder/coder.h"
#include "log/logger.h"
#include "io_block.h"
#include "pbgz_types.h"
#include "utils/timer.h"

class PbgzManager {
public:
    PbgzManager() {
        totalReadLen = 0;
        totalWriteLen = 0;
    }

    PbgzManager& operator=(const PbgzManager& item) = delete;
    PbgzManager(const PbgzManager& item) = delete;

    ~PbgzManager() {

    }

    void exitProc(int errorCode, const char* errorMessage);

    static PbgzManager& getInstance();

    std::string getVersion();

    void updateReadDataLen(RoughIOBlock* blockPtr);
    
    void updateWriteDataLen(RoughIOBlock* blockPtr);

    void printHeadInfo(PbgzParameter& para);

    void printFileType(BlockType blockType);

    void printTailInfo(Timer costTime, PbgzParameter& para);
private:
    void updateDataInfo();

    std::vector<std::pair<std::string, bool>> outfiles; 

    uint64_t totalReadLen;

    uint64_t totalWriteLen;
};

void pbgzExitProc(int errorCode, const char* errorMessage);

void coderLog(int logLevel, const char* logMessage);
