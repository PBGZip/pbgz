#pragma once

#include <vector>
#include <string>
#include <string.h>

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

    std::vector<char> getVersionAsArray();

    void updateReadDataLen(RoughIOBlock* blockPtr);
    
    void updateWriteDataLen(RoughIOBlock* blockPtr);

    void printHeadInfo(PbgzParameter& para);

    void printFileType(BlockType blockType);

    void printTailInfo(Timer costTime, bool isPrintRatio);

    void printBufferContent(uint8_t* buffer, uint32_t bufferLen) {
        char temp[2048 + 1] = {0};
        bufferLen = std::min<uint32_t>(bufferLen, 2048);
        memcpy(temp, (char*)buffer, bufferLen);
        temp[2048] = 0;
        LOG_DEBUG("%s", temp);
    }
private:
    void updateDataInfo();

    std::vector<std::pair<std::string, bool>> outfiles; 

    uint64_t totalReadLen;

    uint64_t totalWriteLen;
};

void pbgzExitProc(int errorCode, const char* errorMessage);

void coderLog(int logLevel, const char* logMessage);

int32_t powerof2Proximal(int32_t i);
