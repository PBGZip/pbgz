#pragma once

#include <stdint.h>
#include <map>



class PbgzIndex {

public:

    void setBlockPosition(uint32_t blockId, int64_t position) {
        blockPositions[blockId] = position;
    }

    int32_t loadIndex() {
        return 0;
    }

    int32_t initIndex() {
        return 0;
    }

    int32_t packIndex() {
        return 0;
    }

private:
    int64_t refeBeginPos;
    int64_t refeEndPos;
    
    std::map<uint32_t, int64_t> blockPositions;

};