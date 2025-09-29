#pragma once

#include <json/json.h>

#include "io_block.h"

class Actuator {
public:
    virtual int32_t decompress() = 0;
    virtual int32_t compress() = 0;

    Actuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr): inBlockPtr(inPtr), outBlockPtr(outPtr) {};
    
    virtual ~Actuator() {
        inBlockPtr = nullptr;
        outBlockPtr = nullptr;
    }

protected: 
    RoughIOBlock* inBlockPtr;
    RoughIOBlock* outBlockPtr;  
    Json::Value meta;
};