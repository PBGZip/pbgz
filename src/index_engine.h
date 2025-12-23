#pragma once

#include "pbgz_types.h"
#include "decompress_engine.h"

class IndexEngine : public DecompressEngine {
public:
    IndexEngine(PbgzParameter& para) : DecompressEngine(para) {
    }

    ~IndexEngine() {
    }

private:
    
};
