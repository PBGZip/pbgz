#pragma once

#include "pbgz_engine.h"
#include "actuator.h"


class DecompressEngine : public PbgzEngine {
public:
    DecompressEngine(PbgzParameter& para) : PbgzEngine(para) {

    }

protected:
    BlockReader* createBlockReader();

    BlockWriter* createBlockWriter();

    virtual void readBlocks(BlockReader* blockReader);

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*);

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*);

private:
    bool initRefGene(PbgzBlockReader* blockReader);

    bool initRefeIndex(PbgzBlockReader* blockReader);

    void readBlockByPostition(BlockReader* blockReader);

    bool unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta);
};