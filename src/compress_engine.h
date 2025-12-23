#pragma once

#include "pbgz_engine.h"
#include "actuator.h"
#include <mutex>

class CompressEngine : public PbgzEngine {
public:
    CompressEngine(PbgzParameter& para) : PbgzEngine(para) {
        
    }

    virtual ~CompressEngine();

    virtual int32_t init();

protected:
    virtual BlockReader* createBlockReader() override;
    
    virtual BlockWriter* createBlockWriter() override;

    virtual int32_t engineStartPreProc() override;

    virtual int32_t engineStartAfterProc() override;

    virtual Actuator* actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr);

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock* outBlockPtr);

    virtual bool isPrintRatio() {
        return true;
    }

    virtual void setDataBlockPosition(uint32_t blockId);

    void startWriteIndexTask();

    virtual int32_t beforeCoderProc();

private:
    bool initReference();

    int64_t packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef = true);

    uint32_t calcPackRefeBlockSize();

private:
    std::map<uint32_t, std::vector<int64_t>> blockRefePos;
    std::map<int64_t, uint32_t> blockRefeIndex;
    BlockingQueue<RoughIOBlock*> indexBlockQueue;
    BlockingQueue<RoughIOBlock*> freeIndexBlockQueue;  
};
