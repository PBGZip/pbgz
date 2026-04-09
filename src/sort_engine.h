#pragma once

#include "pbgz_engine.h"
#include "utils/timer.h"

class SortEngine : public PbgzEngine {

public:
    using PbgzEngine::PbgzEngine;

    virtual ~SortEngine();

protected:
    virtual void printHeadInfo() override;

    virtual void printTailInfo(Timer&) override;

    virtual int32_t startEnginePreProc() override;

    virtual int32_t startWorkPreProc() override;


    virtual int32_t startReadPreProc() override;


    virtual int32_t startEnginePostProc() override;

    virtual BlockReader* createBlockReader() override;

    virtual BlockWriter* createBlockWriter() override;

    virtual void releaseBlockReader(BlockReader* &blockReader) override;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) override;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) override;

    virtual int32_t actuatorProc(Actuator*, RoughIOBlock*, RoughIOBlock*) override;
};