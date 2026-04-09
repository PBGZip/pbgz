#include "sort_engine.h"

SortEngine::~SortEngine() {

}

void SortEngine::printHeadInfo() {

}

void SortEngine::printTailInfo(Timer&) {

}

int32_t SortEngine::startEnginePreProc() {
    return 0;
}



int32_t SortEngine::startWorkPreProc() {
    return 0;
}

int32_t SortEngine::startReadPreProc() {
    return 0;
}

int32_t SortEngine::startEnginePostProc() {
    return 0;
}

BlockReader* SortEngine::createBlockReader() {
    return nullptr;
}

BlockWriter* SortEngine::createBlockWriter() {
    return nullptr;
}

void SortEngine::releaseBlockReader(BlockReader* &blockReader) {

}

void SortEngine::releaseBlockWriter(BlockWriter* &blockWriter) {

}

Actuator* SortEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    return nullptr;
}

int32_t SortEngine::actuatorProc(Actuator*, RoughIOBlock*, RoughIOBlock*) {
    return 0;
}
