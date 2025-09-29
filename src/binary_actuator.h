#include "actuator.h"

class BinaryActuator : public Actuator {
public:
    BinaryActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr): Actuator(inPtr, outPtr) {}
    virtual ~BinaryActuator() {}
    int32_t decompress() override ;
    int32_t compress() override ;
};
