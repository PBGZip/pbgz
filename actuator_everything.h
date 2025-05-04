#ifndef _ACTUATOR_EVERYTHING_H_
#define _ACTUATOR_EVERYTHING_H_

#include "block.h"

/* 通用传动装置 */

class actuator_everything
{
public:
    actuator_everything(const block_rough_ptr bptr_in, block_rough_ptr bptr_out);
    virtual ~actuator_everything();

    bool compress();

    bool decompress();

private:
    block_rough_ptr indata;
    block_rough_ptr outdata;
};

#endif