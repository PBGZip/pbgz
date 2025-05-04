#ifndef _CODER_H_
#define _CODER_H_

#include "coder_io.h"
#include <stdint.h>

class coder
{
public:
    virtual ~coder() {}
    virtual int32_t decode_line([[maybe_unused]] uint8_t *dst, 
                                [[maybe_unused]] int32_t len, 
                                [[maybe_unused]] uint8_t split_ch = '\n', 
                                [[maybe_unused]] bool need2hold = false) { return 0; }
    virtual int32_t decode_line([[maybe_unused]] uint8_t *dst, 
                                [[maybe_unused]] int32_t len, 
                                [[maybe_unused]] uint8_t *rely = nullptr, 
                                [[maybe_unused]] uint8_t split_ch = '\n', 
                                [[maybe_unused]] bool need2hold = false)
    {
        return 0;
    }

    virtual void set_level(int32_t level)
    {
        io->set_level(level);
    }

protected:
    coder_io *io;
};

#endif