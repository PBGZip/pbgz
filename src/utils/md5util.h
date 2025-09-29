#pragma once

#include <string>
#include <stdint.h>

#include "md5sum.h"

inline void calcMd5sum(std::string& md5, const uint8_t* data, uint32_t dataLen) {
    MD5_CONTEXT md5Contex;
    md5_init(&md5Contex);
    md5_write(&md5Contex, const_cast<uint8_t*>(data), dataLen);
    md5_final(&md5Contex);
    md5 = md5Contex.hexstr();
}