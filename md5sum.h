#ifndef _MD5SUM_H_
#define _MD5SUM_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <string>

#undef BIG_ENDIAN_HOST
typedef unsigned int u32;

typedef struct
{
    u32 A, B, C, D; /* chaining variables */
    u32 nblocks;
    unsigned char buf[64];
    int count;

    std::string hexstr()
    {
        if (hexval.empty())
        {
            int len = (16 << 1) + 1;
            char str[len];
            snprintf(str, sizeof(str),
                     "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x",
                     buf[0], buf[1],  buf[2],  buf[3],    buf[4],   buf[5],    buf[6],    buf[7],
                     buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15]);
            hexval = std::string(str, len);
        }
        return hexval;
    }
    std::string hexval;
} MD5_CONTEXT;

void md5_init(MD5_CONTEXT *ctx);

void md5_write(MD5_CONTEXT *hd, unsigned char *inbuf, size_t inlen);

void md5_final(MD5_CONTEXT *hd);

#endif