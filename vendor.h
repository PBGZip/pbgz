#ifndef _VERSION_H_
#define _VERSION_H_

#include <iostream>
#include "singleton.h"

#define ZIP_SUFFIX ".pbgz"
#define GZIP_SUFFIX ".gz"

/* 流的标识，如果后面流格式有修改，可以把改为其他字样  */
#define STREAM_ID "PBSM"

/* pbgz zip格式文件的标识id */
#define FILE_ID "PBGZ"

/*  块大小如果设置太小，可能会导致一个块不能包含一个完整的fastq 4行，
  * 从而导致fastq块被当作二进制压缩，导致压缩率降低 
  */
#define BLOCK_SIZE (268435456) /*  block size: fastq/ binary*/
// #define BLOCK_SIZE (20971520) /*  block size: fastq/ binary*/

/* mapping上的阀值 */
#define MAPPED_THRESHOLD_GEN2 2

#define MAJOR 1
#define MINOR 0
#define PATCH 0

class version
{
private:
    std::string v;

public:
    version(/* args */)
    {
        v = std::to_string(MAJOR);
        v += ".";
        v += std::to_string(MINOR);
        v += ".";
        v += std::to_string(PATCH);
    }
    ~version(){}

    const std::string ver() const { return this->v; }
};

typedef singleton<version> ver;

#endif