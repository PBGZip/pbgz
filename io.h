#ifndef _IO_H_
#define _IO_H_

#include <isa-l/igzip_lib.h>
#include "block.h"
#include "manager.h"
#include <condition_variable>
#include <thread>
#include <htslib/bgzf.h>
#include "wait_event.h"

/*  支持以对应格式（如gz) 读取或者压缩数据 */

enum iomode
{
    mr,
    mw,
    mwgz
}; /* 读|写|写gz，读模式时文件类型由io自己判断并调用对应的接口 */

class io
{
public:
    io(const std::string &file, const iomode &mode, bool rw2cache = true);

    virtual ~io();

    /*  读一个完整的块，如读一个完整的fastq，最大读len数据，返回实际读到的长度，块的类型 */
    virtual int64_t read_one_block(block_rough &block);

    /* 读指定长度的数据，返回实际读到的长度，为0表示结束（内部会处理出错的情形）*/
    virtual int64_t read(uint8_t *buffer, int64_t len);

    /*  写len长度数据，返回实际写的数据长度 */
    virtual int64_t write(const uint8_t *buffer, int64_t len, bool eof);

    /*  写操作时刷新文件io */
    virtual bool flush();

    /*  写头信息*/
    virtual int64_t write_header(const uint8_t *header, int64_t len);

    /*  初始化块切割时cache长度 */
    virtual void init_bsplit_cache(const uint32_t &len);

    /*  读模式时获取文件的size */
    virtual int64_t get_filesize();

    /* 将文件指针偏移到离头offset的位置 */
    virtual void fseek2pos(int64_t offset);

private:
    /*  根据文件类型初始化文件的读写 */
    virtual bool initialize();

    /* 内部读接口 */
    virtual int64_t read_inner(uint8_t *buffer, int64_t len);

    /*  初始化后该函数会指向文件格式对应的read file函数 */
    int64_t (io::*read_fun)(uint8_t *, int64_t);

    /*  读文件，读len长度，返回实际读到的长度 */
    virtual int64_t read_file(uint8_t *buffer, int64_t len);

    /*  读gz格式文件，读len长度，返回实际读到的长度 */
    virtual int64_t read_gzfile(uint8_t *buffer, int64_t len);

    /*  读gz格式文件，读len长度，返回实际读到的长度，指令加速方式读 */
    virtual int64_t read_gzfile_fast(uint8_t *buffer, int64_t len);

    /* 内部写接口 */
    virtual void write_inner(const uint8_t *buf, int64_t size, bool eof);

    /*  初始化后该函数会指向文件格式对应的write file函数 */
    void (io::*write_fun)(const uint8_t *, int64_t);

    /* 写文件 */
    virtual void write_file(const uint8_t *buf, int64_t size);

    /* 写gz文件 */
    virtual void write_gzfile(const uint8_t *buf, int64_t size);

    /* 写gz文件，快速 */
    virtual void write_gzfile_fast(const uint8_t *buf, int64_t size);

private:
    FILE *fp;
    std::string file_name;
    int64_t file_size;
    bool eof;

    /* 块做切割时的缓存buffer */
    uint8_t *cache_bsplit;
    int64_t cache_bsplit_capacity;
    int64_t cache_bsplit_current;
    iomode m;
    int32_t parallel;

    bool is_gzipfile;
    BGZF *fpGZ;

    int64_t len_firstread; /* 第一次读的长度，用于预估gzip文件的压缩率 */

    /* 支持硬件指令集加速 */
    bool support_simd;

    /*  设计ring buffer给读写加缓存，防止调用者不断地读写小文件影响性能 */
private:
    uint8_t *ringbuffer;         /* ring buffer，类似linux内核的kfifo */
    int64_t ringbuffer_capacity; 
    int64_t ringbuffer_read_pos;
    int64_t ringbuffer_write_pos;
    std::thread *ringbuffer_guard; /* 守护线程 */
    wait_event put_event;
    wait_event get_event;
    bool ringbuffer_enable;
    bool ringbuffer_exit; /* 已经完成了，退出 */

    /* 从ringbuffer里拿数据 */
    virtual int64_t ringbuffer_get(uint8_t *buffer, int64_t len);

    /* 往ring buffer写数据 */
    virtual void ringbuffer_put(const uint8_t *buf, int64_t size);

    /* ring buffer剩余空间长度 */
    int64_t ringbuffer_space_left() const;

    /* ring buffer当前数据长度 */
    int64_t ringbuffer_len() const;

    /* ring buffer是否为空 */
    bool ringbuffer_empty() const;

    /* ring buffer 是否已满 */
    bool ringbuffer_full() const;

    /*  for intel gz */
    int64_t isal_remain_gzlen;       /* gz文件剩余未读的长度 */
    struct inflate_state isal_state; /* 解压 */
    struct isal_gzip_header isal_header;
    uint32_t isal_extra_each;   /*  每次补充的gz原始数据的长度 */
    uint8_t *isal_buffer_extra; /* 用于需要补充数据时的临时使用的buffer */
    uint32_t isal_extra_len;
    uint8_t *isal_buffer_in; /* 存读进来的gzip格式的原始数据 */
    int64_t isal_in_len;
};

#endif