#ifndef _ENGINE_H_
#define _ENGINE_H_

#include "manager.h"
#include "block.h"
#include "io.h"
#include "blocking_collection.h"
#include <list>
#include "pbgz_file.h"
#include "actuator_fastq.h"
#include "actuator_everything.h"
#include "reference.h"

#define check_bcstatus(s)                                                       \
    if (status != block_collection::BlockingCollectionStatus::Ok)               \
    {                                                                           \
        check_exit(false, ERR_INTERNEL, "%s: %d", s, static_cast<int>(status)); \
    }

class engine
{
    typedef block_collection::BlockingCollection<block_rough_ptr> bcollect;
    typedef bcollect *bcollect_ptr;
    typedef block_collection::BlockingCollectionStatus bcstat;

    /* block id, pair< pair<block id, out buffer>, pair<offset, handler len>> */
    typedef block_collection::BlockingCollection<std::pair<std::pair<int64_t, uint8_t *>, std::pair<int64_t, int64_t>>> bcollect_refe;
    typedef bcollect_refe *bcollect_refe_ptr;

public:
    engine();
    virtual ~engine();

    /*  启动任务 */
    bool start();

    /*  根据文件类型获取需要设置的cache buffer的长度*/
    uint32_t get_cachelen(const filetype &ft);

private:
    /* 压缩时获取原数据 */
    const block_rough_ptr encoder_readdata(io &io);

    /* 解压时获取原数据 */
    const block_rough_ptr decoder_readdata();

    /* 压缩后数据落磁盘 */
    void encoder_task_outdata();

    /* 解压后数据落磁盘 */
    void decoder_task_outdata();

    /* 写一个压缩后的数据块到文件句柄中 */
    void write_encoded_block(const block_rough_ptr &block);

    /* 压缩数据 */
    bool encoder(const block_rough_ptr &bptr_in, block_rough_ptr &bptr_out);

    /* 解压数据 */
    bool decoder(const block_rough_ptr &bptr_in, block_rough_ptr &bptr_out);

    /* 压缩完成进行收尾工作 */
    void encoder_close();

    /* 解压完成进行收尾工作 */
    void decoder_close();

private:
    info conf;
    filetype ftype;
    int64_t file_size;

    /* pbgz 文件对应的句柄 */
    pbgz_file *fp_pbgz;

    /* 解压时输出io */
    io *io2write;

    int64_t block_cnt; /* 当前block的计数 */

    /*  原始数据存储队列，长度为1，读线程解析数据往里塞，压缩或者解压线程从里取 */
    bcollect_ptr block_input;

    /*  原始数据处理线程 */
    std::thread *block_input_thread;

    /* 原始数据的缓存队列，该队列用于并发控制 */
    bcollect_ptr block_input_pool;

    /* 输出数据的缓存队列，该队列用于并发控制 */
    bcollect_ptr block_output_pool;

    /* 输出数据写磁盘前的缓存队列，块的输出顺序是不定的，需要缓存排序写磁盘 */
    bcollect_ptr block_output_sort;

    /*  输出数据线程 */
    std::thread *block_output_thread;

    /*  任务线程池 */
    std::vector<thread *> task_pool;

    /* 分配的block池 */
    std::vector<block_rough_ptr> block_pool;

    /* 参考基因组 */
    reference *refgene;

    /* 文件的全局meta信息 */
    Json::Value file_meta;

private:

    /*  保存参考基因组 */
    int64_t reference_pack(int64_t &max_block_len, int64_t &total_enclen);

    /*  解压参考基因组 */
    void reference_unpack();

    /* 解压时初始化参考基因组 */
    bool decoder_initialize_reference();

    /* 标识是从压缩文件里解压的参考基因组 */
    bool is_unpack_refe;

    /* 当前读长度 */
    int64_t currlen_readed;

    /* 当前写长度 */
    int64_t currlen_writed;
};

#endif