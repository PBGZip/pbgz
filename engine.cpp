#include "engine.h"
#include "io.h"
#include "file_op.h"
#include "md5sum.h"
#include "coder/coder_json.h"
#include "coder/coder_ppmd.h"
#include "coder/coder_bwt_cm.h"
#include "coder/coder_fc.h"

engine::engine()
{
    block_cnt = 0;
    ftype = TYPE_UNKNOW;
    conf = manage::instance().get_zipinfo();
    fp_pbgz = nullptr;
    file_size = 0;
    block_input_thread = nullptr;
    block_output_thread = nullptr;
    io2write = nullptr;
    refgene = nullptr;
    block_pool.clear();
    currlen_readed = currlen_writed = 0;
}

engine::~engine()
{
    if (fp_pbgz)
        delete fp_pbgz;
    if (refgene)
        delete refgene;
}

void __debug_test_coder__()
{
    const int src_len = 12 << 20;
    uint8_t *src;
    safe_alloc(src_len, uint8_t, src);
    printf("src_len %u\n", src_len);
    for(int n = 0; n < src_len; n++)
        src[n] = n;
    uint32_t out_len = src_len << 2;
    uint8_t *out;
    safe_alloc(out_len, uint8_t, out);
    uint8_t *decout;
    safe_alloc(out_len, uint8_t, decout);
    coder_io match_io(out, out_len);

    // coder_bwt_cm match_cm(&match_io);
    coder_fc match_cm(&match_io);

    match_cm.encode_line((uint8_t *)src, src_len);
    {
        FILE *fp = fopen("fakedata", "wb");
        fwrite(src, src_len, 1, fp);
        fflush(fp);
        fclose(fp);
    }
    match_cm.encode_flush();
    
    printf("match_io len %u\n", match_io.data_len);

    coder_io match_io1(match_io.data, match_io.data_len);

    // coder_bwt_cm match_cm1(&match_io1);
    coder_fc match_cm1(&match_io1);

    int declen = match_cm1.decode_line((uint8_t *)decout, out_len, UINT8_MAX, false);
    check_exit(declen == src_len, ERR_INTERNEL, "error 1");

    {
        FILE *fp = fopen("fakedata.deocde", "wb");
        fwrite(decout, declen, 1, fp);
        fflush(fp);
        fclose(fp);
    }
    _Exit(0);
}

 /*  启动任务 */
bool engine::start()
{
    fprintf(stderr, "\033[37mpbgz version => %d.%d.%d\033[0m\n\n", MAJOR, MINOR, PATCH);
    fprintf(stderr, "Parallel set: %d\n", manage::instance().get_zipinfo().get_threads());
    
    if (conf.get_mode() == ZIP) /*  压缩模式 */
    {
        int64_t n, size_origin = 0, size_dest = 0;
        bcstat status;
        timer cost_ms(true);
        block_rough_ptr bptr_in, bptr_out;
        
        blocktype btype = TYPE_UNKNOW;
        ftype= judge_file_type(conf.get_infile());
        if (block_type_isfastq(ftype) && !conf.get_refgene().empty())
            safe_new_class(reference, manage::instance().get_zipinfo().get_refgene(), refgene);
        if (refgene && !refgene->make_index())
            refgene = nullptr;
        fprintf(stderr, "File type: %s\n", get_typename(ftype).get());

        fp_pbgz = new pbgz_file(conf.get_outfile(), false);
        check_exit(fp_pbgz, ERR_MEM_NOENOUGH, "pbgz file fp initialize failed: no enough memory");

        /* 读数据任务并发为1 */
        safe_new_class(bcollect, 1, block_input);

        /* 压缩数据队列长度为并发数 */
        safe_new_class(bcollect, conf.get_threads(), block_input_pool);
        n = block_input_pool->bounded_capacity();

        /* 压缩后输出数据长度设置并发数两倍，因为数据落磁盘需要时间 */
        safe_new_class(bcollect, (n << 1), block_output_pool);
        safe_new_class(bcollect, (n << 1), block_output_sort);

        if ((block_output_pool->bounded_capacity() != ((block_input_pool->bounded_capacity()) << 1)) || 
            (block_output_pool->bounded_capacity() != block_output_sort->bounded_capacity()))
            check_exit(false, ERR_INTERNEL, "buffer pool capacity check fail: output  %ld, input %ld", 
                block_output_pool->bounded_capacity(), block_input_pool->bounded_capacity());

        io read_io(conf.get_infile(), mr, true);
        read_io.init_bsplit_cache(get_cachelen(ftype));
        file_size = read_io.get_filesize();

        /* 启动写数据线程，将压缩后输出落磁盘 */
        encoder_task_outdata();
      
        for(;;) {
        
            bptr_in = encoder_readdata(read_io);
            // fprintf(stderr, "read data get block_cnt: %d, readed len %ld\n", block_cnt, bptr_in->current_len);

            if (bptr_in->current_len == 0) { /*  待压缩原始数据已经读完，塞结束标识，压缩线程拿到结束标识后会退出*/
                for (n = 0; n < std::min(std::size_t(conf.get_threads()), task_pool.size()); n++) {
                    status = block_input->add(nullptr);
                    check_bcstatus("block input add finish flag failed with errno");
                }
                // block_input->complete_adding();
                break;
            }
            size_origin += bptr_in->current_len;
            status = block_input->add(bptr_in);
            check_bcstatus("block input add data failed with errno");

            if (block_cnt <= conf.get_threads()) {
                task_pool.push_back(new  std::thread([&]() {

                    bcstat status;
                    block_rough_ptr bptr_in, bptr_out;

                    for (;;) {
                        /* 拿原始数据 */
                        status = block_input->take(bptr_in);
                        if (!bptr_in) {/* 拿到结束标识退出 */
                            // fprintf(stderr, "compress exit one\n");
                            break;
                        }
                        // fprintf(stderr, "compress get block_cnt: %d, readed len %ld\n", block_cnt, bptr_in->current_len);
                        check_bcstatus("block input take data failed with errno");
                        
                        /* 拿输出buffer空间*/
                        status = block_output_pool->take(bptr_out);
                        bptr_out->reset();
                        check_bcstatus("block output pool take data failed with errno");

                        /* 开始压缩 */
                        check_exit(encoder(bptr_in, bptr_out), ERR_INTERNEL, "encoder block %ld failed", bptr_in->block_id);

                        /* 压缩完成 ，将block塞回队列池 */
                        bptr_in->reset();
                        status = block_input_pool->add(bptr_in);
                        check_bcstatus("block input pool add data failed with errno");
                    }
                }));
                check_exit(task_pool.back(), ERR_MEM_NOENOUGH, "io task thread initialize failed: no enough memory");
            }
        }

        for (const auto &p : task_pool) {
            p->join();
            delete p;
        }
        status = block_output_sort->add(nullptr); /* 数据已经读完了，告知write线程写完buffer里数据后可以退出了 */
        check_bcstatus("block output write add finish flag failed with errno");
        block_output_thread->join();
        delete block_output_thread;

        for(auto&b: block_pool)
            delete b;
        delete block_output_pool;
        delete block_output_sort;
        delete block_input;
        delete block_input_pool;

        encoder_close(); /* 收尾工作，写meta或reference信息 */
        size_dest = file_length(conf.get_outfile());

        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\033[0m\n", "from/to", currlen_readed, currlen_writed);
        fprintf(stderr, "\nCompress finish, cost %lum%lus.\n", (cost_ms.elapsed() / 1000) / 60, (cost_ms.elapsed() / 1000) % 60);
        fprintf(stderr, "Total size_dest size %ld bytes, compressed to %ld bytes, ratio %0.2f%\n", size_origin, size_dest, (size_dest * 1.0) * 100 / size_origin);
    }

    else if (conf.get_mode() == UNZIP)
    { /* 解压模式 */
        int64_t n;
        bcstat status;
        timer cost_ms(true);
        block_rough_ptr bptr_in, bptr_out;

        fp_pbgz = new pbgz_file(conf.get_infile(), true);
        check_exit(fp_pbgz, ERR_MEM_NOENOUGH, "pbgz file fp initialize failed: no enough memory");

        fp_pbgz->get_file_meta(file_meta);
        // cerr << file_meta << endl;
        fprintf(stderr, "File type: %s\n", file_meta["origin_filetype"].asString().c_str());

        /* 读数据任务并发为1 */
        safe_new_class(bcollect, 1, block_input);

        /* 设置解压数据队列长度为并发数 */
        safe_new_class(bcollect, conf.get_threads(), block_input_pool);
        n = block_input_pool->bounded_capacity();

        /* 解压后输出数据长度设置并发数两倍，因为数据落磁盘需要时间 */
        safe_new_class(bcollect, (n << 1), block_output_pool);
        safe_new_class(bcollect, (n << 1), block_output_sort);

        if ((block_output_pool->bounded_capacity() != ((block_input_pool->bounded_capacity()) << 1)) || 
            (block_output_pool->bounded_capacity() != block_output_sort->bounded_capacity()))
            check_exit(false, ERR_INTERNEL, "buffer pool capacity check fail: output  %ld, input %ld", 
                block_output_pool->bounded_capacity(), block_input_pool->bounded_capacity());

        if (file_meta["refe"].isObject()) /* 有使用参考基因组压缩 */
            check_exit(decoder_initialize_reference(), ERR_INTERNEL, "initialize reference failed");

        io2write = new io(conf.get_outfile(), conf.get_decompress2gz() ? mwgz : mw, false); // 读写cache不稳定，先关闭
        check_exit(io2write, ERR_MEM_NOENOUGH, "io write initialize failed: no enough memory");

        /* 启动写数据线程，将解压后输出落磁盘 */
        decoder_task_outdata();
      
        for(;;) {
        
            bptr_in = decoder_readdata();
            // fprintf(stderr, "read data get block_cnt: %d, readed len %ld\n", block_cnt, bptr_in->current_len);

            if (bptr_in->current_len == 0) { /*  待压缩原始数据已经读完，塞结束标识，压缩线程拿到结束标识后会退出 */
                for (n = 0; n < std::min(std::size_t(conf.get_threads()), task_pool.size()); n++) {
                    status = block_input->add(nullptr);
                    check_bcstatus("block input add finish flag failed with errno");
                }
                // block_input->complete_adding();
                break;
            }
            status = block_input->add(bptr_in);
            check_bcstatus("block input add data failed with errno");

            if (block_cnt <= conf.get_threads()) {
                task_pool.push_back(new  std::thread([&]() {

                    bcstat status;
                    block_rough_ptr bptr_in, bptr_out;

                    for (;;) {
                        /* 拿原始数据 */
                        status = block_input->take(bptr_in);
                        if (!bptr_in) {/* 拿到结束标识退出 */
                            // fprintf(stderr, "decompress exit one\n");
                            break;
                        }
                        // fprintf(stderr, "decompress get block_cnt: %d, readed len %ld\n", block_cnt, bptr_in->current_len);
                        check_bcstatus("block input take data failed with errno");
                        
                        /* 拿输出buffer空间*/
                        status = block_output_pool->take(bptr_out);
                        bptr_out->reset();
                        check_bcstatus("block output pool take data failed with errno");

                        /* 开始解压 */
                        check_exit(decoder(bptr_in, bptr_out), ERR_INTERNEL, "decoder block %ld failed", bptr_in->block_id);

                        /* 解压完成 ，将block塞回队列池 */
                        bptr_in->reset();
                        status = block_input_pool->add(bptr_in);
                        check_bcstatus("block input pool add data failed with errno");
                    }
                }));
                check_exit(task_pool.back(), ERR_MEM_NOENOUGH, "io task thread initialize failed: no enough memory");
            }
        }

        for (const auto &p : task_pool) {
            p->join();
            delete p;
        }
        status = block_output_sort->add(nullptr); /* 数据已经读完了，告知write线程写完buffer里数据后可以退出了 */
        check_bcstatus("block output write add finish flag failed with errno");
        block_output_thread->join();
        delete block_output_thread;

        for(auto&b: block_pool)
            delete b;
        delete block_output_pool;
        delete block_output_sort;
        delete block_input;
        delete block_input_pool;

        decoder_close(); /* 收尾工作 */
        delete io2write;
        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\033[0m\n", "from/to", currlen_readed, currlen_writed);
        fprintf(stderr, "\nDecompress finish, cost %lum%lus.\n", (cost_ms.elapsed() / 1000) / 60, (cost_ms.elapsed() / 1000) % 60);

        return true;
    }

    else
    {
        check(false, false, "mode is invalid");
    }

    return true;
}

/* 压缩时获取原数据 */
const block_rough_ptr engine::encoder_readdata(io &io)
{
    int64_t len_readed, n;
    block_rough_ptr bptr_in, bptr_out;
    bcstat status;

    /* 当小于并发队列长度时分配block空间 */
    if (block_cnt < conf.get_threads()) {
        safe_new_class(block_rough, BLOCK_SIZE, bptr_in);
        status = block_input_pool->add(bptr_in);
        check_bcstatus("block output add data failed with errno");
        block_pool.push_back(bptr_in);

        for (n = 0; n < 2; n++) {
            safe_new_class(block_rough, BLOCK_SIZE, bptr_out);
            status = block_output_pool->add(bptr_out);
            check_bcstatus("block output add data failed with errno");
            block_pool.push_back(bptr_out);
        }
    }

    /* 从缓存队列里取空闲的block，没有空闲block buffer时会阻塞等待 */
    status = block_input_pool->take(bptr_in);
    check_bcstatus("block input pool take data failed with errno");

    len_readed = io.read_one_block(*bptr_in);
    if (len_readed > 0)
        bptr_in->block_id = ++block_cnt;
    bptr_in->current_len = len_readed;
    currlen_readed += len_readed;
    // if (bptr_in->block_id == 1)
        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
    return bptr_in;
}

/* 写一个压缩后的数据块到文件句柄中 */
void engine::write_encoded_block(const block_rough_ptr &block)
{
    nstype stream_type = (block_type_isfastq(block->btype)) ? NST_FASTQ_DATA_BLOCK : NST_OTHER_DATA_BLOCK;
    
    /* 写当前block的meta压缩后的流 */
    fp_pbgz->write_one_stream(block->buffer + block->current_len - block->meta_encoded_len, block->meta_encoded_len, NST_META_BLOCK);

    /* 写当前block的数据流 */
    fp_pbgz->write_one_stream(block->buffer, block->current_len - block->meta_encoded_len, stream_type);
}

/* 压缩后数据落磁盘 */
void engine::encoder_task_outdata()
{
    if (block_output_thread)
        return;
    block_output_thread = new  std::thread([&]() {
        bcstat status;
        int64_t block2write = 1;
        block_rough_ptr bptr;
        std::list<block_rough_ptr> sort_cache;

        for (;;) {
            status = block_output_sort->take(bptr);
            check_bcstatus("block output write take data failed with errno");

            if (!bptr) { /* 不会有新的数据被压缩了，处理完缓存的数据就可以退出了 */
                while (!sort_cache.empty()) {
                    write_encoded_block(sort_cache.front());
                    currlen_writed = fp_pbgz->get_offset();
                    fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
                    sort_cache.pop_front();
                }
                return;

            } else { 
                /* 排序将block id小的放前面 */
                sort_cache.push_back(bptr);
                sort_cache.sort([](const block_rough_ptr &p1, const block_rough_ptr &p2 ){
                    return p1->block_id <= p2->block_id;
                });

                /* 循环查找下一个写的block，如果没找到则退出，找到了则继续查找下一个 */
                while (!sort_cache.empty()) {
                    const block_rough_ptr &current = sort_cache.front();
                    if (current->block_id == block2write) {
                        write_encoded_block(current);
                        currlen_writed = fp_pbgz->get_offset();
                        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
     
                        /* 当前encoder buffer已经用完了，丢入pool给压缩线程重新使用 */
                        status = block_output_pool->add(current);
                        check_bcstatus("block output pool add data failed with errno");

                        sort_cache.pop_front();
                        block2write++;
                    } else 
                        break;
                }
            }
        }
    });
    check_exit(block_output_thread, ERR_MEM_NOENOUGH, "io block output thread initialize failed: no enough memory");
}

/* 压缩数据 */
bool engine::encoder(const block_rough_ptr &bptr_in, block_rough_ptr &bptr_out)
{
    bool res = false;
    bcstat status;

    // fprintf(stderr, "Block type ------[%ld][%s]------\n", bptr_in->block_id, get_typename(bptr_in->btype).get());

    if (block_type_isfastq(bptr_in->btype))
    {
        actuator_fastq actuator(bptr_in, bptr_out, refgene);
        if (actuator.analyze_fastq())
            res = actuator.compress();
        else {
            actuator_everything actuator(bptr_in, bptr_out);
            res = actuator.compress();
        }
    } else {
        actuator_everything actuator(bptr_in, bptr_out);
        res = actuator.compress();
    }
    check(res, false, "compress false: block id: %ld", bptr_out->block_id);

    // fprintf(stderr, "Block type ------[%ld][%s]------done\n", bptr_in->block_id, get_typename(bptr_in->btype).get());

    /* 将压缩后数据写入磁盘前的缓存队列，块的输出顺序是不定的，需要缓存排序写磁盘 */
    status = block_output_sort->add(bptr_out);
    check_bcstatus("block output write add data failed with errno");
    return true; 
}

/* 解压数据 */
bool engine::decoder(const block_rough_ptr &bptr_in, block_rough_ptr &bptr_out)
{
    bool res = false;
    bcstat status;
    int64_t sid_meta, sid_data, enclen_meta;
    nstype stype;
    pbgz_stream_header header_meta(bptr_in->buffer);

    /* 检测block 的meta流 */
    check(header_meta.check_mark(), false, "header check failed: meta id is not: %s", STREAM_ID);
    sid_meta = header_meta.get_subid();
    enclen_meta = header_meta.get_datalen();
    stype = (nstype)(header_meta.get_type());

    check(stype == NST_META_BLOCK, false,
          "header check failed: type is not %d, current is %d", NST_META_BLOCK, stype);
    
    /* 检测block 的data流 */
    pbgz_stream_header header_data(bptr_in->buffer + header_meta.get_bufferlen() + enclen_meta);
    check(header_data.check_mark(), false, "header check failed: data id is not: %s", STREAM_ID);
    sid_data = header_data.get_subid();
    stype = (nstype)(header_data.get_type());
    check(sid_data == sid_meta + 1, false, "pair stream check failed: data id  (%ld) != meta id + 1 (%ld)", sid_data, sid_meta + 1);

    // fprintf(stderr, "enclen_meta: %ld, encdata len: %ld, stype %d\n", enclen_meta, header_data.get_datalen(), stype);

    if (NST_FASTQ_DATA_BLOCK == stype) {
        actuator_fastq actuator(bptr_in, bptr_out, refgene);
        res = actuator.decompress();
    } else if(NST_OTHER_DATA_BLOCK == stype) {
        actuator_everything actuator(bptr_in, bptr_out);
        res = actuator.decompress();
    } else {
        check_exit(false, ERR_INTERNEL, "undefined decompress type.");
    }
    check(res, false, "decompress false: block id: %ld", bptr_out->block_id);

      /* 将压缩后数据写入磁盘前的缓存队列，块的输出顺序是不定的，需要缓存排序写磁盘 */
    status = block_output_sort->add(bptr_out);
    check_bcstatus("block output write add data failed with errno");
    return true;
}

/* 压缩完成进行收尾工作 */
void engine::encoder_close()
{
    Json::Value refe_meta;
    coder_json cjson;
    std::string file_meta_enc, ft, ni_name, fasta_name;
    int64_t offset1, offset2, block_cnt;
    int64_t max_reflen, total_enclen;
    std::shared_ptr<uint8_t> type;
    type = get_typename(ftype);
    MD5_CONTEXT md5;

    /*  如果有使能pack参考基因组, 则需要将有效的参考基因组pack进压缩文件*/
    if (refgene) {
        file_abspath_filename(refgene->get_fasta_name(), fasta_name);
        refe_meta["squash_len"] = (Json::Value::Int64)(refgene->get_squashlen());
        refe_meta["fasta_name"] = fasta_name;
        refe_meta["fasta_len"] = (Json::Value::Int64)(file_length(refgene->get_fasta_name()));
        refe_meta["fasta_md5"] = refgene->get_fasta_md5();
        file_abspath_filename(refgene->get_ni_path(), ni_name);
        refe_meta["ni_name"] = ni_name; /* 包含了md5信息，用于解压校验 */

        if (!conf.get_nopackref()) {
            timer cost_ms(true);
            refe_meta["file_offset"] = (Json::Value::Int64)(fp_pbgz->get_offset()); /* 记录pack的reference的偏移 */
            block_cnt = reference_pack(max_reflen, total_enclen);
            // md5_init(&md5);
            // md5_write(&md5, (uint8_t *)(refgene->get_squash()), refgene->get_squashlen());
            // md5_final(&md5);
            fprintf(stderr, "reference: origin len %lu -> packed len %lu, cost ms %lu\n",
                    refgene->get_squashlen(), total_enclen, cost_ms.elapsed());
            refe_meta["max_block_len"] = (Json::Value::Int64)(max_reflen);
            refe_meta["blocks"] = (Json::Value::Int64)(block_cnt);
            // refe_meta["md5"] = md5.hexstr();
        }
    }
    file_meta["refe"] = refe_meta;

    /* 先将没有写完的数据都写到磁盘 */
    fp_pbgz->flush();

    /* 写文件对应的全局meta信息流，写两份 */
    file_meta["origin_filename"] = conf.get_infile();
    file_meta["origin_filetype"] = std::string((char *)(type.get()), strlen((char *)(type.get())));
    file_meta["origin_filesize"] = (Json::Value::Int64)(file_size);
    file_meta["block_size"] = (Json::Value::Int)BLOCK_SIZE;
    file_meta["compress_cmd"] = conf.get_cmd();
    file_meta["stream_offsets"] = fp_pbgz->get_stream_offset();
    file_meta["stream_count"] = fp_pbgz->get_subid() + 2; /* 两条文件meta流offset记录在第一条流中 */
    
    cjson.encoder(file_meta, file_meta_enc);

    offset1 = fp_pbgz->get_offset();
    fp_pbgz->write_one_stream((uint8_t *)(file_meta_enc.c_str()), file_meta_enc.length(), NST_META_FILE);
    offset2 = fp_pbgz->get_offset();
    fp_pbgz->write_one_stream((uint8_t *)(file_meta_enc.c_str()), file_meta_enc.length(), NST_META_FILE);

    /* 将全局meta流的地址更新至第一条流 */
    fp_pbgz->set_meta_offset(offset1, offset2);

    /* 再刷一次磁盘将刚才写的数据落盘*/
    fp_pbgz->flush();
}

/* 解压完成进行收尾工作 */
void engine::decoder_close()
{
    io2write->flush();
}

/* 解压时获取原数据 */
const block_rough_ptr engine::decoder_readdata()
{
    int64_t len_readed, n;
    block_rough_ptr bptr_in, bptr_out;
    bcstat status; 
    int64_t offset_refe = 0;
    const int32_t block_size = file_meta["block_size"].asInt();
    if (file_meta["refe"].isObject() && file_meta["refe"]["file_offset"].isInt64())
        offset_refe = file_meta["refe"]["file_offset"].asInt64();

    /* 当小于并发队列长度时分配block空间 */
    if (block_cnt < conf.get_threads()) {
        safe_new_class(block_rough, block_size, bptr_in);
        status = block_input_pool->add(bptr_in);
        check_bcstatus("block output add data failed with errno");
        block_pool.push_back(bptr_in);

        for (n = 0; n < 2; n++) {
            safe_new_class(block_rough, block_size, bptr_out);
            status = block_output_pool->add(bptr_out);
            check_bcstatus("block output add data failed with errno");
            block_pool.push_back(bptr_out);
        }
    }

    /* 从缓存队列里取空闲的block，没有空闲block buffer时会阻塞等待 */
    status = block_input_pool->take(bptr_in);
    check_bcstatus("block input pool take data failed with errno");

    if (offset_refe > 0 && fp_pbgz->get_offset() == offset_refe) /* 已经解析到reference了 */
        len_readed = 0;
    else {
        len_readed = fp_pbgz->read_pair_stream(*bptr_in);
        if (len_readed > 0)
            bptr_in->block_id = ++block_cnt;
    }
    bptr_in->current_len = len_readed;
    currlen_readed += len_readed;
    // if (bptr_in->block_id == 1)
        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
    return bptr_in;
}

/* 解压后数据落磁盘 */
void engine::decoder_task_outdata()
{
    if (block_output_thread)
        return;
    block_output_thread = new  std::thread([&]() {
        bcstat status;
        int64_t block2write = 1;
        block_rough_ptr bptr;
        std::list<block_rough_ptr> sort_cache;

        for (;;) {
            status = block_output_sort->take(bptr);
            check_bcstatus("block output write take data failed with errno");

            if (!bptr) { /* 不会有新的数据被解压了，处理完缓存的数据就可以退出了 */
                while (!sort_cache.empty()) {
                    const block_rough_ptr &block = sort_cache.front();
                    io2write->write(block->buffer, block->current_len, false);
                    currlen_writed += block->current_len;
                    fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
                    io2write->flush();
                    sort_cache.pop_front();
                }
                return;

            } else { 
                /* 排序将block id小的放前面 */
                sort_cache.push_back(bptr);
                sort_cache.sort([](const block_rough_ptr &p1, const block_rough_ptr &p2 ){
                    return p1->block_id <= p2->block_id;
                });

                /* 循环查找下一个写的block，如果没找到则退出，找到了则继续查找下一个 */
                while (!sort_cache.empty()) {
                    const block_rough_ptr &current = sort_cache.front();
                    if (current->block_id == block2write) {
                        io2write->write(current->buffer, current->current_len, false);
                        currlen_writed += current->current_len;
                        fprintf(stderr, "\033[37m%-s ---[%ld/%ld]---\r", "from/to", currlen_readed, currlen_writed);
                        // io2write->flush();
     
                        /* 当前encoder buffer已经用完了，丢入pool给压缩线程重新使用 */
                        status = block_output_pool->add(current);
                        check_bcstatus("block output pool add data failed with errno");

                        sort_cache.pop_front();
                        block2write++;
                    } else 
                        break;
                }
            }
        }
    });
    check_exit(block_output_thread, ERR_MEM_NOENOUGH, "io block output thread initialize failed: no enough memory");
}

/*  根据文件类型获取需要设置的cache buffer的长度 */
uint32_t engine::get_cachelen(const filetype &ft)
{
    uint32_t len = (8 << 20);
    switch (ft)
    {
    case FASTQ_GEN2:
    case FASTQ_GEN2_GZIP:
        len = (GENE2_MAX_BASE) << 2;
        break;
    case FASTQ_GEN3:
    case FASTQ_GEN3_GZIP:
        len = (GENE3_MAX_BASE) << 2;
        break;
    case PBGZFILE:
    case GZIP:
    case BINARY:
    case BINARY_GZIP:
        len = 0;
        break;
    default:
        check_exit(false, ERR_INTERNEL, "can not set cache len");
        break;
    }
    return len;
}

/*  保存参考基因组 */
int64_t engine::reference_pack(int64_t &max_block_len, int64_t &total_enclen)
{
    std::mutex m;
    uint8_t *output;
    int64_t each, current, remain, total;
    int64_t n, pcnt = conf.get_threads();
    int64_t block = 0, offset = 0;
    bcollect_refe_ptr input_pool;
    bcollect_refe_ptr input_block;
    std::pair<std::pair<int64_t, uint8_t *>, std::pair<int64_t, int64_t>> refe2do;
    pbgz_file *fp2write = fp_pbgz;
    std::vector<std::thread> tpools;
    reference *refe = refgene;

    total = refgene->get_squashlen();
    each = (16 << 20);
    remain = total;
    max_block_len = total_enclen = 0;
    safe_alloc((pcnt * each), uint8_t, output);
    safe_new_class(bcollect_refe, pcnt, input_pool);
    safe_new_class(bcollect_refe, 1, input_block);

    for (;;)
    {
        current = std::min(each, remain);

        if (block < pcnt)
        {
            input_pool->add(std::make_pair(std::make_pair(block, output + (block * each)), std::make_pair(offset, current)));

            tpools.push_back(std::thread([&input_pool, &output, &input_block, &m, &fp2write, &refe, &each, &max_block_len, &total_enclen]()
                                         {
                                             for (;;)
                                             {
                                                 Json::Value meta;
                                                 coder_json cmeta;
                                                 std::string meta_string;
                                                 int64_t curr_block_len, plen;
                                                 std::pair<std::pair<int64_t, uint8_t *>, std::pair<int64_t, int64_t>> refe2do;

                                                 input_pool->take(refe2do);
                                                 plen = refe2do.second.second;
                                                 if (plen == 0)
                                                     return;
                                                 const uint8_t *p = refe->get_squash() + refe2do.second.first;
                                                 refe->reference_squash_sanitize(refe2do.second.first, refe2do.second.second);
                       
                                                 coder_io refe_io(refe2do.first.second, each);
                                                 coder_ppmd cppmd(&refe_io);
                                                 cppmd.encode(p, plen);
                                                 cppmd.encode_flush();
                                                //  coder_bwt_cm match_cm(&refe_io);
                                                //  match_cm.encode_line(p, plen);
                                                //  match_cm.encode_flush();

                                                 meta["block"] = (Json::Value::Int)refe2do.first.first;
                                                 meta["offset"] = (Json::Value::Int64)refe2do.second.first;
                                                 meta["srclen"] = (Json::Value::Int)plen;
                                                 meta["dstlen"] = (Json::Value::Int)refe_io.data_len;
                                                 meta["coder"] = refe_io.meta;
                                                 cmeta.encoder(meta, meta_string); /*  压缩block meta */
                                                
                                                 m.lock();
                                                 curr_block_len = meta_string.length() + refe_io.data_len;
                                                 check_exit(curr_block_len < plen, ERR_INTERNEL, "reference block pack failed");
                                                 /* 写当前block的meta压缩后的流 */
                                                 fp2write->write_one_stream((const uint8_t *)(meta_string.c_str()), meta_string.length(), NST_META_REFE);
                                                 /* 写当前block的数据流 */
                                                 fp2write->write_one_stream(refe_io.data, refe_io.data_len, NST_DATA_REFE);
                                                 max_block_len = (curr_block_len > max_block_len) ? curr_block_len : max_block_len;
                                                 total_enclen += curr_block_len;
                                                m.unlock();

                                                 input_block->add(refe2do);
                                             }
                                         }));
        }
        else
        {
            input_block->take(refe2do);
            refe2do.first.first = block;
            refe2do.second.first = offset;
            refe2do.second.second = current;
            input_pool->add(refe2do);
        }

        block++;
        offset += current;
        remain -= current;
        if (remain <= 0)
            break;
    }
    for (auto &t : tpools)
    {
        input_block->take(refe2do);
        refe2do.second.second = 0;
        input_pool->add(refe2do);
    }
    for (auto &t : tpools)
    {
        if (t.joinable())
            t.join();
    }
    delete input_pool;
    delete input_block;
    free(output);
    return block;
}

/*  解压参考基因组 */
void engine::reference_unpack()
{
    Json::Value meta_refe;
    int64_t len_readed, len_header, len;
    int64_t refe_squash_len, refe_offset, max_len;
    int64_t n, pcnt = conf.get_threads(), blocks;
    std::string fasta_name, md5_packed;
    uint8_t *refe_squash, *p;
    block_rough_ptr block[pcnt];
    block_rough_ptr current;
    bcollect_ptr input_pool;
    bcollect_ptr input_block;
    std::vector<std::thread> tpools;
    MD5_CONTEXT md5;
    pbgz_stream_header header;

    len_header = header.get_bufferlen();
    meta_refe = file_meta["refe"];
    refe_squash_len = meta_refe["squash_len"].asInt64();
    fasta_name = meta_refe["fasta_name"].asString();
    refe_offset = meta_refe["file_offset"].asInt64();
    max_len = meta_refe["max_block_len"].asInt64();
    max_len += (len_header << 1);
    blocks = meta_refe["blocks"].asInt64();
    md5_packed = meta_refe["md5"].asString();
    for (n = 0; n < pcnt; n++)
        safe_new_class(block_rough, max_len, block[n]);
    safe_new_class(bcollect, pcnt, input_pool);
    safe_new_class(bcollect, 1, input_block);

    safe_new_class(reference, fasta_name, refgene);
    refe_squash = refgene->initialize_squash_2(refe_squash_len);

    /* 起一个关闭ringbuffer的文件操作，因为需要fseek读取文件对应位置的pack的reference流 */
    io fileio_nokb(conf.get_infile(), mr, false);
    fileio_nokb.fseek2pos(refe_offset);
    // printf("\trefe_offset %ld\n", refe_offset);

    for (n = 0; n < blocks; n++)
    {
        if (n < pcnt)
        {
            {
                 /* reference block meta stream header */
                 p = block[n]->buffer;
                block[n]->buffer_size = fileio_nokb.read(p, len_header);
                pbgz_stream_header header_meta(p);
                p += len_header;
                /* reference block meta stream data */
                len = header_meta.get_datalen();
                block[n]->buffer_size += fileio_nokb.read(p, len);
                p += len;

                /* reference block data stream header */
                block[n]->buffer_size += fileio_nokb.read(p, len_header);
                pbgz_stream_header header_data(p);
                p += len_header;
                /* reference block data stream data */
                len = header_data.get_datalen();
                block[n]->buffer_size += fileio_nokb.read(p, len);
            }

            input_pool->add(block[n]);

            tpools.push_back(std::thread([&input_pool, &input_block, &refe_squash, &n]()
                                         {
                                             for (;;)
                                             {
                                                 nstype stype;
                                                 coder_json cmeta;
                                                 Json::Value block_meta;
                                                 int32_t src_len, dst_len, len;
                                                 int64_t sid_meta, sid_data, enclen_meta, offset;
                                                 block_rough_ptr curr_block;

                                                 input_pool->take(curr_block);
                                                 if (curr_block == nullptr)
                                                     return;

                                                 pbgz_stream_header header_meta(curr_block->buffer);

                                                 /* 检测block 的meta流 */
                                                 check_exit(header_meta.check_mark(), ERR_INTERNEL, "header check failed: meta id is not: %s", STREAM_ID);
                                                 sid_meta = header_meta.get_subid();
                                                 enclen_meta = header_meta.get_datalen();
                                                 stype = (nstype)(header_meta.get_type());
                                                 check_exit(stype == NST_META_REFE, ERR_INTERNEL,
                                                            "header check failed: type is not %d, current is %d", NST_META_REFE, stype);

                                                 /* 检测block 的data流 */
                                                 pbgz_stream_header header_data(curr_block->buffer + header_meta.get_bufferlen() + enclen_meta);
                                                 check_exit(header_data.check_mark(), ERR_INTERNEL, "header check failed: data id is not: %s", STREAM_ID);
                                                 sid_data = header_data.get_subid();
                                                 stype = (nstype)(header_data.get_type());
                                                 check_exit(stype == NST_DATA_REFE, ERR_INTERNEL,
                                                            "header check failed: type is not %d, current is %d", NST_DATA_REFE, stype);
                                                 check_exit(sid_data == sid_meta + 1, ERR_INTERNEL, "pair stream check failed: data id  (%ld) != meta id + 1 (%ld)", sid_data, sid_meta + 1);
  
                                                 /* 首先解析出当前block的meta信息，获取对应行的流信息和编码器等信息 */
                                                 cmeta.decoder(curr_block->buffer + header_meta.get_bufferlen(), header_meta.get_datalen(), block_meta);
                                                 if (block_meta["coder"]["magic"].asString() == "coder_ppmd")
                                                //  if (block_meta["coder"]["magic"].asString() == "coder_bwt_cm")
                                                 {
                                                     offset = block_meta["offset"].asInt64();
                                                     src_len = block_meta["srclen"].asInt();
                                                     dst_len = block_meta["dstlen"].asInt();

                                                     /* 指向reference当前块待解压的数据位置 */
                                                     curr_block->current_len = header_meta.get_bufferlen() + header_meta.get_datalen() + header_meta.get_bufferlen();

                                                     check_exit(curr_block->get_remain() == dst_len, ERR_INTERNEL,
                                                                "reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                                                block_meta["block"].asInt(), dst_len, curr_block->get_remain());
                                                    // std::cerr << block_meta << std::endl;

                                                    // {
                                                    //     coder_io data_io(curr_block->get_curr(), curr_block->get_remain());
                                                    //     coder_bwt_cm match_cm(&data_io);
                                                    //     len = match_cm.decode_line(refe_squash + offset, src_len, UINT8_MAX, false);
                                                    //     check_exit(src_len == len, ERR_INTERNEL,
                                                    //                "reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                                    //                block_meta["block"].asInt(), src_len, len);
                                                    // }

                                                    {
                                                        coder_io refe_io(curr_block->get_curr(), curr_block->get_remain());
                                                        refe_io.meta = block_meta["coder"];
                                                        coder_ppmd cppmd(&refe_io);
                                                        len = cppmd.decode(refe_squash + offset, src_len);
                                                        check_exit(src_len == len, ERR_INTERNEL,
                                                                   "reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                                                   block_meta["block"].asInt(), src_len, len);
                                                    }

                                                 }
                                                 else
                                                     check_exit(false, ERR_INTERNEL, "undefined coder %s", block_meta["coder"]["magic"].asString().c_str());

                                                 input_block->add(curr_block);
                                             }
                                         }));
        }
        else
        {
            input_block->take(current);
            {
                 /* reference block meta stream header */
                p = current->buffer;
                current->buffer_size = fileio_nokb.read(p, len_header);
                pbgz_stream_header header_meta(p);
                p += len_header;
                /* reference block meta stream data */
                len = header_meta.get_datalen();
                current->buffer_size += fileio_nokb.read(p, len);
                p += len;

                /* reference block data stream header */
                current->buffer_size += fileio_nokb.read(p, len_header);
                pbgz_stream_header header_data(p);
                p += len_header;
                /* reference block data stream data */
                len = header_data.get_datalen();
                current->buffer_size += fileio_nokb.read(p, len);
            }
            input_pool->add(current);
        }
    }

    for (auto &t : tpools)
    {
        input_block->take(current);
        input_pool->add(nullptr);
    }
    for (auto &t : tpools)
    {
        if (t.joinable())
            t.join();
    }

    // md5_init(&md5);
    // md5_write(&md5, (uint8_t *)(refgene->get_squash()), refgene->get_squashlen());
    // md5_final(&md5);
    // check_exit(md5.hexstr() == md5_packed, ERR_INTERNEL,
    //            "reference unpack failed in check md5: should be %s, actual is %s", md5_packed.c_str(), md5.hexstr().c_str());

    for (n = 0; n < pcnt; n++)
        delete block[n];
    delete input_pool;
    delete input_block;
}

 /* 解压时初始化参考基因组 */
bool engine::decoder_initialize_reference()
{
    int64_t file_len1, file_len2;
    Json::Value meta_refe;
    std::string fasta_name1, fasta_name2;
    std::string ni_name1, ni_name2;
    
    meta_refe = file_meta["refe"];
    fasta_name2 = meta_refe["fasta_name"].asString();
    file_len2 = meta_refe["fasta_len"].asInt64();
    ni_name2 = meta_refe["ni_name"].asString();

    if (meta_refe["file_offset"].isInt64())
        /*  压缩包里有pack reference，直接unpack*/
        reference_unpack();
    else
    {
        fasta_name1 = conf.get_refgene();
        if (fasta_name1.empty())
        { /* 没有指定reference文件 */
            fprintf(stderr, "need to specify the following FASTA file:\n\n");
            fprintf(stderr, "\t%-12s : %s\n", "File Name", fasta_name2.c_str());
            fprintf(stderr, "\t%-12s : %ld\n", "File Length", file_len2);
            fprintf(stderr, "\t%-12s : %s\n", "File MD5", meta_refe["fasta_md5"].asString().c_str());
            check_exit(false, ERR_INTERNEL, "reference file needs to be specified to complete decompression");
        }
        file_len1 = file_length(fasta_name1);
        file_abspath_filename(fasta_name1, fasta_name1);

        /* check whether fasta is matched */
        check_exit(fasta_name1 == fasta_name2, ERR_INTERNEL,
                   "initialize reference failed: used fasta %s, should be %s", fasta_name1.c_str(), fasta_name2.c_str());
        check_exit(file_len1 == file_len2, ERR_INTERNEL,
                   "initialize reference failed: used fasta file len %ld, should be %ld", file_len1, file_len2);

        reference refe_check(conf.get_refgene());
        refe_check.ni_from_reference(ni_name1);
        file_abspath_filename(ni_name1, ni_name1);
        check_exit(ni_name1 == ni_name2, ERR_INTERNEL,
                   "initialize reference failed: used ni file %s, should be %s", ni_name1.c_str(), ni_name2.c_str());
        /* matched, do make index */
        safe_new_class(reference, conf.get_refgene(), refgene);
        check_exit(refgene->initialize_squash_1(), ERR_INTERNEL, "initialize reference failed");
    }
    return true;
}
