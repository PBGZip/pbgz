#include "platform_compat.h"

#include "io.h"

io::io(const std::string &file, const iomode &mode, bool rw2cache)
{
    this->file_name = file;
    this->m = mode;
    this->parallel = manage::instance().get_zipinfo().get_threads();
    this->cache_bsplit_capacity = 0;
    this->cache_bsplit_current = 0;
    this->file_size = 0;
    this->eof = false;
    this->len_firstread = 0;
    this->is_gzipfile = false;
    this->fp = nullptr;
    this->fpGZ = nullptr;
    this->cache_bsplit = nullptr;
    this->isal_buffer_in = nullptr;
    this->isal_buffer_in = nullptr;
    this->isal_buffer_extra = nullptr;
    this->ringbuffer_guard = nullptr;
    this->ringbuffer = nullptr;
    this->ringbuffer_exit = false;
    this->ringbuffer_enable = rw2cache;
    this->ringbuffer_read_pos = 0;
    this->ringbuffer_write_pos = 0;
    this->support_simd = manage::instance().get_zipinfo().get_support_simd();

    check_exit(!file_name.empty(), ERR_FILE_READ, "file name is empty");
    check_exit(initialize(), ERR_INTERNEL, "io initialize fail.");
}

io::~io()
{
    if (m == mw || m == mwgz)
        flush();
    if (cache_bsplit)
        free(cache_bsplit);
    if (fp)
        fclose(fp);
    if (fpGZ)
        bgzf_close(fpGZ);
    if (ringbuffer_guard)
    {
        if (!ringbuffer_exit) /* 把环形缓存buffer中数据刷入磁盘，并结束guard线程 */
            ringbuffer_exit = true;
        get_event.wakeup();
        put_event.wakeup();
        ringbuffer_guard->join();
        delete ringbuffer_guard;
        ringbuffer_guard = nullptr;
    }
    if (ringbuffer)
        free(ringbuffer);
    if (isal_buffer_in)
        free(isal_buffer_in);
    if (isal_buffer_extra)
        free(isal_buffer_extra);
}

/*  根据文件类型初始化文件的读写 */
bool io::initialize()
{
    if (mr == m)
    { /* read */
        this->fp = fopen64(this->file_name.c_str(), "rb");
        check_exit(fp, ERR_FILE_READ, "open failed: %s", file_name.c_str());
        check_exit(fseeko64(fp, 0, SEEK_END) == 0, ERR_FILE_READ, "fseek failed: %s", file_name.c_str());
        this->file_size = ftello64(fp);
        rewind(fp);

        if (this->file_size > 3)
        { /* 判断是否是gz格式 */
            uint8_t tmp[3];
            fread(tmp, 3, 1, fp);
            this->is_gzipfile = is_gzip(tmp, 3);
        }
        fclose(fp);
        this->fp = nullptr;

        if (is_gzipfile)
        {
            if (support_simd)
            { /*  intel gz加速读 */
                int64_t len2read;
                fp = fopen64(this->file_name.c_str(), "rb");
                check_exit(fp, ERR_FILE_READ, "open for speed up read failed: %s", file_name.c_str());

                /*  intel gz解压一段gz数据时可能因为原始数据少了导致解压不了，所以需要每次补充一段原数据再去解压 */
                /* 建议 extra <= 20K，太小读文件io会影响性能，太大会触发realloc */
                isal_extra_each = 2048;
                isal_extra_len = (2 << 20);
                isal_in_len = (128 << 20);
                safe_alloc(isal_in_len, uint8_t, isal_buffer_in);
                safe_alloc(isal_extra_len, uint8_t, isal_buffer_extra);
                isal_remain_gzlen = this->file_size;
                isal_gzip_header_init(&isal_header);
                isal_inflate_init(&isal_state);
                isal_state.crc_flag = ISAL_GZIP_NO_HDR_VER;
                isal_state.next_in = this->isal_buffer_in;

                len2read = std::max<int64_t>((64 << 20), isal_in_len);
                len2read = std::min<int64_t>(isal_remain_gzlen, len2read);
                isal_state.avail_in = read_file(isal_state.next_in, len2read);
                len_firstread = isal_state.avail_in;
                isal_remain_gzlen -= isal_state.avail_in;
                check_exit(isal_read_gzip_header(&isal_state, &isal_header) == 0,
                           ERR_FILE_FORMAT, "found invalid gzip header in file %s", file_name.c_str());

                this->read_fun = &io::read_gzfile_fast;
            }
            else
            {                                              /* 普通gz读 */
                fpGZ = bgzf_open(file_name.c_str(), "rb"); /* 初始化对应的open */
                check(fpGZ, false, "bgzf open for read fail");
                if (fpGZ->is_gzip)
                    bgzf_mt(fpGZ, this->parallel, 256);

                this->read_fun = &io::read_gzfile; /* 初始化对应的read */
            }
        }
        else
        {
            fp = fopen64(this->file_name.c_str(), "rb");
            check_exit(fp, ERR_FILE_READ, "open for read failed: %s", file_name.c_str());

            this->read_fun = &io::read_file;
        }
    }
    else if (mwgz == m)
    { /* write gz */
        /*  即使支持硬件加速也不使用intel gz做压缩，因为测试结果是压缩时相比并发gz没有性能优势 */

        fpGZ = bgzf_open(file_name.c_str(), "wb");
        check(fpGZ, false, "bgzf open for write fail");

        fpGZ->compress_level = 4; /*  默认为4，有需要可以设置*/
        check(fpGZ->compress_level >= 0 && fpGZ->compress_level <= 9, false, "gz  compress level %d is invalid, should [0, 9]", fpGZ->compress_level);
        bgzf_mt(fpGZ, this->parallel, 256);

        this->write_fun = &io::write_gzfile;
    }
    else if (mw == m)
    { /* write */
        this->fp = fopen64(this->file_name.c_str(), "wb");
        check_exit(fp, ERR_FILE_WRITE, "open for write failed: %s", file_name.c_str());

        this->write_fun = &io::write_file;
    }
    else
    {
        check_exit(false, ERR_INTERNEL, "io initialize fail");
    }

    /* 初始化ring buffer */
    if (!ringbuffer_enable)
        return true;
    this->ringbuffer_capacity = std::max(BLOCK_SIZE, (256 << 20));
    this->ringbuffer_capacity = powerof2_proximal(this->ringbuffer_capacity);
    safe_alloc(ringbuffer_capacity, uint8_t, ringbuffer);

    this->ringbuffer_guard = new std::thread([&]()
                                             {
        int64_t remain, l, n, thresold = ringbuffer_capacity / 3; 
        
        if (m == mr) { /* 通过对应的文件read接口，从文件中往ring buffer补充数据 */
            for (;;) { 
                while(!ringbuffer_exit && ringbuffer_full()) /* ring buffer已满，需要等待把数据读走腾出空间 */
                    put_event.wait();
                remain = ringbuffer_space_left(); /* 全部补满 */

                l = min(remain, ringbuffer_capacity - (ringbuffer_write_pos & (ringbuffer_capacity - 1)));
                n = (this->*read_fun)(ringbuffer + (ringbuffer_write_pos & (ringbuffer_capacity - 1)), l);
                ringbuffer_write_pos += n;
                if (n == 0 || n != l) {
                    ringbuffer_exit = true;
                    get_event.wakeup();
                    return;
                }
                if (l < remain) {
                    n = (this->*read_fun)(ringbuffer, remain - l);
                    ringbuffer_write_pos += n;
                    if (n == 0 || n != remain - l) {
                        ringbuffer_exit = true;
                        get_event.wakeup();
                        return;
                    }
                }
                get_event.wakeup();
            }
        
        } else { /*  将buffer塞入ringbuffer，再由ring buffer调用对应的文件结果将buffer处理后写磁盘（譬如gz压缩) */
            for (;;) {
                while(!ringbuffer_exit && (ringbuffer_len() < thresold)) /* 累够阀值才写 */
                    get_event.wait();

                n = ringbuffer_len(); /* 全部写完 */
                l = min(n, ringbuffer_capacity - (ringbuffer_read_pos & (ringbuffer_capacity - 1)));
                (this->*write_fun)(ringbuffer + (ringbuffer_read_pos & (ringbuffer_capacity - 1)), l);
                (this->*write_fun)(ringbuffer, n - l);
                ringbuffer_read_pos += n;

                put_event.wakeup(); /* 如果往ring buffer在等待，则唤醒它继续写 */

                if (ringbuffer_exit && ringbuffer_empty())
                    break;
            }
        } });
    check_exit(ringbuffer_guard, ERR_MEM_NOENOUGH, "io guard thread initialize failed: no enough memory");

    return true;
}

/*  读一个完整的块，如读一个完整的fastq，最大读len数据，返回实际读到的长度，返回为0时表示已经读完了*/
int64_t io::read_one_block(block_rough &block)
{
    uint32_t type = BINARY;
    uint8_t *buffer;
    bool eof = false;
    int64_t len_remain, len_need_cache, lcnt, len;
    int64_t len_readed = 0, maxbase_len;

    block.reset();
    buffer = block.buffer;
    len = block.buffer_size;
    len_readed = block.current_len = 0;
    std::vector<uint32_t> &npos = block.npos;

    /*  step 1: 先把len长度的数据全部读出来 */
    if (cache_bsplit_current > 0)
    {
        check_exit(cache_bsplit_current <= block.buffer_size,
                   ERR_INTERNEL, "check failed in: cache size %ld <= %u", cache_bsplit_current, block.buffer_size);
        memcpy(block.buffer, this->cache_bsplit, this->cache_bsplit_current);
        len_readed += this->cache_bsplit_current;
        this->cache_bsplit_current = 0;
    }
    if (len_readed < block.buffer_size)
    {
        len = this->read(block.buffer + len_readed, block.buffer_size - len_readed);
        eof = len < (block.buffer_size - len_readed);
        len_readed += len;
    }

    /*  step 2: 再做块类型分析，先试图分析当前块是否是fastq */
    if ((maxbase_len = is_fastq(block.buffer, len_readed, npos, eof)) > 0)
    { /*  格式一：判断当前块是否fastq */
        len_remain = len_readed - npos[((npos.size() >> 2) << 2) - 1] - 1;

        if (len_remain)
        { /* 将本次读到的未完整fastq块之外的数据缓存到cache buffer */
            len_need_cache = len_remain << 2;
            if (this->cache_bsplit_capacity < len_need_cache)
            { /* 保证能缓存一个完整的4行的fastq小块 */
                safe_realloc(this->cache_bsplit_capacity, uint8_t, this->cache_bsplit, len_need_cache);
                this->cache_bsplit_capacity = len_need_cache;
            }
            memcpy(this->cache_bsplit, buffer + len_readed - len_remain, len_remain);
            this->cache_bsplit_current = len_remain;
            this->eof = false; /*  缓存buffer里有数据，不能算结束*/
        }
        if (maxbase_len > GENE3_MAX_BASE)
            type = BINARY;
        else if (maxbase_len > GENE2_MAX_BASE)
            type = FASTQ_GEN3;
        else
            type = FASTQ_GEN2;
        if (is_gzipfile)
            type |= GZIP;
        block.max_line_len = maxbase_len;
        block.btype = (blocktype)(type);
        block.current_len = len_readed - len_remain;
        lcnt = (npos.size() & 0x3); /*  行数是否以4对齐 */
        for (int i = 0; i < lcnt; i++)
            npos.pop_back();
        return block.current_len;

        /*  格式二：xxx */

        /* 未识别出格式: 整个块压缩 */
    }
    else
    {
        if (npos.size() >= 4) /*  已经解析完了4行但是没有被识别为fastq，那么就当作二进制类型了 */
            type = BINARY;
        if (is_gzipfile)
            type |= GZIP;
        block.btype = (blocktype)(type);
        block.current_len = len_readed;
        return block.current_len;
    }
}

/* 读指定长度的数据，返回实际读到的长度，为0表示结束（内部会处理出错的情形）*/
int64_t io::read(uint8_t *buffer, int64_t len)
{
    int64_t len_readed = 0, len_fromfile = 0;
    if (this->cache_bsplit_current)
    {
        /*  检查cache size 与 read size */
        check_exit(this->cache_bsplit_current <= len, ERR_INTERNEL, "innter error: cache size %ld > read size %ld", cache_bsplit_current, len);

        memcpy(buffer, this->cache_bsplit, this->cache_bsplit_current);
        len_readed += this->cache_bsplit_current;
        len_fromfile = read_inner(buffer + len_readed, len - len_readed);
        this->cache_bsplit_current = 0;
    }
    else
        len_fromfile = read_inner(buffer, len);

    check_exit(len_fromfile >= 0, ERR_FILE_READ, "read %s failed, get len %ld", file_name.c_str(), len_fromfile);
    len_readed += len_fromfile;
    return len_readed;
}

/*  读文件，读len长度，返回实际读到的长度 */
int64_t io::read_file(uint8_t *buffer, int64_t len)
{
    int64_t pos = ftello64(fp);
    if (fread(buffer, len, 1, this->fp) != 1)
    {
        if (feof(fp)) /*  遇到文件结束时，文件指针偏移重新读*/
        {
            check_exit(fseeko64(fp, pos, SEEK_SET) == 0, ERR_FILE_READ, "fseek failed: %s", file_name.c_str());

            int64_t remain = this->file_size - pos;
            if (remain == 0)
                return remain;
            check_exit(fread(buffer, remain, 1, this->fp) == 1, ERR_INTERNEL, "read remain data failed, len %ld, file %s ", remain, file_name.c_str());
            return remain;
        }
        else
        {
            perror("read error");
            manage::instance().exit(ERR_FILE_READ);
        }
    }
    return len;
}

/*  读gz格式文件，读len长度，返回实际读到的长度 */
int64_t io::read_gzfile(uint8_t *buffer, int64_t len)
{
    int64_t len_readed;
    check_exit((len_readed = bgzf_read(fpGZ, buffer, len)) >= 0, ERR_FILE_READ, "read gz file error: %s", file_name.c_str());
    return reinterpret_cast<int64_t>(len_readed);
}

/*  读gz格式文件，读len长度，返回实际读到的长度，指令加速方式读 */
int64_t io::read_gzfile_fast(uint8_t *buffer, int64_t len)
{
    int64_t len_unzipped, isal_res;
    uint8_t *isal_in, *isal_out;
    uint32_t isal_inlen, isal_outlen;
    int64_t len2read, len_in2cache;

    len_unzipped = 0;
    isal_res = ISAL_DECOMP_OK;
    isal_state.next_out = buffer;
    isal_state.avail_out = len;

    for (;;)
    {
        if (len_unzipped >= len || !isal_state.next_in) /*  文件已经解完或者已经解压到需要的长度 */
            break;

        /*  解压当前读到的gz原数据，如果出错则每次补充 isal_extra_len 原数据再解压 */
        for (;;)
        {
            /*  记录当前解压位置信息，包括原始数据的位置的长度，输出数据的位置和长度 */
            isal_in = isal_state.next_in;
            isal_inlen = isal_state.avail_in;
            isal_out = isal_state.next_out;
            isal_outlen = isal_state.avail_out;

            isal_res = isal_inflate(&isal_state);
            if (ISAL_DECOMP_OK == isal_res) /* 解压gz ok，退出*/
                break;
            check_exit(ISAL_END_INPUT == isal_res && 0 != this->isal_remain_gzlen,
                       ERR_FILE_FORMAT, "invalid gz file, errno %d", isal_res);
            memset(&this->isal_header, 0, sizeof(struct isal_gzip_header));
            isal_gzip_header_init(&(this->isal_header));
            memset(&this->isal_state, 0, sizeof(struct inflate_state));
            isal_inflate_init(&(this->isal_state));
            isal_state.crc_flag = ISAL_GZIP_NO_HDR_VER;

            /* 填充本次解压的gz原数据buffer: 未解压完的原数据+从文件读isal_extra_each原数据 */
            len2read = std::min<int64_t>(this->isal_remain_gzlen, isal_extra_each); /* 本次需要读的gz原始数据长度 */
            len_in2cache = isal_inlen + len2read;                                   /*  需要缓存的gz原始数据长度 */
            safe_realloc(isal_extra_len, uint8_t, isal_buffer_extra, len_in2cache);
            memcpy(isal_buffer_extra, isal_in, isal_inlen);
            read_file(isal_buffer_extra + isal_inlen, len2read);
            isal_remain_gzlen -= len2read;

            safe_realloc(isal_in_len, uint8_t, isal_buffer_in, len_in2cache);
            memcpy(isal_buffer_in, isal_buffer_extra, len_in2cache);

            /*  指向新的数据位置继续解压 */
            isal_state.avail_in = len_in2cache;
            isal_state.next_in = isal_buffer_in;
            isal_state.avail_out = isal_outlen;
            isal_state.next_out = isal_out;
        }

        len_unzipped = isal_state.next_out - buffer;
        if (len_unzipped != len)
        {
            bool filedone = true;
            if (isal_state.avail_in == 0)
            { /* buffer里没有原始gz数据了,读文件拿数据 */
                if (isal_remain_gzlen > 0)
                {
                    isal_state.next_in = isal_buffer_in;
                    len2read = std::min<int64_t>(this->isal_remain_gzlen, isal_in_len);
                    isal_state.avail_in = read_file(isal_state.next_in, len2read);
                    isal_remain_gzlen -= len2read;
                    filedone = (isal_state.avail_in == 0);
                }
            }
            else
            {
                for (;;)
                {
                    /*  记录当前解压位置信息，包括原始数据的位置的长度，输出数据的位置和长度 */
                    isal_in = isal_state.next_in;
                    isal_inlen = isal_state.avail_in;
                    isal_out = isal_state.next_out;
                    isal_outlen = isal_state.avail_out;

                    /* 当遇到有些特别特殊的，有多个头的gz文件时，intel gz需要重置header */
                    isal_inflate_reset(&isal_state);
                    isal_res = isal_read_gzip_header(&isal_state, &isal_header);
                    if (ISAL_DECOMP_OK == isal_res) /* 解压gz ok，退出*/
                        break;
                    check_exit(ISAL_END_INPUT == isal_res && 0 != this->isal_remain_gzlen,
                               ERR_FILE_FORMAT, "invalid gz header, errno %d", isal_res);

                    memset(&this->isal_header, 0, sizeof(struct isal_gzip_header));
                    isal_gzip_header_init(&(this->isal_header));
                    memset(&this->isal_state, 0, sizeof(struct inflate_state));
                    isal_inflate_init(&(this->isal_state));
                    isal_state.crc_flag = ISAL_GZIP_NO_HDR_VER;

                    /* 填充本次解压的gz原数据buffer: 未解压完的原数据+从文件读isal_extra_each原数据 */
                    len2read = std::min<int64_t>(this->isal_remain_gzlen, isal_extra_each); /* 本次需要读的gz原始数据长度 */
                    len_in2cache = isal_inlen + len2read;                                   /*  需要缓存的gz原始数据长度 */
                    safe_realloc(isal_extra_len, uint8_t, isal_buffer_extra, len_in2cache);
                    memcpy(isal_buffer_extra, isal_in, isal_inlen);
                    read_file(isal_buffer_extra + isal_inlen, len2read);
                    isal_remain_gzlen -= len2read;

                    safe_realloc(isal_in_len, uint8_t, isal_buffer_in, len_in2cache);
                    memcpy(isal_buffer_in, isal_buffer_extra, len_in2cache);

                    /*  指向新的数据位置继续解压 */
                    isal_state.avail_in = len_in2cache;
                    isal_state.next_in = isal_buffer_in;
                    isal_state.avail_out = isal_outlen;
                    isal_state.next_out = isal_out;
                }
                filedone = (isal_res != ISAL_DECOMP_OK);
            }

            if (filedone)
            {
                isal_state.next_in = nullptr;
                return len_unzipped;
            }
        }
    }
    return len_unzipped;
}

/* 写文件 */
void io::write_file(const uint8_t *buf, int64_t size)
{
    if (size <= 0)
        return;

    if (fwrite(buf, size, 1, fp) != 1)
    {
        perror("write error");
        manage::instance().exit(ERR_FILE_WRITE);
    }
}

/* 写gz文件 */
void io::write_gzfile(const uint8_t *buf, int64_t size)
{
    if (size <= 0)
        return;
    if (bgzf_write(fpGZ, buf, size) != size)
        manage::instance().exit(ERR_FILE_WRITE, "write gz file error");
}

/* 写gz文件，快速 */
void io::write_gzfile_fast(const uint8_t *buf, int64_t size)
{
    if (size <= 0)
        return;
    check_exit(false, ERR_INTERNEL, "undefined");
}

/*  写len长度数据，返回实际写的数据长度 */
int64_t io::write(const uint8_t *buffer, int64_t len, bool eof)
{
    write_inner(buffer, len, eof);
    return len;
}

/*  写操作时刷新文件io */
bool io::flush()
{
    bool ok;
    if (ringbuffer_guard)
    {
        if (!ringbuffer_exit) /* 把环形缓存buffer中数据刷入磁盘，并结束guard线程 */
            ringbuffer_exit = true;
        get_event.wakeup();
        put_event.wakeup();
        ringbuffer_guard->join();
        delete ringbuffer_guard;
        ringbuffer_guard = nullptr;
        ringbuffer_enable = false;
    }
    if (fp)
        check(fflush(fp) >= 0, false, "fflush file error");
    if (fpGZ)
        check(bgzf_flush(fpGZ) >= 0, false, "bgzf flush error");
    return true;
}

/*  写头信息 */
int64_t io::write_header(const uint8_t *header, int64_t len)
{
    return 0;
}

/* 将文件指针偏移到离头offset的位置 */
void io::fseek2pos(int64_t offset)
{
    /* 当ring buffer工作时，不支持fseek操作，因为需要把ring buffer清空，否则读到的可能是ringbuffer中的数据 */
    check_exit(!ringbuffer_enable, ERR_INTERNEL, "fseek don't support in current mode");
    check_exit(!is_gzipfile, ERR_INTERNEL, "gzip file don't support fseek!");
    check_exit(fseeko64(fp, offset, SEEK_SET) == 0, ERR_INTERNEL, "fseek to pos %ld failed", offset);
}

/*  初始化块切割时cache长度 */
void io::init_bsplit_cache(const uint32_t &len)
{
    this->cache_bsplit_capacity = len;
    safe_alloc(this->cache_bsplit_capacity, uint8_t, this->cache_bsplit);
}

/*  读模式时获取文件的size */
int64_t io::get_filesize()
{
    return this->file_size;
}

/* ring buffer */
/* 从ringbuffer里拿数据，返回0表示读结束，内部会处理读出错的情形 */
int64_t io::ringbuffer_get(uint8_t *buffer, int64_t len)
{
    uint8_t *p = buffer;
    int64_t l, n, needed = len;

    for (;;)
    {
        while (ringbuffer_empty() && !ringbuffer_exit) /* ring buffer为空且没有退出则等待数据 */
            get_event.wait();

        n = min(needed, ringbuffer_len());
        l = min(n, ringbuffer_capacity - (ringbuffer_read_pos & (ringbuffer_capacity - 1)));
        memcpy(p, ringbuffer + (ringbuffer_read_pos & (ringbuffer_capacity - 1)), l);
        memcpy(p + l, ringbuffer, n - l);
        ringbuffer_read_pos += n;
        needed -= n;
        p += n;

        put_event.wakeup(); /* 通知gurard线程已经读走一部分数据了，看是否需要往ring buffer中补充数据 */

        if (needed == 0 || (ringbuffer_exit && ringbuffer_empty()))
            break;
    }
    return len - needed;
}

/* 往ring buffer写数据 */
void io::ringbuffer_put(const uint8_t *buf, int64_t size)
{
    int64_t l, n, needed = size;

    for (;;)
    {
        while (ringbuffer_full()) /* ring buffer已满等待guard线程处理数据腾出空间 */
            put_event.wait();

        n = min(needed, ringbuffer_space_left());
        l = min(n, ringbuffer_capacity - (ringbuffer_write_pos & (ringbuffer_capacity - 1)));
        memcpy(ringbuffer + (ringbuffer_write_pos & (ringbuffer_capacity - 1)), buf, l);
        memcpy(ringbuffer, buf + l, n - l);

        ringbuffer_write_pos += n;
        needed -= n;

        get_event.wakeup(); /* 通知guard线程有数据来了，可以检查是否需要处理数据了，譬如将ring buffer中已有数据写磁盘 */
        if (needed == 0)
            break;
    }
}

/* 剩余空间长度 */
int64_t io::ringbuffer_space_left() const
{
    return ringbuffer_capacity - ringbuffer_write_pos + ringbuffer_read_pos;
}

/* ring buffer是否为空 */
bool io::ringbuffer_empty() const
{
    return ringbuffer_read_pos == ringbuffer_write_pos;
}

/* ring buffer 是否已满 */
bool io::ringbuffer_full() const
{
    return ringbuffer_capacity == ringbuffer_write_pos - ringbuffer_read_pos;
}

/* ring buffer当前数据长度 */
int64_t io::ringbuffer_len() const
{
    return ringbuffer_write_pos - ringbuffer_read_pos;
}

/* 内部读接口 */
int64_t io::read_inner(uint8_t *buffer, int64_t len)
{
    if (!ringbuffer_enable)
        return (this->*read_fun)(buffer, len);
    return ringbuffer_get(buffer, len);
}

/* 内部写接口 */
void io::write_inner(const uint8_t *buf, int64_t size, bool eof)
{
    if (!ringbuffer_enable)
        return (this->*write_fun)(buf, size);
    else
    {
        if (eof)
            ringbuffer_exit = true;
        return ringbuffer_put(buf, size);
    }
}