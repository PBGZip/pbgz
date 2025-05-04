#include "pbgz_file.h"
#include "coder/coder_json.h"
#include "city.h"

pbgz_file::pbgz_file(const std::string file_name, bool readonly)
{
    std::string mid_str = file_name, md5;
    this->file = file_name;
    this->readonly = readonly;
    this->file_offset = 0;
    this->streams_readed = 0;
    this->file_meta_offset = UINT64_MAX;
    this->meta_stream_offset.clear();

    fileio = new io(file_name, (readonly ? mr : mw), false);
    head = new nshead();
    check_exit(fileio && head, ERR_MEM_NOENOUGH,
               "pbgz file initialize failed: no enough memory");

    if (!readonly)
    { /* 写模式，压缩 */
        file_calc_md5(file_name, 1024, md5);
        mid_str += md5;
        mid2write = (uint32_t)CityHash64((const char *)mid_str.c_str(), mid_str.length());

        /* 写文件标识 */
        fileio->write((uint8_t *)(id_file.c_str()), id_file.length(), false);
        file_offset += id_file.length();

        /* 写版本号 */
        const char version[3] = {MAJOR, MINOR, PATCH};
        fileio->write((uint8_t *)version, 3, false);
        file_offset += 3;

        stream_start = file_offset;

        /* 第一条流类型为NST_META_FILE，记录了整个文件全局meta流的位置，有2份，一份用来备份 */
        sid2write = 0;
        uint8_t buffer[16];
        write_one_stream(buffer, 16, NST_META_FILE);
    }
    else
    { /* 读模式 ，解压 */

        /* 检查文件标识 */
        uint32_t id_file_len = id_file.size();
        uint8_t id_fromfile[id_file_len];
        check_exit(fileio->read(id_fromfile, id_file_len) == id_file_len, ERR_FILE_FORMAT,
                   "It's not pbgz format: %s", file.c_str());
        file_offset += id_file_len;

        /* 检查版本号 */
        uint8_t version[3];
        check_exit(fileio->read(version, 3) == 3, ERR_FILE_FORMAT,
                   "Bad pbgz format: cann't get version info: %s", file.c_str());
        // fprintf(stderr, "parsed  ---[version]---: %d.%d.%d\n", version[0], version[1], version[2]);
        if ((version[0] > MAJOR) ||
            (version[0] == MAJOR && version[1] > MINOR) ||
            (version[0] == MAJOR && version[1] == MINOR && version[2] > PATCH))
            check_exit(false, ERR_VERSION_MATCH,
                       "version is too old, at least: %d.%d.%d", version[0], version[1], version[2]);
        file_offset += 3;
    }
}

pbgz_file::~pbgz_file()
{
    if (fileio)
    {
        fileio->flush();
        delete fileio;
    }
    if (head)
        delete head;
}

/* 刷磁盘 */
void pbgz_file::flush()
{
    if (fileio)
        fileio->flush();
}

/* 设置主备两条meta流的文件偏移地址 */
void pbgz_file::set_meta_offset(int64_t offset_1, int64_t offset_2)
{
    int64_t meta_offset = stream_start + head->get_bufferlen();
    fileio->fseek2pos(meta_offset);
    fileio->write((uint8_t *)(&offset_1), 8, false);
    fileio->write((uint8_t *)(&offset_2), 8, true);
}

/* 解析全局的文件meta流 */
void pbgz_file::get_file_meta(Json::Value &file_meta)
{
    uint8_t buffer[16];
    uint8_t *data;
    nstype stream_type;
    int64_t len, offset_1, offset_2;
    int64_t len_filemeta_1, len_filemeta_2;
    coder_json cjson;
    Json::Value meta_2;

    /* 起一个关闭ringbuffer的文件操作，因为需要fseek读取文件对应位置的file meta流 */
    io fileio_nokb(file, mr, false);
    check_exit(readonly, ERR_INTERNEL, "parse file meta must work on readonly mode");

    check_exit(this->file_offset == (id_file.size() + 3), ERR_INTERNEL,
               "cannot get file meta info, current file offset %ld", file_offset);

    /* 读第一条流的头 */
    fileio_nokb.fseek2pos(file_offset);
    len = fileio_nokb.read(head->get_buffer(), head->get_bufferlen());
    check_exit(len == head->get_bufferlen(), ERR_FILE_FORMAT,
               "read first stream header failed, read len: %d, expected len %d", len, head->get_bufferlen());

    /* 校验该条流数据长度是否为 16 */
    check_exit(head->get_datalen() == 16, ERR_FILE_FORMAT,
               "check first stream data len failed: expected len %d, current is %d", 16, head->get_datalen());

    /* 读第一条流的数据部分 */
    len = fileio_nokb.read(buffer, head->get_datalen());
    check_exit(len == 16, ERR_FILE_FORMAT,
               "read first stream data failed, read len: %d, expected len %d", len, head->get_datalen());

    /* 读第一条file meta流，并检查合法性 */
    offset_1 = *((int64_t *)(buffer));
    offset_2 = *((int64_t *)(buffer + 8));
    fileio_nokb.fseek2pos(offset_1);
    len = fileio_nokb.read(head->get_buffer(), head->get_bufferlen());
    check_exit(len == head->get_bufferlen(), ERR_FILE_FORMAT,
               "read first file meta stream header failed, read len: %d, expected len %d", len, head->get_bufferlen());

    stream_type = (nstype)head->get_type();
    check_exit(stream_type == NST_META_FILE, ERR_FILE_FORMAT,
               "check first file meta stream header failed, expected type %d, current is %d", NST_META_FILE, stream_type);

    len_filemeta_1 = offset_2 - offset_1 - head->get_bufferlen();
    check_exit(head->get_datalen() == len_filemeta_1, ERR_FILE_FORMAT,
               "check first file meta stream data len failed: expected len %ld, current is %ld", len_filemeta_1, head->get_datalen());
    safe_alloc(len_filemeta_1, uint8_t, data);
    len = fileio_nokb.read(data, len_filemeta_1);
    check_exit(len == len_filemeta_1, ERR_FILE_FORMAT,
               "read first file meta stream data failed, read len: %d, expected len %d", len, len_filemeta_1);
    cjson.decoder(data, len_filemeta_1, file_meta);

    /* 读和解析第二条file meta流，并检查合法性 */
    fileio_nokb.fseek2pos(offset_2);
    len = fileio_nokb.read(head->get_buffer(), head->get_bufferlen());
    check_exit(len == head->get_bufferlen(), ERR_FILE_FORMAT,
               "read second file meta stream header failed, read len: %d, expected len %d", len, head->get_bufferlen());

    stream_type = (nstype)head->get_type();
    check_exit(stream_type == NST_META_FILE, ERR_FILE_FORMAT,
               "check second file meta stream header failed, expected type %d, current is %d", NST_META_FILE, stream_type);

    len_filemeta_2 = offset_2 - offset_1 - head->get_bufferlen();
    check_exit(head->get_datalen() == len_filemeta_2, ERR_FILE_FORMAT,
               "check second file meta stream data len failed: expected len %ld, current is %ld", len_filemeta_2, head->get_datalen());
    check_exit(len_filemeta_1 == len_filemeta_2, ERR_FILE_FORMAT,
               "check file meta stream data len failed, first len %ld != second len %ld", len_filemeta_1, len_filemeta_2);
    len = fileio_nokb.read(data, len_filemeta_2);
    check_exit(len == len_filemeta_2, ERR_FILE_FORMAT,
               "read second file meta stream data failed, read len: %d, expected len %d", len, len_filemeta_2);
    cjson.decoder(data, len_filemeta_2, meta_2);

    check_exit(fileio_nokb.read(buffer, 1) == 0, ERR_FILE_FORMAT,
               "check second file meta stream failed at end pos");
    check_exit(file_meta == meta_2, ERR_FILE_FORMAT, "check two file meta streams data failed: first != second");

    /*  预留：如果第一条流和第二条不想等，其中一条流可能是正常，或者更坏情形，即使file meta信息丢失了，还可以根据
      * 每条流做恢复
      */
    free(data);

    file_meta_offset = offset_1;
    meta_stream_offset = file_meta["stream_offsets"];

    /* 将第一条流读走，让文件指针指向数据块的地址，这里不做校验了前面已经做过了 */
    fileio->read(head->get_buffer(), head->get_bufferlen());
    fileio->read(buffer, head->get_datalen());
    file_offset += head->get_bufferlen() + head->get_datalen();
    check_exit(file_offset == meta_stream_offset[0]["e"].asInt(), ERR_FILE_FORMAT,
               "pbgz format error: expect offset %ld, current is %ld", file_offset, meta_stream_offset[0]["e"].asInt64());
    this->streams_readed++;
}

/* 获取文件当前的offset */
const int64_t pbgz_file::get_offset() const
{
    return this->file_offset;
}

/* 获取当前的sub id */
const int64_t pbgz_file::get_subid() const
{
    return this->sid2write;
}

/* 得到所有流的offset */
const Json::Value pbgz_file::get_stream_offset()
{
    return this->meta_stream_offset;
}

/* 从pbgz文件中读一对流：data meta + data，如果读出错了内部会处理，返回读到的流的长度 */
int64_t pbgz_file::read_pair_stream(block_rough &block)
{
    int64_t start, end, len2read;
    if (file_meta_offset == file_offset)
        return 0;
    start = meta_stream_offset[(Json::Value::ArrayIndex)(this->streams_readed)]["s"].asInt64();this->streams_readed++;
    end = meta_stream_offset[(Json::Value::ArrayIndex)(this->streams_readed)]["e"].asInt64();this->streams_readed++;
    len2read = end - start;
    fileio->read(block.buffer, len2read);
    this->file_offset += len2read;
    return len2read;
}

/* 写一条流，如果写出错内部会处理 */
void pbgz_file::write_one_stream(const uint8_t *buffer, uint32_t buffer_len, nstype stream_type, bool eof)
{
    Json::Value offset;
    offset["s"] = this->get_offset();
    check_exit(buffer_len != 0, ERR_INTERNEL, "write stream %ld failed: current stream data len is 0", sid2write);
    head->set(mid2write, sid2write++, stream_type, buffer_len);
    fileio->write(head->get_buffer(), head->get_bufferlen(), eof);
    fileio->write(buffer, buffer_len, eof);
    file_offset += head->get_bufferlen() + buffer_len;
    offset["e"] = this->get_offset();
    this->meta_stream_offset.append(offset);
}