#include "platform_compat.h"

#include "reference.h"
#include "file_op.h"
#include "vendor.h"
#include <algorithm>
#include "io.h"
#include "cfgpath/cfgpath.h"
#include "wait_event.h"
#include "pbgz_file.h"
#include "md5sum.h"
#include "coder/coder_json.h"
#include "actg.h"
#include "blake2.h"
#include "city.h"
#include <fstream>

reference::reference(const std::string &fasta)
{
    mbar = nullptr;
    reference_squash = nullptr;
    reference_squashlen = 0;
    hash_bucket_cnt = nullptr;
    hash_table_buffer = nullptr;
    reference_squash_matched = nullptr;
    reference_squash_matchedlen = 0;
    reference_file = fasta;
    parallel = manage::instance().get_zipinfo().get_threads();
    check_exit((basegroup_step & 0x3) == 0,
               ERR_INTERNEL, "invalid base group step %d", basegroup_step);
    /* 首先它必须为奇数，其次只处理了31的情形，改了mapping和reference的建表逻辑也得适配，当然这个长度没有必要去调整 */
    check_exit(basegroup_len == 31,
               ERR_INTERNEL, "invalid base group len %d", basegroup_len);
    check_exit(basegroup_len <= basegroup_step,
               ERR_INTERNEL, "reference basegroup step %d must be less than basegroup len %d", basegroup_step, basegroup_len);

    /* 检查fasta最大大小 need do */
}

reference::~reference()
{
    if (reference_squash)
        free(reference_squash);
    if (hash_table_buffer)
        free(hash_table_buffer);
    if (hash_bucket_cnt)
        free(hash_bucket_cnt);
    if (reference_squash_matched)
        free(reference_squash_matched);
    if (mbar)
        delete mbar;
}

/* 初始化参考基因squash buffer, 通过.ni索引文件 */
bool reference::initialize_squash_1()
{
    uint8_t buffer[4096], *p;
    Json::Value meta_ni;
    coder_json cmeta;
    std::string name_ni, path_ni;
    const std::string id_file = FILE_ID;
    int64_t len_data, len, len_meta, offset = 0;

    ni_from_reference(path_ni);
    std::cerr << "Reference index path: " << path_ni << "\n";

    /* 检查.ni文件的有效性，有效则直接从.ni文件读取 */
    if (!ni_valid(path_ni))
    {
        std::cerr << "do ni make..." << std::endl;
        check_warning(ni_make(path_ni), false, "make ni file failed: %s", path_ni.c_str());
        check_warning(ni_valid(path_ni), false, "check ni file failed after maked");
    }

    /* 初始化 squash buffer */
    io ni2read(path_ni, mr, false);
    len = ni2read.read(buffer, id_file.length());
    offset += len;
    check_warning(len == id_file.length() && !memcmp(buffer, id_file.c_str(), id_file.length()),
                  false, "ni file format error in file magic header");

    len = ni2read.read(buffer, 3);
    offset += len;
    check_warning(len == 3, false, "ni file version check failed");
    len = ni2read.read(buffer, 2);
    offset += len;
    check_warning(len == 2 && !memcmp(buffer, "ni", 2), false, "ni file format error in magic");

    len = ni2read.read((uint8_t *)(&len_data), sizeof(len_data));
    offset += len;
    check_warning(len == sizeof(len_data), false, "ni file check failed in meta offset");

    reference_squashlen = len_data;
    safe_alloc(reference_squashlen, uint8_t, reference_squash);

    len = ni2read.read(reference_squash, reference_squashlen);
    offset += len;
    check_warning(len == reference_squashlen, false, "ni file check falied in data");

    p = buffer;
    len_meta = ni2read.get_filesize() - offset;
    if (len_meta > 4096)
        safe_alloc(len_meta, uint8_t, p);
    len = ni2read.read(p, len_meta);
    check_warning(len == len_meta, false, "ni file check failed in meta");
    cmeta.decoder(p, len_meta, meta_ni);
    check_warning(ni2read.read(p, 1) == 0, false, "ni file check failed in eof");
    if (p != buffer)
        free(p);

    return true;
}

/* 初始化参考基因squash buffer, 通过stream */
uint8_t * reference::initialize_squash_2(int64_t squash_len)
{
    reference_squashlen = squash_len;
    safe_alloc(reference_squashlen, uint8_t, reference_squash);
    return reference_squash;
}

/* .ni文件是否有效 */
bool reference::ni_valid(const std::string &ni)
{
    Json::Value niconf;
    coder_json cmeta;
    uint8_t *buffer;
    int32_t n;
    int64_t file_size;
    char cfgdir[MAX_PATH];
    std::string conf, refename, ni_name;

    /*  检查.conf文件是否存在 */
    get_user_config_folder(cfgdir, sizeof(cfgdir), FILE_ID);
    if (cfgdir[0] == 0)
        create_dir(cfgdir);
    conf = cfgdir;
    conf += ".conf";
    if (!file_exists(conf) || !file_exists(ni))
        return false;

    /* 读conf文件 */
    io conf2read(conf, mr, false);
    file_size = conf2read.get_filesize();
    safe_alloc(file_size, uint8_t, buffer);
    check_exit(conf2read.read(buffer, file_size) == file_size,
               ERR_FILE_READ, "conf file read failed: %s", conf.c_str());
    cmeta.decoder(buffer, file_size, ref2ni_cache);
    free(buffer);
    file_abspath_filename(reference_file, refename);
    if (!ref2ni_cache[refename.c_str()].isArray()) /* 这里存数组是考虑同名的reference，但是内容不相同的情形 */
        return false;

    io ni2read(ni, mr, false);
    file_size = ni2read.get_filesize();
    file_abspath_filename(ni, ni_name);
    for (n = 0; n < ref2ni_cache[refename.c_str()].size(); n++)
    {
        niconf = ref2ni_cache[refename.c_str()][n];
        if (!niconf["ni_name"].isNull() && niconf["ni_name"] == ni_name)
        {
            /* conf中有记录当前ni文件信息*/
            /* 检查最后修改时间 */
            if (file_last_mtime(ni) != niconf["ni_mtime"].asInt64())
                return false;
            /* 检查文件长度 */
            if (file_size != niconf["ni_fsize"].asInt64())
                return false;
            return true;
        }
    }
    return false;
}

/* 创建.ni索引文件 */
bool reference::ni_make(const std::string &ni)
{
    timer tcost(true);
    Json::Value meta;
    int64_t n, refelen, left, docnt;
    int64_t offset = 0, wlen = 0;
    uint8_t *buffer, *s, *e;
    // uint8_t *squash_buffer;
    uint8_t last[4], ch;
    uint32_t last_len = 0;
    uint32_t buffer_len, len, len_actg, len_squash, l, l4align;
    const std::string id_file = FILE_ID;
    const char version[3] = {MAJOR, MINOR, PATCH};
    MD5_CONTEXT md5;
    coder_json cmeta;
    std::string meta2write;
    std::string refe_md5, refename;
    std::string refe = reference_file;

    /* 检查是否有程序已经在make ni */
    int32_t fd;
    bool done = false, skip = false;
    std::string file_lock = ni;
    file_lock += ".lock";
    std::thread *progress = nullptr;
    for (;;)
    {
        if ((fd = acquire_lock(file_lock)) < 0)
        {
            skip = true; /* 检测到有程序在make该ni时，当前程序直接等待ni文件制作完成直接使用 */
            if (!progress)
            {
                progress = new std::thread([&done]()
                                           {
                                               const std::string prompt = "......";
                                               int32_t cnt = 0, mask = powerof2_proximal(prompt.length()) - 1;
                                               for (;;)
                                               {
                                                   if (done)
                                                       break;
                                                   ++cnt;
                                                   cnt &= mask;
                                                   fprintf(stderr, "\33[2K\ranother program is making ni file, waiting %s", std::string(prompt.c_str(), cnt).c_str());
                                                   usleep(500000);
                                               }
                                           });
            }
            usleep(10000);
            continue;
        }
        done = true;
        if (progress)
            progress->join();
        if (skip)
        {
            fprintf(stderr, "\n");
            return true;
        }
        break;
    }

    std::thread calc_md5 = std::thread([&refe, &refe_md5, &refelen]()
                                       {
                                           uint8_t *buffer;
                                           uint32_t buffer_len = 4096, len;
                                           refelen = 0;
                                           MD5_CONTEXT md5;
                                           safe_alloc(buffer_len, uint8_t, buffer);
                                           md5_init(&md5);
                                           io io2read(refe, mr, false);
                                           for (;;)
                                           {
                                               len = io2read.read(buffer, buffer_len);
                                               refelen += len;
                                               md5_write(&md5, buffer, len);
                                               if (len != buffer_len)
                                                   break;
                                           }
                                           md5_final(&md5);
                                           refe_md5 = md5.hexstr();
                                       });

    buffer_len = (5 << 20);
    safe_alloc(buffer_len, uint8_t, buffer);
    len_actg = (4 << 20);
    // squash_buffer = buffer + len_actg;
    // io ref2read(reference_file, mr, false); /*  单线程处理时不开起读写cache，会影响性能 */
    io ni2write(ni, mw, false);
    md5_init(&md5);

    ni2write.write((uint8_t *)(id_file.c_str()), id_file.length(), false);
    ni2write.write((uint8_t *)version, 3, false);
    ni2write.write((uint8_t *)"ni", 2, false);
    offset += id_file.length() + 3 + 2;
    ni2write.write((uint8_t *)(&offset), sizeof(offset), false);

    /* 串行处理就行，基本就是读文件的时间, base_squash时间可以忽略 */
    std::ifstream file(reference_file.c_str());
    std::string line;

    uint8_t cache_actg[4];
    uint32_t cache_len = 0;
    uint8_t *squash_buffer;
    
    uint32_t squash_bufferlen = 1024 >> 2;
    safe_alloc(1024, uint8_t, squash_buffer);

    for (;;)
    {
        if (!std::getline(file, line))
            break;
        if (line[0] == '>')
            continue;
        
        if (cache_len + line.length() < 4) {
            memcpy(cache_actg + cache_len, line.c_str(), line.length());
            cache_len += line.length();
        } else {
            /* 先处理cache */
            docnt = 4 - cache_len;
            memcpy(cache_actg + cache_len, line.c_str(), docnt);
            actg_squash((uint8_t *)cache_actg, 4, &ch);
            ni2write.write(&ch, 1, false);
            wlen += 1;
            md5_write(&md5, &ch, 1);
            cache_len = 0;

            left = line.length() - docnt;
            if (left < 4) {
                /* 剩余的不够cache */
                memcpy(cache_actg + cache_len, line.c_str() + docnt, left);
                cache_len += left;
            } else {
                /* 剩余的够cache */
                l4align = left >> 2 << 2;
                safe_realloc(squash_bufferlen, uint8_t, squash_buffer, l4align >> 2);
                len_squash = actg_squash((const uint8_t *)(line.c_str() + docnt), l4align, squash_buffer);
                ni2write.write(squash_buffer, len_squash, false);
                wlen += len_squash;
                md5_write(&md5, squash_buffer, len_squash);

                /* 处理未处理的未以4对齐的字节 */
                docnt += l4align;
                left = line.length() - docnt;
                if (left > 0) {
                    memcpy(cache_actg + cache_len, line.c_str() + docnt, left);
                    cache_len += left;
                }
            }
        }
    }

    ni2write.fseek2pos(offset);
    ni2write.write((uint8_t *)(&wlen), sizeof(wlen), false);
    ni2write.fseek2pos(offset + sizeof(offset) + wlen);

    md5_final(&md5);
    calc_md5.join();

    file_abspath_filename(reference_file, refename);
    meta["refe_name"] = refename;
    meta["refe_len"] = Json::Value::UInt64(refelen);
    meta["refe_orgfile_md5"] = refe_md5; /* reference文件的原始md5，即直接读未解开 */
    meta["ni_data_md5"] = md5.hexstr();  /* ni文件的数据内容的md5 */
    fasta_len = refelen;
    fasta_md5 = refe_md5;

    cmeta.encoder(meta, meta2write);
    ni2write.write((uint8_t *)(meta2write.c_str()), meta2write.length(), true);
    ni2write.flush();

    { /* 将当前ni信息写conf */
        Json::Value niconf, nicurr;
        std::string conf, ni_name, out;
        int64_t file_size;
        char cfgdir[MAX_PATH];

        get_user_config_folder(cfgdir, sizeof(cfgdir), FILE_ID);
        if (cfgdir[0] == 0)
            create_dir(cfgdir);
        conf = cfgdir;
        conf += ".conf";

        io ni2read(ni, mr, false);
        file_abspath_filename(ni, ni_name);
        nicurr["ni_name"] = ni_name;
        nicurr["ni_mtime"] = (Json::Value::Int64)(file_last_mtime(ni));
        nicurr["ni_fsize"] = (Json::Value::Int64)(ni2read.get_filesize());

        if (file_exists(conf))
        {
            io conf2read(conf, mr, false);
            file_size = conf2read.get_filesize();
            safe_realloc(buffer_len, uint8_t, buffer, file_size);
            check_exit(conf2read.read(buffer, file_size) == file_size,
                        ERR_FILE_READ, "conf file read failed: %s", conf.c_str());
            cmeta.decoder(buffer, file_size, ref2ni_cache);

            if (!ref2ni_cache[refename.c_str()].isArray()) /* 没有conf信息则新增 */
                ref2ni_cache[refename.c_str()].append(nicurr);
            else
            {
                for (n = 0; n < ref2ni_cache[refename.c_str()].size(); n++)
                {
                    niconf = ref2ni_cache[refename.c_str()][(int32_t)n];
                    if (!niconf["ni_name"].isNull() && niconf["ni_name"] == ni_name)
                    { /* 存在则更新 */
                        ref2ni_cache[refename.c_str()][(int32_t)n]["ni_mtime"] = nicurr["ni_mtime"];
                        ref2ni_cache[refename.c_str()][(int32_t)n]["ni_fsize"] = nicurr["ni_fsize"];
                        break;
                    }
                }
                if (n >= ref2ni_cache[refename.c_str()].size())
                    ref2ni_cache[refename.c_str()].append(nicurr);
            }
        }
        else
            ref2ni_cache[refename.c_str()].append(nicurr);

        io conf2write(conf, mw, false);
        cmeta.encoder(ref2ni_cache, out);
        check_exit(conf2write.write((uint8_t *)(out.c_str()), out.length(), true) == out.length(),
                    ERR_FILE_WRITE, "conf file write failed: %s", conf.c_str());
}

    free(buffer);
    free(squash_buffer);
    return true;
}

// bool reference::ni_make(const std::string &ni)
// {
//     timer tcost(true);
//     Json::Value meta;
//     int64_t n, refelen;
//     int64_t offset = 0, wlen = 0;
//     uint8_t *buffer, *s, *e;
//     uint8_t *squash_buffer;
//     uint8_t last[4], ch;
//     uint32_t last_len = 0;
//     uint32_t buffer_len, len, len_actg, len_squash, l, l4align;
//     const std::string id_file = FILE_ID;
//     const char version[3] = {MAJOR, MINOR, PATCH};
//     MD5_CONTEXT md5;
//     coder_json cmeta;
//     std::string meta2write;
//     std::string refe_md5, refename;
//     std::string refe = reference_file;

//     /* 检查是否有程序已经在make ni */
//     int32_t fd;
//     bool done = false, skip = false;
//     std::string file_lock = ni;
//     file_lock += ".lock";
//     std::thread *progress = nullptr;
//     for (;;)
//     {
//         if ((fd = acquire_lock(file_lock)) < 0)
//         {
//             skip = true; /* 检测到有程序在make该ni时，当前程序直接等待ni文件制作完成直接使用 */
//             if (!progress)
//             {
//                 progress = new std::thread([&done]()
//                                            {
//                                                const std::string prompt = "......";
//                                                int32_t cnt = 0, mask = powerof2_proximal(prompt.length()) - 1;
//                                                for (;;)
//                                                {
//                                                    if (done)
//                                                        break;
//                                                    ++cnt;
//                                                    cnt &= mask;
//                                                    fprintf(stderr, "\33[2K\ranother program is making ni file, waiting %s", std::string(prompt.c_str(), cnt).c_str());
//                                                    usleep(500000);
//                                                }
//                                            });
//             }
//             usleep(10000);
//             continue;
//         }
//         done = true;
//         if (progress)
//             progress->join();
//         if (skip)
//         {
//             fprintf(stderr, "\n");
//             return true;
//         }
//         break;
//     }

//     std::thread calc_md5 = std::thread([&refe, &refe_md5, &refelen]()
//                                        {
//                                            uint8_t *buffer;
//                                            uint32_t buffer_len = 4096, len;
//                                            refelen = 0;
//                                            MD5_CONTEXT md5;
//                                            safe_alloc(buffer_len, uint8_t, buffer);
//                                            md5_init(&md5);
//                                            io io2read(refe, mr, false);
//                                            for (;;)
//                                            {
//                                                len = io2read.read(buffer, buffer_len);
//                                                refelen += len;
//                                                md5_write(&md5, buffer, len);
//                                                if (len != buffer_len)
//                                                    break;
//                                            }
//                                            md5_final(&md5);
//                                            refe_md5 = md5.hexstr();
//                                        });

//     buffer_len = (5 << 20);
//     safe_alloc(buffer_len, uint8_t, buffer);
//     len_actg = (4 << 20);
//     squash_buffer = buffer + len_actg;
//     io ref2read(reference_file, mr, false); /*  单线程处理时不开起读写cache，会影响性能 */
//     io ni2write(ni, mw, false);
//     md5_init(&md5);

//     ni2write.write((uint8_t *)(id_file.c_str()), id_file.length(), false);
//     ni2write.write((uint8_t *)version, 3, false);
//     ni2write.write((uint8_t *)"ni", 2, false);
//     offset += id_file.length() + 3 + 2;
//     ni2write.write((uint8_t *)(&offset), sizeof(offset), false);

//     /* 串行处理就行，基本就是读文件的时间, base_squash时间可以忽略 */
//     for (;;)
//     {
//         len = ref2read.read(buffer, len_actg);

//         s = buffer;
//         for (n = 0; n < len; n++)
//         {
//             if (*(buffer + n) == '\n' || n + 1 == len) /* 是新的一行或者本次读结束时 */
//             {
//                 if ((*s) == '>') {
//                       s = (*(buffer + n) == '\n') ? (buffer + n + 1) : nullptr;
//                     continue; /* 非碱基开头头丢弃 */
//                 }
//                 // switch (*s)
//                 // {
//                 // case 'a':
//                 //     break;
//                 // case 'c':
//                 //     break;
//                 // case 't':
//                 //     break;
//                 // case 'g':
//                 //     break;
//                 // case 'n':
//                 //     break;
//                 // case 'A':
//                 //     break;
//                 // case 'C':
//                 //     break;
//                 // case 'T':
//                 //     break;
//                 // case 'G':
//                 //     break;
//                 // case 'N':
//                 //     break;
//                 // default:
//                 //     s = (*(buffer + n) == '\n') ? (buffer + n + 1) : nullptr;
//                 //     continue; /* 非碱基开头头丢弃 */
//                 //     break;
//                 // }

//                 e = buffer + n + (*(buffer + n) != '\n'); /* 当为本次读结束处理时需要补1 */
//                 if (last_len > 0)
//                 { /* 补齐4个碱基 */
//                     l = e - s;
//                     if (l + last_len < 4)
//                     { /* 需要考虑本行很短的特殊情形 */
//                         memcpy(last + last_len, s, l);
//                         last_len += l;
//                         s = e + 1;
//                         continue;
//                     }
//                     else if (l + last_len == 4)
//                     {
//                         memcpy(last + last_len, s, l);
//                         actg_squash((uint8_t *)last, 4, &ch);
//                         ni2write.write(&ch, 1, false);
//                         wlen += 1;
//                         md5_write(&md5, &ch, 1);
//                         s = e + 1;
//                         continue;
//                     }
//                     else
//                     {
//                         memcpy(last + last_len, s, 4 - last_len);
//                         actg_squash((uint8_t *)last, 4, &ch);
//                         ni2write.write(&ch, 1, false);
//                         wlen += 1;
//                         md5_write(&md5, &ch, 1);
//                         s += 4 - last_len;
//                     }
//                 }
//                 l = e - s; /* 当前行未处理的碱基长度 */
//                 l4align = (l >> 2) << 2;
//                 if (l4align > 0)
//                 { /* 4对齐编码碱基 */
//                     len_squash = actg_squash(s, l4align, squash_buffer);
//                     ni2write.write(squash_buffer, len_squash, false);
//                     wlen += len_squash;
//                     md5_write(&md5, squash_buffer, len_squash);
//                 }
//                 /* 保存本行未编码的碱基 */
//                 last_len = l - l4align;
//                 if (last_len > 0)
//                     memcpy(last, s + l4align, last_len);

//                 s = (*(buffer + n) == '\n') ? (buffer + n + 1) : nullptr;
//             }
//         }

//         if (len != len_actg)
//         {
//             if (last_len > 0)
//             { /* 如果整个参考基因组都处理完了，把尾巴未对齐部分处理下 */
//                 const std::string actg2align = "ACTG";
//                 memcpy(last + last_len, actg2align.c_str(), 4 - last_len);
//                 actg_squash((uint8_t *)last, 4, &ch);
//                 ni2write.write(&ch, 1, true);
//                 wlen += 1;
//                 md5_write(&md5, &ch, 1);
//             }
//             ni2write.fseek2pos(offset);
//             ni2write.write((uint8_t *)(&wlen), sizeof(wlen), false);
//             ni2write.fseek2pos(offset + sizeof(offset) + wlen);
//             break;
//         }
//     }
//     md5_final(&md5);
//     calc_md5.join();

//     file_abspath_filename(reference_file, refename);
//     meta["refe_name"] = refename;
//     meta["refe_len"] = Json::Value::UInt64(refelen);
//     meta["refe_orgfile_md5"] = refe_md5; /* reference文件的原始md5，即直接读未解开 */
//     meta["ni_data_md5"] = md5.hexstr();  /* ni文件的数据内容的md5 */
//     fasta_len = refelen;
//     fasta_md5 = refe_md5;

//     cmeta.encoder(meta, meta2write);
//     ni2write.write((uint8_t *)(meta2write.c_str()), meta2write.length(), true);
//     ni2write.flush();

//     { /* 将当前ni信息写conf */
//         Json::Value niconf, nicurr;
//         std::string conf, ni_name, out;
//         int64_t file_size;
//         char cfgdir[MAX_PATH];

//         get_user_config_folder(cfgdir, sizeof(cfgdir), FILE_ID);
//         if (cfgdir[0] == 0)
//             create_dir(cfgdir);
//         conf = cfgdir;
//         conf += ".conf";

//         io ni2read(ni, mr, false);
//         file_abspath_filename(ni, ni_name);
//         nicurr["ni_name"] = ni_name;
//         nicurr["ni_mtime"] = (Json::Value::Int64)(file_last_mtime(ni));
//         nicurr["ni_fsize"] = (Json::Value::Int64)(ni2read.get_filesize());

//         if (file_exists(conf))
//         {
//             io conf2read(conf, mr, false);
//             file_size = conf2read.get_filesize();
//             safe_realloc(buffer_len, uint8_t, buffer, file_size);
//             check_exit(conf2read.read(buffer, file_size) == file_size,
//                        ERR_FILE_READ, "conf file read failed: %s", conf.c_str());
//             cmeta.decoder(buffer, file_size, ref2ni_cache);

//             if (!ref2ni_cache[refename.c_str()].isArray()) /* 没有conf信息则新增 */
//                 ref2ni_cache[refename.c_str()].append(nicurr);
//             else
//             {
//                 for (n = 0; n < ref2ni_cache[refename.c_str()].size(); n++)
//                 {
//                     niconf = ref2ni_cache[refename.c_str()][(int32_t)n];
//                     if (!niconf["ni_name"].isNull() && niconf["ni_name"] == ni_name)
//                     { /* 存在则更新 */
//                         ref2ni_cache[refename.c_str()][(int32_t)n]["ni_mtime"] = nicurr["ni_mtime"];
//                         ref2ni_cache[refename.c_str()][(int32_t)n]["ni_fsize"] = nicurr["ni_fsize"];
//                         break;
//                     }
//                 }
//                 if (n >= ref2ni_cache[refename.c_str()].size())
//                     ref2ni_cache[refename.c_str()].append(nicurr);
//             }
//         }
//         else
//             ref2ni_cache[refename.c_str()].append(nicurr);

//         io conf2write(conf, mw, false);
//         cmeta.encoder(ref2ni_cache, out);
//         check_exit(conf2write.write((uint8_t *)(out.c_str()), out.length(), true) == out.length(),
//                    ERR_FILE_WRITE, "conf file write failed: %s", conf.c_str());
//     }

//     free(buffer);
//     return true;
// }

/* 得到reference对应的ni文件名 */
void reference::ni_from_reference(std::string &ni)
{
    FILE *fp;
    int64_t file_size, read_each = 1024;
    int64_t offset, offset_each, len, n;
    char cfgdir[MAX_PATH];
    uint8_t buffer[(read_each << 2) + sizeof(file_size)];
    std::string name_ni;
    MD5_CONTEXT md5;

    md5_init(&md5);
    file_abspath_filename(reference_file, name_ni);
    get_user_data_folder(cfgdir, sizeof(cfgdir), FILE_ID);
    if (cfgdir[0] == 0)
        create_dir(cfgdir);
    path_ni = cfgdir;

    fp = fopen(reference_file.c_str(), "rb");
    check_exit(fp, ERR_FILE_READ, "reference file read failed: %s", reference_file.c_str());
    check_exit(fseeko64(fp, 0, SEEK_END) == 0, ERR_FILE_READ, "fseek failed: %s", reference_file.c_str());
    file_size = ftello64(fp);
    rewind(fp);

    offset_each = file_size >> 2;
    if (offset_each < read_each)
    {
        uint8_t buffer[file_size + sizeof(file_size)];
        check_exit(fread(buffer, file_size, 1, fp) == 1, ERR_FILE_READ, "reference file read failed: %s", reference_file.c_str());
        *((int64_t *)(buffer + file_size)) = file_size;
        md5_write(&md5, buffer, file_size + sizeof(file_size));
    }
    else
    {
        offset = len = 0;
        for (n = 0; n < 4; n++)
        {
            check_exit(fseeko64(fp, offset, SEEK_SET) == 0, ERR_FILE_READ, "reference file fseek failed: %s", reference_file.c_str());
            check_exit(fread(buffer + len, read_each, 1, fp) == 1, ERR_FILE_READ, "reference file read failed: %s", reference_file.c_str());
            offset += offset_each;
            len += read_each;
        }
        *((int64_t *)(buffer + len)) = file_size;
        md5_write(&md5, buffer, (read_each << 2) + sizeof(file_size));
    }
    fclose(fp);

    md5_final(&md5);
    name_ni += ".";
    name_ni += std::string(md5.hexstr().c_str(), 8); /* 取md5的前8位 */
    name_ni += ".ni";
    path_ni += name_ni;
    ni = path_ni;
}

/* 对外接口，创建索引表 */
bool reference::make_index()
{
    timer cost_ms(true);
    int64_t support_max = ((int64_t)2 << 30) * basegroup_step;
    bg_hash *hash;
    htable hash_table;
    std::string tips;
    /* 记录每个bucket下一个hash值写的位置，该位置为相对hash_buff起始位置的偏移，这样是为了并发写 */
    uint32_t *hash_bucket_curpos;

    safe_alloc(hash_buckets, uint32_t, hash_bucket_curpos);
    check_warning(initialize_squash_1(), false, "initialize reference failed");
    if ((reference_squashlen << 2) >= support_max)
        check_exit(false, ERR_INTERNEL, "Reference max support %lu(M), current size %ld(M)\n",
                   support_max / 1024 / 1024, (reference_squashlen << 2) / 1024 / 1024);
    else
        fprintf(stderr, "Reference max support %lu(M), current size %ld(M)\n",
                support_max / 1024 / 1024, (reference_squashlen << 2) / 1024 / 1024);

    cost_ms.reset();
    make_step1_fetchBG(hash);
    // printf("\n\t >>> step 1 cost ms: %d\n", cost_ms.elapsed()); cost_ms.reset();
    make_step2_calcHT(hash_table);
    // printf("\n\t >>> step 2 cost ms: %d\n", cost_ms.elapsed()); cost_ms.reset();
    make_step3_initHT(hash_table, hash_bucket_curpos);
    // printf("\n\t >>> step 3 cost ms: %d\n", cost_ms.elapsed()); cost_ms.reset();
    make_step4_buildHT(hash, hash_bucket_curpos);
    //   dump_hash_table();
    // printf("\n\t >>> step 4 cost ms: %d\n", cost_ms.elapsed()); cost_ms.reset();
    make_step5_sortHT();
    //  dump_hash_table();
    // printf("\n\t >>> step 5 cost ms: %d\n", cost_ms.elapsed()); cost_ms.reset();

    tips = "elapsed ms: ";
    tips += std::to_string((int64_t)(cost_ms.elapsed()));
    mbar->done(tips);
   
    free(hash);
    free(hash_bucket_curpos);

    reference_squash_matchedlen = reference_squashlen; /* 一个字节表示一个squash字节是否有matched */
    safe_alloc(reference_squash_matchedlen, uint8_t, reference_squash_matched);
    return true;
}

/* 制作索引表：从参考基因组中拿碱基*/
void reference::make_step1_fetchBG(bg_hash *&hash)
{
    uint8_t *p; 
    uint32_t *hb_cnt;
    bg_hash *bhash;
    std::vector<std::thread> tpools;
    int64_t each, current, remain, total, offset_start;
    int64_t n, pcnt = this->parallel;
    int64_t *pnn = nullptr;
    spinlock *slocks;

    const uint32_t bg_step = basegroup_step;
    const uint32_t bg_len = basegroup_len;
    const uint32_t *actg_stretch_tab = actg_stretch;
    const int32_t hbuckets = hash_buckets;
    const int32_t hmask = hash_mask;

    total = (reference_squashlen << 2) / basegroup_step;
    each = total / pcnt;
    remain = (total - (each * pcnt));
    p = this->reference_squash;
    safe_alloc(total, bg_hash, hash);
    safe_alloc(hash_buckets, uint32_t, hash_bucket_cnt);
    safe_alloc(hash_buckets, spinlock, slocks);
    hb_cnt = hash_bucket_cnt;
    bhash = hash;
    bg_hash *bhash_start = hash;

    for (n = 0; n < pcnt; n++)
    {
        current = (n + 1 == pcnt) ? (each + remain) : each; /* 当前处理的key数 */
        if ((n + 1) == pcnt) {
            mbar_total = 0;
            mbar_current = current;
            mbar = new guard_bar(mbar_current, &mbar_total, "Building index from reference");
            mbar->start();
            pnn = &mbar_total;
        }

        offset_start = bhash - bhash_start;

        tpools.push_back(std::thread([p, current, offset_start, hbuckets, hmask, &bhash_start,
                                      &bg_step, &bg_len, &actg_stretch_tab, &hb_cnt, &slocks, pnn]()
                                     {
                                         uint32_t len_hash = 0, len_bucket = 0;
                                         int64_t n = 0, m = 0;
                                         uint8_t current_cnt;
                                         uint8_t *s = p;
                                         uint32_t kpos = bg_len >> 1;
                                         int64_t align4 = bg_len >> 2 << 2;
                                         uint32_t hash32, curr_bucket;
                                         int64_t offset = offset_start, pos;
                                         int64_t *pnnn = (pnn) ? pnn : (&n);
                                         uint64_t xsquash;

                                         const uint32_t len_bgs = (bg_len >> 2) + ((bg_len & 0x3) ? 1 : 0);
                                         const char actg4[4] = {'A', 'C', 'T', 'G'};
                                         char actg_bg[bg_len + 1];
                                         char actg_bg_pair[bg_len + 1]; /* 互补碱基*/
                                         char actg_bgs[len_bgs];    /* squash base group */
                                         char actg_bgs_pair[len_bgs];
                                         
                                         for (*pnnn = 0; *pnnn < current; *pnnn = (*pnnn) + 1)
                                         {
                                             /* step 1 : get base group */
                                             for (m = 0; m < align4; m += 4)
                                                 *((uint32_t *)(actg_bg + m)) = actg_stretch_tab[*s++];
                                             switch (bg_len - align4)
                                             {
                                             case 0:
                                                 break;
                                             case 1:
                                                 m = align4;
                                                 *(actg_bg + m) = actg4[(((*s) >> 6) & 0x3)];
                                                 s++;
                                                 break;
                                             case 2:
                                                 m = align4;
                                                 *(actg_bg + m++) = actg4[(((*s) >> 6) & 0x3)];
                                                 *(actg_bg + m) = actg4[(((*s) >> 4) & 0x3)];
                                                 s++;
                                                 break;
                                             case 3:
                                                 m = align4;
                                                 *(actg_bg + m++) = actg4[(((*s) >> 6) & 0x3)];
                                                 *(actg_bg + m++) = actg4[(((*s) >> 4) & 0x3)];
                                                 *(actg_bg + m) = actg4[(((*s) >> 2) & 0x3)];
                                                 s++;
                                                 break;
                                             default:
                                                 break;
                                             }

                                             /* step 2: get base group pair */
                                             actg_pair((uint8_t *)actg_bg_pair, (uint8_t *)actg_bg, bg_len);

                                             /* step 3:  save current base group with direction*/
                                             if (actg_bg[kpos] < actg_bg_pair[kpos])
                                             {
                                                 actg_squash((uint8_t *)actg_bg_pair, bg_len, (uint8_t *)actg_bgs_pair);
                                                 xsquash = *((uint64_t *)(actg_bgs_pair));
                                                 xsquash &= 0xFCFFFFFFFFFFFFFF;
                                                 hash32 = (uint32_t) CityHash64((const char *)(&xsquash), len_bgs);
                                                 pos = offset | 0x80000000;
                                             }
                                             else
                                             {
                                                 actg_squash((uint8_t *)actg_bg, bg_len, (uint8_t *)actg_bgs);
                                                 xsquash = *((uint64_t *)(actg_bgs));
                                                 xsquash &= 0xFCFFFFFFFFFFFFFF;
                                                 hash32 = (uint32_t)CityHash64((const char *)(&xsquash), len_bgs);
                                                 pos = offset;
                                             }
                                             curr_bucket = hash32 & hmask;
                                             (bhash_start + offset)->hash_bucket = curr_bucket;
                                             (bhash_start + offset)->basegroup_pos = pos;
                                             offset++;

                                            /* step 3: 更新对应bucket的计数 */
                                            //  if (hb_cnt[curr_bucket] < 16)
                                             {
                                                 spinlock &sl = slocks[curr_bucket];
                                                 spin_lock(&sl);
                                                 hb_cnt[curr_bucket]++;
                                                 spin_unlock(&sl);
                                             }
                                         }
                                     }));

        p += (current * (basegroup_step >> 2)); /* 这里限定了basegroup_step必须为4的整数，有需要可以修改 */
        bhash += current;
    }

    for (n = 0; n < pcnt; n++)
    {
        if (tpools[n].joinable())
            tpools[n].join();
    }
    free(slocks);
}

/* 制作索引表：计算hash table的size */
void reference::make_step2_calcHT(htable &hash_table)
{
    uint32_t *p;
    int64_t each, current, remain, total;
    uint32_t n, pcnt = this->parallel;
    std::vector<std::thread> tpools;

    /* 第一段存hash buffer的内容和总长度，第二段存hash butcket的长度 */
    hash_table.resize(pcnt);

    total = hash_buckets;
    each = hash_buckets / pcnt;
    remain = (total - (each * pcnt));
    p = hash_bucket_cnt;

    for (n = 0; n < pcnt; n++)
    {
        current = (n + 1 == pcnt) ? (each + remain) : each;

        std::pair<std::pair<uint32_t *, uint32_t>, uint32_t> &hash_buffer = hash_table[n];
        tpools.push_back(std::thread([p, current, &hash_buffer]()
                                     {
                                         int64_t len_hash = 0, len_bucket = 0;
                                         int64_t n = 0, m = 0;
                                         uint32_t current_cnt;

                                         for (n = 0; n < current; n++)
                                         {
                                             current_cnt = *(p + n);
                                             len_bucket++;
                                             /* 如果当前bucket中有多个hash值，那么hash buffer第一个存指针，该指针指向当前bucket的hash值对应的buffer */
                                             len_hash += (current_cnt <= 1) ? 1 : (current_cnt + 1);
                                         }
                                         hash_buffer.first.second = len_hash;
                                         hash_buffer.second = len_bucket;
                                     }));

        p += current;
    }

    for (n = 0; n < pcnt; n++)
    {
        if (tpools[n].joinable())
            tpools[n].join();
    }
}

/* 制作索引表：初始化hash table */
void reference::make_step3_initHT(const htable &hash_table, uint32_t *&hash_bucket_curpos)
{
    uint32_t hash_bufflen = 0;
    uint32_t hash_bucketlen = 0;
    uint32_t *phbucket;
    uint32_t *pbucket_next;
    uint32_t *phbuff, *phbuff_conflict, *phbuff_start;
    int64_t each, current, remain, total;
    uint32_t n, pcnt = this->parallel;
    std::vector<std::thread> tpools;

    for (n = 0; n < hash_table.size(); n++)
    {
        hash_bufflen += hash_table[n].first.second;
        hash_bucketlen += hash_table[n].second;
    }
    safe_alloc(hash_bufflen, uint32_t, hash_table_buffer);
    
    total = hash_buckets;
    each = hash_buckets / pcnt;
    remain = (total - (each * pcnt));
    phbucket = hash_bucket_cnt;
    phbuff = hash_table_buffer;
    phbuff_start = hash_table_buffer;
    phbuff_conflict = hash_table_buffer + hash_bucketlen;
    pbucket_next = hash_bucket_curpos;

    for (n = 0; n < pcnt; n++)
    {
        current = (n + 1 == pcnt) ? (each + remain) : each;
        tpools.push_back(std::thread([phbucket, phbuff, phbuff_conflict, phbuff_start, pbucket_next, &current]()
                                     {
                                         uint32_t n = 0;
                                         uint32_t current_cnt;
                                         uint32_t *p = phbuff;
                                         uint32_t *pc = phbuff_conflict;
                                         uint32_t *pn = pbucket_next;

                                         for (n = 0; n < current; n++)
                                         {
                                             current_cnt = *(phbucket + n);
                                             switch (current_cnt)
                                             {
                                             case 0:
                                                 pn++;
                                                 p++;
                                                 break;
                                             case 1:
                                                 *pn++ = p - phbuff_start; /* 没有hash冲突时，下一个位置直接写bucket */
                                                 p++;
                                                 break;
                                             default:
                                                 *pn++ = pc - phbuff_start;
                                                 *p = pc - phbuff_start; /*  hash冲突时bucket第一个元素存hash冲突buffer相对hash buffer起始位置的偏移 */
                                                 pc += current_cnt;      /*  当前bucket有current_cnt个hash值，故做current_cnt个偏移 */
                                                 p++;
                                                 break;
                                             }
                                         }
                                     }));

        pbucket_next += current;
        phbucket += current;
        phbuff += hash_table[n].second;
        phbuff_conflict += hash_table[n].first.second - hash_table[n].second;
    }

    for (n = 0; n < pcnt; n++)
    {
        if (tpools[n].joinable())
            tpools[n].join();
    }
}

/* 制作索引表：将碱基的位置信息存入hash table*/
void reference::make_step4_buildHT(const bg_hash *hash, uint32_t *&hash_bucket_curpos)
{
    uint32_t *phbuff_start, id;
    std::vector<std::thread> tpools;
    uint32_t *hb_cnt = hash_bucket_cnt;
    int64_t each, current, remain, total;
    int64_t n, pcnt = this->parallel;
    uint32_t *hb_curpos = hash_bucket_curpos;

    total = (reference_squashlen << 2) / basegroup_step;
    each = total / pcnt;
    remain = (total - (each * pcnt));
    const bg_hash *h = (bg_hash *)hash;
    phbuff_start = hash_table_buffer;

    for (n = 0; n < pcnt; n++)
    {
        id = n;
        tpools.push_back(std::thread([h, total, &phbuff_start, &hb_curpos, &hb_cnt, id, pcnt]()
                                     {
                                         timer cost_ms(true);
                                         uint32_t n = 0, bucket, pos;

                                         for (n = 0; n < total; n++)
                                         {
                                             bucket = (h + n)->hash_bucket;
                                             if ((bucket % pcnt) == id)
                                             {
                                                 switch (hb_cnt[bucket])
                                                 {
                                                 case 0:
                                                     break;
                                                 case 1:
                                                     pos = hb_curpos[bucket];
                                                     *(phbuff_start + pos) = (h + n)->basegroup_pos;
                                                     break;
                                                 default:
                                                     pos = hb_curpos[bucket];
                                                     *(phbuff_start + pos) = (h + n)->basegroup_pos;
                                                     hb_curpos[bucket]++; /* 该bucket指向下一个位置 */
                                                     break;
                                                 }
                                             }
                                         }
                                        //  printf("\n\t id %u cost ms: %u\n\n", id, cost_ms.elapsed());
                                     }));
    }

    for (n = 0; n < pcnt; n++)
    {
        if (tpools[n].joinable())
            tpools[n].join();
    }
}

/* 制作索引表：将同一个bucket中的位置排序，保证并发生成的索引表都一致 */
void reference::make_step5_sortHT()
{
    uint32_t *phb_cnt, *phb_cnt_start;
    uint32_t *phb_buff, *phb_buff_start, offset = 0;
    int64_t each, current, remain, total;
    int64_t n, pcnt = this->parallel;
    std::vector<std::thread> tpools;

    total = hash_buckets;
    each = total / pcnt;
    remain = (total - (each * pcnt));
    phb_buff_start = phb_buff = hash_table_buffer;
    phb_cnt_start = phb_cnt = hash_bucket_cnt;

    for (n = 0; n < pcnt; n++)
    {
        current = (n + 1 == pcnt) ? (each + remain) : each;
        tpools.push_back(std::thread([phb_buff, current, phb_buff_start, &phb_cnt_start, offset]()
                                     {
                                         uint32_t current_cnt;
                                         uint32_t n, m, *p;
                                         uint32_t *pstart = phb_buff_start;

                                         for (n = 0; n < current; n++)
                                         {
                                             current_cnt = *(phb_cnt_start + offset + n);
                                             if (current_cnt > 1)
                                             {
                                                 if (current_cnt < 255) {
	                                                 p = pstart + phb_buff[n];
	                                                 std::sort(p, p + current_cnt, [](const uint32_t &p1, const uint32_t &p2) -> bool
	                                                           { return p1 < p2; }); // 修复可能core问题，相等返回true会越界，当元素个数>16(_S_threshold)时选择快速排序，<=16个则选择插入排序(对象少时快排性能不理想)
                                                 }
                                                *(phb_cnt_start + offset + n) = std::min((uint32_t)16, current_cnt);
                                             }
                                         }
                                     }));
        offset += current;
        phb_buff += current;
    }

    for (n = 0; n < pcnt; n++)
    {
        if (tpools[n].joinable())
            tpools[n].join();
    }
}

/* dump hash table to file */
void reference::dump_hash_table()
{
    uint32_t n, m, *p;
    FILE *fp;

    fp = fopen("hash_table", "wb");
    for (n = 0; n < hash_buckets; n++)
    {
        if (hash_bucket_cnt[n] <= 0)
            continue;
        fprintf(fp, "bucket %llu : ", n);
        p = (hash_bucket_cnt[n] > 1) ? (hash_table_buffer + hash_table_buffer[n]) : (hash_table_buffer + n);
            
        for (m = 0; m < hash_bucket_cnt[n]; m++)
            fprintf(fp, "%llu,", *(p + m));
        fprintf(fp, "\n");
    }
    fflush(fp);
    fclose(fp);
}

/* 获取建索引时取的碱基的长度 */
const uint32_t reference::get_bglen() const
{
    return this->basegroup_len;
}

/* 获取建索引时取的碱基的步长 */
const uint32_t reference::get_bgstep() const
{
    return this->basegroup_step;
}

/* 查询位置 */
const uint32_t *reference::query_pos(const uint32_t &hash, uint32_t &len)
{
    uint32_t hbucket;

    hbucket = hash & hash_mask;
    len = hash_bucket_cnt[hbucket];
    return (len == 1) ? (hash_table_buffer + hbucket) : (hash_table_buffer + hash_table_buffer[hbucket]);
}

/* 更新参考基因组匹配信息，以一个squash字节为单位更新，这样可以并发更新 */
void reference::update_matched(uint64_t actg_pos, uint32_t actg_match_len)
{
    uint64_t spos_start = actg_pos >> 2;
    uint64_t spos_end = (actg_pos + actg_match_len - 1) >> 2;
    memset(reference_squash_matched + spos_start, 1, spos_end - spos_start + 1);
}

/* 得到reference squash buffer */
const uint8_t *reference::get_squash() const
{
    return this->reference_squash;
}

/* 得到reference squash buffer的长度 */
const int64_t reference::get_squashlen() const
{
    return this->reference_squashlen;
}

/* 获取fasta文件名 */
const std::string reference::get_fasta_name() const
{
    return this->reference_file;
}

/* 获取fasta文件内容长度 */
const int64_t reference::get_fasta_len() const
{
    return this->fasta_len;
}

/* 获取fasta文件内容md5 */
const std::string reference::get_fasta_md5() const
{
    return this->fasta_md5;
}

/* 获取ni文件路径 */
const std::string reference::get_ni_path() const
{
    return this->path_ni;
}

/* 将没有matched上的reference清零 */
void reference::reference_squash_sanitize(int64_t start_squash_pos, int64_t len)
{
    uint64_t n, e = start_squash_pos + len;
    for (n = start_squash_pos; n < e; n++)
        *(reference_squash + n) = (*(reference_squash_matched + n)) ? (*(reference_squash + n)) : 0;
}

/* 获取指定位置对应长度的actg碱基 */
void reference::get_stretch_actg(uint8_t *out, uint32_t out_len, uint64_t actg_pos)
{
    uint32_t n, len_need, offset = 0;
    uint64_t squash_pos = actg_pos >> 2;
    const char actg4[4] = {'A', 'C', 'T', 'G'};
    uint8_t *p, ch;

    p = reference_squash + squash_pos;
    ch = *p;

    /* 左边不对齐 */
    switch (actg_pos & 0x3) 
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    case 2:
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    case 3:
        *(out + offset++) = actg4[(ch) & 0x3];
        p++;
        break;
    default:
        break;
    }
    if (offset == out_len)
        return;

    /* 对齐部分 */
    len_need = (out_len - offset) >> 2;
    for (n = 0; n < len_need; n++) {
        *((uint32_t *)(out + offset)) = actg_stretch[*p];
        offset += 4;
        p++;
    }
    if (offset == out_len)
        return;

    /* 右边未对齐部分 */
    len_need = out_len - offset;
    ch = *p;
    switch (len_need & 0x3) 
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        break;
    case 2:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        break;
    case 3:
        *(out + offset++) = actg4[(ch >> 6) & 0x3];
        *(out + offset++) = actg4[(ch >> 4) & 0x3];
        *(out + offset++) = actg4[(ch >> 2) & 0x3];
        break;
    default:
        break;
    }
}

/*  根据每个字节末尾的2个bits，转换成actg */
void reference::get_actg_from2bits(const uint8_t *src_2bits, uint32_t src_2blits_len, uint8_t *dst_actg)
{
    uint8_t *s = (uint8_t *)src_2bits, *p;
    uint32_t n, align4 = src_2blits_len >> 2 << 2, offset = 0;
    const uint8_t actg4[4] = {'A', 'C', 'T', 'G'};

    for (n = 0; n < align4; n += 4) {
        p = s + n;
        *((uint32_t *)(dst_actg + offset)) = actg_stretch[((*p) << 6) | ((*(p + 1)) << 4) | ((*(p + 2)) << 2) | (*(p + 3))];
        offset += 4;
    }

    if (n == src_2blits_len)
        return;

    for (; n < src_2blits_len; n++)
        *(dst_actg + offset++) = actg4[*(src_2bits + n)];
}

/* 获取指定位置对应长度的squash碱基，即2个bits 放到一个字符的末尾*/
void reference::get_stretch_2bits1char(uint8_t *out, uint32_t out_len, uint64_t actg_pos)
{
    uint32_t n, len_need, offset = 0;
    uint64_t squash_pos = actg_pos >> 2;
    uint8_t *p, ch;

    p = reference_squash + squash_pos;
    ch = *p;

    /* 左边不对齐 */
    switch (actg_pos & 0x3) 
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = (ch >> 4) & 0x3;
        *(out + offset++) = (ch >> 2) & 0x3;
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    case 2:
        *(out + offset++) = (ch >> 2) & 0x3;
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    case 3:
        *(out + offset++) = (ch) & 0x3;
        p++;
        break;
    default:
        break;
    }
    if (offset == out_len)
        return;

    /* 对齐部分 */
    len_need = (out_len - offset) >> 2;
    for (n = 0; n < len_need; n++) {
        *((uint32_t *)(out + offset)) = actg_stretch_2bits[*p];
        offset += 4;
        p++;
    }
    if (offset == out_len)
        return;

    /* 右边未对齐部分 */
    len_need = out_len - offset;
    ch = *p;
    switch (len_need & 0x3) 
    {
    case 0:
        break;
    case 1:
        *(out + offset++) = (ch >> 6) & 0x3;
        break;
    case 2:
        *(out + offset++) = (ch >> 6) & 0x3;
        *(out + offset++) = (ch >> 4) & 0x3;
        break;
    case 3:
        *(out + offset++) = (ch >> 6) & 0x3;
        *(out + offset++) = (ch >> 4) & 0x3;
        *(out + offset++) = (ch >> 2) & 0x3;
        break;
    default:
        break;
    }
}
