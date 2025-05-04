#ifndef _FILE_OP_H_
#define _FILE_OP_H_

#include <sys/stat.h>
#if defined(__aarch64__)
#include <sys/uio.h>
#else
    // Conditionally include sys/io.h only on Linux platforms
    #if defined(__linux__) || defined(__linux) || defined(linux)
    #include <sys/io.h>
    #endif
#endif
#include <unistd.h>
#include <iostream>
#include <limits.h>
#include <memory.h>

/* 文件操作 */

/*  文件或文件夹是否存在*/
bool file_exists(const std::string &file);

/*  文件或文件夹是否有读权限*/
bool file_readable(const std::string &file);

/*  文件或文件夹是否有写权限*/
bool file_writeable(const std::string &file);

/* 路径是否是目录 */
bool is_directory(const std::string &path);

/* 路径是否是文件 */
bool is_file(const std::string &path);

/* 返回文件的绝对路径对应的目录*/ 
void file_abspath_dir(const std::string &file, std::string &get);

/* 返回文件的绝对路径对应的文件名*/ 
void file_abspath_filename(const std::string &file, std::string &get);

/* 判断文件是否以后缀结尾 */
bool file_suffix_check(const std::string &file, const std::string &suffix);

/* 文件追加后缀，如果已经带了需要追加的后缀则不再追加 */ 
void file_suffix_add(const std::string &file, const std::string &suffix, std::string &get);

/* 文件删除后缀 ，如果文件后缀不是需要删除的后缀则返回原文件，如果是则一直递归删除 */
void file_suffix_del(const std::string &file, const std::string &suffix, std::string &get);

/* 创建目录 */
void create_dir(const char *dir2create);

/* 获取文件锁 */
int32_t acquire_lock(std::string filelock);

/* 释放文件锁 */
void release_lock(int lockFd, std::string filelock);

/* 获取文件最后的修改时间 */
int64_t file_last_mtime(const std::string &path);

/* 获取文件长度 */
int64_t file_length(const std::string &file_name);

/* 读文件指定长度内容计算md5 */
void file_calc_md5(const std::string &file_name, int32_t len, std::string &md5);

#endif
