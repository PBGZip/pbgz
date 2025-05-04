#include "platform_compat.h"

#include "file_op.h"
#include "manager.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include "cfgpath/cfgpath.h"
#include "block.h"

/*  文件或文件夹是否存在*/
bool file_exists(const std::string &file)
{
    return (access(file.c_str(), F_OK) == 0);
}

/*  文件或文件夹是否有读权限*/
bool file_readable(const std::string &file)
{
    return (access(file.c_str(), R_OK) == 0);
}

/*  文件或文件夹是否有写权限*/
bool file_writeable(const std::string &file)
{
    return (access(file.c_str(), W_OK) == 0);
}

/* 路径是否是目录 */
bool is_directory(const std::string &path)
{
    struct stat s_buf;
    int result;
    if (0 != stat(path.c_str(), &s_buf))
        return false;
    return S_ISDIR(s_buf.st_mode);
}

/* 路径是否是文件 */
bool is_file(const std::string &path)
{
    struct stat s_buf;
    int result;
    if (0 != stat(path.c_str(), &s_buf))
        return false;
    return S_ISREG(s_buf.st_mode);
}

/* 返回文件的绝对路径对应的目录*/
void file_abspath_dir(const std::string &file, std::string &get)
{
    char abs_path[PATH_MAX];
    memset(abs_path, 0, PATH_MAX);
    realpath(file.c_str(), abs_path);
    std::string file_path(abs_path, strlen(abs_path));
    if (is_directory(abs_path))
    {
        get = abs_path;
        get += "/";
    }
    else
        get = file_path.substr(0, file_path.find_last_of("\\/") + 1);
}

/* 返回文件的绝对路径对应的文件名*/
void file_abspath_filename(const std::string &file, std::string &get)
{
    char abs_path[PATH_MAX];
    memset(abs_path, 0, PATH_MAX);
    realpath(file.c_str(), abs_path);
    std::string file_path(abs_path, strlen(abs_path));
    int pos = file_path.find_last_of("\\/");
    get.assign(file_path, pos + 1, file_path.size() - pos - 1);
}

/* 判断文件是否以后缀结尾 */
bool file_suffix_check(const std::string &file, const std::string &suffix)
{
    if (suffix.empty())
        return true;
    if (suffix.length() > file.length())
        return false;
    return (!strncmp(file.c_str() + file.length() - suffix.length(), suffix.c_str(), suffix.length()));
}

/* 文件追加后缀，如果已经带了需要追加的后缀则不再追加 */
void file_suffix_add(const std::string &file, const std::string &suffix, std::string &get)
{
    get = file;
    if (suffix.empty() || file_suffix_check(file, suffix))
        return;
    get += suffix;
}

/* 文件删除后缀 ，如果文件后缀不是需要删除的后缀则返回原文件，如果是则一直递归删除 */
void file_suffix_del(const std::string &file, const std::string &suffix, std::string &get)
{
    get = file;
    if (suffix.empty() || !file_suffix_check(file, suffix))
        return;
    while (file_suffix_check(get, suffix))
        get.assign(get.c_str(), get.length() - suffix.length());
}

/* 创建目录 */
void create_dir(const char *dir2create)
{
    char dname[MAX_PATH];
    strcpy(dname, dir2create);
    int32_t i, len = strlen(dname);
    for (i = 1; i < len; i++)
    {
        if (dname[i] == '/')
        {
            dname[i] = 0;
            if (!file_exists(dname))
                check_exit(mkdir(dname, 0755) != -1, ERR_FILE_WRITE, "create dir failed: %s", dir2create);
            dname[i] = '/';
        }
    }
}

/* 获取文件锁 */
int32_t acquire_lock(std::string filelock)
{
    int32_t lockFd;

    if ((lockFd = open(filelock.c_str(), O_CREAT | O_RDWR, 0666)) < 0)
        return -1;

    if (flock(lockFd, LOCK_EX | LOCK_NB) < 0)
    {
        close(lockFd);
        return -1;
    }

    return lockFd;
}

/* 释放文件锁 */
void release_lock(int lockFd, std::string filelock)
{
    flock(lockFd, LOCK_UN);
    close(lockFd);
    unlink(filelock.c_str());
}

/* Get file's last modification time */
int64_t file_last_mtime(const std::string &path)
{
    struct stat s_buf;
    if (0 != stat(path.c_str(), &s_buf))
        return INT64_MAX;

#if defined(__APPLE__) || defined(__MACH__)
    return s_buf.st_mtimespec.tv_sec;
#else
    return s_buf.st_mtim.tv_sec;
#endif
}

/* 获取文件长度 */
int64_t file_length(const std::string &file_name)
{
    FILE *fp;
    int64_t file_size;

    fp = fopen64(file_name.c_str(), "rb");
    check_exit(fp, ERR_FILE_READ, "open failed: %s", file_name.c_str());
    check_exit(fseeko64(fp, 0, SEEK_END) == 0, ERR_FILE_READ, "fseek failed: %s", file_name.c_str());
    file_size = ftello64(fp);
    fclose(fp);

    return file_size;
}

/* 读文件指定长度内容计算md5 */
void file_calc_md5(const std::string &file_name, int32_t len, std::string &md5)
{
    FILE *fp;
    uint8_t *buff;
    
    safe_alloc(len, uint8_t, buff);
    fp = fopen64(file_name.c_str(), "rb");
    check_exit(fp, ERR_FILE_READ, "open failed: %s", file_name.c_str());

    fread(buff, len, 1, fp);
    calc_md5(md5, buff, len);

    fclose(fp);
    free(buff);
}