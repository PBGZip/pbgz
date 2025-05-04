#include "manager.h"

manager::manager()
{
}

manager::~manager()
{
}

void manager::set_zipinfo(const info &val)
{
    zipinfo = val;
}

const info &manager::get_zipinfo()
 {
    return this->zipinfo;
 }

void manager::exit(int code, const std::string &info, int err)
{
    if (code != ERR_NO)
    {
        /* 非正常的输出文件要清除 */
        for (auto &curfile : outfiles)
        {
            const std::string &file = curfile.first;
            if (!curfile.second && file_exists(file))
                unlink(file.c_str());
        }
    }

    if (info.size())
        fprintf(err ? stderr : stdout, "pbgz exit : %s\n", info.c_str());

    _Exit(code);
}

/* 将n转换为最接近2的幂次方的整数 */
int32_t powerof2_proximal(int32_t i)
{
    i |= (i >> 1);
    i |= (i >> 2);
    i |= (i >> 4);
    i |= (i >> 8);
    i |= (i >> 16);
    return i - (i >> 1);
}