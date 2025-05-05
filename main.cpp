#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include "manager.h"
#include "vendor.h"
#include "engine.h"

typedef struct
{
    char arg_short;
    std::string arg_long;
    std::string arg_describe;
} arg;
typedef std::vector<arg> ziparg;

static void get_arg(ziparg &a)
{
    a.push_back({'d', "decompress", "<file.pbgz> specify file to decompress"});
    a.push_back({'z', "gz", "decompress to gz format"});
    a.push_back({'o', "outfile", "<outfile> specify output filename for compress/decompress"});
    a.push_back({'O', "outdir", "<outdir> specify output directory for compress/decompress"});
    a.push_back({'f', "force", "force overwrite of output file. use for compress/decompress"});
    a.push_back({'r', "reference", "<FASTA> specify reference gene file (it's not Mandatory). use for compress/decompress"});
    a.push_back({'n', "refunpack", "unpack reference to pbgz file, so reference gene is needed when decompress"});
    a.push_back({'t', "threads", "<number> specify number of threads, default is CPUS. use for compress/decompress"});
    a.push_back({'l', "level", "<1-3> specify compress level, default is 2, 1 is fast, 3 is best"});
    a.push_back({'e', "remove", "if compress succeed, will remove origin file, else this option is invalid. use for compress"});
    a.push_back({'h', "help", "give help info"});
}

static void usage(int err, const ziparg &arg)
{
#define ALIGN " "
#define CMD_ARG(f, s, d, p) \
    fprintf(f, "%2s-%c,%s--%-12s%s\n", ALIGN, s, ALIGN, d, p)

    FILE *fp = err ? stderr : stdout;
    fprintf(fp, "pbgz: %s\n", ver::instance().ver().c_str());
    fprintf(fp, "Usage: pbgz [OPTION] [FILE]\n\n");
    fprintf(fp, "Mandatory arguments to long options are mandatory for short options too.\n\n");

    for (auto &a : arg)
        CMD_ARG(fp, a.arg_short, a.arg_long.c_str(), a.arg_describe.c_str());

    fprintf(fp, "\nTo compress\n  pbgz human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
    fprintf(fp, "To decompress\n  pbgz -d human.fq.gz.pbgz\n\n");

    manage::instance().exit(ERR_BAD_ARGS, "", err);
}

static bool cmd_parse(int arg_n, char **arg_s, const ziparg &arg, info &zip_info)
{
    std::string reference, cmd;
    std::string infile, outfile, outdir;
    std::string outpath;
    int d = -1, z = -1, f = -1, n = -1, t = -1, l = -1, h = -1, e = -1;

    for (int i = 0; i < arg_n; i++)
    {
        cmd += " ";
        cmd += arg_s[i];
    }
    zip_info.set_cmd(cmd);

    if (arg_n == 1)
    {
        usage(1, arg);
        return true;
    }

#define check_arg(c, info)                                  \
    if (!(c))                                               \
    {                                                       \
        fprintf(stderr, "error: %s: %s\n", info, arg_s[i]); \
        return false;                                       \
    }

    for (int i = 1; i < arg_n; i++)
    {
        char *p = arg_s[i];
        int len = strlen(p);

        if (*p == '-')
        { // parse as arg
            char opt_short = '\0';
            std::string opt_long;

            check_arg(len >= 2, "invalid option");
            if (len == 2)
                opt_short = *(++p);
            else
            {
                check_arg(*(++p) == '-', "invalid option");
                opt_long.assign(++p, len - 2);
            }

            if (opt_short == 'z' || opt_long == "gz")
            {
                check_arg(z == -1, "exists multi option");
                z = 1;
            }
            else if (opt_short == 'f' || opt_long == "force")
            {
                check_arg(f == -1, "exists multi option");
                f = 1;
            }
            else if (opt_short == 'n' || opt_long == "refunpack")
            {
                check_arg(n == -1, "exists multi option");
                n = 1;
            }
            else if (opt_short == 'e' || opt_long == "remove")
            {
                check_arg(e == -1, "exists multi option");
                e = 1;
            }
            else if (opt_short == 'h' || opt_long == "help")
            {
                check_arg(h == -1, "exists multi option");
                h = 1;
            }
            else if (opt_short == 'd' || opt_long == "decompress")
            { /* 开始解析需要带参数的命令行 */
                check_arg(d == -1, "exists multi option");
                d = 1;
                check_arg(infile.empty(), "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                infile = std::string(arg_s[++i]);
                check(file_exists(infile), false, "file isn't exists: %s", infile.c_str());
                check(file_readable(infile), false, "file does not have read permission: %s", infile.c_str());
            }
            else if (opt_short == 'r' || opt_long == "reference")
            {
                check_arg(reference.empty(), "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                std::string refgene = std::string(arg_s[++i]);
                check(file_exists(refgene), false, "file isn't exists: %s", refgene.c_str());
                check(file_readable(refgene), false, "file does not have read permission: %s", refgene.c_str());
                reference.assign(refgene);
            }
            else if (opt_short == 't' || opt_long == "threads")
            {
                check_arg(t == -1, "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                t = atoi(arg_s[++i]);
                t = std::min(static_cast<long int>(t), sysconf(_SC_NPROCESSORS_ONLN));
            }
            else if (opt_short == 'l' || opt_long == "level")
            {
                check_arg(l == -1, "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                l = atoi(arg_s[++i]);
                check(l >= 1 && l <= 3, false, "level %d is invalid, should in [%d, %d]", l, 1, 3);
            }
            else if (opt_short == 'o' || opt_long == "outfile")
            {
                check_arg(outfile.empty(), "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                outfile = std::string(arg_s[++i]);
            }
            else if (opt_short == 'O' || opt_long == "outdir")
            {
                check_arg(outdir.empty(), "exists multi option");
                check_arg(i + 1 < arg_n, "missing parameter");
                outdir = std::string(arg_s[++i]);
            }
            else
            {
                check_arg(false, "invalid option");
            }
        }
        else
        { // parse as input file, compress
            check_arg(infile.empty(), "exists multi option");
            infile = std::string(arg_s[i]);
            check(file_exists(infile), false, "file isn't exists: %s", infile.c_str());
            check(file_readable(infile), false, "file does not have read permission: %s", infile.c_str());
        }
    }

    check(!infile.empty(), false, "no input file is specfied");
    if (d == 1)
    {
        check(file_suffix_check(infile, ZIP_SUFFIX), false, "decompress file %s don't suffix with %s", infile.c_str(), ZIP_SUFFIX);
        check(l == -1, false, "it's a compress option: 'l'");
        check(n == -1, false, "it's a compress option: 'n'");
    }
    else
    {
        check(z == -1, false, "it's a decompress option: 'z'");
    }

    /* 检查输出文件的权限，及是否存在，是否有强制覆盖选项 */

    /*  得到输出文件的文件路径*/
    if (!outfile.empty() && !outdir.empty())
    { /* 如果同时指定了输出文件和输出路径 ，则报错 */
        check(false, false, "conflicting option: '-o' and '-O'");
    }
    else if (outfile.empty() && outdir.empty())
    { /*  没有指定输出文件名或者输出路径时，默认输出到原文件同路径下*/
        std::string outfile_path, outfile_name;
        file_abspath_filename(infile, outfile_name);
        file_abspath_dir(infile, outfile_path);

        check(file_exists(outfile_path), false, "directory doesn't exist: %s", outfile_path.c_str());
        check(file_writeable(outfile_path), false, "directory does not have write permission: %s", outfile_path.c_str());

        outpath = outfile_path;
        outpath += outfile_name;
        if (d != 1) /* 加压时追加ZIP_SUFFIX */
            file_suffix_add(outpath, ZIP_SUFFIX, outpath);
        else
        {
            file_suffix_del(outpath, ZIP_SUFFIX, outpath);
            if (z == 1)
                file_suffix_add(outpath, GZIP_SUFFIX, outpath);
        }
    }
    else if (!outdir.empty())
    { /* 指定了输出路径*/
        check(file_exists(outdir), false, "directory isn't exists: %s", outdir.c_str());
        check(file_writeable(outdir), false, "directory does not have write permission: %s", outdir.c_str());

        std::string outfile_path, outfile_name;
        file_abspath_filename(infile, outfile_name);
        file_abspath_dir(outdir, outfile_path);
        outpath = outfile_path;
        outpath += outfile_name;

        if (d != 1) /* 加压时追加ZIP_SUFFIX */
            file_suffix_add(outpath, ZIP_SUFFIX, outpath);
        else
        {
            file_suffix_del(outpath, ZIP_SUFFIX, outpath);
            if (z == 1)
                file_suffix_add(outpath, GZIP_SUFFIX, outpath);
        }
    }
    else if (!outfile.empty())
    { /* 指定了输出文件的文件路径*/
        /*  检查输出文件名是否合法*/
        if (d != 1) /* 加压时输出文件名必须以ZIP_SUFFIX结尾 */
        {
            check(file_suffix_check(outfile, ZIP_SUFFIX), false, "compressed output file %s don't suffix with %s", outfile.c_str(), ZIP_SUFFIX);
        }
        else
        {
            check(!file_suffix_check(outfile, ZIP_SUFFIX), false, "decompressed output file %s can not suffix with %s", outfile.c_str(), ZIP_SUFFIX);
            if (z == 1)
            {
                check(file_suffix_check(outfile, GZIP_SUFFIX), false, "decompressed output file %s to gz must suffix with %s", outfile.c_str(), GZIP_SUFFIX);
            }
            else
            {
                check(!file_suffix_check(outfile, GZIP_SUFFIX), false, "decompressed output file %s suffix with %s, but no option '-z'", outfile.c_str(), GZIP_SUFFIX);
            }
        }

        /*  检查输出文件路径是否有写权限*/
        std::string outfile_path, outfile_name;
        file_abspath_dir(outfile, outfile_path);
        check(file_exists(outfile_path), false, "directory isn't exists: %s", outfile_path.c_str());
        check(file_writeable(outfile_path), false, "directory does not have write permission: %s", outfile_path.c_str());

        file_abspath_filename(outfile, outfile_name);
        outpath = outfile_path;
        outpath += outfile_name;
    }

    /* 得到输出文件后，检查文件是否存在，是否有写权限 */
    check(!outpath.empty(), false, "cannot parse output file path");
    if (file_exists(outpath))
    {
        check(f == 1, false, "%s is already exists, use '-f' to overwrite", outpath.c_str());
        check(file_writeable(outpath), false, "file does not have write permission: %s", outpath.c_str());
    }

    zip_info.set_mode((d == 1) ? UNZIP : ZIP);
    zip_info.set_decompres2gz(z == 1);
    zip_info.set_infile(infile);
    zip_info.set_outfile(outpath);
    zip_info.set_force(f == 1);
    zip_info.set_refgene(reference);
    zip_info.set_nopackref(n == 1);
    zip_info.set_threads((t == -1) ? (sysconf(_SC_NPROCESSORS_ONLN)) : t);
    zip_info.set_clevel((l == -1) ? 0 : l);
    zip_info.set_remove(e == 1);
    zip_info.dump();
    return true;
}

int main(int argc, char **argv)
{
    ziparg arg;
    info zinfo;

    // 参数获取
    get_arg(arg);

    // 命令行解析
    if (!cmd_parse(argc, argv, arg, zinfo))
        manage::instance().exit(ERR_BAD_ARGS);
    else
        manage::instance().set_zipinfo(zinfo);

    // 启动引擎
    engine().start();

    return 0;
}