#include <iostream>
#include <getopt.h>
#include <filesystem>

#include "pbgz_file.h"
#include "log/logger.h"
#include "pbgz_engine.h"
#include "pbgz_types.h"
#include "pbgz_manager.h"
#include "utils/path_util.h"
#include "config_manager.h"


typedef struct
{
    char argShort;
    std::string argLong;
    int hasFlag;
    std::string argDescribe;
} Arg;

typedef std::vector<Arg> PbgzArg;

static PbgzArg pbgzArgs = 
{
    {'d', "decompress", no_argument, "<file.pbgz> specify file to decompress"},
    {'z', "gz", no_argument, "decompress to gz format"},
    {'o', "outfile", required_argument, "<outfile> specify output filename for compress/decompress, - means output to pipe"},
    {'O', "outdir", required_argument, "<outdir> specify output directory for compress/decompress"},
    {'f', "force", no_argument, "force overwrite of output file. use for compress/decompress"},
    {'r', "reference", required_argument, "<FASTA> specify reference gene file (it's not Mandatory). use for compress/decompress"},
    {'n', "refunpack", no_argument, "unpack reference to pbgz file, so reference gene is needed when decompress"},
    {'t', "threads", required_argument, "<number> specify number of threads, default is CPUS. use for compress/decompress"},
    {'l', "level", required_argument, "<1-9> specify compress level, default is 5, 1 is fast, 9 is best"},
    {'e', "remove", no_argument, "if compress succeed, will remove origin file, else this option is invalid. use for compress"},
    {'h', "help",  no_argument, "show help"},
    {'v', "version", no_argument, "show version"},
    {'g', "loglevel", required_argument, "<1-6> sepcify log level, default 6. debug(1), info(2), warning(3), error(4), fatal(5), off(6)"},
    {'G', "logfile", required_argument, "sepcify log file"},
};


void printUsage() {
#define ALIGN " "
#define CMD_ARG(f, s, d, p) \
    fprintf(f, "%2s-%c,%s--%-12s%s\n", ALIGN, s, ALIGN, d, p)

    FILE *fp = stdout;
    fprintf(fp, "pbgz: %s\n", PbgzManager::getInstance().getVersion().c_str());
    fprintf(fp, "Usage: pbgz [OPTION] [FILE]\n\n");
    fprintf(fp, "Mandatory arguments to long options are mandatory for short options too.\n\n");

    for (auto &a : pbgzArgs)
        CMD_ARG(fp, a.argShort, a.argLong.c_str(), a.argDescribe.c_str());

    fprintf(fp, "\nTo compress\n  pbgz human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
    fprintf(fp, "To decompress\n  pbgz -d human.fq.gz.pbgz\n\n");
}

void printVersion() {
    fprintf(stdout, "pbgz: %s\n", PbgzManager::getInstance().getVersion().c_str());
}

int reconstrucParamters(PbgzParameter& para) {
    if (para.inputFile.empty()) {
        para.inputFile = "/dev/stdin";
    }

    if ("-" == para.outputFile ) {
        para.outputFile = "/dev/stdout";
    }

    para.inputFile = PathUtil::getAbsPath(para.inputFile);
    if (!para.outputDir.empty() && !para.outputFile.empty()) {  // 输出目录和文件同时指定
        fprintf(stderr, "Output directory and filename cannot be specified simultaneously.\n");
        return -1;
    } else if(!para.outputDir.empty()) {   // 只指定输出目录场景
        if (para.inputFile == "/dev/stdin") {
            fprintf(stderr, "Specify an output filename, or do not specify both the filename and the directory.\n");
            return -1;
        } else {
            // 输出文件名为输出文件名加后缀(.pbgz)
            std::string inFileName = PathUtil::getFileName(para.inputFile);
            if (inFileName.empty()) {
                fprintf(stderr, "%s not exits or no permission to access.\n", para.inputFile.c_str());
                return -1;
            }

            std::string outPath = PathUtil::getFilePath(para.outputDir);
            if (outPath.empty()) {
                // 目录不存在就尝试创建此目录
                outPath = PathUtil::createDir(para.outputDir);
                if (outPath.empty()) {
                    fprintf(stderr, "%s not exits or no permission to access.\n", para.outputDir.c_str());
                    return -1;
                }
            }
            para.outputFile = outPath + inFileName + ".pbgz";
        }
    } else if (!para.outputFile.empty()) {     // 指定文件输出场景
        // 指定了输出文件，转成绝对路径
        para.outputFile = PathUtil::getAbsPath(para.outputFile);
    } else {
        // 输出目录和文件均未指定
        if (para.inputFile == "/dev/stdin") {
            para.outputFile = "/dev/stdout";
        } else {
            if (!PathUtil::fileExists(para.inputFile)) {
                fprintf(stderr, "%s not exits or no permission to access.\n", para.inputFile.c_str());
                return -1;
            }
            // 压缩场景，在出入文件同目录下输出,输出文件名加后缀(.pbgz)
            if (!para.isDecompress) {
                para.outputFile = PathUtil::getFilePath(para.inputFile) + PathUtil::getFileName(para.inputFile) + ".pbgz";
            } else {  
                // 解压场景，去掉文件的pbgz后缀
                std::string name = PathUtil::getFileName(para.inputFile);
                if (name.length() <= 5) {
                    fprintf(stderr, "input file %s invalid.", name.c_str());
                    return -1;
                }
                if (!PathUtil::suffixCheck(name, ".pbgz")) {
                    fprintf(stdout, "Decompress file must be pbgz format.\n");
                    return -1;
                }
                name.resize(name.length() - 5);
                para.outputFile = PathUtil::getFilePath(para.inputFile) + name;
            }
        }
    }

    // 输入了参考基因，转换成绝对路径
    if (!para.referenceGenic.empty()) {
        std::string abspath = PathUtil::getAbsPath(para.referenceGenic);
        if (abspath.empty()) {
            fprintf(stderr, "%s not exits or no permission to access.\n", para.inputFile.c_str());
            return -1;
        }
        para.referenceGenic = abspath;
    }

    // 如果文件已经存在，并且配置了强制覆盖，则删除输出文件
    if (PathUtil::fileExists(para.outputFile) && para.isOverwriteOutFile) {
        PathUtil::removeFile(para.outputFile);
    }

    ConfigManager::getInstance().init(para);
    return 0;
}

int checkParameters(PbgzParameter& para) {
    if (para.isDecompress) {   // 解压场景，输入文件名为pbgz结尾
        if ("/dev/stdin" != para.inputFile) {
            if (PathUtil::isDir(para.inputFile)) {
                fprintf(stdout, "Decompress file cannot be directory.\n");
                return -1;
            }
            
            if (!PathUtil::suffixCheck(para.inputFile, "pbgz")) {
                fprintf(stdout, "Decompress file must be pbgz format.\n");
                return -1;
            }
        }
    } else {   // 压缩场景，输出文件名以pbgz结尾
        if ("/dev/stdout" != para.outputFile){
            if (PathUtil::isDir(para.outputFile)) {
                fprintf(stdout, "Output file cannot be directory.\n");
                return -1;
            }

            if (!PathUtil::suffixCheck(para.outputFile, "pbgz")) {
                fprintf(stdout, "Output file must be pbgz format.\n");
                return -1;
            }
        }
    }

    if (para.outputFile != "/dev/stdout") {
        // 输出文件已经存在时,提醒使用-f强制覆盖输出
        if (PathUtil::fileExists(para.outputFile) && !para.isOverwriteOutFile) {
            fprintf(stdout, "%s exits, use -f to force overwrite.\n", para.outputFile.c_str());
            return -1;
        }

        // 指定文件并且覆盖时，检查是否有写的权限
        if (PathUtil::fileExists(para.outputFile)) {
            if (!PathUtil::fileWriteable(para.outputFile)) {
                fprintf(stdout,"no permisstion to write %s.\n", para.outputFile.c_str());
                return -1;
            }
        } else {  // 文件不存在, 则需要检查文件所在的目录是否有权限
            // 到了这里输出文件已经转换成了绝对路径，必然是带文件名的
            char separator = std::filesystem::path::preferred_separator;
            int pos = para.outputFile.find_last_of(separator);
            if (pos == -1) {
                fprintf(stderr, "invalid file path, %s.\n", para.outputFile.c_str());
                return -1;
            }

            std::string path = para.outputFile.substr(0, pos + 1);
            if (!PathUtil::fileExists(path)) {
                fprintf(stderr, "%s not exists.\n", para.outputFile.c_str());
                return -1;
            }

            if (!PathUtil::fileWriteable(path)) {
                fprintf(stderr, "no permssion to access %s.\n", para.outputFile.c_str());
                return -1;
            }
        }
    }

    if (para.inputFile != "/dev/stdin") {
        if (!PathUtil::fileExists(para.inputFile)) {
            fprintf(stderr, "%s not exits.\n", para.inputFile.c_str());
            return -1;
        }

        if (!PathUtil::fileReadble(para.inputFile)) {
            fprintf(stderr, "no permission to access %s.\n", para.inputFile.c_str());
            return -1;
        }
    }

    if (para.compressLevel < 1 || para.compressLevel > 9) {
        fprintf(stderr, "Compress level is invalid.\n");
        return -1;
    }

    if (para.logLevel < 1 || para.logLevel > 6) {
         fprintf(stderr, "Log level is invalid.\n");
        return -1;
    }

    return 0;
}

int main(int argc, char** argv) {
    if (argc == 1) {
        printUsage();
        return 0;
    }

    std::string argOption = "";
    option longopts[pbgzArgs.size() + 1];
    for (uint32_t i = 0; i < pbgzArgs.size(); ++i) {
        longopts[i] = {pbgzArgs[i].argLong.c_str(), pbgzArgs[i].hasFlag, nullptr, pbgzArgs[i].argShort};
        argOption += pbgzArgs[i].argShort;
        switch (longopts[i].has_arg) 
        {
        case optional_argument:
            argOption += "::";
            break;
        case required_argument:
            argOption += ":";
            break;
        case no_argument:
        default:
            break;
        }
    }
    longopts[pbgzArgs.size()] =  {nullptr, 0, nullptr, 0};

    int opt = 0;
    PbgzParameter parameter;
    while((opt = getopt_long(argc, argv, argOption.c_str(), longopts, NULL)) != -1) {
        switch (opt) {
        case 'h': {
            parameter.showHelp = true;
            break;
        }
        case 'v' : {
            parameter.showVersion = true;
            break;
        }
        case 'd' : {
            parameter.isDecompress = true;
            break;
        }
        case 'z':{
            parameter.isDecToGZ = true;
            break;
        }
        case 'o': {
            parameter.outputFile = optarg;
            break;
        }
        case 'O':{
            parameter.outputDir = optarg;
            break;
        }
        case 'f':{
            parameter.isOverwriteOutFile = true;
            break;
        }
        case 'r':{
            parameter.referenceGenic = optarg;
            break;
        }
        case 'n': {
            parameter.isUnpackRef = true;
            break;
        }
        case 't': {
            parameter.threadNum = atoi(optarg);
            break;
        }
        case 'l': {
            parameter.compressLevel = atoi(optarg);
            break;
        }
        case 'e': {
            parameter.isRemoveOriginFile = true;
            break;
        }
        case 'g': {
            parameter.logLevel = atoi(optarg);
            break;
        }
         case 'G': {
            parameter.logFile = optarg;
            break;
        }
        case '?':
        default:
            fprintf(stdout, "Use -h or --help to show usage");
            return 0;
        }
    }

    if (optind < argc) {
        if (argc - optind > 1) {
            fprintf(stdout, "Only one file is allowed");
            return 0;
        }
        parameter.inputFile = argv[optind];
    }

    if (parameter.showHelp || parameter.showVersion) {
        if (argc > 2) {
            fprintf(stdout, "Invalid parameter: exits multi option");
            return -1;
        }
    }

    if (parameter.showHelp) {
        printUsage();
        return 0;
    }

    if (parameter.showVersion) {
        printVersion();
        return 0;
    }

    if (reconstrucParamters(parameter) != 0) {
        LOG_ERROR("Invalid parameter");
        return -1;
    }

    if (checkParameters(parameter) != 0) {
        LOG_ERROR("Invalid parameter");
        return -1;
    }

    PbgzEngine engin(parameter);
    if(0 != engin.init()) {
        LOG_ERROR("Pbgz Engin init failed");
        return -1;
    }

    if (0 != engin.start()) {
        LOG_ERROR("Pbgz Engin start failed");
        return -1;
    }
    
    return 0;
}

