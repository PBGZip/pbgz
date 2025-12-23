#include <iostream>
#include <getopt.h>
#include <filesystem>

#include "pbgz_file.h"
#include "log/logger.h"
#include "pbgz_engine.h"
#include "pbgz_types.h"
#include "pbgz_manager.h"
#include "utils/path_util.h"
#include "utils/memory_util.h"
#include "config_manager.h"
#include "decompress_engine.h"
#include "compress_engine.h"
#include "index_engine.h"


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
    {'g', "loglevel", required_argument, "<1-6> sepcify log level, default 6. debug(1), info(2), warning(3), error(4), fatal(5), off(6)"},
    {'G', "logfile", required_argument, "sepcify log file"},
    {'p', "position", required_argument, "sepecify the reference gene posision"},
    {'i', "index", no_argument, "make index"},
};

class CommandProc {
protected:
    PbgzParameter parameter;
    PbgzEngine* engine;

public:
    CommandProc(PbgzParameter& para) {
        parameter = para;
        engine = nullptr;
    }

    virtual ~CommandProc() {
        MemoryUtil::safeDeleteClass(engine);
    }

    int32_t reconstruct() {
        if (0 != reconstructPorc()) {
            return -1;
        }

        return afterConstructProc();
    }

    virtual int32_t reconstructPorc() {
        if (parameter.inputFile.empty()) {
            parameter.inputFile = "/dev/stdin";
        }
        if ("-" == parameter.outputFile ) {
            parameter.outputFile = "/dev/stdout";
        }
        parameter.inputFile = PathUtil::getAbsPath(parameter.inputFile);

        if (!parameter.outputDir.empty() && !parameter.outputFile.empty()) {  // Output directory and file specified simultaneously
            fprintf(stderr, "Output directory and filename cannot be specified simultaneously.\n");
            return -1;
        } else if(!parameter.outputDir.empty()) {   // Only output directory specified
            if (parameter.inputFile == "/dev/stdin") {
                fprintf(stderr, "Specify an output filename, or do not specify both the filename and the directory.\n");
                return -1;
            } else {
                // Output filename is input filename with suffix (.pbgz)
                std::string inFileName = PathUtil::getFileName(parameter.inputFile);
                if (inFileName.empty()) {
                    fprintf(stderr, "%s not exits or no permission to access.\n", parameter.inputFile.c_str());
                    return -1;
                }

                std::string outPath = PathUtil::getFilePath(parameter.outputDir);
                if (outPath.empty()) {
                    // Try to create directory if it doesn't exist
                    outPath = PathUtil::createDir(parameter.outputDir);
                    if (outPath.empty()) {
                        fprintf(stderr, "%s not exits or no permission to access.\n", parameter.outputDir.c_str());
                        return -1;
                    }
                }
                parameter.outputFile = outPath + inFileName + ".pbgz";
            }
        } else if (!parameter.outputFile.empty()) {     // File output specified
            // Convert output file to absolute path
            parameter.outputFile = PathUtil::getAbsPath(parameter.outputFile);
        } else {  // Neither output directory nor file specified
            if (parameter.inputFile == "/dev/stdin") {
                parameter.outputFile = "/dev/stdout";
            } else {
                if (!PathUtil::fileExists(parameter.inputFile)) {
                    fprintf(stderr, "%s not exits or no permission to access.\n", parameter.inputFile.c_str());
                    return -1;
                }
            }
        }

        // Reference gene input, convert to absolute path
        if (!parameter.referenceGenic.empty()) {
            std::string abspath = PathUtil::getAbsPath(parameter.referenceGenic);
            if (abspath.empty()) {
                fprintf(stderr, "%s not exits or no permission to access.\n", parameter.inputFile.c_str());
                return -1;
            }
            parameter.referenceGenic = abspath;
        }

        return 0;
    }

    int32_t afterConstructProc() {
        // If file exists and force overwrite is configured, delete output file
        if (PathUtil::fileExists(parameter.outputFile) && parameter.isOverwriteOutFile) {
            PathUtil::removeFile(parameter.outputFile);
        }
        ConfigManager::getInstance().init(parameter);
        return 0;
    }

    int32_t check() {
        return checkProc();
    }

    virtual int32_t checkProc() {
        if (parameter.outputFile != "/dev/stdout") {
            // Output file exists, remind to use -f to force overwrite
            if (PathUtil::fileExists(parameter.outputFile) && !parameter.isOverwriteOutFile) {
                fprintf(stdout, "%s exits, use -f to force overwrite.\n", parameter.outputFile.c_str());
                return -1;
            }

            // File specified and overwrite, check write permission
            if (PathUtil::fileExists(parameter.outputFile)) {
                if (!PathUtil::fileWriteable(parameter.outputFile)) {
                    fprintf(stdout,"no permisstion to write %s.\n", parameter.outputFile.c_str());
                    return -1;
                }
            } else {  // File doesn't exist, need to check directory permission
                // At this point, output file has been converted to absolute path, must include filename
                char separator = std::filesystem::path::preferred_separator;
                int pos = parameter.outputFile.find_last_of(separator);
                if (pos == -1) {
                    fprintf(stderr, "invalid file path, %s.\n", parameter.outputFile.c_str());
                    return -1;
                }

                std::string path = parameter.outputFile.substr(0, pos + 1);
                if (!PathUtil::fileExists(path)) {
                    fprintf(stderr, "%s not exists.\n", parameter.outputFile.c_str());
                    return -1;
                }

                if (!PathUtil::fileWriteable(path)) {
                    fprintf(stderr, "no permssion to access %s.\n", parameter.outputFile.c_str());
                    return -1;
                }
            }
        }

        if (parameter.inputFile != "/dev/stdin") {
            if (!PathUtil::fileExists(parameter.inputFile)) {
                fprintf(stderr, "%s not exits.\n", parameter.inputFile.c_str());
                return -1;
            }

            if (!PathUtil::fileReadble(parameter.inputFile)) {
                fprintf(stderr, "no permission to access %s.\n", parameter.inputFile.c_str());
                return -1;
            }
        }

        if (parameter.compressLevel < 1 || parameter.compressLevel > 9) {
            fprintf(stderr, "Compress level is invalid.\n");
            return -1;
        }

        if (parameter.logLevel < 1 || parameter.logLevel > 6) {
            fprintf(stderr, "Log level is invalid.\n");
            return -1;
        }

        if (!parameter.referenceGenic.empty()) {
            if (!PathUtil::fileExists(parameter.referenceGenic)) {
                fprintf(stderr, "Reference genic(%s) is not exists.\n", parameter.referenceGenic.c_str());
                return -1;
            }
        }

        return 0;
    }

    virtual int32_t startEngine() {
        return 0;
    }
};

class CompressCmdProc : public CommandProc {
public:
    CompressCmdProc(PbgzParameter& para) : CommandProc(para) {}
    
    int32_t reconstructPorc() override { 
        if (0 != CommandProc::reconstructPorc()) {
            return -1;
        }

        if (parameter.outputFile.empty() && parameter.outputDir.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + PathUtil::getFileName(parameter.inputFile) + ".pbgz";
        }
        return 0;
    }
    
    virtual int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        if ("/dev/stdout" != parameter.outputFile){
            if (PathUtil::isDir(parameter.outputFile)) {
                fprintf(stdout, "Output file cannot be directory.\n");
                return -1;
            }
            if (!PathUtil::suffixCheck(parameter.outputFile, "pbgz")) {
                fprintf(stdout, "Output file must be pbgz format.\n");
                return -1;
            }
        }
        if (!parameter.refeGenePos.empty()) {
            fprintf(stdout, "Reference position is not support in compress scenario.\n");
            return -1;
        }
        return 0;
    }
    
    int32_t startEngine() override {
        // 实现启动压缩引擎的逻辑
        engine = new CompressEngine(parameter);
        if (engine->init() != 0) {
            LOG_ERROR("Compress engine init error.");
            return -1;
        }

        int32_t ret = engine->start();
        return ret;
    }
};

class DecompressCmdProc : public CommandProc {
public:
    DecompressCmdProc(PbgzParameter& para) : CommandProc(para) {}
    
    int32_t reconstructPorc() override {
        if (0 != CommandProc::reconstructPorc()) {
            return -1;
        }
        if (parameter.outputFile.empty() && parameter.outputDir.empty()) { 
            // Decompression scenario, remove pbgz suffix from file
            std::string name = PathUtil::getFileName(parameter.inputFile);
            if (name.length() <= 5) {
                fprintf(stderr, "input file %s invalid.", name.c_str());
                return -1;
            }
            if (!PathUtil::suffixCheck(name, ".pbgz")) {
                fprintf(stdout, "Decompress file must be pbgz format.\n");
                return -1;
            }
            name.resize(name.length() - 5);
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + name;
        }

        return 0;
    }
    
    int32_t checkProc() override { 
        if (0 != CommandProc::checkProc()) {
            return -1;
        }
        if ("/dev/stdin" != parameter.inputFile) {
            if (PathUtil::isDir(parameter.inputFile)) {
                fprintf(stdout, "Decompress file cannot be directory.\n");
                return -1;
            }
            
            if (!PathUtil::suffixCheck(parameter.inputFile, "pbgz")) {
                fprintf(stdout, "Decompress file must be pbgz format.\n");
                return -1;
            }
        }
        return 0;
    }
    

    int32_t startEngine() override {
        // 实现启动解压缩引擎的逻辑
        engine = new DecompressEngine(parameter);
        if (engine->init() != 0) {
            return -1;
        }
        
        int32_t ret = engine->start();
        return ret;
    }
};

class IndexCmdProc : public CommandProc {
public:
    IndexCmdProc(PbgzParameter& para) : CommandProc(para) {}
    
    int32_t reconstructPorc() override {
        if (0 != CommandProc::reconstructPorc()) {
            return -1;
        }
        return 0;
    }
    
    int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        return -1;
    }
    
    int32_t startEngine() override {
        // 实现启动索引引擎的逻辑
        engine = new IndexEngine(parameter);
        if (engine->init() != 0) {
            return -1;
        }
        return engine->start();
    }
};

using CommandHandler = std::function<CommandProc*(PbgzParameter&)>;

// 二级命令结构
struct SubCommand {
    std::string name;               // 命令名称（如 "compress"）
    std::string description;        // 命令描述
    std::vector<char> args;          // 该命令支持的选项
    CommandHandler handler;
};

static std::vector<SubCommand> subCommands = { 
    {   
        "compress", 
        "Compress file to pbgz format",
        {'h', 'v', 'o', 'O', 'f', 'r', 'n', 't', 'l', 'e', 'i', 'g', 'G', 'h'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<CompressCmdProc>(para);
        },
    },
    {   
        "decompress", 
        "Decompress file from pbgz file",
        {'h', 'v', 'o', 'O', 'f', 'r', 't', 'e', 'p', 'g', 'G', 'h'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<DecompressCmdProc>(para);
        },
    },
    {   
        "index", 
        "Create index file from pbgz file",
        {'h', 'v', 'o', 'O', 'f', 't', 'g', 'G', 'h'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<IndexCmdProc>(para);
        },
    },
};

void printUsage(const std::string& subCommandName = "") {
#define ALIGN " "
#define CMD_ARG(f, s, d, p) \
    fprintf(f, "%2s-%c,%s--%-12s%s\n", ALIGN, s, ALIGN, d, p)

    FILE *fp = stdout;
    fprintf(fp, "pbgz: %s\n", PbgzManager::getInstance().getVersion().c_str());
    
    if (subCommandName.empty()) {
        // 显示通用帮助信息
        fprintf(fp, "Usage: pbgz subCommand [FILE] [OPTION]\n\n");
        fprintf(fp, "Available subcommands:\n");
        for (auto& cmd : subCommands) {
            fprintf(fp, "  %-12s %s\n", cmd.name.c_str(), cmd.description.c_str());
        }
        fprintf(fp, "\nUse 'pbgz help <subcommand>' for detailed information on a specific subcommand.\n");
        
        fprintf(fp, "\nExamples:\n");
        fprintf(fp, "  pbgz compress human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
        fprintf(fp, "  pbgz decompress human.fq.gz.pbgz\n");
        fprintf(fp, "  pbgz index human.fq.gz.pbgz\n\n");
    } else {
        // 显示特定子命令的帮助信息
        SubCommand* targetCmd = nullptr;
        for (auto& cmd : subCommands) {
            if (cmd.name == subCommandName) {
                targetCmd = &cmd;
                break;
            }
        }
        
        if (!targetCmd) {
            fprintf(fp, "Error: Unknown subcommand '%s'\n\n", subCommandName.c_str());
            fprintf(fp, "Available subcommands:\n");
            for (auto& cmd : subCommands) {
                fprintf(fp, "  %-12s %s\n", cmd.name.c_str(), cmd.description.c_str());
            }
            fprintf(fp, "\nUse 'pbgz help' to see all available subcommands.\n\n");
            return;
        }
        
        fprintf(fp, "Usage: pbgz %s [FILE] [OPTION]\n\n", subCommandName.c_str());
        fprintf(fp, "Description: %s\n\n", targetCmd->description.c_str());
        fprintf(fp, "Options for %s:\n", subCommandName.c_str());
        
        // 显示该子命令支持的选项
        for (auto &a : pbgzArgs) {
            bool supported = false;
            for (char supportedArg : targetCmd->args) {
                if (a.argShort == supportedArg) {
                    supported = true;
                    break;
                }
            }
            if (supported) {
                CMD_ARG(fp, a.argShort, a.argLong.c_str(), a.argDescribe.c_str());
            }
        }
        
        // 显示示例
        if (subCommandName == "compress") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz compress human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
            fprintf(fp, "  pbgz compress input.fastq -O /output/directory/ -t 8 -l 6\n\n");
        } else if (subCommandName == "decompress") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz decompress human.fq.gz.pbgz\n");
            fprintf(fp, "  pbgz decompress input.pbgz -o output.fastq -p chr1:1000-2000\n\n");
        } else if (subCommandName == "index") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz index human.fq.gz.pbgz\n");
            fprintf(fp, "  pbgz index input.pbgz -O /output/directory/\n\n");
        }
    }
}

void printVersion() {
    fprintf(stdout, "pbgz: %s\n", PbgzManager::getInstance().getVersion().c_str());
}

int main(int argc, char** argv) {
    if (argc == 1) {
        printUsage();
        return 0;
    }

    // 检查第一个参数是否是 help 或 version
    std::string firstArg = argv[1];
    if (firstArg == "help") {
        if (argc >= 3) {
            // 显示特定子命令的帮助信息
            printUsage(argv[2]);
        } else {
            // 显示通用帮助信息
            printUsage();
        }
        return 0;
    }

    if (firstArg == "version") {
        printVersion();
        return 0;
    }

    bool foundSubCommand = false;
    SubCommand* selectedSubCommand = nullptr;
    
    for (auto& subCmd : subCommands) {
        if (firstArg == subCmd.name) {
            foundSubCommand = true;
            selectedSubCommand = &subCmd;
            break;
        }
    }

    if (!foundSubCommand) {
        printUsage();
        return 0;
    }

    // 使用子命令的新逻辑
    PbgzParameter parameter;

    // 构建选项字符串，只包含子命令支持的选项
    std::string argOption = "";
    option longopts[pbgzArgs.size() + 1];
    int optIndex = 0;
    
    for (uint32_t i = 0; i < pbgzArgs.size(); ++i) {
        // 检查这个选项是否被当前子命令支持
        bool supported = false;
        for (char supportedArg : selectedSubCommand->args) {
            if (pbgzArgs[i].argShort == supportedArg) {
                supported = true;
                break;
            }
        }
        
        if (supported) {
            longopts[optIndex] = {pbgzArgs[i].argLong.c_str(), pbgzArgs[i].hasFlag, nullptr, pbgzArgs[i].argShort};
            argOption += pbgzArgs[i].argShort;
            switch (pbgzArgs[i].hasFlag) 
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
            optIndex++;
        }
    }
    longopts[optIndex] = {nullptr, 0, nullptr, 0};

    // 重置optind以跳过子命令名称
    optind = 2;
    
    int opt = 0;
    while((opt = getopt_long(argc, argv, argOption.c_str(), longopts, NULL)) != -1) {
        switch (opt) {
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
        case 'p':{
            parameter.refeGenePos = optarg;
            break;
        }
        case 'i':
            parameter.isMakeIndex = true;
            break;
        case '?':
        case 'h':
            parameter.showHelp = true;
            break;
        default:
            fprintf(stdout, "Use pbgz help %s to show usage.\n\n", selectedSubCommand->name.c_str());
            return 0;
        }
    }

    if (parameter.showHelp) {
        printUsage(selectedSubCommand->name);
        return 0;
    }

    // 处理剩余的参数（输入文件等）
    if (optind < argc) {
        if (argc - optind > 1) {
            fprintf(stdout, "Only one file is allowed");
            return 0;
        }
        parameter.inputFile = argv[optind];
    }

    // 创建并执行子命令处理器
    CommandProc* processor = selectedSubCommand->handler(parameter);
    if (!processor) {
        LOG_ERROR("Failed to create command processor");
        return -1;
    }

    int result = 0;
    if (processor->reconstruct() != 0) {
        LOG_ERROR("Command reconstruction failed");
        return -1;
    } 
    
    if (processor->check() != 0) {
        LOG_ERROR("Command parameter check failed");
        return -1;
    } 
    
    if (processor->startEngine() != 0) {
        LOG_ERROR("Command start failed");
        return -1;
    }

    MemoryUtil::safeDeleteClass(processor);
    return 0;
}
