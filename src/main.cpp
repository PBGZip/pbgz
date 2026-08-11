/*
 * main.cpp - Source file for pbgz project
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
#include "sort_engine.h"


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
    {'N', "ni", required_argument, "<NI> use an index built by 'pbgz makeref' instead of re-reading the FASTA. -r is still required and decides which reference is used"},
    {'n', "refunpack", no_argument, "unpack reference to pbgz file, so reference gene is needed when decompress"},
    {'t', "threads", required_argument, "<number> specify number of threads, default is CPUS. use for compress/decompress"},
    {'l', "level", required_argument, "<1-9> specify compress level, default is 5, 1 is fast, 9 is best"},
    {'e', "remove", no_argument, "if compress succeed, will remove origin file, else this option is invalid. use for compress"},
    {'h', "help",  no_argument, "show help"},
    {'g', "loglevel", required_argument, "<1-6> sepcify log level, default 6. debug(1), info(2), warning(3), error(4), fatal(5), off(6)"},
    {'G', "logfile", required_argument, "sepcify log file"},
    {'p', "position", required_argument, "sepecify the reference gene posision"},
    {'i', "index", no_argument, "make index"},
    {'s', "stat ", no_argument, "show statistics data"},
    {'v', "verbose", no_argument, "show codec selection per field and compression speed (MB/s)"},
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
            parameter.inputFile = STDIN;
        }
        if ("-" == parameter.outputFile ) {
            parameter.outputFile = STDOUT;
        }
        parameter.inputFile = PathUtil::getAbsPath(parameter.inputFile);

        if (!parameter.outputDir.empty() && !parameter.outputFile.empty()) {  // Output directory and file specified simultaneously
            fprintf(stderr, "Output directory and filename cannot be specified simultaneously.\n");
            return -1;
        } else if(!parameter.outputDir.empty()) {   // Only output directory specified
            if (parameter.inputFile == STDIN) {
                fprintf(stderr, "Specify an output filename, or do not specify both the filename and the directory.\n");
                return -1;
            } else {
                // Output filename is input filename with suffix (.pbgz)
                std::string inFileName = PathUtil::getFileName(parameter.inputFile);
                if (inFileName.empty()) {
                    fprintf(stderr, "%s not exits or no permission to access.\n", parameter.inputFile.c_str());
                    return -1;
                }

                // Append a '/' if the input directory path doesn't end with one
                if (!(parameter.outputDir.back() == '/' || parameter.outputDir.back() == '\\')) {
                    std::filesystem::path outputDir = parameter.outputDir;
                    outputDir = outputDir / "";
                    parameter.outputDir = outputDir.string();
                }

                std::string outPath = PathUtil::getFilePath(parameter.outputDir);
                if (!PathUtil::fileExists(outPath)) {
                    // Try to create directory if it doesn't exist
                    outPath = PathUtil::createDir(parameter.outputDir);
                    if (outPath.empty()) {
                        fprintf(stderr, "%s not exits or no permission to access.\n", parameter.outputDir.c_str());
                        return -1;
                    }
                }
            }
        } else if (!parameter.outputFile.empty()) {     // File output specified
            // Convert output file to absolute path
            parameter.outputFile = PathUtil::getAbsPath(parameter.outputFile);
        } else {  // Neither output directory nor file specified
            if (parameter.inputFile == STDIN) {
                parameter.outputFile = STDOUT;
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
        ConfigManager::getInstance().init(parameter);
        return 0;
    }

    int32_t check() {
        return checkProc();
    }

    virtual int32_t checkProc() {
        if (parameter.outputFile != STDOUT) {
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
                    fprintf(stderr, "%s not exists.\n", path.c_str());
                    return -1;
                }

                if (!PathUtil::fileWriteable(path)) {
                    fprintf(stderr, "no permssion to access %s.\n", parameter.outputFile.c_str());
                    return -1;
                }
            }
        }

        if (parameter.inputFile != STDIN) {
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
            if (Reference::isNiFile(parameter.referenceGenic)) {
                fprintf(stderr, "%s is a pbgz index, pass it with --ni and give the FASTA to -r.\n",
                        parameter.referenceGenic.c_str());
                return -1;
            }
        }

        /* 索引只是加速手段, -r 才决定用哪个参考; 索引不可信时还得靠 -r 退回去重建。 */
        if (!parameter.niIndexFile.empty()) {
            if (parameter.referenceGenic.empty()) {
                fprintf(stderr, "--ni requires -r: the FASTA it was built from decides which reference is used.\n");
                return -1;
            }
            if (!PathUtil::fileExists(parameter.niIndexFile)) {
                fprintf(stderr, "Index file(%s) is not exists.\n", parameter.niIndexFile.c_str());
                return -1;
            }
        }

        return 0;
    }

    int32_t afterCheck() {
        return afterCheckProc();
    }

    virtual int32_t afterCheckProc() {
        // If file exists and force overwrite is configured, delete output file
        if (PathUtil::fileExists(parameter.outputFile) && parameter.isOverwriteOutFile) {
            PathUtil::removeFile(parameter.outputFile);
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
            std::string inputFile = PathUtil::getFileName(parameter.inputFile);
            if (inputFile.length() > 3 && PathUtil::suffixCheck(inputFile, ".gz")) {
                inputFile = inputFile.substr(0, inputFile.length() - 3);
            }
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + inputFile + ".pbgz";
        }

        if (!parameter.outputDir.empty() && parameter.outputFile.empty()) {
            std::string inputFile = PathUtil::getFileName(parameter.inputFile);
            if (inputFile.length() > 3 && PathUtil::suffixCheck(inputFile, ".gz")) {
                inputFile = inputFile.substr(0, inputFile.length() - 3);
            }
            parameter.outputFile = PathUtil::getFilePath(parameter.outputDir) + inputFile + ".pbgz";
        }
        return 0;
    }

    virtual int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        if (STDOUT != parameter.outputFile){
            if (PathUtil::isDir(parameter.outputFile)) {
                fprintf(stdout, "Output file cannot be directory.\n");
                return -1;
            }
            if (!PathUtil::suffixCheck(parameter.outputFile, "pbgz")) {
                fprintf(stdout, "Output file must be pbgz format.\n");
                return -1;
            }
        }
        return 0;
    }

    int32_t startEngine() override {
        // Implement compression engine startup logic
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
        if (parameter.outputFile.empty()) {
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
            if (parameter.outputDir.empty()) {
                parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + name;
            } else {
                parameter.outputFile = PathUtil::getAbsPath(parameter.outputDir) + name;
            }

            if (parameter.isDecToGZ) {
                parameter.outputFile = parameter.outputFile + ".gz";
            }
        }

        return 0;
    }

    int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }
        if (STDIN != parameter.inputFile) {
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
        // Implement decompression engine startup logic
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

        if (!parameter.outputFile.empty()) {
            return 0;
        }

        std::string inputFile = PathUtil::getFileName(parameter.inputFile);
        if (inputFile.length() <= 5) {
            fprintf(stderr, "input file %s invalid.", inputFile.c_str());
            return -1;
        }
        if (!PathUtil::suffixCheck(inputFile, ".pbgz")) {
            fprintf(stdout, "Decompress file must be pbgz format.\n");
            return -1;
        }

        inputFile = inputFile + ".pbgzi";
        if (parameter.isDecToGZ) {
            inputFile = inputFile + ".gz";
        }

        if (parameter.outputFile.empty() && parameter.outputDir.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + inputFile;
        } else if (!parameter.outputDir.empty() && parameter.outputFile.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.outputDir) + inputFile;
        }

        return 0;
    }

    int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        if (STDIN != parameter.inputFile) {
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
        // Implement index engine startup logic
        engine = new IndexEngine(parameter);
        if (engine->init() != 0) {
            return -1;
        }
        return engine->start();
    }
};

class SortCmdProc : public CommandProc {
public:
    SortCmdProc(PbgzParameter& para) : CommandProc(para) {}

    int32_t reconstructPorc() override {
        if (0 != CommandProc::reconstructPorc()) {
            return -1;
        }

        if (!parameter.outputFile.empty()) {
            return 0;
        }

        std::string inputFile = PathUtil::getFileName(parameter.inputFile);
        if (inputFile.length() > 3 && PathUtil::suffixCheck(inputFile, ".gz")) {
            inputFile = inputFile.substr(0, inputFile.length() - 3);
        }

        // Split the string into two parts based on the position of the last . in inputFile
        size_t lastDotPos = inputFile.find_last_of('.');
        if (lastDotPos != std::string::npos && lastDotPos > 0) {
            std::string fileNamePrefix = inputFile.substr(0, lastDotPos);
            std::string  fileExtension = inputFile.substr(lastDotPos);
            inputFile = fileNamePrefix + ".sort" + fileExtension;
        } else {
            inputFile = inputFile + ".sort";
        }

        if (parameter.isDecToGZ) {
            inputFile = inputFile + ".gz";
        }

        if (parameter.outputFile.empty() && parameter.outputDir.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + inputFile;
        } else if (!parameter.outputDir.empty() && parameter.outputFile.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.outputDir) + inputFile;
        }

        return 0;
    }

    int32_t checkProc() override {
        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        return 0;
    }

    int32_t startEngine() override {
        // Implement index engine startup logic
        engine = new SortEngine(parameter);
        if (engine->init() != 0) {
            return -1;
        }
        return engine->start();
    }
} ;

/*
 * 参考索引只在用户显式要求时生成, 且落在用户指定的路径; 之后压缩要复用, 也得由用户
 * 显式把它交给 -r。不做隐式缓存: 悄悄留下的索引一旦被人改动, 后续压缩会静默压错,
 * 而用户手上没有任何线索。显式指定意味着用户认得这个文件, 也认这份风险。
 */
class MakeRefCmdProc : public CommandProc {
public:
    MakeRefCmdProc(PbgzParameter& para) : CommandProc(para) {}

    int32_t reconstructPorc() override {
        if (0 != CommandProc::reconstructPorc()) {
            return -1;
        }

        if (!parameter.outputFile.empty()) {
            return 0;
        }

        std::string inputFile = PathUtil::getFileName(parameter.inputFile);
        if (inputFile.empty()) {
            fprintf(stderr, "input file %s invalid.\n", parameter.inputFile.c_str());
            return -1;
        }
        inputFile += ".ni";

        if (parameter.outputDir.empty()) {
            parameter.outputFile = PathUtil::getFilePath(parameter.inputFile) + inputFile;
        } else {
            parameter.outputFile = PathUtil::getFilePath(parameter.outputDir) + inputFile;
        }
        return 0;
    }

    int32_t checkProc() override {
        if (STDIN == parameter.inputFile) {
            fprintf(stdout, "Reference file must be a regular file, not a pipe.\n");
            return -1;
        }

        if (0 != CommandProc::checkProc()) {
            return -1;
        }

        if (PathUtil::isDir(parameter.inputFile)) {
            fprintf(stdout, "Reference file cannot be directory.\n");
            return -1;
        }

        if (Reference::isNiFile(parameter.inputFile)) {
            fprintf(stdout, "%s is already a reference index.\n", parameter.inputFile.c_str());
            return -1;
        }
        return 0;
    }

    int32_t startEngine() override {
        Reference reference(parameter.inputFile, parameter.threadNum);
        if (!reference.makeNiFile(parameter.outputFile)) {
            fprintf(stderr, "Make reference index %s failed.\n", parameter.outputFile.c_str());
            return -1;
        }
        fprintf(stderr, "Reference index created: %s\n", parameter.outputFile.c_str());
        return 0;
    }
};

using CommandHandler = std::function<CommandProc*(PbgzParameter&)>;
// Subcommand structure
struct SubCommand {
    std::string name;               // Subcommand name, e.g., "compress"
    std::string description;        // Subcommand description
    std::vector<char> args;          // Options supported by subcommand
    CommandHandler handler;
};

static std::vector<SubCommand> subCommands = {
    {
        "compress",
        "Compress file to pbgz format",
        {'h', 'v', 'o', 'O', 'f', 'r', 'N', 'n', 't', 'l', 'e', 'i', 'g', 'G', 's'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<CompressCmdProc>(para);
        },
    },
    {
        "decompress",
        "Decompress file from pbgz file",
        {'h', 'v', 'o', 'O', 'f', 'r', 'N', 't', 'e', 'p', 'g', 'G', 'z', 'l'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<DecompressCmdProc>(para);
        },
    },
    {
        "index",
        "Create index file from pbgz file",
        {'h', 'v', 'o', 'O', 'f', 't', 'g', 'G'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<IndexCmdProc>(para);
        },
    },
    {
        "sort",
        "Sort sam file",
        {'h', 'v', 'o', 'O', 'f', 't', 'g', 'G', 'z'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<SortCmdProc>(para);
        },
    },
    {
        "makeref",
        "Build a reusable index from a reference FASTA",
        {'h', 'o', 'O', 'f', 't', 'g', 'G'},
        [](PbgzParameter& para) {
            return MemoryUtil::safeNewClass<MakeRefCmdProc>(para);
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
        // Display general help information
        fprintf(fp, "Usage: pbgz subCommand [FILE] [OPTION]\n\n");
        fprintf(fp, "Available subcommands:\n");
        for (auto& cmd : subCommands) {
            fprintf(fp, "  %-12s %s\n", cmd.name.c_str(), cmd.description.c_str());
        }
        fprintf(fp, "\nUse 'pbgz help <subcommand>' for detailed information on a specific subcommand.\n");

        fprintf(fp, "\nExamples:\n");
        fprintf(fp, "  pbgz compress human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
        fprintf(fp, "  pbgz decompress human.fq.gz.pbgz\n");
        fprintf(fp, "  pbgz index human.sam.gz.pbgz\n\n");
    } else {
        // Display specific subcommand help information
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

        // Display options supported by subcommand
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

        // Display examples
        if (subCommandName == "compress") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz compress human.fq.gz -o /path/human.fq.gz.pbgz -r /path/ucsc.hg19.fa\n");
            fprintf(fp, "  pbgz compress input.fastq -O /output -t 8 -l 6\n\n");
        } else if (subCommandName == "decompress") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz decompress human.fq.gz.pbgz\n");
            fprintf(fp, "  pbgz decompress input.pbgz -O /output \n\n");
        } else if (subCommandName == "index") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz index human.sam.gz.pbgz\n");
            fprintf(fp, "  pbgz index input.pbgz -O /output\n\n");
        } else if (subCommandName == "sort") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz sort human.sam -O /output\n");
            fprintf(fp, "  pbgz sort human.sam.gz -O /output -z\n\n");
        } else if (subCommandName == "makeref") {
            fprintf(fp, "\nExample:\n");
            fprintf(fp, "  pbgz makeref ucsc.hg19.fa -o /path/ucsc.hg19.ni\n");
            fprintf(fp, "  pbgz compress input.sam -o out.pbgz -r ucsc.hg19.fa --ni /path/ucsc.hg19.ni\n\n");
            fprintf(fp, "The index is never created or picked up implicitly, and it never replaces -r:\n");
            fprintf(fp, "-r names the reference, --ni only saves parsing it again. The index must carry\n");
            fprintf(fp, "that same file name and pass its own checksum, otherwise pbgz warns and builds\n");
            fprintf(fp, "the index from the FASTA, so the reference in use is always the one you chose.\n\n");
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

    // Check if the first argument is help or version
    std::string firstArg = argv[1];
    if (firstArg == "help") {
        if (argc >= 3) {
            // Display specific subcommand help information
            printUsage(argv[2]);
        } else {
            // Display general help information
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

    // Use command line argument processing logic
    PbgzParameter parameter;

    // Build option string containing only options supported by current subcommand
    std::string argOption = "";
    option longopts[pbgzArgs.size() + 1];
    int optIndex = 0;

    for (uint32_t i = 0; i < pbgzArgs.size(); ++i) {
        // Check if subcommand option is supported by current subcommand
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

    // Reset optind to skip program name and subcommand
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
        case 'N':{
            parameter.niIndexFile = optarg;
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
        case 's':
            parameter.showStat = true;
            break;
        case 'v':
            parameter.verbose = true;
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

    // Process remaining arguments (such as input files, etc.)
    if (optind < argc) {
        if (argc - optind > 1) {
            fprintf(stdout, "Only one file is allowed");
            return 0;
        }
        parameter.inputFile = argv[optind];
    }

    // Create and execute command processor
    CommandProc* processor = selectedSubCommand->handler(parameter);
    if (!processor) {
        LOG_ERROR("Failed to create command processor");
        return -1;
    }

    if (processor->reconstruct() != 0) {
        LOG_ERROR("Command reconstruction failed");
        return -1;
    }

    if (processor->check() != 0) {
        LOG_ERROR("Command parameter check failed");
        return -1;
    }

    if (processor->afterCheck() != 0) {
        LOG_ERROR("Command parameter after check processs failed");
        return -1;
    }

    int32_t ret = processor->startEngine();
    if (ret != 0) {
        LOG_ERROR("Command start failed");
        /* 异常退出：清理并删除已创建的输出文件（如缺参考基因导致的失败）。 */
        PbgzManager::getInstance().exitProc(ret, "pbgz start failed");
    }

    MemoryUtil::safeDeleteClass(processor);
    return 0;
}
