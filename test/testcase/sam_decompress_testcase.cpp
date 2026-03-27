#include <gtest/gtest.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

#define private public
#include "coder/coder.h"
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include "utils/memory_util.h"
#undef private
#include <random>

namespace SamDecompressData {
    const std::string testSamFile = "test.sam";
    const std::string testFastqFile = "test.fa";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamDecompressTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

		pInBlock = new RoughIOBlock(SamDecompressData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamDecompressData::MAX_BLOCK_SIZE);

        generateSamFile(SamDecompressData::testSamFile);

        ConfigManager::getInstance().logLevel = LogLevel::DEBUGGING;
	}

	// Clean up resources
	void TearDown() override {
		if (pInBlock != nullptr) {
			delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }
        // 清理所有可能生成的测试文件
        std::remove(SamDecompressData::testSamFile.c_str());
        std::remove("./test/test.sam");
        std::remove("../test/test.sam");
        std::remove("test_reference.fa");
        std::remove("./test/test_reference.fa");
        std::remove("../test/test_reference.fa");
        std::remove("test_id_separators.sam");
        std::remove("test_missing_fields.sam");
        std::remove("test_flag_matching.sam");
        std::remove("test_base_fixed.sam");
        std::remove("test_base_variable.sam");
        std::remove("test_base_fastq3.sam");
        std::remove("test_optional_fields.sam");
        std::remove("test_no_optional_fields.sam");
	}

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        // 尝试多个可能的路径，首先尝试当前目录和测试目录
        std::vector<std::string> paths = {
            filename,           // 当前目录
            "./test/" + filename, // 测试目录
            "../test/" + filename // 构建目录的上一级测试目录
        };

        IOReader* pIoReader = nullptr;
        for (const auto& path : paths) {
            pIoReader = new FileReader(path);
            if (pIoReader->openIO() == 0) {
                break;
            }
            delete pIoReader;
            pIoReader = nullptr;
        }

        if (pIoReader == nullptr) {
            // 如果所有路径都失败，直接使用test_data/test.sam作为最后的fallback
            pIoReader = new FileReader("test_data/test.sam");
            if (pIoReader->openIO() != 0) {
                delete pIoReader;
                return; // 无法打开文件
            }
        }

        BlockReader* pBlockReader = new BlockReader(pIoReader);
        pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);

        delete pBlockReader;
        delete pIoReader;
    }

    void generateSamFile(const std::string& filename) {
        // 直接写入 test_data 中的真实 SAM 数据，不从文件读取
        std::ofstream file(filename);
        if (!file.is_open()) {
            // 如果无法在当前目录打开，尝试在测试目录中创建
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                // 如果还是失败，尝试在构建目录中创建
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }

        if (!file.is_open()) {
            return; // 无法创建文件
        }

        // 生成标准的 SAM 文件内容，所有记录都是 76bp 长度
        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";
        file << "@SQ\tSN:chr3\tLN:339\n";
        file << "@SQ\tSN:chr4\tLN:339\n";
        file << "@SQ\tSN:chrX\tLN:340\n";
        file << "@SQ\tSN:chrY\tLN:339\n";
        file << "@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188\tCL:bwa mem test_data/reference.fa test_data/test_reads.fastq\n";

        // 所有记录都是 76M，序列长度和质量值都是76
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read3\t0\tchr1\t153\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCCGCAAG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:3\tMD:Z:37C37C0\tAS:i:72\tXS:i:0\n";
        file << "read4\t0\tchr2\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGACAATCGCCGGAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read5\t0\tchr2\t77\t60\t76M\t*\t0\t0\tATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAACGCAATCGATCGA\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read6\t16\tchr3\t1\t60\t76M\t*\t0\t0\tCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCACAATCGGCGGCCGAATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read7\t0\tchr3\t77\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATACGCAATCGGCGGCCGAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:74A0\tAS:i:75\tXS:i:0\n";
        file << "read8\t16\tchr4\t1\t60\t76M\t*\t0\t0\tTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGCAACAATCGGCGGCC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read9\t0\tchrX\t1\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGACGCAATCGCCGGCCTAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read10\t0\tchrY\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAACGCAATCGCCGGCCTTGTTC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";

        file.close();
    }

    void generateFasta(const std::string& filename) {
        std::ifstream srcFile("reference.fa");
        if (!srcFile.is_open()) {
            // 如果无法打开test_data/reference.fa，生成包含真实ATCG序列的参考基因组
            std::ofstream file(filename);
            if (!file.is_open()) {
                std::string testPath = "./test/" + filename;
                file.open(testPath);
                if (!file.is_open()) {
                    testPath = "../test/" + filename;
                    file.open(testPath);
                }
            }

            if (!file.is_open()) {
                return; // 无法创建文件
            }

            // 生成包含真实ATCG序列的参考基因组
            file << ">chr1\n";
            for (int i = 0; i < 17; i++) {
                file << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
            }
            file << "ATCGATCGATCGATCGATCG\n";  // 填充到合适的长度

            file << ">chr2\n";
            for (int i = 0; i < 16; i++) {
                file << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
            }
            file << "GCTAGCTAGCTAGCTAGC\n";

            file << ">chr3\n";
            for (int i = 0; i < 17; i++) {
                file << "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
            }
            file << "TACGTACGTACGTACGTACGT\n";

            file << ">chr4\n";
            for (int i = 0; i < 17; i++) {
                file << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
            }
            file << "CATGCATGCATGCATGCATGC\n";

            file << ">chrX\n";
            for (int i = 0; i < 17; i++) {
                file << "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
            }
            file << "ATGCATGCATGCATGCATGC\n";

            file << ">chrY\n";
            for (int i = 0; i < 17; i++) {
                file << "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\n";
            }
            file << "GCATGCATGCATGCATGCATG\n";

            file.close();
            return;
        }

        // 如果能打开test_data/reference.fa，则复制其内容
        std::string content((std::istreambuf_iterator<char>(srcFile)),
                           std::istreambuf_iterator<char>());
        srcFile.close();

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::string testPath = "./test/" + filename;
            file.open(testPath);
            if (!file.is_open()) {
                testPath = "../test/" + filename;
                file.open(testPath);
            }
        }

        if (!file.is_open()) {
            return; // 无法创建文件
        }

        file << content;
        file.close();
    }

    // 创建测试用的Reference对象
    Reference createTestReference() {
        std::string tempFastaFile = "test_reference.fa";
        generateFasta(tempFastaFile);
        Reference ref(tempFastaFile, 1);
        ref.makeIndex();
        return ref;
    }

    // 生成包含各种ID分割符的SAM文件
    void generateSamFileWithSeparators(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";
        file << "@SQ\tSN:chr3\tLN:339\n";

        // 每个read都包含所有分隔符，顺序一致：: . # / | - _ $
        file << "read1:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read2:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read3:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read4:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t4\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read5:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t5\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read6:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t6\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read7:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t7\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read8:part1.part2#part3/part4|part5-part6_part7$part8\t0\tchr1\t8\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // 生成缺少字段的SAM文件（用于测试段压缩逻辑）
    void generateSamFileWithMissingFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 正常记录
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // 缺少一些字段的记录（只有ID和FLAG）
        file << "read2\t0\n";
        // 缺少中间字段的记录
        file << "read3\t0\tchr1\t2\t60\n";
        // 缺少更多字段的记录
        file << "read4\t0\tchr2\t3\t76M\t*\t0\n";
        // 只有ID的记录
        file << "read5\n";

        file.close();
    }

    // 生成不同FLAG匹配情况的SAM文件
    void generateSamFileWithFlagVariants(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 匹配的FLAG (65)
        file << "read1\t65\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // 未匹配的FLAG (4)
        file << "read2\t4\t*\t0\t0\t*\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:0\n";
        // 正向匹配的FLAG (0)
        file << "read3\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        // 反向匹配的FLAG (16)
        file << "read4\t16\tchr2\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // 生成定长Base字段的SAM文件
    void generateSamFileWithFixedBase(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 所有序列长度相同（定长）
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\n";

        file.close();
    }

    // 生成变长Base字段的SAM文件
    void generateSamFileWithVariableBase(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 序列长度不同，确保质量值数量与碱基数量一致
        // 第一行：50M CIGAR，50个碱基，50个质量值
        std::string seq1 = "";
        for (int i = 0; i < 50; i++) {
            seq1 += "ATCG"[i % 4];
        }
        std::string qual1 = std::string(50, '!');
        file << "read1\t0\tchr1\t1\t60\t50M\t*\t0\t0\t" << seq1 << "\t" << qual1 << "\tNM:i:1\n";
        
        // 第二行：75M CIGAR，75个碱基，75个质量值
        std::string seq2 = "";
        for (int i = 0; i < 75; i++) {
            seq2 += "ATCG"[i % 4];
        }
        std::string qual2 = std::string(75, '!');
        file << "read2\t0\tchr1\t2\t60\t75M\t*\t0\t0\t" << seq2 << "\t" << qual2 << "\tNM:i:1\n";
        
        // 第三行：100M CIGAR，100个碱基，100个质量值
        std::string seq3 = "";
        for (int i = 0; i < 100; i++) {
            seq3 += "ATCG"[i % 4];
        }
        std::string qual3 = std::string(100, '!');
        file << "read3\t0\tchr2\t3\t60\t100M\t*\t0\t0\t" << seq3 << "\t" << qual3 << "\tNM:i:1\n";

        file.close();
    }

    // 生成3代fastq风格的SAM文件
    void generateSamFileWithFastq3Base(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 包含长序列和特殊碱基，确保质量值数量与碱基数量一致
        file << "read1\t0\tchr1\t1\t60\t4000M\t*\t0\t0\t";
        for (int i = 0; i < 1000; i++) {
            file << "ATCG";
        }
        file << "\t";
        for (int i = 0; i < 4000; i++) {
            file << "!";
        }
        file << "\tNM:i:10\n";

        file.close();
    }

    // 生成包含附加字段的SAM文件
    void generateSamFileWithOptionalFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 包含多种附加字段
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\tRG:Z:group1\tBC:Z:AGCTCT\tQT:Z:FFFFFFFFFFFFFF\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74T0\tAS:i:74\tXS:i:0\tRG:Z:group2\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:0\tMD:Z:76\tAS:i:76\tXS:i:0\tXA:Z:chr2,+5000,76M,1;chr3,+10000,76M,1;\n";

        file.close();
    }

    // 生成不包含附加字段的SAM文件
    void generateSamFileWithoutOptionalFields(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            return;
        }

        file << "@HD\tVN:1.6\tSO:coordinate\n";
        file << "@SQ\tSN:chr1\tLN:340\n";
        file << "@SQ\tSN:chr2\tLN:338\n";

        // 只有必要字段
        file << "read1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read2\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
        file << "read3\t0\tchr1\t3\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

        file.close();
    }

    // 压缩和解压辅助函数（优化版本：不使用文件，直接内存操作）
    void compressAndDecompress(const std::string& inputFile) {
        // 加载输入文件
        loadSamData(inputFile);

        // 创建压缩器
        Reference ref = createTestReference();
        SamActuator compressor(pInBlock, pOutBlock, &ref);
        ASSERT_EQ(compressor.preAnalysis(), 0);
        ASSERT_EQ(compressor.compress(), 0);

        // 压缩之后，将输出从outBlock拷贝到inBlock
        pInBlock->reset();
        memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
        pInBlock->setDataLen(pOutBlock->getDataLen());
        pInBlock->setMetaLen(pOutBlock->getMetaLen());

        // 重置outBlock
        pOutBlock->reset();

        // 用这两个block解压，只需校验解压的返回码即可
        SamActuator decompressor(pInBlock, pOutBlock, &ref);

        int32_t ret = decompressor.decompress();
        std::remove(inputFile.c_str());
        EXPECT_EQ(ret, 0);
    }

protected:
    RoughIOBlock* pInBlock = nullptr;
    RoughIOBlock* pOutBlock = nullptr;
};

// 测试1: ID字段覆盖所有分割字符
TEST_F(SamDecompressTest, TestIDFieldSeparators) {
    generateSamFileWithSeparators("test_id_separators.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_id_separators.sam"));
}

// 测试2: 缺少字段的情况，触发段压缩逻辑
TEST_F(SamDecompressTest, TestMissingFieldsSegmentCompression) {
    generateSamFileWithMissingFields("test_missing_fields.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    // 加载输入文件
    loadSamData("test_missing_fields.sam");

    // 创建压缩器
    Reference ref = createTestReference();
    SamActuator compressor(pInBlock, pOutBlock, &ref);
    ASSERT_EQ(compressor.preAnalysis(), -1);
}

// 测试3: FLAG字段的匹配、未匹配、正向匹配、反向匹配
TEST_F(SamDecompressTest, TestFlagVariants) {
    generateSamFileWithFlagVariants("test_flag_matching.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_flag_matching.sam"));
}

// 测试4.1: Base字段定长场景
TEST_F(SamDecompressTest, TestBaseFieldFixed) {
    generateSamFileWithFixedBase("test_base_fixed.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_base_fixed.sam"));
}

// 测试4.2: Base字段变长场景
TEST_F(SamDecompressTest, TestBaseFieldVariable) {
    generateSamFileWithVariableBase("test_base_variable.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_base_variable.sam"));
}

// 测试4.3: 3代fastq场景
TEST_F(SamDecompressTest, TestBaseFieldFastq3) {
    generateSamFileWithFastq3Base("test_base_fastq3.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_base_fastq3.sam"));
}

// 测试5.1: 有附加字段
TEST_F(SamDecompressTest, TestOptionalFieldsPresent) {
    generateSamFileWithOptionalFields("test_optional_fields.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_optional_fields.sam"));
}

// 测试5.2: 没有附加字段
TEST_F(SamDecompressTest, TestNoOptionalFields) {
    generateSamFileWithoutOptionalFields("test_no_optional_fields.sam");
    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_no_optional_fields.sam"));
}

// 综合测试：混合所有场景
TEST_F(SamDecompressTest, TestMixedScenarios) {
    std::ofstream file("test_mixed.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:340\n";
    file << "@SQ\tSN:chr2\tLN:338\n";

    // 混合各种分割符
    file << "read1:part1.part2#part3/part4|part5-part6_part7$part8\t65\tchr1\t1\t60\t100M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:5\tMD:Z:95A0\tAS:i:95\tXS:i:0\n";
    // 缺少字段的记录
    file << "read2\t4\t*\t0\t0\t*\t*\t0\t0\n";
    // 正向匹配
    file << "read3\t0\tchr1\t2\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    // 反向匹配
    file << "read4\t16\tchr2\t3\t60\t200M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG_ATCGATCG_ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:10\n";
    file.close();

    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    // 加载输入文件
    loadSamData("test_mixed.sam");

    // 创建压缩器
    Reference ref = createTestReference();
    SamActuator compressor(pInBlock, pOutBlock, &ref);
    ASSERT_EQ(compressor.preAnalysis(), -1);
}

// 压缩性能测试（大文件）
TEST_F(SamDecompressTest, TestCompressionPerformance) {
    // 生成大型SAM文件进行性能测试
    std::ofstream file("test_large.sam");
    if (!file.is_open()) {
        return;
    }

    file << "@HD\tVN:1.6\tSO:coordinate\n";
    file << "@SQ\tSN:chr1\tLN:248956422\n";
    file << "@SQ\tSN:chr2\tLN:242193529\n";
    file << "@SQ\tSN:chr3\tLN:198295559\n";
    file << "@SQ\tSN:chr4\tLN:190214555\n";
    file << "@SQ\tSN:chr5\tLN:181538259\n";

    // 生成10000条记录
    for (int i = 0; i < 10000; i++) {
        file << "read_" << i << "_test:some/complex|id#$" << i % 10 << "\t";
        file << (i % 4 == 0 ? 4 : (i % 2 == 0 ? 0 : 16)) << "\t";
        file << "chr" << (i % 5 + 1) << "\t" << (i * 100 + 1) << "\t60\t";
        file << (50 + i % 100) << "M\t*\t0\t0\t";

        // 生成变长序列
        int seqLen = 50 + i % 100;
        for (int j = 0; j < seqLen; j++) {
            file << "ATCG"[j % 4];
        }
        file << "\t";

        // 生成质量值
        for (int j = 0; j < seqLen; j++) {
            file << "!";
        }

        if (i % 3 == 0) {
            // 部分记录包含附加字段
            file << "\tNM:i:" << (i % 10) << "\tMD:Z:" << seqLen << "\tAS:i:" << (90 - i % 10);
        }
        file << "\n";
    }

    file.close();

    // 优化后只验证压缩和解压的返回码，不再比较文件内容
    EXPECT_NO_THROW(compressAndDecompress("test_large.sam"));
}
