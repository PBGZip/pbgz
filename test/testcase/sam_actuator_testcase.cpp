#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

#define private public
#include "sam_actuator.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#undef private
#include <random>

namespace SamTestData {
    const std::string testSamFile = "test.sam";
    const std::string testFastqFile = "test.fa";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class SamActuatorTest : public ::testing::Test {
public:
    // Prepare data objects for testing
	void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

		pInBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(SamTestData::MAX_BLOCK_SIZE);

        generateSamFile(SamTestData::testSamFile);

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
        std::remove(SamTestData::testSamFile.c_str());
        std::remove("./test/test.sam");
        std::remove("../test/test.sam");
        std::remove("test_reference.fa");
        std::remove("./test/test_reference.fa");
        std::remove("../test/test_reference.fa");
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
        file << "read1_11.1\t0\tchr1\t1\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read3_11.2\t0\tchr1\t153\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCCGCAAG\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!T!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:3\tMD:Z:37C37C0\tAS:i:72\tXS:i:0\n";
        file << "read4_11.3\t0\tchr2\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGACAATCGCCGGAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!F!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read5_11.4\t0\tchr2\t77\t60\t76M\t*\t0\t0\tATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAAATTGCAACGCAATCGATCGA\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:75A0\tAS:i:75\tXS:i:0\n";
        file << "read6_11.5\t16\tchr3\t1\t60\t76M\t*\t0\t0\tCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCGTTAGGCCACAATCGGCGGCCGAATC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read7_11.6\t0\tchr3\t77\t60\t76M\t*\t0\t0\tATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATACGCAATCGGCGGCCGAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:1\tMD:Z:74A0\tAS:i:75\tXS:i:0\n";
        file << "read8_11.7\t16\tchr4\t1\t60\t76M\t*\t0\t0\tTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGGCAACAATCGGCGGCC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read9_11.8\t0\tchrX\t1\t60\t76M\t*\t0\t0\tGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGACGCAATCGCCGGCCTAAT\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";
        file << "read10_11.9\t16\tchrY\t1\t60\t76M\t*\t0\t0\tAATTAAATTTAAATTTCCGGAAATTAAATTTAAATTTCCGGAAATTAAATTTAAACGCAATCGCCGGCCTTGTTC\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tNM:i:2\tMD:Z:74A0\tAS:i:74\tXS:i:0\n";

        file.close();
    }

    void generateFasta(const std::string& filename) {
        std::ifstream srcFile("test_data/reference.fa");
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
                file << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
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

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
    
    // 映射数据，用于测试
    std::map<uint32_t, uint16_t> mappedFlag;
    std::map<uint32_t, uint64_t> mappedPos;
};

TEST_F(SamActuatorTest, testPreAnalysis) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    EXPECT_EQ(actuator.headEndLine, 8);  // 头部有8行 (@HD, 6个@SQ, 1个@PG, 可能还有其他头部行)
    EXPECT_EQ(actuator.contentPos.size(), 10);  // 有10条比对记录
}

TEST_F(SamActuatorTest, testCompressQuality) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressQuality method with field index 10 (QUAL field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressQuality(10, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressChrName) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressChrName method with field index 2 (RNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressChrName(2, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressWithoutRef) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressWithoutRef
    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testCompressWithRef) {
    loadSamData(SamTestData::testSamFile);
    Reference refGene = createTestReference();
    SamActuator actuator(pInBlock, pOutBlock, &refGene);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressWithRef
    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testCompressBaseWithoutRef) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressBaseWithoutRef method with field index 9 (SEQ field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressBaseWithoutRef(9, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressBaseWithRef) {
    loadSamData(SamTestData::testSamFile);
    Reference refGene = createTestReference();
    SamActuator actuator(pInBlock, pOutBlock, &refGene);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressBaseWithRef method with field index 9 (SEQ field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressBaseWithRef(9, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressIdFieldSplit) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressIdFieldSplit method with field index 0 (QNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressIdFieldSplit(fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressIdFieldInAll) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressIdFieldInAll method with field index 0 (QNAME field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressIdFieldInAll(fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testCompressRegularField) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Pre-analysis
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compressRegularField method with field index 6 (PNEXT field)
    uint32_t fieldSrcLen = 0;
    Json::Value fieldMeta;
    result = actuator.compressRegularField(6, fieldSrcLen, fieldMeta);
    EXPECT_GT(result, 0);
}

TEST_F(SamActuatorTest, testNotifyFlag) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    // Test getNotifyFlag
    bool hasData = actuator.getNotifyFlag();
    EXPECT_TRUE(hasData) << "Should have data after loading SAM file";
}

TEST_F(SamActuatorTest, testSetReference) {
    loadSamData(SamTestData::testSamFile);
    SamActuator actuator(pInBlock, pOutBlock);
    
    Reference refGene = createTestReference();
    
    // Test setReference
    actuator.setReference(&refGene);
    
    // Verify reference is set by trying to compress
    int32_t result = actuator.preAnalysis();
    EXPECT_EQ(result, 0);
    
    result = actuator.compress();
    EXPECT_EQ(result, 0);
}

TEST_F(SamActuatorTest, testDecompress) {
    loadSamData(SamTestData::testSamFile);
    
    // 打印压缩前的原始内容
    std::cout << "=== 压缩前的原始内容 ===" << std::endl;
    std::cout << "原始数据长度: " << pInBlock->getDataLen() << std::endl;
    std::cout << "原始内容:" << std::endl;
    std::cout.write(reinterpret_cast<const char*>(pInBlock->getBuffer()), pInBlock->getDataLen());
    std::cout << std::endl << std::endl;
    
    // 创建一个用于压缩的SamActuator对象
    SamActuator compressor(pInBlock, pOutBlock);
    
    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    pInBlock->reset();
    
    // 将压缩后的输出Block内容拷贝到新的input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());
    
    // 重置输出Block
    pOutBlock->reset();
    
    // 创建用于解压缩的SamActuator对象
    SamActuator decompressor(pInBlock, pOutBlock);
    
    // 解压缩不需要preAnalysis，直接调用decompress
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);
    
    // 打印解压后的内容
    std::cout << "=== 解压后的内容 ===" << std::endl;
    std::cout << "解压后数据长度: " << pOutBlock->getDataLen() << std::endl;
    std::cout << "解压后内容:" << std::endl;
    std::cout.write(reinterpret_cast<const char*>(pOutBlock->getBuffer()), pOutBlock->getDataLen());
    std::cout << std::endl << std::endl;
    
    // 基本检查：确保解压缩产生了数据
    EXPECT_GT(pOutBlock->getDataLen(), 0);
}

TEST_F(SamActuatorTest, testDecompressWithRef) {
    loadSamData(SamTestData::testSamFile);
    
    // 打印压缩前的原始内容
    std::cout << "=== 压缩前的原始内容 ===" << std::endl;
    std::cout << "原始数据长度: " << pInBlock->getDataLen() << std::endl;
    std::cout << "原始内容:" << std::endl;
    std::cout.write(reinterpret_cast<const char*>(pInBlock->getBuffer()), pInBlock->getDataLen());
    std::cout << std::endl << std::endl;

    Reference reference = createTestReference();
    
    // 创建一个用于压缩的SamActuator对象
    SamActuator compressor(pInBlock, pOutBlock, &reference);
    
    // Pre-analysis for compression
    int32_t result = compressor.preAnalysis();
    EXPECT_EQ(result, 0);
    
    // Test compress
    result = compressor.compress();
    EXPECT_EQ(result, 0);

    pInBlock->reset();
    
    // 将压缩后的输出Block内容拷贝到新的input Block
    memcpy(pInBlock->getBuffer(), pOutBlock->getBuffer(), pOutBlock->getDataLen() + pOutBlock->getMetaLen());
    pInBlock->setDataLen(pOutBlock->getDataLen());
    pInBlock->setMetaLen(pOutBlock->getMetaLen());
    pInBlock->setBlockType(pOutBlock->getBlockType());
    
    // 重置输出Block
    pOutBlock->reset();
    
    // 创建用于解压缩的SamActuator对象
    SamActuator decompressor(pInBlock, pOutBlock, &reference);
    
    // 解压缩不需要preAnalysis，直接调用decompress
    result = decompressor.decompress();
    EXPECT_EQ(result, 0);
    
    // 打印解压后的内容
    std::cout << "=== 解压后的内容 ===" << std::endl;
    std::cout << "解压后数据长度: " << pOutBlock->getDataLen() << std::endl;
    std::cout << "解压后内容:" << std::endl;
    std::cout.write(reinterpret_cast<const char*>(pOutBlock->getBuffer()), pOutBlock->getDataLen());
    std::cout << std::endl << std::endl;
    
    // 基本检查：确保解压缩产生了数据
    EXPECT_GT(pOutBlock->getDataLen(), 0);
}