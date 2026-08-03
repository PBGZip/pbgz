

#include "pbgz_index.h"
#include "pbgz_testcase_util.h"
#include "sam_info.h"
#include <fstream>
#include <filesystem>


#define protected public
#define private public
#include "index_engine.h"
#undef protected
#undef private


// MockSortEngine - Shadow SortEngine's queues with MockBlockingQueue
class MockIndexEngine : public IndexEngine {
public:
    MockIndexEngine(const PbgzParameter& para, uint32_t blockSz = 8 << 20)
        : IndexEngine(para) {
        freeOutputPool = std::make_unique<MockBlockingQueue>(blockSz);
        outputDataPool = std::make_unique<MockBlockingQueue>(blockSz);
    }

    int32_t init(const std::string& outputFile) {
        /*
         * IndexEngine::startEnginePostProc 一进门就要求 ioReader 是个 FileReader
         * ——索引只对可 seek 的输入才有意义。这个门是后加的，而本 mock 从来只设
         * ioWriter，于是 postProc 在空指针检查处直接 return 0：
         * 断言 result == 0 照样成立，却什么都没写，检查输出的用例才会失败。
         * 这里补一个真实存在的输入文件，让用例走到它真正要验证的那段逻辑。
         */
        inputFileName = outputFile + ".in";
        { std::ofstream seed(inputFileName); seed << "seed"; }
        ioReader = MemoryUtil::safeNewClass<FileReader>(inputFileName);
        ioReader->openIO();

        ioWriter = MemoryUtil::safeNewClass<FileWriter>(outputFile);
        ioWriter->openIO();
        ((MockBlockingQueue*)freeOutputPool.get())->setIOWriter(ioWriter);
        ((MockBlockingQueue*)outputDataPool.get())->setIOWriter(ioWriter);
        return 0;
    }

    void finish() {
        ioWriter->closeIO();
    }

    virtual ~MockIndexEngine(){
        MemoryUtil::safeDeleteClass(ioWriter);
        MemoryUtil::safeDeleteClass(ioReader);
        if (!inputFileName.empty()) {
            std::filesystem::remove(inputFileName);
        }
    }

private:
    std::string inputFileName;
};

class IndexEngineTest : public ::testing::Test {
public:
    void SetUp() override {
        tempOutputFile = "/tmp/test_index_output.txt";
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().addChromosomeInfo("chr1", 1000000);
        SamInfo::getInstance().addChromosomeInfo("chr2", 1000000);
        SamInfo::getInstance().calculateChromosomePositions();
    }

    void TearDown() override {
        SamIndex::getInstance().clear();
        SamInfo::getInstance().clearChromosomeInfo();
        std::filesystem::remove(tempOutputFile);
    }

    std::string tempOutputFile;
};

TEST_F(IndexEngineTest, StartEnginePostProc_EmptyIndex_Succeeds) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 1;

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();
}

TEST_F(IndexEngineTest, StartEnginePostProc_WithSingleBlock_Succeeds) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 1;

    // Add index data for block 0
    SamIndex::getInstance().addSamIndex(1, 1000, 5, 0);
    SamIndex::getInstance().addSamIndex(2, 2000, 3, 0);

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();

    // Verify output file contains expected data
    std::ifstream file(tempOutputFile);
    ASSERT_TRUE(file.is_open());

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    EXPECT_EQ(2, lines.size());
    // Lines should be sorted by blockId first (all block 0), then chrIndex, then position
    EXPECT_EQ("1\t1000\t5\t0\t0", lines[0]);
    EXPECT_EQ("2\t2000\t3\t0\t0", lines[1]);
}

TEST_F(IndexEngineTest, StartEnginePostProc_WithMultipleBlocks_Succeeds) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 3;

    // Add index data for multiple blocks
    SamIndex::getInstance().addSamIndex(1, 1000, 5, 0);
    SamIndex::getInstance().addSamIndex(2, 2000, 3, 0);
    SamIndex::getInstance().addSamIndex(1, 5000, 8, 1);
    SamIndex::getInstance().addSamIndex(2, 6000, 2, 1);
    SamIndex::getInstance().addSamIndex(3, 7000, 15, 2);

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();

    // Verify output file contains expected data in correct order
    std::ifstream file(tempOutputFile);
    ASSERT_TRUE(file.is_open());

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    EXPECT_EQ(5, lines.size());
    // Sorted: blockId, then chrIndex, then referenceMapPos
    EXPECT_EQ("1\t1000\t5\t0\t0", lines[0]);
    EXPECT_EQ("2\t2000\t3\t0\t0", lines[1]);
    EXPECT_EQ("1\t5000\t8\t1\t0", lines[2]);
    EXPECT_EQ("2\t6000\t2\t1\t0", lines[3]);
    EXPECT_EQ("3\t7000\t15\t2\t0", lines[4]);
}

TEST_F(IndexEngineTest, StartEnginePostProc_SameBlockDifferentChromosomes_SortingWorks) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 1;

    // Add index data for block 0 with different chromosomes (out of order)
    SamIndex::getInstance().addSamIndex(3, 3000, 4, 0);
    SamIndex::getInstance().addSamIndex(1, 1000, 5, 0);
    SamIndex::getInstance().addSamIndex(2, 2000, 3, 0);

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();

    // Verify output is sorted
    std::ifstream file(tempOutputFile);
    ASSERT_TRUE(file.is_open());

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    EXPECT_EQ(3, lines.size());
    EXPECT_EQ("1\t1000\t5\t0\t0", lines[0]);
    EXPECT_EQ("2\t2000\t3\t0\t0", lines[1]);
    EXPECT_EQ("3\t3000\t4\t0\t0", lines[2]);
}

TEST_F(IndexEngineTest, StartEnginePostProc_SameBlockSameChromosome_SortingWorks) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 1;

    // Add index data for block 0 with same chromosome, different positions (out of order)
    SamIndex::getInstance().addSamIndex(1, 3000, 4, 0);
    SamIndex::getInstance().addSamIndex(1, 1000, 5, 0);
    SamIndex::getInstance().addSamIndex(1, 2000, 3, 0);

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();

    // Verify output is sorted by position
    std::ifstream file(tempOutputFile);
    ASSERT_TRUE(file.is_open());

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }
    file.close();

    EXPECT_EQ(3, lines.size());
    EXPECT_EQ("1\t1000\t5\t0\t0", lines[0]);
    EXPECT_EQ("1\t2000\t3\t0\t0", lines[1]);
    EXPECT_EQ("1\t3000\t4\t0\t0", lines[2]);
}

TEST_F(IndexEngineTest, StartEnginePostProc_OutputWrittenToFile_Success) {
    PbgzParameter para;
    para.compressLevel = 1;
    para.threadNum = 1;

    MockIndexEngine engine(para);
    engine.init(tempOutputFile);
    engine.blockCount = 1;

    SamIndex::getInstance().addSamIndex(1, 1000, 5, 0);
    SamIndex::getInstance().addSamIndex(2, 2000, 3, 0);

    int32_t result = engine.startEnginePostProc();
    EXPECT_EQ(0, result);

    engine.finish();

    // Verify output file was written with correct content
    std::ifstream file(tempOutputFile);
    ASSERT_TRUE(file.is_open());

    std::string content;
    std::string line;
    while (std::getline(file, line)) {
        if (!content.empty()) {
            content += "\n";
        }
        content += line;
    }
    file.close();

    // Content should contain the expected format with tabs and newlines
    EXPECT_NE(std::string::npos, content.find("1\t1000\t5\t0"));
    EXPECT_NE(std::string::npos, content.find("2\t2000\t3\t0"));
}

