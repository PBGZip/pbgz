#include "gtest/gtest.h"
#include "io_wrapper.h"

TEST(IOReaderTest, IOReaderRead) {
    std::string fileName = "io_writer.txt";
    (void)remove(fileName.c_str());
    // First construct a file
    FILE* pFile = fopen(fileName.c_str(), "w+");
    std::string testTxt = "this is a test file";
    fwrite(testTxt.c_str(), testTxt.length(), 1, pFile);
    fflush(pFile);
    fclose(pFile);

    FileReader ioReader(fileName);
    ioReader.openIO();
    char* pReadBuf = (char*)malloc(testTxt.size());
    ioReader.readIO(pReadBuf, testTxt.length());
    ioReader.closeIO();
    EXPECT_EQ(memcmp(testTxt.c_str(), pReadBuf, testTxt.length()), 0);
    
    // Clean up
    free(pReadBuf);
    remove(fileName.c_str());
}

TEST(IOWriterTest, IOWriterWrite) {
    std::string fileName = "io_writer.txt";
    (void)remove(fileName.c_str());
    
    FileWriter ioWriter(fileName);
    ioWriter.openIO();
    std::string testTxt = "this is a test file";
    ioWriter.writeIO(testTxt.c_str(), testTxt.length());
    ioWriter.closeIO();

    FileReader ioReader(fileName);
    ioReader.openIO();
    EXPECT_EQ(memcmp(testTxt.c_str(), ioReader.getAt(0), testTxt.length()), 0);
    ioReader.closeIO();
    
    // Clean up
    remove(fileName.c_str());
}

TEST(GzFileReaderTest, GzFileReaderRead) {
    std::string fileName = "gz_reader.txt";
    std::string gzFileName = "gz_reader.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter ioWriter(fileName);
    ioWriter.openIO();
    std::string testTxt = "this is a test file for gzip";
    ioWriter.writeIO(testTxt.c_str(), testTxt.length());
    ioWriter.closeIO();

    system("gzip gz_reader.txt -f");

    GzFileReader gzReader(gzFileName, 1);
    EXPECT_EQ(gzReader.openIO(), 0);

    char buffer[1024];
    uint32_t readLen = gzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    // Clean up
    gzReader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

TEST(GzFileReaderTest, FastGzFileReaderRead) {
    std::string fileName = "gz_reader.txt";
    std::string gzFileName = "gz_reader.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter ioWriter(fileName);
    ioWriter.openIO();
    std::string testTxt = "this is a test file for Fast Gzip";
    ioWriter.writeIO(testTxt.c_str(), testTxt.length());
    ioWriter.closeIO();

    system("gzip gz_reader.txt -f");

    FastGzFileReader gzReader(gzFileName);
    EXPECT_EQ(gzReader.openIO(), 0);
    
    char buffer[1024];
    uint32_t readLen = gzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    // Clean up
    gzReader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

TEST(GzFileWireterTest, GzFileWriterWrite) {
    std::string gzFileName = "gz_writer.txt.gz";
    (void)remove(gzFileName.c_str());
    
    // 创建GzFileWriter并写入数据
    GzFileWriter gzWriter(gzFileName, 1);
    EXPECT_EQ(gzWriter.openIO(), 0);
    std::string testTxt = "this is a test file for gzip writer";
    size_t writeLen = gzWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    gzWriter.closeIO();

    // 验证gzip文件是否创建成功
    FILE* gzFile = fopen(gzFileName.c_str(), "rb");
    EXPECT_NE(gzFile, nullptr);
    if (gzFile) {
        fclose(gzFile);
    }

    // 使用GzFileReader读取并验证数据
    GzFileReader gzReader(gzFileName, 1);
    EXPECT_EQ(gzReader.openIO(), 0);
    
    char buffer[1024];
    uint32_t readLen = gzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    gzReader.closeIO();
    
    // Clean up
    remove(gzFileName.c_str());
}

TEST(PipeReaderTest, PipeReaderRead) {
    // 创建测试数据文件
    std::string fileName = "pipe_test.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for pipe reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // 通过管道将数据发送给PipeReader
    std::string cmd = "cat " + fileName;
    FILE* pipe = popen(cmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // 保存原始的stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // 将管道输出重定向到stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // 使用PipeReader读取数据
    PipeReader pipeReader;
    EXPECT_EQ(pipeReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = pipeReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    pipeReader.closeIO();
    
    // 恢复原始的stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // 关闭管道
    // pclose(pipe);
    
    // Clean up
    remove(fileName.c_str());
}

TEST(PipeWriterTest, PipeWriterWrite) {
    // 创建临时文件来捕获PipeWriter的输出
    std::string fileName = "pipe_output.txt";
    (void)remove(fileName.c_str());
    
    // 保存原始的stdout
    int saved_stdout = dup(STDOUT_FILENO);
    
    // 重定向stdout到文件
    int fd = open(fileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    EXPECT_GE(fd, 0);
    dup2(fd, STDOUT_FILENO);
    close(fd);
    
    // 使用PipeWriter写入数据
    PipeWriter pipeWriter;
    EXPECT_EQ(pipeWriter.openIO(), 0);
    
    std::string testTxt = "this is a test file for pipe writer";
    size_t writeLen = pipeWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    
    pipeWriter.closeIO();
    
    // 恢复原始的stdout
    dup2(saved_stdout, STDOUT_FILENO);
    close(saved_stdout);
    
    // 验证文件内容
    FileReader reader(fileName);
    reader.openIO();
    EXPECT_EQ(memcmp(testTxt.c_str(), reader.getAt(0), testTxt.length()), 0);
    reader.closeIO();
    
    // Clean up
    remove(fileName.c_str());
}

TEST(PileGzReaderTest, PileGzReaderRead) {
    // 创建测试gzip文件
    std::string fileName = "pile_gz_test.txt";
    std::string gzFileName = "pile_gz_test.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for pile gzip reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // 压缩文件
    std::string cmd = "gzip " + fileName + " -f";
    system(cmd.c_str());

    // 通过管道将gzip数据发送给PileGzReader
    std::string pipeCmd = "cat " + gzFileName;
    FILE* pipe = popen(pipeCmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // 保存原始的stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // 将管道输出重定向到stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // 使用PileGzReader读取数据
    GzPipeReader pileGzReader(1);
    EXPECT_EQ(pileGzReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = pileGzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    pileGzReader.closeIO();
    
    // 恢复原始的stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // 关闭管道
    pclose(pipe);
    
    // Clean up
    remove(gzFileName.c_str());
}

TEST(PileGzWriterTest, PileGzWriterWrite) {
    std::string outputFileName = "pile_gz_output.txt.gz";
    (void)remove(outputFileName.c_str());
    
    // 保存原始的stdout
    int saved_stdout = dup(STDOUT_FILENO);
    
    // 重定向stdout到文件
    int fd = open(outputFileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    EXPECT_GE(fd, 0);
    dup2(fd, STDOUT_FILENO);
    close(fd);
    
    // 使用PileGzWriter写入数据
    GzPipeWriter pileGzWriter(1);
    EXPECT_EQ(pileGzWriter.openIO(), 0);
    
    std::string testTxt = "this is a test file for pile gzip writer";
    size_t writeLen = pileGzWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    
    pileGzWriter.closeIO();
    
    // 恢复原始的stdout
    dup2(saved_stdout, STDOUT_FILENO);
    close(saved_stdout);
    
    // 验证gzip文件是否创建成功
    FILE* gzFile = fopen(outputFileName.c_str(), "rb");
    EXPECT_NE(gzFile, nullptr);
    if (gzFile) {
        fclose(gzFile);
    }

    // 使用GzFileReader读取并验证数据
    GzFileReader gzReader(outputFileName, 1);
    EXPECT_EQ(gzReader.openIO(), 0);
    
    char buffer[1024];
    uint32_t readLen = gzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    gzReader.closeIO();
    
    // Clean up
    remove(outputFileName.c_str());
}

TEST(FastPileGzReaderTest, FastPileGzReaderRead) {
    // 创建测试gzip文件
    std::string fileName = "fast_pile_gz_test.txt";
    std::string gzFileName = "fast_pile_gz_test.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for fast pile gzip reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // 压缩文件
    std::string cmd = "gzip " + fileName + " -f";
    system(cmd.c_str());

    // 通过管道将gzip数据发送给FastPileGzReader
    std::string pipeCmd = "cat " + gzFileName;
    FILE* pipe = popen(pipeCmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // 保存原始的stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // 将管道输出重定向到stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // 使用FastPileGzReader读取数据
    FastGzPipeReader fastPileGzReader;
    EXPECT_EQ(fastPileGzReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = fastPileGzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    fastPileGzReader.closeIO();
    
    // 恢复原始的stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // 关闭管道
    pclose(pipe);
    
    // Clean up
    remove(gzFileName.c_str());
}
