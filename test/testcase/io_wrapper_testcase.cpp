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

// FileReader readLine测试 - 基本功能
TEST(FileReaderTest, ReadLineBasic) {
    std::string fileName = "readline_basic.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    writer.writeIO("Line 1\nLine 2\nLine 3\n", 21);
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // FileReader返回的是行内容长度（不包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 1" = 6个字符
    EXPECT_EQ(line, "Line 1");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 2" = 6个字符
    EXPECT_EQ(line, "Line 2");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 3" = 6个字符
    EXPECT_EQ(line, "Line 3");
    
    // 读取EOF
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine测试 - 超长行（超过8KB缓冲区）
TEST(FileReaderTest, ReadLineLongLine) {
    std::string fileName = "readline_long.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // 创建一个超过8KB的长行
    std::string longLine(10000, 'X');  // 10KB长行
    longLine += "\n";
    writer.writeIO(longLine.c_str(), longLine.length());
    
    // 添加一个短行
    writer.writeIO("Short\n", 6);
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // FileReader返回的是行内容长度（不包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10000);  // 10000个'X'（不包括换行符）
    EXPECT_EQ(line.length(), 10000);
    EXPECT_EQ(line[0], 'X');
    EXPECT_EQ(line[9999], 'X');
    
    // 读取短行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 5);  // "Short" = 5个字符
    EXPECT_EQ(line, "Short");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine测试 - 跨缓冲区场景
TEST(FileReaderTest, ReadLineCrossBuffer) {
    std::string fileName = "readline_cross_buffer.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // 创建多个行，总长度超过8KB，确保跨缓冲区读取
    std::string content;
    for (int i = 0; i < 100; i++) {
        content += "Line " + std::to_string(i) + ": " + std::string(80, 'A' + (i % 26)) + "\n";
    }
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t totalBytes = 0;
    int lineCount = 0;
    
    // 读取所有行
    while (true) {
        size_t bytesRead = reader.readLine(line);
        if (bytesRead == 0) break;
        
        totalBytes += bytesRead;
        lineCount++;
        
        // 验证行内容格式
        EXPECT_GT(line.length(), 0);
        EXPECT_TRUE(line.find("Line ") == 0);
        
        // 验证行号
        std::string expectedPrefix = "Line " + std::to_string(lineCount - 1) + ": ";
        EXPECT_TRUE(line.find(expectedPrefix) == 0);
    }
    
    // 验证读取了所有行
    EXPECT_EQ(lineCount, 100);
    EXPECT_GT(totalBytes, 8000);  // 确保总字节数超过8KB
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine测试 - Windows换行符
TEST(FileReaderTest, ReadLineWindowsCRLF) {
    std::string fileName = "readline_crlf.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // 写入Windows风格的换行符
    std::string content = "Line 1\r\nLine 2\r\nLine 3\r\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // 读取第一行（应该去掉\r）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 1");
    
    // 读取第二行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 2");
    
    // 读取第三行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 3");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine测试 - 空行和单字符行
TEST(FileReaderTest, ReadLineEdgeCases) {
    std::string fileName = "readline_edge.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // 写入各种边界情况
    std::string content = "\nA\n\nBC\n\n\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // 读取空行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    // 读取单字符行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 1);
    EXPECT_EQ(line, "A");
    
    // 读取空行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    // 读取两字符行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 2);
    EXPECT_EQ(line, "BC");
    
    // 读取连续空行
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// GzFileReader readLine测试 - 基本功能
TEST(GzFileReaderTest, ReadLineBasic) {
    std::string fileName = "gz_readline_basic.txt";
    std::string gzFileName = "gz_readline_basic.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    writer.writeIO("Line 1\nLine 2\nLine 3\n", 21);
    writer.closeIO();

    system("gzip gz_readline_basic.txt -f");

    GzFileReader reader(gzFileName, 1);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // 读取第一行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 1\n" = 7个字符
    EXPECT_EQ(line, "Line 1");
    
    // 读取第二行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 2\n" = 7个字符
    EXPECT_EQ(line, "Line 2");
    
    // 读取第三行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 3\n" = 7个字符
    EXPECT_EQ(line, "Line 3");
    
    // 读取EOF
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// GzFileReader readLine测试 - 超长行和跨缓冲区
TEST(GzFileReaderTest, ReadLineLongAndCrossBuffer) {
    std::string fileName = "gz_readline_long.txt";
    std::string gzFileName = "gz_readline_long.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // 创建超长行和多个行，测试跨缓冲区
    std::string content;
    content += std::string(12000, 'X') + "\n";  // 12KB长行
    for (int i = 0; i < 50; i++) {
        content += "GzLine " + std::to_string(i) + ": " + std::string(100, 'A' + (i % 26)) + "\n";
    }
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    system("gzip gz_readline_long.txt -f");

    GzFileReader reader(gzFileName, 1);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    int lineCount = 0;
    
    // 读取超长行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 12001);  // 12000个'X' + 1个换行符
    EXPECT_EQ(line.length(), 12000);
    lineCount++;
    
    // 读取后续所有行
    while (true) {
        bytesRead = reader.readLine(line);
        if (bytesRead == 0) break;
        
        lineCount++;
        // 验证行内容格式
        EXPECT_TRUE(line.find("GzLine ") == 0);
    }
    
    // 验证读取了所有行
    EXPECT_EQ(lineCount, 51);  // 1个长行 + 50个普通行
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// FastGzFileReader readLine测试
TEST(FastGzFileReaderTest, ReadLineBasic) {
    std::string fileName = "fastgz_readline_basic.txt";
    std::string gzFileName = "fastgz_readline_basic.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    writer.writeIO("FastLine1\nFastLine2\nFastLine3\n", 30);
    writer.closeIO();

    system("gzip fastgz_readline_basic.txt -f");

    FastGzFileReader reader(gzFileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // 读取所有行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine1\n" = 10个字符
    EXPECT_EQ(line, "FastLine1");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine2\n" = 10个字符
    EXPECT_EQ(line, "FastLine2");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine3\n" = 10个字符
    EXPECT_EQ(line, "FastLine3");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// PipeReader readLine测试
TEST(PipeReaderTest, ReadLineBasic) {
    std::string fileName = "pipe_readline.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string content = "PipeLine1\nPipeLine2\n";
    content += std::string(9000, 'Y') + "\n";  // 9KB长行，测试跨缓冲区
    content += "PipeLine3\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    // 通过管道发送数据
    std::string cmd = "cat " + fileName;
    FILE* pipe = popen(cmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    int saved_stdin = dup(STDIN_FILENO);
    dup2(fileno(pipe), STDIN_FILENO);
    
    PipeReader reader;
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // 读取第一行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine1\n" = 10个字符
    EXPECT_EQ(line, "PipeLine1");
    
    // 读取第二行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine2\n" = 10个字符
    EXPECT_EQ(line, "PipeLine2");
    
    // 读取长行（跨缓冲区，包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 9001);  // 9000个'Y' + 1个换行符
    EXPECT_EQ(line.length(), 9000);
    EXPECT_EQ(line[0], 'Y');
    EXPECT_EQ(line[8999], 'Y');
    
    // 读取第三行（包括换行符）
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine3\n" = 10个字符
    EXPECT_EQ(line, "PipeLine3");
    
    reader.closeIO();
    
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    pclose(pipe);
    
    remove(fileName.c_str());
}
