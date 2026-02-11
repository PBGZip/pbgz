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
    
    // Create GzFileWriter and write data
    GzFileWriter gzWriter(gzFileName, 1);
    EXPECT_EQ(gzWriter.openIO(), 0);
    std::string testTxt = "this is a test file for gzip writer";
    size_t writeLen = gzWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    gzWriter.closeIO();

    // Verify gzip file was created successfully
    FILE* gzFile = fopen(gzFileName.c_str(), "rb");
    EXPECT_NE(gzFile, nullptr);
    if (gzFile) {
        fclose(gzFile);
    }

    // Use GzFileReader to read and verify data
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
    // Create test data file
    std::string fileName = "pipe_test.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for pipe reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // Send data to PipeReader through pipe
    std::string cmd = "cat " + fileName;
    FILE* pipe = popen(cmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // Save original stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // Redirect pipe output to stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // Use PipeReader to read data
    PipeReader pipeReader;
    EXPECT_EQ(pipeReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = pipeReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    pipeReader.closeIO();
    
    // Restore original stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // Close pipe
    // pclose(pipe);
    
    // Clean up
    remove(fileName.c_str());
}

TEST(PipeWriterTest, PipeWriterWrite) {
    // Create temporary file to capture PipeWriter output
    std::string fileName = "pipe_output.txt";
    (void)remove(fileName.c_str());
    
    // Save original stdout
    int saved_stdout = dup(STDOUT_FILENO);
    
    // Redirect stdout to file
    int fd = open(fileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    EXPECT_GE(fd, 0);
    dup2(fd, STDOUT_FILENO);
    close(fd);
    
    // Use PipeWriter to write data
    PipeWriter pipeWriter;
    EXPECT_EQ(pipeWriter.openIO(), 0);
    
    std::string testTxt = "this is a test file for pipe writer";
    size_t writeLen = pipeWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    
    pipeWriter.closeIO();
    
    // Restore original stdout
    dup2(saved_stdout, STDOUT_FILENO);
    close(saved_stdout);
    
    // Verify file content
    FileReader reader(fileName);
    reader.openIO();
    EXPECT_EQ(memcmp(testTxt.c_str(), reader.getAt(0), testTxt.length()), 0);
    reader.closeIO();
    
    // Clean up
    remove(fileName.c_str());
}

TEST(PileGzReaderTest, PileGzReaderRead) {
    // Create test gzip file
    std::string fileName = "pile_gz_test.txt";
    std::string gzFileName = "pile_gz_test.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for pile gzip reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // Compress file
    std::string cmd = "gzip " + fileName + " -f";
    system(cmd.c_str());

    // Send gzip data to PileGzReader through pipe
    std::string pipeCmd = "cat " + gzFileName;
    FILE* pipe = popen(pipeCmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // Save original stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // Redirect pipe output to stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // Use PileGzReader to read data
    GzPipeReader pileGzReader(1);
    EXPECT_EQ(pileGzReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = pileGzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    pileGzReader.closeIO();
    
    // Restore original stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // Close pipe
    pclose(pipe);
    
    // Clean up
    remove(gzFileName.c_str());
}

TEST(PileGzWriterTest, PileGzWriterWrite) {
    std::string outputFileName = "pile_gz_output.txt.gz";
    (void)remove(outputFileName.c_str());
    
    // Save original stdout
    int saved_stdout = dup(STDOUT_FILENO);
    
    // Redirect stdout to file
    int fd = open(outputFileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    EXPECT_GE(fd, 0);
    dup2(fd, STDOUT_FILENO);
    close(fd);
    
    // Use PileGzWriter to write data
    GzPipeWriter pileGzWriter(1);
    EXPECT_EQ(pileGzWriter.openIO(), 0);
    
    std::string testTxt = "this is a test file for pile gzip writer";
    size_t writeLen = pileGzWriter.writeIO(testTxt.c_str(), testTxt.length());
    EXPECT_EQ(writeLen, testTxt.length());
    
    pileGzWriter.closeIO();
    
    // Restore original stdout
    dup2(saved_stdout, STDOUT_FILENO);
    close(saved_stdout);
    
    // Verify gzip file was created successfully
    FILE* gzFile = fopen(outputFileName.c_str(), "rb");
    EXPECT_NE(gzFile, nullptr);
    if (gzFile) {
        fclose(gzFile);
    }

    // Use GzFileReader to read and verify data
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
    // Create test gzip file
    std::string fileName = "fast_pile_gz_test.txt";
    std::string gzFileName = "fast_pile_gz_test.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string testTxt = "this is a test file for fast pile gzip reader";
    writer.writeIO(testTxt.c_str(), testTxt.length());
    writer.closeIO();

    // Compress file
    std::string cmd = "gzip " + fileName + " -f";
    system(cmd.c_str());

    // Send gzip data to FastPileGzReader through pipe
    std::string pipeCmd = "cat " + gzFileName;
    FILE* pipe = popen(pipeCmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    // Save original stdin
    int saved_stdin = dup(STDIN_FILENO);
    
    // Redirect pipe output to stdin
    dup2(fileno(pipe), STDIN_FILENO);
    
    // Use FastPileGzReader to read data
    FastGzPipeReader fastPileGzReader;
    EXPECT_EQ(fastPileGzReader.openIO(), 0);
    
    char buffer[1024];
    size_t readLen = fastPileGzReader.readIO(buffer, 1024);
    EXPECT_EQ(readLen, testTxt.size());
    EXPECT_EQ(memcmp(buffer, testTxt.c_str(), testTxt.size()), 0);
    
    fastPileGzReader.closeIO();
    
    // Restore original stdin
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    
    // Close pipe
    pclose(pipe);
    
    // Clean up
    remove(gzFileName.c_str());
}

// FileReader readLine test - basic functionality
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
    
    // FileReader returns line content length (excluding newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 1" = 6 characters
    EXPECT_EQ(line, "Line 1");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 2" = 6 characters
    EXPECT_EQ(line, "Line 2");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);  // "Line 3" = 6 characters
    EXPECT_EQ(line, "Line 3");
    
    // Read EOF
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine test - long line (exceeding 8KB buffer)
TEST(FileReaderTest, ReadLineLongLine) {
    std::string fileName = "readline_long.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // Create a long line exceeding 8KB
    std::string longLine(10000, 'X');  // 10KB long line
    longLine += "\n";
    writer.writeIO(longLine.c_str(), longLine.length());
    
    // Add a short line
    writer.writeIO("Short\n", 6);
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // FileReader returns line content length (excluding newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10000);  // 10000 'X' characters (excluding newline)
    EXPECT_EQ(line.length(), 10000);
    EXPECT_EQ(line[0], 'X');
    EXPECT_EQ(line[9999], 'X');
    
    // Read short line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 5);  // "Short" = 5 characters
    EXPECT_EQ(line, "Short");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine test - cross buffer scenario
TEST(FileReaderTest, ReadLineCrossBuffer) {
    std::string fileName = "readline_cross_buffer.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // Create multiple lines with total length exceeding 8KB to ensure cross-buffer reading
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
    
    // Read all lines
    while (true) {
        size_t bytesRead = reader.readLine(line);
        if (bytesRead == 0) break;
        
        totalBytes += bytesRead;
        lineCount++;
        
        // Verify line content format
        EXPECT_GT(line.length(), 0);
        EXPECT_TRUE(line.find("Line ") == 0);
        
        // Verify line number
        std::string expectedPrefix = "Line " + std::to_string(lineCount - 1) + ": ";
        EXPECT_TRUE(line.find(expectedPrefix) == 0);
    }
    
    // Verify all lines were read
    EXPECT_EQ(lineCount, 100);
    EXPECT_GT(totalBytes, 8000);  // Ensure total bytes exceed 8KB
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine test - Windows line endings
TEST(FileReaderTest, ReadLineWindowsCRLF) {
    std::string fileName = "readline_crlf.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // Write Windows-style line endings
    std::string content = "Line 1\r\nLine 2\r\nLine 3\r\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // Read first line (should remove \r)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 1");
    
    // Read second line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 2");
    
    // Read third line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 6);
    EXPECT_EQ(line, "Line 3");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// FileReader readLine test - empty lines and single character lines
TEST(FileReaderTest, ReadLineEdgeCases) {
    std::string fileName = "readline_edge.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // Write various edge cases
    std::string content = "\nA\n\nBC\n\n\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    FileReader reader(fileName);
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // Read empty line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    // Read single character line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 1);
    EXPECT_EQ(line, "A");
    
    // Read empty line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    // Read two character line
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 2);
    EXPECT_EQ(line, "BC");
    
    // Read consecutive empty lines
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    EXPECT_EQ(line, "");
    
    reader.closeIO();
    remove(fileName.c_str());
}

// GzFileReader readLine test - basic functionality
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
    
    // Read first line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 1\n" = 7 characters
    EXPECT_EQ(line, "Line 1");
    
    // Read second line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 2\n" = 7 characters
    EXPECT_EQ(line, "Line 2");
    
    // Read third line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 7);  // "Line 3\n" = 7 characters
    EXPECT_EQ(line, "Line 3");
    
    // Read EOF
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// GzFileReader readLine test - long lines and cross buffer
TEST(GzFileReaderTest, ReadLineLongAndCrossBuffer) {
    std::string fileName = "gz_readline_long.txt";
    std::string gzFileName = "gz_readline_long.txt.gz";
    (void)remove(fileName.c_str());
    (void)remove(gzFileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    
    // Create long lines and multiple lines to test cross buffer
    std::string content;
    content += std::string(12000, 'X') + "\n";  // 12KB long line
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
    
    // Read long line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 12001);  // 12000 'X' characters + 1 newline
    EXPECT_EQ(line.length(), 12000);
    lineCount++;
    
    // Read all subsequent lines
    while (true) {
        bytesRead = reader.readLine(line);
        if (bytesRead == 0) break;
        
        lineCount++;
        // Verify line content format
        EXPECT_TRUE(line.find("GzLine ") == 0);
    }
    
    // Verify all lines were read
    EXPECT_EQ(lineCount, 51);  // 1 long line + 50 normal lines
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// FastGzFileReader readLine test
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
    
    // Read all lines (including newlines)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine1\n" = 10 characters
    EXPECT_EQ(line, "FastLine1");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine2\n" = 10 characters
    EXPECT_EQ(line, "FastLine2");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "FastLine3\n" = 10 characters
    EXPECT_EQ(line, "FastLine3");
    
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 0);
    
    reader.closeIO();
    remove(gzFileName.c_str());
    remove(fileName.c_str());
}

// PipeReader readLine test
TEST(PipeReaderTest, ReadLineBasic) {
    std::string fileName = "pipe_readline.txt";
    (void)remove(fileName.c_str());
    
    FileWriter writer(fileName);
    writer.openIO();
    std::string content = "PipeLine1\nPipeLine2\n";
    content += std::string(9000, 'Y') + "\n";  // 9KB long line, test cross buffer
    content += "PipeLine3\n";
    writer.writeIO(content.c_str(), content.length());
    writer.closeIO();

    // Send data through pipe
    std::string cmd = "cat " + fileName;
    FILE* pipe = popen(cmd.c_str(), "r");
    EXPECT_NE(pipe, nullptr);
    
    int saved_stdin = dup(STDIN_FILENO);
    dup2(fileno(pipe), STDIN_FILENO);
    
    PipeReader reader;
    EXPECT_EQ(reader.openIO(), 0);
    
    std::string line;
    size_t bytesRead;
    
    // Read first line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine1\n" = 10 characters
    EXPECT_EQ(line, "PipeLine1");
    
    // Read second line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine2\n" = 10 characters
    EXPECT_EQ(line, "PipeLine2");
    
    // Read long line (cross buffer, including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 9001);  // 9000 'Y' characters + 1 newline
    EXPECT_EQ(line.length(), 9000);
    EXPECT_EQ(line[0], 'Y');
    EXPECT_EQ(line[8999], 'Y');
    
    // Read third line (including newline)
    bytesRead = reader.readLine(line);
    EXPECT_EQ(bytesRead, 10);  // "PipeLine3\n" = 10 characters
    EXPECT_EQ(line, "PipeLine3");
    
    reader.closeIO();
    
    dup2(saved_stdin, STDIN_FILENO);
    close(saved_stdin);
    pclose(pipe);
    
    remove(fileName.c_str());
}
