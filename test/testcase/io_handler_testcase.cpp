#include "gtest/gtest.h"
#include "io_wrapper.h"

TEST(IOReaderTest, IOReaderRead) {
    // 先构造一个文件
    std::string fileName = "io_reader.txt";
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
}