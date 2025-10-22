#include "gtest/gtest.h"
#include "utils/path_util.h"

TEST(PathUtil, GetFilePath0) {
    EXPECT_EQ(PathUtil::getFilePath("./io_writer.txt"), "/home/huangmei/workspace/vscode/pbgz_v2/build/test/");
}

TEST(PathUtil, GetFileName0) {
    EXPECT_EQ(PathUtil::getFileName("io_writer.txt"), "io_writer.txt");
}

TEST(PathUtil, GetAbsPath0) {
    EXPECT_EQ(PathUtil::getAbsPath("./io_writer.txt"), "/home/huangmei/workspace/vscode/pbgz_v2/build/test/io_writer.txt");
}

TEST(PathUtil, GetAbsPath1) {
    EXPECT_EQ(PathUtil::getAbsPath("/home/huangmei/workspace/vscode/pbgz_v2/build/test/io_writer.txt"), 
              "/home/huangmei/workspace/vscode/pbgz_v2/build/test/io_writer.txt");
}

