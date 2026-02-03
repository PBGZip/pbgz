#include "gtest/gtest.h"
#include "utils/path_util.h"
#include <unistd.h>
#include <limits.h>

// 获取当前二进制所在目录
std::string getBinaryDir() {
    char path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (len != -1) {
        path[len] = '\0';
        std::string binaryPath(path);
        return binaryPath.substr(0, binaryPath.find_last_of('/'));
    }
    return "";
}

TEST(PathUtil, GetFileName0) {
    EXPECT_EQ(PathUtil::getFileName("io_writer.txt"), "io_writer.txt");
}

TEST(PathUtil, GetAbsPath0) {
    std::string binaryDir = getBinaryDir();
    std::string expectedPath = binaryDir + "/io_writer.txt";
    EXPECT_EQ(PathUtil::getAbsPath("./io_writer.txt"), expectedPath);
}

TEST(PathUtil, GetAbsPath1) {
    std::string binaryDir = getBinaryDir();
    std::string expectedPath = binaryDir + "/io_writer.txt";
    EXPECT_EQ(PathUtil::getAbsPath(expectedPath), expectedPath);
}

TEST(PathUtil, GetAbsPathWithBaseDir) {
    std::string binaryDir = getBinaryDir();
    std::string expectedPath = binaryDir + "/io_writer.txt";
    EXPECT_EQ(PathUtil::getAbsPath("./io_writer.txt", binaryDir), expectedPath);
}
