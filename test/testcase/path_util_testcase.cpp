/*
 * path_util_testcase.cpp - Test cases for path utility functionality
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

#include "gtest/gtest.h"
#include "utils/path_util.h"
#include <unistd.h>
#include <limits.h>

// Get the directory where the current binary is located
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

// TEST(PathUtil, GetAbsPath0) {
//     std::string binaryDir = getBinaryDir();
//     std::string expectedPath = binaryDir + "/io_writer.txt";
//     EXPECT_EQ(PathUtil::getAbsPath("./io_writer.txt", expectedPath), expectedPath);
// }

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

// TEST(PathUtil, GetFileNameFromGz) {
//     EXPECT_EQ(PathUtil::getFileNameFromGz("../../data/GCA_000001405.29_GRCh38.p14_genomic.fna.gz"), "GCA_000001405.29_GRCh38.p14_genomic.fna");
// }