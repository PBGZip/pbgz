/*
 * path_util.h - Header file for pbgz project
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

#pragma once

#include <string>

namespace PathUtil {
    /// @brief Check if filename ends with suffix
    /// @param fileName  File name
    /// @param suffix    Suffix
    /// @return bool     true/false
    bool suffixCheck(const std::string& fileName, const std::string& suffix);

    /// @brief Check if file exists
    /// @param fileName  File name
    /// @return bool     true: exists / false: not exists
    bool fileExists(const std::string& fileName);

    /// @brief Check if file is readable
    /// @param fileName 
    /// @return  bool  true/false
    bool fileReadble(const std::string& fileName);

    /// @brief Check if file is writable
    /// @param fileName 
    /// @return  bool  true/false
    bool fileWriteable(const std::string& fileName);

    /// @brief Get file name without directory from full path
    /// @param fullFileName 
    /// @return File name without directory
    std::string getFileName(const std::string& fullFileName);

    /// @brief Get directory path from full file name
    /// @param fullFileName 
    /// @return Directory where file is located
    std::string getFilePath(const std::string& fullFileName);

    /// @brief Get absolute path of file
    /// @param fileName 
    /// @return Absolute path of file
    std::string getAbsPath(const std::string& fileName);

    /// @brief Get absolute path of file relative to base directory
    /// @param fileName 
    /// @param baseDir Base directory path
    /// @return Absolute path of file
    std::string getAbsPath(const std::string& fileName, const std::string& baseDir);

    /// @brief Create directory
    /// @param dir 
    /// @return Absolute path of created directory
    std::string createDir(const std::string& dir);

    /// @brief Check if path is a directory
    /// @param filename 
    /// @return  bool  true/false
    bool isDir(const std::string& filename);

    /// @brief Check if path is a regular file
    /// @param filename 
    /// @return bool  true/false
    bool isFile(const std::string& filename);

    /// @brief Remove file
    /// @param fileName 
    /// @return bool  true/false
    bool removeFile(const std::string& fileName);

    /// @brief Get home directory
    /// @return Absolute path of home directory
    std::string getHomePath();

    /// @brief Get file last modification time
    /// @return 
    int64_t getFileMtime(const std::string& fileName);

    /// @brief Get file size
    /// @param fileName 
    /// @return File size
    int64_t getFileSize(const std::string& fileName);

    bool isGzFile(const std::string& fileName);

    std::string getFileNameFromGz(const std::string& fileName);
}
