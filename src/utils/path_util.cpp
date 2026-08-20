/*
 * path_util.cpp - Source file for pbgz project
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

#include <cstring>
#include <filesystem>
#include <system_error>
#include <fstream>
#include <pwd.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "path_util.h"
#include "log/logger.h"
#include <vector>


namespace PathUtil {
    bool suffixCheck(const std::string& fileName,const std::string& suffix) {
        if (suffix.empty())
            return true;
        if (suffix.length() > fileName.length())
            return false;
        return (!strncmp(fileName.c_str() + fileName.length() - suffix.length(), suffix.c_str(), suffix.length()));
    }

    bool fileExists(const std::string& fileName) {
        return std::filesystem::exists(fileName);
    }

    bool fileReadble(const std::string& fileName) {
        std::error_code ec;
        auto perms = std::filesystem::status(fileName, ec).permissions();
        return (perms & std::filesystem::perms::owner_read) != std::filesystem::perms::none && !ec;
    }

    bool fileWriteable(const std::string& fileName) {
        std::error_code ec;
        auto perms = std::filesystem::status(fileName, ec).permissions();
        return (perms & std::filesystem::perms::owner_write) != std::filesystem::perms::none && !ec;
    }

    bool isFile(const std::string& filename) {
        return std::filesystem::is_regular_file(filename);
    }

    bool isDir(const std::string& filename) {
        return std::filesystem::is_directory(filename);
    }

    std::string getFileName(const std::string& fullFileName) {
        return std::filesystem::path(fullFileName).filename().string();
    }

    std::string getFilePath(const std::string& fullFileName) {
        std::string&& absFileName = std::filesystem::absolute(fullFileName).lexically_normal();
        return std::filesystem::path(absFileName).parent_path().string() + std::filesystem::path::preferred_separator;
    }

    std::string getAbsPath(const std::string& fileName) {
        std::filesystem::path fileNamePath(fileName);
        if (fileNamePath.is_absolute()) {
            return fileNamePath.lexically_normal();
        }
        return std::filesystem::absolute(fileName).lexically_normal();
    }

    std::string getAbsPath(const std::string& fileName, const std::string& baseDir) {
        std::filesystem::path fileNamePath(fileName);
        // If path is absolute, use it directly
        if (fileNamePath.is_absolute()) {
            return fileNamePath.lexically_normal();
        }
        
        // Otherwise, use base directory
        return std::filesystem::path((baseDir + std::filesystem::path::preferred_separator + fileName)).lexically_normal();
    }

    std::string createDir(const std::string& dir) {
        if (dir.empty()) {
            return "";
        }
        std::error_code ec;
        if (std::filesystem::exists(dir) && !std::filesystem::is_directory(dir)) {
            return "";
        }
        if (std::filesystem::create_directories(dir, ec)) {
            return getAbsPath(dir);
        }
        return "";
    }

    bool removeFile(const std::string& fileName) {
        std::error_code ec;
        return std::filesystem::remove(fileName, ec) && !ec;
    }

    std::string getHomePath() {
#if defined(_WIN32)
        char path[MAX_PATH];
        if (SUCCEEDED(SHGetFolderPath(NULL, CSIDL_PERSONAL, NULL, 0, path))) {
            return std::string(path);
        }
#else 
        const char* home = getenv("HOME");
        if (home != nullptr) {
            return std::string(home);
        }

        struct passwd *pwd = getpwuid(getuid());
        if (pwd != nullptr) {
            return std::string(pwd->pw_dir);
        }

        return "";
#endif
    }

    int64_t getFileMtime(const std::string& fileName) {
        struct stat st;
        if (0 != stat(fileName.c_str(), &st)) {
            return INT64_MAX;
        }
#if defined(__APPLE__) || defined(__MACH__) 
        return st.st_mtimespec.tv_sec;
#else
        return st.st_mtim.tv_sec;
#endif
    }

    int64_t getFileSize(const std::string& fileName) {
        struct stat st;
        if (0 != stat(fileName.c_str(), &st)) {
            return -1;
        }
        return st.st_size;
    }

    bool isGzFile(const std::string& fileName) {
        // Open file in binary mode for reading
        FILE* fp = fopen(fileName.c_str(), "rb");
        if (fp == nullptr) {
            LOG_ERROR("File %s open failed.", fileName.c_str());
            return false;
        }
        
        // Get file status to check file size
        struct stat st;
        if (0 != stat(fileName.c_str(), &st)) {
            fclose(fp);
            return false;
        }
        
        int64_t fileSize = st.st_size;
        // Check if file has at least 3 bytes for gzip header magic number
        if (fileSize > 3) {
            uint8_t tmpBuf[3] = {0};
            if (fread(tmpBuf, 3, 1, fp) != 1) {
                fclose(fp);
                return false;
            }
            fclose(fp);
            // Check for gzip magic bytes: 0x1F 0x8B
            return (tmpBuf[0] == 0x1F && tmpBuf[1] == 0x8B);
        }
        fclose(fp);
        return false;
    }

    std::string getFileNameFromGz(const std::string& gzFileName) {
        std::string originalName;
        do {
            std::ifstream file(gzFileName, std::ios::binary);
            if (!file.is_open()) {
                LOG_ERROR("Cannot open file:%s", gzFileName.c_str());
                break;
            }
            
            // Read gzip header
            unsigned char header[10];
            file.read(reinterpret_cast<char*>(header), 10);
            if (file.gcount() != 10) {
                LOG_ERROR("Not a valid gz file: %s", gzFileName.c_str());
                break;
            }
            
            // Check gzip magic bytes
            if (header[0] != 0x1F || header[1] != 0x8B) {
                LOG_ERROR("Not a valid gz file: %s", gzFileName.c_str());
                break;
            }
            
            // Check compression method (must be 8 = DEFLATE)
            if (header[2] != 8) {
                LOG_ERROR("Not supported compress format for gz file: %s", gzFileName.c_str());
                break;
            }
            
            // Get flags
            uint8_t flags = header[3];
            
            // Skip modification time (4 bytes), extra flags (1 byte), operating system (1 byte)
            // File pointer is already after the 10th byte
            
            // If extra field exists
            if (flags & 0x04) {
                uint8_t xlen[2];
                file.read(reinterpret_cast<char*>(xlen), 2);
                if (file.gcount() != 2) {
                    LOG_ERROR("Get extend info failed for gz file: %s", gzFileName.c_str());
                    break;
                }
                uint16_t extra_len = xlen[0] | (xlen[1] << 8);
                // Skip extra field
                file.seekg(extra_len, std::ios::cur);
            }
            
            // Read original filename (if exists)
            if (flags & 0x08) {
                std::vector<char> name;
                char ch;
                while (file.get(ch)) {
                    if (ch == '\0') {
                        break;
                    }
                    name.push_back(ch);
                }
                
                if (file.eof()) {
                    LOG_ERROR("Not a valid gz file: %s", gzFileName.c_str());
                    break;
                }
                
                originalName.assign(name.begin(), name.end());
                break;
            } else {
                LOG_ERROR("Cannot find file name in gz file: %s", gzFileName.c_str());
                break;
            }
        } while(0);

        
        if (originalName.empty()) {
            std::string gzName = getFileName(gzFileName);
            // Remove .gz extension
            if (gzName.length() > 3 && gzName.substr(gzName.length() - 3) == ".gz") {
                originalName = gzName.substr(0, gzName.length() - 3);
            } 
            // Remove .tgz extension and add .tar
            else if (gzName.length() > 4 && gzName.substr(gzName.length() - 4) == ".tgz") {
                originalName = gzName.substr(0, gzName.length() - 4) + ".tar";
            } 
            else {
                originalName = gzName;
            }
        }

        return originalName;
    }
}



