
#include <cstring>
#include <filesystem>
#include <system_error>

#include "path_util.h"

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
        return std::filesystem::path(fullFileName).parent_path().string() + std::filesystem::path::preferred_separator;
    }

    std::string getAbsPath(const std::string& fileName) {
        return std::filesystem::absolute(fileName).lexically_normal();
    }

    std::string getAbsPath(const std::string& fileName, const std::string& baseDir) {
        std::filesystem::path fileNamePath(fileName);
        // 如果路径是绝对的，直接使用
        if (fileNamePath.is_absolute()) {
            return fileNamePath.lexically_normal();
        }
        
        // 否则基于基准目录
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
        if (std::filesystem::create_directory(dir, ec)) {
            return getAbsPath(dir);
        }
        return "";
    }

    bool removeFile(const std::string& fileName) {
        std::error_code ec;
        return std::filesystem::remove(fileName, ec) && !ec;
    }

}
