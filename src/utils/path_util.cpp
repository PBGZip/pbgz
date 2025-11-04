
#include <cstring>
#include <filesystem>
#include <system_error>
#include <pwd.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "path_util.h"
#include "log/logger.h"


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
        if (std::filesystem::create_directory(dir, ec)) {
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

    bool isGzFile(std::string& fileName) {
        FILE* fp = fopen64(fileName.c_str(), "rb");
        if (fp == nullptr) {
            LOG_ERROR("File %s open failed.", fileName.c_str());
            return false;
        }
        int64_t fileSize = PathUtil::getFileSize(fileName);
        if (fileSize > 3) {
            uint8_t tmpBuf[3] = {0};
            fread(tmpBuf, 3, 1, fp);
            return (tmpBuf[0] == 0x1F && tmpBuf[1] == 0x8B);
        }
        return false;
    }
}



