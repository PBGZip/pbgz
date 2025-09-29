#pragma once

#include <string>

namespace PathUtil {
    /// @brief 判断文件名是否以后缀结尾
    /// @param fileName  文件名
    /// @param suffix    后缀
    /// @return bool     true/false
    bool suffixCheck(const std::string& fileName, const std::string& suffix);

    /// @brief 判断文件是否存在
    /// @param fileName  文件名
    /// @return bool     true:存在 / false:不存在
    bool fileExists(const std::string& fileName);

    /// @brief 判断文件是否可读
    /// @param fileName 
    /// @return  bool  true/false
    bool fileReadble(const std::string& fileName);

    /// @brief 判断文件是否可写
    /// @param fileName 
    /// @return  bool  true/false
    bool fileWriteable(const std::string& fileName);

    /// @brief 根据文件名获取文件的名字，不带目录
    /// @param fullFileName 
    /// @return 不带目录的文件名
    std::string getFileName(const std::string& fullFileName);

    /// @brief 根据文件名获取文件的目录
    /// @param fullFileName 
    /// @return 文件所在的目录
    std::string getFilePath(const std::string& fullFileName);

    /// @brief 获取文件的绝对路径
    /// @param fileName 
    /// @return 文件的决定路径
    std::string getAbsPath(const std::string& fileName);

    /// @brief 创建目录
    /// @param dir 
    /// @return 创建的目录名称的绝对路径
    std::string createDir(const std::string& dir);

    /// @brief 判断文件是否是目录
    /// @param filename 
    /// @return  bool  true/false
    bool isDir(const std::string& filename);

    /// @brief 判断文件是否是普通文件
    /// @param filename 
    /// @return bool  true/false
    bool isFile(const std::string& filename);

    /// @brief 删除文件
    /// @param fileName 
    /// @return bool  true/false
    bool removeFile(const std::string& fileName);
}

