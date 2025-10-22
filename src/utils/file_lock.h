#pragma once 

#include <string>
#include <fcntl.h>
#include <unistd.h>

#include <chrono>
#include <thread>
#include <sys/file.h>

class FileLock {
private:
    int fd = -1;
    std::string filename;
    bool isLocked = false;

public:
    enum LockType {
        SHARED_LOCK,    // 共享锁（读锁）
        EXCLUSIVE_LOCK  // 独占锁（写锁）
    };

    FileLock(const std::string& file) : filename(file) {}

    // 阻塞式获取锁（无限等待）
    bool lock(LockType type = EXCLUSIVE_LOCK) {
        if (isLocked) {
            return true; // 已经锁定
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        int operation = (type == SHARED_LOCK) ? LOCK_SH : LOCK_EX;
        
        // 阻塞式获取锁
        if (flock(fd, operation) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        isLocked = true;
        return true;
    }

    // 带超时的阻塞锁
    bool lock_for(std::chrono::milliseconds timeout, LockType type = EXCLUSIVE_LOCK) {
        auto start = std::chrono::steady_clock::now();
        
        while (std::chrono::steady_clock::now() - start < timeout) {
            if (try_lock(type)) {
                return true;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        // 最后一次尝试
        return try_lock(type);
    }

    // 非阻塞尝试获取锁
    bool try_lock(LockType type = EXCLUSIVE_LOCK) {
        if (isLocked) {
            return true;
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        int operation = (type == SHARED_LOCK) ? LOCK_SH : LOCK_EX;
        operation |= LOCK_NB; // 非阻塞模式
        
        if (flock(fd, operation) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        isLocked = true;
        return true;
    }

    void unlock() {
        if (isLocked && fd != -1) {
            flock(fd, LOCK_UN);
            close(fd);
            fd = -1;
            isLocked = false;
        }
    }

    bool locked() const {
        return isLocked;
    }

    ~FileLock() {
        unlock();
    }

    // 删除拷贝操作
    FileLock(const FileLock&) = delete;
    FileLock& operator=(const FileLock&) = delete;
};


class FcntlFileLock {
private:
    int fd = -1;
    std::string filename;
    bool is_locked = false;
    struct flock lockStruct;

public:
    enum LockType {
        READ_LOCK,   // 共享读锁
        WRITE_LOCK   // 独占写锁
    };

    FcntlFileLock(const std::string& file) : filename(file) {
        lockStruct.l_whence = SEEK_SET;
        lockStruct.l_start = 0;
        lockStruct.l_len = 0; // 0 表示锁定到文件末尾
    }

    // 阻塞锁定整个文件
    bool lock(LockType type) {
        return lock_region(0, 0, type); // 长度为0表示整个文件
    }

    // 阻塞锁定文件区域
    bool lock_region(off_t start, off_t length, LockType type) {
        if (is_locked) {
            return true;
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        lockStruct.l_type = (type == READ_LOCK) ? F_RDLCK : F_WRLCK;
        lockStruct.l_start = start;
        lockStruct.l_len = length;

        // F_SETLKW - 阻塞等待锁
        if (fcntl(fd, F_SETLKW, &lockStruct) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        is_locked = true;
        return true;
    }

    // 带超时的区域锁定
    bool lock_region_for(off_t start, off_t length, LockType type, 
                        std::chrono::milliseconds timeout) {
        auto start_time = std::chrono::steady_clock::now();
        
        while (std::chrono::steady_clock::now() - start_time < timeout) {
            if (try_lock_region(start, length, type)) {
                return true;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        
        return try_lock_region(start, length, type);
    }

    // 非阻塞尝试区域锁定
    bool try_lock_region(off_t start, off_t length, LockType type) {
        if (is_locked) {
            return true;
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        lockStruct.l_type = (type == READ_LOCK) ? F_RDLCK : F_WRLCK;
        lockStruct.l_start = start;
        lockStruct.l_len = length;

        // F_SETLK - 非阻塞尝试
        if (fcntl(fd, F_SETLK, &lockStruct) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        is_locked = true;
        return true;
    }

    void unlock() {
        if (is_locked && fd != -1) {
            lockStruct.l_type = F_UNLCK;
            fcntl(fd, F_SETLK, &lockStruct);
            close(fd);
            fd = -1;
            is_locked = false;
        }
    }

    ~FcntlFileLock() {
        unlock();
    }

    FcntlFileLock(const FcntlFileLock&) = delete;
    FcntlFileLock& operator=(const FcntlFileLock&) = delete;
};