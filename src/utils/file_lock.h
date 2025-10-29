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
        SHARED_LOCK,    // Shared lock (read lock)
        EXCLUSIVE_LOCK  // Exclusive lock (write lock)
    };

    FileLock(const std::string& file) : filename(file) {}

    // Blocking lock acquisition (infinite wait)
    bool lock(LockType type = EXCLUSIVE_LOCK) {
        if (isLocked) {
            return true; // Already locked
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        int operation = (type == SHARED_LOCK) ? LOCK_SH : LOCK_EX;
        
        // Blocking lock acquisition
        if (flock(fd, operation) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        isLocked = true;
        return true;
    }

    // Blocking lock with timeout
    bool lock_for(std::chrono::milliseconds timeout, LockType type = EXCLUSIVE_LOCK) {
        auto start = std::chrono::steady_clock::now();
        
        while (std::chrono::steady_clock::now() - start < timeout) {
            if (try_lock(type)) {
                return true;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        // Final attempt
        return try_lock(type);
    }

    // Non-blocking try lock acquisition
    bool try_lock(LockType type = EXCLUSIVE_LOCK) {
        if (isLocked) {
            return true;
        }

        fd = open(filename.c_str(), O_RDWR | O_CREAT, 0644);
        if (fd == -1) {
            return false;
        }

        int operation = (type == SHARED_LOCK) ? LOCK_SH : LOCK_EX;
        operation |= LOCK_NB; // Non-blocking mode
        
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

    // Delete copy operations
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
        READ_LOCK,   // Shared read lock
        WRITE_LOCK   // Exclusive write lock
    };

    FcntlFileLock(const std::string& file) : filename(file) {
        lockStruct.l_whence = SEEK_SET;
        lockStruct.l_start = 0;
        lockStruct.l_len = 0; // 0 means lock to end of file
    }

    // Block lock entire file
    bool lock(LockType type) {
        return lock_region(0, 0, type); // Length 0 means entire file
    }

    // Block lock file region
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

        // F_SETLKW - Blocking wait for lock
        if (fcntl(fd, F_SETLKW, &lockStruct) == -1) {
            close(fd);
            fd = -1;
            return false;
        }

        is_locked = true;
        return true;
    }

    // Region lock with timeout
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

    // Non-blocking try region lock
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

        // F_SETLK - Non-blocking try
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