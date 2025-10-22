#pragma once
#include <mutex>
#include <condition_variable>

class WaitEvent {
private:
    std::mutex mtx;
    std::condition_variable cv;
    bool event_occurred = false;

public:
    void wait() {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this]() { return event_occurred; });
    }

    bool wait_for(std::chrono::milliseconds timeout) {
        std::unique_lock<std::mutex> lock(mtx);
        return cv.wait_for(lock, timeout, [this]() { return event_occurred; });
    }

    bool wait_until(std::chrono::steady_clock::time_point deadline) {
        std::unique_lock<std::mutex> lock(mtx);
        return cv.wait_until(lock, deadline, [this]() { return event_occurred; });
    }

    void start() {
        {
            std::lock_guard<std::mutex> lock(mtx);
            event_occurred = true;
        }
        cv.notify_all();
    }

    void reset() {
        std::lock_guard<std::mutex> lock(mtx);
        event_occurred = false;
    }
};