#pragma once
#include <cstdint>
#include <chrono>

class Timer {
public:
    /**
     * @brief Construct a new Timer object
     * @param autostart If true, timer starts upon construction
     */
    explicit Timer(bool autostart = false) {
        if (autostart) {
            reset();
        }
    }

    /**
     * @brief Reset the timer to current time
     */
    void reset() {
        _start = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Get elapsed milliseconds since reset
     * @return uint32_t Elapsed time in milliseconds
     */
    uint32_t elapsed() const {
        auto now = std::chrono::high_resolution_clock::now();
        return static_cast<uint32_t>(
            std::chrono::duration_cast<std::chrono::milliseconds>(now - _start).count()
        );
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> _start;
};
