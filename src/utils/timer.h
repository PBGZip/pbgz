/*
 * timer.h - Header file for pbgz project
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
