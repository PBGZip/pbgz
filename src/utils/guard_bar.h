/*
 * guard_bar.h - Header file for pbgz project
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

#include <stdint.h>
#include <string>
#include <thread>
#include <future>
#include "unistd.h"

#include "wait.h"

class GuardBar {
public:
    GuardBar(const int64_t& totalSize, int64_t& currSize, const std::string& barInfo) {
        isExit = false;
        info = barInfo;
        const std::string tips = barInfo;
        bool &extFlag = isExit;
        std::future<void>& event = future;
        guard = std::thread([&currSize, totalSize, tips, &event, &extFlag](){
            fflush(stderr);
            event.wait();
            do {
                int64_t current = std::min((currSize * 100)/totalSize, (int64_t)100);
                fprintf(stderr, "\033[?25l%s --[%lld%%]--\r", tips.c_str(), (long long int)current);
                fflush(stderr);
                usleep(20000);
                if (extFlag) {
                    break;
                }
            } while(true);
        });
    }

    void start() {
        promise.set_value(); 
    }

    void done(std::string& append) {
        isExit = true;
        if (guard.joinable()) {
            guard.join();
        }

        fprintf(stderr, "\33[?25l%-s --[%d%%]-- %s\r", info.c_str(), 100, append.c_str());
        fprintf(stderr, "\n\033[?25h");
        fflush(stderr);
    }

    virtual ~GuardBar() { }

private:
    bool isExit;
    std::string info;
    std::thread guard;
    // Create promise-future pair
    std::promise<void> promise;
    std::future<void> future = promise.get_future();
};
