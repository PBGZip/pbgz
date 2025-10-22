#pragma once

#include <stdint.h>
#include <string>
#include <thread>
#include <future>

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
    // 创建 promise-future 对
    std::promise<void> promise;
    std::future<void> future = promise.get_future();
};