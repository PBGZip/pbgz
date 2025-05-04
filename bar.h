#ifndef _BAR_H_
#define _BAR_H_

#include "wait_event.h"
#include <iostream>
#include <thread>

/* 进度条 */

/* 自动读值更新进度直到退出，没加锁，不支持多线程 */
class guard_bar
{
public:
    guard_bar(const int64_t &total_size, int64_t *curr_size, std::string bar_info)
    {
        this->exit = false;
        this->info = bar_info;
        const std::string tips = this->info;
        wait_event &event = this->wait;
        bool &eflag = this->exit;
        this->guard = new std::thread([curr_size, total_size, tips, &event, &eflag]()
                                      {
                                          int64_t current;
                                          fflush(stderr);
                                          event.wait();
                                          do
                                          {
                                              current = std::min((((*curr_size) * 100) / total_size), (int64_t)100);
                                              fprintf(stderr, "\033[?25l%s --[%lld%%]--\r", tips.c_str(), current);
                                              fflush(stderr);
                                              usleep(20000);
                                              if (eflag)
                                                  break;
                                          } while (true);
                                      });
    }

    void start()
    {
        this->wait.wakeup();
    }

    void done(std::string append = "")
    {
        this->exit = true;
        if (this->guard->joinable())
            this->guard->join();
        fprintf(stderr, "\033[?25l%-s --[%d%%]-- %s\r", this->info.c_str(), 100, append.c_str());
        fprintf(stderr, "\n\033[?25h");
        fflush(stderr);
    }

    virtual ~guard_bar()
    {
        delete this->guard;
    }

private:
    std::string info;
    std::thread *guard;
    bool exit;
    wait_event wait;
};

#endif
