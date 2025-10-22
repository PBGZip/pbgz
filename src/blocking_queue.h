
#pragma once

#include <queue>
#include <mutex>
#include <stdint.h>
#include <condition_variable>

template<typename T>
class BlockingQueue {

public:
    BlockingQueue(const BlockingQueue&) = delete;
    BlockingQueue& operator=(const BlockingQueue&) = delete;

    BlockingQueue() : maxQueueSize(0) {}

    BlockingQueue(uint32_t maxSize) : maxQueueSize(maxSize) { }

    /*
     * 将元素加入队列,如果队列满就阻塞
     **/
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        conditonVar.wait(lock, [this]{return (0 == maxQueueSize) ? true : (dataQueue.size() < maxQueueSize);});
        dataQueue.push(std::move(item));
        conditonVar.notify_one();
    }

    /*
     * 将元素加入队列, 不带阻塞
     **/
    void pushForce(const T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        dataQueue.push(std::move(item));
        conditonVar.notify_one();
    }

    /*
     * 从队列中弹出一个元素,如果队列为空就阻塞
     **/
    T get() {
        std::unique_lock<std::mutex> lock(mutex);
        conditonVar.wait(lock, [this]{return !this->dataQueue.empty();});
        T value = std::move(dataQueue.front());
        dataQueue.pop();
        conditonVar.notify_one();
        return value;
    }

    /**
     * 返回队列是否为空
     */
    bool empty() {
        std::lock_guard<std::mutex> lock(mutex);
        return dataQueue.empty();
    }

    void clear() {
        std::lock_guard<std::mutex> lock(mutex);
        while(!dataQueue.empty()) {
            dataQueue.pop();
        }
    }

    /*
     * 返回队列中元素数个
     * */
    size_t size() {
        std::lock_guard<std::mutex> lock(mutex);
        return dataQueue.size();
    }

    void setCapility(uint32_t capility) {
        maxQueueSize = capility;
    }

    uint32_t getCapility() {
        return maxQueueSize;
    }

private:
    uint32_t maxQueueSize;
    mutable std::mutex mutex;
    mutable std::condition_variable conditonVar;
    std::queue<T> dataQueue;
}; 
