
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
     * Push element to queue, block if queue is full
     **/
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        conditonVar.wait(lock, [this]{return (0 == maxQueueSize) ? true : (dataQueue.size() < maxQueueSize);});
        dataQueue.push(std::move(item));
        conditonVar.notify_one();
    }

    /*
     * Push element to queue without blocking
     **/
    void pushForce(const T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        dataQueue.push(std::move(item));
        conditonVar.notify_one();
    }

    /*
     * Pop element from queue, block if queue is empty
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
     * Return whether queue is empty
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
     * Return number of elements in queue
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
