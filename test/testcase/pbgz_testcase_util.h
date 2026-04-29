#pragma once

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <thread>
#include <mutex>
#include "coder/coder.h"
#include "blocking_queue.h"
#include "io_block.h"
#include "utils/memory_util.h"
#include "io_wrapper.h"
#include <gmock/gmock.h>

// MockBlockingQueue - Simple wrapper to intercept queue operations
class MockBlockingQueue : public BlockingQueue<RoughIOBlock*> {
public:
    MockBlockingQueue() : BlockingQueue<RoughIOBlock*>(), blockSize(8 << 20), ioWriter(nullptr) {}

    MockBlockingQueue(uint32_t maxSize) : BlockingQueue<RoughIOBlock*>(maxSize),
                                          blockSize(8 << 20),
                                          ioWriter(nullptr) {}
    MockBlockingQueue(uint32_t maxSize, IOWriter* writer) {
        blockSize = maxSize;
        ioWriter = writer;
    }

    void setBlockSize(uint32_t size) {
        blockSize = size;
    }

    void setIOWriter(IOWriter* writer) {
        ioWriter = writer;
    }

    // Override get() to create blocks
    virtual RoughIOBlock* get() override {
        return MemoryUtil::safeNewClass<RoughIOBlock>(blockSize);
    }

    // Override push()
    virtual void push(RoughIOBlock* const &item) override {
        ioWriter->writeIO(item->getBuffer(), item->getDataLen());
        RoughIOBlock* dd = const_cast<RoughIOBlock*>(item);
        MemoryUtil::safeDeleteClass(dd);
    }

    virtual void pushForce(RoughIOBlock* const &item) override {
        push(item);
    }

private:
    uint32_t blockSize;
    IOWriter* ioWriter;
};


