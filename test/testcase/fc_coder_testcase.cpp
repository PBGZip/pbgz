/*
 * fc_coder_testcase.cpp - Test cases for FC coder
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

#include <gtest/gtest.h>
#include <fstream>

#include "coder_fc.h"
#include "utils/memory_util.h"

class FcCoderTest : public ::testing::Test {
public:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
     }

    void TearDown() override { 

    }
};

TEST_F(FcCoderTest, FcDecodeText) {
    uint32_t length = 2 << 20; 
    uint8_t* pOut =static_cast<uint8_t*>(malloc(length));
    coder_io*  fcIo = new coder_io(pOut, length);
    coder_fc*  fcCoder = new coder_fc(fcIo);
    char* text = "1212121212112121212112121212121221212121212121211212121";
    fcCoder->encode_line((uint8_t*)text, strlen(text));
    EXPECT_EQ(19, fcIo->data_len);
    delete fcCoder;
    delete fcIo;
    free(pOut);
}

TEST_F(FcCoderTest, FcDecodeRandom) {
    uint32_t length = 10 << 20; 
    uint8_t* pOut =static_cast<uint8_t*>(malloc(length));
    coder_io*  fcIo = new coder_io(pOut, length);
    coder_fc*  fcCoder = new coder_fc(fcIo);
    std::ifstream urandom("/dev/urandom", std::ios::binary);
    uint32_t count = 1024;
    std::vector<unsigned char> buffer(count);
    urandom.read(reinterpret_cast<char*>(buffer.data()), count);
    if (!urandom) {
        throw std::runtime_error("Read data from /dev/urandom failed.");
    }
    fcCoder->encode_line(buffer.data(), count);
    EXPECT_LT(1024, fcIo->data_len);
    delete fcCoder;
    delete fcIo;
    free(pOut);
}