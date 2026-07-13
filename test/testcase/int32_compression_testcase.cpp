/*
 * int32_compression_testcase.cpp - Test cases for int32 array compression comparison
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
#include <string>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
#include <sstream>

#include "coder/coder_bwt_cm.h"
#include "coder/coder_io.h"
#include "pbgz_errno.h"
#include "utils/memory_util.h"

class Int32CompressionTest : public ::testing::Test 
{
protected:
    void SetUp() override 
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    }

    void TearDown() override 
    {
    }
    
    std::vector<int32_t> generate_test_int32_array(size_t size, int32_t min_val, int32_t max_val) {
        std::vector<int32_t> result;
        result.reserve(size);
        for (size_t i = 0; i < size; ++i) {
            result.push_back(min_val + (rand() % (max_val - min_val + 1)));
        }
        return result;
    }
    
    std::string int32_array_to_string(const std::vector<int32_t>& data, const std::string& delimiter = ",") {
        std::stringstream ss;
        for (size_t i = 0; i < data.size(); ++i) {
            if (i > 0) {
                ss << delimiter;
            }
            ss << data[i];
        }
        return ss.str();
    }
    
    std::vector<int32_t> string_to_int32_array(const std::string& str, const std::string& delimiter = ",") {
        std::vector<int32_t> result;
        std::string token;
        std::stringstream ss(str);
        
        while (std::getline(ss, token, delimiter[0])) {
            if (!token.empty()) {
                result.push_back(std::stoi(token));
            }
        }
        return result;
    }
    
    double calculate_compression_ratio(size_t original_size, size_t compressed_size) {
        if (compressed_size == 0) return 0.0;
        return (double)compressed_size / original_size;
    }
};

TEST_F(Int32CompressionTest, CompareNumericAndStringCompression) {
    const size_t array_size = 10000;
    const int32_t min_val = 1;
    const int32_t max_val = 1000;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Numeric compression:\n");
    printf("  Original size: %zu bytes\n", numeric_original_size);
    printf("  Compressed size: %zu bytes\n", numeric_compressed_size);
    printf("  Compression ratio: %.4f\n", numeric_ratio);
    printf("  Compression savings: %.2f%%\n", (1.0 - numeric_ratio) * 100.0);
    
    printf("\nString compression:\n");
    printf("  Original size: %zu bytes\n", string_original_size);
    printf("  Compressed size: %zu bytes\n", string_compressed_size);
    printf("  Compression ratio: %.4f\n", string_ratio);
    printf("  Compression savings: %.2f%%\n", (1.0 - string_ratio) * 100.0);
    
    printf("\nComparison:\n");
    printf("  Original size ratio (string/numeric): %.2f\n", (double)string_original_size / numeric_original_size);
    printf("  Compressed size ratio (string/numeric): %.2f\n", (double)string_compressed_size / numeric_compressed_size);
    printf("  Difference in compression ratio: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
    EXPECT_LE(numeric_ratio, 1.0);
    EXPECT_LE(string_ratio, 1.0);
}

TEST_F(Int32CompressionTest, CompareSmallRangeNumbers) {
    const size_t array_size = 5000;
    const int32_t min_val = 1;
    const int32_t max_val = 10;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Small range numbers (1-10):\n");
    printf("Numeric compression ratio: %.4f (%.2f%% savings)\n", numeric_ratio, (1.0 - numeric_ratio) * 100.0);
    printf("String compression ratio: %.4f (%.2f%% savings)\n", string_ratio, (1.0 - string_ratio) * 100.0);
    printf("Difference: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
}

TEST_F(Int32CompressionTest, CompareLargeRangeNumbers) {
    const size_t array_size = 5000;
    const int32_t min_val = 1;
    const int32_t max_val = 1000000;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Large range numbers (1-1000000):\n");
    printf("Numeric compression ratio: %.4f (%.2f%% savings)\n", numeric_ratio, (1.0 - numeric_ratio) * 100.0);
    printf("String compression ratio: %.4f (%.2f%% savings)\n", string_ratio, (1.0 - string_ratio) * 100.0);
    printf("Difference: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
}

TEST_F(Int32CompressionTest, CompareSequentialNumbers) {
    const size_t array_size = 10000;
    
    std::vector<int32_t> test_data;
    test_data.reserve(array_size);
    for (size_t i = 0; i < array_size; ++i) {
        test_data.push_back(static_cast<int32_t>(i));
    }
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Sequential numbers (0-%zu):\n", array_size);
    printf("Numeric compression ratio: %.4f (%.2f%% savings)\n", numeric_ratio, (1.0 - numeric_ratio) * 100.0);
    printf("String compression ratio: %.4f (%.2f%% savings)\n", string_ratio, (1.0 - string_ratio) * 100.0);
    printf("Difference: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
}

TEST_F(Int32CompressionTest, CompareRepetitiveNumbers) {
    const size_t array_size = 10000;
    
    std::vector<int32_t> test_data;
    test_data.reserve(array_size);
    const int32_t pattern[] = {1, 5, 10, 100, 1000};
    
    for (size_t i = 0; i < array_size; ++i) {
        test_data.push_back(pattern[i % 5]);
    }
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Repetitive numbers:\n");
    printf("Numeric compression ratio: %.4f (%.2f%% savings)\n", numeric_ratio, (1.0 - numeric_ratio) * 100.0);
    printf("String compression ratio: %.4f (%.2f%% savings)\n", string_ratio, (1.0 - string_ratio) * 100.0);
    printf("Difference: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
}

TEST_F(Int32CompressionTest, CompareDifferentDelimiters) {
    const size_t array_size = 5000;
    const int32_t min_val = 1;
    const int32_t max_val = 100;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    std::vector<std::string> delimiters = {",", ";", " ", "|", "\n"};
    
    printf("Testing different delimiters:\n");
    for (const auto& delimiter : delimiters) {
        uint32_t buffer_size = 2 << 20; 
        std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
        std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
        std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
        
        std::string test_string = int32_array_to_string(test_data, delimiter);
        size_t string_original_size = test_string.size();
        
        string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
        string_encoder->encode_flush();
        size_t string_compressed_size = string_io->data_len;
        
        double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
        
        printf("Delimiter '%s': ratio=%.4f (%.2f%% savings)\n", 
               delimiter.c_str(), string_ratio, (1.0 - string_ratio) * 100.0);
        
        EXPECT_GT(string_compressed_size, 0);
    }
}

TEST_F(Int32CompressionTest, CompareMixedPositiveNegativeNumbers) {
    const size_t array_size = 5000;
    const int32_t min_val = -1000;
    const int32_t max_val = 1000;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    size_t numeric_compressed_size = numeric_io->data_len;
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    size_t string_compressed_size = string_io->data_len;
    
    double numeric_ratio = calculate_compression_ratio(numeric_original_size, numeric_compressed_size);
    double string_ratio = calculate_compression_ratio(string_original_size, string_compressed_size);
    
    printf("Mixed positive/negative numbers (-1000 to 1000):\n");
    printf("Numeric compression ratio: %.4f (%.2f%% savings)\n", numeric_ratio, (1.0 - numeric_ratio) * 100.0);
    printf("String compression ratio: %.4f (%.2f%% savings)\n", string_ratio, (1.0 - string_ratio) * 100.0);
    printf("Difference: %.4f\n", string_ratio - numeric_ratio);
    
    EXPECT_GT(numeric_compressed_size, 0);
    EXPECT_GT(string_compressed_size, 0);
}

TEST_F(Int32CompressionTest, VerifyCompressionCorrectness) {
    const size_t array_size = 1000;
    const int32_t min_val = 1;
    const int32_t max_val = 100;
    
    std::vector<int32_t> test_data = generate_test_int32_array(array_size, min_val, max_val);
    
    uint32_t buffer_size = 2 << 20; 
    std::unique_ptr<uint8_t[]> numeric_output(new uint8_t[buffer_size]);
    std::unique_ptr<uint8_t[]> string_output(new uint8_t[buffer_size]);
    
    std::unique_ptr<coder_io> numeric_io(new coder_io(numeric_output.get(), buffer_size));
    std::unique_ptr<coder_io> string_io(new coder_io(string_output.get(), buffer_size));
    
    std::unique_ptr<coder_bwt_cm> numeric_encoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_encoder(new coder_bwt_cm(string_io.get()));
    
    size_t numeric_original_size = test_data.size() * sizeof(int32_t);
    numeric_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), numeric_original_size);
    numeric_encoder->encode_flush();
    
    std::string test_string = int32_array_to_string(test_data, ",");
    size_t string_original_size = test_string.size();
    string_encoder->encode_line(reinterpret_cast<const uint8_t*>(test_string.data()), string_original_size);
    string_encoder->encode_flush();
    
    numeric_io->data_len = 0;
    string_io->data_len = 0;
    
    std::unique_ptr<coder_bwt_cm> numeric_decoder(new coder_bwt_cm(numeric_io.get()));
    std::unique_ptr<coder_bwt_cm> string_decoder(new coder_bwt_cm(string_io.get()));
    
    std::vector<uint8_t> numeric_decoded(numeric_original_size);
    int32_t numeric_decoded_len = numeric_decoder->decode_line(numeric_decoded.data(), numeric_original_size);
    
    std::vector<uint8_t> string_decoded(string_original_size);
    int32_t string_decoded_len = string_decoder->decode_line(string_decoded.data(), string_original_size);
    
    EXPECT_EQ(numeric_decoded_len, numeric_original_size);
    EXPECT_EQ(string_decoded_len, string_original_size);
    
    EXPECT_EQ(memcmp(numeric_decoded.data(), test_data.data(), numeric_original_size), 0);
    
    std::string decoded_string(string_decoded.begin(), string_decoded.begin() + string_decoded_len);
    EXPECT_EQ(decoded_string, test_string);
}