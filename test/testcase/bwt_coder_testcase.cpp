/*
 * bwt_coder_testcase.cpp - Test cases for BWT coder
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

#include "coder/coder_bwt_cm.h"
#include "coder/coder_io.h"
#include "pbgz_errno.h"
#include "utils/memory_util.h"

class BwtCoderTest : public ::testing::Test 
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
};


// Test basic encode/decode cycle
TEST_F(BwtCoderTest, EncodeDecodeBasic) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test data
    const std::string test_data = "ATCGATCGATCGATCG";
    std::vector<uint8_t> input(test_data.begin(), test_data.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(input.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, input.size());
    EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), test_data);
}

// Test single character
TEST_F(BwtCoderTest, SingleCharacter) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test single character
    const uint8_t single_char = 'A';
    encoder.encode_line(&single_char, 1);
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    uint8_t output;
    int32_t decoded_len = decoder.decode_line(&output, 1);
    
    EXPECT_EQ(decoded_len, 1);
    EXPECT_EQ(output, single_char);
    
}

// Test repeated characters
TEST_F(BwtCoderTest, RepeatedCharacters) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test repeated characters
    const std::string repeated_data = "AAAAAAAAAAAAAAAAAAAA";
    std::vector<uint8_t> input(repeated_data.begin(), repeated_data.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(input.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, input.size());
    EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), repeated_data);
}

// Test random data
TEST_F(BwtCoderTest, RandomData) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Generate random test data
    std::vector<uint8_t> random_data;
    for (int i = 0; i < 1000; ++i) {
        random_data.push_back(static_cast<uint8_t>(rand() % 256));
    }
    
    // Encode
    encoder.encode_line(random_data.data(), random_data.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(random_data.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, random_data.size());
    EXPECT_EQ(memcmp(output.data(), random_data.data(), random_data.size()), 0);
}

// Test DNA sequence data (typical use case)
TEST_F(BwtCoderTest, DNASequence) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test DNA sequence
    const std::string dna_sequence = "ATGCGATCGTAGCTAGCTAGCGATCGATCGTACGATCGATCGATCGATCGATCGATCG";
    std::vector<uint8_t> input(dna_sequence.begin(), dna_sequence.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(input.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, input.size());
    EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), dna_sequence);
}

// Test split character functionality
TEST_F(BwtCoderTest, SplitCharacter) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test data with split characters
    const std::string test_data = "A,T,C,G,A,T,C,G";
    const uint8_t split_char = ',';
    std::vector<uint8_t> input(test_data.begin(), test_data.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode with split character
    std::vector<uint8_t> output(input.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size(), split_char);
    
    EXPECT_GT(decoded_len, 0);
    EXPECT_LE(decoded_len, input.size());
    
    // Verify output contains split character
    bool found_split = false;
    for (int i = 0; i < decoded_len; ++i) {
        if (output[i] == split_char) {
            found_split = true;
            break;
        }
    }
    EXPECT_TRUE(found_split);
    
}

// Test multiple encode calls
TEST_F(BwtCoderTest, MultipleEncodeCalls) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test multiple small encode calls
    const std::string chunk1 = "ATCGATCG";
    const std::string chunk2 = "GCTAGCTA";
    const std::string chunk3 = "TTTTAAAA";
    
    // Encode in chunks
    encoder.encode_line(reinterpret_cast<const uint8_t*>(chunk1.data()), chunk1.size());
    encoder.encode_line(reinterpret_cast<const uint8_t*>(chunk2.data()), chunk2.size());
    encoder.encode_line(reinterpret_cast<const uint8_t*>(chunk3.data()), chunk3.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode all at once
    std::string combined = chunk1 + chunk2 + chunk3;
    std::vector<uint8_t> output(combined.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, combined.size());
    EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), combined);
    
}

// Test large data (multiple blocks)
TEST_F(BwtCoderTest, LargeData) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Generate large test data (should span multiple blocks)
    std::vector<uint8_t> large_data;
    const std::string pattern = "ATCGATCGGCTAGCTAG";
    for (int i = 0; i < 100000; ++i) {
        for (char c : pattern) {
            large_data.push_back(static_cast<uint8_t>(c));
        }
    }
    
    // Encode
    encoder.encode_line(large_data.data(), large_data.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(large_data.size());
    int32_t decoded_len = 0;
    
    // Decode in chunks if needed
    int32_t chunk_size = 1024;
    while (decoded_len < large_data.size()) {
        int32_t chunk_len = decoder.decode_line(
            output.data() + decoded_len, 
            std::min(chunk_size, static_cast<int>(large_data.size() - decoded_len))
        );
        
        if (chunk_len == 0) break; // End of data
        decoded_len += chunk_len;
    }
    
    EXPECT_EQ(decoded_len, large_data.size());
    EXPECT_EQ(memcmp(output.data(), large_data.data(), large_data.size()), 0);
    
}

// Test decode with insufficient buffer
TEST_F(BwtCoderTest, DecodeInsufficientBuffer) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test data
    const std::string test_data = "ATCGATCGATCG";
    std::vector<uint8_t> input(test_data.begin(), test_data.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode with insufficient buffer
    std::vector<uint8_t> small_output(4); // Smaller than input
    int32_t decoded_len = decoder.decode_line(small_output.data(), small_output.size());
    
    EXPECT_EQ(decoded_len, small_output.size());
    EXPECT_EQ(std::string(small_output.begin(), small_output.end()), test_data.substr(0, 4));
    
}

// Test decode input length tracking
TEST_F(BwtCoderTest, DecodeInputLength) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test data
    const std::string test_data = "ATCGATCG";
    std::vector<uint8_t> input(test_data.begin(), test_data.end());
    
    // Encode
    encoder.encode_line(input.data(), input.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Check input length before decoding
    int32_t initial_inlen = decoder.decode_inlen();
    
    // Decode
    std::vector<uint8_t> output(input.size());
    decoder.decode_line(output.data(), output.size());
    
    // Check input length after decoding
    int32_t final_inlen = decoder.decode_inlen();
    
    EXPECT_GT(final_inlen, initial_inlen);
}

// Test memory cleanup
TEST_F(BwtCoderTest, MemoryCleanup) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    
    // Test that encoder and decoder clean up properly when destructed
    {
        coder_bwt_cm encoder(coderIo.get());
        const std::string test_data = "ATCGATCG";
        encoder.encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), test_data.size());
        encoder.encode_flush();
    }
    
    // Should not crash when encoder goes out of scope
    SUCCEED();
}

// Test edge cases with special characters
TEST_F(BwtCoderTest, SpecialCharacters) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test data with all possible byte values
    std::vector<uint8_t> special_data;
    for (int i = 0; i < 256; ++i) {
        special_data.push_back(static_cast<uint8_t>(i));
    }
    
    // Encode
    encoder.encode_line(special_data.data(), special_data.size());
    encoder.encode_flush();
    
    // Setup for decoding
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    
    // Decode
    std::vector<uint8_t> output(special_data.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, special_data.size());
    EXPECT_EQ(memcmp(output.data(), special_data.data(), special_data.size()), 0);
}

// Test compression mode detection
TEST_F(BwtCoderTest, CompressionModeDetection) 
{
    uint32_t length = 2 << 20; 
    
    // Test with compression mode
    std::unique_ptr<uint8_t[]> pOut1(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo1(new coder_io(pOut1.get(), length));
    coderIo1->m = coder_io::MENC;
    EXPECT_NO_THROW({
        coder_bwt_cm encoder(coderIo1.get());
    });
    
    // Test with decompression mode
    std::unique_ptr<uint8_t[]> pOut2(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo2(new coder_io(pOut2.get(), length));
    coderIo2->m = coder_io::MDEC;
    EXPECT_NO_THROW({
        coder_bwt_cm decoder(coderIo2.get());
    });
    
    // Test with unset mode (should default appropriately)
    std::unique_ptr<uint8_t[]> pOut3(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo3(new coder_io(pOut3.get(), length));
    coderIo3->m = coder_io::MUNSET;
    EXPECT_NO_THROW({
        coder_bwt_cm coder(coderIo3.get());
    });
}

// Test metadata handling
TEST_F(BwtCoderTest, MetadataHandling) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    // Test with various metadata settings
    coderIo->meta["level"] = 6;
    coderIo->meta["magic"] = "BWT_CM";
    coderIo->meta["custom_param"] = "test_value";
    
    EXPECT_NO_THROW({
        coder_bwt_cm encoder(coderIo.get());
    });
    
    // Test encoding with metadata
    coder_bwt_cm encoder(coderIo.get());
    const std::string test_data = "ATCGATCGATCG";
    encoder.encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), test_data.size());
    encoder.encode_flush();
    
    // Verify metadata is preserved
    EXPECT_EQ(coderIo->get_level(), 6);
    EXPECT_EQ(coderIo->get_magic(), "coder_bwt_cm");
}

// Test very large single block
TEST_F(BwtCoderTest, VeryLargeSingleBlock) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Generate data close to block size limit
    std::vector<uint8_t> large_block(32000 - 100, 'A'); // Just under typical block size
    for (size_t i = 0; i < large_block.size(); i += 100) {
        large_block[i] = static_cast<uint8_t>('A' + (i / 100) % 4); // Add some variation
    }
    
    // Encode
    EXPECT_NO_THROW({
        encoder.encode_line(large_block.data(), large_block.size());
        encoder.encode_flush();
    });
    
    // Decode
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    std::vector<uint8_t> output(large_block.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, large_block.size());
    EXPECT_EQ(memcmp(output.data(), large_block.data(), large_block.size()), 0);
}

// Test performance with repetitive data
TEST_F(BwtCoderTest, RepetitiveDataCompression) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());
    
    // Test highly repetitive data (should compress well)
    std::vector<uint8_t> repetitive_data(10000, 'A');
    
    // Encode
    encoder.encode_line(repetitive_data.data(), repetitive_data.size());
    encoder.encode_flush();
    
    // Decode
    coderIo->data_len = 0; // Reset data length for decoding
    coder_bwt_cm decoder(coderIo.get());
    std::vector<uint8_t> output(repetitive_data.size());
    int32_t decoded_len = decoder.decode_line(output.data(), output.size());
    
    EXPECT_EQ(decoded_len, repetitive_data.size());
    EXPECT_EQ(memcmp(output.data(), repetitive_data.data(), repetitive_data.size()), 0);
}

// Test mixed compression levels in same session
TEST_F(BwtCoderTest, MixedCompressionLevels) 
{
    const std::string test_data = "ATCGATCGATCGATCG";
    
    // Test with different compression levels
    for (int level = 1; level <= 9; ++level) {
        uint32_t length = 2 << 20; 
        std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
        std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
        coderIo->meta["level"] = level;
        
        coder_bwt_cm encoder(coderIo.get());
        std::vector<uint8_t> input(test_data.begin(), test_data.end());
        
        encoder.encode_line(input.data(), input.size());
        encoder.encode_flush();
        
        // Decode and verify
        coderIo->data_len = 0; // Reset data length for decoding
        coder_bwt_cm decoder(coderIo.get());
        std::vector<uint8_t> output(input.size());
        int32_t decoded_len = decoder.decode_line(output.data(), output.size());
        
        EXPECT_EQ(decoded_len, input.size());
        EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), test_data);
    }
}


// Test mixed compression levels in same session
TEST_F(BwtCoderTest, LargeRepeateSingleChar) 
{
    uint32_t length = 2 << 20; 
    std::unique_ptr<uint8_t[]> pOut(new uint8_t[length]);
    std::unique_ptr<coder_io> coderIo(new coder_io(pOut.get(), length));
    coder_bwt_cm encoder(coderIo.get());

    for (int i = 0; i < 2000; i++) {
        int random_num = (rand() % 10) + 1; // Generate random number from 1 to 10
        std::string test_data = std::to_string(random_num) + ":";
        encoder.encode_line(reinterpret_cast<const uint8_t*>(test_data.data()), test_data.length());
    }
    encoder.encode_flush();
}
   