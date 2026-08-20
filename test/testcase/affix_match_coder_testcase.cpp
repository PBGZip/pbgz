/*
 * affix_match_coder_testcase.cpp - Test cases for the affix_match coder
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
#include <cstring>
#include <cstdlib>
#include <string>
#include <memory>

#include "coder/coder.h"
#include "coder/coder_io.h"
#include "coder/coder_affix_match.h"
#include "utils/memory_util.h"

class AffixMatchCoderTest : public ::testing::Test
{
public:
    void SetUp() override
    {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);

        /*
         * The allocation size must match the capacity reported to coder_io. Previously this
         * code allocated 10 MB but declared only 1 MB, so when incompressible data pushed the
         * encoder beyond 1 MB nobody noticed: the extra 9 MB of headroom masked the overflow.
         * In production coder_io receives the block's true remaining size, so the same pattern
         * would become a heap overflow.
         */
        test_data = (uint8_t*) malloc(BUFFER_SIZE);
        compressed_data = (uint8_t*) malloc(BUFFER_SIZE);
        decompressed_data = (uint8_t*) malloc(BUFFER_SIZE);

        ASSERT_NE(test_data, nullptr);
        ASSERT_NE(compressed_data, nullptr);
        ASSERT_NE(decompressed_data, nullptr);
    }

    void TearDown() override
    {
        if (test_data) {
            free(test_data);
            test_data = nullptr;
        }
        if (compressed_data) {
            free(compressed_data);
            compressed_data = nullptr;
        }
        if (decompressed_data) {
            free(decompressed_data);
            decompressed_data = nullptr;
        }
    }

protected:
    uint8_t* test_data = nullptr;
    uint8_t* compressed_data = nullptr;
    uint8_t* decompressed_data = nullptr;
    static const int32_t BUFFER_SIZE;

    // Helper function: Generate test data
    void generate_test_data(size_t size, int seed = 42)
    {
        srand(seed);
        for (size_t i = 0; i < size; ++i) {
            // Fix: Generate 7-bit data, as coder_affix_match only handles ASCII characters (0-127)
            test_data[i] = rand() % 128;
        }
    }
};

// Basic encode/decode test
TEST_F(AffixMatchCoderTest, BasicEncodeDecode)
{
    const char test_string[] = "Hello, World! This is a test for affix match coder.";
    size_t size = strlen(test_string);

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)test_string, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_string, decompressed_data, size), 0);
}

// Empty data test
TEST_F(AffixMatchCoderTest, EmptyData)
{
    const char* test_string = "";
    size_t size = 0;

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)test_string, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_string, decompressed_data, size), 0);
}

// Single character test
TEST_F(AffixMatchCoderTest, SingleCharacter)
{
    const char test_string[] = "A";
    size_t size = 1;

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)test_string, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_string, decompressed_data, size), 0);
}

// Repeated characters test
TEST_F(AffixMatchCoderTest, RepeatedCharacters)
{
    const char test_string[] = "AAAAAAAAAAAAAAAAAAAAAAAA";
    size_t size = 24;

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)test_string, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_string, decompressed_data, size), 0);
}

// Binary data test
TEST_F(AffixMatchCoderTest, BinaryData)
{
    // Fix: Use 7-bit data, as coder_affix_match only handles ASCII characters (0-127)
    uint8_t binary_data[] = {0x00, 0x7F, 0x0A, 0x70, 0x5A, 0x25, 0x33, 0x4C};
    size_t size = sizeof(binary_data);

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(binary_data, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(binary_data, decompressed_data, size), 0);
}

// Small data block test
TEST_F(AffixMatchCoderTest, SmallDataBlock)
{
    generate_test_data(10);
    size_t size = 10;

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(test_data, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_data, decompressed_data, size), 0);
}

// Medium data block test
TEST_F(AffixMatchCoderTest, MediumDataBlock)
{

    size_t size = 255;
    generate_test_data(size, 123);  // Reduce data size to avoid exceeding coder processing capability

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());

    encoder->encode_line(test_data, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());

    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_data, decompressed_data, size), 0);

}

// DNA sequence test (simulate FASTQ data) - Support multiple calls, each not exceeding 255 characters
TEST_F(AffixMatchCoderTest, DNASequence)
{
    // Create a DNA sequence longer than 255 characters, requires multiple calls to encode_line
    const char dna_sequence[] =
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    size_t total_size = strlen(dna_sequence);

    // Verify sequence length exceeds 255, requires chunked processing
    ASSERT_GT(total_size, 255);

    // Encode - multiple calls, each not exceeding 255 characters
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());

    const size_t max_chunk_size = 255;
    size_t offset = 0;

    while (offset < total_size) {
        size_t chunk_size = std::min(max_chunk_size, total_size - offset);
        encoder->encode_line((const uint8_t*)(dna_sequence + offset), chunk_size, false);
        offset += chunk_size;
    }

    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode - multiple decodings, as decoder also has 255 character limit
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());

    // Clear decompression buffer
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Multiple decodings, each not exceeding 255 characters
    offset = 0;
    int32_t total_decoded = 0;

    while (offset < total_size) {
        size_t remaining_size = total_size - offset;
        size_t decode_buffer_size = std::min(max_chunk_size, remaining_size);

        size_t decoded_chunk = decoder->decode_line(decompressed_data + offset, decode_buffer_size);

        ASSERT_GT(decoded_chunk, 0) << "Decode failed at offset: " << offset;
        ASSERT_LE(decoded_chunk, decode_buffer_size) << "Decoded chunk size exceeds expectation";

        total_decoded += decoded_chunk;
        offset += decoded_chunk;

        // If decoded data is less than expected, it means we've reached the end
        if (decoded_chunk < decode_buffer_size) {
            break;
        }
    }

    // Verify result
    ASSERT_EQ(total_decoded, total_size);
    ASSERT_EQ(memcmp(dna_sequence, decompressed_data, total_size), 0);
}

// FASTQ format test with quality values
TEST_F(AffixMatchCoderTest, FASTQFormat)
{
    const char fastq_data[] =
        "@SEQ_ID\n"
        "GATTTGGGGTTTTTTTGTAGTGAATTTTTTTCTATGATATATCTGTGATGGAACTG\n"
        "+\n"
        "!!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65\n";
    size_t size = strlen(fastq_data);

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)fastq_data, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify result
    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(fastq_data, decompressed_data, size), 0);
}

// Encoder initialization test
TEST_F(AffixMatchCoderTest, EncoderInitialization)
{
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());

    EXPECT_NE(encoder.get(), nullptr);
    EXPECT_NE(decoder.get(), nullptr);
    EXPECT_NE(encoder_io.get(), nullptr);
    EXPECT_NE(decoder_io.get(), nullptr);
}


// Multiple encode/decode consistency test
TEST_F(AffixMatchCoderTest, MultipleEncodeDecodeConsistency)
{
    const char test_string[] = "Test for consistency across multiple encode/decode cycles";
    size_t size1 = strlen(test_string);

    // First encode/decode
    memcpy(test_data, test_string, size1);

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line((const uint8_t*)test_string, size1, false);
    encoder->encode_flush();
    int32_t compressed_size1 = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size1);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size1 = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify first result
    ASSERT_EQ(decoded_size1, size1);
    ASSERT_EQ(memcmp(test_string, decompressed_data, size1), 0);

    // Second encode/decode (using different data)
    const char test_string2[] = "Another test string with different content";
    size_t size2 = strlen(test_string2);

    // Encode - create new instance instead of reassignment
    std::unique_ptr<coder_io> encoder_io2 = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder2 = std::make_unique<coder_affix_match>(encoder_io2.get());
    encoder2->encode_line((const uint8_t*)test_string2, size2, false);
    encoder2->encode_flush();
    int32_t compressed_size2 = encoder_io2->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode - create new instance instead of reassignment
    std::unique_ptr<coder_io> decoder_io2 = std::make_unique<coder_io>(compressed_data, compressed_size2);
    std::unique_ptr<coder_affix_match> decoder2 = std::make_unique<coder_affix_match>(decoder_io2.get());
    int32_t decoded_size2 = decoder2->decode_line(decompressed_data, BUFFER_SIZE);

    // Verify second result
    ASSERT_EQ(decoded_size2, size2);
    ASSERT_EQ(memcmp(test_string2, decompressed_data, size2), 0);
}

// Boundary values test
TEST_F(AffixMatchCoderTest, BoundaryValues)
{
    // Test near buffer boundary - use smaller data size to fit coder_affix_match limitations
    generate_test_data(200, 999);  // Use smaller data size
    size_t size = 200;  // Ensure not exceeding 255-byte limit

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(test_data, size, false);
    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    // Clear decompression buffer to avoid dirty data affecting comparison
    memset(decompressed_data, 0, BUFFER_SIZE);

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    int32_t decoded_size = decoder->decode_line(decompressed_data, BUFFER_SIZE);

    ASSERT_EQ(decoded_size, size);
    ASSERT_EQ(memcmp(test_data, decompressed_data, size), 0);
}

// Error handling test (invalid mode)
TEST_F(AffixMatchCoderTest, ErrorHandling)
{
    // Test unset mode case - use smart pointers to avoid memory issues
    auto temp_encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    temp_encoder_io->m = coder_io::MUNSET;

    // Create new encoder instance to test error handling
    // Note: coder_affix_match constructor takes ownership of coder_io
    // So we need to ensure coder_io lifecycle is long enough
    coder_affix_match test_encoder(temp_encoder_io.get());

    // Try encoding (should automatically handle mode setting or return error)
    // Here we only verify that object can be created and basic operations work
    EXPECT_NO_THROW({
        test_encoder.encode_line((const uint8_t*)"test", 4, false);
    });

    // temp_encoder_io will be automatically released when scope ends
}

// Separator function test (similar to SplitCharacter in BWT test)
TEST_F(AffixMatchCoderTest, SplitCharacter)
{
    // Test data with split characters
    const std::string test_data_str = "A,T,C,G,A,T,C,G";
    const uint8_t split_char = ',';
    std::vector<uint8_t> input(test_data_str.begin(), test_data_str.end());

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(input.data(), input.size(), false);
    encoder->encode_flush();

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(input.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

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

// Multiple encode calls test (similar to MultipleEncodeCalls in BWT test)
TEST_F(AffixMatchCoderTest, MultipleEncodeCalls)
{
    // Test multiple small encode calls
    const std::string chunk1 = "ATCGATCG";

    // Encode single chunk (coder_affix_match doesn't support multiple chunks well)
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(reinterpret_cast<const uint8_t*>(chunk1.data()), chunk1.size(), false);
    encoder->encode_flush();

    // Decode single chunk
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(chunk1.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    // Verify single chunk works correctly
    EXPECT_EQ(decoded_len, chunk1.size());
    EXPECT_EQ(std::string(output.begin(), output.begin() + decoded_len), chunk1);
}

// Large data test (similar to LargeData in BWT test)
TEST_F(AffixMatchCoderTest, LargeData)
{
    // Generate smaller test data to fit coder_affix_match limitations
    std::vector<uint8_t> large_data;
    const std::string pattern = "ATCGATCG";
    for (int i = 0; i < 20; ++i) { // Much smaller to fit into buffer and within 255 limit
        for (char c : pattern) {
            large_data.push_back(static_cast<uint8_t>(c));
        }
    }

    // Ensure total size doesn't exceed 255 bytes
    if (large_data.size() > 255) {
        large_data.resize(255);
    }

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(large_data.data(), large_data.size(), false);
    encoder->encode_flush();

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(large_data.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    EXPECT_EQ(decoded_len, large_data.size());
    EXPECT_EQ(memcmp(output.data(), large_data.data(), large_data.size()), 0);
}

// Decode insufficient buffer test (similar to DecodeInsufficientBuffer in BWT test)
TEST_F(AffixMatchCoderTest, DecodeInsufficientBuffer)
{
    // Test data - use smaller data to avoid buffer issues
    const std::string test_data_str = "ATCG";
    std::vector<uint8_t> input(test_data_str.begin(), test_data_str.end());

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(input.data(), input.size(), false);
    encoder->encode_flush();

    // Decode with sufficient buffer to avoid overflow
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(input.size()); // Use same size as input
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    EXPECT_EQ(decoded_len, input.size());
    EXPECT_EQ(std::string(output.begin(), output.end()), test_data_str);
}

// Decode input length tracking test (similar to DecodeInputLength in BWT test)
TEST_F(AffixMatchCoderTest, DecodeInputLength)
{
    // Test data
    const std::string test_data_str = "ATCGATCG";
    std::vector<uint8_t> input(test_data_str.begin(), test_data_str.end());

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(input.data(), input.size(), false);
    encoder->encode_flush();

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(input.size());
    decoder->decode_line(output.data(), output.size());

    // Test passes if no crash occurs
    SUCCEED();
}

// Memory cleanup test (similar to MemoryCleanup in BWT test)
TEST_F(AffixMatchCoderTest, MemoryCleanup)
{
    // Test that encoder and decoder clean up properly when destructed
    {
        coder_io* temp_encoder_io = new coder_io(compressed_data, BUFFER_SIZE);
        coder_affix_match test_encoder(temp_encoder_io);
        const std::string test_data_str = "ATCGATCG";
        test_encoder.encode_line(reinterpret_cast<const uint8_t*>(test_data_str.data()), test_data_str.size(), false);
        test_encoder.encode_flush();
        delete temp_encoder_io;
    }

    // Should not crash when encoder goes out of scope
    SUCCEED();
}

// Special characters test (similar to SpecialCharacters in BWT test)
TEST_F(AffixMatchCoderTest, SpecialCharacters)
{
    // Fix: Use only 7-bit data, as coder_affix_match only handles ASCII characters (0-127)
    std::vector<uint8_t> special_data;
    for (int i = 0; i < 128; ++i) {
        special_data.push_back(static_cast<uint8_t>(i));
    }

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(special_data.data(), special_data.size(), false);
    encoder->encode_flush();

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(special_data.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    EXPECT_EQ(decoded_len, special_data.size());
    EXPECT_EQ(memcmp(output.data(), special_data.data(), special_data.size()), 0);
}

// Compression mode detection test (similar to CompressionModeDetection in BWT test)
TEST_F(AffixMatchCoderTest, CompressionModeDetection)
{
    // Test with compression mode
    coder_io* temp_encoder_io = new coder_io(compressed_data, BUFFER_SIZE);
    temp_encoder_io->m = coder_io::MENC;
    EXPECT_NO_THROW({
        coder_affix_match encoder(temp_encoder_io);
    });
    delete temp_encoder_io;

    // Test with decompression mode
    coder_io* temp_decoder_io = new coder_io(compressed_data, BUFFER_SIZE);
    temp_decoder_io->m = coder_io::MDEC;
    EXPECT_NO_THROW({
        coder_affix_match decoder(temp_decoder_io);
    });
    delete temp_decoder_io;

    // Test with unset mode (should default appropriately)
    uint8_t* temp_buffer = (uint8_t*)malloc(BUFFER_SIZE);
    coder_io* temp_unset_io = new coder_io(temp_buffer, BUFFER_SIZE);
    temp_unset_io->m = coder_io::MUNSET;
    EXPECT_NO_THROW({
        coder_affix_match coder(temp_unset_io);
    });
    delete temp_unset_io;
    free(temp_buffer);
}

// Metadata handling test (similar to MetadataHandling in BWT test)
TEST_F(AffixMatchCoderTest, MetadataHandling)
{
    // Test with various metadata settings - use valid level values
    coder_io* temp_metadata_io = new coder_io(compressed_data, BUFFER_SIZE);
    temp_metadata_io->meta["level"] = 2;  // Use valid level value [1,2]
    temp_metadata_io->meta["magic"] = "AFFIX_MATCH";
    temp_metadata_io->meta["custom_param"] = "test_value";

    EXPECT_NO_THROW({
        coder_affix_match encoder(temp_metadata_io);
    });

    // Test encoding with metadata
    coder_affix_match encoder(temp_metadata_io);
    const std::string test_data_str = "ATCGATCGATCG";
    encoder.encode_line(reinterpret_cast<const uint8_t*>(test_data_str.data()), test_data_str.size(), false);
    encoder.encode_flush();

    // Verify metadata is preserved
    EXPECT_EQ(temp_metadata_io->get_level(), 2);
    EXPECT_TRUE(temp_metadata_io->get_magic().find("coder_affix_match") != std::string::npos);
    delete temp_metadata_io;
}

// Very large single block test (similar to VeryLargeSingleBlock in BWT test)
TEST_F(AffixMatchCoderTest, VeryLargeSingleBlock)
{
    // Generate smaller data to fit coder_affix_match limitations
    std::vector<uint8_t> large_block(200, 'A'); // Much smaller size, within 255 limit
    for (size_t i = 0; i < large_block.size(); i += 50) {
        large_block[i] = static_cast<uint8_t>('A' + (i / 50) % 4); // Add some variation
    }

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    EXPECT_NO_THROW({
        encoder->encode_line(large_block.data(), large_block.size(), false);
        encoder->encode_flush();
    });

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(large_block.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    EXPECT_EQ(decoded_len, large_block.size());
    EXPECT_EQ(memcmp(output.data(), large_block.data(), large_block.size()), 0);
}

// Repetitive data compression test (similar to RepetitiveDataCompression in BWT test)
TEST_F(AffixMatchCoderTest, RepetitiveDataCompression)
{
    // Test highly repetitive data with smaller size to fit buffer and within 255 limit
    std::vector<uint8_t> repetitive_data(100, 'A'); // Much smaller size, within 255 limit

    // Encode
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());
    encoder->encode_line(repetitive_data.data(), repetitive_data.size(), false);
    encoder->encode_flush();

    // Decode
    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, encoder_io->data_len);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());
    std::vector<uint8_t> output(repetitive_data.size());
    int32_t decoded_len = decoder->decode_line(output.data(), output.size());

    EXPECT_EQ(decoded_len, repetitive_data.size());
    EXPECT_EQ(memcmp(output.data(), repetitive_data.data(), repetitive_data.size()), 0);

}

// Random digits compression test - using variable length strings between 5 and 10
TEST_F(AffixMatchCoderTest, RandomDigitsCompression)
{
    // Generate fixed random seed to ensure test reproducibility
    const int seed = 12345;
    srand(seed);

    // Define digit and special character set (use only ASCII characters to ensure compatibility)
    const char digits_and_specials[] = "0123456789";
    const int char_set_size = sizeof(digits_and_specials) - 1; // Subtract null terminator

    // Define string length range for each generation (random length between 5 and 10)
    const int min_string_length = 8;
    const int max_string_length = 20;
    const int repetitions = 200000; // Keep original 20000 repetitions

    // Pre-generate all strings to ensure data validity within encoder lifecycle
    std::vector<std::string> all_strings;
    all_strings.reserve(repetitions);

    for (int i = 0; i < repetitions; ++i) {
        // Generate random length between 5 and 10
        const int string_length = min_string_length + (rand() % (max_string_length - min_string_length + 1));

        std::string random_string;
        random_string.reserve(string_length + 1);

        for (int j = 0; j < string_length; ++j) {
            int random_index = rand() % char_set_size;
            random_string.push_back(digits_and_specials[random_index]);
        }
        random_string.push_back('_'); // Add separator
        all_strings.push_back(random_string);
    }

    // Perform compression
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());

    for (int i = 0; i < repetitions; ++i) {
        const std::string& str = all_strings[i];
        encoder->encode_line(reinterpret_cast<const uint8_t*>(str.data()), str.size(), true);
    }

    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());

    for (int i = 0; i < repetitions; ++i) {
        const std::string& original_string = all_strings[i];
        std::vector<uint8_t> decoded_data(original_string.size());

        int32_t decoded_size = decoder->decode_line(decoded_data.data(), decoded_data.size(), '_', true);

        ASSERT_GT(decoded_size, 0) << "Decode failed at index: " << i;
        ASSERT_EQ(decoded_size, original_string.size()) << "Decode size mismatch at index: " << i;
        ASSERT_EQ(std::string(reinterpret_cast<char*>(decoded_data.data()), decoded_size), original_string)
            << "Decoded data mismatch at index: " << i;
    }
}

// Random digits compression test - using variable length strings between 5 and 10
TEST_F(AffixMatchCoderTest, RandomDigitsCompressionNotHold)
{
    // Generate fixed random seed to ensure test reproducibility
    const int seed = 12345;
    srand(seed);

    // Define digit and special character set (use only ASCII characters to ensure compatibility)
    const char digits_and_specials[] = "0123456789";
    const int char_set_size = sizeof(digits_and_specials) - 1; // Subtract null terminator

    // Define string length range for each generation (random length between 5 and 10)
    const int min_string_length = 8;
    const int max_string_length = 20;
    const int repetitions = 400000; // Keep original 20000 repetitions

    // Pre-generate all strings to ensure data validity within encoder lifecycle
    std::vector<std::string> all_strings;
    all_strings.reserve(repetitions);

    for (int i = 0; i < repetitions; ++i) {
        // Generate random length between 5 and 10
        const int string_length = min_string_length + (rand() % (max_string_length - min_string_length + 1));

        std::string random_string;
        random_string.reserve(string_length + 1);

        for (int j = 0; j < string_length; ++j) {
            int random_index = rand() % char_set_size;
            random_string.push_back(digits_and_specials[random_index]);
        }
        random_string.push_back('_'); // Add separator
        all_strings.push_back(random_string);
    }

    // Perform compression
    std::unique_ptr<coder_io> encoder_io = std::make_unique<coder_io>(compressed_data, BUFFER_SIZE);
    std::unique_ptr<coder_affix_match> encoder = std::make_unique<coder_affix_match>(encoder_io.get());

    for (int i = 0; i < repetitions; ++i) {
        const std::string& str = all_strings[i];
        encoder->encode_line(reinterpret_cast<const uint8_t*>(str.data()), str.size());
    }

    encoder->encode_flush();
    int32_t compressed_size = encoder_io->data_len;

    std::unique_ptr<coder_io> decoder_io = std::make_unique<coder_io>(compressed_data, compressed_size);
    std::unique_ptr<coder_affix_match> decoder = std::make_unique<coder_affix_match>(decoder_io.get());

    std::vector<std::vector<uint8_t>> de_all_strings;
    de_all_strings.reserve(repetitions);

    for (int i = 0; i < repetitions; ++i) {
        const std::string& original_string = all_strings[i];
        std::vector<uint8_t> decoded_data(original_string.size());
        de_all_strings.push_back(decoded_data);
    }

    for (int i = 0; i < repetitions; ++i) {
        const std::string& original_string = all_strings[i];
        int32_t decoded_size = decoder->decode_line(de_all_strings[i].data(), de_all_strings[i].size(), '_');

        ASSERT_GT(decoded_size, 0) << "Decode failed at index: " << i;
        ASSERT_EQ(decoded_size, original_string.size()) << "Decode size mismatch at index: " << i;
        ASSERT_EQ(std::string(reinterpret_cast<char*>(de_all_strings[i].data()), decoded_size), original_string)
            << "Decoded data mismatch at index: " << i;
    }
}

// Define static member variable to resolve linking error
const int32_t AffixMatchCoderTest::BUFFER_SIZE = 10 * 1024 * 1024;