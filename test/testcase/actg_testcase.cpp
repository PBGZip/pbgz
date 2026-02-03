#include <gtest/gtest.h>
#include <memory>
#include <cstring>
#include <algorithm>
#include <random>
#include <chrono>

// Include the header with the functions to test
#include "actg.h"

// Note: actgStretch and actgStretch2bits constants are defined in reference.cpp
// For testing purposes, we'll use direct computation instead of these lookup tables

class ActgTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test data
        std::random_device rd;
        gen = std::mt19937(rd());
        dis = std::uniform_int_distribution<>(0, 3);
    }

    void TearDown() override {
        // Cleanup if needed
    }

    // Helper function to generate random ACTG sequence
    std::string generateRandomACTG(int length) {
        std::string seq;
        const char bases[] = {'A', 'C', 'T', 'G'};
        for (int i = 0; i < length; ++i) {
            seq += bases[dis(gen)];
        }
        return seq;
    }

    // Helper function to convert ACTG char to 2-bit encoding
    uint8_t actgTo2Bit(char c) {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'T': return 2;
            case 'G': return 3;
            default: return 0;
        }
    }

    // Helper function to convert 2-bit encoding to ACTG char
    char bit2ToActg(uint8_t b) {
        switch (b & 0x3) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'T';
            case 3: return 'G';
            default: return 'A';
        }
    }

    std::mt19937 gen;
    std::uniform_int_distribution<> dis;
};

// Test actgSquash function
TEST_F(ActgTest, ActgSquash) {
    // Test case 1: Perfect 4-byte alignment
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'};
        uint8_t output[2];
        int64_t result = actgSquash(input, 8, output);
        
        EXPECT_EQ(result, 2);
        // ACTG -> 00011011, so "ACTGACTG" should become 0x1B1B
        EXPECT_EQ(output[0], 0x1B);
        EXPECT_EQ(output[1], 0x1B);
    }

    // Test case 2: Non-4-byte alignment (1 extra byte)
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A'};
        uint8_t output[2];
        int64_t result = actgSquash(input, 5, output);
        
        EXPECT_EQ(result, 2);
        EXPECT_EQ(output[0], 0x1B);
        // "A" should be in the high 2 bits: 00xxxxxx = 0x00
        EXPECT_EQ(output[1], 0x00);
    }

    // Test case 3: Non-4-byte alignment (2 extra bytes)
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A', 'C'};
        uint8_t output[2];
        int64_t result = actgSquash(input, 6, output);
        
        EXPECT_EQ(result, 2);
        EXPECT_EQ(output[0], 0x1B);
        // "AC" -> 00010000 = 0x10
        EXPECT_EQ(output[1], 0x10);
    }

    // Test case 4: Non-4-byte alignment (3 extra bytes)
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T'};
        uint8_t output[2];
        int64_t result = actgSquash(input, 7, output);
        
        EXPECT_EQ(result, 2);
        EXPECT_EQ(output[0], 0x1B);
        // "ACT" -> 00011000 = 0x18
        EXPECT_EQ(output[1], 0x18);
    }

    // Test case 5: Empty input
    {
        uint8_t output[1];
        int64_t result = actgSquash(nullptr, 0, output);
        EXPECT_EQ(result, 0);
    }

    // Test case 6: Large random data
    {
        std::string seq = generateRandomACTG(1000);
        uint8_t output[250];
        int64_t result = actgSquash((const uint8_t*)seq.c_str(), 1000, output);
        EXPECT_EQ(result, 250);
        
        // Verify first few bytes manually
        for (int i = 0; i < 4; ++i) {
            uint8_t expected = 0;
            for (int j = 0; j < 4; ++j) {
                expected |= (actgTo2Bit(seq[i*4 + j]) << (6 - j*2));
            }
            EXPECT_EQ(output[i], expected);
        }
    }

    {
        const uint8_t input[] = {'N', 'N', 'N', 'N', 'N', 'N', 'N'};
        uint8_t output[2];
        int64_t result = actgSquash(input, 7, output);
        
        EXPECT_EQ(result, 2);
        EXPECT_EQ(output[0], 0xFF);
        // NNN -> 0x11111100
        EXPECT_EQ(output[1], 0xFC);
    }
}

// Test actgPair function
TEST_F(ActgTest, ActgPair) {
    // Test case 1: Basic complement (reverse order due to reverse writing)
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G'};
        uint8_t output[4];
        actgPair(output, input, 4);
        
        // Due to reverse writing: output[0] gets complement of input[3], etc.
        // input[3]='G'→'C', input[2]='T'→'A', input[1]='C'→'G', input[0]='A'→'T'
        EXPECT_EQ(output[0], 'C');  // G's complement
        EXPECT_EQ(output[1], 'A');  // T's complement
        EXPECT_EQ(output[2], 'G');  // C's complement
        EXPECT_EQ(output[3], 'T');  // A's complement
    }

    // Test case 2: Longer sequence (reverse order due to reverse writing)
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'};
        uint8_t output[8];
        actgPair(output, input, 8);
        
        // Expected output in reverse order:
        // input[7]='G'→'C', input[6]='T'→'A', input[5]='C'→'G', input[4]='A'→'T'
        // input[3]='G'→'C', input[2]='T'→'A', input[1]='C'→'G', input[0]='A'→'T'
        const uint8_t expected[] = {'C', 'A', 'G', 'T', 'C', 'A', 'G', 'T'};
        for (int i = 0; i < 8; ++i) {
            EXPECT_EQ(output[i], expected[i]);
        }
    }

    // Test case 3: Random sequence (verify reverse order complement)
    {
        std::string seq = generateRandomACTG(100);
        uint8_t input[100];
        uint8_t output[100];
        
        for (int i = 0; i < 100; ++i) {
            input[i] = seq[i];
        }
        
        actgPair(output, input, 100);
        
        // Verify complement relationship in reverse order
        // output[i] should be complement of input[99-i]
        for (int i = 0; i < 100; ++i) {
            char expected;
            switch (seq[99 - i]) {  // Reverse index due to reverse writing
                case 'A': expected = 'T'; break;
                case 'C': expected = 'G'; break;
                case 'T': expected = 'A'; break;
                case 'G': expected = 'C'; break;
                default: expected = seq[99 - i]; break;
            }
            EXPECT_EQ(output[i], expected);
        }
    }

    // Test case 4: Empty input
    {
        uint8_t output[1];
        actgPair(output, nullptr, 0);
        // Should not crash
    }
}

// Test getDiffCnt function
TEST_F(ActgTest, GetDiffCnt) {
    // Test case 1: Identical strings
    {
        const uint8_t s1[] = {'A', 'C', 'T', 'G'};
        const uint8_t s2[] = {'A', 'C', 'T', 'G'};
        uint32_t result = getDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 0);
    }

    // Test case 2: Completely different strings
    {
        const uint8_t s1[] = {'A', 'C', 'T', 'G'};
        const uint8_t s2[] = {'T', 'G', 'A', 'C'};
        uint32_t result = getDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 4);
    }

    // Test case 3: Partially different
    {
        const uint8_t s1[] = {'A', 'C', 'T', 'G'};
        const uint8_t s2[] = {'A', 'G', 'T', 'C'};
        uint32_t result = getDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 2);
    }

    // Test case 4: Empty strings
    {
        uint32_t result = getDiffCnt(nullptr, nullptr, 0);
        EXPECT_EQ(result, 0);
    }
}

// Test actgSquashDiffCnt function
TEST_F(ActgTest, ActgSquashDiffCnt) {
    // Test case 1: Safety checks
    {
        uint32_t result = actgSquashDiffCnt(nullptr, nullptr, 10);
        EXPECT_EQ(result, 0);
        
        uint8_t dummy[1] = {0};
        result = actgSquashDiffCnt(dummy, nullptr, 10);
        EXPECT_EQ(result, 0);
        
        result = actgSquashDiffCnt(nullptr, dummy, 10);
        EXPECT_EQ(result, 0);
        
        result = actgSquashDiffCnt(dummy, dummy, 0);
        EXPECT_EQ(result, 0);
    }

    // Test case 2: Identical squash data
    {
        const uint8_t s1[] = {0x1B, 0x1B, 0x1B, 0x1B};
        const uint8_t s2[] = {0x1B, 0x1B, 0x1B, 0x1B};
        uint32_t result = actgSquashDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 0);
    }

    // Test case 3: Different squash data
    {
        const uint8_t s1[] = {0x1B, 0x1B, 0x1B, 0x1B}; // ACTGACTGACTGACTG
        const uint8_t s2[] = {0xE4, 0xE4, 0xE4, 0xE4}; // TGACTGACTGACTGAC
        uint32_t result = actgSquashDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 16); // All 4 bytes * 4 bases each = 16 differences
    }

    // Test case 4: Partial differences
    {
        const uint8_t s1[] = {0x1B, 0x1B, 0x1B, 0x1B};
        const uint8_t s2[] = {0x1B, 0xE4, 0x1B, 0xE4};
        uint32_t result = actgSquashDiffCnt(s1, s2, 4);
        EXPECT_EQ(result, 8); // 2 bytes * 4 bases each = 8 differences
    }

    // Test case 5: Large data test
    {
        const int size = 1000;
        uint8_t s1[size];
        uint8_t s2[size];
        
        for (int i = 0; i < size; ++i) {
            s1[i] = i % 256;
            s2[i] = (i + 1) % 256;
        }
        
        uint32_t result = actgSquashDiffCnt(s1, s2, size);
        EXPECT_GT(result, 0);
        EXPECT_LE(result, size * 4);
    }
}

// Test actgStretchMapping function
TEST_F(ActgTest, ActgStretchMapping) {
    // Test case 1: Identical sequences (should produce all zeros)
    {
        const uint8_t base[] = {0x1B, 0x1B};
        const uint8_t refe[] = {0x1B, 0x1B};
        uint8_t output[16];
        uint32_t result = actgStretchMapping(base, refe, 2, output);
        
        // 2 bytes input → 8 bytes output (function processes 2-byte blocks)
        EXPECT_EQ(result, 8);
        for (int i = 0; i < 8; ++i) {
            EXPECT_EQ(output[i], 0);
        }
    }

    // Test case 2: Different sequences
    {
        const uint8_t base[] = {0x1B, 0xE4}; // ACTG, TGAC
        const uint8_t refe[] = {0xE4, 0x1B}; // TGAC, ACTG
        uint8_t output[16];
        uint32_t result = actgStretchMapping(base, refe, 2, output);
        
        // 2 bytes input → 8 bytes output
        EXPECT_EQ(result, 8);
        // Should have non-zero values for differences
        bool hasNonZero = false;
        for (int i = 0; i < 8; ++i) {
            if (output[i] != 0) {
                hasNonZero = true;
                break;
            }
        }
        EXPECT_TRUE(hasNonZero);
    }

    // Test case 3: Single byte difference
    {
        const uint8_t base[] = {0x1B, 0x1B};
        const uint8_t refe[] = {0x1B, 0xE4};
        uint8_t output[16];
        uint32_t result = actgStretchMapping(base, refe, 2, output);
        
        // 2 bytes input → 8 bytes output
        EXPECT_EQ(result, 8);
        // All 8 bytes should be processed (function processes 2-byte blocks as a unit)
        bool hasNonZero = false;
        for (int i = 0; i < 8; ++i) {
            if (output[i] != 0) {
                hasNonZero = true;
                break;
            }
        }
        EXPECT_TRUE(hasNonZero);
    }

    // Test case 4: Odd length
    {
        const uint8_t base[] = {0x1B, 0x1B, 0x1B};
        const uint8_t refe[] = {0x1B, 0xE4, 0x1B};
        uint8_t output[16];
        uint32_t result = actgStretchMapping(base, refe, 3, output);
        
        // 3 bytes input: first 2 bytes processed as one block (8 bytes), 
        // remaining 1 byte processed separately (4 bytes)
        EXPECT_GT(result, 0);
        EXPECT_LE(result, 12); // 8 + 4 = 12 bytes maximum
    }
}

// Test actgStretchMappingXor function
TEST_F(ActgTest, ActgStretchMappingXor) {
    // Test case 1: Identical sequences (should produce all zeros)
    {
        const uint8_t base[] = {0x1B, 0x1B};
        const uint8_t refe[] = {0x1B, 0x1B};
        uint8_t output[16];
        uint32_t result = actgStretchMappingXor(base, refe, 2, output);
        
        // 2 bytes input → 8 bytes output (function processes 2-byte blocks)
        EXPECT_EQ(result, 8);
        for (int i = 0; i < 8; ++i) {
            EXPECT_EQ(output[i], 0);
        }
    }

    // Test case 2: Different sequences
    {
        const uint8_t base[] = {0x1B, 0xE4};
        const uint8_t refe[] = {0xE4, 0x1B};
        uint8_t output[16];
        uint32_t result = actgStretchMappingXor(base, refe, 2, output);
        
        // 2 bytes input → 8 bytes output
        EXPECT_EQ(result, 8);
        // Should have non-zero values
        bool hasNonZero = false;
        for (int i = 0; i < 8; ++i) {
            if (output[i] != 0) {
                hasNonZero = true;
                break;
            }
        }
        EXPECT_TRUE(hasNonZero);
    }

    // Test case 3: Large data test
    {
        const int size = 100;
        uint8_t base[size];
        uint8_t refe[size];
        uint8_t output[400];  // 100 bytes → 96/8=12 blocks × 32 + 4 bytes = 400 bytes max
        
        for (int i = 0; i < size; ++i) {
            base[i] = i % 256;
            refe[i] = (i + 128) % 256;
        }
        
        uint32_t result = actgStretchMappingXor(base, refe, size, output);
        // 100 bytes: 96 bytes (12 blocks of 8) × 32 + 4 bytes (2 blocks of 2) = 400 bytes
        EXPECT_EQ(result, 400);
    }

    // Test case 4: Various lengths
    {
        for (int len = 1; len <= 10; ++len) {
            uint8_t base[len];
            uint8_t refe[len];
            uint8_t output[40];  // Enough for 10 bytes max
            
            for (int i = 0; i < len; ++i) {
                base[i] = i % 256;
                refe[i] = (i + 1) % 256;
            }
            
            uint32_t result = actgStretchMappingXor(base, refe, len, output);
            EXPECT_GT(result, 0);
            
            // Calculate expected output based on function logic
            uint32_t expected_max;
            if (len <= 7) {
                expected_max = len * 4;  // Each byte processed separately
            } else if (len == 8) {
                expected_max = 32;  // One 8-byte block
            } else {
                expected_max = (len / 8) * 32 + (len % 8) * 4;
            }
            EXPECT_EQ(result, expected_max);
        }
    }
}

// Test actgEncode function
TEST_F(ActgTest, ActgEncode) {
    // Test case 1: Basic encoding
    {
        const uint8_t input[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'};
        uint8_t output[8];
        actgEncode(input, output, 8);
        
        // A=0, C=1, T=2, G=3, so should be 0x01230123
        EXPECT_EQ(output[0], 0);
        EXPECT_EQ(output[1], 1);
        EXPECT_EQ(output[2], 2);
        EXPECT_EQ(output[3], 3);
        EXPECT_EQ(output[4], 0);
        EXPECT_EQ(output[5], 1);
        EXPECT_EQ(output[6], 2);
        EXPECT_EQ(output[7], 3);
    }

    // Test case 2: Large data
    {
        const int size = 1000;
        uint8_t input[size];
        uint8_t output[size];
        
        for (int i = 0; i < size; ++i) {
            input[i] = "ACTG"[i % 4];
        }
        
        actgEncode(input, output, size);
        
        for (int i = 0; i < size; ++i) {
            EXPECT_EQ(output[i], i % 4);
        }
    }

    // Test case 3: Empty input
    {
        uint8_t output[1];
        actgEncode(nullptr, output, 0);
        // Should not crash
    }
}

// Test actgXor function
TEST_F(ActgTest, ActgXor) {
    // Test case 1: Basic XOR
    {
        const uint8_t x1[] = {0x12, 0x34, 0x56, 0x78};
        const uint8_t x2[] = {0xAB, 0xCD, 0xEF, 0x90};
        uint8_t output[4];
        actgXor(x1, x2, output, 4);
        
        EXPECT_EQ(output[0], 0x12 ^ 0xAB);
        EXPECT_EQ(output[1], 0x34 ^ 0xCD);
        EXPECT_EQ(output[2], 0x56 ^ 0xEF);
        EXPECT_EQ(output[3], 0x78 ^ 0x90);
    }

    // Test case 2: Identical inputs (should be all zeros)
    {
        const uint8_t x1[] = {0x12, 0x34, 0x56, 0x78};
        uint8_t output[4];
        actgXor(x1, x1, output, 4);
        
        for (int i = 0; i < 4; ++i) {
            EXPECT_EQ(output[i], 0);
        }
    }

    // Test case 3: Various lengths
    {
        for (int len = 1; len <= 20; ++len) {
            uint8_t x1[len];
            uint8_t x2[len];
            uint8_t output[len];
            
            for (int i = 0; i < len; ++i) {
                x1[i] = i % 256;
                x2[i] = (i * 2) % 256;
            }
            
            actgXor(x1, x2, output, len);
            
            for (int i = 0; i < len; ++i) {
                EXPECT_EQ(output[i], x1[i] ^ x2[i]);
            }
        }
    }

    // Test case 4: Large data test
    {
        const int size = 1000;
        uint8_t x1[size];
        uint8_t x2[size];
        uint8_t output[size];
        
        for (int i = 0; i < size; ++i) {
            x1[i] = i % 256;
            x2[i] = (i + 100) % 256;
        }
        
        actgXor(x1, x2, output, size);
        
        for (int i = 0; i < size; ++i) {
            EXPECT_EQ(output[i], x1[i] ^ x2[i]);
        }
    }
}

// Test integration scenarios
TEST_F(ActgTest, IntegrationScenarios) {
    // Test complement and mapping
    {
        const uint8_t original[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'};
        uint8_t complement[8];
        uint8_t squash[2];
        uint8_t squashComp[2];
        uint8_t mapping[16];
        
        // Get complement
        actgPair(complement, original, 8);
        
        // Squash both
        actgSquash(original, 8, squash);
        actgSquash(complement, 8, squashComp);
        
        // They should be different
        bool different = false;
        for (int i = 0; i < 2; ++i) {
            if (squash[i] != squashComp[i]) {
                different = true;
                break;
            }
        }
        EXPECT_TRUE(different);
        
        // Test mapping
        uint32_t mapLen = actgStretchMapping(squash, squashComp, 2, mapping);
        EXPECT_EQ(mapLen, 8);
        
        // Should have differences
        bool hasDiff = false;
        for (int i = 0; i < 8; ++i) {
            if (mapping[i] != 0) {
                hasDiff = true;
                break;
            }
        }
        EXPECT_TRUE(hasDiff);
    }
    
    // Test squash round-trip verification
    {
        const uint8_t original[] = {'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G'};
        uint8_t squash[2];
        
        // Convert to squash
        int64_t squashLen = actgSquash(original, 8, squash);
        EXPECT_EQ(squashLen, 2);
        EXPECT_EQ(squash[0], 0x1B); // ACTG -> 00011011
        EXPECT_EQ(squash[1], 0x1B); // ACTG -> 00011011
    }
}

// Performance test
TEST_F(ActgTest, PerformanceTest) {
    const int size = 10000;
    uint8_t data[size];
    uint8_t output[size];
    
    // Generate test data
    for (int i = 0; i < size; ++i) {
        data[i] = "ACTG"[i % 4];
    }
    
    // Test actgSquash performance
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; ++i) {
        actgSquash(data, size, output);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Should complete reasonably fast (adjust threshold as needed)
    EXPECT_LT(duration.count(), 100000); // Less than 100ms for 1000 iterations
    
    // Test actgPair performance
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; ++i) {
        actgPair(output, data, size);
    }
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    EXPECT_LT(duration.count(), 100000); // Less than 100ms for 1000 iterations
}

// Edge cases and error handling
TEST_F(ActgTest, EdgeCases) {
    // Test with null pointers (should not crash)
    {
        uint8_t dummy[1] = {0};
        EXPECT_NO_THROW(actgSquash(nullptr, 0, dummy));
        EXPECT_NO_THROW(actgPair(dummy, nullptr, 0));
        EXPECT_EQ(getDiffCnt(nullptr, nullptr, 0), 0);
        EXPECT_EQ(actgSquashDiffCnt(nullptr, nullptr, 0), 0);
        EXPECT_NO_THROW(actgEncode(nullptr, dummy, 0));
        EXPECT_NO_THROW(actgXor(nullptr, nullptr, dummy, 0));
    }
    
    // Test with maximum values
    {
        const int maxSize = 100000;
        uint8_t* largeInput = new uint8_t[maxSize];
        uint8_t* largeOutput = new uint8_t[maxSize];
        
        for (int i = 0; i < maxSize; ++i) {
            largeInput[i] = i % 256;
        }
        
        EXPECT_NO_THROW(actgSquash(largeInput, maxSize, largeOutput));
        EXPECT_NO_THROW(actgPair(largeOutput, largeInput, maxSize));
        EXPECT_NO_THROW(actgXor(largeInput, largeOutput, largeOutput, maxSize));
        
        delete[] largeInput;
        delete[] largeOutput;
    }
}
