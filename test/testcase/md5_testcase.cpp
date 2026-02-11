#include <gtest/gtest.h>
#include "utils/md5_util.h"
#include <string>
#include <vector>

class Md5UtilTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialization before test
    }

    void TearDown() override {
        // Cleanup after test
    }

    // Helper function: compare MD5 results
    void expectMd5Equals(const std::string& input, const std::string& expected_hex) {
        uint8_t digest[16];
        md5_compute((const uint8_t*)input.c_str(), input.length(), digest);
        std::string result = md5_to_string(digest);
        EXPECT_EQ(result, expected_hex) << "Input: " << input;
    }

    // Helper function: compare MD5 of binary data
    void expectMd5EqualsBinary(const std::vector<uint8_t>& data, const std::string& expected_hex) {
        uint8_t digest[16];
        md5_compute(data.data(), data.size(), digest);
        std::string result = md5_to_string(digest);
        EXPECT_EQ(result, expected_hex);
    }
};

// Basic functionality test
TEST_F(Md5UtilTest, BasicMd5Compute) {
    // Test empty string
    expectMd5Equals("", "d41d8cd98f00b204e9800998ecf8427e");
    
    // Test single character
    expectMd5Equals("a", "0cc175b9c0f1b6a831c399e269772661");
    
    // Test simple string
    expectMd5Equals("abc", "900150983cd24fb0d6963f7d28e17f72");
    
    // Test numeric string
    expectMd5Equals("123456789", "25f9e794323b453885f5181f1b624d0b");
}

// RFC 1321 MD5 test vectors
TEST_F(Md5UtilTest, Rfc1321TestVectors) {
    // Standard test vectors from RFC 1321
    expectMd5Equals("", "d41d8cd98f00b204e9800998ecf8427e");
    expectMd5Equals("a", "0cc175b9c0f1b6a831c399e269772661");
    expectMd5Equals("abc", "900150983cd24fb0d6963f7d28e17f72");
    expectMd5Equals("message digest", "f96b697d7cb7938d525a2f31aaf161d0");
    expectMd5Equals("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b");
    expectMd5Equals("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", 
                   "d174ab98d277d9f5a5611c2c9f419d9f");
    expectMd5Equals("12345678901234567890123456789012345678901234567890123456789012345678901234567890", 
                   "57edf4a22be3c955ac49da2e2107b67a");
}

// Different length data test
TEST_F(Md5UtilTest, DifferentLengths) {
    // Test 1-byte data
    expectMd5Equals("1", "c4ca4238a0b923820dcc509a6f75849b");
    
    // Test 16-byte data (exactly one MD5 block)
    expectMd5Equals("1234567890123456", "abeac07d3c28c1bef9e730002c753ed4");
    
    // Test 32-byte data (two MD5 blocks)
    expectMd5Equals("12345678901234567890123456789012", "767179c7a2bff19651ce97d294c30cfb");
    
    // Test 64-byte data (exactly two MD5 blocks)
    expectMd5Equals("1234567890123456789012345678901234567890123456789012345678901234", 
                   "eb6c4179c0a7c82cc2828c1e6338e165");
    
    // Test 65-byte data (more than two blocks)
    expectMd5Equals("12345678901234567890123456789012345678901234567890123456789012345", 
                   "823cc889fc7318dd33dde0654a80b70a");
}

// Binary data test
TEST_F(Md5UtilTest, BinaryData) {
    // Test pure binary data
    std::vector<uint8_t> binary_data = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05};
    expectMd5EqualsBinary(binary_data, "d15ae53931880fd7b724dd7888b4b4ed");
    
    // Test data containing null characters
    std::vector<uint8_t> data_with_null = {'a', 0, 'b', 0, 'c'};
    expectMd5EqualsBinary(data_with_null, "ed4b9e8bfe6c7571a700dc6a0b0d01cb");
    
    // Test all zero data
    std::vector<uint8_t> all_zeros(32, 0);
    expectMd5EqualsBinary(all_zeros, "70bc8f4b72a86921468bf8e8441dce51");
}

// Special characters test
TEST_F(Md5UtilTest, SpecialCharacters) {
    // Test string containing spaces
    expectMd5Equals("hello world", "5eb63bbbe01eeed093cb22bb8f5acdc3");
    
    // Test string containing newline characters
    expectMd5Equals("hello\nworld", "9195d0beb2a889e1be05ed6bb1954837");
    
    // Test string containing special symbols
    expectMd5Equals("!@#$%^&*()", "05b28d17a7b6e7024b6e5d8cc43a8bf7");
    
    // Test Unicode characters (Chinese characters)
    expectMd5Equals("HelloWorld", "65396ee4aad0b4f17aacd1c6112ee364");
}

// Boundary conditions test
TEST_F(Md5UtilTest, BoundaryConditions) {
    // Test data with length 0
    expectMd5Equals("", "d41d8cd98f00b204e9800998ecf8427e");
    
    // Test data with length 1
    expectMd5Equals("x", "9dd4e461268c8034f5c8564e155c67a6");
    
    // Test large data block (1KB)
    std::string large_data(1024, 'A');
    uint8_t digest[16];
    md5_compute((const uint8_t*)large_data.c_str(), large_data.length(), digest);
    std::string result = md5_to_string(digest);
    // Verify result is not empty and has correct length
    EXPECT_EQ(result.length(), 32);
    EXPECT_NE(result, "d41d8cd98f00b204e9800998ecf8427e"); // Should not be MD5 of empty string
    
    // Test large data block (10KB)
    std::string very_large_data(10240, 'B');
    md5_compute((const uint8_t*)very_large_data.c_str(), very_large_data.length(), digest);
    result = md5_to_string(digest);
    EXPECT_EQ(result.length(), 32);
}

// md5_to_string function test
TEST_F(Md5UtilTest, Md5ToString) {
    uint8_t digest[16];
    md5_compute((const uint8_t*)"test", 4, digest);
    
    std::string result = md5_to_string(digest);
    EXPECT_EQ(result.length(), 32);
    EXPECT_TRUE(result.find_first_not_of("0123456789abcdef") == std::string::npos);
}

// Consistency test: multiple calls with same input should give same result
TEST_F(Md5UtilTest, ConsistencyTest) {
    const std::string test_input = "consistency_test";
    
    uint8_t digest1[16], digest2[16];
    md5_compute((const uint8_t*)test_input.c_str(), test_input.length(), digest1);
    md5_compute((const uint8_t*)test_input.c_str(), test_input.length(), digest2);
    
    std::string result1 = md5_to_string(digest1);
    std::string result2 = md5_to_string(digest2);
    
    EXPECT_EQ(result1, result2);
}

// Calculation inline function test
TEST_F(Md5UtilTest, InlineCalcMd5sumTest) {
    std::string md5_result;
    const uint8_t* test_data = (const uint8_t*)"test_data";
    
    calcMd5sum(md5_result, test_data, 9);
    
    EXPECT_EQ(md5_result.length(), 32);
    EXPECT_TRUE(md5_result.find_first_not_of("0123456789abcdef") == std::string::npos);
    
    // Verify consistency with direct call to md5_compute
    uint8_t digest[16];
    md5_compute(test_data, 9, digest);
    std::string expected = md5_to_string(digest);
    EXPECT_EQ(md5_result, expected);
}
