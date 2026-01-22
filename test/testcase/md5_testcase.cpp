#include <gtest/gtest.h>
#include "utils/md5_util.h"
#include <string>
#include <vector>

class Md5UtilTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 测试前的初始化
    }

    void TearDown() override {
        // 测试后的清理
    }

    // 辅助函数：比较MD5结果
    void expectMd5Equals(const std::string& input, const std::string& expected_hex) {
        uint8_t digest[16];
        md5_compute((const uint8_t*)input.c_str(), input.length(), digest);
        std::string result = md5_to_string(digest);
        EXPECT_EQ(result, expected_hex) << "Input: " << input;
    }

    // 辅助函数：比较二进制数据的MD5
    void expectMd5EqualsBinary(const std::vector<uint8_t>& data, const std::string& expected_hex) {
        uint8_t digest[16];
        md5_compute(data.data(), data.size(), digest);
        std::string result = md5_to_string(digest);
        EXPECT_EQ(result, expected_hex);
    }
};

// 基本功能测试
TEST_F(Md5UtilTest, BasicMd5Compute) {
    // 测试空字符串
    expectMd5Equals("", "d41d8cd98f00b204e9800998ecf8427e");
    
    // 测试单字符
    expectMd5Equals("a", "0cc175b9c0f1b6a831c399e269772661");
    
    // 测试简单字符串
    expectMd5Equals("abc", "900150983cd24fb0d6963f7d28e17f72");
    
    // 测试数字字符串
    expectMd5Equals("123456789", "25f9e794323b453885f5181f1b624d0b");
}

// RFC 1321 MD5测试向量
TEST_F(Md5UtilTest, Rfc1321TestVectors) {
    // RFC 1321中的标准测试向量
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

// 不同长度数据测试
TEST_F(Md5UtilTest, DifferentLengths) {
    // 测试1字节数据
    expectMd5Equals("1", "c4ca4238a0b923820dcc509a6f75849b");
    
    // 测试16字节数据（正好一个MD5块）
    expectMd5Equals("1234567890123456", "abeac07d3c28c1bef9e730002c753ed4");
    
    // 测试32字节数据（两个MD5块）
    expectMd5Equals("12345678901234567890123456789012", "767179c7a2bff19651ce97d294c30cfb");
    
    // 测试64字节数据（正好两个MD5块）
    expectMd5Equals("1234567890123456789012345678901234567890123456789012345678901234", 
                   "eb6c4179c0a7c82cc2828c1e6338e165");
    
    // 测试65字节数据（超过两个块）
    expectMd5Equals("12345678901234567890123456789012345678901234567890123456789012345", 
                   "823cc889fc7318dd33dde0654a80b70a");
}

// 二进制数据测试
TEST_F(Md5UtilTest, BinaryData) {
    // 测试纯二进制数据
    std::vector<uint8_t> binary_data = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05};
    expectMd5EqualsBinary(binary_data, "d15ae53931880fd7b724dd7888b4b4ed");
    
    // 测试包含null字符的数据
    std::vector<uint8_t> data_with_null = {'a', 0, 'b', 0, 'c'};
    expectMd5EqualsBinary(data_with_null, "ed4b9e8bfe6c7571a700dc6a0b0d01cb");
    
    // 测试全零数据
    std::vector<uint8_t> all_zeros(32, 0);
    expectMd5EqualsBinary(all_zeros, "70bc8f4b72a86921468bf8e8441dce51");
}

// 特殊字符测试
TEST_F(Md5UtilTest, SpecialCharacters) {
    // 测试包含空格的字符串
    expectMd5Equals("hello world", "5eb63bbbe01eeed093cb22bb8f5acdc3");
    
    // 测试包含换行符的字符串
    expectMd5Equals("hello\nworld", "9195d0beb2a889e1be05ed6bb1954837");
    
    // 测试包含特殊符号的字符串
    expectMd5Equals("!@#$%^&*()", "05b28d17a7b6e7024b6e5d8cc43a8bf7");
    
    // 测试Unicode字符（中文字符）
    expectMd5Equals("你好世界", "65396ee4aad0b4f17aacd1c6112ee364");
}

// 边界条件测试
TEST_F(Md5UtilTest, BoundaryConditions) {
    // 测试长度为0的数据
    expectMd5Equals("", "d41d8cd98f00b204e9800998ecf8427e");
    
    // 测试长度为1的数据
    expectMd5Equals("x", "9dd4e461268c8034f5c8564e155c67a6");
    
    // 测试大块数据（1KB）
    std::string large_data(1024, 'A');
    uint8_t digest[16];
    md5_compute((const uint8_t*)large_data.c_str(), large_data.length(), digest);
    std::string result = md5_to_string(digest);
    // 验证结果不为空且长度正确
    EXPECT_EQ(result.length(), 32);
    EXPECT_NE(result, "d41d8cd98f00b204e9800998ecf8427e"); // 不应该是空字符串的MD5
    
    // 测试大块数据（10KB）
    std::string very_large_data(10240, 'B');
    md5_compute((const uint8_t*)very_large_data.c_str(), very_large_data.length(), digest);
    result = md5_to_string(digest);
    EXPECT_EQ(result.length(), 32);
}

// md5_to_string函数测试
TEST_F(Md5UtilTest, Md5ToString) {
    uint8_t digest[16];
    md5_compute((const uint8_t*)"test", 4, digest);
    
    std::string result = md5_to_string(digest);
    EXPECT_EQ(result.length(), 32);
    EXPECT_TRUE(result.find_first_not_of("0123456789abcdef") == std::string::npos);
}

// 一致性测试：多次调用相同输入应该得到相同结果
TEST_F(Md5UtilTest, ConsistencyTest) {
    const std::string test_input = "consistency_test";
    
    uint8_t digest1[16], digest2[16];
    md5_compute((const uint8_t*)test_input.c_str(), test_input.length(), digest1);
    md5_compute((const uint8_t*)test_input.c_str(), test_input.length(), digest2);
    
    std::string result1 = md5_to_string(digest1);
    std::string result2 = md5_to_string(digest2);
    
    EXPECT_EQ(result1, result2);
}

// 计算内联函数测试
TEST_F(Md5UtilTest, InlineCalcMd5sumTest) {
    std::string md5_result;
    const uint8_t* test_data = (const uint8_t*)"test_data";
    
    calcMd5sum(md5_result, test_data, 9);
    
    EXPECT_EQ(md5_result.length(), 32);
    EXPECT_TRUE(md5_result.find_first_not_of("0123456789abcdef") == std::string::npos);
    
    // 验证与直接调用md5_compute的结果一致
    uint8_t digest[16];
    md5_compute(test_data, 9, digest);
    std::string expected = md5_to_string(digest);
    EXPECT_EQ(md5_result, expected);
}
