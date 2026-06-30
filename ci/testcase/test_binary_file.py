"""
测试用例：BinaryFileTest - 基本二进制文件压缩解压测试

测试场景：
- 测试pbgz对基本的二进制文件进行压缩和解压功能
- 验证压缩率和数据完整性
- 基础功能测试，无特殊参数

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩命令成功执行，生成压缩文件
- 解压命令成功执行，恢复原始文件
- 解压后的文件与原始文件MD5完全一致
- 压缩率合理（应显示具体数值）
- 文件内容完整性验证通过
"""


import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class BinaryFileTest(PBGZTestCase):
    def __init__(self):
        super().__init__("BinaryFileTest")
        self.source_file = "test_data.bin"
        self.compressed_file = "test_data.bin.pbgz"
        self.decompressed_file = "test_data.bin.dec"
    
    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成一个非fastq非sam格式的二进制文件
        binary_data = bytearray(random.getrandbits(8) for _ in range(1024 * 1024))  # 1MB随机二进制数据
        self.original_hash = hashlib.md5(binary_data).hexdigest()
        
        with open(self.source_file, 'wb') as f:
            f.write(binary_data)
        
        
        # 添加压缩测试命令
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # 添加解压测试命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
    def cleanup_test_data(self):
        """清理测试用例创建的临时文件"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filepath in files_to_clean:
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                print(f"Warning: Failed to remove {filepath}: {e}")
    
    def verify_expected_results(self) -> bool:
        # 验证压缩文件是否生成
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        # 验证解压文件是否生成
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        # 验证压缩比是否合理
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        # 打印获取到的执行时间和压缩比
        exec_time = self.get_execution_time()
        
        
        # 验证压缩比是否有效
        if compression_ratio <= 0:
            print(f"Error: Invalid compression ratio: {compression_ratio:.2f}%")
            return False
        
        # 对比原文件和解压文件是否一致
        with open(self.source_file, 'rb') as f1, open(self.decompressed_file, 'rb') as f2:
            original_data = f1.read()
            decompressed_data = f2.read()
        
        # 使用哈希值对比，避免内存对比问题
        decompressed_hash = hashlib.md5(decompressed_data).hexdigest()
        original_hash = hashlib.md5(original_data).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = BinaryFileTest()
    result = test_case.execute()
    test_case.print_results()
    test_case.save_to_json()
    
    if result:
        print("\\nTest passed successfully!")
    else:
        print("\\nTest failed!")
        sys.exit(1)