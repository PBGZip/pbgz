"""
测试用例：BinaryPipelineTest - 二进制文件多级管道处理测试

测试场景：
- 测试pbgz在管道中进行多级处理的能力
- 验证管道中的中间处理步骤
- 测试复杂管道场景下的稳定性

使用命令：
1. cat {source_file} | ./bin/pbgz compress | cat | ./bin/pbgz decompress -o {decompressed_file}

预期结果：
- 多级管道处理成功
- 压缩和解压在管道中正常工作
- 中间的cat命令集成正常
- 最终解压文件与原始文件一致

技术说明：
- 测试复杂的管道组合
- 验证pbgz在管道链中的稳定性
- 测试错误处理和流管理
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class BinaryPipelineTest(PBGZTestCase):
    def __init__(self):

        super().__init__("BinaryPipelineTest")
        self.source_file = "pipeline_test_data.bin"
        self.compressed_file = "pipeline_test_data.bin.pbgz"
        self.decompressed_file = "pipeline_test_data.bin.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成一个非fastq非sam格式的二进制文件
        binary_data = bytearray(random.getrandbits(8) for _ in range(1024 * 1024))  # 1MB随机二进制数据
        self.original_hash = hashlib.md5(binary_data).hexdigest()
        
        with open(self.source_file, 'wb') as f:
            f.write(binary_data)
        
        
        # 添加测试命令：管道输入，文件输出
        compress_command = f"cat {self.source_file} | ./bin/pbgz compress -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        # 添加解压测试命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        """清理测试用例创建的临时文件"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
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
    test_case = BinaryPipelineTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # 为调试结果保存JSON文件
        if not result:
            test_case.save_to_json()
        sys.exit(1)
