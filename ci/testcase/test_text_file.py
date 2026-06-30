"""
测试用例：TextFileTest - 文本文件压缩解压测试

测试场景：
- 测试pbgz对普通文本文件的压缩和解压
- 验证文本内容的完整性保留
- 测试多样化字符（字母、数字、符号）的处理

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

参数说明：
- 测试文件包含1000行文本，每行包含多样化字符
- 字符包括：大小写字母、数字、特殊符号
- 每行长度约为100字符

预期结果：
- 压缩成功，生成压缩文件
- 解压成功，恢复原始文本文件
- 解压后文件的MD5与原始文件完全一致
- 所有字符和行结构正确保留
- 压缩率合理（因为文本有重复模式）

技术说明：
- 测试最常用的文本文件类型
- 验证pbgz处理通用文本数据的能力
- 文本文件压缩是基础功能验证
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class TextFileTest(PBGZTestCase):
    def __init__(self):

        super().__init__("TextFileTest")
        self.source_file = "text_file_data.txt"
        self.compressed_file = "text_file_data.txt.pbgz"
        self.decompressed_file = "text_file_data.txt.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成文本文件内容
        full_text = []
        num_lines = 1000  # 生成1000行文本
        
        for i in range(num_lines):
            line = f"This is line {i} of text file data for compression testing. "
            line += "Each line contains various characters: ABCDEF123456!@#$%^&*()"
            full_text.append(line)
        
        self.original_hash = hashlib.md5('\n'.join(full_text).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in full_text:
                f.write(line + '\n')
        
        
        # 添加压缩测试命令
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # 添加解压测试命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
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
        
        
        # 对比原文件和解压文件是否一致
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read().encode('utf-8')
            decompressed_content = f2.read().encode('utf-8')
        
        decompressed_hash = hashlib.md5(decompressed_content).hexdigest()
        original_hash = hashlib.md5(original_content).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = TextFileTest()
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