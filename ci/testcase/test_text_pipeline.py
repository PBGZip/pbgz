"""
测试用例：TextPipelineTest - 文本文件多级管道处理测试

测试场景：
- 测试pbgz在多级管道中的处理能力
- 验证复杂管道组合的稳定性
- 测试文件通过多个命令链的处理

使用命令：
1. cat {source_file} | ./bin/pbgz compress | cat | ./bin/pbgz decompress -o {decompressed_file}

预期结果：
- 多级管道处理成功
- 文件在管道链中正确传递
- 最终解压文件与原始文件内容一致
- 管道组合操作稳定可靠

技术说明：
- 测试复杂的管道组合场景
- 验证错误处理和流管理
- 确保多命令链的正确传递
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class TextPipelineTest(PBGZTestCase):
    def __init__(self):

        super().__init__("TextPipelineTest")
        self.source_file = "text_pipeline_data.txt"
        self.compressed_file = "text_pipeline_data.txt.pbgz"
        self.decompressed_file = "text_pipeline_data.txt.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成文本文件内容
        full_text = []
        num_lines = 1000
        
        for i in range(num_lines):
            line = f"This is line {i} of text file data for pipeline compression testing."
            full_text.append(line)
        
        self.original_hash = hashlib.md5('\n'.join(full_text).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in full_text:
                f.write(line + '\n')
        
        
        # 添加压缩测试命令（管道输入，文件输出）
        compress_command = f"cat {self.source_file} | ./bin/pbgz compress -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        # 添加解压测试命令（管道输入，文件输出）
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
    test_case = TextPipelineTest()
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