"""
测试用例：FastqInvalidBaseTest - FASTQ文件非法碱基处理测试

测试场景：
- 测试FASTQ文件中包含非法碱基（非ATCGN字符）的情况
- 验证pbgz对非法碱基的容错能力和处理策略
- 测试数据格式验证和异常处理

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

参数说明：
- FASTQ标准碱基为A、T、C、G、N，测试非法字符Z等
- 验证pbgz对异常数据的处理方式
- 测试格式严格性与容错性的平衡

预期结果：
- 压缩成功，处理非法碱基数据
- 解压成功，保留原始数据（包括非法字符）
- MD5校验完全一致
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqInvalidBaseTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqInvalidBaseTest")
        self.source_file = "invalid_base_test.fq"
        self.compressed_file = "invalid_base_test.fq.pbgz"
        self.decompressed_file = "invalid_base_test.fq.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成包含无效碱基的FASTQ文件
        fastq_content = []
        num_reads = 50
        
        for i in range(num_reads):
            # ID行
            seq_id = f"@INVALID{i:06d} 1/1"
            
            # 序列内容：在随机序列中插入一些无效字符
            valid_bases = ''.join(random.choices('ATGCN', k=45))
            # 在序列中插入无效字符（非ATGCN）
            if i % 3 == 0:
                seq = valid_bases + 'XZ'  # 插入X和Z
            elif i % 3 == 1:
                seq = valid_bases + 'KLM'  # 插入KLM
            else:
                seq = valid_bases + 'OPQR'  # 插入OPQR
            
            # + 行
            separator = '+'
            
            # 质量分数行
            quality_length = len(seq)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, seq, separator, quality])
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        
        # 添加压缩测试命令
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}"
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
        
        
        # 对比原文件和解压文件是否一致
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqInvalidBaseTest()
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