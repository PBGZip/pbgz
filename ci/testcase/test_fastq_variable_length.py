"""
测试用例：FastqVariableLengthTest - 变长序列FASTQ文件压缩解压测试

测试场景：
- 测试包含不同长度序列的FASTQ文件处理能力
- 验证pbgz对多样长度序列的适应性
- 测试真实测序数据的特点（序列长度不一）

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

参数说明：
- 测序数据通常包含不同长度的序列片段
- 变长序列是真实生物数据的典型特征
- 测试pbgz对序列多样性处理的稳定性

预期结果：
- 压缩成功，有效压缩变长序列
- 解压成功，所有序列完整恢复
- MD5校验完全一致
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqVariableLengthTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqVariableLengthTest")
        self.source_file = "variable_length_test.fq"
        self.compressed_file = "variable_length_test.fq.pbgz"
        self.decompressed_file = "variable_length_test.fq.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        fastq_content = []
        num_reads = 200
        
        for i in range(num_reads):
            seq_id = f"@VAR{i:06d} 1/1"
            
            bases = ''.join(random.choices('ATGCN', k=random.randint(20, 200)))
            
            separator = '+'
            
            quality_length = len(bases)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, bases, separator, quality])
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        

        
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        exec_time = self.get_execution_time()
        
        
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
    test_case = FastqVariableLengthTest()
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