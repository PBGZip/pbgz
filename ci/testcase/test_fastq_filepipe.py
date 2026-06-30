"""
测试用例：FastqFilepipeTest - FASTQ文件管道压缩解压测试

测试场景：
- 测试FASTQ格式文件通过管道方式进行压缩和解压
- 验证FASTQ格式与标准Unix管道的兼容性
- 测试cat命令与pbgz的组合使用

使用命令：
1. cat {source_file} | ./bin/pbgz compress -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

预期结果：
- 管道压缩解压成功
- FASTQ文件格式完整保留（@、+分隔符）
- 解压后文件的MD5与原始文件完全一致
- 管道模式功能正常工作
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqFilepipeTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqFilepipeTest")
        self.source_file = "filepipe_test.fq"
        self.compressed_file = "filepipe_test.fq.pbgz"
        self.decompressed_file = "filepipe_test.fq.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成正常的FASTQ文件数据
        fastq_content = []
        num_reads = 200  # 生成200条序列
        
        for i in range(num_reads):
            # ID行：@开头的序列标识符
            seq_id = f"@FILEPIPE{i:06d} 1/1"
            
            # 序列内容：只包含ATGCN的随机序列
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            # + 行：分隔符
            separator = '+'
            
            # 质量分数行：生成与序列长度匹配的质量分数（使用常见的Phred质量分数范围）
            quality_length = len(bases)
            # 使用ASCII 33-73范围内的字符（对应Phred质量分数0-40）
            quality = ''.join(random.choice('\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, bases, separator, quality])
        
        full_content = 'n'.join(fastq_content) + 'n'  # 正常FASTQ文件以换行符结尾
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        

        
        # 添加压缩测试命令：文件输入，管道输出
        # -o - 参数将压缩数据输出到stdout，> 重定向到文件
        compress_command = f"./bin/pbgz compress {self.source_file} -o - > {self.compressed_file}"
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
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        # 使用哈希值对比
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqFilepipeTest()
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
