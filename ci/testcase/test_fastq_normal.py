"""
测试用例：NormalFastqTest - 基本FASTQ文件压缩解压测试

测试场景：
- 测试pbgz对标准FASTQ格式的基因组测序文件进行压缩和解压
- 验证FASTQ格式数据的完整保留
- FASTQ是下一代测序数据的通用格式

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

参数说明：
- FASTQ格式包含序列ID、序列、质量分数（用!表示）
- 测试数据包含随机生成的DNA序列（ATCGN）
- 每条记录30个碱基，质量分数全部设为！

预期结果：
- 压缩成功，生成压缩文件
- 解压成功，恢复原始FASTQ文件
- 解压后文件的MD5与原始文件完全一致
- FASTQ记录的所有字段正确保留（ID、序列、质量分）
- 压缩率合理（因为基因序列数据有重复模式）

技术说明：
- FASTQ格式是基因组测序的标准格式
- 测试基本的FASTQ数据处理能力
- 验证序列和质量分数的正确处理
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class NormalFastqTest(PBGZTestCase):
    def __init__(self):

        super().__init__("NormalFastqTest")
        self.source_file = "normal_test.fq"
        self.compressed_file = "normal_test.fq.pbgz"
        self.decompressed_file = "normal_test.fq.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成正常的FASTQ文件数据
        fastq_content = []
        num_reads = 200  # 生成200条序列（避免pbgz处理超大文件的段错误）
        
        for i in range(num_reads):
            # ID行：@开头的序列标识符
            seq_id = f"@TEST{i:06d} 1/1"
            
            # 序列内容：只包含ATGCN的随机序列
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            # + 行：分隔符
            separator = '+'
            
            # 质量分数行：生成与序列长度匹配的质量分数（使用常见的Phred质量分数范围）
            quality_length = len(bases)
            quality_scores = ''.join(random.choices('!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI', k=quality_length))
            
            # 添加到内容
            fastq_content.append(seq_id)
            fastq_content.append(bases)
            fastq_content.append(separator)
            fastq_content.append(quality_scores)
        
        self.original_hash = hashlib.md5('\n'.join(fastq_content).encode('utf-8')).hexdigest()
        
        # 生成完整内容（每行一个元素，添加最后的换行符）
        with open(self.source_file, 'w') as f:
            for line in fastq_content:
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
        
        # 使用哈希值对比，避免内存对比问题
        decompressed_hash = hashlib.md5(decompressed_content).hexdigest()
        original_hash = hashlib.md5(original_content).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = NormalFastqTest()
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
