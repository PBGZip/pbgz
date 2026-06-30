"""
测试用例：SamMismatchedCigarBaseLengthTest - 碱基字段长度和CIGAR不匹配测试

测试场景：
- 测试pbgz处理SAM文件中碱基字段长度和CIGAR不匹配的情况
- 验证对于这类 不一致性，pbgz仍能正常压缩和解压
- 确保数据完整性不受影响

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，即使存在碱基长度和CIGAR不匹配的记录
- 解压成功，恢复原始SAM文件
- 解压后文件的MD5与原始文件完全一致
- 所有不一致记录的原始信息被完整保留
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamMismatchedCigarBaseLengthTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMismatchedCigarBaseLengthTest")
        self.source_file = "sam_mismatched_cigar.sam"
        self.compressed_file = "sam_mismatched_cigar.sam.pbgz"
        self.decompressed_file = "sam_mismatched_cigar.sam.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成包含碱基字段长度和CIGAR不匹配的SAM文件
        sam_content = []
        
        # SAM文件头部
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # 添加正常reads作为对照 - 增加数量以达到32M数据量
        for i in range(20000):
            qname = f"NORMAL_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = "150M"
            rnext = "*"
            pnext = 0
            tlen = 0
            # 正常长度与CIGAR匹配
            seq = ''.join(random.choices('ATCGN', k=150))
            qual = '!' * 150
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # 添加碱基长度大于CIGAR的reads - 增加数量
        for i in range(20000):
            qname = f"SEQ_LONGER_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 2000
            mapq = 60
            cigar = "100M"  # CIGAR表示100个碱基
            rnext = "*"
            pnext = 0
            tlen = 0
            # 实际序列长度大于CIGAR (150 vs 100)
            seq = ''.join(random.choices('ATCGN', k=150))
            qual = '!' * 150
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # 添加碱基长度小于CIGAR的reads - 增加数量
        for i in range(20000):
            qname = f"SEQ_SHORTER_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 3000
            mapq = 60
            cigar = "150M"  # CIGAR表示150个碱基
            rnext = "*"
            pnext = 0
            tlen = 0
            # 实际序列长度小于CIGAR (100 vs 150)
            seq = ''.join(random.choices('ATCGN', k=100))
            qual = '!' * 100
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # 添加复杂CIGAR与序列长度不匹配的reads - 增加数量
        for i in range(20000):
            qname = f"COMPLEX_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 4000
            mapq = 60
            cigar = "50M10I50M10D50M"  # 复杂CIGAR
            rnext = "*"
            pnext = 0
            tlen = 0
            # 序列长度与CIGAR明显不匹配
            seq = ''.join(random.choices('ATCGN', k=200))
            qual = '!' * 200
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
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
        # 验证压缩和解压命令是否成功
        compress_success = decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # 先检查decompress，避免"decompress"包含"compress"的问题
            if "decompress" in command and "pbgz" in command:
                decompress_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        if not compress_success or not decompress_success:
            print(f"压缩或解压命令失败: compress_success={compress_success}, decompress_success={decompress_success}")
            return False
        
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
        
        exec_time = self.get_execution_time()
        print(f"Verification: Compression ratio = {compression_ratio:.2f}%, Execution time = {exec_time:.2f}s")
        
        # 对比原文件和解压文件是否一致 - 对于不匹配的case，应该保留原始数据
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
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
    test_case = SamMismatchedCigarBaseLengthTest()
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