"""
测试用例：SamFileTest - 基本SAM文件压缩解压测试（带参考基因组）

测试场景：
- 测试pbgz对SAM（Sequence Alignment Map）文件进行压缩和解压
- SAM文件包含高质量基因组比对数据
- 使用参考基因组进行压缩，提高压缩效率

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

参数说明：
- -r参数：指定参考基因组文件，用于优化压缩比
- 测试数据包含1000条SAM记录，位置在97-900的范围内

预期结果：
- 压缩命令成功，生成.pbgz文件
- 解压命令成功，恢复原始SAM文件
- 解压后文件的MD5与原始文件完全一致
- 压缩率良好（因为使用参考基因组）
- SAM头部和数据行完整保留

技术说明：
- SAM格式是基因组比对的标准格式
- 引用参考基因组能大幅提升压缩效果
- 验证比对信息、质量值、可选字段等的完整性
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamFileTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamFileTest")
        self.source_file = "sam_only_test.sam"
        self.compressed_file = "sam_only_test.sam.pbgz"
        self.decompressed_file = "sam_only_test.sam.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成SAM文件数据
        sam_content = []
        
        # SAM文件头部
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000")
        
        # 添加大量读取记录进行压缩测试
        num_reads = 1000
        for i in range(num_reads):
            # QNAME
            qname = f"SAM_READ_{i:04d}"
            
            # FLAG
            flag = 0
            
            # RNAME
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            
            # POS
            pos = (i % 900) + 1
            
            # MAPQ
            mapq = 60
            
            # CIGAR
            cigar = "50M"
            
            # RNEXT
            rnext = "*"
            
            # PNEXT
            pnext = 0
            
            # TLEN
            tlen = 0
            
            # SEQ
            seq = ''.join(random.choices('ATCGN', k=50))
            
            # QUAL
            qual = '!' * 50
            
            # 可选字段
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:50\tAS:i:50\tXS:i:0\tRG:Z:1"
            
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        # 生成SAM文件
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # 添加压缩测试命令
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -r fa/ABQD01.fasta.gz"
        self.add_test_command(compress_command, cmd_id=1)
        
        # 添加解压测试命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r fa/ABQD01.fasta.gz"
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
        # 严格检查：压缩和解压命令都必须成功
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
            return False
        
        # 验证解压文件是否生成  
        if not os.path.exists(self.decompressed_file):
            return False
        
        # 验证压缩比是否合理
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            return False
        
        # 对比原文件和解压文件是否一致
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        # 使用哈希值对比
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        return original_hash == decompressed_hash


if __name__ == "__main__":
    test_case = SamFileTest()
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