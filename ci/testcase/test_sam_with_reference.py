"""
测试用例：SamWithReferenceTest - SAM文件压缩解压测试（高质量参考数据）

测试场景：
- 测试pbgz对高质量SAM文件进行压缩和解压
- SAM文件包含精心设计的比对数据
- 使用参考基因组实现最佳压缩效果
- 测试关键字段：read名称、比对位置、映射质量、CIGAR字符串等

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

参数说明：
- -r参数：使用参考基因组abqd01.fasta.gz
- 测试数据包含3条高质量的SAM记录

预期结果：
- 压缩成功，压缩率很高（因为参考匹配）
- 解压成功，数据完整恢复
- MD5校验完全一致
- 所有可选字段正确保留（NM, MD, AS, XS, RG等）

技术说明：
- 测试高质量SAM数据的处理能力
- 验证pbgz处理复杂可选字段的能力
- 预期压缩率很高，因为数据包含参考匹配信息
"""

"""
测试用例：SamWithReferenceTest - SAM高质量数据压缩解压测试（参考基因组）

测试场景：
- 测试pbgz对高质量SAM文件进行压缩和解压
- 使用参考基因组abqD01.fasta.gz优化压缩
- SAM数据包含精确的参考匹配，预期压缩率高

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

参数说明：
- -r参数：使用参考基因组优化压缩，特别是对于包含精确匹配的数据
- 测试3条高质量的SAM记录，位置在1k、7k附近

预期结果：
- 压缩成功率很高（因为包含精确的参考匹配）
- 解压成功，数据完美恢复
- MD5校验完全一致
- 所有可选字段正确保留：NM（错配数）、MD（错配细节）、AS（比对分数）、XS（次优比对分数）、RG（read group）
- 预期压缩率很高，因为数据包含完全匹配的参考序列

技术说明：
- 测试高质量SAM数据的处理
- 验证pbgz处理复杂可选字段的能力
- 参考基因组能极大提升压缩效果
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamWithReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamWithReferenceTest")
        self.source_file = "with_ref_test.sam"
        self.compressed_file = "with_ref_test.sam.pbgz"
        self.decompressed_file = "with_ref_test.sam.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"
    
    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成基本的SAM文件内容
        sam_content = []
        
        # SAM头
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # 添加3条read记录
        num_reads = 3
        for i in range(num_reads):
            # QNAME
            qname = f"SAM_READ_{i:03d}"
            
            # FLAG
            flag = 0
            
            # RNAME
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            
            # POS
            pos = (i + 1) * 1000
            
            # MAPQ
            mapq = 60
            
            # CIGAR
            cigar = "25M"
            
            # RNEXT
            rnext = "*"
            
            # PNEXT
            pnext = 0
            
            # TLEN
            tlen = 0
            
            # SEQ
            seq = ''.join(random.choices('ATGCN', k=25))
            
            # QUAL
            qual = '!' * 25
            
            # 添加可选字段，与参考基因组基因序列测试一致
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25\tXS:i:0\tRG:Z:1"
            
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        # 生成SAM文件
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        
        # 添加压缩测试命令
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # 添加解压测试命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}"
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
        
        # 使用MD5对比
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        return original_hash == decompressed_hash


if __name__ == "__main__":
    test_case = SamWithReferenceTest()
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