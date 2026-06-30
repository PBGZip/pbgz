"""
测试用例：SamMissingChromosomeInHeaderTest - SAM文件缺少染色体头部定义测试

测试场景：
- 测试SAM文件中实际数据包含的染色体没有在@SQ header中定义的情况
- 验证pbgz对不完整SAM格式的容错能力
- 测试染色体引用完整性检查的边界处理

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，pbgz能够处理缺少染色体定义的数据
- 解压成功，所有数据完整保留
- MD5校验完全一致
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamMissingChromosomeInHeaderTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMissingChromosomeInHeaderTest")
        self.source_file = "sam_missing_chr.sam"
        self.compressed_file = "sam_missing_chr.sam.pbgz"
        self.decompressed_file = "sam_missing_chr.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        
        # Header中只定义部分染色体
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")
        sam_content.append("@SQ\tSN:chr1\tLN:1000000")  # 定义了chr1
        sam_content.append("@SQ\tSN:chr2\tLN:500000")   # 定义了chr2
        sam_content.append("@SQ\tSN:chr3\tLN:2000000")  # 定义了chr3
        # 注意：这里没有定义chr4
        
        # 添加包含已定义染色体的数据
        for pos in [100, 200, 300]:
            qname = f"CHR1_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        for pos in [50, 150]:
            qname = f"CHR2_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr2", pos)
            sam_content.append(read_line)
        
        # 关键点：添加包含未在header中定义的染色体的数据（chr4）
        for pos in [1000, 1200]:
            qname = f"CHR4_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr4", pos)  # chr4未在header中定义
            sam_content.append(read_line)
        
        # 再添加一些已定义染色体的数据
        for pos in [400, 500]:
            qname = f"CHR3_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr3", pos)
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # 压缩
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        
        # 解压
        self.add_test_command(f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}")

    def _create_sam_line(self, qname, rname, pos):
        """创建SAM记录行"""
        flag = 0
        mapq = 60
        cigar = "30M"
        rnext = "*"
        pnext = 0
        tlen = 0
        seq = ''.join(random.choices('ATCGN', k=30))
        qual = '!' * 30
        return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30"

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass

    def verify_expected_results(self) -> bool:
        # 检查命令执行情况
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # 基于命令类型精确判断
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "compress":
                    compress_success = success
                elif main_cmd == "decompress":
                    decompress_success = success
        
        # 验证基本功能
        if not compress_success or not decompress_success:
            print(f"压缩成功: {compress_success}, 解压成功: {decompress_success}")
            return False
        
        # 验证压缩文件存在
        if not os.path.exists(self.compressed_file):
            print("压缩文件不存在")
            return False
        
        # 验证解压文件存在
        if not os.path.exists(self.decompressed_file):
            print("解压文件不存在")
            return False
        
        # 验证数据完整性
        with open(self.source_file, 'r') as f:
            original_lines = f.readlines()
        
        with open(self.decompressed_file, 'r') as f:
            decompressed_lines = f.readlines()
        
        if len(original_lines) != len(decompressed_lines):
            print(f"行数不匹配: 原始{len(original_lines)}行 vs 解压{len(decompressed_lines)}行")
            return False
        
        # MD5校验
        with open(self.source_file, 'rb') as f:
            original_md5 = hashlib.md5(f.read()).hexdigest()
        
        with open(self.decompressed_file, 'rb') as f:
            decompressed_md5 = hashlib.md5(f.read()).hexdigest()
        
        if original_md5 != decompressed_md5:
            print(f"MD5不匹配: {original_md5} vs {decompressed_md5}")
            return False
        
        # 验证关键点：包含未定义染色体的记录
        has_un_defined_chr = False
        for line in decompressed_lines:
            if line.startswith('@'):  # 跳过header
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3 and fields[2] == "chr4":
                has_un_defined_chr = True
                break
        
        if not has_un_defined_chr:
            print("解压文件中缺少未定义染色体的记录(chr4)")
            return False
        
        print("SAM文件包含未在header中定义的染色体时也能正常压缩解压")
        return True


if __name__ == "__main__":
    test_case = SamMissingChromosomeInHeaderTest()
    result = test_case.execute()
    test_case.print_results()
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)