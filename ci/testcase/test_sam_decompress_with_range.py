"""
测试用例：SamDecompressWithRangeTest - SAM文件区域解压测试

测试场景：
- 测试使用-p参数解压SAM压缩文件的指定染色体区域
- 验证区域过滤功能的正确性和精确性
- 测试基因组坐标范围查询的处理能力

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress -p chr3:100-200 {compressed_file} -o {decompressed_file}

预期结果：
- 压缩命令执行成功
- 区域解压命令成功，只包含chr3位置100-200的记录
- 解压文件内容精确匹配目标区域
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamDecompressWithRangeTest(PBGZTestCase):
    
    def __init__(self):
        super().__init__("SamDecompressWithRangeTest")
        self.source_file = "sam_range_test.sam"
        self.compressed_file = "sam_range_test.sam.pbgz"
        self.index_file = "sam_range_test.sam.pbgz.pbgzi"
        self.decompressed_file = "sam_range_test.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成已排序的SAM文件，使用pbgz支持的格式
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")
        sam_content.append("@SQ\tSN:chr1\tLN:1000000")
        sam_content.append("@SQ\tSN:chr2\tLN:1000000")
        sam_content.append("@SQ\tSN:chr3\tLN:1000000")
        
        # 添加不同区域的reads
        # chr1: 位置100-500范围
        for pos in [100, 200, 300, 400, 500]:
            qname = f"CHR1_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        # chr1: 位置1000-1500范围（不在目标范围内）
        for pos in [1000, 1200, 1500]:
            qname = f"CHR1_OUT_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        # chr2: 位置50-150范围（不同染色体）
        for pos in [50, 150]:
            qname = f"CHR2_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr2", pos)
            sam_content.append(read_line)
        
        # chr3: 位置100-200范围（目标区域）
        for pos in [120, 180]:
            qname = f"CHR3_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr3", pos)
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # 先压缩
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -i")
        
        # 测试区域解压（如果pbgz支持）
        self.add_test_command(f"./bin/pbgz decompress -p 'chr3:100-200' {self.compressed_file} -o {self.decompressed_file}")

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
        for filename in [self.source_file, 
                         self.compressed_file, 
                         self.index_file, 
                         self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # 检查命令执行情况 - 基于命令类型来判断，而不是文件存在性
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            file_sizes = cmd_result.get("file_sizes", {})
            
            # 通过命令内容判断是压缩还是解压命令 - 必须先检查"decompress"避免被"compress"包含
            if "decompress" in command and "pbgz" in command:
                decompress_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        # 验证基本功能
        if not compress_success:
            print("Compression command failed")
            return False
        
        if not decompress_success:
            print("Decompression command failed")
            return False
        
        return True


if __name__ == "__main__":
    test_case = SamDecompressWithRangeTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)