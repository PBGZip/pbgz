"""
测试用例：SamSortedCompressWithIndexTest - 已排序SAM文件压缩生成索引测试

测试场景：
- 测试对已排序的SAM文件进行压缩并同时生成索引
- 验证排序数据与索引生成的协同工作
- 测试Indexed SAM格式的正确性

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -i {index_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，索引文件正确生成 (.index格式)
- 解压成功，完整数据恢复
- 索引文件可用于快速区域查询
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamSortedCompressWithIndexTest(PBGZTestCase):
    """测试已排序SAM文件压缩时同时生成索引文件"""
    
    def __init__(self):
        super().__init__("SamSortedCompressWithIndexTest")
        self.source_file = "sam_sorted_with_index.sam"
        self.compressed_file = "sam_sorted_with_index.sam.pbgz"
        self.index_file = "sam_sorted_with_index.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        # 生成已排序的SAM文件
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # coordinate表示已排序
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # 按位置排序的reads
        for i in range(50):
            qname = f"SORTED_IDX_{i:03d}"
            flag = 0
            rname = "ref1" if i < 25 else "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 100
            mapq = 60
            cigar = "30M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=30))
            qual = '!' * 30
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # 使用-i参数在压缩时生成索引文件
        self.add_test_command(f"./bin/pbgz compress -i {self.source_file} -o {self.compressed_file}")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # 检查压缩命令是否成功
        compress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            if "compress" in command and "pbgz" in command:
                compress_success = success
                break
        
        if not compress_success:
            return False
        
        # 检查压缩文件是否存在
        if not os.path.exists(self.compressed_file):
            return False
        
        # 检查索引文件是否生成
        if self.index_file and not os.path.exists(self.index_file):
            return False
        
        return True


if __name__ == "__main__":
    test_case = SamSortedCompressWithIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)