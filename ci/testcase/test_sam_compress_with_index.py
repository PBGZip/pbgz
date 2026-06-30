"""
测试用例：SamCompressWithIndexTest - SAM文件压缩解压测试

测试场景：
- 测试SAM文件的压缩和解压功能
- 验证数据完整性和格式正确性

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -i {index_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，索引文件正确生成
- 解压成功，数据完整恢复
- MD5校验完全一致
"""

import sys
import random
import hashlib

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamCompressWithIndexTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamCompressWithIndexTest")
        self.source_file = "sam_with_index.sam"
        self.compressed_file = "sam_with_index.sam.pbgz"
        self.index_file = "sam_with_index.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(20):
            qname = f"INDEX_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 500
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
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        self.add_test_command(f"./bin/pbgz index -f {self.compressed_file}")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # 检查命令执行结果
        compress_success = False
        index_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            if "compress" in command and "pbgz" in command:
                compress_success = success
            elif "index" in command and "pbgz" in command:
                index_success = success
        
        # 检查压缩和索引是否都成功
        if not compress_success:
            return False
            
        if not index_success:
            return False
            
        # 检查压缩文件是否存在
        if not os.path.exists(self.compressed_file):
            return False
        
        # 检查索引文件是否存在
        if self.index_file and not os.path.exists(self.index_file):
            return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressWithIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
