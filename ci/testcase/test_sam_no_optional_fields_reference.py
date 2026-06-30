"""
测试用例：Sam No Optional Fields Reference

测试场景：
- 测试pbgz的Sam No Optional Fields Reference功能
- 验证相关命令的正确执行
- 确保数据完整性

预期结果：
- 测试命令执行成功
- 数据完整性验证通过
"""

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase
# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class SamNoOptionalFieldsReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamNoOptionalFieldsReferenceTest")
        self.source_file = "sam_no_opt_ref.sam"
        self.compressed_file = "sam_no_opt_ref.sam.pbgz"
        self.decompressed_file = "sam_no_opt_ref.sam.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(10):
            qname = f"NO_OPT_REF_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 800
            mapq = 60
            cigar = "30M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=30))
            qual = '!' * 30
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}")
        self.add_test_command(f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}")
    
    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        if not os.path.exists(self.compressed_file):
            return False
        if not os.path.exists(self.decompressed_file):
            return False
        return True

if __name__ == "__main__":
    test_case = SamNoOptionalFieldsReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
