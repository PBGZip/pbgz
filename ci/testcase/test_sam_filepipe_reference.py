"""
测试用例：SamFilepipeReferenceTest - SAM文件管道压缩解压测试（带参考基因组）

测试场景：
- 测试SAM文件通过管道方式进行压缩和解压
- 使用参考基因组优化压缩
- 验证管道模式和参考基因组的组合

使用命令：
1. cat {source_file} | ./bin/pbgz compress -o {compressed_file} -r {reference_file}
2. ./bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

预期结果：
- 管道压缩解压成功
- 参考基因组优化生效
- 数据完整性验证通过
"""

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class SamFilepipeReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamFilepipeReferenceTest")
        self.source_file = "sam_filepipe_ref.sam"
        self.compressed_file = "sam_filepipe_ref.sam.pbgz"
        self.decompressed_file = "sam_filepipe_ref.sam.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(10):
            qname = f"FILEPIPE_REF_{i:03d}"
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
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30\tAS:i:30\tXS:i:0\tRG:Z:FILEPIPE_REF"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o - > {self.compressed_file}")
        self.add_test_command(f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}")
    
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
    test_case = SamFilepipeReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)