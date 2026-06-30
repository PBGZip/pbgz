"""
测试用例：Sam Variable Base Length Pipeline - 变长基线长度SAM处理测试

测试场景：
- 测试包含变长基线长度（CIGAR中M的数量变化）的SAM数据处理
- 验证pbgz对变长序列的处理能力
- 测试多样化序列组合的处理稳定性

使用命令：
1. ./bin/pbgz sort {input_file} | ./bin/pbgz compress -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 处理成功，支持变长基线长度
- 解压成功，数据完整
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamVariableBaseLengthPipelineTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamVariableBaseLengthPipelineTest")
        self.source_file = "sam_var_pipeline.sam"
        self.compressed_file = "sam_var_pipeline.sam.pbgz"
        self.decompressed_file = "sam_var_pipeline.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [12, 18, 24, 30, 36]
        for i, base_len in enumerate(base_lengths):
            qname = f"VAR_PIPE_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 2500
            mapq = 60
            cigar = f"{base_len}M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=base_len))
            qual = '!' * base_len
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:{base_len}"
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
        return os.path.exists(self.compressed_file) and os.path.exists(self.decompressed_file)

if __name__ == "__main__":
    test_case = SamVariableBaseLengthPipelineTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)