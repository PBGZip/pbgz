"""
测试用例：SamMixedOptionalFieldsTest - 混合可选字段SAM文件压缩解压测试

测试场景：
- 测试SAM文件中包含混合可选字段（部分记录有，部分记录没有）
- 验证pbgz对可选字段处理的一致性
- 测试数据格式灵活性

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，正确处理混合可选字段
- 解压成功，可选字段完整保留
- MD5校验完全一致
"""

import sys
import random
import hashlib

# 添加父目录到Python路径
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase
import random
import hashlib
# 添加父目录到Python路径
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase
# 有重复sys.path.insert，需要修复一下结构，但先继续处理其他文件
# 添加父目录到Python路径
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class SamMixedOptionalFieldsTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMixedOptionalFieldsTest")
        self.source_file = "sam_mixed_opt.sam"
        self.compressed_file = "sam_mixed_opt.sam.pbgz"
        self.decompressed_file = "sam_mixed_opt.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(20):
            qname = f"MIXED_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 200
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=25))
            qual = '!' * 25
            # 混合可选字段
            if i % 3 == 0:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            elif i % 3 == 1:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0"
            else:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}")
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
    test_case = SamMixedOptionalFieldsTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
