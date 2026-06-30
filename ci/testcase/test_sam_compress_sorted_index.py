"""
测试用例：SamCompressSortedIndexTest - 已排序SAM文件索引创建测试

测试场景：
- 测试对已排序的SAM文件创建索引文件
- 验证索引文件的正确生成和格式
- 测试已排序文件的索引处理（SO:coordinate）

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz index -f {compressed_file}

预期结果：
- 压缩命令成功执行
- 索引命令成功执行
- 索引文件存在且不为空
- 索引文件格式正确：每行按tab分割，有5个字段

技术说明：
- 已排序的SAM文件可以创建索引
- 索引文件格式为tab分隔的5字段格式
- -f参数强制覆盖现有索引文件
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamCompressSortedIndexTest(PBGZTestCase):
    """测试已排序SAM文件的索引创建"""
    
    def __init__(self):
        super().__init__("SamCompressSortedIndexTest")
        self.source_file = "sam_sorted.sam"
        self.compressed_file = "sam_sorted.sam.pbgz"
        self.index_file = "sam_sorted.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        # 生成已排序的SAM文件
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # 注意这里是coordinate（已排序）
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        
        # 按位置排序的reads
        for i in range(10):
            qname = f"SORTED_IDX_{i:03d}"
            flag = 0
            rname = "ref1"
            pos = (i + 1) * 100  # 位置递增，已排序
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
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}", cmd_id=1)
        self.add_test_command(f"./bin/pbgz index -f {self.compressed_file}", cmd_id=2)

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # 检查命令执行结果
        compress_success = index_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # 先检查"index"避免被"compress"包含
            if "index" in command and "pbgz" in command:
                index_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        # 已排序SAM文件：压缩和索引都应该成功
        if not compress_success:
            return False
            
        # 索引应该成功（因为是已排序文件）
        if not index_success:
            return False
        
        # 验证索引文件内容不为空
        if not os.path.exists(self.index_file):
            print(f"索引文件 {self.index_file} 未生成")
            return False
        
        if os.path.getsize(self.index_file) == 0:
            print(f"索引文件 {self.index_file} 为空")
            return False
        
        # 验证索引文件内容格式：每行按照tab分割，有5个字段
        with open(self.index_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # 跳过空行
                    continue
                
                fields = line.split('\t')
                if len(fields) != 5:
                    print(f"索引文件第{line_num}行字段数量不正确：期望5个，实际{len(fields)}个")
                    print(f"  内容: {line}")
                    return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressSortedIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)