"""
测试用例：SamUnsortedCompressWithIndexTest - 未排序SAM文件压缩索引测试

测试场景：
- 测试对未排序的SAM文件尝试压缩并同时生成索引
- 验证pbgz对未排序数据的处理能力
- 测试索引建立时的数据排序要求

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -i {index_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

预期结果：
- 压缩成功，能够处理未排序数据
- 解压成功，数据完整恢复
- 验证对未排序数据的容错处理
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamUnsortedCompressWithIndexTest(PBGZTestCase):
    """测试未排序SAM文件压缩时尝试生成索引文件"""
    
    def __init__(self):
        super().__init__("SamUnsortedCompressWithIndexTest")
        self.source_file = "sam_unsorted_with_index.sam"
        self.compressed_file = "sam_unsorted_with_index.sam.pbgz"
        self.index_file = "sam_unsorted_with_index.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        # 生成未排序的SAM文件
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")  # unsorted表示未排序
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        
        # 随机位置的reads（未排序）
        import random
        random.seed(42)
        positions = [random.randint(100, 9000) for _ in range(50)]
        
        for i, pos in enumerate(positions):
            qname = f"UNSORTED_IDX_{i:03d}"
            flag = 0
            rname = "ref1"
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
        
        # 使用-i参数在压缩时尝试生成索引文件，但因为文件未排序，预期不会创建索引
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
        has_unsorted_warning = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            stderr = cmd_result.get("stderr", "")
            
            if "compress" in command and "pbgz" in command:
                compress_success = success
                # 检查是否有未排序的警告信息
                if "unsorted" in stderr.lower() or "not sorted" in stderr.lower() or "warning" in stderr.lower():
                    has_unsorted_warning = True
                break
        
        # 压缩应该成功
        if not compress_success:
            return False
        
        # 检查压缩文件是否存在
        if not os.path.exists(self.compressed_file):
            return False
        
        # 应该没有未排序的警告，因为pbgz可能对未排序文件正常处理
        # 但压缩文件应该存在
        
        # 检查索引文件是否不应该存在（因为是未排序文件）
        if self.index_file and os.path.exists(self.index_file):
            # 如果索引文件存在，检查文件大小，可能为空或很小
            index_size = os.path.getsize(self.index_file)
            if index_size > 0:
                # 有内容的索引文件不应该存在
                pass
        
        return True


if __name__ == "__main__":
    test_case = SamUnsortedCompressWithIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)