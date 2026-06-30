"""
测试用例：SamHeaderOnlyTest - 纯Header SAM文件压缩解压测试

测试场景：
- 测试只包含header、没有实际记录的SAM文件处理
- 验证pbgz对空数据体SAM文件的处理能力
- 测试文件边界情况和最小数据集处理

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file}

参数说明：
- SAM header以@开头的元数据行
- 测试只有header、无sequence records的极端情况
- 验证最小数据集的压缩解压稳定性

预期结果：
- 压缩成功，处理只有header的文件
- 解压成功，header完整保留
- MD5校验完全一致
"""

import os
import sys
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class SamHeaderOnlyTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamHeaderOnlyTest")
        self.source_file = "sam_header_only.sam"
        self.compressed_file = "sam_header_only.sam.pbgz"
        self.decompressed_file = "sam_header_only.sam.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成只有头部的SAM文件
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
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
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        if original_content == decompressed_content:
            
            return True
        else:
            print(f"✗ Header-only SAM file test failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = SamHeaderOnlyTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)