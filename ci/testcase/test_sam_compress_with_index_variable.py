"""
测试用例：SamCompressWithIndexVariableTest - 有序变长SAM文件索引生成测试

测试场景：
- 测试包含变长序列的有序SAM文件能够成功生成索引文件
- 验证有序数据与索引创建的协同工作
- 测试变量长度序列的索引处理

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file}
2. ./bin/pbgz index -f {compressed_file}

预期结果：
- 压缩成功
- 索引文件成功生成且不为空
- 索引文件格式正确（每行5个tab分隔字段）
"""

import sys
import random
import hashlib

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamCompressWithIndexVariableTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamCompressWithIndexVariableTest")
        self.source_file = "sam_index_var.sam"
        self.compressed_file = "sam_index_var.sam.pbgz"
        self.index_file = "sam_index_var.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # 修改为有序
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [15, 20, 25, 30, 35]
        for i, base_len in enumerate(base_lengths):
            qname = f"INDEX_VAR_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 700  # 位置递增，确保有序
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
        
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        # 有序文件应该能够成功生成索引文件
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
            
            # 基于命令类型精确判断，避免字符串包含问题
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "index":
                    index_success = success
                elif main_cmd == "compress":
                    compress_success = success
        
        # 压缩必须成功
        if not compress_success:
            print("Compression command failed")
            return False
        
        # 压缩文件必须存在
        if not os.path.exists(self.compressed_file):
            print(f"Compressed file {self.compressed_file} not found")
            return False
        
        # 对于有序SAM文件，我们期望索引命令成功
        if not index_success:
            print("Index command failed - ordered SAM should be able to generate index")
            return False
        
        # 索引文件必须存在且不为空
        if not os.path.exists(self.index_file):
            print(f"Index file {self.index_file} not created")
            return False
        
        if os.path.getsize(self.index_file) == 0:
            print(f"Index file {self.index_file} is empty")
            return False
        
        # 验证索引文件格式：每行5个tab分隔字段
        try:
            with open(self.index_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) != 5:
                        print(f"Index file line {line_num} has incorrect field count: expected 5, got {len(fields)}")
                        print(f"  Content: {line}")
                        return False
        except Exception as e:
            print(f"Error reading index file: {e}")
            return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressWithIndexVariableTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
