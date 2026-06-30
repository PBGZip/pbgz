"""
测试用例：SamPipelineTest - SAM文件管道压缩解压测试

测试场景：
- 测试SAM文件通过Unix管道方式进行压缩解压处理
- 验证SAM格式数据与标准Unix工具的兼容性
- 测试SAM格式验证在管道中的正确执行

使用命令：
1. cat {source_file} | ./bin/pbgz compress -o {compressed_file}
2. ./bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

预期结果：
- 管道压缩解压成功
- SAM格式完整保留（@header、11列数据）
- 解压后文件MD5与原始文件完全一致
- 管道模式功能正常工作
"""

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class SamPipelineTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamPipelineTest")
        self.source_file = "sam_pipe_test.sam"
        self.compressed_file = "sam_pipe_test.sam.pbgz"
        self.decompressed_file = "sam_pipe_test.sam.dec"

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 生成基本的SAM文件内容
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # 添加3条read记录
        num_reads = 3
        for i in range(num_reads):
            qname = f"SAM_PIPE_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATGCN', k=25))
            qual = '!' * 25
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25\tXS:i:0\tRG:Z:PIPE"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # 添加压缩和解压命令（使用管道）
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o - > {self.compressed_file}")
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
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            
            return True
        else:
            print(f"✗ Pipeline test failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = SamPipelineTest()
    result = test_case.execute()
    test_case.print_results()
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)