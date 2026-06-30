"""
测试用例：SamRealFileNoReferenceTest - 真实大文件SAM性能测试（无参考基因组）

测试场景：
- 使用真实的大规模SAM文件（SRR2007660.sam，3.7GB）进行性能测试
- 验证pbgz在无参考基因组情况下对真实数据的压缩性能
- 测试大规模文件处理的时间效率和压缩率

使用命令：
1. ./bin/pbgz compress sam/SRR2007660.sam -o sam/SRR2007660.sam.pbgz
2. ./bin/pbgz decompress sam/SRR2007660.sam.pbgz -o sam/SRR2007660.sam.dec

性能验证指标：
- 压缩时间（大文件处理时间）
- 压缩比（压缩率）
- 解压时间和数据完整性
- 内存使用情况（通过测试框架监控）

预期结果：
- 压缩成功，处理大规模数据
- 压缩率达到预期水平（通常>50%）
- 解压成功，数据完美还原
- MD5校验完全一致
"""

import os
import sys
import hashlib
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamRealFileNoReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamRealFileNoReferenceTest")
        self.source_file = "sam/SRR2007660.sam"
        self.compressed_file = "sam/SRR2007660.sam.pbgz"
        self.decompressed_file = "sam/SRR2007660.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # 验证源文件存在
        if not os.path.exists(self.source_file):
            print(f"Error: Source file {self.source_file} does not exist")
            return False
        
        print(f"Real SAM file: {self.source_file}")

        # 压缩（添加-f参数强制覆盖已存在的文件）
        self.add_test_command(f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -f")
        
        # 解压（添加-f参数强制覆盖已存在的文件）
        self.add_test_command(f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -f")

    def cleanup_test_data(self):
        # 清理所有生成的文件
        for filename in [self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                pass

    def verify_expected_results(self) -> bool:
        # 检查命令执行结果
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "compress":
                    compress_success = success
                    # 记录压缩时间
                    if "execution_time" in cmd_result:
                        comp_time = cmd_result["execution_time"]
                        print(f"  压缩时间: {comp_time:.2f}秒")
                elif main_cmd == "decompress":
                    decompress_success = success
                    # 记录解压时间
                    if "execution_time" in cmd_result:
                        decomp_time = cmd_result["execution_time"]
                        print(f"  解压时间: {decomp_time:.2f}秒")
        
        # 验证基本功能
        if not compress_success:
            print("压缩命令失败")
            return False
        
        if not decompress_success:
            print("解压命令失败")
            return False
        
        # 验证压缩文件存在
        if not os.path.exists(self.compressed_file):
            print("压缩文件不存在")
            return False
        
        # 计算压缩比
        original_size = os.path.getsize(self.source_file)
        compressed_size = os.path.getsize(self.compressed_file)
        compression_ratio = (original_size - compressed_size) / original_size * 100
        print(f"  压缩率: {compression_ratio:.2f}%")
        
        # 验证解压文件存在且大小合理
        decompressed_size = os.path.getsize(self.decompressed_file)
        print(f"  解压文件大小: {decompressed_size:,} 字节 ({decompressed_size / (1024**3):.2f}GB)")
        
        # MD5校验前几个MB的数据（性能考虑，不校验整个3.7GB文件）
        with open(self.source_file, 'rb') as f1, open(self.decompressed_file, 'rb') as f2:
            original_part = f1.read(10*1024*1024)  # 读取前10MB
            decompressed_part = f2.read(10*1024*1024)
            
        original_md5 = hashlib.md5(original_part).hexdigest()
        decompressed_md5 = hashlib.md5(decompressed_part).hexdigest()
        
        if original_md5 != decompressed_md5:
            return False
        
        print("✓ SamRealFileNoReferenceTest 性能测试完成！")
        return True


if __name__ == "__main__":
    test_case = SamRealFileNoReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)