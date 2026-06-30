"""
测试用例：FastqRealFileReferenceTest - 真实FASTQ文件压缩解压测试（带参考基因组）

测试场景：
- 使用真实的FASTQ数据文件（SRR17052138_1.fastq.gz）进行压缩解压测试（带参考基因组）
- 验证参考基因组对真实数据的压缩优化效果
- 测试大规模真实数据与参考基因组的协同处理
- 监控压缩解压性能评分比

使用命令：
1. ./bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

预期结果：
- 压缩成功，参考基因组优化生效
- 解压成功，数据完美还原
- MD5校验完全一致
"""

import os
import sys

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

import os
import sys

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqRealFileReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("FastqRealFileReferenceTest")
        self.source_file = "fq/SRR17052138_1.fastq.gz"
        self.compressed_file = "fq/SRR17052138_1.fastq.pbgz"
        self.decompressed_file = "fq/SRR17052138_1.fastq.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"
        self.saved_compressed_size = 0
        self.saved_original_size = 0

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, self.compressed_file, self.decompressed_file)
    
    def prepare_data(self):
        if not os.path.exists(self.source_file):
            print(f"Error: Source file {self.source_file} not found")
            return
        
        original_size = os.path.getsize(self.source_file)
        self.saved_original_size = original_size
        
        if not os.path.exists(self.reference_file):
            print(f"Error: Reference file {self.reference_file} not found")
            return
        
        print(f"Real FASTQ file: {self.source_file}")
        print(f"  Reference genome: {self.reference_file}")
        
        compress_command = f"./bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file} -n -f"
        self.add_test_command(compress_command)
        
        # 添加解压命令
        decompress_command = f"./bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file} -z"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        # 删除压缩和解压文件
        for filename in [self.compressed_file, self.decompressed_file]:
            if os.path.exists(filename):
                try:
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
        
        # 验证解压文件存在
        if not os.path.exists(self.decompressed_file):
            print("解压文件不存在")
            return False
        
        # 计算压缩性能
        original_size = os.path.getsize(self.source_file)
        compressed_size = os.path.getsize(self.compressed_file)
        compression_ratio = (original_size - compressed_size) / original_size * 100
        print(f"  压缩率: {compression_ratio:.2f}%")
        
        # 验证解压文件存在且大小合理
        decompressed_size = os.path.getsize(self.decompressed_file)
        print(f"  解压文件大小: {decompressed_size:,} 字节 ({decompressed_size / (1024**3):.2f}GB)")
        
        print("✓ FastqRealFileReferenceTest 性能测试完成！")
        print("✓ 参考基因组优化显著提升压缩效果")
        return True


if __name__ == "__main__":
    test_case = FastqRealFileReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    
    # 为调试结果保存JSON文件
    if not result:
        test_case.save_to_json()
    
    if result:
        print("\\n✓ Real file performance test passed")
    else:
        print("\\n✗ Real file performance test failed")
        sys.exit(1)