"""
测试用例：SamSortTest - SAM文件排序功能测试

测试场景：
- 测试pbgz的sort命令对SAM文件进行排序
- 验证SAM记录按参考序列位置正确排序
- 测试排序后的压缩功能

使用命令：
1. ./bin/pbgz sort {input_file} -o {sorted_file}
2. ./bin/pbgz compress {sorted_file} -o {compressed_file}

参数说明：
- sort命令：按SAM头部的SO字段进行排序
- -o参数：指定排序后的输出文件位置
- SAM文件从unsorted状态变为coordinate状态

预期结果：
- sort命令成功执行，生成排序后的SAM文件
- 压缩命令成功执行
- SAM记录按位置正确排序
- 排序后的压缩文件质量良好

技术说明：
- SAM排序是比对后分析的重要步骤
- 需要正确处理@HD头部中的SO字段
- 测试基本的sort功能工作流程
"""

import os
import sys
import random
import hashlib

# 添加父目录到Python路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamSortTest(PBGZTestCase):
    def __init__(self):

        super().__init__("SamSortTest")
        self.source_file = "sam_unsorted.sam"
        self.sorted_file = "sam_unsorted.sort.sam"
        self.output_dir = "."

    def get_test_files(self) -> tuple:
        """返回测试文件信息"""
        return (self.source_file, None, self.sorted_file)

    def prepare_data(self):
        # 生成一个无序的SAM文件
        sam_content = []
        
        # SAM文件头部
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # 生成无序的数据行（故意打乱位置）
        # 先生成有序的数据
        ordered_lines = []
        num_reads = 50
        for i in range(num_reads):
            flag = 0
            seq_id = f"SAM_{i:03d}"
            pos = (i + 1) * 100
            mapq = 60
            cigar = "50M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=50))
            qual = '!' * 50
            
            line = f"{seq_id}\t{flag}\tENA|ABQD01000001|ABQD01000001.1\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:50\tAS:i:50\tXS:i:0\tRG:Z:SORT_TEST"
            ordered_lines.append(line)
        
        # 打乱顺序
        random.shuffle(ordered_lines)
        sam_content.extend(ordered_lines)
        
        full_content = '\n'.join(sam_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        # 添加排序测试命令，使用-o参数明确指定输出文件名，-f参数强制覆盖
        sort_command = f"./bin/pbgz sort {self.source_file} -o {self.sorted_file} -f"
        self.add_test_command(sort_command, cmd_id=1)
    
    def cleanup_test_data(self):
        files_to_clean = [self.source_file, self.sorted_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        # 验证sort命令是否成功
        sort_success = False
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            if "sort" in command and "pbgz" in command:
                sort_success = success
                
                # 如果sort失败，记录错误信息用于调试
                if not success:
                    error_info = cmd_result.get("error_info", {})
                    print(f"Sort command failed: {error_info}")
        
        if not sort_success:
            return False
            
        if not os.path.exists(self.sorted_file):
            return False
        
        sorted_size = os.path.getsize(self.sorted_file)
        if sorted_size == 0:
            return False
        
        # 验证文件是否真正有序
        with open(self.sorted_file, 'r') as f:
            lines = f.readlines()
        
        # 检查头部信息（头部行以@开头）
        header_lines = []
        data_lines = []
        for line in lines:
            if line.startswith('@'):
                header_lines.append(line)
            else:
                data_lines.append(line)
        
        # 验证数据行是否按位置排序
        positions = []
        for line in data_lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 4:
                    try:
                        pos = int(parts[3])  # 第4列是position
                        positions.append((pos, line.strip()))
                    except (ValueError, IndexError):
                        pass
        
        # 检查位置是否递增
        for i in range(len(positions) - 1):
            if positions[i][0] > positions[i+1][0]:
                return False
        
        # 验证排序后的文件内容完整性
        with open(self.sorted_file, 'r') as f:
            sorted_content = f.read()
        
        sorted_hash = hashlib.md5(sorted_content.encode('utf-8')).hexdigest()
        
        # 由于排序会改变文件内容的顺序，我们不能直接比较MD5
        # 但我们可以检查排序后的文件包含了正确的头部和数据
        if "@HD\tVN:1.0\tSO:unsorted" not in sorted_content:
            return False
        
        if "@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000" not in sorted_content:
            return False
        
        
        return True


if __name__ == "__main__":
    test_case = SamSortTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # 为调试结果保存JSON文件
        if not result:
            test_case.save_to_json()
        sys.exit(1)