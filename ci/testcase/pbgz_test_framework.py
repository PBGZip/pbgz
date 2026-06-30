import os
import sys
import subprocess
import hashlib
import json
import re
import time
import shutil
from typing import List, Dict, Any, Optional
from abc import ABC, abstractmethod
import tempfile


class PBGZTestCase(ABC):
    """PBGZ压缩/解压缩测试用例的抽象基类
    
    提供测试框架的基础功能，包括：
    - 命令执行和结果收集
    - 性能指标（执行时间、压缩率）
    - 错误处理和报告
    - 文件清理
    """
    
    def __init__(self, test_name: str):
        self.test_name = test_name
        self.test_commands = []
        self.test_results = {
            "test_name": test_name,
            "passed": True,
            "error": None,
            "commands": []
        }
        self.execution_metrics = {
            "total_time": 0.0,
            "compression_ratios": [],
            "average_compression_ratio": 0.0
        }
        self.temp_files = []
        
    @abstractmethod
    def prepare_data(self):
        """准备测试数据
        
        子类必须实现此方法来：
        - 创建/生成测试文件
        - 设置测试命令
        - 定义压缩/解压缩流程
        """
        pass
    
    @abstractmethod
    def verify_expected_results(self) -> bool:
        """验证测试结果
        
        子类必须实现此方法来验证：
        - MD5校验
        - 预期的压缩率
        - 其他自定义验证逻辑
        
        Returns:
            bool: 如果所有验证通过返回True，否则返回False
        """
        pass
    
    def cleanup_test_data(self):
        """清理测试数据
        
        子类应该重写此方法来清理特定的测试数据文件。
        这是推荐的做法，确保测试用例清理自己创建的所有临时文件。
        """
        pass
    
    def get_test_files(self) -> tuple:
        """获取测试文件信息
        
        子类可以重写此方法来返回测试中使用的文件信息。
        
        Returns:
            tuple: (source_file, compressed_file, decompressed_file)
                   - source_file: 原始源文件名
                   - compressed_file: 压缩后文件名
                   - decompressed_file: 解压后文件名（如果有的话）
                   如果任何文件不存在，返回None
        """
        return (None, None, None)
    
    def _cleanup_framework_files(self):
        """清理框架内部产生的临时文件
        
        注意：这个方法只清理框架自身产生的文件，不清理测试用例的临时文件。
        测试用例的临时文件应该通过重写cleanup_test_data()方法来清理。
        
        关于JSON结果文件的清理策略：
        - 测试成功的用例：清理对应的 *_results.json 文件（结果正确，无需保留）
        - 测试失败的用例：保留对应的 *_results.json 文件（用于调试分析）
        - 测试套件总结文件 test_suite_results.json：保留（包含所有测试的汇总信息）
        """
        # 根据测试结果决定是否清理JSON结果文件
        json_file = f"{self.test_name}_results.json"
        if self.test_results.get("passed", False):
            # 测试成功，清理JSON结果文件
            try:
                if os.path.exists(json_file):
                    os.remove(json_file)
            except Exception as e:
                print(f"Warning: Failed to remove framework JSON file {json_file}: {e}")
        else:
            # 测试失败，保留JSON结果文件用于调试
            pass
        
        # 清理框架内部维护的临时文件列表
        # 此列表由框架内部根据需要添加
        for filepath in self.temp_files:
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                print(f"Warning: Failed to remove framework temp file {filepath}: {e}")
    
    def _get_file_md5(self, filepath: str) -> Optional[str]:
        """计算文件的MD5哈希值
        
        Args:
            filepath: 文件路径
            
        Returns:
            str: MD5哈希值（小写十六进制），如果失败返回None
        """
        try:
            md5_hash = hashlib.md5()
            with open(filepath, 'rb') as f:
                for chunk in iter(lambda: f.read(8192), b''):
                    md5_hash.update(chunk)
            return md5_hash.hexdigest()
        except Exception as e:
            print(f"Error calculating MD5 for {filepath}: {e}")
            return None
    
    def _calculate_compression_ratio(self, original_size: int, compressed_size: int) -> float:
        """计算压缩比
        
        Args:
            original_size: 原始文件大小
            compressed_size: 压缩后文件大小
            
        Returns:
            float: 压缩比（压缩后大小/原始大小的百分比，越小越好）
        """
        if original_size == 0:
            return 0.0
        return (compressed_size / original_size) * 100
    
    def _run_command(self, command: str) -> Dict[str, Any]:
        """执行命令并返回结果
        
        Args:
            command: 要执行的命令
            
        Returns:
            dict: 包含执行结果的字典
        """
        result = {
            "command": command,
            "success": False,
            "exit_code": None,
            "stdout": "",
            "stderr": "",
            "execution_time": 0.0,
            "error_info": {},
            "file_sizes": {},
            "md5_hashes": {}
        }
        
        # 获取测试文件信息
        test_files = self.get_test_files()
        
        # 直接使用get_test_files提供的文件信息
        if test_files and test_files[0]:
            # 构建文件列表用于检查
            files_to_check = [test_files[0]]
            if test_files[1]:  # 如果有压缩文件，也检查
                files_to_check.append(test_files[1])
            if test_files[2]:  # 如果有解压文件，也检查
                files_to_check.append(test_files[2])
            pre_file_sizes = get_file_sizes(files_to_check)
        else:
            pre_file_sizes = {}
        
        start_time = time.time()
        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate()
            exit_code = process.returncode
            
            result["exit_code"] = exit_code
            result["stdout"] = stdout
            result["stderr"] = stderr
            
            # 获取命令执行后的文件大小
            if test_files and test_files[0]:  # 使用子类提供的方法
                files_to_check_post = [
                    f for f in [test_files[0], test_files[1], test_files[2]] if f
                ]
                post_file_sizes = get_file_sizes(files_to_check_post)
            else:
                post_file_sizes = {}
            
            result["file_sizes"] = post_file_sizes
            
            # 专门记录输出文件信息
            output_files = {}
            if "compress" in command and "pbgz" in command:
                # 压缩命令的输出文件
                _, compressed_file, _ = test_files
                if compressed_file and compressed_file in post_file_sizes:
                    output_files[compressed_file] = post_file_sizes[compressed_file]
                    result["output_file"] = compressed_file
                    result["output_file_size"] = post_file_sizes[compressed_file]
                    
            elif "decompress" in command and "pbgz" in command:
                # 解压命令的输出文件
                _, _, decompressed_file = test_files
                if decompressed_file and decompressed_file in post_file_sizes:
                    output_files[decompressed_file] = post_file_sizes[decompressed_file]
                    result["output_file"] = decompressed_file
                    result["output_file_size"] = post_file_sizes[decompressed_file]
                    
            elif "sort" in command and "pbgz" in command:
                # sort命令的输出文件（通常是 .sort.sam 格式）
                # 需要从命令中提取输出文件名
                if "-o" in command:
                    # 解析 -o 后面的文件名
                    parts = command.split()
                    for i, part in enumerate(parts):
                        if part == "-o" and i + 1 < len(parts):
                            output_file = parts[i + 1]
                            if output_file in post_file_sizes:
                                output_files[output_file] = post_file_sizes[output_file]
                                result["output_file"] = output_file
                                result["output_file_size"] = post_file_sizes[output_file]
                                break
                    
            elif "index" in command and "pbgz" in command:
                # 索引命令生成的索引文件
                source_file, compressed_file, _ = test_files
                if compressed_file:
                    index_file = compressed_file + ".index"
                    if os.path.exists(index_file):
                        output_files[index_file] = os.path.getsize(index_file)
                        result["output_file"] = index_file
                        result["output_file_size"] = output_files[index_file]
            
            result["output_files"] = output_files
            
            # 计算压缩率和收集MD5哈希
            if "compress" in command and "pbgz" in command:
                try:
                    # 使用子类提供的接口获取文件信息
                    source_file, compressed_file, _ = self.get_test_files()
                    
                    if source_file and compressed_file:
                        original_size = pre_file_sizes.get(source_file)
                        compressed_size = post_file_sizes.get(compressed_file)
                        
                        if original_size and compressed_size:
                            compression_ratio = self._calculate_compression_ratio(original_size, compressed_size)
                            result["compression_ratio"] = compression_ratio
                            
                            # 收集原始文件的MD5哈希
                            filepath = os.path.join(os.getcwd(), source_file)
                            if os.path.exists(filepath):
                                original_md5 = self._get_file_md5(filepath)
                                if original_md5:
                                    result["md5_hashes"][source_file] = original_md5
                except Exception as e:
                    pass  # 静默失败，不影响主要逻辑
            
            # 如果是解压命令，收集解压文件的MD5哈希
            if "decompress" in command and "pbgz" in command:
                for filename, size in post_file_sizes.items():
                    if not filename.endswith('.pbgz') and size > 0:  # 跳过.pbgz文件
                        filepath = os.path.join(os.getcwd(), filename)
                        if os.path.exists(filepath):
                            decompressed_md5 = self._get_file_md5(filepath)
                            if decompressed_md5:
                                result["md5_hashes"][filename] = decompressed_md5
            
            # 检查命令是否成功
            if exit_code == 0:
                result["success"] = True
            else:
                result["success"] = False
                result["error_info"] = self._parse_error(stderr)
                    
        except Exception as e:
            result["error_info"] = {
                "error_type": type(e).__name__,
                "error_message": str(e)
            }
        finally:
            result["execution_time"] = time.time() - start_time
        
        return result
    
    def _parse_error(self, stderr: str) -> Dict[str, Any]:
        """解析错误信息
        
        Args:
            stderr: 错误输出
            
        Returns:
            dict: 解析后的错误信息
        """
        error_info = {
            "error_type": "Unknown error",
            "error_message": "",
            "error_details": []
        }
        
        if not stderr:
            return error_info
        
        try:
            # 尝试解析常见的错误格式
            # 例如："error: file not found" 或 "Error: Invalid argument"
            error_pattern = r'error:\s*([^,\n]+)'
            match = re.search(error_pattern, stderr, re.IGNORECASE)
            if match:
                error_info["error_type"] = match.group(1).strip()
            
            # 提取第一行错误消息
            clean_stderr = stderr.strip()
            first_line = clean_stderr.split('\n')[0] if clean_stderr else ""
            error_info["error_message"] = first_line
            error_info["error_details"] = clean_stderr.split('\n')[:3]  # 最多3行详细信息
            
        except:
            pass
        
        return error_info
    
    def add_test_command(self, command: str, cmd_id: int = 0, expected_to_fail: bool = False):
        """添加测试命令
        
        Args:
            command: 要执行的命令
            cmd_id: 命令ID，用于排序执行，ID越小越先执行
            expected_to_fail: 是否预期命令失败（如对未排序文件进行索引）
        """
        self.test_commands.append({
            'command': command,
            'cmd_id': cmd_id,
            'expected_to_fail': expected_to_fail
        })
    
    def add_cleanup_command(self, command: str):
        """添加清理命令（已废弃，使用自动清理）"""
        pass
    
    def get_compression_rate(self) -> Optional[float]:
        """获取平均压缩比
        
        Returns:
            float: 平均压缩比，如果没有压缩数据返回None
        """
        if self.execution_metrics["compression_ratios"]:
            return self.execution_metrics["average_compression_ratio"]
        return None
    
    def get_execution_time(self) -> float:
        """获取总执行时间
        
        Returns:
            float: 总执行时间（秒）
        """
        return self.execution_metrics["total_time"]
    
    def _verify_md5_consistency(self, original_file: str, decompressed_file: str) -> bool:
        """验证解压后的文件与原始文件的MD5是否一致
        
        Args:
            original_file: 原始文件名
            decompressed_file: 解压后的文件名
            
        Returns:
            bool: 如果MD5一致返回True，否则返回False
        """
        files_to_clean = []
        try:
            # 获取原始文件的MD5
            for cmd_result in self.test_results["commands"]:
                if original_file in cmd_result.get("md5_hashes", {}):
                    original_md5 = cmd_result["md5_hashes"][original_file]
                    
                    # 获取解压文件的MD5
                    if decompressed_file in cmd_result.get("md5_hashes", {}):
                        decompressed_md5 = cmd_result["md5_hashes"][decompressed_file]
                        
                        # 比较MD5
                        if original_md5 != decompressed_md5:
                            print(f"MD5 mismatch: {original_file} ({original_md5}) != {decompressed_file} ({decompressed_md5})")
                            return False
                        else:
                            return True
            
            return False
        except Exception as e:
            print(f"Error during MD5 verification: {e}")
            return False
    
    def execute(self) -> bool:
        """执行测试用例
        
        Returns:
            bool: 测试是否通过
        """
        try:
            self.prepare_data()
            
            # 按cmd_id排序执行命令
            sorted_commands = sorted(self.test_commands, key=lambda x: x.get('cmd_id', 0))
            
            for cmd_info in sorted_commands:
                result = self._run_command(cmd_info['command'])
                # 记录命令ID和执行结果
                result['cmd_id'] = cmd_info.get('cmd_id', 0)
                self.test_results["commands"].append(result)
                
                # 框架只负责记录结果，不判断成功失败
                execution_time = result.get("execution_time", 0)
                self.execution_metrics["total_time"] += execution_time
                
                # 记录压缩比（仅对压缩命令）
                compression_ratio = result.get("compression_ratio")
                if compression_ratio is not None:
                    self.execution_metrics["compression_ratios"].append(compression_ratio)
            
            # 计算平均压缩比
            if self.execution_metrics["compression_ratios"]:
                self.execution_metrics["average_compression_ratio"] = sum(self.execution_metrics["compression_ratios"]) / len(self.execution_metrics["compression_ratios"])
            
            # 调用子类的验证方法，由子类判断成功失败
            self.test_results["passed"] = self.verify_expected_results()
            
        except Exception as e:
            print(f"Test execution error: {e}")
            self.test_results["passed"] = False
            self.test_results["error"] = str(e)
        finally:
            self.cleanup_test_data()
            self._cleanup_framework_files()
        
        return self.test_results["passed"]
    
    def _check_decompression_consistency(self) -> bool:
        """检查解压缩是否一致
        
        遍历所有命令结果，查找解压失败的证据：
        1. 解压命令返回的stderr包含"Warning: block(0) process failed"
        2. 解压文件大小为0
        3. 解压文件MD5与原始文件不匹配
        
        Returns:
            bool: 如果解压缩一致返回True，否则返回False
        """
        for cmd_result in self.test_results["commands"]:
            if "decompress" in cmd_result["command"]:
                stderr = cmd_result.get("stderr", "")
                
                # 检查解压错误标志
                if "Warning: block(0) process failed" in stderr or "MD5 mismatch" in stderr:
                    print("Decompression failed: MD5 mismatch detected in stderr")
                    return False
                
                # 检查解压文件大小和MD5
                md5_hashes = cmd_result.get("md5_hashes", {})
                if md5_hashes:
                    # 查找原始文件和解压文件的MD5
                    original_md5 = None
                    decompressed_md5 = None
                    
                    for filename, md5_hash in md5_hashes.items():
                        if filename.endswith('.dec') or filename.endswith('_dec') or (not filename.endswith('.pbgz') and 'dec' in filename):
                            decompressed_md5 = md5_hash
                        elif not filename.endswith('.pbgz'):
                            original_md5 = md5_hash
                    
                    # 对比MD5
                    if original_md5 and decompressed_md5:
                        if original_md5 != decompressed_md5:
                            print("Decompression failed: MD5 mismatch!")
                            return False
                
                # 检查解压文件大小 - 只检查解压产生的文件（.dec结尾的文件）
                decompressed_files = [
                    (filename, size) for filename, size in cmd_result.get("file_sizes", {}).items()
                    if filename.endswith('.dec') or filename.endswith('_dec') or 
                       (not filename.endswith('.pbgz') and 'dec' in filename)
                ]
                
                decompressed_size = sum(size for _, size in decompressed_files)
                if decompressed_size == 0:
                    print("Decompression failed: Empty decompressed file detected")
                    return False
        
        return True
    
    def print_results(self):
        """打印测试结果"""
        passed = self.test_results.get("passed", False)
        status = "✅ PASS" if passed else "❌ FAIL"
        
        exec_time = self.execution_metrics["total_time"]
        
        if self.execution_metrics["compression_ratios"]:
            avg_compression = self.execution_metrics["average_compression_ratio"]
            print(f"{status} | {self.test_name} | {exec_time:.2f}s | {avg_compression:.2f}%")
        else:
            print(f"{status} | {self.test_name} | {exec_time:.2f}s | N/A")
    
    def save_to_json(self, filename: Optional[str] = None):
        """保存测试结果到JSON文件
        
        Args:
            filename: 输出文件名，如果为None则使用默认命名
        """
        if filename is None:
            filename = f"{self.test_name}_results.json"
        
        # 添加时间戳
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        self.test_results["timestamp"] = timestamp
        self.test_results["execution_metrics"] = self.execution_metrics
        
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            # 只在测试失败时打印保存消息（因为失败的测试会保留JSON文件）
            if not self.test_results.get("passed", True):
                print(f"Results saved to {filename} (for debugging)")
        except Exception as e:
            print(f"Error saving results to JSON: {e}")


def run_test_suite(test_cases: List[PBGZTestCase]):
    """运行测试套件
    
    Args:
        test_cases: 测试用例列表
    """
    results = []
    
    print(f"开始执行测试套件，共 {len(test_cases)} 个测试用例\n")
    
    for test_case in test_cases:
        print(f"运行测试用例: {test_case.test_name}")
        passed = test_case.execute()
        test_case.print_results()
        print()  # 用空行分隔用例
        
        # 只为失败的测试保存JSON（保留用于调试）
        if not passed:
            test_case.save_to_json()
            
        results.append({
            "test_name": test_case.test_name,
            "passed": passed
        })
    
    # 打印总结
    print(f"测试完成")
    passed_count = sum(1 for r in results if r["passed"])
    print(f"通过: {passed_count}/{len(test_cases)} ({passed_count/len(test_cases)*100:.1f}%)")
    
    # 保存测试套件结果
    suite_results = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_tests": len(test_cases),
        "passed": passed_count,
        "failed": len(test_cases) - passed_count,
        "results": results
    }
    
    try:
        with open("test_suite_results.json", 'w') as f:
            json.dump(suite_results, f, indent=2)
    except Exception as e:
        print(f"Error saving suite results: {e}")
    
    # 智能清理JSON结果文件：成功测试的JSON文件将被清理，失败的保留用于调试
    failed_tests = [r["test_name"] for r in results if not r["passed"]]
    success_cleaned, failure_kept = cleanup_all_json_results(failed_tests)


def get_file_sizes(files_to_check: list = None) -> Dict[str, int]:
    """获取指定文件的大小
    
    Args:
        files_to_check: 要检查的文件名列表
        
    Returns:
        dict: 文件名到文件大小的映射
    """
    sizes = {}
    if files_to_check is None:
        return sizes
    
    try:
        for filename in files_to_check:
            if filename and os.path.exists(filename):
                sizes[filename] = os.path.getsize(filename)
    except:
        pass
    return sizes


def cleanup_all_json_results(failed_tests=None):
    """清理所有测试框架生成的JSON结果文件
    
    此函数清理：
    1. 测试成功的用例的 *_results.json 文件
    2. 保留测试失败的用例的 *_results.json 文件（用于调试）
    3. 保留测试套件结果 test_suite_results.json 文件
    
    Args:
        failed_tests: 失败的测试用例名称列表，对应的JSON将被保留
    
    Returns:
        tuple: (清理的成功JSON数量, 保留的失败JSON数量)
    """
    import glob
    import json
    
    if failed_tests is None:
        failed_tests = []
    
    # 尝试从test_suite_results.json中读取失败的测试信息
    try:
        with open("test_suite_results.json", 'r') as f:
            suite_results = json.load(f)
            for result in suite_results.get("results", []):
                if not result.get("passed", False):
                    test_name = result.get("test_name", "")
                    if test_name not in failed_tests:
                        failed_tests.append(test_name)
    except:
        pass
    
    # 清理成功的测试JSON，保留失败的测试JSON
    json_pattern = '*_results.json'
    success_cleaned = 0
    failure_kept = 0
    
    for filepath in glob.glob(json_pattern):
        try:
            if os.path.exists(filepath):
                # 提取测试名称（去掉_results.json后缀）
                test_name = filepath.replace("_results.json", "")
                
                if test_name in failed_tests:
                    # 失败的测试，保留JSON文件
                    failure_kept += 1
                else:
                    # 成功的测试，清理JSON文件
                    os.remove(filepath)
                    success_cleaned += 1
        except Exception as e:
            print(f"Warning: Failed to process {filepath}: {e}")
    
    return (success_cleaned, failure_kept)