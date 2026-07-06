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
    """Abstract base class for PBGZ compression/decompression test cases

    Provides basic functionality for the test framework, including:
    - Command execution and result collection
    - Performance metrics (execution time, compression ratio)
    - Error handling and reporting
    - File cleanup
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
        """Prepare test data"""
        pass

    @abstractmethod
    def verify_expected_results(self) -> bool:
        """Verify test results

        Returns:
            bool: Return True if all verifications pass, otherwise return False
        """
        pass

    def cleanup_test_data(self):
        """Clean up test data"""
        pass

    def get_test_files(self) -> tuple:
        """Get test file information

        Returns:
            tuple: (source_file, compressed_file, decompressed_file)
        """
        return (None, None, None)

    def _cleanup_framework_files(self):
        """Clean up framework files"""
        json_file = f"{self.test_name}_results.json"
        if self.test_results.get("passed", False):
            try:
                if os.path.exists(json_file):
                    os.remove(json_file)
            except Exception as e:
                print(f"Warning: Failed to remove framework JSON file {json_file}: {e}")

        for filepath in self.temp_files:
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                print(f"Warning: Failed to remove framework temp file {filepath}: {e}")

    def _get_file_md5(self, filepath: str) -> Optional[str]:
        """Calculate file MD5"""
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
        """Calculate compression ratio"""
        if original_size == 0:
            return 0.0
        return (compressed_size / original_size) * 100

    def _run_command(self, command: str) -> Dict[str, Any]:
        """Execute command and return result"""
        result = {
            "command": command,
            "success": False,
            "exit_code": None,
            "stdout": "",
            "stderr": "",
            "execution_time": 0.0,
            "error_info": {},
            "file_sizes": {},
            "md5_hashes": {},
            "compression_ratio": None
        }

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

            if exit_code == 0:
                result["success"] = True
                # Try to parse compression ratio from command output
                result["compression_ratio"] = self._parse_compression_ratio(stdout + stderr)
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
        """Parse error information"""
        error_info = {
            "error_type": "Unknown error",
            "error_message": "",
            "error_details": []
        }

        if not stderr:
            return error_info

        try:
            error_pattern = r'error:\s*([^\,\n]+)'
            match = re.search(error_pattern, stderr, re.IGNORECASE)
            if match:
                error_info["error_type"] = match.group(1).strip()

            clean_stderr = stderr.strip()
            first_line = clean_stderr.split('\n')[0] if clean_stderr else ""
            error_info["error_message"] = first_line
            error_info["error_details"] = clean_stderr.split('\n')[:3]

        except:
            pass

        return error_info

    def _parse_compression_ratio(self, output: str) -> Optional[float]:
        """Parse compression ratio from command output"""
        try:
            # Look for pattern like "ratio 100.15%" or "ratio: 100.15%"
            ratio_pattern = r'ratio\s*[:\s]*([0-9]*\.?[0-9]+)%'
            match = re.search(ratio_pattern, output)
            if match:
                ratio = float(match.group(1))
                self.execution_metrics["compression_ratios"].append(ratio)
                return ratio
            return None
        except Exception as e:
            return None

    def add_test_command(self, command: str, cmd_id: int = 0, expected_to_fail: bool = False):
        """Add test command"""
        self.test_commands.append({
            'command': command,
            'cmd_id': cmd_id,
            'expected_to_fail': expected_to_fail
        })

    def add_cleanup_command(self, command: str):
        """Add cleanup command (deprecated)"""
        pass

    def get_compression_rate(self) -> Optional[float]:
        """Get average compression ratio"""
        if self.execution_metrics["compression_ratios"]:
            return self.execution_metrics["average_compression_ratio"]
        return None

    def get_execution_time(self) -> float:
        """Get total execution time"""
        return self.execution_metrics["total_time"]

    def execute(self) -> bool:
        """Execute test case"""
        try:
            self.prepare_data()

            sorted_commands = sorted(self.test_commands, key=lambda x: x.get('cmd_id', 0))

            for cmd_info in sorted_commands:
                result = self._run_command(cmd_info['command'])
                result['cmd_id'] = cmd_info.get('cmd_id', 0)
                self.test_results["commands"].append(result)

                execution_time = result.get("execution_time", 0)
                self.execution_metrics["total_time"] += execution_time

            if self.execution_metrics["compression_ratios"]:
                self.execution_metrics["average_compression_ratio"] = sum(self.execution_metrics["compression_ratios"]) / len(self.execution_metrics["compression_ratios"])

            self.test_results["passed"] = self.verify_expected_results()

        except Exception as e:
            print(f"Test execution error: {e}")
            self.test_results["passed"] = False
            self.test_results["error"] = str(e)
        finally:
            self.cleanup_test_data()
            self._cleanup_framework_files()

        return self.test_results["passed"]

    def print_results(self):
        """Print test results"""
        passed = self.test_results.get("passed", False)
        status = "PASS" if passed else "FAIL"

        exec_time = self.execution_metrics["total_time"]

        if self.execution_metrics["compression_ratios"]:
            avg_compression = self.execution_metrics["average_compression_ratio"]
            print(f"{status} | {self.test_name} | {exec_time:.2f}s | {avg_compression:.2f}%")
        else:
            print(f"{status} | {self.test_name} | {exec_time:.2f}s | N/A")

    def save_to_json(self, filename: Optional[str] = None):
        """Save test results to JSON file"""
        if filename is None:
            filename = f"{self.test_name}_results.json"

        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        self.test_results["timestamp"] = timestamp
        self.test_results["execution_metrics"] = self.execution_metrics

        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            if not self.test_results.get("passed", True):
                print(f"Results saved to {filename} (for debugging)")
        except Exception as e:
            print(f"Error saving results to JSON: {e}")


def run_test_suite(test_cases: List[PBGZTestCase]):
    """Run test suite"""
    results = []

    print(f"Starting test suite execution, total {len(test_cases)} test cases\n")

    for test_case in test_cases:
        print(f"Running test case: {test_case.test_name}")
        passed = test_case.execute()
        test_case.print_results()
        print()

        if not passed:
            test_case.save_to_json()

        results.append({
            "test_name": test_case.test_name,
            "passed": passed
        })

    print(f"Test completed")
    passed_count = sum(1 for r in results if r["passed"])
    print(f"Passed: {passed_count}/{len(test_cases)} ({passed_count/len(test_cases)*100:.1f}%)")

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

    failed_tests = [r["test_name"] for r in results if not r["passed"]]
    success_cleaned, failure_kept = cleanup_all_json_results(failed_tests)


def get_file_sizes(files_to_check: list = None) -> Dict[str, int]:
    """Get file sizes"""
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
    """Clean up JSON result files"""
    import glob
    import json

    if failed_tests is None:
        failed_tests = []

    json_pattern = '*_results.json'
    success_cleaned = 0
    failure_kept = 0

    for filepath in glob.glob(json_pattern):
        try:
            if os.path.exists(filepath):
                test_name = filepath.replace("_results.json", "")

                if test_name in failed_tests:
                    failure_kept += 1
                else:
                    os.remove(filepath)
                    success_cleaned += 1
        except Exception as e:
            print(f"Warning: Failed to process {filepath}: {e}")

    return (success_cleaned, failure_kept)