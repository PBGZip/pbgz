import os
import sys
import inspect
import importlib.util
from testcase.pbgz_test_framework import PBGZTestCase, run_test_suite

def find_testcase_subclasses(testcase_dir="testcase"):
    """查找所有继承自PBGZTestCase的子类"""
    testcase_classes = []
    
    if not os.path.exists(testcase_dir):
        print(f"Testcase directory not found: {testcase_dir}")
        return testcase_classes
    
    # 遍历testcase目录下的所有Python文件
    for filename in os.listdir(testcase_dir):
        if filename.endswith('.py') and filename != '__init__.py' and filename != 'pbgz_test_framework.py':
            filepath = os.path.join(testcase_dir, filename)
            
            try:
                # 动态导入模块
                module_name = filename[:-3]  # 去掉.py后缀
                file_path = os.path.join(os.getcwd(), filepath)
                spec = importlib.util.spec_from_file_location(module_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                # 查找模块中所有的类
                for name, obj in inspect.getmembers(module, inspect.isclass):
                    # 检查是否是PBGZTestCase的子类，且不是PBGZTestCase本身
                    if issubclass(obj, PBGZTestCase) and obj != PBGZTestCase:
                        # 确保该类是在当前模块中定义的，避免重复
                        if obj.__module__ == module_name or obj.__module__.startswith(module_name):
                            testcase_classes.append(obj)
                            
            except Exception as e:
                print(f"Error loading test case from {filename}: {e}")
    
    return testcase_classes

def run_all_testcases():
    """运行所有测试用例"""
    print("=" * 60)
    print("开始扫描测试用例...")
    print("=" * 60)
    
    # 查找所有测试用例类
    testcase_classes = find_testcase_subclasses()
    
    if not testcase_classes:
        print("No test cases found!")
        return
    
    print(f"Found {len(testcase_classes)} test case(s):")
    for tc_class in testcase_classes:
        print(f"  - {tc_class.__name__}")
    print("=" * 60)
    
    # 实例化所有测试用例
    test_cases = []
    for tc_class in testcase_classes:
        try:
            test_case = tc_class()
            test_cases.append(test_case)
        except Exception as e:
            print(f"Error instantiating {tc_class.__name__}: {e}")
    
    # 运行所有测试用例
    if test_cases:
        run_test_suite(test_cases)

if __name__ == "__main__":
    # 清理之前的JSON文件（如果存在）
    from testcase.pbgz_test_framework import cleanup_all_json_results
    cleanup_all_json_results()
    
    run_all_testcases()