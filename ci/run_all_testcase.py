import os
import sys
import inspect
import importlib.util

from testcase.pbgz_test_framework import PBGZTestCase, run_test_suite

def find_testcase_subclasses(testcase_dir="testcase"):
    """Find all subclasses inheriting from PBGZTestCase"""
    testcase_classes = []
    
    if not os.path.exists(testcase_dir):
        print(f"Test case directory not found: {testcase_dir}")
        return testcase_classes
    
    # Traverse all Python files in testcase directory
    for filename in os.listdir(testcase_dir):
        if filename.endswith('.py') and filename != '__init__.py' and filename != 'pbgz_test_framework.py':
            filepath = os.path.join(testcase_dir, filename)
            try:
                # Dynamically import module
                module_name = filename[:-3]  # Remove .py extension
                file_path = os.path.join(os.getcwd(), filepath)
                spec = importlib.util.spec_from_file_location(module_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                # Find all classes in module
                for name, obj in inspect.getmembers(module, inspect.isclass):
                    # Check if it is a subclass of PBGZTestCase and not PBGZTestCase itself
                    if issubclass(obj, PBGZTestCase) and obj != PBGZTestCase:
                        # Ensure the class is defined in the current module to avoid duplication
                        if obj.__module__ == module_name or obj.__module__.startswith(module_name):
                            testcase_classes.append(obj)
            except Exception as e:
                print(f"Error loading test case from {filename}: {e}")
    
    return testcase_classes

def run_all_testcases():
    """Run all test cases"""
    print("=" * 60)
    print("Start scanning test cases...")
    print("=" * 60)
    
    # Find testcase classes
    testcase_classes = find_testcase_subclasses()
    
    if not testcase_classes:
        print("No test cases found!")
        return
    
    print(f"Found {len(testcase_classes)} test case(s):")
    for tc_class in testcase_classes:
        print(f" - {tc_class.__name__}")
    
    print("=" * 60)
    
    # Instantiate all test cases
    test_cases = []
    for tc_class in testcase_classes:
        try:
            test_case = tc_class()
            test_cases.append(test_case)
        except Exception as e:
            print(f"Error instantiating {tc_class.__name__}: {e}")
    
    # Run all test cases
    if test_cases:
        run_test_suite(test_cases)

if __name__ == "__main__":
    # Clean up previous JSON files (if exist)
    from testcase.pbgz_test_framework import cleanup_all_json_results
    cleanup_all_json_results()
    
    run_all_testcases()
