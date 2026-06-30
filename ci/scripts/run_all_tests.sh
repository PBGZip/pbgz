#!/bin/bash
# Script: run_all_tests.sh
# Description: Run complete test suite including C++ unit tests and Python test cases
# Usage: ./run_all_tests.sh [--pbgz-test-path=PATH] [--python-test-path=PATH]

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default paths
DEFAULT_PBGZ_TEST_PATH="./build-x86_64-gcc-release/test/pbgz_test"
DEFAULT_PYTHON_TEST_PATH="./ci/run_all_testcase.py"

# Parse command line arguments
PBGZ_TEST_PATH="$DEFAULT_PBGZ_TEST_PATH"
PYTHON_TEST_PATH="$DEFAULT_PYTHON_TEST_PATH"

for arg in "$@"; do
    case $arg in
        --pbgz-test-path=*)
            PBGZ_TEST_PATH="${arg#*=}"
            shift
            ;;
        --python-test-path=*)
            PYTHON_TEST_PATH="${arg#*=}"
            shift
            ;;
        *)
            log_warning "Unknown argument: $arg"
            ;;
    esac
done

# Track overall test results
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

log_info "================================================"
log_info "Complete Test Suite Execution"
log_info "================================================"
log_info "Starting at: $(date +'%Y-%m-%d %H:%M:%S')"
log_info "Test environment: $(uname -s) $(uname -m)"
log_info "C++ test path: $PBGZ_TEST_PATH"
log_info "Python test path: $PYTHON_TEST_PATH"
log_info "================================================"
echo ""

# Function to run C++ unit tests
run_cpp_tests() {
    log_info "================================================"
    log_info "Running C++ Unit Tests (GoogleTest)"
    log_info "================================================"
    
    if [ ! -f "$PBGZ_TEST_PATH" ]; then
        log_error "C++ test executable not found: $PBGZ_TEST_PATH"
        log_error "Please build the project first using: ./build-release.sh"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        return 1
    fi
    
    local test_results_file="/tmp/cpp_test_results_$$"
    local test_log_file="/tmp/cpp_test_log_$$"
    
    log_info "Executing: $PBGZ_TEST_PATH"
    
    # Capture test results
    if "$PBGZ_TEST_PATH" > "$test_results_file" 2>"$test_log_file"; then
        log_success "C++ unit tests PASSED"
        
        # Parse test results
        if grep -q "tests from" "$test_results_file"; then
            local test_count=$(grep -o "[0-9]* tests from" "$test_results_file" | grep -o "[0-9]*" | head -1)
            local passed_count=$(grep -o "[0-9]* test" "$test_results_file" | head -1)
            
            if [ -n "$test_count" ]; then
                TOTAL_TESTS=$((TOTAL_TESTS + test_count))
                PASSED_TESTS=$((PASSED_TESTS + test_count))
                log_success "Ran $test_count C++ unit tests, all passed"
            fi
        fi
        
        # Show summary if available
        if grep -q "test" "$test_results_file"; then
            log_info "C++ test summary:"
            tail -10 "$test_results_file" | while IFS= read -r line; do
                if [[ "$line" =~ test ]] || [[ "$line" =~ PASSED ]] || [[ "$line" =~ FAILED ]]; then
                    echo "  $line"
                fi
            done
        fi
    else
        log_error "C++ unit tests FAILED"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        
        # Show errors
        if [ -s "$test_log_file" ]; then
            log_error "Error output:"
            cat "$test_log_file"
        fi
    fi
    
    # Cleanup
    rm -f "$test_results_file" "$test_log_file"
    
    log_info "================================================"
    echo ""
}

# Function to run Python test cases
run_python_tests() {
    log_info "================================================"
    log_info "Running Python Test Cases"
    log_info "================================================"
    
    # Check Python availability
    if ! command -v python3 &> /dev/null; then
        log_error "python3 not found. Please install Python 3."
        FAILED_TESTS=$((FAILED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        return 1
    fi
    
    local python_version=$(python3 --version 2>&1)
    log_info "Python version: $python_version"
    
    if [ ! -f "$PYTHON_TEST_PATH" ]; then
        log_error "Python test script not found: $PYTHON_TEST_PATH"
        log_error "Current directory: $(pwd)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        return 1
    fi
    
    log_info "Executing: python3 $PYTHON_TEST_PATH"
    
    # Run Python tests and capture output
    local python_output
    local python_exit_code
    
    python_output=$(python3 "$PYTHON_TEST_PATH" 2>&1) && python_exit_code=0 || python_exit_code=$?
    
    # Parse Python test results
    local python_test_count=0
    local python_passed_count=0
    local python_failed_count=0
    
    if [ $python_exit_code -eq 0 ]; then
        log_success "Python test cases PASSED"
        
        # Try to parse test summary from output
        if echo "$python_output" | grep -q "Found.*test"; then
            local found_tests=$(echo "$python_output" | grep "Found.*test" | head -1)
            if [[ "$found_tests" =~ ([0-9]+) ]]; then
                python_test_count=${BASH_REMATCH[1]}
            fi
        fi
        
        if [ -n "$python_output" ]; then
            log_info "Python test output:"
            echo "$python_output" | head -20
        fi
        
        TOTAL_TESTS=$((TOTAL_TESTS + python_test_count))
        PASSED_TESTS=$((PASSED_TESTS + python_test_count))
    else
        log_error "Python test cases FAILED (exit code: $python_exit_code)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        
        if [ -n "$python_output" ]; then
            log_error "Python test error output:"
            echo "$python_output"
        fi
    fi
    
    log_info "================================================"
    echo ""
}

# Generate final report
generate_report() {
    log_info "================================================"
    log_info "Final Test Report"
    log_info "================================================"
    log_info "Total test suites: $((TOTAL_TESTS > 0 ? 2 : 0))"
    log_info "Total tests executed: $TOTAL_TESTS"
    log_info "Tests passed: $PASSED_TESTS"
    log_info "Tests failed: $FAILED_TESTS"
    
    if [ $TOTAL_TESTS -gt 0 ]; then
        local pass_rate=$((PASSED_TESTS * 100 / TOTAL_TESTS))
        log_info "Pass rate: ${pass_rate}%"
    fi
    
    log_info "Completion time: $(date +'%Y-%m-%d %H:%M:%S')"
    log_info "================================================"
    
    if [ $FAILED_TESTS -eq 0 ] && [ $TOTAL_TESTS -gt 0 ]; then
        log_success "All tests PASSED! 🎉"
        return 0
    elif [ $TOTAL_TESTS -eq 0 ]; then
        log_warning "No tests were executed"
        return 1
    else
        log_error "Some tests FAILED"
        return 1
    fi
}

# Main execution
main() {
    cd "$SCRIPT_DIR"
    
    # Run C++ tests
    run_cpp_tests
    
    # Run Python tests
    run_python_tests
    
    # Generate report
    generate_report
}

# Execute main function
main "$@"