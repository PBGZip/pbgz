#!/bin/bash
# Script: record_build_results.sh
# Description: Record CI build and test results for scheduled builds
# Usage: ./record_build_results.sh --build-status=STATUS --test-status=STATUS --test-results=FILE

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

log_error() {
        echo -e "${RED}[ERROR]${NC} $1"
}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default values
BUILD_STATUS=""
TEST_STATUS=""
BUILD_START_TIME=""
BUILD_END_TIME=""
TEST_RESULTS_FILE=""
PY_TEST_RESULTS_FILE=""
BUILD_LOG_FILE=""
REPORT_OUTPUT_FILE="$SCRIPT_DIR/build_results_report.json"
REPORT_TEXT_FILE="$SCRIPT_DIR/build_results_report.txt"

# Parse command line arguments
for arg in "$@"; do
      case $arg in
          --build-status=*)
          BUILD_STATUS="${arg#*=}"
          shift
          ;;
          --test-status=*)
          TEST_STATUS="${arg#*=}"
          shift
          ;;
          --test-results=*)
          TEST_RESULTS_FILE="${arg#*=}"
          shift
          ;;
          --py-test-results=*)
          PY_TEST_RESULTS_FILE="${arg#*=}"
          shift
          ;;
          --build-start=*)
          BUILD_START_TIME="${arg#*=}"
          shift
          ;;
          --build-end=*)
          BUILD_END_TIME="${arg#*=}"
          shift
          ;;
          --build-log=*)
          BUILD_LOG_FILE="${arg#*=}"
          shift
          ;;
          --output=*)
          REPORT_OUTPUT_FILE="${arg#*=}"
          shift
          ;;
          *)
          log_info "Unknown argument: $arg"
          ;;
      esac
done

# Set default timestamps if not provided
if [ -z "$BUILD_START_TIME" ]; then
      BUILD_START_TIME=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
fi

if [ -z "$BUILD_END_TIME" ]; then
      BUILD_END_TIME=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
fi

# Function to parse test results from log file
parse_test_results() {
      local log_file="$1"
      if [ ! -f "$log_file" ]; then
          echo "0"
          echo "0"
          echo "0"
          return
      fi
      
      # Try to parse gtest output format
      if grep -q "tests.*from" "$log_file" 2>/dev/null; then
          total_tests=$(grep -o "[0-9]* tests from" "$log_file" 2>/dev/null | head -1 | grep -o "[0-9]*" || echo "0")
          passed_tests=$(grep -o "[0-9]* test" "$log_file" 2>/dev/null | head -1 | grep -o "[0-9]*" || echo "0")
          failed_tests=$(grep -o "[0-9]* failures" "$log_file" 2>/dev/null | head -1 | grep -o "[0-9]*" || echo "0")
      else
          total_tests="0"
          passed_tests="0"
          failed_tests="0"
      fi
      
      echo "$total_tests"
      echo "$passed_tests"
      echo "$failed_tests"
}

# Function to parse Python test results
parse_python_test_results() {
      local log_file="$1"
      if [ ! -f "$log_file" ]; then
          echo "0"
          echo "0"
          echo "0"
          return
      fi
      
      # Try to parse Python test output
      if grep -q "Found.*test" "$log_file" 2>/dev/null; then
          total_tests=$(grep "Found.*test" "$log_file" 2>/dev/null | head -1 | grep -o "[0-9]*" || echo "0")
          passed_tests=$(grep -q "PASSED" "$log_file" 2>/dev/null && echo "$total_tests" || echo "0")
          failed_tests=$(grep -q "FAILED" "$log_file" 2>/dev/null && echo "0" || echo "0")
      else
          total_tests="0"
          passed_tests="0"
          failed_tests="0"
      fi
      
      echo "$total_tests"
      echo "$passed_tests"
      echo "$failed_tests"
}

# Function to extract failed test cases
extract_failed_tests() {
      local log_file="$1"
      
      if [ ! -f "$log_file" ]; then
          echo "No test log file available"
          return
      fi
      
      echo "Failed test cases:"
      
      # Try to extract gtest failures
      if grep -q "[  FAILED  ]" "$log_file" 2>/dev/null; then
          grep "[  FAILED  ]" "$log_file" 2>/dev/null | sed 's/\[  FAILED  \]//' | sed 's/^/  - /' || echo "  None detected"
      else
          echo "  None detected in C++ tests"
      fi
      
      # Try to extract Python test failures
      if grep -q "Error\|FAILED" "$log_file" 2>/dev/null; then
          echo "Python test issues:"
          grep -E "Error|FAILED" "$log_file" 2>/dev/null | head -10 | sed 's/^/  - /' || echo "  None detected"
      fi
}

# Function to generate JSON report
generate_json_report() {
      log_info "Generating JSON report: $REPORT_OUTPUT_FILE"
      
      # Parse C++ test results
      if [ -n "$BUILD_LOG_FILE" ] && [ -f "$BUILD_LOG_FILE" ]; then
          cpp_total=$(parse_test_results "$BUILD_LOG_FILE" | head -1)
          cpp_passed=$(parse_test_results "$BUILD_LOG_FILE" | tail -2 | head -1)
          cpp_failed=$(parse_test_results "$BUILD_LOG_FILE" | tail -1)
      else
          cpp_total="0"
          cpp_passed="0"
          cpp_failed="0"
      fi
      
      # Parse Python test results
      if [ -n "$PY_TEST_RESULTS_FILE" ] && [ -f "$PY_TEST_RESULTS_FILE" ]; then
          py_total=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | head -1)
          py_passed=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | tail -2 | head -1)
          py_failed=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | tail -1)
      else
          py_total="0"
          py_passed="0"
          py_failed="0"
      fi
      
      # Calculate overall statistics
      overall_total=$((cpp_total + py_total))
      overall_passed=$((cpp_passed + py_passed))
      overall_failed=$((cpp_failed + py_failed))
      
      # Calculate success rate
      if [ "$overall_total" -gt 0 ]; then
          success_rate=$((overall_passed * 100 / overall_total))
      else
          success_rate="N/A"
      fi
      
      # Extract failed tests
      failed_tests="None"
      
      if [ -n "$BUILD_LOG_FILE" ] && [ -f "$BUILD_LOG_FILE" ]; then
          failed_tests=$(extract_failed_tests "$BUILD_LOG_FILE")
      fi
      
      if [ -n "$PY_TEST_RESULTS_FILE" ] && [ -f "$PY_TEST_RESULTS_FILE" ]; then
          py_failed_cases=$(extract_failed_tests "$PY_TEST_RESULTS_FILE")
          if [ "$py_failed_cases" != "None detected" ]; then
              failed_tests="$failed_tests"$'\n'"$py_failed_cases"
          fi
      fi
      
      # Create JSON report
      cat > "$REPORT_OUTPUT_FILE" << EOF
{
      "build_info": {
      "build_start_time": "$BUILD_START_TIME",
      "build_end_time": "$BUILD_END_TIME",
      "build_status": "$BUILD_STATUS",
      "build_duration_secs": "dynamic",
      "build_type": "Scheduled Build"
      },
      "test_results": {
      "cpp_tests": {
          "total": $cpp_total,
          "passed": $cpp_passed,
          "failed": $cpp_failed
      },
      "python_tests": {
          "total": $py_total,
          "passed": $py_passed,
          "failed": $py_failed
      },
      "overall": {
          "total": $overall_total,
          "passed": $overall_passed,
          "failed": $overall_failed,
          "success_rate": "$success_rate"
      }
      },
      "failed_tests": {
      "details": $(echo "$failed_tests" | jq -Rs .)
      },
      "metadata": {
      "report_generated_at": "$BUILD_END_TIME",
      "report_version": "1.0",
      "ci_system": "Gradle-based CI"
      }
}
EOF

      log_success "JSON report generated successfully"
}

# Function to generate human-readable text report
generate_text_report() {
      log_info "Generating text report: $REPORT_TEXT_FILE"
      
      # Parse test results for text report
      if [ -n "$BUILD_LOG_FILE" ] && [ -f "$BUILD_LOG_FILE" ]; then
          cpp_total=$(parse_test_results "$BUILD_LOG_FILE" | head -1)
          cpp_passed=$(parse_test_results "$BUILD_LOG_FILE" | tail -2 | head -1)
          cpp_failed=$(parse_test_results "$BUILD_LOG_FILE" | tail -1)
      else
          cpp_total="0"
          cpp_passed="0"
          cpp_failed="0"
      fi
      
      if [ -n "$PY_TEST_RESULTS_FILE" ] && [ -f "$PY_TEST_RESULTS_FILE" ]; then
          py_total=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | head -1)
          py_passed=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | tail -2 | head -1)
          py_failed=$(parse_python_test_results "$PY_TEST_RESULTS_FILE" | tail -1)
      else
          py_total="0"
          py_passed="0"
          py_failed="0"
      fi
      
      overall_total=$((cpp_total + py_total))
      overall_passed=$((cpp_passed + py_passed))
      overall_failed=$((cpp_failed + py_failed))
      
      if [ "$overall_total" -gt 0 ]; then
          success_rate=$((overall_passed * 100 / overall_total))
      else
          success_rate="N/A"
      fi
      
      # Create text report
      cat > "$REPORT_TEXT_FILE" << EOF
================================================================================
                      CI BUILD RESULTS REPORT
================================================================================

Build Information:
- Build Type: Scheduled Build
- Start Time: $BUILD_START_TIME
- End Time: $BUILD_END_TIME
- Build Status: $BUILD_STATUS

Test Results Summary:
================================================================================

C++ Unit Tests:
- Total: $cpp_total
- Passed: $cpp_passed
- Failed: $cpp_failed

Python Test Cases:
- Total: $py_total
- Passed: $py_passed
- Failed: $py_failed

Overall Results:
- Total Tests: $overall_total
- Passed: $overall_passed
- Failed: $overall_failed
- Success Rate: ${success_rate}%

================================================================================

Failed Test Details:
$(extract_failed_tests "$BUILD_LOG_FILE")

================================================================================

Report generated at: $BUILD_END_TIME
Report version: 1.0
================================================================================
EOF

      log_success "Text report generated successfully"
}

# Main execution
main() {
      log_info "================================================"
      log_info "CI Build Results Recording"
      log_info "================================================"
      
      # Validate required parameters
      if [ -z "$BUILD_STATUS" ]; then
          log_error "BUILD_STATUS not specified"
          log_info "Usage: $0 --build-status=STATUS --test-status=STATUS [options]"
          exit 1
      fi
      
      if [ -z "$TEST_STATUS" ]; then
          log_error "TEST_STATUS not specified"
          exit 1
      fi
      
      log_info "Build Status: $BUILD_STATUS"
      log_info "Test Status: $TEST_STATUS"
      log_info "Start Time: $BUILD_START_TIME"
      log_info "End Time: $BUILD_END_TIME"
      
      # Generate reports
      generate_json_report
      generate_text_report
      
      log_info "================================================"
      log_success "Build results recording completed"
      log_info "JSON report: $REPORT_OUTPUT_FILE"
      log_info "Text report: $REPORT_TEXT_FILE"
      log_info "================================================"
      
      # Display summary
      echo ""
      cat "$REPORT_TEXT_FILE"
}

# Execute main function
main "$@"