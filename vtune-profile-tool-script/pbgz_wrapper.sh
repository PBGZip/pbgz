#!/bin/bash

# PBGZ VTune Profiling Wrapper Script
# This script wraps the pbgz command for VTune profiling
# Usage: ./pbgz_wrapper.sh <source_directory> [output_file]

set -e

# Script directory and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${PROJECT_ROOT}/vtune-profiles-result"
PBGZ_BIN="${PROJECT_ROOT}/release-release-profile/bin/pbgz"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 <source_directory> [output_file] [pbgz_options]

Arguments:
    source_directory    Path to the directory to compress (required)
    output_file         Output file name (optional, default: test.pbgz)
    pbgz_options       Additional options to pass to pbgz (optional)

Options:
    -h, --help          Show this help message and exit

Examples:
    $0 ./mydata
    $0 ./mydata my_archive.pbgz
    $0 ./mydata test.pbgz -t 4 -v
    $0 --help

Environment Variables:
    PBGZ_BIN           Path to pbgz executable (default: ../release-release-profile/bin/pbgz)

VTune Analysis:
    This script runs threading analysis (concurrency and parallelism) by default.
    Results are saved in: ../vtune-profiles-result/

EOF
}

# Function to setup VTune environment
setup_vtune_environment() {
    print_info "Setting up VTune environment..."
    
    # Check if vtune command is already available
    if command -v vtune &> /dev/null; then
        print_info "VTune command is already available"
        return 0
    fi
    
    # Try to find and source VTune environment script
    local vtune_env_scripts=(
        "/opt/intel/oneapi/vtune/latest/vtune-vars.sh"
        "/opt/intel/vtune_amplifier/vtune-vars.sh"
        "/opt/intel/vtune/vtune-vars.sh"
    )
    
    for vtune_script in "${vtune_env_scripts[@]}"; do
        if [[ -f "$vtune_script" ]]; then
            print_info "Found VTune environment script: $vtune_script"
            if source "$vtune_script"; then
                print_info "Successfully sourced VTune environment script"
                return 0
            else
                print_warn "Failed to source VTune environment script: $vtune_script"
            fi
        fi
    done
    
    print_warn "VTune environment script not found in common locations"
    print_warn "Please install Intel VTune Profiler or set up the environment manually"
    return 1
}

# Function to check dependencies
check_dependencies() {
    print_info "Checking dependencies..."
    
    # Check if pbgz exists
    if [[ ! -f "$PBGZ_BIN" ]]; then
        print_error "pbgz binary not found at: $PBGZ_BIN"
        print_error "Please build pbgz first or set PBGZ_BIN environment variable"
        exit 1
    fi
    
    # Check if pbgz is executable
    if [[ ! -x "$PBGZ_BIN" ]]; then
        print_error "pbgz binary is not executable: $PBGZ_BIN"
        exit 1
    fi
    
    # Setup VTune environment
    setup_vtune_environment
    
    # Check if vtune is available after environment setup
    if ! command -v vtune &> /dev/null; then
        print_warn "vtune command still not available after environment setup"
        print_warn "VTune profiling will be skipped"
        return 1
    else
        print_info "VTune is available for profiling"
        return 0
    fi
}

# Function to validate input
validate_input() {
    local source_dir="$1"
    
    if [[ -z "$source_dir" ]]; then
        print_error "Source directory is required"
        show_usage
        exit 1
    fi
    
    if [[ ! -d "$source_dir" ]]; then
        print_error "Source directory does not exist: $source_dir"
        exit 1
    fi
    
    if [[ ! -r "$source_dir" ]]; then
        print_error "Source directory is not readable: $source_dir"
        exit 1
    fi
    
    # Get absolute path
    SOURCE_DIR="$(realpath "$source_dir")"
    print_info "Source directory: $SOURCE_DIR"
}

# Function to create results directory
create_results_dir() {
    if [[ ! -d "$RESULTS_DIR" ]]; then
        mkdir -p "$RESULTS_DIR"
        print_info "Created results directory: $RESULTS_DIR"
    fi
}

# Function to run pbgz command
run_pbgz() {
    local output_file="$1"
    shift
    local pbgz_opts="$@"
    
    print_info "Starting pbgz compression..."
    
    # Use pipe compression method
    print_info "Command: tar -cf - \"$SOURCE_DIR\" | \"$PBGZ_BIN\" $pbgz_opts -o \"$output_file\""
    
    # Create a temporary script for VTune to profile
    local temp_script="${RESULTS_DIR}/pbgz_temp_script.sh"
    cat > "$temp_script" << EOF
#!/bin/bash
# Temporary script for VTune profiling
set -e
cd "$SOURCE_DIR"
tar -cf - . | "$PBGZ_BIN" $pbgz_opts -o "$PWD/$output_file"
EOF
    
    chmod +x "$temp_script"
    
    # Run the command
    cd "$RESULTS_DIR"
    if tar -cf - "$SOURCE_DIR" | "$PBGZ_BIN" $pbgz_opts -o "$output_file"; then
        print_info "pbgz compression completed successfully"
        print_info "Output file: $RESULTS_DIR/$output_file"
        
        # Show file size info
        local original_size=$(du -sh "$SOURCE_DIR" | cut -f1)
        local compressed_size=$(du -sh "$output_file" | cut -f1)
        print_info "Original size: $original_size"
        print_info "Compressed size: $compressed_size"
    else
        print_error "pbgz compression failed"
        exit 1
    fi
    
    # Clean up temporary script
    rm -f "$temp_script"
}

# Function to run VTune profiling
run_vtune_profiling() {
    local output_file="$1"
    shift
    local pbgz_opts="$@"
    
    if ! command -v vtune &> /dev/null; then
        print_warn "VTune not available, skipping profiling"
        return 0
    fi
    
    print_info "Starting VTune profiling..."
    
    # Create profiling results directory
    local profiling_dir="${RESULTS_DIR}/vtune_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$profiling_dir"
    
    # Run different types of analysis
    local analyses=( "threading" )
    
    for analysis in "${analyses[@]}"; do
        print_info "Running $analysis analysis..."
        
        local analysis_result_dir="${profiling_dir}/${analysis}"
        mkdir -p "$analysis_result_dir"
        
        # Create a script for this analysis
        local analysis_script="${analysis_result_dir}/run_analysis.sh"
        cat > "$analysis_script" << EOF
#!/bin/bash
cd "$SOURCE_DIR"
temp_tar="\${PWD}/temp_\${output_file%.pbgz}.tar"
tar -cf "\$temp_tar" -C "$SOURCE_DIR" .
"$PBGZ_BIN" $pbgz_opts -o "\$PWD/$output_file" "\$temp_tar"
rm -f "\$temp_tar"
EOF
        
        chmod +x "$analysis_script"
        
        # Run VTune analysis
        if vtune -collect "$analysis" \
                -result-dir "$analysis_result_dir" \
                -- "$analysis_script"; then
            print_info "$analysis analysis completed successfully"
            
            # Generate HTML report
            local report_file="${analysis_result_dir}/${analysis}_report.html"
            vtune -report "$analysis" \
                    -result-dir "$analysis_result_dir" \
                    -format html \
                    -output "$report_file" || true
            
            print_info "Report generated: $report_file"
        else
            print_warn "$analysis analysis failed"
        fi
        
        # Clean up analysis script
        rm -f "$analysis_script"
    done
    
    print_info "VTune profiling completed"
    print_info "Results directory: $profiling_dir"
    print_info "To view results in GUI: vtune-gui $profiling_dir/hotspots"
}

# Function to show summary
show_summary() {
    local output_file="$1"
    
    print_info "=== Summary ==="
    print_info "Source directory: $SOURCE_DIR"
    print_info "Output file: $RESULTS_DIR/$output_file"
    
    if [[ -f "$RESULTS_DIR/$output_file" ]]; then
        local file_size=$(du -sh "$RESULTS_DIR/$output_file" | cut -f1)
        print_info "File size: $file_size"
    fi
    
    if [[ -d "$RESULTS_DIR" ]] && [[ "$(ls -A "$RESULTS_DIR" 2>/dev/null)" ]]; then
        print_info "Results directory: $RESULTS_DIR"
        print_info "Contents:"
        ls -la "$RESULTS_DIR"
    fi
    
    print_info "=== End Summary ==="
}

# Function to list available analyses
list_analyses() {
    cat << EOF
Available VTune Analyses:
1. hotspots          - Hotspots analysis (CPU bottlenecks)
   Description: Identifies CPU hotspots and functions consuming the most CPU time
   
2. threading         - Threading analysis (concurrency and parallelism)
   Description: Analyzes threading performance, concurrency issues, and parallelism
   
3. performance-snapshot - Performance snapshot overview
   Description: Provides a quick overview of system performance metrics
   
4. memory-access     - Memory access analysis
   Description: Analyzes memory access patterns, cache misses, and memory bottlenecks
   
5. hpc-performance   - HPC performance analysis
   Description: High-performance computing specific performance analysis

Usage:
   By default, only 'threading' analysis is run (most commonly used for concurrency)
   To run all analyses, modify the 'analyses' array in the run_vtune_profiling function
   
EOF
}

# Main function
main() {
    # Handle help flag first
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        show_usage
        exit 0
    fi
    
    local source_dir="$1"
    local output_file="${2:-test.pbgz}"
    shift 2
    local pbgz_opts="$@"
    
    print_info "=== PBGZ VTune Profiling Wrapper ==="
    
    # Validate input
    validate_input "$source_dir"
    
    # Check dependencies
    check_dependencies
    
    # Create results directory
    create_results_dir
    
    # Run pbgz
    run_pbgz "$output_file" $pbgz_opts
    
    # Run VTune profiling
    run_vtune_profiling "$output_file" $pbgz_opts
    
    # Show summary
    show_summary "$output_file"
    
    print_info "=== Process Completed ==="
}

# Run main function with all arguments
main "$@"
