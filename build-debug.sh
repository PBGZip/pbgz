#!/bin/bash
set -e

# Get script directory
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Default compiler and jobs
COMPILER="gcc"
JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --compiler=*)
      COMPILER="${1#*=}"
      shift
      ;;
    --jobs=*)
      JOBS="${1#*=}"
      shift
      ;;
    -j)
      if [[ $# -gt 1 ]]; then
        JOBS="$2"
        shift 2
      else
        echo "Error: -j requires a number"
        exit 1
      fi
      ;;
    --help)
      echo "Usage: $0 [--compiler=gcc|clang|icc] [--jobs=N|-j N]"
      echo "Build Debug version of PBGZ"
      echo ""
      echo "Options:"
      echo "  --compiler=NAME    Select compiler (gcc, clang, icc)"
      echo "  --jobs=N, -j N     Number of parallel build jobs (default: auto)"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo "Use --help for usage information"
      exit 1
      ;;
  esac
done

# Determine architecture
ARCH=$(uname -m)
if [[ "$ARCH" == "x86_64" ]]; then
    ARCH="x86_64"
elif [[ "$ARCH" == "arm64" || "$ARCH" == "aarch64" ]]; then
    ARCH="arm64"
else
    echo "Unknown architecture: $ARCH"
    exit 1
fi

# Check dependencies
echo "Checking dependencies..."
source "$ROOT_DIR/scripts/check_dependencies.sh"

# Set compiler
case $COMPILER in
  gcc)
    export CC=gcc
    export CXX=g++
    ;;
  clang)
    export CC=clang
    export CXX=clang++
    ;;
  icc)
    export CC=icc
    export CXX=icpc
    ;;
  *)
    echo "Unsupported compiler: $COMPILER"
    echo "Supported compilers: gcc, clang, icc"
    exit 1
    ;;
esac

# Check if compiler exists
if ! command -v $CC &> /dev/null; then
    echo "Error: Compiler $CC not found"
    exit 1
fi

# Create build directory
BUILD_DIR="$ROOT_DIR/build-$ARCH-$COMPILER-debug"
INSTALL_DIR="$ROOT_DIR/release-debug"

echo "Creating build directory: $BUILD_DIR"
mkdir -p "$BUILD_DIR"

echo "Configuring CMake project..."
cd "$BUILD_DIR"
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DENABLE_SANITIZER=ON \
      "$ROOT_DIR"

echo "Building project with $JOBS parallel jobs..."
cmake --build . -- VERBOSE=1 -j"$JOBS"

echo "Installing to $INSTALL_DIR..."
mkdir -p "$INSTALL_DIR"
cmake --install .

echo "Build complete. Executable located at: $INSTALL_DIR/bin/pbgz"