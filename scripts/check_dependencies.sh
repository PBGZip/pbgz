#!/bin/bash

# Script to check third-party library dependencies

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
THIRDPARTY_DIR="$ROOT_DIR/3rd_party/release"
MISSING_LIBS=0

echo "Checking third-party library dependencies..."

# Check required library files
declare -a REQUIRED_LIBS=(
    "intel/isa-l/lib/libisal.a"
    "intel/isa-l/include/isa-l.h"
    "libdeflate/lib/libdeflate.a"
    "libdeflate/include/libdeflate.h"
    "hts/lib/libhts.a"
    "hts/include/htslib/hts.h"
    "jsoncpp/lib/libjsoncpp.a"
    "jsoncpp/include/json/json.h"
    "zstd/lib/libzstd.a"
    "zstd/include/zstd.h"
    "bzip2/lib/libbz2.a"
    "zlib/lib/libz.a"
)

for lib in "${REQUIRED_LIBS[@]}"; do
    if [ ! -e "$THIRDPARTY_DIR/$lib" ]; then
        echo "Missing library file: $THIRDPARTY_DIR/$lib"
        MISSING_LIBS=1
    fi
done

if [ $MISSING_LIBS -eq 1 ]; then
    echo "Error: Some dependency libraries are missing."
    echo "Please run '$ROOT_DIR/3rdparty/build.intel.sh' to build required third-party libraries."
    exit 1
else
    echo "All dependency libraries verified."
fi