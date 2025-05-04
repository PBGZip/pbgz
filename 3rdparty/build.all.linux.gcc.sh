#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status.


ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
pushd "$ROOT_DIR" || exit 1

# All third-party libraries are considered stable and robust code, 
# we do not need to debug these codes, and the performance of the 
# code in release mode is better

buildtype=release 

rm -rf ../3rdlib/$buildtype

# Build intel. This script must be placed first because it includes 
# tools like autoconf, which may be required by the libraries built later.
./build.intel.sh $buildtype

./build.zstd.sh $buildtype

./build.bzip2.sh $buildtype

./build.libdeflate.sh $buildtype

./build.htslib.sh $buildtype

./build.jsoncpp.sh $buildtype

./build.zlib.sh $buildtype



popd || exit 1