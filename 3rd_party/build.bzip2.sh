#!/bin/bash

set -e # Add this line to make the script exit immediately if any command fails

buildtype=release

# build type
 if [[ $1 ]] && [[ $1 == "debug" ]] ; then
    buildtype=debug
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Modified to use explicit error exit code
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# build path
buildpath=`pwd`/../3rd_party/$buildtype
if [ ! -d "$buildpath" ]; then
  mkdir -p $buildpath
fi
echo "----[build path: $buildpath]----"

libname=bzip2-1.0.8
root=`pwd`/$libname
rm -rf $libname || true # Add || true to ignore rm errors
tar -xvf $libname.tar.gz
cd $libname

# start to build
 if [ $buildtype == "release" ] ; then
  # bzip2 is linked into htslib's libhts.so, so the static lib must be PIC
  # (mirrors the -O2 baseline flags in bzip2's own Makefile + -fPIC).
  make -j $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) CFLAGS="-Wall -Winline -O2 -g -D_FILE_OFFSET_BITS=64 -fPIC" # Using nproc instead of cat /proc/cpuinfo | grep "physical id" | wc -l for better CPU count detection
  make install PREFIX=$buildpath/bzip2
else
  echo "----[error type: build type $buildtype not defined]----"
  exit 1 # Modified to use explicit error exit code
fi

rm -rf $root || true # Add || true to ignore rm errors

