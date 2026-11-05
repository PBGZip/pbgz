#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

buildtype=release

# build type
 if [[ $1 ]] && [[ $1 == "debug" ]] ; then
    buildtype=debug
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Added explicit error code
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# build path
buildpath=`pwd`/../3rd_party/$buildtype
if [ ! -d "$buildpath" ]; then
  mkdir -p $buildpath
fi
echo "----[build path: $buildpath]----"

libname=libdeflate-1.26
root=`pwd`/$libname
rm -rf $libname || true # Ignore errors from rm command
tar -xvf $libname.tar.gz
cd $libname

# start to build
# libdeflate >= 1.19 only ships a CMake build (the top-level Makefile was dropped),
# so the old `make && make install` flow no longer applies. Configure CMake to build
# only the static lib and install into the same layout the rest of the project expects
# (release/libdeflate/{include/libdeflate.h,lib/libdeflate.a}).
 if [ $buildtype == "release" ] ; then
  mkdir -p build
  cd build
  cmake -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=$buildpath/libdeflate \
        -DCMAKE_C_FLAGS=-fPIC \
        -DLIBDEFLATE_BUILD_SHARED_LIB=OFF \
        -DLIBDEFLATE_BUILD_TESTS=OFF \
        -DLIBDEFLATE_BUILD_GZIP=OFF \
        ..
  make -j $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) # Using nproc instead of cat /proc/cpuinfo for better portability
  make install
  cd ..
else
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Added explicit error code
fi

rm -rf $root || true # Ignore errors from rm command
