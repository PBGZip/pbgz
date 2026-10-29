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

libname=zlib-1.3.1
root=`pwd`/$libname
rm -rf $libname || true # Ignore errors from rm command
tar -xvf $libname.tar.gz
cd $libname

# start to build
 if [ $buildtype == "release" ] ; then
  # -fPIC alone would *override* zlib's default CFLAGS, silently dropping its
  # -O3 and leaving deflate/crc at unoptimized speed (measured: 17.7 MB/s vs
  # 33.6 MB/s on the same 236MB BAM stream). Keep -O3.
  # On Linux ARMv8, enable the CRC32 extension so zlib uses hardware CRC.
  # x86 builds keep exactly the previous flags.
  ZLIB_CFLAGS="-O3 -fPIC"
  if [ "$(uname -s)" = "Linux" ] && { [ "$(uname -m)" = "aarch64" ] || [ "$(uname -m)" = "arm64" ]; }; then
    ZLIB_CFLAGS="$ZLIB_CFLAGS -march=armv8-a+crc"
  fi
  CFLAGS="$ZLIB_CFLAGS" ./configure --prefix=$buildpath/zlib
  make -j $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2) # Using nproc instead of cat /proc/cpuinfo for better portability
  make install 
else
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Added explicit error code
fi

rm -rf $root || true # Ignore errors from rm command
