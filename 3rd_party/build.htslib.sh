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

libname=htslib-1.12
rm -rf $libname || true # Ignore errors from rm command
tar -xvf $libname.tar.bz2
cd $libname

export PATH=$buildpath/tools/bin:$PATH
export LD_LIBRARY_PATH=$buildpath/tools/lib:$LD_LIBRARY_PATH

export LIBRARY_PATH=$buildpath/libdeflate/lib:$buildpath/bzip2/lib
export C_INCLUDE_PATH=$buildpath/libdeflate/include:$buildpath/bzip2/include

# start to build
 if [ $buildtype == "release" ] ; then
    autoreconf -i
    ./configure CFLAGS=-fPIC --prefix=$buildpath/hts --disable-lzma --with-libdeflate --enable-bz2 --disable-libcurl --disable-s3
    make -j $(nproc) # Using nproc instead of cat /proc/cpuinfo for better portability
    make install 
else
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Added explicit error code
fi


cd ../ && rm -rf $libname || true # Ignore errors from rm command

