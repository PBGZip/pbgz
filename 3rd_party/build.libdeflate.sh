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

libname=libdeflate-1.7
root=`pwd`/$libname
rm -rf $libname || true # Ignore errors from rm command
tar -xvf $libname.tar.gz
cd $libname

# start to build
export INSTALL_PATH=$buildpath/libdeflate
 if [ $buildtype == "release" ] ; then
    make -j $(nproc) # Using nproc instead of cat /proc/cpuinfo for better portability
    make install 
else
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Added explicit error code
fi
unset INSTALL_PATH

rm -rf $root || true # Ignore errors from rm command

