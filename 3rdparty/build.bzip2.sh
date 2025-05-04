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
buildpath=`pwd`/../3rdlib/$buildtype
if [ ! -d "$buildpath" ]; then
  mkdir -p $buildpath
fi
echo "----[build path: $buildpath]----"

libname=bzip2
root=`pwd`/$libname
rm -rf $libname || true # Add || true to ignore rm errors
tar -xvf $libname.tar.gz
cd $libname

# start to build
 if [ $buildtype == "release" ] ; then
  make -j $(nproc) # Using nproc instead of cat /proc/cpuinfo | grep "physical id" | wc -l for better CPU count detection
  make install PREFIX=$buildpath/bzip2
else
  echo "----[error type: build type $buildtype not defined]----"
  exit 1 # Modified to use explicit error exit code
fi

rm -rf $root || true # Add || true to ignore rm errors

