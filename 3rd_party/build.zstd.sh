#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status.

buildtype=release

# build type
if [[ $1 ]] && [[ $1 == "debug" ]] ; then
    buildtype=debug
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Use non-zero exit code for errors
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# install path
installpath=`pwd`/../3rd_party/$buildtype/zstd
rm -rf $installpath || true # Ignore errors from rm
mkdir -p $installpath/lib
mkdir -p $installpath/include

echo "----[install path: $installpath]----"

root=`pwd`
libname=zstd-1.5.0
rm -rf $root/$libname || true # Ignore errors from rm
tar -xvf $libname.tar.gz
cd $libname

# start to build
if [ $buildtype == "release" ] ; then
  # Use nproc for better portability and correctness in getting CPU count
  make -j `nproc`
  # Use cp -f instead of yes|cp
  cp -f $root/$libname/lib/libzstd.a $installpath/lib/
  cp -f $root/$libname/lib/zstd.h $installpath/include/
else
  echo "----[error type: build type $buildtype not defined]----"
  exit 1
fi

cd $root && rm -rf $root/$libname || true # Ignore errors from rm

