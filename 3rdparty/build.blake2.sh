#!/bin/bash

set -e # Exit immediately if any command fails

buildtype=release

# build type
if [[ $1 ]] && [[ $1 == "debug" ]] ; then
    buildtype=debug
    echo "----[error type: build type $buildtype not defined]----"
    exit 1
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# build path
root=`pwd`
buildpath=$root/../3rdlib/$buildtype
if [ ! -d "$buildpath" ]; then
  mkdir -p $buildpath
fi
echo "----[build path: $buildpath]----"

installpath=$buildpath/blake2
rm -rf $installpath || true
mkdir -p $installpath/include
mkdir -p $installpath/lib

libname=BLAKE2-20190724
rm -rf $libname || true
tar -xvf $libname.tar.gz
cd $libname

# Detect architecture
ARCH=$(uname -m)
echo "----[detected architecture: $ARCH]----"

# start to build
if [ $buildtype == "release" ] ; then
  # Use SSE implementation for x86/x86_64 and reference implementation for ARM
  if [[ "$ARCH" == "x86_64" || "$ARCH" == "i686" ]]; then
    echo "----[using SSE implementation]----"
    cd sse
    gcc -O3 -Wall -Wextra -std=c89 -pedantic -Wno-long-long -c blake2b.c
    ar -crv libblake2b.a blake2b.o
  else
    echo "----[using reference implementation]----"
    cd ref
    gcc -O3 -Wall -Wextra -std=c89 -pedantic -Wno-long-long -c blake2b-ref.c
    ar -crv libblake2b.a blake2b-ref.o
  fi
  
  cp -f libblake2b.a $installpath/lib
  cp -f blake2.h $installpath/include
else
  echo "----[error type: build type $buildtype not defined]----"
  exit 1
fi

cd ../ && rm -rf $root/$libname || true