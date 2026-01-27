#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status

buildtype=release

# build type
 if [[ $1 ]] && [[ $1 == "debug" ]] ; then
    buildtype=debug
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# install path
installpath=`pwd`/../3rd_party/$buildtype
if [ ! -d "$installpath" ]; then
  mkdir -p $installpath
fi
echo "----[install path: $installpath]----"

root=`pwd`
libname=jsoncpp-1.9.6
rm -rf $libname || true # Ignore errors from rm command
tar -xvf $libname.tar.gz
cd $libname

# start to build
mkdir -p build/$buildtype
cd build/$buildtype
cmake -DCMAKE_BUILD_TYPE=$buildtype -DBUILD_STATIC_LIBS=ON -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=$installpath/jsoncpp -G "Unix Makefiles" ../..
make -j $(nproc) # Using nproc instead of cat /proc/cpuinfo for better portability
make install

cd $root && rm -rf $root/$libname || true # Ignore errors from rm command

