#!/bin/bash

set -e # Exit script immediately on error

buildtype=release

# build type
if [[ $1 ]] && [[ $1 == "debug" ]]; then
    buildtype=debug
    echo "----[error type: build type $buildtype not defined]----"
    exit 1 # Correct behavior - debug build type not supported, display error and exit
fi

export BUILDTYPE=$buildtype
echo "----[build type: $buildtype]----"

# build path
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

buildpath=$ROOT_DIR/../3rd_party/$buildtype

if [ ! -d "$buildpath" ]; then
  mkdir -p $buildpath
fi
echo "----[build path: $buildpath]----"

cd $ROOT_DIR

libname=intel
rm -rf $libname || true # Add || true to ignore rm errors
tar -jxf $libname.tar.bz2
root=`pwd`/$libname
path=`pwd`
cd $libname

# start to build

rm -frf $buildpath/tools || true

export PATH=$buildpath/tools/bin:$PATH
export LD_LIBRARY_PATH=$buildpath/tools/lib:$LD_LIBRARY_PATH

# 检查当前系统是否为aarch64

if [ "$(uname -m)" == "aarch64" ]; then
    other_configure_flags="--build=aarch64-unknown-linux"
fi

# install m4
cd $path/
echo "----[info: installing m4 1.4.9]----"
#wget -O m4-1.4.9.tar.gz http://ftp.gnu.org/gnu/m4/m4-1.4.9.tar.gz
rm -rf m4-1.4.9 || true
tar -zxf m4-1.4.9.tar.gz
cd m4-1.4.9
./configure --prefix=$buildpath/tools  $other_configure_flags
make && $supper_execute make install 
rm -rf $path/m4-1.4.9 || true;


# install autoconf 2.69
cd $path/
echo "----[info: installing autoconf 2.69]----"
#wget ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
rm -rf autoconf-2.69 || true
tar zxvf autoconf-2.69.tar.gz
cd autoconf-2.69
./configure --prefix=$buildpath/tools $other_configure_flags
make && $supper_execute make install && cd .. && { rm -rf $path/autoconf-2.69 || true; }
  

# install automake 1.17
cd $path/
echo "----[info: installing automake 1.17]----"
# wget http://mirrors.kernel.org/gnu/automake/automake-1.17.tar.gz 
rm -rf automake-1.17 || true
tar zxvf automake-1.17.tar.gz
cd automake-1.17
./configure --prefix=$buildpath/tools  $other_configure_flags
make && $supper_execute make install && cd .. && { rm -rf $path/automake-1.17 || true; }


# install yasm 1.3.0
cd $path/
echo "----[info: installing yasm 1.3.0]----"
#wget http://www.tortall.net/projects/yasm/releases/yasm-1.3.0.tar.gz
rm -rf yasm-1.3.0 || true
tar zxvf yasm-1.3.0.tar.gz
cd yasm-1.3.0
./configure $other_configure_flags --prefix=$buildpath/tools && make && $supper_execute make install && cd .. && { rm -rf $path/yasm-1.3.0 || true; }
  
# install nasm 2.14.02
cd $path/
echo "----[info: installing nasm 2.14.02]----"
#wget https://www.nasm.us/pub/nasm/releasebuilds/2.14.02/nasm-2.14.02.tar.bz2
rm -rf nasm-2.14.02 || true
tar -xf nasm-2.14.02.tar.bz2
cd nasm-2.14.02
./configure $other_configure_flags --prefix=$buildpath/tools && make && $supper_execute make install && cd .. && { rm -rf $path/nasm-2.14.02 || true; }
  

# install libtool 2.4.6
cd $path/
echo "----[info: installing libtool 2.4.6]----"
rm -rf libtool-2.4.6 || true
#wget http://mirrors.ustc.edu.cn/gnu/libtool/libtool-2.4.6.tar.gz
tar zxvf libtool-2.4.6.tar.gz && cd libtool-2.4.6
./configure $other_configure_flags --prefix=$buildpath/tools && make && $supper_execute make install && cd .. && { rm -rf $path/libtool-2.4.6 || true; }


cd $path/intel/isa-l/
echo "----[info: installing isa-l]----"
INSTALL_PATH=$buildpath/intel/isa-l/
mkdir -p $INSTALL_PATH 2> /dev/null

# if MACOSX
if [[ $OSTYPE == darwin* ]]; then
    echo "----[info: isa-l for MACOSX]----"
    mkdir -p m4 2> /dev/null

    # Completely fix configure.ac file - remove unnecessary indentation
    cat > configure.ac.new << 'EOF'
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([libisal],
        [2.26.0],
        [sg.support.isal@intel.com],
        [isa-l],
        [http://01.org/storage-acceleration-library])
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([
    foreign
    1.11
    -Wall
    -Wno-portability
    silent-rules
    tar-pax
    no-dist-gzip
    dist-xz
    subdir-objects
])
EOF

    # Add the rest of the original file (starting from AM_PROG_AS line)
    sed -n '/^AM_PROG_AS/,$p' configure.ac >> configure.ac.new

    # Replace original file
    mv configure.ac.new configure.ac

    # Modify Makefile.am to add ACLOCAL_AMFLAGS
    if ! grep -q "ACLOCAL_AMFLAGS" Makefile.am; then
        sed -i.bak '1s/^/ACLOCAL_AMFLAGS = -I m4\n\n/' Makefile.am
    fi

    # Run autotools toolchain
    $buildpath/tools/bin/libtoolize --force --copy
    $buildpath/tools/bin/aclocal -I m4
    $buildpath/tools/bin/autoheader
    $buildpath/tools/bin/automake --add-missing --copy
    $buildpath/tools/bin/autoconf

    # Choose appropriate configuration based on processor architecture
    ARCH=`uname -m`
    if [[ $ARCH == "arm64" ]]; then
        echo "----[info: isa-l for MACOSX arm64]----"
        # for MAC OSX arm64 version
        ./configure --target=darwin --prefix=$INSTALL_PATH CPPFLAGS="-fPIC -fdeclspec" CFLAGS="-fPIC -fdeclspec" 
    else
        echo "----[info: isa-l for MACOSX intel]----"
        # for MAC OSX intel version
        ./configure --target=darwin  --prefix=$INSTALL_PATH CPPFLAGS="-fPIC" CFLAGS="-fPIC" 
    fi
else
    ./autogen.sh
    ./configure $other_configure_flags --prefix=$INSTALL_PATH CPPFLAGS="-fPIC" CFLAGS="-fPIC"
fi

make clean
make
make install
make clean



rm -rf $root || true # Ignore rm errors

