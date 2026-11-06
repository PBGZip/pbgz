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
source "$ROOT_DIR/prepare-build-env.sh"

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

# Provide an explicit build triplet instead of relying on autoconf's
# config.guess.  On some hosts (notably newer AMD servers) config.guess can
# return an empty result, making every bundled ./configure fail with
# "cannot guess build type; you must specify one".  The compiler reports a
# config.sub-compatible triplet on both x86_64 and aarch64; fall back to the
# common uname values only when it cannot.
build_triplet=$(bash "$ROOT_DIR/detect-build-triplet.sh")
other_configure_flags=""
if [ -n "$build_triplet" ]; then
    other_configure_flags="--build=$build_triplet"
    echo "----[info: using build triplet $build_triplet]----"
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
tar zxvf libtool-2.4.6.tar.gz
cd libtool-2.4.6
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
        ./configure $other_configure_flags --target=darwin --prefix=$INSTALL_PATH CPPFLAGS="-fPIC -fdeclspec" CFLAGS="-O3 ${PBGZ_ARM_MARCH} -fPIC -fdeclspec" 
    else
        echo "----[info: isa-l for MACOSX intel]----"
        # for MAC OSX intel version
        ./configure $other_configure_flags --target=darwin  --prefix=$INSTALL_PATH CPPFLAGS="-fPIC" CFLAGS="-O3 ${PBGZ_ARM_MARCH} -fPIC" 
    fi
else
    # ---- aarch64 non-SVE build workaround -------------------------------------
    # isa-l >= 2.32 always compiles erasure_code/aarch64/gf_nvect_dot_prod_sve.c on
    # aarch64, which requires GCC's SVE intrinsics header <arm_sve.h>. ARMv8.0
    # baseline toolchains (FT-D3000 / Kylin GCC 9.3) ship no SVE support, so that
    # translation unit cannot compile. On such a CPU the SVE kernels are never
    # selected at runtime (ISA-L dispatches on HWCAP_SVE, which is 0), and pbgz
    # only uses ISA-L's isal_inflate (igzip path - no SVE). We therefore drop the
    # SVE C file from the aarch64 source list and satisfy the EC dispatcher /
    # high-level references with stubs that delegate to the generic
    # gf_vect_dot_prod_base implementation. The official library tarball is left
    # untouched; this patch is applied to the extracted build tree.
    ARCH=`uname -m`
    if [[ "$ARCH" == "aarch64" || "$ARCH" == "arm64" ]]; then
        echo "----[info: isa-l aarch64 - applying non-SVE build workaround]----"
        sed -i -e '/erasure_code\/aarch64\/gf_nvect_dot_prod_sve\.c/d' \
               -e '/erasure_code\/aarch64\/gf_nvect_dot_prod_sve\.lo/d' \
               erasure_code/aarch64/Makefile.am
        cat >> erasure_code/aarch64/ec_aarch64_highlevel_func.c << 'EOF'

/*
 * Generated by build.intel.sh (aarch64, no SVE): ISA-L >= 2.32 builds
 * erasure_code/aarch64/gf_nvect_dot_prod_sve.c on aarch64 which requires GCC SVE
 * intrinsics (<arm_sve.h>). ARMv8.0-baseline toolchains have no SVE, so that file
 * is excluded from the build. The SVE/SVE2 n-vect dot-product entry points below
 * are referenced by the EC multibinary dispatcher and this high-level file, but
 * are never reached at runtime (ISA-L gates them behind HWCAP_SVE, which is 0 on
 * this CPU). These stubs delegate to the generic gf_vect_dot_prod_base so the
 * library always links and stays correct even if the SVE path were selected.
 */
#include <stddef.h>

void
gf_vect_dot_prod_base(int len, int vlen, unsigned char *v, unsigned char **src,
                      unsigned char *dest);

static void
gf_vect_dot_prod_base_n(int len, int vlen, unsigned char *gftbls,
                        unsigned char **src, unsigned char **dest, int n)
{
        int d;
        for (d = 0; d < n; d++)
                gf_vect_dot_prod_base(len, vlen,
                                      gftbls + (size_t) d * vlen * 32, src,
                                      dest[d]);
}

void
gf_vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                     unsigned char *dest)
{
        gf_vect_dot_prod_base(len, vlen, gftbls, src, dest);
}

void
gf_2vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 2);
}

void
gf_3vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 3);
}

void
gf_4vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 4);
}

void
gf_5vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 5);
}

void
gf_6vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 6);
}

void
gf_7vect_dot_prod_sve(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 7);
}

void
gf_vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                      unsigned char *dest)
{
        gf_vect_dot_prod_base(len, vlen, gftbls, src, dest);
}

void
gf_2vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 2);
}

void
gf_3vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 3);
}

void
gf_4vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 4);
}

void
gf_5vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 5);
}

void
gf_6vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 6);
}

void
gf_7vect_dot_prod_sve2(int len, int vlen, unsigned char *gftbls, unsigned char **src,
                       unsigned char **dest)
{
        gf_vect_dot_prod_base_n(len, vlen, gftbls, src, dest, 7);
}
EOF
        echo "----[info: isa-l aarch64 - SVE C source excluded, stubs appended]----"
    else
        echo "----[info: isa-l NOT aarch64, no SVE workaround needed]----"
    fi

    ./autogen.sh
    ./configure $other_configure_flags --prefix=$INSTALL_PATH CPPFLAGS="-fPIC" CFLAGS="-O3 ${PBGZ_ARM_MARCH} -fPIC"
fi

make clean
make
make install
make clean



rm -rf $root || true # Ignore rm errors

