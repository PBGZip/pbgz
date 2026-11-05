cmake_minimum_required(VERSION 3.10)
project(pbgz VERSION 2.0.0 LANGUAGES C CXX)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include build options
include(BuildOptions)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Define third-party library paths first (make these available globally)
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rd_party)

# Include third-party library headers globally
include_directories(
    ${THIRDPARTY_DIR}/release/jsoncpp/include
    ${THIRDPARTY_DIR}/release/zstd/include
    ${THIRDPARTY_DIR}/release/hts/include
    ${THIRDPARTY_DIR}/release/intel/isa-l/include
    ${THIRDPARTY_DIR}/release/libdeflate/include
    ${THIRDPARTY_DIR}/release/bzip2/include
    ${THIRDPARTY_DIR}/release/zlib/include
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/release/zstd/lib/libzstd.a
    ${THIRDPARTY_DIR}/release/hts/lib/libhts.a
    ${THIRDPARTY_DIR}/release/intel/isa-l/lib/libisal.a
    ${THIRDPARTY_DIR}/release/libdeflate/lib/libdeflate.a
    ${THIRDPARTY_DIR}/release/bzip2/lib/libbz2.a
    ${THIRDPARTY_DIR}/release/zlib/lib/libz.a
)

add_subdirectory(coder)

# Find all source files
file(GLOB PBGZ_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/log/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/farmhash/src/farmhash.cc
)

set(MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )

# Define main executable
add_executable(pbgz ${MAIN_SOURCES} ${PBGZ_SOURCES} ${PBGZ_OPEN_SOURCES})

# ARMv8.2 DotProd base packing: only coder_factory.cpp (the TU that uses
# packBases) is compiled with the DotProd extension, so the inline UDOT path in
# base_pack.h is used there. The rest of the program keeps the detected -march
# (incl. +nolse), so LSE atomics stay disabled everywhere else.
if(PBGZ_ARCH STREQUAL "arm64" AND PBGZ_CPU_HAS_DOTPROD)
    if(PBGZ_CPU_HAS_FP16)
        set(_pbgz_dot_march "armv8.2-a+dotprod+fp16")
    else()
        set(_pbgz_dot_march "armv8.2-a+dotprod")
    endif()
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/coder_factory.cpp PROPERTIES
        COMPILE_FLAGS "-march=${_pbgz_dot_march}"
        COMPILE_DEFINITIONS "PBGZ_BASEPACK_INLINE_DOT=1")
endif()

target_link_libraries(pbgz PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/release/zstd/lib/libzstd.a
    ${THIRDPARTY_DIR}/release/hts/lib/libhts.a
    ${THIRDPARTY_DIR}/release/intel/isa-l/lib/libisal.a
    ${THIRDPARTY_DIR}/release/libdeflate/lib/libdeflate.a
    ${THIRDPARTY_DIR}/release/bzip2/lib/libbz2.a
    ${THIRDPARTY_DIR}/release/zlib/lib/libz.a
    pthread
    curl
)

if(APPLE) 
    target_link_libraries(pbgz_test PRIVATE  "-framework CoreServices")
endif()


install(TARGETS pbgz DESTINATION bin)

