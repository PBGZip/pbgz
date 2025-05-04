cmake_minimum_required(VERSION 3.10)
project(pbgz VERSION 1.0.0 LANGUAGES C CXX)

# Set C++ standard (add these lines right after project declaration)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")


# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Include build options
include(BuildOptions)

# Define third-party library paths
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rdlib/release)

# Check for required third-party libraries
set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/intel/isa-l/lib/libisal.a
    ${THIRDPARTY_DIR}/libdeflate/lib/libdeflate.a
    ${THIRDPARTY_DIR}/hts/lib/libhts.a
    ${THIRDPARTY_DIR}/blake2/lib/libblake2b.a
    ${THIRDPARTY_DIR}/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/zstd/lib/libzstd.a
    ${THIRDPARTY_DIR}/bzip2/lib/libbz2.a
    ${THIRDPARTY_DIR}/zlib/lib/libz.a
)

foreach(lib ${REQUIRED_LIBS})
    if(NOT EXISTS ${lib})
        message(FATAL_ERROR "Missing required library: ${lib}\nPlease run scripts/build_3rdlibs.sh to build dependencies")
    endif()
endforeach()

# Include subdirectories
add_subdirectory(coder)

# Main project source files
set(MAIN_SOURCES
    main.cpp
    engine.cpp
    analy_type.cpp
    io.cpp
    file_op.cpp
    manager.cpp
    wait_event.cpp
    pbgz_file.cpp
    actuator_fastq.cpp
    actuator_everything.cpp
    reference.cpp
    md5sum.cpp
    city.cpp
)

# Add include directories
include_directories(
    ${THIRDPARTY_DIR}/intel/isa-l/include
    ${THIRDPARTY_DIR}/libdeflate/include
    ${THIRDPARTY_DIR}/hts/include
    ${THIRDPARTY_DIR}/blake2/include
    ${THIRDPARTY_DIR}/jsoncpp/include
    ${THIRDPARTY_DIR}/zstd/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/coder
)

# Define main executable
add_executable(pbgz ${MAIN_SOURCES})

# Link to libraries - Add PRIVATE keyword to match style with CoreServices linkage
target_link_libraries(pbgz PRIVATE
    pbgz_coder
    pthread
    ${THIRDPARTY_DIR}/hts/lib/libhts.a
    ${THIRDPARTY_DIR}/intel/isa-l/lib/libisal.a
    ${THIRDPARTY_DIR}/libdeflate/lib/libdeflate.a
    ${THIRDPARTY_DIR}/bzip2/lib/libbz2.a
    ${THIRDPARTY_DIR}/zlib/lib/libz.a
    ${THIRDPARTY_DIR}/blake2/lib/libblake2b.a
    ${THIRDPARTY_DIR}/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/zstd/lib/libzstd.a
)

# Check if we're on macOS
if(APPLE)
  # Find and link CoreServices framework
  find_library(CORE_SERVICES_FRAMEWORK CoreServices REQUIRED)
  target_link_libraries(pbgz PRIVATE ${CORE_SERVICES_FRAMEWORK})
  
  # You can also add a define to indicate we're using macOS APIs
  target_compile_definitions(pbgz PRIVATE USING_MAC_CORE_SERVICES=1)
endif()

# Installation rules
install(TARGETS pbgz DESTINATION bin)