# GoogleTest
project(pbgz_test)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_BUILD_TYPE "Debug")

# Include build options
include(BuildOptions)

# Set test source directory
set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/testcase)

# Define third-party library paths first (make these available globally)
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rd_party)

# Find all test source files
file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testcase/*.cpp)

# Find all main project source files
file(GLOB PBGZ_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp
                       ${CMAKE_SOURCE_DIR}/src/log/*.cpp
                       ${CMAKE_SOURCE_DIR}/src/utils/*.cpp
                       ${CMAKE_SOURCE_DIR}/src/farmhash/src/farmhash.cc
)

list(REMOVE_ITEM PBGZ_SOURCES ${CMAKE_SOURCE_DIR}/src/main.cpp)

# Add main test file
set(TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/pbgz_test_main.cpp)

# Generate test executable
add_executable(pbgz_test ${TEST_MAIN} ${TEST_SOURCES} ${PBGZ_SOURCES} ${PBGZ_OPEN_SOURCES})

# Include main project header files
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/log)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/jsoncpp/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/zstd/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/hts/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/intel/isa-l/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/libdeflate/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/bzip2/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/zlib/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/googletest/googletest/include)

# Link to libraries - Add PRIVATE keyword to match style with CoreServices linkage
target_link_libraries(pbgz_test PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/release/zstd/lib/libzstd.a
    ${THIRDPARTY_DIR}/release/hts/lib/libhts.a
    ${THIRDPARTY_DIR}/release/intel/isa-l/lib/libisal.a
    ${THIRDPARTY_DIR}/release/libdeflate/lib/libdeflate.a
    ${THIRDPARTY_DIR}/release/bzip2/lib/libbz2.a
    ${THIRDPARTY_DIR}/release/zlib/lib/libz.a
    pthread
)

add_subdirectory(googletest)

# Optional: Add gtest support
target_link_libraries(pbgz_test PRIVATE gtest_main)

# Enable testing
enable_testing()
add_test(NAME pbgz_test COMMAND pbgz_test)
