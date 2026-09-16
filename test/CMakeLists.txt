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

# The model snapshot test does not depend on an external dataset; it is kept as a small
# standalone command-line program so that prior binary compatibility, error fallback, and codec
# state synchronization can be verified before integration into the main compression path.
add_executable(fcv2_model_roundtrip fcv2_model_roundtrip.cpp)

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
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/googletest/googlemock/include)

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
    curl
    lzma
)

if(APPLE) 
    target_link_libraries(pbgz_test PRIVATE  "-framework CoreServices")
endif()

target_link_libraries(fcv2_model_roundtrip PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# Per-record fixed-length round-trip verification for coder_bwt_cm: used to pinpoint the root
# cause of decompression failures when QUAL is encoded via bwt_cm.
add_executable(bwtcm_record_roundtrip bwtcm_record_roundtrip.cpp)
target_link_libraries(bwtcm_record_roundtrip PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# How many mixing models fcv2 should use for QUAL. The mixer's cost is one
# counter update per model per coded bit, so the count sets its speed; this tool
# measures size and throughput per count on real QUAL data (see its file header).
add_executable(qual_model_scan qual_model_scan.cpp)
target_link_libraries(qual_model_scan PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# Which fcv2 context tier (cycleMax / cycleBucket / deltaMax / deltaBucket /
# prevShift) is best for a given QUAL volume. The axis is selected on the command
# line; see the file header for why cycleMax is the axis that had never been
# swept, and for the bwt_cm reference row.
add_executable(qual_tier_scan qual_tier_scan.cpp)
target_link_libraries(qual_tier_scan PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# Which fcv2 context parameters and mixing-model count are best for a given read-length range.
# The product currently decides both by rule (one threshold at 500 bytes: 4 models above it, 6
# below; five hand-picked presets trialled as they are), and this tool measures the real encoded
# size per length bucket and descends to the configuration that wins inside each one. See its
# file header for the search and for how to paste the result back as the per-length rule.
add_executable(fcv2_len_train fcv2_len_train.cpp)
target_link_libraries(fcv2_len_train PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)


add_subdirectory(googletest)

# Optional: Add gtest support
target_link_libraries(pbgz_test PRIVATE gmock gmock_main gtest)

# Enable testing
enable_testing()
add_test(NAME pbgz_test COMMAND pbgz_test)
add_test(NAME fcv2_model_roundtrip COMMAND fcv2_model_roundtrip)
