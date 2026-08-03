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

# 模型快照测试不依赖外部数据集，单独保留为小型命令行程序，便于在接入主压缩路径之前
# 直接验证先验的二进制兼容性、错误回退和编解码状态同步。
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
    "-framework CoreServices"
)

target_link_libraries(fcv2_model_roundtrip PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# 先验快照落盘工具：用真实 SAM 训练并导出模型快照，供评估通用压缩器对快照的压缩效果。
# 仅用于评估，不加入 add_test，也不参与产品功能。
add_executable(fcv2_prior_dump fcv2_prior_dump.cpp)
target_link_libraries(fcv2_prior_dump PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# 先验收益实测工具：对照"无先验整段/无先验分片/有先验分片"三种模式，量出先验在本树
# fcv2 上的实际收益。仅用于评估，不加入 add_test，也不参与产品功能。
add_executable(fcv2_prior_gain fcv2_prior_gain.cpp)
target_link_libraries(fcv2_prior_gain PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# QUAL 列编码器横向实测：补上 coder_bwt_cm 这个当前不在 QualSelector 候选集里的对照。
# 仅用于评估，不加入 add_test，也不参与产品功能。
add_executable(qual_coder_bench qual_coder_bench.cpp)
target_link_libraries(qual_coder_bench PRIVATE
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
