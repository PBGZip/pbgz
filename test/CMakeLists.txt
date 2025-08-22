# GoogleTest
project(pbgz_test)

# 设置测试源码目录
set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/testcase)

# Define third-party library paths first (make these available globally)
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rd_party)

# 查找所有测试源文件
file(GLOB TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testcase/*.cpp)

# 查找所有主项目源文件
file(GLOB PBGZ_SOURCES ${CMAKE_SOURCE_DIR}/code/*.cpp
                       ${CMAKE_SOURCE_DIR}/code/log/*.cpp
)

list(REMOVE_ITEM PBGZ_SOURCES ${CMAKE_SOURCE_DIR}/code/main.cpp)

set(PBGZ_OPEN_SOURCES ${THIRDPARTY_DIR}/farmhash/src/farmhash.cc)

# 添加主测试文件
set(TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/pbgz_test_main.cpp)

# 生成测试可执行文件
add_executable(pbgz_test ${TEST_MAIN} ${TEST_SOURCES} ${PBGZ_SOURCES} ${PBGZ_OPEN_SOURCES})

# 包含主项目头文件
include_directories(${CMAKE_SOURCE_DIR}/code)
include_directories(${CMAKE_SOURCE_DIR}/code/log)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/jsoncpp/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/googletest/googletest/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/farmhash/src)


set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# Link to libraries - Add PRIVATE keyword to match style with CoreServices linkage
target_link_libraries(pbgz_test PRIVATE
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
)

# 可选：添加gtest支持
target_link_libraries(pbgz_test PRIVATE gtest_main)

# 启用测试
enable_testing()
add_test(NAME pbgz_test COMMAND pbgz_test)


