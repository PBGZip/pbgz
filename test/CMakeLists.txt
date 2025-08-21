# GoogleTest
project(pbgz_test)

# 设置测试源码目录
set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/testcase)

# 查找所有测试源文件
file(GLOB TEST_SOURCES ${TEST_SRC_DIR}/*.cpp)

# 添加主测试文件
set(TEST_MAIN ${CMAKE_CURRENT_SOURCE_DIR}/pbgz_test_main.cpp)

# 生成测试可执行文件
add_executable(pbgz_test ${TEST_MAIN} ${TEST_SOURCES})

# 包含主项目头文件
include_directories(${CMAKE_SOURCE_DIR}/code)
include_directories(${CMAKE_SOURCE_DIR}/code/log)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/release/jsoncpp/include)
include_directories(${CMAKE_SOURCE_DIR}/3rd_party/googletest/googletest/include)


set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/jsoncpp/lib/libjsoncpp.a
)

# 可选：添加gtest支持
target_link_libraries(pbgz_test PRIVATE gtest_main)

# 启用测试
enable_testing()
add_test(NAME pbgz_test COMMAND pbgz_test)

