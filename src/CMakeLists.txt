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
    ${THIRDPARTY_DIR}/md5sum
    ${CMAKE_CURRENT_SOURCE_DIR}
)

set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/release/zstd/lib/libzstd.a
)

add_subdirectory(coder)

# Find all source files
file(GLOB PBGZ_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/log/*.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp
)

# Directly compiled open source code
set(PBGZ_OPEN_SOURCES ${THIRDPARTY_DIR}/md5sum/md5sum.cpp)

set(MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp )

# Define main executable
add_executable(pbgz ${MAIN_SOURCES} ${PBGZ_SOURCES} ${PBGZ_OPEN_SOURCES})

# Link to libraries - Add PRIVATE keyword to match style with CoreServices linkage
target_link_libraries(pbgz PRIVATE
    pbgz_coder
    ${THIRDPARTY_DIR}/release/jsoncpp/lib/libjsoncpp.a
    ${THIRDPARTY_DIR}/release/zstd/lib/libzstd.a
    pthread
)

install(TARGETS pbgz DESTINATION bin)

