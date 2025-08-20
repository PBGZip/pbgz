cmake_minimum_required(VERSION 3.10)
project(pbgz VERSION 2.0.0 LANGUAGES C CXX)

# Set C++ standard (add these lines right after project declaration)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add cmake directory to module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include build options
include(BuildOptions)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Define third-party library paths first (make these available globally)
set(THIRDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rd_party/release)

# Include third-party library headers globally
include_directories(
    ${THIRDPARTY_DIR}/jsoncpp/include
)

add_subdirectory(test)
add_subdirectory(${CMAKE_SOURCE_DIR}/3rd_party/googletest)

set(REQUIRED_LIBS
    ${THIRDPARTY_DIR}/jsoncpp/lib/libjsoncpp.a
)

set(MAIN_SOURCES
    code/log/logger.cpp
    code/pbgz_file_handler.cpp
    code/pbgz_file.cpp
    code/main.cpp
)

# Define main executable
add_executable(pbgz ${MAIN_SOURCES})


# Link to libraries - Add PRIVATE keyword to match style with CoreServices linkage
target_link_libraries(pbgz PRIVATE
    ${THIRDPARTY_DIR}/jsoncpp/lib/libjsoncpp.a
)
